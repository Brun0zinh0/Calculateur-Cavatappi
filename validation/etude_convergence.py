"""etude_convergence.py - etude de convergence systematique du moteur Base.py.

Item 4.2 du plan de correction (audit 2026-08). Ce script balaie les
parametres de discretisation du modele TCPA (Base.py ; la version du
moteur effectivement utilisee est enregistree dans le champ meta du JSON
de sortie -- etude d'aout 2026 executee sur 2026.08.21-audit-phase4-14)
sur un protocole bloque de reference et
quantifie l'erreur relative de chaque reglage par rapport au reglage le
plus raffine de son balayage, afin de fixer des reglages de production
recommandes.

Protocole de reference (bloque)
-------------------------------
    - pression cyclique LINEAIRE (defaut projet depuis l'audit phase 2),
    - Pmax = 1.4 MPa, eps (pre-etirement) = 0.8, 3 cycles,
    - debit 10 mL/min, volume 1.50 mL -> periode 18 s, duree totale 54 s.

L'historique de pression est TOUJOURS genere explicitement par
Base.cyclic_pressure_history(dt=...) et passe via pressure_time /
pressure_MPa : le dt de la configuration est IGNORE des qu'une histoire
de pression est fournie, le pas de temps effectif du solveur est la
grille de pression (constat de l'audit 2026-08).

Balayages
---------
    A. n_layers dans {1, 2, 3, 4, 6, 8, 10}   (n_phi=16, dt=0.5)
    B. n_phi    dans {4, 8, 16, 32}           (n_layers=4, dt=0.5)
    C. dt       dans {2, 1, 0.5, 0.25, 0.1}   (n_layers=4, n_phi=16)
    D. integrateur paper_explicit vs exponential (cas commun n_layers=4,
       n_phi=16 ; dt=0.5 et dt=2, tous deux < 2*tau_min ~ 14.96 s,
       garde de stabilite d'Euler explicite de Base._validate_time_step)
    E. section_update_mode='updated' : n_layers {2, 4, 6} x dt {1, 0.5}
    F. mode suspendu (run_suspended_actuation, charge 1 N, rampe de
       pression simple 0 -> 1.4 MPa en 9 s puis maintien 60 s, dt=0.5) :
       n_layers {2, 4, 6}
    G. validation croisee des reglages recommandes (erreur COMBINEE des
       candidats "interactif" et "production" contre une reference tres
       raffinee n_layers=10, n_phi=16, dt=0.1)

Metriques par run
-----------------
    - pic de force du dernier cycle (force_total_mN, max),
    - vallee de force du dernier cycle (force_total_mN, min),
    - pic de couple d'actionnement du dernier cycle (|torque_act_microNm|),
    - temps de calcul mur (s), nombre de pas, residu max,
    - en suspendu : contraction max et finale (%), longueur finale (mm).

Robustesse aux timeouts
-----------------------
Les resultats sont ecrits de maniere INCREMENTALE (un run a la fois,
ecriture atomique de validation/out/convergence_results.json) : le script
peut etre relance autant de fois que necessaire, il ne re-execute que les
runs manquants. Options :

    python etude_convergence.py                # execute tous les runs manquants
    python etude_convergence.py --max-seconds 480   # s'arrete proprement avant le timeout
    python etude_convergence.py --only A_L4 C_dt0.1 # runs choisis
    python etude_convergence.py --list         # etat de la grille
    python etude_convergence.py --report       # tables markdown -> out/convergence_tables.md

Sorties (validation/out/)
-------------------------
    - convergence_results.json   : resultats bruts incrementaux,
    - convergence_tables.md      : tableaux + erreurs relatives par balayage.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))  # racine "alpha V3" pour importer Base
import Base  # noqa: E402

OUT_DIR = HERE / "out"
RESULTS_JSON = OUT_DIR / "convergence_results.json"
TABLES_MD = OUT_DIR / "convergence_tables.md"

# --- Protocole de reference -------------------------------------------------
PMAX_MPA = 1.4
EPS = 0.8
N_CYCLES = 3
FLOW_ML_MIN = 10.0
VOLUME_ML = 1.50
PERIOD_S = 2.0 * 60.0 * VOLUME_ML / FLOW_ML_MIN  # 18 s

# --- Cas suspendu -----------------------------------------------------------
SUSP_LOAD_N = 1.0
SUSP_RAMP_S = 9.0
SUSP_HOLD_S = 60.0
SUSP_DT = 0.5


def build_runs():
    """Construit la grille complete des runs (liste de dicts ordonnee)."""
    runs = []

    def add(run_id, sweep, **params):
        runs.append({"id": run_id, "sweep": sweep, **params})

    # A. n_layers (bloque, section fixe, exponential)
    for nl in (1, 2, 3, 4, 6, 8, 10):
        add(f"A_L{nl}", "A_n_layers", mode="blocked", n_layers=nl, n_phi=16,
            dt=0.5, integration="exponential", section_update_mode="fixed")

    # B. n_phi
    for np_ in (4, 8, 16, 32):
        add(f"B_P{np_}", "B_n_phi", mode="blocked", n_layers=4, n_phi=np_,
            dt=0.5, integration="exponential", section_update_mode="fixed")

    # C. dt (grille de pression = pas de temps effectif)
    for dt in (2.0, 1.0, 0.5, 0.25, 0.1):
        add(f"C_dt{dt:g}", "C_dt", mode="blocked", n_layers=4, n_phi=16,
            dt=dt, integration="exponential", section_update_mode="fixed")

    # D. integrateur explicite de l'article (dt admissible < 2*tau_min ~ 14.96 s)
    for dt in (0.5, 2.0):
        add(f"D_expl_dt{dt:g}", "D_integrateur", mode="blocked", n_layers=4,
            n_phi=16, dt=dt, integration="paper_explicit",
            section_update_mode="fixed")

    # E. section reactualisee (cout ~2x, audit item 3.x)
    for nl in (2, 4, 6):
        for dt in (1.0, 0.5):
            add(f"E_upd_L{nl}_dt{dt:g}", "E_updated", mode="blocked",
                n_layers=nl, n_phi=16, dt=dt, integration="exponential",
                section_update_mode="updated")

    # F. suspendu (charge 1 N, rampe simple)
    for nl in (2, 4, 6):
        add(f"F_susp_L{nl}", "F_suspendu", mode="suspended", n_layers=nl,
            n_phi=16, dt=SUSP_DT, integration="exponential",
            section_update_mode="fixed")

    # G. validation croisee des reglages recommandes : erreurs COMBINEES
    # mesurees contre une reference tres raffinee (et non par parametre isole).
    add("G_ref", "G_reco", mode="blocked", n_layers=10, n_phi=16, dt=0.1,
        integration="exponential", section_update_mode="fixed")
    add("G_interactif", "G_reco", mode="blocked", n_layers=3, n_phi=8, dt=0.5,
        integration="exponential", section_update_mode="fixed")
    add("G_production", "G_reco", mode="blocked", n_layers=6, n_phi=8, dt=0.1,
        integration="exponential", section_update_mode="fixed")

    return runs


# --- Metriques ---------------------------------------------------------------


def _last_cycle_mask(t):
    return np.asarray(t, dtype=float) >= float(t[-1]) - PERIOD_S - 1.0e-9


def blocked_metrics(arr):
    t = np.asarray(arr["time"], dtype=float)
    m = _last_cycle_mask(t)
    force = np.asarray(arr["force_total_mN"], dtype=float)[m]
    couple = np.asarray(arr["torque_act_microNm"], dtype=float)[m]
    return {
        "pic_force_mN": float(np.nanmax(force)),
        "vallee_force_mN": float(np.nanmin(force)),
        "pic_couple_microNm": float(np.nanmax(np.abs(couple))),
        "residu_max_Nmm": float(np.nanmax(np.abs(np.asarray(arr["residual"], dtype=float)))),
        "n_pas": int(len(t)),
    }


def suspended_metrics(arr):
    t = np.asarray(arr["time"], dtype=float)
    force = np.asarray(arr["force_total_mN"], dtype=float)
    couple = np.asarray(arr["torque_act_microNm"], dtype=float)
    contr = np.asarray(arr["free_actuation_percent"], dtype=float)
    return {
        "pic_force_mN": float(np.nanmax(force)),
        "vallee_force_mN": float(np.nanmin(force)),
        "pic_couple_microNm": float(np.nanmax(np.abs(couple))),
        "contraction_max_pct": float(np.nanmax(contr)),
        "contraction_finale_pct": float(contr[-1]),
        "longueur_finale_mm": float(np.asarray(arr["axial_length_mm"], dtype=float)[-1]),
        "residu_max_Nmm": float(np.nanmax(np.abs(np.asarray(arr["residual"], dtype=float)))),
        "n_pas": int(len(t)),
    }


# --- Execution d'un run -------------------------------------------------------


def execute_run(run):
    """Execute un run et retourne (metriques + temps de calcul)."""
    geom = Base.default_geometry_params(section_update_mode=run["section_update_mode"])
    if run["mode"] == "blocked":
        t, p = Base.cyclic_pressure_history(
            n_cycles=N_CYCLES,
            Pmax=PMAX_MPA,
            flow_rate_mL_min=FLOW_ML_MIN,
            volume_mL=VOLUME_ML,
            dt=run["dt"],
            nonlinear=False,  # profil lineaire (defaut projet)
        )
        t0 = time.perf_counter()
        _, arr = Base.run_blocked_actuation(
            eps=EPS,
            Pmax=PMAX_MPA,
            n_cycles=N_CYCLES,
            n_layers=run["n_layers"],
            n_phi=run["n_phi"],
            integration=run["integration"],
            geom=geom,
            pressure_time=t,
            pressure_MPa=p,
        )
        elapsed = time.perf_counter() - t0
        metrics = blocked_metrics(arr)
    elif run["mode"] == "suspended":
        t, p = Base.ramp_hold_pressure_history(
            P_hold=PMAX_MPA,
            ramp_time=SUSP_RAMP_S,
            hold_time=SUSP_HOLD_S,
            dt=run["dt"],
        )
        t0 = time.perf_counter()
        _, arr = Base.run_suspended_actuation(
            eps=EPS,
            load_N=SUSP_LOAD_N,
            n_layers=run["n_layers"],
            n_phi=run["n_phi"],
            integration=run["integration"],
            geom=geom,
            pressure_time=t,
            pressure_MPa=p,
        )
        elapsed = time.perf_counter() - t0
        metrics = suspended_metrics(arr)
    else:
        raise ValueError(f"mode inconnu: {run['mode']}")
    metrics["temps_s"] = round(elapsed, 3)
    return metrics


# --- Persistance incrementale --------------------------------------------------


def load_results():
    if RESULTS_JSON.exists():
        with open(RESULTS_JSON, "r", encoding="utf-8") as f:
            return json.load(f)
    return {
        "meta": {
            "model_version": Base.MODEL_VERSION,
            "protocole": {
                "Pmax_MPa": PMAX_MPA,
                "eps": EPS,
                "n_cycles": N_CYCLES,
                "profil_pression": "lineaire",
                "periode_s": PERIOD_S,
                "suspendu": {
                    "load_N": SUSP_LOAD_N,
                    "rampe_s": SUSP_RAMP_S,
                    "maintien_s": SUSP_HOLD_S,
                    "dt": SUSP_DT,
                },
            },
        },
        "runs": {},
    }


def save_results(results):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tmp = RESULTS_JSON.with_suffix(".json.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, RESULTS_JSON)


# --- Rapport -----------------------------------------------------------------


def _rel_err(value, ref):
    if ref == 0.0:
        return float("nan")
    return 100.0 * (value - ref) / abs(ref)


def _fmt(x, nd=1):
    return f"{x:.{nd}f}"


def make_report(results):
    """Genere les tableaux markdown (valeurs + erreur relative vs le plus raffine)."""
    runs_def = build_runs()
    data = results["runs"]
    lines = [
        "# Tableaux de convergence (bruts)",
        "",
        f"Moteur : Base.py {results['meta']['model_version']} — protocole bloque : "
        f"pression cyclique lineaire, Pmax {PMAX_MPA} MPa, eps {EPS}, {N_CYCLES} cycles "
        f"(periode {PERIOD_S:g} s).",
        "",
    ]
    metric_keys = ("pic_force_mN", "pic_couple_microNm", "vallee_force_mN")
    refs = {
        "A_n_layers": "A_L10",
        "B_n_phi": "B_P32",
        "C_dt": "C_dt0.1",
        "E_updated": "E_upd_L6_dt0.5",
        "F_suspendu": "F_susp_L6",
        "G_reco": "G_ref",
    }
    sweeps = {}
    for run in runs_def:
        sweeps.setdefault(run["sweep"], []).append(run)

    for sweep, sruns in sweeps.items():
        lines.append(f"## Balayage {sweep}")
        lines.append("")
        ref_id = refs.get(sweep)
        ref = data.get(ref_id, {}) if ref_id else {}
        if sweep == "F_suspendu":
            header = ("| run | n_layers | dt (s) | contraction max (%) | err. (%) | "
                      "contraction finale (%) | pic force (mN) | temps (s) |")
            sep = "|---|---|---|---|---|---|---|---|"
            lines += [header, sep]
            for run in sruns:
                r = data.get(run["id"])
                if not r:
                    lines.append(f"| {run['id']} | {run['n_layers']} | {run['dt']} | (non exécuté) | | | | |")
                    continue
                err = _rel_err(r["contraction_max_pct"], ref.get("contraction_max_pct", 0.0)) if ref else float("nan")
                lines.append(
                    f"| {run['id']} | {run['n_layers']} | {run['dt']} | "
                    f"{_fmt(r['contraction_max_pct'], 3)} | {_fmt(err, 3)} | "
                    f"{_fmt(r['contraction_finale_pct'], 3)} | {_fmt(r['pic_force_mN'])} | "
                    f"{_fmt(r['temps_s'])} |"
                )
        else:
            header = ("| run | n_layers | n_phi | dt (s) | integ. | section | "
                      "pic force (mN) | err. (%) | pic couple (uN.m) | err. (%) | "
                      "vallee force (mN) | err. (%) | temps (s) |")
            sep = "|" + "---|" * 13
            lines += [header, sep]
            for run in sruns:
                r = data.get(run["id"])
                base = (f"| {run['id']} | {run['n_layers']} | {run['n_phi']} | {run['dt']:g} | "
                        f"{run['integration']} | {run['section_update_mode']} | ")
                if not r:
                    lines.append(base + "(non exécuté) | | | | | | |")
                    continue
                cells = []
                for key in metric_keys:
                    cells.append(_fmt(r[key], 2))
                    if ref and ref_id != run["id"]:
                        cells.append(_fmt(_rel_err(r[key], ref[key]), 3))
                    elif ref_id == run["id"]:
                        cells.append("réf.")
                    else:
                        cells.append("—")
                cells.append(_fmt(r["temps_s"]))
                lines.append(base + " | ".join(cells) + " |")
        lines.append("")

    # Comparaison integrateurs : paper_explicit vs exponential a dt identique.
    lines.append("## Comparaison intégrateurs (écart paper_explicit vs exponential, dt identique)")
    lines.append("")
    lines.append("| dt (s) | Δ pic force (%) | Δ pic couple (%) | Δ vallée force (%) |")
    lines.append("|---|---|---|---|")
    for dt in (0.5, 2.0):
        expl = data.get(f"D_expl_dt{dt:g}")
        expo = data.get(f"C_dt{dt:g}")
        if expl and expo:
            lines.append(
                f"| {dt:g} | {_fmt(_rel_err(expl['pic_force_mN'], expo['pic_force_mN']), 4)} | "
                f"{_fmt(_rel_err(expl['pic_couple_microNm'], expo['pic_couple_microNm']), 4)} | "
                f"{_fmt(_rel_err(expl['vallee_force_mN'], expo['vallee_force_mN']), 4)} |"
            )
    lines.append("")

    # Comparaison fixed vs updated a reglages identiques (dt=0.5).
    lines.append("## Écart section fixe → réactualisée (dt = 0.5 s, n_phi = 16)")
    lines.append("")
    lines.append("| n_layers | Δ pic force (%) | Δ pic couple (%) | ratio temps updated/fixed |")
    lines.append("|---|---|---|---|")
    for nl in (2, 4, 6):
        fixed = data.get(f"A_L{nl}")
        upd = data.get(f"E_upd_L{nl}_dt0.5")
        if fixed and upd:
            ratio = upd["temps_s"] / fixed["temps_s"] if fixed["temps_s"] > 0 else float("nan")
            lines.append(
                f"| {nl} | {_fmt(_rel_err(upd['pic_force_mN'], fixed['pic_force_mN']), 3)} | "
                f"{_fmt(_rel_err(upd['pic_couple_microNm'], fixed['pic_couple_microNm']), 3)} | "
                f"{_fmt(ratio, 2)} |"
            )
    lines.append("")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_MD.write_text("\n".join(lines), encoding="utf-8")
    return TABLES_MD


# --- Point d'entree -----------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description="Étude de convergence Base.py (item 4.2, audit 2026-08).")
    parser.add_argument("--only", nargs="*", default=None, help="identifiants de runs à exécuter")
    parser.add_argument("--max-seconds", type=float, default=None,
                        help="budget temps : s'arrête proprement une fois dépassé")
    parser.add_argument("--list", action="store_true", help="affiche l'état de la grille")
    parser.add_argument("--report", action="store_true", help="génère out/convergence_tables.md")
    args = parser.parse_args(argv)

    runs = build_runs()
    results = load_results()

    if args.list:
        for run in runs:
            status = "fait" if run["id"] in results["runs"] else "à faire"
            print(f"  {run['id']:<18} [{run['sweep']}] {status}")
        done = sum(1 for r in runs if r["id"] in results["runs"])
        print(f"{done}/{len(runs)} runs faits.")
        return 0

    if args.report:
        path = make_report(results)
        print(f"Tables écrites : {path}")
        return 0

    wanted = set(args.only) if args.only else None
    start = time.perf_counter()
    executed = 0
    for run in runs:
        if wanted is not None and run["id"] not in wanted:
            continue
        if wanted is None and run["id"] in results["runs"]:
            continue  # déjà fait
        if args.max_seconds is not None and time.perf_counter() - start > args.max_seconds:
            print(f"Budget temps atteint ({args.max_seconds:g} s) — arrêt propre.")
            break
        print(f"-> {run['id']} ...", flush=True)
        metrics = execute_run(run)
        results["runs"][run["id"]] = {**{k: v for k, v in run.items() if k != "id"}, **metrics}
        save_results(results)  # flush incrémental : survit aux timeouts
        executed += 1
        print(f"   ok ({metrics['temps_s']:.1f} s) : "
              + ", ".join(f"{k}={v:.4g}" for k, v in metrics.items()
                          if isinstance(v, float) and k != "temps_s"), flush=True)

    done = sum(1 for r in runs if r["id"] in results["runs"])
    print(f"{executed} run(s) exécuté(s) ; grille : {done}/{len(runs)}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
