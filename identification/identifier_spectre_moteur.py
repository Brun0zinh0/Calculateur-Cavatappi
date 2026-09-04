# -*- coding: utf-8 -*-
"""Identification du spectre de Maxwell EN BOUCLE FERMÉE à travers le moteur.

Alpha V4-6 (proposition retenue n° 5 de la contre-expertise). Le spectre
identifié directement sur F(t)/F(0) (identification/identifier_materiau.py,
essais/identifier_spectre.py) suppose F(t)/F(0) = E(t)/E(0) ; réinjecté dans
la chaîne (projection hélicoïdale, partage tube/nylon, relaxation axiale
seule), il produit une relaxation de force ~14 fois plus faible que la mesure
(rapport 7.8.7). Ici les {E_i, eta_i} sont ajustés pour que la force SIMULÉE
par le moteur en maintien bloqué à P = 0 reproduise la force mesurée.

Fonction coût : écart entre F_sim(t)/F_sim(0) et F_mes(t)/F_mes(0) (forme
normalisée — le niveau absolu de précontrainte, biais B1 du rapport, ne pollue
pas le spectre ; contre-expertise, correction 1). ΣE est maintenu à sa valeur
d'entrée : seules les fractions et les temps bougent.

Usage :
    python identifier_spectre_moteur.py essai_10N.csv --fiche fiche.json --eps 1.0
    python identifier_spectre_moteur.py essai_10N.csv --branches 2 --hold 600

Le CSV doit contenir time_s et force_unfiltered_mN (ou force_mN). La fiche
JSON est un dictionnaire de réglages (mêmes clés que parametres.DEFAULT_SETTINGS,
par exemple identification/fiches/muscle_J_eps10.json).

ATTENTION : chaque évaluation du coût est une simulation complète. Avec les
réglages interactifs (L3/nφ8) et un maintien de 600 s, compter 20-40 s par
évaluation et 30-80 évaluations : 15-45 min par muscle. Réduire --hold,
--dt ou le maillage pour dégrossir.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time as _time

import numpy as np
from scipy.optimize import least_squares

APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if APP_DIR not in sys.path:
    sys.path.insert(0, APP_DIR)

import Base  # noqa: E402
import parametres  # noqa: E402


def load_force_csv(path: str):
    """Lit time_s et force_unfiltered_mN (ou force_mN). Délimiteur ',' ';' ou
    tabulation détecté automatiquement, virgule décimale tolérée (export Excel
    en locale française). Erreurs explicites sur fichier vide ou colonne absente
    (contre-expertise C14)."""
    import io

    with open(path, "rb") as fh:
        raw = fh.read()
    if not raw.strip():
        raise ValueError(f"CSV vide : {path}")
    text = raw.decode("utf-8-sig", errors="replace")
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)
    fields = [str(name).strip() for name in (reader.fieldnames or [])]
    reader.fieldnames = fields
    if "time_s" not in fields:
        raise ValueError(f"colonne time_s absente ; colonnes disponibles : {fields}")
    if "force_unfiltered_mN" in fields:
        col = "force_unfiltered_mN"
    elif "force_mN" in fields:
        col = "force_mN"
    else:
        raise ValueError(f"colonne force_unfiltered_mN ou force_mN absente ; colonnes : {fields}")

    def number(value) -> float:
        return float(str(value).strip().replace(",", "."))

    t, f = [], []
    for row in reader:
        if not row.get("time_s") or not row.get(col):
            continue
        t.append(number(row["time_s"]))
        f.append(number(row[col]))
    if len(t) < 2:
        raise ValueError(f"moins de deux lignes numériques exploitables dans {path}")
    return np.asarray(t), np.asarray(f)


def simulate_hold(settings: dict, hold_s: float, dt: float, prestretch_rate_mm_min: float):
    """Précontrainte viscoélastique rapide puis maintien bloqué à P = 0."""
    cfg = parametres.build_config(settings)
    disc = Base.default_discretization(
        n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05)
    )
    model = Base.TCPAMaxwellBlockedModel(
        mat=cfg.mat, geom=cfg.geom, disc=disc, integration="exponential",
        prestrain_reference_mode="viscoelastic_history",
    )
    model.prestretch_to(cfg.eps, strain_rate_mm_min=prestretch_rate_mm_min)
    t_start = model.helix.time
    model.step(0.0, 0.0, h_target=model.h_blocked)
    model.lock_blocked_series_reference()
    i0 = len(model.history) - 1
    t_hold = np.arange(0.0, hold_s + 0.5 * dt, dt)
    model.run_pressure_history(t_start + t_hold, np.zeros_like(t_hold))
    arr = model.history_arrays()
    return np.asarray(arr["time"][i0:]) - t_start, np.asarray(arr["force_mN"][i0:])


def settings_with_spectrum(base: dict, fractions: np.ndarray, taus: np.ndarray, sigma_total: float) -> dict:
    s = dict(base)
    E = fractions * sigma_total
    E0 = sigma_total - float(E.sum())
    s["maxwell_E0_mpa"] = max(E0, 1.0e-6)
    for k in range(3):
        Ek = float(E[k]) if k < len(E) else 0.0
        tk = float(taus[k]) if k < len(taus) else 1.0
        s[f"maxwell_E{k + 1}_mpa"] = Ek
        s[f"maxwell_eta{k + 1}_mpa_s"] = max(Ek * tk, 1.0e-9)
    return s


def identify(
    t_meas, f_meas, base_settings: dict, n_branches: int = 2, hold_s: float = 600.0,
    dt: float = 2.0, prestretch_rate_mm_min: float = 20.0, sigma_total: float | None = None,
    verbose: bool = True, t0: float | None = None, min_points: int = 10,
):
    """t0 : instant (s, temps brut du CSV) du blocage — début de la fenêtre
    ajustée [t0, t0 + hold_s] ; défaut : premier échantillon. Le spectre
    identifié est conditionnel à prestretch_rate_mm_min (relaxation consommée
    pendant la rampe) : défaut 20 mm/min = Base.prestretch_to (C12)."""
    if sigma_total is None:
        sigma_total = sum(float(base_settings[k]) for k in ("maxwell_E0_mpa", "maxwell_E1_mpa", "maxwell_E2_mpa", "maxwell_E3_mpa"))
    t_meas = np.asarray(t_meas, dtype=float)
    f_meas = np.asarray(f_meas, dtype=float)
    if t0 is None:
        t0 = float(t_meas[0])
    t_rel = t_meas - float(t0)
    mask = (t_rel >= 0.0) & (t_rel <= hold_s)
    n_in = int(np.sum(mask))
    if n_in < min_points:
        raise ValueError(
            f"fenêtre de maintien [t0 = {t0:.1f} s, t0 + {hold_s:.0f} s] : {n_in} point(s) du CSV "
            f"(temps brut de {t_meas[0]:.1f} à {t_meas[-1]:.1f} s) ; précisez --t0 ou --hold."
        )
    t_fit = t_rel[mask]
    n_ref = max(1, min(3, int(np.sum(t_fit <= t_fit[0] + 1.0))))
    y_meas = f_meas[mask] / float(np.mean(f_meas[mask][:n_ref]))
    # Paramètres : log-fractions et log-taus (positivité garantie)
    f0 = np.array([0.05, 0.08, 0.05][:n_branches])
    tau0 = np.array([20.0, 400.0, 3000.0][:n_branches])
    x0 = np.concatenate([np.log(f0), np.log(tau0)])
    lower = np.concatenate([np.full(n_branches, np.log(1.0e-4)), np.log([2.0, 60.0, 600.0][:n_branches])])
    upper = np.concatenate([np.full(n_branches, np.log(0.6)), np.log([120.0, 3000.0, 30000.0][:n_branches])])
    evaluations = []

    def residuals(x):
        fr = np.exp(x[:n_branches])
        ta = np.exp(x[n_branches:])
        if fr.sum() >= 0.95:
            return np.full_like(y_meas, 10.0)
        try:
            s = settings_with_spectrum(base_settings, fr, ta, sigma_total)
            s["dt"] = dt
            t_sim, F_sim = simulate_hold(s, hold_s, dt, prestretch_rate_mm_min)
            F_on = np.interp(t_fit, t_sim, F_sim)
            y_sim = F_on / float(np.mean(F_on[:n_ref]))
        except Exception as exc:  # état inadmissible : pénalité
            if verbose:
                print("  evaluation rejetee :", exc)
            return np.full_like(y_meas, 10.0)
        r = y_sim - y_meas
        evaluations.append((float(np.sqrt(np.mean(r * r))), fr.copy(), ta.copy()))
        if verbose:
            print(f"  eval {len(evaluations):3d} : RMS {evaluations[-1][0]:.5f} | fractions {np.round(fr, 4)} | taus {np.round(ta, 1)}", flush=True)
        return r

    t0 = _time.time()
    res = least_squares(residuals, x0, bounds=(lower, upper), xtol=1.0e-4, ftol=1.0e-5, max_nfev=80, diff_step=0.05)
    fr = np.exp(res.x[:n_branches])
    ta = np.exp(res.x[n_branches:])
    E = fr * sigma_total
    out = {
        "success": bool(res.success),
        "rms_normalized": float(np.sqrt(np.mean(res.fun ** 2))),
        "n_evaluations": len(evaluations),
        "elapsed_s": _time.time() - t0,
        "sigma_total_mpa": float(sigma_total),
        "prestretch_rate_mm_min": float(prestretch_rate_mm_min),
        "t0_s": float(t0),
        "hold_s": float(hold_s),
        "E0_mpa": float(sigma_total - E.sum()),
        "E0_over_sigma": float(1.0 - fr.sum()),
        "branches": [
            {"E_mpa": float(E[k]), "tau_s": float(ta[k]), "eta_mpa_s": float(E[k] * ta[k]), "fraction": float(fr[k])}
            for k in range(n_branches)
        ],
    }
    return out


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description="Identification du spectre de Maxwell en boucle fermee a travers le moteur.")
    ap.add_argument("csv", help="Essai de relaxation à P = 0 (time_s, force_unfiltered_mN)")
    ap.add_argument("--fiche", help="Fiche de réglages JSON (géométrie du muscle)")
    ap.add_argument("--eps", type=float, default=None, help="Pré-étirement (surcharge la fiche)")
    ap.add_argument("--branches", type=int, default=2, choices=(1, 2, 3))
    ap.add_argument("--hold", type=float, default=600.0, help="Durée ajustée (s)")
    ap.add_argument("--dt", type=float, default=2.0)
    ap.add_argument("--rate", type=float, default=20.0,
                    help="Vitesse de pre-etirement de l'essai (mm/min) ; defaut 20 = Base.prestretch_to. "
                         "Le spectre identifie est conditionnel a cette vitesse (relaxation consommee en rampe).")
    ap.add_argument("--t0", type=float, default=None,
                    help="Instant du blocage dans le temps brut du CSV (s) ; defaut : premier echantillon")
    ap.add_argument("--sigma-total", type=float, default=None,
                    help="Somme des modules E maintenue (MPa) ; defaut : somme des modules de la fiche")
    ap.add_argument("--out", default=None, help="Fichier JSON de sortie")
    args = ap.parse_args()

    settings = dict(parametres.DEFAULT_SETTINGS)
    if args.fiche:
        with open(args.fiche, encoding="utf-8") as fh:
            fiche = json.load(fh)
        settings.update({k: v for k, v in fiche.items() if k in settings})
    if args.eps is not None:
        settings["eps"] = float(args.eps)
    settings.update(dict(n_layers=3, n_phi=8, pre_steps=24, prestrain_reference_mode="viscoelastic_history",
                         maxwell_anisotropy_mode="axial_test_only", integration="exponential"))
    t, f = load_force_csv(args.csv)
    print(f"Base {Base.MODEL_VERSION} — {len(t)} points, {t[-1]:.0f} s ; eps {settings['eps']} ; {args.branches} branche(s)")
    out = identify(t, f, settings, n_branches=args.branches, hold_s=args.hold, dt=args.dt,
                   prestretch_rate_mm_min=args.rate, sigma_total=args.sigma_total, t0=args.t0)
    print(json.dumps(out, indent=1, ensure_ascii=False))
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, ensure_ascii=False)
        print("ecrit", args.out)


if __name__ == "__main__":
    main()
