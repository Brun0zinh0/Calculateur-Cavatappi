# -*- coding: utf-8 -*-
"""Identification matériau pour les muscles de l'équipe (audit 2026-08, item 3.3).

Pipeline en trois volets, honnête sur ce qui est identifiable :

A) PALIERS — détection des paliers de pression dans un essai bloqué mesuré,
   puis ajustement multi-exponentiel de la force pendant chaque palier :
   -> constantes de temps tau_k et fractions de relaxation, SANS géométrie.
   C'est la partie transférable du spectre. Un palier où la force monte
   (fuite, dérive capteur, raidissement) est signalé, pas forcé.

B) ÉCHELLE — simulation de l'essai complet avec le moteur Base (pression
   mesurée injectée, géométrie de la fiche muscle ou défauts de l'article)
   et ajustement en forme fermée d'un facteur d'échelle unique sur la force
   d'actionnement : scale = <F_exp_act . F_sim_act> / <F_sim_act²>.
   Tant que la fiche géométrique (essai E1) manque, ce facteur absorbe
   l'incertitude de géométrie — les modules absolus ne sont PAS identifiés.

C) COURBURE — verdict « spectre mal paramétré vs limite structurelle » :
   ajustement F_act = a.P + b.P² sur la première montée en pression, mesuré
   et simulé. Si la courbure mesurée reste hors de portée du modèle quelle
   que soit l'échelle, elle est structurelle (dead-band / raidissement en
   pression, hors du cadre de l'article) et doit être documentée comme
   limite du modèle plutôt que compensée par un fit.

Usage :

    python identifier_materiau.py --essai <fichier.csv|.xlsx> [--fiche fiche.json]
        [--mode elastic_tk_reference|viscoelastic_history] [--sans-simulation]

Entrées : CSV (séparateur , ou ; ; virgule décimale acceptée) avec des
colonnes temps/pression/force (unités inférées des en-têtes : s, bar/psi/MPa,
mN/N/g), ou xlsx (première feuille, mêmes conventions). Fiche muscle : JSON
avec les clés de géométrie de `parametres.DEFAULT_SETTINGS` (rout_mm, rin_mm,
nylon_diameter_mm, rho0_mm, alpha0_deg, theta_f_deg, initial_length_mm,
uncoiled_length_mm, eps) — les clés absentes prennent les valeurs de
l'article, avec avertissement (voir fiche_muscle_exemple.json).

Sorties : rapport console, `out/identification_<essai>.json` et
`out/identification_<essai>.png` à côté du présent script.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
APP_DIR = SCRIPT_DIR.parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import Base  # noqa: E402
import parametres  # noqa: E402
import pression  # noqa: E402

OUT_DIR = SCRIPT_DIR / "out"

FICHE_GEOMETRY_KEYS = (
    "rout_mm",
    "rin_mm",
    "nylon_diameter_mm",
    "rho0_mm",
    "alpha0_deg",
    "theta_f_deg",
    "initial_length_mm",
    "uncoiled_length_mm",
    "eps",
)


# ---------------------------------------------------------------------------
# Chargement des données
# ---------------------------------------------------------------------------


def load_fiche(path: Path | None) -> tuple[dict, list[str]]:
    """Retourne (surcharges de réglages, avertissements)."""
    overrides: dict = {}
    warnings: list[str] = []
    if path is None:
        warnings.append(
            "Aucune fiche muscle fournie : géométrie = spécimen de l'article "
            "(essai E1 du plan d'essais manquant). Les modules absolus ne sont "
            "pas identifiables ; seul un facteur d'échelle global est ajusté."
        )
        return overrides, warnings
    with open(path, encoding="utf-8") as fh:
        fiche = json.load(fh)
    for key in FICHE_GEOMETRY_KEYS:
        if key in fiche:
            overrides[key] = float(fiche[key])
        else:
            warnings.append(f"Fiche muscle : clé '{key}' absente, valeur de l'article utilisée.")
    for key in fiche:
        if key not in FICHE_GEOMETRY_KEYS:
            warnings.append(f"Fiche muscle : clé '{key}' ignorée (inconnue).")
    return overrides, warnings


def _columns_from_xlsx(path: Path) -> dict[str, np.ndarray]:
    import pandas as pd

    frame = pd.read_excel(path)
    columns: dict[str, np.ndarray] = {}
    for name in frame.columns:
        values = pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(values).sum() >= 3:
            columns[str(name)] = values
    if not columns:
        raise ValueError("Aucune colonne numérique exploitable dans le xlsx.")
    return columns


def _pick_column(columns: dict[str, np.ndarray], keywords: tuple[str, ...]) -> str | None:
    for name in columns:
        lowered = name.lower()
        if any(word in lowered for word in keywords):
            return name
    return None


def load_essai(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Retourne (t_s, P_MPa, F_mN, notes)."""
    notes: list[str] = []
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        columns = _columns_from_xlsx(path)
    else:
        columns = pression.parse_uploaded_numeric_csv(path.read_bytes())

    time_col = _pick_column(columns, ("time", "temps", "t_s", "t ("))
    pressure_col = _pick_column(columns, ("press", "pression", "bar", "psi", "mpa"))
    force_col = _pick_column(columns, ("force", "load", "charge", "mn", "poids"))
    if pressure_col is None or force_col is None:
        raise ValueError(
            f"Colonnes pression/force introuvables parmi {sorted(columns)} — "
            "renommer les en-têtes ou préciser les unités."
        )
    if time_col is None:
        notes.append(
            "Pas de colonne de temps : index d'échantillon utilisé avec un pas de "
            "1 s PAR HYPOTHÈSE — les constantes de temps identifiées ne sont "
            "fiables qu'à cette hypothèse près (cadence réelle à documenter, essai E3)."
        )
        t = np.arange(len(columns[pressure_col]), dtype=float)
    else:
        t = np.asarray(columns[time_col], dtype=float)
        unit_t = pression.infer_time_unit(time_col)
        if unit_t == "ms":
            t = t / 1000.0
        t = t - t[0]

    unit_p = pression.infer_pressure_unit(pressure_col)
    factor_p = {"MPa": 1.0, "bar": 0.1, "kPa": 1.0e-3, "psi": parametres.PSI_TO_MPA}[unit_p]
    P = np.asarray(columns[pressure_col], dtype=float) * factor_p
    notes.append(f"Pression : colonne '{pressure_col}' lue en {unit_p}.")

    unit_f = pression.infer_force_unit(force_col)
    factor_f = {"mN": 1.0, "N": 1000.0, "g": 9.80665, "kg": 9806.65}[unit_f]
    F = np.asarray(columns[force_col], dtype=float) * factor_f
    notes.append(f"Force : colonne '{force_col}' lue en {unit_f}.")

    valid = np.isfinite(t) & np.isfinite(P) & np.isfinite(F)
    t, P, F = t[valid], P[valid], F[valid]
    keep = np.concatenate(([True], np.diff(t) > 0.0))
    return t[keep], P[keep], F[keep], notes


# ---------------------------------------------------------------------------
# A) Paliers et ajustement multi-exponentiel
# ---------------------------------------------------------------------------


def detect_plateaus(
    t: np.ndarray,
    P: np.ndarray,
    min_duration_s: float = 12.0,
    pressure_tolerance_mpa: float = 0.01,
    min_pressure_mpa: float = 0.02,
) -> list[tuple[int, int]]:
    """Segments [i0, i1] où la pression reste constante à la tolérance près."""
    plateaus: list[tuple[int, int]] = []
    i = 0
    n = len(t)
    while i < n - 1:
        j = i
        while j + 1 < n and abs(P[j + 1] - P[i]) <= pressure_tolerance_mpa:
            j += 1
        if t[j] - t[i] >= min_duration_s and np.median(P[i : j + 1]) >= min_pressure_mpa:
            plateaus.append((i, j))
        i = j + 1 if j > i else i + 1
    return plateaus


def fit_plateau_relaxation(t: np.ndarray, F: np.ndarray) -> dict:
    """Ajuste F(t) = F_inf + A1 exp(-t/tau1) [+ A2 exp(-t/tau2)] sur un palier.

    Retourne taus, amplitudes, fraction de relaxation signée
    (positive = la force décroît, comme un Maxwell ; négative = elle monte).
    """
    from scipy.optimize import curve_fit

    tt = t - t[0]
    F0 = float(F[0])
    span = float(F[0] - F[-1])
    duration = float(tt[-1])
    result: dict = {
        "duration_s": duration,
        "force_start_mN": F0,
        "force_end_mN": float(F[-1]),
        "relative_change": float((F[-1] - F0) / max(abs(F0), 1e-9)),
    }
    if span <= 0.0:
        result["verdict"] = (
            "force croissante ou constante pendant le palier : aucune relaxation de "
            "Maxwell identifiable (fuite, dérive capteur ou raidissement — essai E3 requis)"
        )
        return result

    def one_exp(x, f_inf, a1, tau1):
        return f_inf + a1 * np.exp(-x / tau1)

    def two_exp(x, f_inf, a1, tau1, a2, tau2):
        return f_inf + a1 * np.exp(-x / tau1) + a2 * np.exp(-x / tau2)

    best = None
    for model, p0, bounds, labels in (
        (one_exp, [F[-1], span, duration / 3.0], ([-np.inf, 0.0, 0.1], [np.inf, np.inf, 10.0 * duration]), ("tau1",)),
        (
            two_exp,
            [F[-1], 0.7 * span, duration / 6.0, 0.3 * span, duration],
            ([-np.inf, 0.0, 0.1, 0.0, 0.1], [np.inf, np.inf, 10.0 * duration, np.inf, 10.0 * duration]),
            ("tau1", "tau2"),
        ),
    ):
        try:
            popt, _ = curve_fit(model, tt, F, p0=p0, bounds=bounds, maxfev=20000)
            residual = float(np.sqrt(np.mean((model(tt, *popt) - F) ** 2)))
            if best is None or residual < 0.98 * best["rmse_mN"]:
                taus = sorted(float(x) for x in popt[2::2])
                amplitudes = [float(x) for x in popt[1::2]]
                best = {
                    "rmse_mN": residual,
                    "taus_s": taus,
                    "amplitudes_mN": amplitudes,
                    "force_infinite_mN": float(popt[0]),
                    "relaxation_fraction": float(sum(amplitudes) / max(abs(F0), 1e-9)),
                    "n_terms": len(labels),
                }
        except Exception:
            continue
    if best is None:
        result["verdict"] = "ajustement exponentiel non convergé"
    else:
        result.update(best)
        result["verdict"] = "relaxation identifiée"
    return result


# ---------------------------------------------------------------------------
# B) Simulation et facteur d'échelle
# ---------------------------------------------------------------------------


def simulate_essai(
    t: np.ndarray,
    P: np.ndarray,
    overrides: dict,
    mode: str = "elastic_tk_reference",
    target_dt_s: float = 1.0,
    n_layers: int = 2,
    n_phi: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Simule l'essai bloqué avec la pression mesurée (rééchantillonnée)."""
    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update(
        {
            "n_layers": n_layers,
            "n_phi": n_phi,
            "pre_steps": 12,
            "prestrain_reference_mode": mode,
            "p_max_mpa": min(1.5, float(np.max(P)) + 0.05),
        }
    )
    settings.update(overrides)
    config = parametres.build_config(settings)
    t_ds = np.arange(0.0, float(t[-1]), target_dt_s)
    P_ds = np.clip(np.interp(t_ds, t, P), 0.0, None)
    _, data = Base.run_blocked_actuation(config, pressure_time=t_ds, pressure_MPa=P_ds)
    return np.asarray(data["time"], float), np.asarray(data["force_total_mN"], float)


def closed_form_scale(F_exp_act: np.ndarray, F_sim_act: np.ndarray) -> float:
    denominator = float(np.dot(F_sim_act, F_sim_act))
    return float(np.dot(F_exp_act, F_sim_act) / denominator) if denominator > 0 else float("nan")


# ---------------------------------------------------------------------------
# C) Courbure de la première montée
# ---------------------------------------------------------------------------


def first_ramp_curvature(t: np.ndarray, P: np.ndarray, F_act: np.ndarray) -> dict | None:
    """Ajuste F_act = a.P + b.P² sur la première montée monotone en pression."""
    i_start = int(np.argmax(P > 0.02)) if np.any(P > 0.02) else None
    if i_start is None:
        return None
    i = i_start
    peak = P[i]
    i_end = i
    while i_end + 1 < len(P) and P[i_end + 1] >= peak - 0.005:
        i_end += 1
        peak = max(peak, P[i_end])
        if P[i_end] < 0.8 * peak:
            break
    segment = slice(i_start, i_end + 1)
    Ps, Fs = P[segment], F_act[segment]
    if len(Ps) < 8 or np.ptp(Ps) < 0.1:
        return None
    design = np.column_stack([Ps, Ps**2])
    coeffs, *_ = np.linalg.lstsq(design, Fs, rcond=None)
    a, b = (float(x) for x in coeffs)
    return {
        "a_mN_per_MPa": a,
        "b_mN_per_MPa2": b,
        "curvature_ratio_b_over_a": float(b / a) if abs(a) > 1e-9 else float("nan"),
        "pressure_span_MPa": float(np.ptp(Ps)),
        "n_points": int(len(Ps)),
    }


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def identify(essai_path: Path, fiche_path: Path | None, mode: str, run_simulation: bool = True) -> dict:
    overrides, warnings = load_fiche(fiche_path)
    t, P, F, notes = load_essai(essai_path)
    report: dict = {
        "essai": str(essai_path),
        "fiche": str(fiche_path) if fiche_path else None,
        "mode": mode,
        "notes": notes + warnings,
        "n_samples": int(len(t)),
        "duration_s": float(t[-1]),
        "pressure_max_MPa": float(np.max(P)),
        "preload_mN": float(np.median(F[P < 0.01][:50])) if np.any(P < 0.01) else float(F[0]),
    }

    plateaus = detect_plateaus(t, P)
    report["plateaus"] = []
    for i0, i1 in plateaus:
        fit = fit_plateau_relaxation(t[i0 : i1 + 1], F[i0 : i1 + 1])
        fit["pressure_MPa"] = float(np.median(P[i0 : i1 + 1]))
        fit["t_start_s"] = float(t[i0])
        report["plateaus"].append(fit)

    taus = sorted(
        tau
        for plateau in report["plateaus"]
        if plateau.get("verdict") == "relaxation identifiée"
        for tau in plateau["taus_s"]
    )
    fractions = [
        plateau["relaxation_fraction"]
        for plateau in report["plateaus"]
        if plateau.get("verdict") == "relaxation identifiée"
    ]
    report["identified_taus_s"] = taus
    report["identified_relaxation_fractions"] = fractions
    report["table_a1_taus_s"] = [154.57 / 20.67, 977.79 / 5.98, 11044.83 / 4.75]
    report["table_a1_long_term_fraction"] = 1.0 - 6.36 / 37.76

    F_act_exp = F - report["preload_mN"]
    curvature_exp = first_ramp_curvature(t, P, F_act_exp)
    report["curvature_experimental"] = curvature_exp

    if run_simulation:
        t_sim, F_sim = simulate_essai(t, P, overrides, mode=mode)
        F_sim_act = F_sim - F_sim[0]
        F_sim_on_exp = np.interp(t, t_sim, F_sim_act)
        scale = closed_form_scale(F_act_exp, F_sim_on_exp)
        scaled = scale * F_sim_on_exp
        rmse = float(np.sqrt(np.mean((scaled - F_act_exp) ** 2)))
        correlation = (
            float(np.corrcoef(F_act_exp, F_sim_on_exp)[0, 1]) if np.std(F_sim_on_exp) > 0 else float("nan")
        )
        report["simulation"] = {
            "scale_factor": scale,
            "rmse_scaled_mN": rmse,
            "rmse_relative_to_peak": float(rmse / max(np.max(np.abs(F_act_exp)), 1e-9)),
            "correlation": correlation,
        }
        report["curvature_simulated"] = first_ramp_curvature(t, P, scaled)
        if curvature_exp and report["curvature_simulated"]:
            ratio_exp = curvature_exp["curvature_ratio_b_over_a"]
            ratio_sim = report["curvature_simulated"]["curvature_ratio_b_over_a"]
            structural = abs(ratio_sim) < 0.2 * abs(ratio_exp)
            report["verdict_courbure"] = {
                "b_over_a_experimental_per_MPa": ratio_exp,
                "b_over_a_simulated_per_MPa": ratio_sim,
                "structurel": bool(structural),
                "commentaire": (
                    "La courbure F(P) mesurée est hors de portée du modèle quel que soit "
                    "le facteur d'échelle : limite STRUCTURELLE (dead-band/raidissement en "
                    "pression, hors du cadre de l'article) — à documenter, pas à compenser "
                    "par le spectre."
                    if structural
                    else "La courbure mesurée est comparable à celle du modèle : un ajustement "
                    "de spectre/géométrie peut suffire."
                ),
            }
    return report


def print_report(report: dict) -> None:
    print(f"\n=== Identification matériau — {Path(report['essai']).name} ===")
    for note in report["notes"]:
        print(f"  note : {note}")
    print(
        f"  {report['n_samples']} points, {report['duration_s']:.0f} s, "
        f"P max {report['pressure_max_MPa']:.3f} MPa, précharge {report['preload_mN']:.0f} mN"
    )
    print(f"  paliers détectés : {len(report['plateaus'])}")
    for plateau in report["plateaus"]:
        base = (
            f"    - P={plateau['pressure_MPa']:.3f} MPa, {plateau['duration_s']:.0f} s, "
            f"variation {100 * plateau['relative_change']:+.1f} % : {plateau['verdict']}"
        )
        if plateau.get("verdict") == "relaxation identifiée":
            taus = ", ".join(f"{tau:.1f} s" for tau in plateau["taus_s"])
            base += f" (tau = {taus} ; fraction {100 * plateau['relaxation_fraction']:.1f} %)"
        print(base)
    if report["identified_taus_s"]:
        print(f"  taus identifiés : {[round(x, 1) for x in report['identified_taus_s']]} s "
              f"(Table A1 : {[round(x, 1) for x in report['table_a1_taus_s']]} s)")
        print(f"  fractions de relaxation identifiées : "
              f"{[round(100 * x, 1) for x in report['identified_relaxation_fractions']]} % "
              f"(Table A1 à long terme : {100 * report['table_a1_long_term_fraction']:.0f} %)")
    else:
        print("  aucun palier ne montre de relaxation de Maxwell exploitable "
              "(voir plan d'essais E2/E3 : relaxation uniaxiale propre, maintien sans fuite).")
    if "simulation" in report:
        sim = report["simulation"]
        print(f"  simulation ({report['mode']}) : échelle {sim['scale_factor']:.2f}, "
              f"corrélation {sim['correlation']:.3f}, RMSE {sim['rmse_scaled_mN']:.0f} mN "
              f"({100 * sim['rmse_relative_to_peak']:.1f} % du pic)")
    verdict = report.get("verdict_courbure")
    if verdict:
        print(f"  courbure F(P) : b/a mesuré {verdict['b_over_a_experimental_per_MPa']:+.1f} /MPa "
              f"vs simulé {verdict['b_over_a_simulated_per_MPa']:+.1f} /MPa -> "
              f"{'STRUCTURELLE' if verdict['structurel'] else 'compatible spectre'}")
        print(f"    {verdict['commentaire']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Identification matériau (audit 2026-08, item 3.3).")
    parser.add_argument("--essai", required=True, help="essai bloqué mesuré (.csv ou .xlsx)")
    parser.add_argument("--fiche", default=None, help="fiche géométrique du muscle (.json)")
    parser.add_argument(
        "--mode",
        choices=("elastic_tk_reference", "viscoelastic_history"),
        default="elastic_tk_reference",
    )
    parser.add_argument("--sans-simulation", action="store_true", help="volet A (paliers) uniquement")
    args = parser.parse_args()

    essai_path = Path(args.essai)
    report = identify(
        essai_path,
        Path(args.fiche) if args.fiche else None,
        args.mode,
        run_simulation=not args.sans_simulation,
    )
    print_report(report)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = OUT_DIR / f"identification_{essai_path.stem}.json"
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)
    print(f"  rapport : {out_json}")


if __name__ == "__main__":
    main()
