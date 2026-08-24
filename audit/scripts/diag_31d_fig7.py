# -*- coding: utf-8 -*-
"""Phase 3.1d — figure 7 avec le mode updated corrigé (point milieu complet).

Protocole figure 7 complet (préfixe 0-180 s + pression numérisée recalée),
combinaisons {bias} x {section}, L3p16 dt=1. Compare aux cibles théorie
article et ancrages texte (comme la matrice d'audit).
"""
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parent / "alpha V3"
sys.path.insert(0, str(APP))
VAL = APP / "validation"
sys.path.insert(0, str(VAL))

import Base  # noqa: E402
from parametres import DEFAULT_SETTINGS, build_config  # noqa: E402
import validation_figure7 as v  # noqa: E402

THEORY = {0.8: dict(F=(1758.0, 1806.0), T=(591.0, 733.0)), 1.0: dict(F=(2389.0, 2466.0), T=(702.0, 876.0))}
TEXT = {0.8: dict(F=1824.60, T=877.01), 1.0: dict(F=2424.67, T=873.41)}

print("digitisation + recalage...", flush=True)
targets, pinfo = v.build_targets("rescale")

def run(bias, section, eps, dt=1.0, n_layers=3, n_phi=16):
    s = dict(DEFAULT_SETTINGS)
    n_cycles = int(np.ceil(v.PAPER_T_MAX / (2.0 * 60.0 * v.PAPER_VOLUME_ML / v.PAPER_FLOW_RATE_ML_MIN))) + 1
    s.update(
        {
            "eps": eps, "n_cycles": n_cycles, "use_fixed_duration": False,
            "p_max_mpa": v.PAPER_PMAX_MPA, "dt": dt, "n_layers": n_layers, "n_phi": n_phi,
            "pre_steps": 24, "flow_rate_mL_min": v.PAPER_FLOW_RATE_ML_MIN,
            "volume_mL": v.PAPER_VOLUME_ML, "nonlinear_pressure": True,
            "bias_angle_profile": bias, "section_update_mode": section,
        }
    )
    cfg = build_config(s)
    pt, pm = v.digitized_pressure_history(eps, targets, prefix_dt=dt, figure_dt=0.5 * dt)
    _, data = Base.run_blocked_actuation(cfg, pressure_time=pt, pressure_MPa=pm)
    win = v.window_model_data(data)
    exF = v.cycle_extrema(win["time"], win["force_total_mN"])
    exT = v.cycle_extrema(win["time"], win["torque_act_microNm"])
    return exF, exT

print(f"{'config':>24} {'eps':>4} | {'F pics':>8} {'vs th.':>14} {'vs texte':>9} | {'T pics':>8} {'vs th.':>16} {'vs texte':>9} | {'T vallee':>9}")
for bias in ("paper_linear", "uniform_twist"):
    for section in ("fixed", "updated"):
        for eps in (0.8, 1.0):
            try:
                exF, exT = run(bias, section, eps)
                thF = THEORY[eps]["F"]; thT = THEORY[eps]["T"]
                fp = exF["peak_mean"]; tp = exT["peak_mean"]
                gF = 100.0 * (fp / (0.5 * (thF[0] + thF[1])) - 1.0)
                gT = 100.0 * (tp / (0.5 * (thT[0] + thT[1])) - 1.0)
                gFt = 100.0 * (fp / TEXT[eps]["F"] - 1.0)
                gTt = 100.0 * (tp / TEXT[eps]["T"] - 1.0)
                print(f"{bias + '/' + section:>24} {eps:>4.1f} | {fp:>8.1f} {gF:>+13.1f}% {gFt:>+8.1f}% | "
                      f"{tp:>8.1f} {gT:>+15.1f}% {gTt:>+8.1f}% | {exT['valley_mean']:>9.1f}", flush=True)
            except Exception as exc:
                print(f"{bias + '/' + section:>24} {eps:>4.1f} | ECHEC : {type(exc).__name__}: {exc}", flush=True)

print("\nrappel matrice d'audit (updated NON corrige, pression corr) :")
print("  lin/updated  e0.8 : F -12.6 % vs th., T -64.5 % ; e1.0 : F -20.3 %, T -84.6 %")
print("  arctan/updated e0.8 : F -11.4 %, T -59.2 % ; e1.0 : F -19.4 %, T -78.0 %")
print("  lin/fixed (reference) e0.8 : F -5.1 %, T -12.1 % ; e1.0 : F -12.3 %, T -26.8 %")
