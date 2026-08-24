# -*- coding: utf-8 -*-
"""Phase 3.2 — validation du mode de précontrainte viscoélastique.

A) Training 11 cycles (protocole article, Pmax 1.4, non linéaire) :
   atténuation des pics c1->c11, élastique vs viscoélastique.
   Cibles : article -10.3 % ; prototype contre-expertise -6.4 % ; élastique -1.8 %.
B) Démonstration fig. 11 d'EXP : mode suspendu, précontrainte viscoélastique,
   60 s de relaxation à P = 0 après étirement (equilibrate=False).
C) Arbitrage figure 7 : {précontrainte} x {section} aux deux eps, pression
   numérisée recalée — niveaux absolus des pics ET des vallées.

Résultats obtenus le 21/08/2026 (moteur 2026.08.21-audit-phase32-13) :
  A) élastique -0.9 % ; viscoélastique -5.0 % (article -10.3 %).
  B) élastique -1.775 mm (signe opposé à l'expérience) ; visco +0.193 mm (bon signe).
  C) élast./updated : F -3.0/-9.0 % vs texte, T -7.0/-5.9 % ; vallée F 1448 mm
     à eps=0.8 = vallée théorie article (~1450) -> la simulation publiée n'a
     vraisemblablement pas relaxé la prétension.
     visco/updated : T -2.1/-0.7 % (quasi parfait) mais F -24.6/-29.9 %
     (le spectre Table A1 sur-relaxe, constat M6).
"""
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(APP))
sys.path.insert(0, str(APP / "validation"))

import Base  # noqa: E402
from parametres import DEFAULT_SETTINGS, build_config  # noqa: E402
import validation_figure7 as v  # noqa: E402

print("Base :", Base.MODEL_VERSION)

# ---------------------------------------------------------------- A) training
print("\n=== A) Training 11 cycles (eps=0.8, Pmax=1.4, profil non lineaire) ===")


def training(mode):
    s = dict(DEFAULT_SETTINGS)
    s.update(
        {
            "eps": 0.8, "n_cycles": 11, "use_fixed_duration": False, "p_max_mpa": 1.4,
            "dt": 1.0, "n_layers": 3, "n_phi": 16, "pre_steps": 24,
            "nonlinear_pressure": True, "prestrain_reference_mode": mode,
        }
    )
    cfg = build_config(s)
    _, data = Base.run_blocked_actuation(cfg)
    t = np.asarray(data["time"], float)
    f = np.asarray(data["force_total_mN"], float)
    period = 18.0
    peaks = []
    for c in range(11):
        m = (t >= c * period) & (t < (c + 1) * period)
        if np.any(m):
            peaks.append(float(np.max(f[m])))
    tk = float(f[0])
    return tk, peaks


for mode in ("elastic_tk_reference", "viscoelastic_history"):
    tk, pk = training(mode)
    att = 100.0 * (pk[-1] / pk[0] - 1.0)
    print(f"  [{mode}] F(t_k)={tk:7.1f} mN ; pics c1={pk[0]:7.1f} c11={pk[-1]:7.1f} mN ; "
          f"attenuation c1->c11 = {att:+.1f} %", flush=True)
print("  cibles : article -10.3 % (1041.7->934.4 mN a eps=0.5) ; prototype visco -6.4 % ; elastique -1.8 %")

# ---------------------------------------------------------------- B) fig. 11
print("\n=== B) Demonstration fig. 11 EXP : 60 s de relaxation a P=0 apres etirement (suspendu, 1 N) ===")
for mode in ("elastic_tk_reference", "viscoelastic_history"):
    s = dict(DEFAULT_SETTINGS)
    s.update({"eps": 0.8, "dt": 2.0, "n_layers": 2, "n_phi": 8, "pre_steps": 12,
              "prestrain_reference_mode": mode})
    cfg = build_config(s)
    t = np.arange(0.0, 62.0, 2.0)
    p = np.zeros_like(t)
    _, data = Base.run_suspended_actuation(
        cfg, load_N=1.0, pressure_time=t, pressure_MPa=p,
        equilibrate_load_before_pressure=False,
    )
    L = np.asarray(data["axial_length_mm"], float)
    print(f"  [{mode}] L(0+)={L[0]:.3f} mm -> L(60 s)={L[-1]:.3f} mm ; "
          f"fluage sous poids pendant la relaxation post-etirement : {L[-1]-L[0]:+.3f} mm", flush=True)

# ---------------------------------------------------------------- C) figure 7
print("\n=== C) Figure 7 (pression numerisee recalee, L3p16, dt=1) ===")
targets, _ = v.build_targets("rescale")
TEXT = {0.8: dict(F=1824.60, T=877.01), 1.0: dict(F=2424.67, T=873.41)}


def fig7(mode, section, eps):
    s = dict(DEFAULT_SETTINGS)
    n_cycles = int(np.ceil(v.PAPER_T_MAX / 18.0)) + 1
    s.update(
        {
            "eps": eps, "n_cycles": n_cycles, "use_fixed_duration": False,
            "p_max_mpa": v.PAPER_PMAX_MPA, "dt": 1.0, "n_layers": 3, "n_phi": 16,
            "pre_steps": 24, "flow_rate_mL_min": 10.0, "volume_mL": 1.5,
            "nonlinear_pressure": True, "prestrain_reference_mode": mode,
            "section_update_mode": section,
        }
    )
    cfg = build_config(s)
    pt, pm = v.digitized_pressure_history(eps, targets, prefix_dt=1.0, figure_dt=0.5)
    _, data = Base.run_blocked_actuation(cfg, pressure_time=pt, pressure_MPa=pm)
    win = v.window_model_data(data)
    return (
        v.cycle_extrema(win["time"], win["force_total_mN"]),
        v.cycle_extrema(win["time"], win["torque_act_microNm"]),
    )


print(f"{'precontrainte':>22} {'section':>8} {'eps':>4} | {'F pics':>8} {'vs texte':>9} {'F vallees':>10} | {'T pics':>8} {'vs texte':>9}")
for mode in ("elastic_tk_reference", "viscoelastic_history"):
    for section in ("fixed", "updated"):
        for eps in (0.8, 1.0):
            try:
                exF, exT = fig7(mode, section, eps)
                gFt = 100.0 * (exF["peak_mean"] / TEXT[eps]["F"] - 1.0)
                gTt = 100.0 * (exT["peak_mean"] / TEXT[eps]["T"] - 1.0)
                print(f"{mode:>22} {section:>8} {eps:>4.1f} | {exF['peak_mean']:>8.1f} {gFt:>+8.1f}% "
                      f"{exF['valley_mean']:>10.1f} | {exT['peak_mean']:>8.1f} {gTt:>+8.1f}%", flush=True)
            except Exception as exc:
                print(f"{mode:>22} {section:>8} {eps:>4.1f} | ECHEC : {type(exc).__name__}: {exc}", flush=True)
print("\nreperes : vallees theorie article ~1450 mN (eps=0.8) / ~2050-2100 mN (eps=1.0) ;")
print("pics theorie 1758-1806 / 2389-2466 mN ; ancrages exp 1824.6/2424.7 mN et 877/873 uNm")
