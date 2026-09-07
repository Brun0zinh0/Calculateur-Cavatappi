# -*- coding: utf-8 -*-
"""Phase 3.1 — instrumentation de la dérive du mode section_update_mode='updated'.

Trois runs bloqués, pression cyclique linéaire générée (Pmax 1.4, 11 cycles) :
  A) updated + Maxwell complet      (la dérive observée par la matrice)
  B) updated + quasi élastique      (eta x 1e9 : si la dérive disparaît -> visqueux)
  C) fixed  + Maxwell complet       (référence stable)

Sorties : diag_31_series.csv (par pas), diag_31_cycles.csv (fins de cycle),
diag_31_drift.png, et un résumé imprimé.
"""
import sys
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parent / "alpha V3"
sys.path.insert(0, str(APP))

import Base  # noqa: E402
import parametres  # noqa: E402

OUT = Path(__file__).resolve().parent
EPS = 0.8
PMAX = 1.4
N_CYCLES = 11
DT = 1.0
LAYERS = 3
PHI = 16


def build(settings_extra):
    s = dict(parametres.DEFAULT_SETTINGS)
    s.update(
        {
            "eps": EPS,
            "n_cycles": N_CYCLES,
            "use_fixed_duration": False,
            "p_max_mpa": PMAX,
            "dt": DT,
            "n_layers": LAYERS,
            "n_phi": PHI,
            "pre_steps": 24,
            "nonlinear_pressure": False,
            # v4-16 : cles historiques conservees -> demi-periode exacte de 9 s
            "flow_rate_mL_min": 10.0,
            "volume_mL": 1.5,
        }
    )
    s.update(settings_extra)
    return parametres.build_config(s)


def run_instrumented(tag, settings_extra):
    cfg = build(settings_extra)
    disc = Base.default_discretization(
        n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05)
    )
    model = Base.TCPAMaxwellBlockedModel(
        mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration,
        prestrain_reference_mode=cfg.prestrain_reference_mode,
    )
    model.prestretch_to(cfg.eps)
    t, p = Base.cyclic_pressure_history(
        n_cycles=cfg.n_cycles, Pmax=cfg.Pmax, dt=cfg.dt, nonlinear=False,
        half_period_s=Base._config_half_period(cfg),
    )
    model.step(float(p[0]), 0.0, h_target=model.h_blocked)
    model.lock_blocked_series_reference()
    rows = []

    def record(time, pressure):
        rows.append(
            dict(
                t=time, P=pressure,
                Rin=float(model.R_edges[0]), Rout=float(model.R_edges[-1]),
                area=float(np.pi * (model.R_edges[-1] ** 2 - model.R_edges[0] ** 2)),
                theta_deg_mean=float(np.rad2deg(np.mean(model.theta_layers))),
                theta_deg_max=float(np.rad2deg(np.max(model.theta_layers))),
                rho=float(model.helix.rho), alpha_deg=float(np.rad2deg(model.helix.alpha)),
                F_mN=1000.0 * float(model.history[-1].Ft),
                T_uNm=1000.0 * float(model.history[-1].Tt),
            )
        )

    record(float(t[0]), float(p[0]))
    for k in range(1, len(t)):
        model.step(float(p[k]), float(t[k] - t[k - 1]), h_target=model.h_blocked)
        record(float(t[k]), float(p[k]))
    return rows


VARIANTS = {
    "A_updated_maxwell": {"section_update_mode": "updated"},
    "B_updated_quasielastic": {
        "section_update_mode": "updated",
        "maxwell_eta1_mpa_s": 154.57e9,
        "maxwell_eta2_mpa_s": 977.79e9,
        "maxwell_eta3_mpa_s": 11044.83e9,
    },
    "C_fixed_maxwell": {"section_update_mode": "fixed"},
}

period = 2.0 * 60.0 * 1.5 / 10.0
all_rows = {}
for tag, extra in VARIANTS.items():
    print(f"--- run {tag} ...", flush=True)
    try:
        all_rows[tag] = run_instrumented(tag, extra)
        print(f"    ok, {len(all_rows[tag])} pas", flush=True)
    except Exception as exc:
        print(f"    ECHEC : {type(exc).__name__}: {exc}", flush=True)
        all_rows[tag] = None

import csv

with open(OUT / "diag_31_series.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["variant", "t", "P", "Rin", "Rout", "area", "theta_deg_mean", "theta_deg_max", "rho", "alpha_deg", "F_mN", "T_uNm"])
    for tag, rows in all_rows.items():
        if rows:
            for r in rows:
                w.writerow([tag] + [r[k] for k in ("t", "P", "Rin", "Rout", "area", "theta_deg_mean", "theta_deg_max", "rho", "alpha_deg", "F_mN", "T_uNm")])

print("\n=== fins de cycle (P retour a 0) : derive par variante ===")
summary = {}
for tag, rows in all_rows.items():
    if not rows:
        continue
    print(f"\n[{tag}]")
    print("cycle |   Rin    |   Rout   |  aire   | theta_max | F(mN)  | T(uNm)")
    cyc = []
    for c in range(1, N_CYCLES + 1):
        tc = c * period
        best = min(rows, key=lambda r: abs(r["t"] - tc))
        cyc.append(best)
        print(f"{c:5d} | {best['Rin']:.5f} | {best['Rout']:.5f} | {best['area']:.5f} | {best['theta_deg_max']:9.4f} | {best['F_mN']:7.1f} | {best['T_uNm']:7.1f}")
    summary[tag] = cyc
    r1, rN = cyc[0], cyc[-1]
    print(f"  derive c1->c{N_CYCLES} : Rin {100*(rN['Rin']/r1['Rin']-1):+.2f} %, Rout {100*(rN['Rout']/r1['Rout']-1):+.2f} %, "
          f"aire {100*(rN['area']/r1['area']-1):+.2f} %, theta_max {rN['theta_deg_max']-r1['theta_deg_max']:+.3f} deg, "
          f"F {rN['F_mN']-r1['F_mN']:+.1f} mN, T {rN['T_uNm']-r1['T_uNm']:+.1f} uNm")

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for tag, rows in all_rows.items():
        if not rows:
            continue
        tt = [r["t"] for r in rows]
        axes[0, 0].plot(tt, [r["Rin"] for r in rows], label=tag, lw=1.0)
        axes[0, 1].plot(tt, [r["theta_deg_max"] for r in rows], label=tag, lw=1.0)
        axes[1, 0].plot(tt, [r["F_mN"] for r in rows], label=tag, lw=1.0)
        axes[1, 1].plot(tt, [r["T_uNm"] for r in rows], label=tag, lw=1.0)
    for ax, ttl in zip(axes.flat, ["Rin (mm)", "theta max (deg)", "Force bloquee (mN)", "Couple (uNm)"]):
        ax.set_title(ttl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
        ax.set_xlabel("t (s)")
    fig.savefig(OUT / "diag_31_drift.png", dpi=140)
    print(f"\nfigure : {OUT / 'diag_31_drift.png'}")
except Exception as exc:
    print("plot impossible:", exc)
