# -*- coding: utf-8 -*-
"""Audit numerique independant du mode masse suspendue / relaxation de Base.py."""
import sys
import numpy as np

sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import Base  # noqa: E402

np.set_printoptions(precision=5, suppress=True)

LOAD = 1.0  # N, comme l'article EXP

# ------------------------------------------------------------------
# 1) Run suspendu, eps=0, rampe lineaire 0 -> 1.2 MPa en 60 s puis maintien
# ------------------------------------------------------------------
t = np.arange(0.0, 121.0, 2.0)
p = np.clip(0.02 * t, 0.0, 1.2)

model, arr = Base.run_suspended_actuation(
    load_N=LOAD,
    pressure_time=t,
    pressure_MPa=p,
    eps=0.0,
    n_layers=3,
    n_phi=8,
    pre_steps=4,
    dt=2.0,
    integration="exponential",
)

alpha = arr["alpha_rad"]
rho = arr["rho_mm"]
Ftube = arr["Ftube_axis_N"]
Fny = arr["Fnylon_axis_N"]
Mtube = arr["Mtube_Nmm"]
Mny = arr["Mnylon_Nmm"]
Ttube = arr["Ttube_axis_Nmm"]
Tny = arr["Tnylon_axis_Nmm"]

r1 = Ftube + Fny - LOAD * np.sin(alpha)
r2 = Mtube + Mny + LOAD * rho * np.sin(alpha)
r3 = Ttube + Tny - LOAD * rho * np.cos(alpha)
print("=== Verification eq. (23) EXP sur l'historique suspendu ===")
print("max |F_tube+F_ny - F sin a|  =", np.nanmax(np.abs(r1)), "N")
print("max |M_tube+M_ny + F rho sin a| =", np.nanmax(np.abs(r2)), "N.mm")
print("max |T_tube+T_ny - F rho cos a| =", np.nanmax(np.abs(r3)), "N.mm")

print("\n=== Conventions de contraction ===")
print("reference_axial_length_mm =", arr["reference_axial_length_mm"][0])
print("initial_length (geom)     =", model.geom.initial_length)
print("axial_length t=0          =", arr["axial_length_mm"][0])
print("axial_length final        =", arr["axial_length_mm"][-1])
print("free_contraction_mm final =", arr["free_contraction_mm"][-1])
print("free_actuation_% final    =", arr["free_actuation_percent"][-1])
print("free_actuation_% (norme L_T0 article) =",
      100.0 * arr["free_contraction_mm"][-1] / model.geom.initial_length)
print("signe: contraction >0 attendu sous pression :", arr["free_contraction_mm"][-1] > 0)

print("\n=== Conservation du nombre de spires (doit differer en mode libre) ===")
print("axial_length_mm (fibre, eq.24)      final =", arr["active_axial_length_mm"][-1])
print("axial_length_geometry_mm (N constant) final =", arr["axial_length_geometry_mm"][-1])
dN = (arr["active_centerline_length_mm"] * np.cos(alpha) / (2 * np.pi * rho)) - model.turns
print("turns implicites - turns initiaux : min", dN.min(), "max", dN.max())

print("\n=== Etat apres equilibrage (t=0, P=0) ===")
print("force_mN[0] =", arr["force_mN"][0], "(doit etre ~ load*1000 =", 1000*LOAD, ")")
print("residual max =", np.nanmax(arr["residual"]))
print("suspended_equivalent_settling_time_s =", arr["suspended_equivalent_settling_time_s"])

# ------------------------------------------------------------------
# 2) Transitoire de mise en charge exclu par defaut ?
#    comparaison equilibrate=True vs False a P=0 constant
# ------------------------------------------------------------------
print("\n=== P=0 constant pendant 300 s : le fluage sous poids est-il elimine par defaut ? ===")
t0 = np.arange(0.0, 301.0, 10.0)
p0 = np.zeros_like(t0)
for eq in (True, False):
    m2, a2 = Base.run_suspended_actuation(
        load_N=LOAD, pressure_time=t0, pressure_MPa=p0,
        eps=0.0, n_layers=3, n_phi=8, pre_steps=4, dt=10.0,
        integration="exponential",
        equilibrate_load_before_pressure=eq,
    )
    L = a2["axial_length_mm"]
    print(f"equilibrate={eq}:  L(0)={L[0]:.4f} mm, L(300 s)={L[-1]:.4f} mm, "
          f"derive={L[-1]-L[0]:+.5f} mm, actionnement final={a2['free_actuation_percent'][-1]:+.4f} %")

# ------------------------------------------------------------------
# 3) Effet de la precontrainte par defaut (eps=0.8) en mode suspendu
# ------------------------------------------------------------------
print("\n=== eps=0.8 (defaut interface) vs eps=0 en mode suspendu equilibre ===")
for eps in (0.0, 0.8):
    m3, a3 = Base.run_suspended_actuation(
        load_N=LOAD, pressure_time=t, pressure_MPa=p,
        eps=eps, n_layers=3, n_phi=8, pre_steps=6, dt=2.0,
        integration="exponential",
    )
    print(f"eps={eps}: L(0)={a3['axial_length_mm'][0]:.4f} mm, "
          f"actionnement final={a3['free_actuation_percent'][-1]:+.4f} %")

# ------------------------------------------------------------------
# 4) run_hold_relaxation : nature de la simulation
# ------------------------------------------------------------------
print("\n=== run_hold_relaxation (bloque, rampe+maintien) ===")
m4, a4 = Base.run_hold_relaxation(
    eps=0.5, P_hold=1.0, hold_time=200.0, ramp_time=9.0, dt=2.0,
    n_layers=3, n_phi=8, pre_steps=4,
)
i0 = a4["hold_start_index"]
print("pression pendant maintien: min", a4["pressure_MPa"][i0:].min(),
      "max", a4["pressure_MPa"][i0:].max(), "MPa (constante ?)")
F = a4["force_total_mN"]
print(f"force au debut du maintien = {F[i0]:.2f} mN ; a la fin = {F[-1]:.2f} mN ; "
      f"variation = {F[-1]-F[i0]:+.2f} mN")
print("h (pas) constant ?", np.allclose(a4["h_mm_per_rad"], a4["h_mm_per_rad"][0], rtol=1e-9))
print("-> relaxation de force a pression constante en mode bloque")
