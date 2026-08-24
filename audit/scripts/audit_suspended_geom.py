# -*- coding: utf-8 -*-
"""Audit mode suspendu : equilibre Eq.(23) EXP et conventions de contraction."""
import sys
import numpy as np

BASE_DIR = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
sys.path.insert(0, BASE_DIR)
import Base  # noqa: E402

model, arr = Base.run_suspended_actuation(
    load_N=1.0,
    eps=0.0,
    n_cycles=1,
    Pmax=1.2,
    dt=2.0,
    n_layers=4,
    n_phi=12,
    pre_steps=8,
    flow_rate_mL_min=10.0,
    volume_mL=1.5,
    nonlinear_pressure=False,
    integration="exponential",
)

alpha = arr["alpha_rad"]; rho = arr["rho_mm"]
sin_a, cos_a = np.sin(alpha), np.cos(alpha)
F = 1.0
r1 = np.max(np.abs(arr["Ftube_axis_N"] + arr["Fnylon_axis_N"] - F * sin_a))
r2 = np.max(np.abs(arr["Mtube_Nmm"] + arr["Mnylon_Nmm"] + F * rho * sin_a))
r3 = np.max(np.abs(arr["Ttube_axis_Nmm"] + arr["Tnylon_axis_Nmm"] - F * rho * cos_a))
print(f"max |F_tube+F_ny - F sin b|      = {r1:.3e} N")
print(f"max |M_tube+M_ny + F rho sin b|  = {r2:.3e} N.mm")
print(f"max |T_tube+T_ny - F rho cos b|  = {r3:.3e} N.mm")

L = arr["axial_length_mm"]
c = arr["free_contraction_mm"]
s = arr["free_actuation_strain"]
print(f"L[0]={L[0]:.4f} mm ; reference={arr['reference_axial_length_mm'][0]:.4f} mm")
iP = np.argmax(arr["pressure_MPa"])
print(f"a P max ({arr['pressure_MPa'][iP]:.2f} MPa) : L={L[iP]:.4f} mm ; contraction={c[iP]:+.4f} mm ; strain={s[iP]:+.4f}")
print(f"convention : contraction>0 quand L diminue -> {'OK' if (c[iP] > 0) == (L[iP] < L[0]) else 'INCOHERENT'}")
# verification L = L0 * stretch * sin(alpha)/sin(alpha0) (Eq. 24 EXP)
L_pred = model.geom.initial_length * arr["axial_stretch"] * np.sin(alpha) / np.sin(model.alpha0)
print(f"max |L_active - L0*stretch*sin a/sin a0| = {np.max(np.abs(arr['active_axial_length_mm'] - L_pred)):.3e} mm")
# en mode suspendu le nombre de tours n'est pas conserve : h*2piN != L_active en general
Lg = arr["axial_length_geometry_mm"]
print(f"max |L_geom(2piN h) - L_active| = {np.max(np.abs(Lg - arr['active_axial_length_mm'])):.3e} mm (attendu ~0 : meme formule ?)")
