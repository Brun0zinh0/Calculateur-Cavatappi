# -*- coding: utf-8 -*-
"""Audit numerique independant : transmission de pression et cinematique du tube.

1. Verification Lame : tube epais quasi-isotrope, la solution radiale du code
   doit reproduire sigma_r / sigma_phi de Lame et les CL -dP / 0.
2. Effet du profil d'angle de biais (paper_linear vs uniform_twist=eq.11).
3. Effet de section fixe vs mise a jour.
4. Verification des historiques de pression (periodes, pentes, valeurs).
"""
import sys
import numpy as np

sys.path.insert(
    0,
    r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2",
)
import Base

print("MODEL_VERSION =", Base.MODEL_VERSION)

# ---------------------------------------------------------------- 1. Lame ---
print("\n=== 1. Verification Lame (quasi-isotrope, theta_f = 0) ===")
E = 10.0
nu = 0.3
mat = Base.default_material_params(
    E_axial=E * (1.0 + 1e-9),
    E_radius=E,
    G12=E / (2.0 * (1.0 + nu)),
    nu12=nu,
    nu23=nu,
)
geom = Base.default_geometry_params(theta_f_deg=0.0)
disc = Base.default_discretization(n_layers=8, n_phi=8)
model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration="exponential")

dP = 1.0  # MPa
dv = 0.0
dw = 0.0
C_alg, hist = model._algorithmic_data(0.0)
coeffs = model._solve_radial_constants(dw, dv, dP, 0.0, C_alg, hist)

Rin, Rout = geom.Rin, geom.Rout
A_lame = dP * Rin**2 / (Rout**2 - Rin**2)
B_lame = dP * Rin**2 * Rout**2 / (Rout**2 - Rin**2)

rows = []
for j, R in enumerate(model.R_centers):
    u, du = model._u_du_layer(j, R, coeffs, dv, dw, C_alg)
    eps = np.array([dw, u / R, du, 0.0, 0.0, dv * R])
    sig = C_alg[j] @ eps
    sr_lame = A_lame - B_lame / R**2
    st_lame = A_lame + B_lame / R**2
    rows.append((R, sig[2], sr_lame, sig[1], st_lame))
print("   R      sig_r(code)  sig_r(Lame)  sig_phi(code)  sig_phi(Lame)")
for R, sr_c, sr_l, st_c, st_l in rows:
    print(f"  {R:5.3f}  {sr_c:+10.6f}  {sr_l:+10.6f}   {st_c:+10.6f}   {st_l:+10.6f}")

# CL exactes aux faces
j = 0
Aj, Bj, muj = model._layer_AB_mu(j, C_alg)
c_sig, k_sig = model._radial_stress_axisym_coeff(j, model.R_edges[0], Aj, Bj, muj, dv, dw, 0.0, C_alg)
sr_in = float(c_sig @ coeffs[0]) + k_sig
j = disc.n_layers - 1
Aj, Bj, muj = model._layer_AB_mu(j, C_alg)
c_sig, k_sig = model._radial_stress_axisym_coeff(j, model.R_edges[-1], Aj, Bj, muj, dv, dw, 0.0, C_alg)
sr_out = float(c_sig @ coeffs[-1]) + k_sig
print(f"  CL interne : sigma_r(Rin) = {sr_in:+.8f}  (attendu {-dP:+.1f})")
print(f"  CL externe : sigma_r(Rout) = {sr_out:+.8e}  (attendu 0)")
err = max(abs(sr_c - sr_l) for _, sr_c, sr_l, _, _ in rows)
errt = max(abs(st_c - st_l) for _, _, _, st_c, st_l in rows)
print(f"  Erreur max sigma_r vs Lame : {err:.3e} MPa ; sigma_phi : {errt:.3e} MPa")

# --------------------------------------- 2. Profil d'angle de biais --------
print("\n=== 2. Effet du profil d'angle de biais (defaut vs eq. 11 arctan) ===")
theta_f = 37.91
frac = np.linspace(0, 1, 11)
lin = frac * np.deg2rad(theta_f)
arct = np.arctan(frac * np.tan(np.deg2rad(theta_f)))
print("  r/Rout   theta_lin(deg)  theta_arctan(deg)  ecart(deg)")
for f, a, b in zip(frac, lin, arct):
    print(f"  {f:5.2f}    {np.rad2deg(a):8.3f}       {np.rad2deg(b):8.3f}      {np.rad2deg(b-a):+6.3f}")


def run_case(bias_profile, section_mode, label):
    geom = Base.default_geometry_params(
        bias_angle_profile=bias_profile, section_update_mode=section_mode
    )
    cfg = Base.default_simulation_config(
        geom=geom, n_cycles=1, Pmax=1.3, dt=0.25, n_layers=4, n_phi=24, pre_steps=24,
        nonlinear_pressure=False,
    )
    model, arr = Base.run_blocked_actuation(cfg, eps=0.8)
    fpk = float(np.max(arr["force_mN"]))
    fmin = float(np.min(arr["force_mN"]))
    tpk = float(np.max(arr["torque_microNm"]))
    f_act_pk = float(np.max(arr["force_act_mN"]))
    t_act_pk = float(np.max(arr["torque_act_microNm"]))
    print(
        f"  {label:44s} Fmax={fpk:9.2f} mN  Fmin={fmin:9.2f}  Tmax={tpk:8.2f} uNm  "
        f"dF_act={f_act_pk:8.2f}  dT_act={t_act_pk:8.2f}"
    )
    return fpk, tpk, f_act_pk, t_act_pk


base = run_case("paper_linear", "fixed", "defaut  : paper_linear + section fixe")
alt = run_case("uniform_twist", "fixed", "article : arctan (eq.11) + section fixe")
upd = run_case("paper_linear", "updated", "paper_linear + section mise a jour")
upd2 = run_case("uniform_twist", "updated", "arctan (eq.11) + section mise a jour")

print("\n  Ecarts relatifs vs defaut :")
for name, case in (("arctan/fixe", alt), ("linear/updated", upd), ("arctan/updated", upd2)):
    dF = 100.0 * (case[0] - base[0]) / abs(base[0])
    dT = 100.0 * (case[1] - base[1]) / abs(base[1])
    dFa = 100.0 * (case[2] - base[2]) / abs(base[2])
    dTa = 100.0 * (case[3] - base[3]) / abs(base[3])
    print(f"   {name:16s}: Fmax {dF:+6.2f}%  Tmax {dT:+6.2f}%  dF_act {dFa:+6.2f}%  dT_act {dTa:+6.2f}%")

# --------------------------------------------- 3. Historiques de pression --
print("\n=== 3. Historiques de pression ===")
t, P = Base.cyclic_pressure_history(n_cycles=3, Pmax=1.3, flow_rate_mL_min=10.0,
                                    volume_mL=1.50, dt=0.15, nonlinear=False)
half = 60.0 * 1.50 / 10.0
print(f"  demi-periode attendue 60*V/Q = {half} s ; periode = {2*half} s ; total = {t[-1]} s (attendu {6*half})")
ipk = np.argmax(P)
print(f"  premier pic : P={P[ipk]:.3f} MPa a t={t[ipk]} s (attendu 1.3 a 9 s)")
peaks = t[np.isclose(P, 1.3)]
print(f"  instants des pics P=Pmax : {peaks}")
print(f"  P(t=0)={P[0]}, P(fin)={P[-1]}")
# pente lineaire
mask = (t > 0.1) & (t < 8.9)
slope = np.polyfit(t[mask], P[mask], 1)[0]
print(f"  pente de charge lineaire : {slope:.5f} MPa/s (attendu {1.3/9:.5f})")

tn, Pn = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.3, dt=0.05, nonlinear=True)
i25 = np.argmin(np.abs(tn - 0.25 * half))
i50 = np.argmin(np.abs(tn - 0.5 * half))
print(f"  profil non lineaire : P(25% charge)={Pn[i25]:.4f} MPa (=Pmax*0.25^3.5={1.3*0.25**3.5:.4f}) ;"
      f" P(50%)={Pn[i50]:.4f} (=Pmax*0.5^3.5={1.3*0.5**3.5:.4f})")

t2, P2 = Base.ramp_hold_pressure_history(P_hold=1.3, ramp_time=9.0, hold_time=30.0, dt=0.5,
                                         nonlinear_ramp=False)
print(f"  ramp_hold lineaire : P(9s)={P2[np.isclose(t2,9.0)]} ; P(fin)={P2[-1]} ; total={t2[-1]} s")

# --------------------------------------------- 4. pression.py --------------
print("\n=== 4. pression.py : conversions ===")
sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import pression

# valeurs de reference exactes
psi_exact = 4.4482216152605 / (0.0254**2) / 1e6  # MPa
print(f"  1 psi = {psi_exact:.15f} MPa (reference NIST)")
cols = {"t (s)": np.array([0.0, 1.0, 2.0, 3.0]),
        "P (psi)": np.array([0.0, 100.0, 50.0, 145.0377])}
out = pression.measured_pressure_payload(cols, "t (s)", "P (psi)", "psi", subtract_initial=False)
print(f"  100 psi -> {out['pressure_MPa'][1]:.6f} MPa (attendu {100*psi_exact:.6f})")
print(f"  145.0377 psi -> {out['pressure_MPa'][3]:.6f} MPa (attendu ~1.0)")
cols2 = {"t": np.array([0.0, 1.0, 2.0]), "P (bar)": np.array([2.0, 5.0, 13.0])}
out2 = pression.measured_pressure_payload(cols2, "t", "P (bar)", "bar", subtract_initial=True)
print(f"  bar avec zero initial 2.0 bar : {out2['pressure_MPa']} (attendu [0, 0.3, 1.1])")
try:
    cols3 = {"t": np.array([0.0, 1.0]), "P": np.array([0.0, 16.0])}
    pression.measured_pressure_payload(cols3, "t", "P", "bar", subtract_initial=False)
    print("  ERREUR : plafond 1.5 MPa non applique !")
except ValueError as e:
    print(f"  plafond 1.5 MPa applique : {e}")
print(f"  infer_force_unit('column1') = {pression.infer_force_unit('column1')} (piege 'mn')")
print(f"  infer_pressure_unit('P_kPa') = {pression.infer_pressure_unit('P_kPa')}")
