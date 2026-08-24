# -*- coding: utf-8 -*-
"""Audit indépendant : géométrie hélicoïdale et décomposition force/couple.

Vérifications :
 A. Courbure/torsion de l'hélice (formules cos^2a/rho et sin2a/(2rho)) vs géométrie différentielle numérique.
 B. _new_geometry_from_dw_and_h : conservation du nombre de tours (Eq. 6 BLOCKED) et pas impose (Eq. 9).
 C. Décomposition force/couple (Eq. 3 BLOCKED) vs statique vectorielle 3D indépendante.
 D. Exécution du modèle : cohérence des sorties (équilibre, longueurs, conservation).
 E. Profil d'angle de biais : "paper_linear" vs Eq. (11) de l'article (arctan).
 F. Poissons effectifs : nu(s->phi) vs nu(s->r) (échange apparent de l'Eq. 18).
"""
import sys
import numpy as np

BASE_DIR = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"
sys.path.insert(0, BASE_DIR)

import Base  # noqa: E402

np.set_printoptions(precision=6, suppress=True)
print("MODEL_VERSION =", Base.MODEL_VERSION)

# ---------------------------------------------------------------- A
print("\n=== A. Courbure / torsion d'une hélice (référence indépendante) ===")
rho, alpha = 2.16, np.deg2rad(10.53)
hh = rho * np.tan(alpha)


def helix(t):
    return np.array([rho * np.cos(t), rho * np.sin(t), hh * t])


def deriv(f, t, eps=1e-6):
    return (f(t + eps) - f(t - eps)) / (2 * eps)


t0 = 0.7
r1 = deriv(helix, t0)
r2 = deriv(lambda t: deriv(helix, t), t0)
r3 = deriv(lambda t: deriv(lambda u: deriv(helix, u), t), t0)
kappa_num = np.linalg.norm(np.cross(r1, r2)) / np.linalg.norm(r1) ** 3
tau_num = np.dot(np.cross(r1, r2), r3) / np.linalg.norm(np.cross(r1, r2)) ** 2
kappa_formula = np.cos(alpha) ** 2 / rho
tau_formula = np.sin(2 * alpha) / (2 * rho)
print(f"kappa numerique={kappa_num:.8f}  formule cos^2a/rho={kappa_formula:.8f}  ecart={abs(kappa_num-kappa_formula):.2e}")
print(f"tau   numerique={tau_num:.8f}  formule sin2a/(2rho)={tau_formula:.8f}  ecart={abs(tau_num-tau_formula):.2e}")

# ---------------------------------------------------------------- B
print("\n=== B. _new_geometry_from_dw_and_h : Eq.(6) et Eq.(9) BLOCKED ===")
model = Base.TCPAMaxwellBlockedModel()
rho_old, alpha_old = model.helix.rho, model.helix.alpha
for dw, h_t in [(0.01, model.h0), (0.03, 1.4 * model.h0), (-0.02, model.h0)]:
    rho_new, alpha_new = model._new_geometry_from_dw_and_h(dw, h_t)
    lhs6 = rho_new / np.cos(alpha_new)
    rhs6 = (1 + dw) * rho_old / np.cos(alpha_old)
    lhs9 = rho_new * np.tan(alpha_new)
    print(f"dw={dw:+.3f} h={h_t:.5f} : Eq6 |rho/cos a - (1+dw)rho0/cos a0| = {abs(lhs6-rhs6):.3e} ; "
          f"Eq9 |rho tan a - h| = {abs(lhs9-h_t):.3e}")

# nombre de spires du specimen papier
h0 = 2.16 * np.tan(np.deg2rad(10.53))
N = 32.45 / (2 * np.pi * h0)
print(f"h0 = {h0:.5f} mm/rad ; N = L0/(2 pi h0) = {N:.3f} spires (specimen BLOCKED, L=32.45mm)")

# ---------------------------------------------------------------- C
print("\n=== C. Equilibre Eq.(3) vs statique vectorielle 3D ===")
rng = np.random.default_rng(42)
ok = True
for _ in range(5):
    a = rng.uniform(np.deg2rad(3), np.deg2rad(60))
    p = rng.uniform(1.5, 4.0)
    F = rng.uniform(-3, 3)     # force axiale exterieure (N), le long de z
    T = rng.uniform(-3, 3)     # couple axial exterieur (N.mm), autour de z
    # coupe au point P=(p,0,0) ; tangente de l'helice t=(0,cos a, sin a)
    tvec = np.array([0.0, np.cos(a), np.sin(a)])
    nvec = np.array([-1.0, 0.0, 0.0])
    bvec = np.cross(tvec, nvec)
    # torseur exterieur transporte en P : force F ez appliquee sur l'axe (0,0,z0), couple T ez
    z0 = rng.uniform(-5, 5)
    Fv = np.array([0.0, 0.0, F])
    MP = np.array([0.0, 0.0, T]) + np.cross(np.array([0.0, 0.0, z0]) - np.array([p, 0.0, 0.0]), Fv)
    F_wire = float(Fv @ tvec)          # effort normal du fil
    T_wire = float(MP @ tvec)          # couple de torsion du fil
    M_wire = float(MP @ bvec)          # moment de flexion (binormale)
    # Eq. (3) de l'article
    F_eq = F * np.sin(a)
    M_eq = T * np.cos(a) - F * p * np.sin(a)
    T_eq = T * np.sin(a) + F * p * np.cos(a)
    e = max(abs(F_wire - F_eq), abs(M_wire - M_eq), abs(T_wire - T_eq))
    ok &= e < 1e-12
    print(f"a={np.rad2deg(a):5.1f}deg rho={p:.2f} F={F:+.2f} T={T:+.2f} : ecart max Eq3 vs vecteurs = {e:.2e}")
print("Eq.(3) verifiee par statique independante :", ok)

# inversion telle que codee (_trial_state) : Ft=(A)/sin, Tt=(B+Ft rho sin)/cos, residu=3e eq.
A_, B_, C_ = 1.234, -0.567, 0.0
a, p = np.deg2rad(25.0), 2.5
Ft = A_ / np.sin(a)
Tt = (B_ + Ft * p * np.sin(a)) / np.cos(a)
C_ = Tt * np.sin(a) + Ft * p * np.cos(a)  # valeur qui annule le residu
# re-projection : les trois equations doivent etre satisfaites
chk = [abs(A_ - Ft * np.sin(a)), abs(B_ - (Tt * np.cos(a) - Ft * p * np.sin(a))), abs(C_ - (Tt * np.sin(a) + Ft * p * np.cos(a)))]
print("Inversion du code (Ft,Tt) coherente avec les 3 equations :", max(chk) < 1e-12)

# ---------------------------------------------------------------- D
print("\n=== D. Execution modele bloque (petit maillage) : coherence des sorties ===")
disc = Base.default_discretization(n_layers=4, n_phi=12, pre_steps=20, dw_bracket=(-0.05, 0.05))
m = Base.TCPAMaxwellBlockedModel(disc=disc, integration="exponential")
m.prestretch_to(0.8)
t, P = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.3, flow_rate_mL_min=10.0, volume_mL=1.5, dt=1.0, nonlinear=False)
m.step(float(P[0]), 0.0, h_target=m.h_blocked)
m.run_pressure_history(m.helix.time + t, P)
arr = m.history_arrays()

sin_a = np.sin(arr["alpha_rad"]); cos_a = np.cos(arr["alpha_rad"]); rho_a = arr["rho_mm"]
Ftube = arr["Ftube_axis_N"]; Fny = arr["Fnylon_axis_N"]
Mtube = arr["Mtube_Nmm"]; Mny = arr["Mnylon_Nmm"]
Ttube = arr["Ttube_axis_Nmm"]; Tny = arr["Tnylon_axis_Nmm"]
Ft = arr["force_N"]; Tt = arr["torque_Nmm"]
r1 = np.max(np.abs(Ftube + Fny - Ft * sin_a))
r2 = np.max(np.abs(Mtube + Mny - (Tt * cos_a - Ft * rho_a * sin_a)))
r3 = np.max(np.abs(Ttube + Tny - (Tt * sin_a + Ft * rho_a * cos_a)))
print(f"max |F_tube+F_ny - Ft sin a|            = {r1:.3e} N")
print(f"max |M_tube+M_ny - (Tt cos a - Ft rho sin a)| = {r2:.3e} N.mm")
print(f"max |T_tube+T_ny - (Tt sin a + Ft rho cos a)| = {r3:.3e} N.mm  (= residu du solveur)")

# blocage : h et longueur axiale constants apres precontrainte
h_arr = arr["h_mm_per_rad"]
i0 = disc.pre_steps  # fin de precontrainte
print(f"h pendant l'actionnement : min={h_arr[i0:].min():.8f}, max={h_arr[i0:].max():.8f} (h_blocked={m.h_blocked:.8f})")
La = arr["active_axial_length_mm"]; Lg = arr["axial_length_geometry_mm"]
print(f"max |L_active(stretch) - L_active(geometrie)| = {np.max(np.abs(La - (Lg - arr['uncoiled_deformed_length_mm']))):.3e} mm")
print(f"L axiale active en actionnement : min={La[i0:].min():.5f} max={La[i0:].max():.5f} "
      f"(attendu (1+0.8)*32.45={1.8*32.45:.3f})")
# conservation du nombre de tours : rho/cos alpha suit le stretch cumule
lw = rho_a / cos_a
lw_pred = (m.geom.rho0 / np.cos(m.alpha0)) * arr["axial_stretch"]
print(f"max |rho/cos a - (rho0/cos a0)*stretch| = {np.max(np.abs(lw - lw_pred)):.3e} mm  (Eq.6 cumulee)")

# conventions de sortie
arr2 = {k: v.copy() for k, v in arr.items()}
Base.add_corrected_output_conventions(arr2)
print("torque_act[0] =", arr2["torque_act_signed_microNm"][0], " force_act[0] =", arr2["force_act_mN"][0])

# ---------------------------------------------------------------- E
print("\n=== E. Profil d'angle de biais : paper_linear vs Eq.(11) arctan ===")
theta_f = np.deg2rad(37.91)
frac = np.linspace(0.4, 1.0, 7)  # Rin/Rout=0.4 pour la geometrie par defaut
lin = frac * np.rad2deg(theta_f)
arct = np.rad2deg(np.arctan(frac * np.tan(theta_f)))
for f, a_, b_ in zip(frac, lin, arct):
    print(f"R/Rout={f:.2f} : lineaire={a_:6.2f} deg ; arctan(article Eq.11)={b_:6.2f} deg ; ecart={b_-a_:+.2f} deg")

# impact sur une simulation complete
def run_case(profile):
    geom = Base.default_geometry_params(bias_angle_profile=profile)
    disc = Base.default_discretization(n_layers=4, n_phi=12, pre_steps=20, dw_bracket=(-0.05, 0.05))
    mm = Base.TCPAMaxwellBlockedModel(geom=geom, disc=disc, integration="exponential")
    mm.prestretch_to(0.8)
    t, P = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.3, flow_rate_mL_min=10.0, volume_mL=1.5, dt=1.0, nonlinear=False)
    mm.step(float(P[0]), 0.0, h_target=mm.h_blocked)
    i0 = len(mm.history) - 1
    mm.run_pressure_history(mm.helix.time + t, P)
    a = mm.history_arrays()
    return {k: v[i0:] for k, v in a.items()}

try:
    r_lin = run_case("paper_linear")
    r_arc = run_case("uniform_twist")
    fl, fa = r_lin["force_mN"], r_arc["force_mN"]
    tl, ta = r_lin["torque_microNm"], r_arc["torque_microNm"]
    print(f"Pic force  : lineaire={fl.max():8.1f} mN ; arctan={fa.max():8.1f} mN ; ecart={100*(fa.max()-fl.max())/fl.max():+.2f} %")
    print(f"Pic couple : lineaire={tl.max():8.1f} uNm ; arctan={ta.max():8.1f} uNm ; ecart={100*(ta.max()-tl.max())/tl.max():+.2f} %")
    print(f"Force  a P=0 fin : lineaire={fl[-1]:8.1f} ; arctan={fa[-1]:8.1f} mN")
    da = fa - fl
    print(f"ecart RMS force sur le cycle = {np.sqrt(np.mean(da**2)):.1f} mN ; ecart max = {np.max(np.abs(da)):.1f} mN")
    dt_ = ta - tl
    print(f"ecart RMS couple sur le cycle = {np.sqrt(np.mean(dt_**2)):.1f} uNm ; ecart max = {np.max(np.abs(dt_)):.1f} uNm")
except Exception as e:
    print("comparaison profils impossible :", e)

# ---------------------------------------------------------------- F
print("\n=== F. Poissons effectifs : nu(s->phi) vs nu(s->r) apres rotation ===")
C_loc = Base.ti_stiffness_from_paper(37.76, 8.82, 7.24, 0.205, 0.422)
for th_deg in [0.0, 10.0, 20.0, 30.0, 37.91]:
    Cb = Base.rotate_stiffness_bias(C_loc, np.deg2rad(th_deg))
    EL, v12b, v13b, v14b = Base.effective_poissons(Cb)
    # v12b = -eps_phi/eps_s (couplage s->phi), v13b = -eps_r/eps_s (couplage s->r)
    print(f"theta={th_deg:5.2f} deg : EL={EL:7.3f} MPa ; nu(s->phi)={v12b:+.4f} (applique a eps_r) ; "
          f"nu(s->r)={v13b:+.4f} (applique a eps_phi) ; nu(s->sphi)={v14b:+.4f}")
print("Le code applique v12b a eps_r et v13b a eps_phi, litteralement comme Eq.(18)/(19) de l'article.")
