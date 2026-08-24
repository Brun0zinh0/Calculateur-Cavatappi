# -*- coding: utf-8 -*-
"""Audit independant : raideur TI, rotation de Voigt, Poissons effectifs, profil de biais."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

np.set_printoptions(precision=5, suppress=True, linewidth=140)

E_ax, E_rad, G12, nu12, nu23 = 37.76, 8.82, 7.24, 0.205, 0.422

# ---------------------------------------------------------------- 1. raideur TI
def C_from_compliance(E1, E2, G12, nu12, nu23):
    """Inversion directe de la compliance TI standard (axe 1 = fibre), Voigt engineering."""
    S = np.zeros((6, 6))
    S[0, 0] = 1.0 / E1
    S[1, 1] = S[2, 2] = 1.0 / E2
    S[0, 1] = S[1, 0] = S[0, 2] = S[2, 0] = -nu12 / E1
    S[1, 2] = S[2, 1] = -nu23 / E2
    S[3, 3] = 2.0 * (1.0 + nu23) / E2   # gamma_23 = tau/G23
    S[4, 4] = S[5, 5] = 1.0 / G12
    return np.linalg.inv(S)

C_code = Base.ti_stiffness_from_paper(E_ax, E_rad, G12, nu12, nu23)
C_ref = C_from_compliance(E_ax, E_rad, G12, nu12, nu23)
print("=== 1. ti_stiffness_from_paper vs inversion compliance standard ===")
print("ecart max |C_code - C_ref| =", np.max(np.abs(C_code - C_ref)))
print("C_code diag :", np.diag(C_code))
print("C11,C12,C22,C23 =", C_code[0,0], C_code[0,1], C_code[1,1], C_code[1,2])
print("C44 code =", C_code[3,3], " vs G23=E2/(2(1+nu23)) =", E_rad/(2*(1+nu23)))
print("valeurs propres :", np.linalg.eigvalsh(C_code))
C_code2 = Base.ti_stiffness_from_paper(31.24, E_rad, G12, nu12, nu23)
C_ref2 = C_from_compliance(31.24, E_rad, G12, nu12, nu23)
print("E_axial=31.24 : ecart max =", np.max(np.abs(C_code2 - C_ref2)),
      " vp min =", np.linalg.eigvalsh(C_code2).min())

# ---------------------------------------------------------------- 2. rotation
def rot_closed_form(C, th):
    """Formules standard (Jones) de rotation autour de l'axe 3, angle +theta."""
    m, n = np.cos(th), np.sin(th)
    C11, C12, C13, C22, C23, C33 = C[0,0], C[0,1], C[0,2], C[1,1], C[1,2], C[2,2]
    C44, C55, C66 = C[3,3], C[4,4], C[5,5]
    Cb = np.zeros((6,6))
    Cb[0,0] = C11*m**4 + 2*(C12 + 2*C66)*m**2*n**2 + C22*n**4
    Cb[0,1] = (C11 + C22 - 4*C66)*m**2*n**2 + C12*(m**4 + n**4)
    Cb[0,2] = C13*m**2 + C23*n**2
    Cb[0,5] = (C11 - C12 - 2*C66)*m**3*n + (C12 - C22 + 2*C66)*m*n**3
    Cb[1,1] = C11*n**4 + 2*(C12 + 2*C66)*m**2*n**2 + C22*m**4
    Cb[1,2] = C13*n**2 + C23*m**2
    Cb[1,5] = (C11 - C12 - 2*C66)*m*n**3 + (C12 - C22 + 2*C66)*m**3*n
    Cb[2,2] = C33
    Cb[2,5] = (C13 - C23)*m*n
    Cb[5,5] = (C11 + C22 - 2*C12)*m**2*n**2 + C66*(m**2 - n**2)**2
    Cb[3,3] = C44*m**2 + C55*n**2
    Cb[4,4] = C55*m**2 + C44*n**2
    Cb[3,4] = (C55 - C44)*m*n
    for i in range(6):
        for j in range(i):
            Cb[i,j] = Cb[j,i]
    return Cb

def rot_paper_printed(C, th):
    """(A.3) telle qu'imprimee : 'C16' (=0 dans A.1) dans Cb11, C44 dans Cb22."""
    m, n = np.cos(th), np.sin(th)
    C11, C12, C22 = C[0,0], C[0,1], C[1,1]
    C44 = C[3,3]
    C16 = 0.0
    Cb = rot_closed_form(C, th).copy()
    Cb[0,0] = C11*m**4 + 2*(C12 + C16)*m**2*n**2 + C22*n**4
    Cb[1,1] = C11*n**4 + 2*(C12 + 2*C44)*m**2*n**2 + C22*m**4
    Cb[1,0] = Cb[0,1]
    return Cb

th = np.deg2rad(37.91)
Cb_code_p = Base.rotate_stiffness_bias(C_code, th)
Cb_code_m = Base.rotate_stiffness_bias(C_code, -th)
Cb_std = rot_closed_form(C_code, th)
Cb_paper = rot_paper_printed(C_code, th)
print("\n=== 2. rotation : code (tenseur) vs forme fermee standard (+theta) ===")
print("ecart max |code(+th) - std(+th)| =", np.max(np.abs(Cb_code_p - Cb_std)))
print("ecart max |code(-th) - std(+th)| =", np.max(np.abs(Cb_code_m - Cb_std)))
print("C16b, C26b, C36b code(+th) :", Cb_code_p[0,5], Cb_code_p[1,5], Cb_code_p[2,5])
print("ecart max |code(+th) - A.3 imprimee| =", np.max(np.abs(Cb_code_p - Cb_paper)))
print("  detail C11b : code", Cb_code_p[0,0], " A.3 imprimee", Cb_paper[0,0])
print("  detail C22b : code", Cb_code_p[1,1], " A.3 imprimee", Cb_paper[1,1])

def T_sigma(th):
    m, n = np.cos(th), np.sin(th)
    return np.array([
        [m*m, n*n, 0, 0, 0, -2*m*n],
        [n*n, m*m, 0, 0, 0,  2*m*n],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, m,  n, 0],
        [0, 0, 0, -n, m, 0],
        [m*n, -m*n, 0, 0, 0, m*m - n*n],
    ])

def T_eps(th):
    m, n = np.cos(th), np.sin(th)
    return np.array([
        [m*m, n*n, 0, 0, 0, -m*n],
        [n*n, m*m, 0, 0, 0,  m*n],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, m,  n, 0],
        [0, 0, 0, -n, m, 0],
        [2*m*n, -2*m*n, 0, 0, 0, m*m - n*n],
    ])

Cb_transform = T_sigma(th) @ C_code @ np.linalg.inv(T_eps(th))
print("verif Reuter : ecart |code(+th) - T_sig C T_eps^-1| =", np.max(np.abs(Cb_code_p - Cb_transform)))
rep = Base.stiffness_sanity_report()
print("sanity report Base :", rep)

# ---------------------------------------------------------------- 3. Poissons effectifs
print("\n=== 3. effective_poissons ===")
EL, v12b, v13b, v14b = Base.effective_poissons(Cb_code_p)
print(f"theta=37.91 deg : EL={EL:.4f}  v12b(s->phi)={v12b:.4f}  v13b(s->r)={v13b:.4f}  v14b(s->sphi)={v14b:.4f}")
S6 = np.linalg.inv(Cb_code_p)
eps_uni = S6 @ np.array([1.0, 0, 0, 0, 0, 0])
print("via S6 : EL =", 1/eps_uni[0], " nu(s->phi) =", -eps_uni[1]/eps_uni[0],
      " nu(s->r) =", -eps_uni[2]/eps_uni[0], " nu(s->sphi) =", -eps_uni[5]/eps_uni[0])
print("couplages residuels eps4, eps5 :", eps_uni[3], eps_uni[4])
for deg in (10, 20, 30, 37.91):
    Cb = Base.rotate_stiffness_bias(C_code, np.deg2rad(deg))
    _, a, b, c = Base.effective_poissons(Cb)
    print(f"  theta={deg:6.2f} deg : nu(s->phi)={a:.4f}  nu(s->r)={b:.4f}  ecart rel={(a-b)/b*100:6.1f}%  nu(s->sphi)={c:.4f}")

# ---------------------------------------------------------------- 4. profil de biais
print("\n=== 4. profil radial de l'angle de biais (n_layers=4, Rin=0.4, Rout=1.0) ===")
theta_f = 37.91
edges = np.linspace(0.4, 1.0, 5)
centers = 0.5*(edges[:-1] + edges[1:])
lin = centers/1.0*theta_f
arct = np.rad2deg(np.arctan(centers/1.0*np.tan(np.deg2rad(theta_f))))
print("R_centres        :", centers)
print("lineaire (code)  :", lin)
print("arctan  (article):", arct)
print("ecart (deg)      :", lin - arct)

print("\n--- run_blocked_actuation : paper_linear vs uniform_twist (=Eq.11 article) ---")
res = {}
for prof in ("paper_linear", "uniform_twist"):
    geom = Base.default_geometry_params(bias_angle_profile=prof)
    model, arr = Base.run_blocked_actuation(eps=0.8, n_cycles=2, Pmax=1.3, dt=0.5,
                                            n_layers=4, n_phi=24, pre_steps=24,
                                            geom=geom)
    res[prof] = arr
    print(f"{prof:14s} : F_max={np.nanmax(arr['force_total_mN']):9.2f} mN  "
          f"F_min={np.nanmin(arr['force_total_mN']):9.2f} mN  "
          f"T_act_max={np.nanmax(arr['torque_act_microNm']):9.2f} uNm  "
          f"residu_max={np.nanmax(np.abs(arr['residual'])):.2e}")
a, b = res["paper_linear"], res["uniform_twist"]
n = min(len(a["force_total_mN"]), len(b["force_total_mN"]))
dF = np.nanmax(np.abs(a["force_total_mN"][:n] - b["force_total_mN"][:n]))
dT = np.nanmax(np.abs(a["torque_act_microNm"][:n] - b["torque_act_microNm"][:n]))
Fmax = np.nanmax(np.abs(b["force_total_mN"][:n]))
Tmax = np.nanmax(np.abs(b["torque_act_microNm"][:n]))
print(f"ecart max force  : {dF:.2f} mN  ({100*dF/Fmax:.2f} % du max)")
print(f"ecart max couple : {dT:.2f} uNm ({100*dT/Tmax:.2f} % du max)")

for prof in ("paper_linear", "uniform_twist"):
    geom = Base.default_geometry_params(bias_angle_profile=prof)
    model, arr = Base.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=1.0,
                                            n_layers=12, n_phi=24, pre_steps=12,
                                            geom=geom)
    print(f"n12 {prof:14s} : F_max={np.nanmax(arr['force_total_mN']):9.2f} mN  "
          f"T_act_max={np.nanmax(arr['torque_act_microNm']):9.2f} uNm")

# ---------------------------------------------------------------- 5. echange v12b/v13b : impact
print("\n=== 5. impact numerique de l'echange v12b <-> v13b (patch temporaire) ===")
orig = Base.effective_poissons
def swapped(Cbar):
    EL, va, vb, vc = orig(Cbar)
    return EL, vb, va, vc
Base.effective_poissons = swapped
try:
    model, arr_sw = Base.run_blocked_actuation(eps=0.8, n_cycles=2, Pmax=1.3, dt=0.5,
                                               n_layers=4, n_phi=24, pre_steps=24)
finally:
    Base.effective_poissons = orig
model, arr_ref = Base.run_blocked_actuation(eps=0.8, n_cycles=2, Pmax=1.3, dt=0.5,
                                            n_layers=4, n_phi=24, pre_steps=24)
n = min(len(arr_sw["force_total_mN"]), len(arr_ref["force_total_mN"]))
dF = np.nanmax(np.abs(arr_sw["force_total_mN"][:n] - arr_ref["force_total_mN"][:n]))
dT = np.nanmax(np.abs(arr_sw["torque_act_microNm"][:n] - arr_ref["torque_act_microNm"][:n]))
print(f"swap Poisson : dF_max={dF:.3f} mN ({100*dF/np.nanmax(np.abs(arr_ref['force_total_mN'][:n])):.3f} %)  "
      f"dT_max={dT:.3f} uNm ({100*dT/np.nanmax(np.abs(arr_ref['torque_act_microNm'][:n])):.3f} %)")
