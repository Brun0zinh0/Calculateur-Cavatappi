# -*- coding: utf-8 -*-
"""Audit independant : raideur TI, rotation de Voigt, profil de biais, Poisson effectifs."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

np.set_printoptions(precision=6, suppress=True, linewidth=140)

E1, E2, G12, nu12, nu23 = 37.76, 8.82, 7.24, 0.205, 0.422

# ---------------------------------------------------------------------------
# 1. Raideur TI exacte par inversion de la souplesse (reference independante)
# ---------------------------------------------------------------------------
S = np.zeros((6, 6))
S[0, 0] = 1.0 / E1
S[1, 1] = S[2, 2] = 1.0 / E2
S[0, 1] = S[1, 0] = S[0, 2] = S[2, 0] = -nu12 / E1
S[1, 2] = S[2, 1] = -nu23 / E2
G23 = E2 / (2.0 * (1.0 + nu23))
S[3, 3] = 1.0 / G23
S[4, 4] = S[5, 5] = 1.0 / G12
C_exact = np.linalg.inv(S)

C_code = Base.ti_stiffness_from_paper(E1, E2, G12, nu12, nu23)
print("=== 1. ti_stiffness_from_paper vs inversion souplesse exacte ===")
print("C_code =\n", C_code)
print("ecart max |C_code - C_exact| =", np.max(np.abs(C_code - C_exact)))
print("valeurs propres C_code :", np.linalg.eigvalsh(C_code))
print("C44 = (C22-C23)/2 =", C_code[3, 3], " vs G23 exact =", G23)

# Verification des 4 identites imprimees (A.2), avec C exact :
C11, C12, C22, C23 = C_exact[0, 0], C_exact[0, 1], C_exact[1, 1], C_exact[1, 2]
nu21 = nu12 * E2 / E1
print("\nIdentites A.2 evaluees avec la raideur exacte :")
print("l1: C11 - (E_axial + 2 nu12 C12)      =", C11 - (E1 + 2 * nu12 * C12))
print("l2: C12 - nu12 (C22 + C23)            =", C12 - nu12 * (C22 + C23))
print("l3: C12 - (C11 nu21 + C12 nu23)       =", C12 - (C11 * nu21 + C12 * nu23))
print("l4: C22 - (C12 nu21 + C23 nu23 + E2)  =", C22 - (C12 * nu21 + C23 * nu23 + E2))

# ---------------------------------------------------------------------------
# 2. Rotation : code (einsum) vs formules standard fermees (A.3 corrigee)
# ---------------------------------------------------------------------------
def rot_standard(C, th):
    """Rotation standard d'une raideur orthotrope autour de l'axe 3 (Voigt ing.),
    formules fermees classiques (laminates), m=cos, n=sin."""
    m, n = np.cos(th), np.sin(th)
    C11, C12, C13 = C[0, 0], C[0, 1], C[0, 2]
    C22, C23, C33 = C[1, 1], C[1, 2], C[2, 2]
    C44, C55, C66 = C[3, 3], C[4, 4], C[5, 5]
    R = np.zeros((6, 6))
    R[0, 0] = C11 * m**4 + 2 * (C12 + 2 * C66) * m**2 * n**2 + C22 * n**4
    R[0, 1] = R[1, 0] = (C11 + C22 - 4 * C66) * m**2 * n**2 + C12 * (m**4 + n**4)
    R[0, 2] = R[2, 0] = C13 * m**2 + C23 * n**2
    R[0, 5] = R[5, 0] = -C22 * m * n**3 + C11 * m**3 * n - (C12 + 2 * C66) * m * n * (m**2 - n**2)
    R[1, 1] = C11 * n**4 + 2 * (C12 + 2 * C66) * m**2 * n**2 + C22 * m**4
    R[1, 2] = R[2, 1] = C13 * n**2 + C23 * m**2
    R[1, 5] = R[5, 1] = -C22 * m**3 * n + C11 * m * n**3 + (C12 + 2 * C66) * m * n * (m**2 - n**2)
    R[2, 2] = C33
    R[2, 5] = R[5, 2] = (C13 - C23) * m * n
    R[5, 5] = (C11 + C22 - 2 * C12) * m**2 * n**2 + C66 * (m**2 - n**2) ** 2
    # bloc cisaillement transverse 4-5
    R[3, 3] = C44 * m**2 + C55 * n**2
    R[4, 4] = C44 * n**2 + C55 * m**2
    R[3, 4] = R[4, 3] = (C55 - C44) * m * n
    return R

def rot_A3_printed(C, th):
    """A.3 telle qu'imprimee : C16 (=0) dans Cbar11, C44=(C22-C23)/2 dans Cbar22."""
    m, n = np.cos(th), np.sin(th)
    C11, C12, C13 = C[0, 0], C[0, 1], C[0, 2]
    C22, C23, C33 = C[1, 1], C[1, 2], C[2, 2]
    C44, C66 = C[3, 3], C[5, 5]
    C16 = 0.0
    R = np.zeros((6, 6))
    R[0, 0] = C11 * m**4 + 2 * (C12 + C16) * m**2 * n**2 + C22 * n**4
    R[1, 1] = C11 * n**4 + 2 * (C12 + 2 * C44) * m**2 * n**2 + C22 * m**4
    return R

print("\n=== 2. rotate_stiffness_bias vs formules standard ===")
for th_deg in (10.0, 20.0, 37.91, 55.0):
    th = np.deg2rad(th_deg)
    Cb_code = Base.rotate_stiffness_bias(C_code, th)
    Cb_std = rot_standard(C_code, th)
    print(f"theta={th_deg:6.2f} deg : ecart max code vs standard =", np.max(np.abs(Cb_code - Cb_std)))

th = np.deg2rad(37.91)
Cb_code = Base.rotate_stiffness_bias(C_code, th)
Cb_std = rot_standard(C_code, th)
Cb_prn = rot_A3_printed(C_code, th)
print("\nA theta=37.91 deg :")
print("Cbar11 code/standard/imprimee :", Cb_code[0, 0], Cb_std[0, 0], Cb_prn[0, 0])
print("Cbar22 code/standard/imprimee :", Cb_code[1, 1], Cb_std[1, 1], Cb_prn[1, 1])
print("Cbar16 code (signe !)          :", Cb_code[0, 5], " standard :", Cb_std[0, 5])
print("Cbar26 code                    :", Cb_code[1, 5], " standard :", Cb_std[1, 5])
print("Cbar36 code                    :", Cb_code[2, 5], " standard :", Cb_std[2, 5])
print("Cbar (code) complet :\n", Cb_code)

# roundtrip Voigt
Cv = Base.tensor_to_voigt(Base.voigt_to_tensor(C_code))
print("\nRoundtrip voigt_to_tensor/tensor_to_voigt, ecart :", np.max(np.abs(Cv - C_code)))

# invariance isotrope
rep = Base.stiffness_sanity_report(0.7)
print("Sanity report Base :", rep)

# ---------------------------------------------------------------------------
# 3. Poisson effectifs : code vs inversion complete 6x6, et question du swap
# ---------------------------------------------------------------------------
print("\n=== 3. effective_poissons ===")
for th_deg in (10.0, 20.0, 37.91):
    th = np.deg2rad(th_deg)
    Cb = Base.rotate_stiffness_bias(C_code, th)
    EL, v12b, v13b, v14b = Base.effective_poissons(Cb)
    # reference : inversion 6x6 complete puis lecture des memes rapports
    S6 = np.linalg.inv(Cb)
    ELs = 1.0 / S6[0, 0]
    eps = S6 @ np.array([ELs, 0, 0, 0, 0, 0.0])
    print(f"theta={th_deg:6.2f} : EL={EL:9.4f}  vbar12(phi)={v12b:8.5f}  vbar13(r)={v13b:8.5f}  vbar14(sphi)={v14b:8.5f}")
    print(f"             ref 6x6 : EL={ELs:9.4f}  -eps_phi={-eps[1]:8.5f}  -eps_r={-eps[2]:8.5f}  -gam_sphi={-eps[5]:8.5f}")
    print(f"             ecart swap |v12b - v13b| = {abs(v12b - v13b):.5f}  (rel {abs(v12b - v13b)/max(abs(v12b), abs(v13b)):.2%})")

# ---------------------------------------------------------------------------
# 4. Profil radial de l'angle de biais : lineaire en theta vs Eq. (11)
# ---------------------------------------------------------------------------
print("\n=== 4. Profil de biais : paper_linear (defaut) vs Eq. (11) arctan ===")
theta_f = np.deg2rad(37.91)
Rin, Rout = 0.4, 1.0
n_layers = 8
edges = np.linspace(Rin, Rout, n_layers + 1)
centers = 0.5 * (edges[:-1] + edges[1:])
frac = centers / Rout
th_lin = frac * theta_f                      # mode 'paper_linear' (defaut du code)
th_eq11 = np.arctan(frac * np.tan(theta_f))  # Eq. (11) de l'article = mode 'uniform_twist'
print("R/Rout    theta_lin(deg)  theta_Eq11(deg)  ecart(deg)")
for f, a, b in zip(frac, np.rad2deg(th_lin), np.rad2deg(th_eq11)):
    print(f"{f:6.3f}   {a:12.3f}   {b:13.3f}   {b - a:9.3f}")

# impact sur la raideur tournee et les coefficients du BVP radial
print("\nImpact sur Cbar et (A, B, mu) par couche :")
print("R/Rout   dC11/C11   dC16/C16   dEL/EL     dA/A      dmu/mu")
for f, ta, tb in zip(frac, th_lin, th_eq11):
    Ca = Base.rotate_stiffness_bias(C_code, float(ta))
    Cb = Base.rotate_stiffness_bias(C_code, float(tb))
    ELa = Base.effective_poissons(Ca)[0]
    ELb = Base.effective_poissons(Cb)[0]
    def ABmu(C):
        C12, C13 = C[0, 1], C[0, 2]
        C22, C26 = C[1, 1], C[1, 5]
        C33, C36 = C[2, 2], C[2, 5]
        mu = np.sqrt(C22 / C33)
        A = (C26 - 2 * C36) / (4 * C33 - C22)
        B = (C12 - C13) / (C33 - C22)
        return A, B, mu
    Aa, Ba, mua = ABmu(Ca)
    Ab, Bb, mub = ABmu(Cb)
    print(f"{f:6.3f}  {(Cb[0,0]-Ca[0,0])/Ca[0,0]:+8.2%}  {(Cb[0,5]-Ca[0,5])/abs(Ca[0,5]) if Ca[0,5] != 0 else float('nan'):+8.2%}"
          f"  {(ELb-ELa)/ELa:+8.2%}  {(Ab-Aa)/abs(Aa):+8.2%}  {(mub-mua)/mua:+8.2%}")

# ---------------------------------------------------------------------------
# 5. Impact global des deux profils sur une simulation bloquee courte
# ---------------------------------------------------------------------------
print("\n=== 5. Simulation bloquee courte : paper_linear vs uniform_twist ===")
results = {}
for profile in ("paper_linear", "uniform_twist"):
    geom = Base.default_geometry_params()
    geom.bias_angle_profile = profile
    model, arr = Base.run_blocked_actuation(
        eps=0.8, n_cycles=1, Pmax=1.4, dt=0.5,
        n_layers=8, n_phi=24, pre_steps=24,
        integration="exponential", geom=geom,
    )
    fpk = float(np.max(arr["force_mN"]))
    tpk = float(np.max(arr["torque_microNm"]))
    f0 = float(arr["force_mN"][0])
    t0 = float(arr["torque_microNm"][0])
    results[profile] = (f0, fpk, tpk)
    print(f"{profile:14s} : F(0)={f0:9.2f} mN  F_pic={fpk:9.2f} mN  T_pic={tpk:9.2f} uNm")
a, b = results["paper_linear"], results["uniform_twist"]
print(f"Ecarts relatifs (Eq.11 vs defaut) : F0 {(b[0]-a[0])/a[0]:+.2%}  Fpic {(b[1]-a[1])/a[1]:+.2%}  Tpic {(b[2]-a[2])/a[2]:+.2%}")

# ---------------------------------------------------------------------------
# 6. Signe de la rotation : orientation de la fibre pour theta > 0
# ---------------------------------------------------------------------------
print("\n=== 6. Convention de signe de la rotation ===")
th = np.deg2rad(30.0)
q = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
e1_local_in_global = q @ np.array([1.0, 0, 0])
print("Direction fibre (axe local 1) exprimee dans (s, phi, r) :", e1_local_in_global)
Cb30 = Base.rotate_stiffness_bias(C_code, th)
print("Cbar16 =", Cb30[0, 5], "; Cbar26 =", Cb30[1, 5], "; Cbar36 =", Cb30[2, 5])
m, n = np.cos(th), np.sin(th)
C16_printed = -C22 * m * n**3 + C11 * m**3 * n - (C12 + 2 * G12) * m * n * (m**2 - n**2)
print("Cbar16 formule imprimee A.3 (n=+sin) =", C16_printed)
