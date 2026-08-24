# -*- coding: utf-8 -*-
"""Verification independante du solveur radial (Pipes & Hubert par couches).

On compare la solution du code (u = C1 R^mu + C2 R^-mu + A dv R^2 + B dw R,
continuite + CL en pression) a une resolution par differences finies de
l'equilibre radial axisymetrique :

    d(sigma_r)/dR + (sigma_r - sigma_phi)/R = 0
    sigma = Cbar @ [dw, u/R, u', 0, 0, dv R]

pour une couche uniforme anisotrope (angle de biais 20 deg partout),
avec sigma_r(Rin) = -dP et sigma_r(Rout) = 0.
"""
import sys
import numpy as np

sys.path.insert(
    0,
    r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2",
)
import Base

mat = Base.default_material_params()
geom = Base.default_geometry_params()
disc = Base.default_discretization(n_layers=6, n_phi=8)
model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc)

theta = np.deg2rad(20.0)
Cbar = Base.rotate_stiffness_bias(model.C_local_total, theta)
C_uniform = [Cbar.copy() for _ in range(disc.n_layers)]

Rin, Rout = geom.Rin, geom.Rout
dP = 1.0

for dv, dw, label in ((0.0, 0.0, "dv=0, dw=0"), (0.02, -0.01, "dv=0.02, dw=-0.01")):
    hist = np.zeros((disc.n_layers, disc.n_phi, 6))
    coeffs = model._solve_radial_constants(dw, dv, dP, 0.0, C_uniform, hist)

    # ---------------- differences finies independantes ----------------
    N = 4001
    R = np.linspace(Rin, Rout, N)
    hgrid = R[1] - R[0]
    C22b, C23b, C33b = Cbar[1, 1], Cbar[1, 2], Cbar[2, 2]
    C12b, C13b = Cbar[0, 1], Cbar[0, 2]
    C26b, C36b = Cbar[1, 5], Cbar[2, 5]
    # sigma_r  = C13 dw + C23 u/R + C33 u' + C36 dv R
    # sigma_ph = C12 dw + C22 u/R + C23 u' + C26 dv R
    # ODE: C33 u'' + C33 u'/R - C22 u/R^2 + (C13-C12) dw/R + (2 C36 - C26) dv = 0
    main = np.zeros((N, N))
    rhs = np.zeros(N)
    for i in range(1, N - 1):
        Ri = R[i]
        main[i, i - 1] = C33b / hgrid**2 - C33b / (2 * hgrid * Ri)
        main[i, i] = -2 * C33b / hgrid**2 - C22b / Ri**2
        main[i, i + 1] = C33b / hgrid**2 + C33b / (2 * hgrid * Ri)
        rhs[i] = -(C13b - C12b) * dw / Ri - (2 * C36b - C26b) * dv
    # CL sigma_r(Rin) = -dP (difference decentree ordre 2)
    main[0, 0] = C23b / Rin + C33b * (-3.0) / (2 * hgrid)
    main[0, 1] = C33b * 4.0 / (2 * hgrid)
    main[0, 2] = C33b * (-1.0) / (2 * hgrid)
    rhs[0] = -dP - C13b * dw - C36b * dv * Rin
    # CL sigma_r(Rout) = 0
    main[N - 1, N - 1] = C23b / Rout + C33b * (3.0) / (2 * hgrid)
    main[N - 1, N - 2] = C33b * (-4.0) / (2 * hgrid)
    main[N - 1, N - 3] = C33b * (1.0) / (2 * hgrid)
    rhs[N - 1] = -C13b * dw - C36b * dv * Rout
    u_fd = np.linalg.solve(main, rhs)

    print(f"\n=== Cas {label} ===")
    print("   R       u(code)        u(FD)         sigma_r(code)  sigma_r(FD)")
    max_du = 0.0
    max_ds = 0.0
    for j, Rc in enumerate(model.R_centers):
        u_c, du_c = model._u_du_layer(j, Rc, coeffs, dv, dw, C_uniform)
        eps = np.array([dw, u_c / Rc, du_c, 0.0, 0.0, dv * Rc])
        s_c = Cbar @ eps
        i = np.argmin(np.abs(R - Rc))
        ui = u_fd[i]
        dui = (u_fd[i + 1] - u_fd[i - 1]) / (2 * hgrid)
        eps_fd = np.array([dw, ui / R[i], dui, 0.0, 0.0, dv * R[i]])
        s_fd = Cbar @ eps_fd
        print(f"  {Rc:5.3f}  {u_c:+.8f}  {ui:+.8f}   {s_c[2]:+.6f}     {s_fd[2]:+.6f}")
        max_du = max(max_du, abs(u_c - ui))
        max_ds = max(max_ds, abs(s_c[2] - s_fd[2]))
    print(f"  ecart max |u| : {max_du:.3e} mm ; |sigma_r| : {max_ds:.3e} MPa")

    # CL du code aux faces
    Aj, Bj, muj = model._layer_AB_mu(0, C_uniform)
    c_sig, k_sig = model._radial_stress_axisym_coeff(0, Rin, Aj, Bj, muj, dv, dw, 0.0, C_uniform)
    sr_in = float(c_sig @ coeffs[0]) + k_sig
    Aj, Bj, muj = model._layer_AB_mu(disc.n_layers - 1, C_uniform)
    c_sig, k_sig = model._radial_stress_axisym_coeff(
        disc.n_layers - 1, Rout, Aj, Bj, muj, dv, dw, 0.0, C_uniform
    )
    sr_out = float(c_sig @ coeffs[-1]) + k_sig
    print(f"  CL code : sigma_r(Rin)={sr_in:+.8f} (attendu {-dP:+.1f}) ; sigma_r(Rout)={sr_out:+.2e}")

# ------- continuite aux interfaces avec les couches reelles (angles varies) -
print("\n=== Continuite inter-couches (profil par defaut, angles varies) ===")
C_real = model.C_total
dv, dw = 0.015, -0.008
hist = np.zeros((disc.n_layers, disc.n_phi, 6))
coeffs = model._solve_radial_constants(dw, dv, dP, 0.0, C_real, hist)
for j in range(disc.n_layers - 1):
    Rb = model.R_edges[j + 1]
    u_j, du_j = model._u_du_layer(j, Rb, coeffs, dv, dw, C_real)
    u_k, du_k = model._u_du_layer(j + 1, Rb, coeffs, dv, dw, C_real)
    eps_j = np.array([dw, u_j / Rb, du_j, 0.0, 0.0, dv * Rb])
    eps_k = np.array([dw, u_k / Rb, du_k, 0.0, 0.0, dv * Rb])
    s_j = (C_real[j] @ eps_j)[2]
    s_k = (C_real[j + 1] @ eps_k)[2]
    print(f"  interface R={Rb:.3f} : du={u_j - u_k:+.2e} mm ; dsigma_r={s_j - s_k:+.2e} MPa")
Aj, Bj, muj = model._layer_AB_mu(0, C_real)
c_sig, k_sig = model._radial_stress_axisym_coeff(0, Rin, Aj, Bj, muj, dv, dw, 0.0, C_real)
print(f"  sigma_r(Rin) = {float(c_sig @ coeffs[0]) + k_sig:+.8f} (attendu -1)")
Aj, Bj, muj = model._layer_AB_mu(disc.n_layers - 1, C_real)
c_sig, k_sig = model._radial_stress_axisym_coeff(disc.n_layers - 1, Rout, Aj, Bj, muj, dv, dw, 0.0, C_real)
print(f"  sigma_r(Rout) = {float(c_sig @ coeffs[-1]) + k_sig:+.2e} (attendu 0)")
