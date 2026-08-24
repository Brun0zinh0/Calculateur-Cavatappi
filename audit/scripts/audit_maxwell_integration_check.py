# -*- coding: utf-8 -*-
"""Audit independant : integration temporelle des branches de Maxwell (Base.py alpha V2).

Tests :
 T1  Relaxation pure (de=0) : integrateur exponentiel vs solution analytique exp(-t/tau),
     avec dt >> tau1 (cas raide).
 T2  Relaxation pure : Euler explicite (paper_explicit) vs (1 - dt/tau)^n et convergence vers l'exact.
 T3  Charge a vitesse de deformation constante : exponentiel vs solution analytique exacte
     sigma(t) = (Ci @ edot) * tau * (1 - exp(-t/tau)).
 T4  Coherence tangente/historique : _algorithmic_data vs _update_maxwell_branches
     (les deux schemas) : dsigma_total == tangent @ de + history.
 T5  Decomposition sans double comptage : C0 + sum(Ci) == C_total par couche,
     modes paper_equal et axial_test_only.
 T6  Valeurs Table A1 et E_axial par defaut ; garde de stabilite dt < 2 tau_min.
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base  # noqa: E402

np.set_printoptions(precision=6, suppress=False)


def make_model(integration, anisotropy="paper_equal"):
    mat = Base.default_material_params(maxwell_anisotropy_mode=anisotropy)
    geom = Base.default_geometry_params()
    disc = Base.default_discretization(n_layers=3, n_phi=8, pre_steps=4)
    return Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration=integration)


model = make_model("exponential")
E = Base._maxwell_E(model.mat.maxwell)
eta = Base._maxwell_eta(model.mat.maxwell)
tau = eta / E
print("tau_i [s] =", tau)
print("E_total (E0+E1+E2+E3) =", Base._maxwell_E_total(model.mat.maxwell))
print("E_axial defaut Base.py =", model.mat.E_axial)

# ------------------------------------------------------------------ T1
print("\n--- T1: relaxation pure, integrateur exponentiel, dt=10 s (tau1=7.48 s) ---")
sig_seed = np.array([1.0, 0.5, -0.3, 0.0, 0.0, 0.2])
for ib in range(model.n_maxwell):
    model.sigma_i[ib, 0, 0] = sig_seed
de0 = np.zeros(6)
dt, nstep = 10.0, 12
for _ in range(nstep):
    model.sigma_i[:, 0, 0] = model._update_maxwell_branches(0, 0, de0, dt)
T = dt * nstep
ok = True
for ib in range(model.n_maxwell):
    exact = sig_seed * np.exp(-T / tau[ib])
    err = np.max(np.abs(model.sigma_i[ib, 0, 0] - exact))
    rel = err / max(np.max(np.abs(exact)), 1e-300)
    print(f"branche {ib+1}: tau={tau[ib]:9.2f} s  err_abs={err:.3e}  err_rel={rel:.3e}")
    ok &= rel < 1e-12
print("T1", "OK (machine precision)" if ok else "ECHEC")

# ------------------------------------------------------------------ T2
print("\n--- T2: relaxation pure, Euler explicite vs (1-dt/tau)^n ---")
for dt in (1.0, 5.0, 14.0):  # 2*tau_min = 14.96 s
    m2 = make_model("paper_explicit")
    for ib in range(m2.n_maxwell):
        m2.sigma_i[ib, 0, 0] = sig_seed
    n2 = max(2, int(round(60.0 / dt)))
    for _ in range(n2):
        m2.sigma_i[:, 0, 0] = m2._update_maxwell_branches(0, 0, de0, dt)
    T2 = dt * n2
    for ib in range(m2.n_maxwell):
        pred = sig_seed * (1.0 - dt / tau[ib]) ** n2
        exact = sig_seed * np.exp(-T2 / tau[ib])
        err_scheme = np.max(np.abs(m2.sigma_i[ib, 0, 0] - pred))
        err_vs_exact = np.max(np.abs(m2.sigma_i[ib, 0, 0] - exact))
        print(f"dt={dt:5.1f}  branche {ib+1}: |code-(1-dt/tau)^n|={err_scheme:.2e}  |code-exact|={err_vs_exact:.3e}")

# ------------------------------------------------------------------ T3
print("\n--- T3: charge a taux constant, exponentiel vs analytique exacte ---")
m3 = make_model("exponential")
edot = np.array([1e-3, 0.0, 0.0, 0.0, 0.0, 5e-4])  # /s
dt, nstep = 2.0, 50
for _ in range(nstep):
    de = edot * dt
    m3.sigma_i[:, 0, 0] = m3._update_maxwell_branches(0, 0, de, dt)
T3 = dt * nstep
okT3 = True
for ib in range(m3.n_maxwell):
    Ci = m3.Ci[0][ib]
    exact = (Ci @ edot) * tau[ib] * (1.0 - np.exp(-T3 / tau[ib]))
    err = np.max(np.abs(m3.sigma_i[ib, 0, 0] - exact))
    rel = err / max(np.max(np.abs(exact)), 1e-300)
    print(f"branche {ib+1}: err_rel={rel:.3e}")
    okT3 &= rel < 1e-10
print("T3", "OK (exact pour taux constant par pas)" if okT3 else "ECHEC")

# ------------------------------------------------------------------ T4
print("\n--- T4: coherence tangente/historique <-> mise a jour des branches ---")
rng = np.random.default_rng(0)
for integ in ("paper_explicit", "exponential"):
    m4 = make_model(integ)
    m4.sigma_i = rng.normal(size=m4.sigma_i.shape) * 0.5
    dt = 1.5
    tangents, history = m4._algorithmic_data(dt)
    worst = 0.0
    for j in range(m4.disc.n_layers):
        for k in (0, 3):
            de = rng.normal(size=6) * 1e-3
            new_branches = m4._update_maxwell_branches(j, k, de, dt)
            dsig_branches = new_branches.sum(axis=0) - m4.sigma_i[:, j, k].sum(axis=0)
            dsig_total = m4.C0[j] @ de + dsig_branches
            pred = tangents[j] @ de + history[j, k]
            worst = max(worst, float(np.max(np.abs(dsig_total - pred))))
    print(f"{integ:15s}: ecart max |dsigma - (C_alg de + hist)| = {worst:.3e} MPa")

# ------------------------------------------------------------------ T5
print("\n--- T5: pas de double comptage C0 + sum Ci == C_total ---")
for mode in ("paper_equal", "axial_test_only"):
    m5 = make_model("exponential", anisotropy=mode)
    worst = 0.0
    for j in range(m5.disc.n_layers):
        S = m5.C0[j] + sum(m5.Ci[j])
        worst = max(worst, float(np.max(np.abs(S - m5.C_total[j]))))
    print(f"{mode:15s}: ecart max ||C0 + sum Ci - C_total|| = {worst:.3e} MPa")
    # module axial d'equilibre (relaxe) de la couche 0 (theta le plus faible)
    ELtot = Base.effective_poissons(m5.C_total[0])[0]
    try:
        EL0 = Base.effective_poissons(m5.C0[0])[0]
        print(f"  E_L instantane couche0 = {ELtot:.3f} MPa ; E_L relaxe (C0) = {EL0:.3f} MPa")
    except np.linalg.LinAlgError:
        print("  C0 singulier pour effective_poissons")

# theta = 0 : verifier module axial local
C_local = Base.ti_stiffness_from_paper(37.76, 8.82, 7.24, 0.205, 0.422)
EL_local = Base.effective_poissons(C_local)[0]
print(f"E_L local (theta=0) avec E_axial=37.76 : {EL_local:.4f} MPa (attendu 37.76)")

# ------------------------------------------------------------------ T6
print("\n--- T6: garde de stabilite explicite ---")
tau_min = float(np.min(tau))
print(f"2*tau_min = {2*tau_min:.4f} s")
m6 = make_model("paper_explicit")
for dt_try in (14.9, 14.96, 15.0):
    try:
        m6._validate_time_step(dt_try)
        print(f"dt={dt_try}: accepte")
    except ValueError as e:
        print(f"dt={dt_try}: rejete ({str(e)[:60]}...)")
m6b = make_model("exponential")
m6b._validate_time_step(1e6)
print("exponential: dt=1e6 s accepte (inconditionnellement stable)")

# ------------------------------------------------------------------ bonus : precision (1-e^-x)/x pour x minuscule
x = 1e-13
naive = (1.0 - np.exp(-x)) / x
print(f"\nfacteur (1-e^-x)/x pour x=1e-13 : {naive!r} (exact ~ 1 - x/2) -> erreur relative {abs(naive-1.0):.2e}")
