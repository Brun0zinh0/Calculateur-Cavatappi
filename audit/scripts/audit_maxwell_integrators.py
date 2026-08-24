# -*- coding: utf-8 -*-
"""Audit independant : integrateurs Maxwell de Base.py vs solutions analytiques.

1. Verification scalaire independante des deux formules (recopiees du code) :
   - relaxation pure d'une branche : sigma(t) = sigma0 * exp(-t/tau)
   - rampe de deformation a taux constant : sigma(t) = eta*rate*(1-exp(-t/tau))
2. Verification directe de TCPAMaxwellBlockedModel._update_maxwell_branches
   (theta_f = 0 pour eviter la rotation) contre les memes solutions.
3. Coherence tangente algorithmique <-> mise a jour des branches.
4. Equivalence matricielle C_i @ inv(eta_i) = (E_i/eta_i) I en mode paper_equal.
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

E1, eta1 = 20.67, 154.57
tau1 = eta1 / E1
print(f"tau1 = {tau1:.6f} s  (tau2={977.79/5.98:.2f} s, tau3={11044.83/4.75:.2f} s)")

# ---------------------------------------------------------------- 1. scalaire
def euler_step(s, E, eta, de, dt):
    # Base.py l.645 : ds = Ci@de - dt*rate*sigma ; rate = E/eta
    return s + E * de - dt * (E / eta) * s

def expo_step(s, E, eta, de, dt):
    # Base.py l.649-654
    tau = eta / E
    a = np.exp(-dt / tau)
    factor = tau / dt * (1.0 - a) if dt > 0 else 1.0
    return a * s + factor * (E * de)

print("\n--- 1a. Relaxation pure (de=0), sigma0=1 MPa, branche 1 ---")
for dt in (0.25, 2.0, 7.478, 20.0):
    T = 60.0
    n = max(1, int(round(T / dt)))
    T = n * dt
    se, sx = 1.0, 1.0
    for _ in range(n):
        se = euler_step(se, E1, eta1, 0.0, dt)
        sx = expo_step(sx, E1, eta1, 0.0, dt)
    exact = np.exp(-T / tau1)
    print(f"dt={dt:7.3f}  T={T:6.1f}  exact={exact:.6e}  expo={sx:.6e} (err {abs(sx-exact):.2e})"
          f"  euler={se:.6e} (err {abs(se-exact):.2e})")

print("\n--- 1b. Rampe eps a taux constant rate=1e-3/s, branche 1 ---")
rate = 1.0e-3
for dt in (0.25, 2.0):
    T = 40.0
    n = int(round(T / dt))
    se, sx = 0.0, 0.0
    for _ in range(n):
        se = euler_step(se, E1, eta1, rate * dt, dt)
        sx = expo_step(sx, E1, eta1, rate * dt, dt)
    exact = eta1 * rate * (1.0 - np.exp(-T / tau1))
    print(f"dt={dt:5.2f}  exact={exact:.8f}  expo={sx:.8f} (err {abs(sx-exact):.2e})"
          f"  euler={se:.8f} (err {abs(se-exact):.2e})")

# ------------------------------------------------- 2. via Base.py directement
print("\n--- 2. _update_maxwell_branches de Base.py (theta_f=0, sans rotation) ---")
geom = Base.default_geometry_params(theta_f_deg=0.0)
mat = Base.default_material_params()
disc = Base.default_discretization(n_layers=2, n_phi=4)

for integ in ("exponential", "paper_explicit"):
    model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration=integ)
    j, k = 0, 0
    dt = 0.25
    n = 240  # 60 s
    # relaxation pure : sigma_i initiale axiale = 1 dans chaque branche
    model.sigma_i[:] = 0.0
    for ib in range(model.n_maxwell):
        model.sigma_i[ib, j, k, 0] = 1.0
    de0 = np.zeros(6)
    for _ in range(n):
        model.sigma_i[:, j, k] = model._update_maxwell_branches(j, k, de0, dt)
    rates = Base._maxwell_rates(mat.maxwell)
    print(f"[{integ}] relaxation 60 s, dt=0.25 :")
    for ib in range(model.n_maxwell):
        exact = np.exp(-60.0 * rates[ib])
        got = model.sigma_i[ib, j, k, 0]
        print(f"   branche {ib+1}: code={got:.8e}  analytique={exact:.8e}  err_rel={abs(got-exact)/exact:.2e}")

# rampe axiale : de = [rate*dt,0,0,0,0,0], comparaison analytique matricielle
model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration="exponential")
j, k = 0, 0
dt, T = 0.25, 40.0
n = int(T / dt)
model.sigma_i[:] = 0.0
de = np.array([rate * dt, 0, 0, 0, 0, 0])
for _ in range(n):
    model.sigma_i[:, j, k] = model._update_maxwell_branches(j, k, de, dt)
rates = Base._maxwell_rates(mat.maxwell)
print("[exponential] rampe axiale 40 s (solution exacte matricielle tau*(1-exp)*Ci@edot) :")
for ib in range(model.n_maxwell):
    tau = 1.0 / rates[ib]
    exact_vec = tau * (1.0 - np.exp(-T / tau)) * (model.Ci[j][ib] @ (de / dt))
    got_vec = model.sigma_i[ib, j, k]
    err = np.max(np.abs(got_vec - exact_vec)) / max(np.max(np.abs(exact_vec)), 1e-30)
    print(f"   branche {ib+1}: err_rel_max = {err:.2e}   sigma_axial code={got_vec[0]:.6f} exact={exact_vec[0]:.6f}")

# ------------------------------------- 3. coherence tangente algorithmique
print("\n--- 3. Coherence _algorithmic_data <-> _update_maxwell_branches ---")
for integ in ("exponential", "paper_explicit"):
    model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration=integ)
    rng = np.random.default_rng(0)
    model.sigma_i[:] = rng.normal(size=model.sigma_i.shape)
    dt = 0.5
    de = rng.normal(scale=1e-3, size=6)
    tangents, history = model._algorithmic_data(dt)
    j, k = 1, 2
    pred = tangents[j] @ de + history[j, k] + model.C0[j] @ de * 0  # tangente inclut C0
    # delta reel = somme des increments de branches + C0@de
    new_branches = model._update_maxwell_branches(j, k, de, dt)
    real = (new_branches - model.sigma_i[:, j, k]).sum(axis=0) + model.C0[j] @ de
    print(f"[{integ}] max|pred - reel| = {np.max(np.abs(pred - real)):.3e}")

# ------------------------------------- 4. equivalence matricielle paper_equal
print("\n--- 4. C_i @ inv(eta_i) = (E_i/eta_i) I sous paper_equal (theta=25 deg) ---")
Cbar = Base.rotate_stiffness_bias(model.C_local_total, np.deg2rad(25.0))
Etot = Base._maxwell_E_total(mat.maxwell)
Es = Base._maxwell_E(mat.maxwell)
etas = Base._maxwell_eta(mat.maxwell)
for ib in range(3):
    Ci_b = (Es[ib] / Etot) * Cbar
    eta_b = (etas[ib] / Etot) * Cbar
    M = Ci_b @ np.linalg.inv(eta_b)
    err = np.max(np.abs(M - (Es[ib] / etas[ib]) * np.eye(6)))
    print(f"   branche {ib+1}: max|Ci@inv(eta_i) - (E_i/eta_i) I| = {err:.3e}")

# ------------------------------------- 5. garde de stabilite explicite
print("\n--- 5. dt limite explicite : 2*tau_min =", 2.0 * min(etas / Es), "s")
try:
    model = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration="paper_explicit")
    model._validate_time_step(15.0)
    print("   dt=15 s ACCEPTE (inattendu)")
except ValueError as e:
    print("   dt=15 s rejete :", e)
model2 = Base.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc, integration="exponential")
model2._validate_time_step(15.0)
print("   dt=15 s accepte en exponentiel (inconditionnellement stable) : OK")
