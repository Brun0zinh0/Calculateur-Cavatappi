# -*- coding: utf-8 -*-
"""Audit independant : viscoelasticite de Maxwell generalisee dans Base.py.

Verifie :
 1. tau_i implicites de la Table A1 ;
 2. relaxation analytique d'une branche seule vs integrateur exponentiel DU CODE ;
 3. rampe a taux constant : exactitude nodale de l'exponentiel, erreur O(dt) d'Euler ;
 4. coherence _algorithmic_data (tangente + historique) vs _update_maxwell_branches ;
 5. absence de double comptage : C0 + somme(Ci) == C_total par couche ;
 6. garde de stabilite d'Euler explicite (dt >= 2 tau_min doit lever) ;
 7. degenerescence propre du pas dt = 0 (reponse elastique instantanee).
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

ok = True

def check(name, cond, detail=""):
    global ok
    status = "OK " if cond else "FAIL"
    if not cond:
        ok = False
    print(f"[{status}] {name} {detail}")

# ---------------------------------------------------------------- 1. tau_i
mw = Base.default_maxwell_tensile_params()
E = np.array([mw.E1, mw.E2, mw.E3]); eta = np.array([mw.eta1, mw.eta2, mw.eta3])
tau = eta / E
print("tau_i [s] =", tau)
check("tau1 ~ 7.478 s", abs(tau[0] - 154.57/20.67) < 1e-12)

# Modele minimal pour acceder aux methodes internes
disc = Base.default_discretization(n_layers=2, n_phi=4, pre_steps=2, dw_bracket=(-0.05, 0.05))
model = Base.TCPAMaxwellBlockedModel(disc=disc, integration="exponential")

# ------------------------- 2. relaxation d'une branche seule (code reel)
# On impose un etat de contrainte de branche puis on relaxe a deformation nulle.
j, k = 0, 0
sig0 = np.array([1.0, 0.3, -0.2, 0.0, 0.0, 0.5])
for ib in range(model.n_maxwell):
    model.sigma_i[ib, j, k] = sig0
de0 = np.zeros(6)
dt = 3.7
n_steps = 40
sim = model.sigma_i[:, j, k].copy()
for _ in range(n_steps):
    # _update_maxwell_branches lit self.sigma_i : on le met a jour a la main
    model.sigma_i[:, j, k] = model._update_maxwell_branches(j, k, de0, dt)
t_tot = n_steps * dt
for ib in range(model.n_maxwell):
    analytic = sig0 * np.exp(-t_tot / tau[ib])
    err = np.max(np.abs(model.sigma_i[ib, j, k] - analytic))
    check(f"relaxation exacte branche {ib+1} (t={t_tot:g}s)", err < 1e-12, f"err={err:.3e}")

# ------------------------- 3. rampe a taux constant, une branche scalaire
# d sigma/dt = E deps/dt - sigma/tau ; eps(t) = r t -> sigma(t) = E r tau (1 - exp(-t/tau))
E1, tau1 = 20.67, tau[0]
r_eps = 1.0e-3
def exp_update(s, deps, dt_, E_, tau_):
    a = np.exp(-dt_ / tau_)
    f = tau_ / dt_ * (1.0 - a) if dt_ > 0 else 1.0
    return a * s + f * E_ * deps

def euler_update(s, deps, dt_, E_, tau_):
    return s + E_ * deps - dt_ * s / tau_

T = 30.0
for n in (10, 100):
    dt_ = T / n
    s_exp, s_eul = 0.0, 0.0
    for _ in range(n):
        s_exp = exp_update(s_exp, r_eps * dt_, dt_, E1, tau1)
        s_eul = euler_update(s_eul, r_eps * dt_, dt_, E1, tau1)
    s_ana = E1 * r_eps * tau1 * (1.0 - np.exp(-T / tau1))
    print(f"  rampe n={n:4d}: exp={s_exp:.10f}  euler={s_eul:.10f}  analytique={s_ana:.10f}")
    check(f"exponentiel exact aux noeuds (n={n})", abs(s_exp - s_ana) < 1e-12,
          f"err={abs(s_exp-s_ana):.2e}")
# Euler converge en O(dt)
errs = []
for n in (10, 20, 40, 80):
    dt_ = T / n
    s_eul = 0.0
    for _ in range(n):
        s_eul = euler_update(s_eul, r_eps * dt_, dt_, E1, tau1)
    errs.append(abs(s_eul - s_ana))
ratios = [errs[i] / errs[i+1] for i in range(3)]
check("Euler ordre 1 (err/2 quand dt/2)", all(1.7 < r < 2.3 for r in ratios),
      f"ratios={['%.2f' % r for r in ratios]}")

# ------------------------- 4. coherence tangente/historique vs update
model2 = Base.TCPAMaxwellBlockedModel(disc=disc, integration="exponential")
rng = np.random.default_rng(0)
model2.sigma_i = rng.normal(size=model2.sigma_i.shape)
de = rng.normal(scale=1e-3, size=6)
for dt_ in (0.0, 0.25, 5.0, 500.0):
    tangents, history = model2._algorithmic_data(dt_)
    jj, kk = 1, 2
    new_branches = model2._update_maxwell_branches(jj, kk, de, dt_)
    total_pred = tangents[jj] @ de + history[jj, kk]
    total_real = model2.C0[jj] @ de + (new_branches.sum(axis=0) - model2.sigma_i[:, jj, kk].sum(axis=0))
    err = np.max(np.abs(total_pred - total_real))
    check(f"tangente+historique == somme branches (dt={dt_:g})", err < 1e-11, f"err={err:.3e}")

# meme test pour paper_explicit
model3 = Base.TCPAMaxwellBlockedModel(disc=disc, integration="paper_explicit")
model3.sigma_i = rng.normal(size=model3.sigma_i.shape)
for dt_ in (0.0, 0.25, 5.0):
    tangents, history = model3._algorithmic_data(dt_)
    jj, kk = 0, 1
    new_branches = model3._update_maxwell_branches(jj, kk, de, dt_)
    total_pred = tangents[jj] @ de + history[jj, kk]
    total_real = model3.C0[jj] @ de + (new_branches.sum(axis=0) - model3.sigma_i[:, jj, kk].sum(axis=0))
    err = np.max(np.abs(total_pred - total_real))
    check(f"paper_explicit coherent (dt={dt_:g})", err < 1e-11, f"err={err:.3e}")
    # forme Eq.(1) de l'article : dsigma_i = Ci de - dt (Ei/eta_i) sigma_i
    for ib in range(3):
        ds_paper = model3.Ci[jj][ib] @ de - dt_ * (E[ib]/eta[ib]) * model3.sigma_i[ib, jj, kk]
        err2 = np.max(np.abs(new_branches[ib] - model3.sigma_i[ib, jj, kk] - ds_paper))
        check(f"  Eq.(1) article, branche {ib+1} (dt={dt_:g})", err2 < 1e-13, f"err={err2:.3e}")

# ------------------------- 5. pas de double comptage
for jj in range(disc.n_layers):
    Csum = model2.C0[jj] + sum(model2.Ci[jj])
    err = np.max(np.abs(Csum - model2.C_total[jj]))
    check(f"C0+somme(Ci) == C_total couche {jj}", err < 1e-10, f"err={err:.3e}")
# fractions
Etot = mw.E0 + E.sum()
check("E_axial defaut == somme Maxwell 37.76", abs(model2.mat.E_axial - 37.76) < 1e-12
      and abs(Etot - 37.76) < 1e-12)
err = np.max(np.abs(model2.C0[0] - (mw.E0/Etot) * model2.C_total[0]))
check("C0 = (E0/Etot) C_E (hyp. ii)", err < 1e-10, f"err={err:.3e}")

# tangente a dt=0 == C_total (reponse instantanee)
tangents0, hist0 = model2._algorithmic_data(0.0)
err = max(np.max(np.abs(tangents0[jj] - model2.C_total[jj])) for jj in range(disc.n_layers))
check("tangente(dt=0) == C_E instantane", err < 1e-10, f"err={err:.3e}")
check("historique(dt=0) == 0", float(np.max(np.abs(hist0))) == 0.0)

# grande dt : tangente -> C0 (branche d'equilibre seule)
tangents_inf, _ = model2._algorithmic_data(1.0e9)
err = np.max(np.abs(tangents_inf[0] - model2.C0[0]))
check("tangente(dt->inf) -> C0 (equilibre)", err < 1e-6, f"err={err:.3e}")

# ------------------------- 6. garde de stabilite Euler
tau_min = tau.min()
try:
    model3._validate_time_step(2.0 * tau_min + 1e-9)
    check("garde Euler dt>=2*tau_min leve", False)
except ValueError:
    check("garde Euler dt>=2*tau_min leve", True, f"(2*tau_min={2*tau_min:.4f}s)")
model3._validate_time_step(2.0 * tau_min - 1e-6)   # ne doit pas lever
check("garde Euler dt<2*tau_min passe", True)
# l'exponentiel n'est pas limite
model2._validate_time_step(1.0e6)
check("exponentiel sans limite de dt", True)

# ------------------------- 7. dt = 0 : reponse elastique instantanee
model4 = Base.TCPAMaxwellBlockedModel(disc=disc, integration="exponential")
new_b = model4._update_maxwell_branches(0, 0, de, 0.0)
for ib in range(3):
    err = np.max(np.abs(new_b[ib] - model4.Ci[0][ib] @ de))
    check(f"dt=0 branche {ib+1}: dsigma = Ci de", err < 1e-14, f"err={err:.3e}")

print()
print("RESULTAT GLOBAL :", "TOUT OK" if ok else "AU MOINS UN ECHEC")
