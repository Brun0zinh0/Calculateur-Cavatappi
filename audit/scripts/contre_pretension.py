# -*- coding: utf-8 -*-
"""Contre-expertise : precontrainte elastique (schema du code) vs precontrainte
viscoelastique (schema de l'article) — effet sur les niveaux absolus de force
des 11 premiers cycles (conditions figure 7 : eps=0.8, 10 mL/min, 1.5 mL).

Aucun fichier du projet n'est modifie : la variante B pilote step() directement
(sans le drapeau _building_reference_state), ce qui reproduit le schema de
l'article ou la loi viscoelastique court des la phase d'elongation.
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

EPS = 0.8
NCYC = 11
PMAX = 1.4
DT = 0.5

# v4-16 : profil par vitesse de pression ; le couple debit/volume historique
# est passe explicitement pour conserver la demi-periode de 9 s de ce script.
cfg = Base.default_simulation_config(eps=EPS, n_cycles=NCYC, Pmax=PMAX, dt=DT,
                                     flow_rate_mL_min=10.0, volume_mL=1.5)

# ---------- estimation analytique (rampe a taux constant, duree T) ----------
mw = cfg.mat.maxwell
E = np.array([mw.E1, mw.E2, mw.E3])
eta = np.array([mw.eta1, mw.eta2, mw.eta3])
tau = eta / E
E0 = mw.E0
Etot = E0 + E.sum()
L0 = cfg.geom.initial_length
T_pre = 60.0 * EPS * L0 / 20.0  # 20 mm/min
Eeff = E * (tau / T_pre) * (1.0 - np.exp(-T_pre / tau))
print("=== Estimation analytique (rampe uniaxiale, T = %.1f s) ===" % T_pre)
print("tau_i [s]           :", np.round(tau, 2))
print("E_i instantane [MPa]:", E, " E0 =", E0, " somme =", Etot)
print("E_i effectif  [MPa] :", np.round(Eeff, 3))
ratio_tk = (E0 + Eeff.sum()) / Etot
print("Module effectif a t_k = %.2f MPa  -> ratio vs instantane = %.3f" % (E0 + Eeff.sum(), ratio_tk))
for t_extra in (60.0, 90.0, 180.0, 1980.0):
    rel = Eeff * np.exp(-t_extra / tau)
    print("  apres +%6.0f s de relaxation : module tube = %.2f MPa (ratio %.3f ; equilibre = %.3f)"
          % (t_extra, E0 + rel.sum(), (E0 + rel.sum()) / Etot, E0 / Etot))

# ---------- variante A : schema du code (reference elastique t_k) ----------
model_a, arr_a = Base.run_blocked_actuation(config=cfg)
print("\n=== Variante A : reference elastique (code, run_blocked_actuation) ===")
print("max|sigma_i| apres precontrainte-plus-pas0 :", float(np.max(np.abs(model_a.sigma_i))) if model_a.sigma_i.size else 0.0)

# ---------- variante B : precontrainte viscoelastique (schema article) ----------
disc = Base.default_discretization(n_layers=cfg.n_layers, n_phi=cfg.n_phi,
                                   pre_steps=cfg.pre_steps, dw_bracket=(-0.05, 0.05))
m = Base.TCPAMaxwellBlockedModel(mat=cfg.mat, geom=cfg.geom, disc=disc,
                                 integration=cfg.integration)
h_end = (1.0 + EPS) * m.h0
Lact0 = 2.0 * np.pi * m.turns * m.h0
total_time = 60.0 * EPS * Lact0 / 20.0
dtp = total_time / disc.pre_steps
for h in np.linspace(m.h0, h_end, disc.pre_steps + 1)[1:]:
    m.step(0.0, dtp, h_target=h)   # branches ACTIVES : loi viscoelastique pendant l'elongation
m.h_blocked = h_end
sig_i_tk = float(np.max(np.abs(m.sigma_i)))
F_tk_B = float(m.history[-1].Ft)
t_start = m.helix.time
t_local, pressure = Base._prepare_actuation_history(cfg.n_cycles, cfg.Pmax, cfg.dt, None, None,
                                                    Base._config_half_period(cfg),
                                                    cfg.nonlinear_pressure)
m.step(float(pressure[0]), 0.0, h_target=m.h_blocked)
m.lock_blocked_series_reference()
i0 = len(m.history) - 1
m.run_pressure_history(t_start + t_local, pressure)
full = m.history_arrays()
arr_b = {k: v[i0:].copy() for k, v in full.items()}
arr_b["time"] = arr_b["time"] - t_start
Base.add_corrected_output_conventions(arr_b)
print("\n=== Variante B : precontrainte viscoelastique (schema article, branches actives) ===")
print("max|sigma_i| a t_k = %.4f MPa ; force a t_k = %.1f mN" % (sig_i_tk, 1000 * F_tk_B))

# ---------- comparaison ----------
period = 2.0 * Base._config_half_period(cfg)  # 18 s


def cycle_stats(arr):
    t = np.asarray(arr["time"])
    F = np.asarray(arr["force_mN"])
    Ftube = np.asarray(arr.get("force_tube_mN", np.full_like(F, np.nan)))
    Fny = np.asarray(arr.get("force_nylon_mN", np.full_like(F, np.nan)))
    rows = []
    for c in range(NCYC):
        m_ = (t >= c * period - 1e-9) & (t <= (c + 1) * period + 1e-9)
        if not m_.any():
            continue
        rows.append((c + 1, F[m_].max(), F[m_].min()))
    return rows, F, Ftube, Fny


rows_a, Fa, Fta, Fna = cycle_stats(arr_a)
rows_b, Fb, Ftb, Fnb = cycle_stats(arr_b)
print("\nDecomposition a t_k (1er echantillon) :")
print("  A : F = %7.1f mN (tube %7.1f, nylon %7.1f)" % (Fa[0], Fta[0], Fna[0]))
print("  B : F = %7.1f mN (tube %7.1f, nylon %7.1f)" % (Fb[0], Ftb[0], Fnb[0]))

print("\ncycle |   pic A   |   pic B   | B-A pic |  vallee A |  vallee B | B-A vallee")
for (ca, pa, va), (cb, pb, vb) in zip(rows_a, rows_b):
    print("  %2d  | %8.1f | %8.1f | %+7.1f | %8.1f | %8.1f | %+7.1f"
          % (ca, pa, pb, pb - pa, va, vb, vb - va))

da = rows_a[0][2] - rows_a[-1][2]
db = rows_b[0][2] - rows_b[-1][2]
print("\nDerive des vallees cycle 1 -> %d : A = %.1f mN ; B = %.1f mN" % (NCYC, da, db))
dpa = rows_a[0][1] - rows_a[-1][1]
dpb = rows_b[0][1] - rows_b[-1][1]
print("Derive des pics    cycle 1 -> %d : A = %.1f mN ; B = %.1f mN" % (NCYC, dpa, dpb))
print("\nResidu max A = %.2e ; B = %.2e N.mm"
      % (max(abs(h.residual) for h in model_a.history), max(abs(h.residual) for h in m.history)))
print("\nRappel article, fig. 7(a) eps=0.8 (11e cycle) : F crete 1824.6 mN, vallees exp ~1300-1400 mN.")
print("Rappel article, fig. 3 (eps=0.5) : F crete 1041.7 -> 934.4 mN sur les 10 cycles de training.")
