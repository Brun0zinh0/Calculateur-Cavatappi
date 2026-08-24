# -*- coding: utf-8 -*-
"""Mises sous pression avec 200 g (muscles B, I, J) vs modele Base.py.

Signal mesure = force (N). Deux confrontations :
  1) actionnement bloque avec precharge : run_blocked_actuation pilote par P(t) mesuree,
     comparaison des increments de force Delta F(P).
  2) actionnement suspendu (masse 200 g) : run_suspended_actuation -> course predite
     (pas de deplacement mesure : prediction seulement).
Cadence d echantillonnage inconnue pour B/I/J -> hypothese dt = 0.2 s (comme le maintien).
"""
import sys, time
sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import Base

BASE = r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/expérimentale/Sacha"
OUT = r"C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad"
DT = 0.2
POIDS = 0.200 * 9.81  # 1.962 N

def load(name):
    if name.endswith(".csv"):
        df = pd.read_csv(BASE + "/" + name, sep=";", decimal=",", usecols=[0, 1])
    else:
        df = pd.read_excel(BASE + "/" + name, header=0)
    P = df.iloc[:, 0].astype(float).values
    F = df.iloc[:, 1].astype(float).values
    return P, F

muscles = {
    "B": "Mise sous pression muscle B 200 g.xlsx",
    "I": "Mise sous pression muscle I 200 g.csv",
    "J": "Mise sous pression muscle J 200 g.csv",
}

results = {}
for key, fname in muscles.items():
    P, F = load(fname)
    t = np.arange(len(P)) * DT
    # retire la queue ou la masse est posee / relachee (F < 0.5 N)
    valid = F > 0.5
    last = len(F) if valid.all() else np.argmax(~valid)
    t_u, P_u, F_u = t[:last], P[:last], F[:last]
    Fb = float(np.median(F_u[:10]))
    dF_exp = float(F_u.max() - Fb)
    # modele bloque pilote par la pression mesuree
    t_hist = t_u.copy(); t_hist[0] = 0.0
    tt0 = time.time()
    _, arrB = Base.run_blocked_actuation(pressure_time=t_hist, pressure_MPa=P_u / 10.0, dt=DT)
    dF_mod = float(arrB["force_act_mN"].max()) / 1000.0
    results[key] = dict(t=t, P=P, F=F, t_u=t_u, Fb=Fb, dF_exp=dF_exp,
                        arrB=arrB, dF_mod=dF_mod)
    print(f"Muscle {key}: baseline={Fb:.3f} N (poids {POIDS:.3f}), Pmax={P_u.max():.2f} bar, "
          f"dF_exp={dF_exp:.3f} N, dF_modele_bloque={dF_mod:.3f} N "
          f"(ratio exp/mod={dF_exp/dF_mod:.1f}) [{time.time()-tt0:.0f}s]")

# ------- run suspendu 200 g pour J (prediction de course) -------
P, F = load(muscles["J"])
t = np.arange(len(P)) * DT
valid = F > 0.5
last = len(F) if valid.all() else np.argmax(~valid)
t_hist = t[:last].copy(); t_hist[0] = 0.0
tt0 = time.time()
model_s, arrS = Base.run_suspended_actuation(load_N=POIDS, pressure_time=t_hist,
                                             pressure_MPa=P[:last] / 10.0, dt=DT)
print(f"\nRun suspendu (J, 200 g) : {time.time()-tt0:.0f} s")
print(f"  longueur de reference = {arrS['reference_axial_length_mm'][0]:.2f} mm")
print(f"  contraction max predite = {np.max(arrS['free_contraction_mm']):.3f} mm "
      f"({np.max(arrS['free_actuation_percent']):.2f} %)")
print(f"  contraction finale = {arrS['free_contraction_mm'][-1]:.3f} mm")

# ------- figure -------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
colors = {"B": "tab:blue", "I": "tab:green", "J": "tab:red"}

ax = axes[0, 0]
for key, r in results.items():
    ax.plot(r["t"], r["F"], color=colors[key], lw=1, label=f"{key} : F exp")
ax.axhline(POIDS, color="k", ls=":", lw=1, label="poids 200 g = 1.962 N")
ax.set_xlabel("t (s, dt suppose 0.2 s)"); ax.set_ylabel("F (N)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("Force mesuree : plateau = poids de la masse (test isotone)", fontsize=10)

ax = axes[0, 1]
for key, r in results.items():
    ax.plot(r["P"][:len(r['t_u'])], r["F"][:len(r['t_u'])] - r["Fb"], color=colors[key], lw=1,
            label=f"{key} exp (dF)")
arr0 = results["J"]["arrB"]
ax.plot(10 * arr0["pressure_MPa"], arr0["force_act_mN"] / 1000.0, "k--", lw=1.5,
        label="modele bloque (defauts)")
ax.set_xlabel("P (bar)"); ax.set_ylabel("Delta F (N)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("Increment de force vs pression : exp vs modele bloque", fontsize=10)

ax = axes[1, 0]
for key, r in results.items():
    tu = r["t_u"]
    ax.plot(tu, (r["F"][:len(tu)] - r["Fb"]) * 1000, color=colors[key], lw=1, label=f"{key} exp")
    ax.plot(r["arrB"]["time"], r["arrB"]["force_act_mN"], color=colors[key], ls="--", lw=1,
            label=f"{key} modele bloque")
ax.set_xlabel("t (s)"); ax.set_ylabel("Delta F (mN)")
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)
ax.set_title("Delta F(t) : exp vs modele bloque pilote par P(t) mesuree", fontsize=10)

ax = axes[1, 1]
ax.plot(arrS["time"], arrS["free_contraction_mm"], "r-", lw=1.4, label="contraction predite (mm)")
ax.set_xlabel("t (s)"); ax.set_ylabel("contraction (mm)", color="r")
ax2 = ax.twinx()
ax2.plot(arrS["time"], 10 * arrS["pressure_MPa"], "b:", lw=1, label="P (bar)")
ax2.set_ylabel("P (bar)", color="b")
ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=0.3)
ax.set_title("Modele suspendu 200 g (J) : course predite - NON mesuree", fontsize=10)

fig.suptitle("Mises sous pression 200 g (B, I, J) vs modele (parametres article par defaut)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT + "/fig3_mise_sous_pression_vs_modele.png", dpi=130)
print("saved fig3_mise_sous_pression_vs_modele.png")

# ------- hysteresis exp (J possede une descente complete) -------
print("\nHysteresis (muscle J, branche montee vs descente) :")
P, F = load(muscles["J"])
Fb = float(np.median(F[:10]))
imax = int(np.argmax(P))
for pq in [2.0, 4.0, 6.0]:
    iu = np.argmin(np.abs(P[:imax] - pq))
    idn = imax + np.argmin(np.abs(P[imax:imax+85] - pq))
    print(f"  P={pq:.0f} bar : dF montee = {F[iu]-Fb:+.3f} N ; descente = {F[idn]-Fb:+.3f} N")
