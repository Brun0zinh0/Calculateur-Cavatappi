# -*- coding: utf-8 -*-
"""Maintien sous pression D 20 min vs modele Base.py.

Deux runs modele :
  A) run_hold_relaxation a P constant = P_exp(0) -> relaxation viscoelastique pure
  B) run_blocked_actuation pilote par la pression MESUREE (qui decroit de 6.39 a 4.03 bar)
Comparaison en force normalisee F(t)/F(debut de maintien).
Fit multi-exponentiel de la decroissance experimentale -> constantes de temps.
"""
import sys, time
sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import Base

BASE = r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/expérimentale/Sacha"
OUT = r"C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad"

# ---------------- donnees experimentales ----------------
df = pd.read_excel(BASE + "/Maintient sous pression D 20 mn.xlsx", header=0)
P_exp = df.iloc[:, 0].astype(float).values      # bar
F_exp = df.iloc[:, 1].astype(float).values      # N
dt_exp = 0.2
t_exp = np.arange(len(P_exp)) * dt_exp

# ---------------- fit exponentiel de l experience ----------------
# F(t) = A0 + A1 exp(-t/tau1) + A2 exp(-t/tau2)
def biexp(t, A0, A1, tau1, A2, tau2):
    return A0 + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

p0 = [1.0, 0.1, 100.0, 0.15, 1000.0]
popt, _ = curve_fit(biexp, t_exp, F_exp, p0=p0, maxfev=50000)
F_fit = biexp(t_exp, *popt)
rmse_fit = float(np.sqrt(np.mean((F_fit - F_exp) ** 2)))
taus_exp = sorted([popt[2], popt[4]])
print(f"Fit bi-exponentiel exp: A0={popt[0]:.3f} N, A1={popt[1]:.3f} N (tau={popt[2]:.0f} s), "
      f"A2={popt[3]:.3f} N (tau={popt[4]:.0f} s), RMSE={rmse_fit*1000:.1f} mN")

# meme fit sur la pression (la pression decroit aussi !)
poptP, _ = curve_fit(biexp, t_exp, P_exp, p0=[4.0, 1.0, 100.0, 1.5, 1000.0], maxfev=50000)
print(f"Fit bi-exponentiel pression: P_inf={poptP[0]:.2f} bar, taus=({poptP[2]:.0f}, {poptP[4]:.0f}) s")

# correlation F vs P
A = np.vstack([P_exp, np.ones_like(P_exp)]).T
coef, res, *_ = np.linalg.lstsq(A, F_exp, rcond=None)
F_lin = A @ coef
r2 = 1 - np.sum((F_exp - F_lin) ** 2) / np.sum((F_exp - F_exp.mean()) ** 2)
print(f"Regression F = a*P + b : a={coef[0]:.4f} N/bar, b={coef[1]:.3f} N, R2={r2:.4f}")

# ---------------- modele A : relaxation pure a P constant ----------------
P0_MPa = P_exp[0] / 10.0
t0 = time.time()
modelA, arrA = Base.run_hold_relaxation(P_hold=P0_MPa, hold_time=1500.0, ramp_time=9.0, dt=2.0)
print(f"run A (P constant {P0_MPa:.3f} MPa) : {time.time()-t0:.0f} s, {len(arrA['time'])} pas")

# ---------------- modele B : pression mesuree ----------------
stride = 10                                   # dt modele = 2.0 s
t_mod = t_exp[::stride].copy()
P_mod = P_exp[::stride] / 10.0                # MPa
if t_mod[0] == 0.0:
    t_mod[0] = 0.0
t0 = time.time()
modelB, arrB = Base.run_blocked_actuation(pressure_time=t_mod, pressure_MPa=P_mod)
print(f"run B (pression mesuree) : {time.time()-t0:.0f} s, {len(arrB['time'])} pas")

# ---------------- normalisation et metriques ----------------
# reference = premier point du maintien pour l exp,
# pour A : fin de rampe ; pour B : t=0 (pression appliquee en echelon)
i0A = int(arrA["hold_start_index"])
tA = arrA["time"] - arrA["time"][i0A]
FA = arrA["force_total_mN"] / 1000.0
FA_n = FA / FA[i0A]

tB = arrB["time"]
FB = arrB["force_total_mN"] / 1000.0
FB_n = FB / FB[0]

Fexp_n = F_exp / F_exp[0]

def decay_at(t_arr, y, t_q):
    i = np.searchsorted(t_arr, t_q)
    i = min(i, len(y) - 1)
    return 100.0 * (1.0 - y[i])

print("\nDecroissance relative de la force (en % de F au debut du maintien)")
print(f"{'t (s)':>8} {'exp':>8} {'mod A (P cst)':>14} {'mod B (P mesuree)':>18}")
for tq in [30, 60, 120, 300, 600, 900, 1200, 1500]:
    print(f"{tq:8d} {decay_at(t_exp, Fexp_n, tq):8.2f} {decay_at(tA[i0A:], FA_n[i0A:], tq):14.2f} {decay_at(tB, FB_n, tq):18.2f}")

# RMSE des courbes normalisees (interp modele B sur t_exp)
FB_i = np.interp(t_exp, tB, FB_n)
rmseB = float(np.sqrt(np.mean((FB_i - Fexp_n) ** 2)))
FA_i = np.interp(t_exp, tA[i0A:], FA_n[i0A:])
rmseA = float(np.sqrt(np.mean((FA_i - Fexp_n) ** 2)))
print(f"\nRMSE force normalisee: modele A vs exp = {100*rmseA:.2f} %, modele B vs exp = {100*rmseB:.2f} %")
print(f"Decroissance totale a 1500 s : exp {decay_at(t_exp, Fexp_n, 1500):.1f} % ; "
      f"mod A {decay_at(tA[i0A:], FA_n[i0A:], 1500):.1f} % ; mod B {decay_at(tB, FB_n, 1500):.1f} %")
print(f"Chute de pression exp sur le maintien : {100*(1-P_exp[-1]/P_exp[0]):.1f} %")
print(f"Force totale absolue: exp debut {F_exp[0]:.3f} N ; modele (eps=0.8 defaut) {FB[0]:.3f} N")

# ---------------- figure ----------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

ax = axes[0, 0]
ax.plot(t_exp, F_exp, "k-", lw=1, label="exp F (N)")
ax.plot(t_exp, F_fit, "g--", lw=1, label=f"fit 2 exp (tau={taus_exp[0]:.0f}, {taus_exp[1]:.0f} s)")
ax.set_xlabel("t (s)"); ax.set_ylabel("F (N)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("Maintien D : force mesuree", fontsize=10)
ax2 = ax.twinx(); ax2.plot(t_exp, P_exp, "r-", lw=0.8, alpha=0.6); ax2.set_ylabel("P (bar)", color="r")

ax = axes[0, 1]
ax.plot(P_exp, F_exp, "b.", ms=1, label="exp")
ax.plot(P_exp, F_lin, "r-", lw=1, label=f"lineaire a={coef[0]:.3f} N/bar, R2={r2:.3f}")
ax.set_xlabel("P (bar)"); ax.set_ylabel("F (N)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("F vs P pendant le maintien (la pression derive)", fontsize=10)

ax = axes[1, 0]
ax.plot(t_exp, 100 * Fexp_n, "k-", lw=1.2, label="exp")
ax.plot(tA[i0A:], 100 * FA_n[i0A:], "b--", lw=1.2, label="modele A : P constant 0.639 MPa")
ax.plot(tB, 100 * FB_n, "r-.", lw=1.2, label="modele B : P(t) mesuree")
ax.set_xlabel("t depuis debut du maintien (s)"); ax.set_ylabel("F / F0  (%)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title("Force normalisee : exp vs modele (defauts article)", fontsize=10)

ax = axes[1, 1]
ax.semilogx(t_exp[1:], 100 * (1 - Fexp_n[1:]), "k-", lw=1.2, label="exp")
ax.semilogx(tA[i0A + 1:] - tA[i0A], 100 * (1 - FA_n[i0A + 1:]), "b--", lw=1.2, label="modele A")
ax.semilogx(tB[1:], 100 * (1 - FB_n[1:]), "r-.", lw=1.2, label="modele B")
for tau in [7.5, 163.5, 2325.2]:
    ax.axvline(tau, color="gray", ls=":", lw=0.8)
ax.text(7.5, 1, "tau1", fontsize=7); ax.text(163.5, 1, "tau2", fontsize=7); ax.text(2325, 1, "tau3", fontsize=7)
ax.set_xlabel("t (s, log)"); ax.set_ylabel("chute de force (%)")
ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
ax.set_title("Chute relative (echelle log) et taus de Maxwell du modele", fontsize=10)

fig.suptitle("Maintien sous pression muscle D (20 min) vs modele bloque Maxwell (parametres article)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT + "/fig2_maintien_vs_modele.png", dpi=130)
print("saved fig2_maintien_vs_modele.png")
