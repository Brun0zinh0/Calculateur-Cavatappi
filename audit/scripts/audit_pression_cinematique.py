# -*- coding: utf-8 -*-
"""Audit 'Transmission de la pression et cinematique' - verifications numeriques independantes.

1. Ecart angulaire entre le profil par defaut 'paper_linear' (theta_j = (R/Rout)*theta_f)
   et l'Eq. (11) de l'article BLOCKED (theta_j = arctan((R/Rout)*tan(theta_f))).
2. Impact sur la simulation bloquee par defaut (force/couple pics).
3. Verification du profil cyclique de pression (periode, amplitude, evenements).
4. Verification des conversions d'unites de pression.py.
5. Verification Eq.6/Eq.9 : conservation de rho*tan(alpha) et du nombre de spires.
"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

# ---------------------------------------------------------------- 1. angles
theta_f = np.deg2rad(37.91)
frac = np.linspace(0.4, 1.0, 7)  # Rin/Rout=0.4 -> 1.0
lin = frac * theta_f
tanlaw = np.arctan(frac * np.tan(theta_f))
print("=== 1. Profil d'angle de biais ===")
print("R/Rout   lineaire(deg)   Eq11(deg)   ecart(deg)")
for f, a, b in zip(frac, np.rad2deg(lin), np.rad2deg(tanlaw)):
    print(f"{f:5.2f}    {a:9.3f}    {b:9.3f}    {a-b:8.3f}")

# ---------------------------------------------------------------- 2. impact simulation
def run(profile):
    cfg = Base.default_simulation_config(
        geom=Base.default_geometry_params(bias_angle_profile=profile),
        n_cycles=2,
        dt=0.25,
    )
    model, arr = Base.run_blocked_actuation(cfg)
    return arr

print("\n=== 2. Impact sur la simulation bloquee (defauts, 2 cycles, P triangulaire 1.3 MPa) ===")
res = {}
for profile in ("paper_linear", "uniform_twist"):
    arr = run(profile)
    m = arr["time"] >= arr["time"][-1] / 2.0  # dernier cycle
    res[profile] = dict(
        Fpic=float(np.nanmax(arr["force_total_mN"][m])),
        Fmin=float(np.nanmin(arr["force_total_mN"][m])),
        Tpic=float(np.nanmax(np.abs(arr["torque_total_signed_microNm"][m]))),
        Fact=float(np.nanmax(arr["force_act_mN"])),
        Tact=float(np.nanmax(np.abs(arr["torque_act_microNm"]))),
    )
    print(profile, {k: round(v, 2) for k, v in res[profile].items()})
for k in res["paper_linear"]:
    a, b = res["paper_linear"][k], res["uniform_twist"][k]
    if abs(b) > 1e-12:
        print(f"  ecart relatif {k}: {100.0*(a-b)/abs(b):+.2f} %")

# ---------------------------------------------------------------- 3. profil cyclique
print("\n=== 3. Profil cyclique ===")
t, P = Base.cyclic_pressure_history(n_cycles=3, Pmax=1.4, flow_rate_mL_min=10.0, volume_mL=1.50, dt=0.25, nonlinear=False)
print("demi-periode attendue 60*1.5/10 = 9 s ; periode 18 s ; total attendu 54 s")
print("t final:", t[-1], "  P final:", P[-1], "  Pmax atteint:", P.max())
ipks = np.where((P[1:-1] >= P[:-2]) & (P[1:-1] >= P[2:]) & (P[1:-1] > 1.3))[0] + 1
print("instants des pics:", t[ipks])
# pente lineaire attendue Pmax/9 :
i9 = np.argmin(np.abs(t - 4.5))
print("P(4.5 s) =", P[i9], " attendu", 1.4 * 4.5 / 9.0)

t2, P2 = Base.cyclic_pressure_history(n_cycles=1, Pmax=1.4, flow_rate_mL_min=10.0, volume_mL=1.22, dt=0.25, nonlinear=False)
print("volume 1.22 mL -> demi-periode:", 60.0 * 1.22 / 10.0, "s ; t_pic =", t2[np.argmax(P2)])

# ---------------------------------------------------------------- 4. unites
print("\n=== 4. Conversions pression.py ===")
psi_exact = 4.4482216152605 / (0.0254**2) / 1e6  # lbf/in2 en MPa
print("1 psi =", psi_exact, "MPa ; code: 0.006894757293168361 ; ecart:", psi_exact - 0.006894757293168361)
print("1 bar = 0.1 MPa (code 0.1) ; 1 kPa = 0.001 MPa (code 0.001)")
print("poids de 1 g =", 1e-3 * 9.80665 * 1e3, "mN (code 9.80665)")

# ---------------------------------------------------------------- 5. cinematique helice
print("\n=== 5. Cinematique : conservation de rho tan(alpha) et Eq.6 ===")
cfg = Base.default_simulation_config(n_cycles=1)
model, arr = Base.run_blocked_actuation(cfg)
rho = arr["rho_mm"]
alpha = np.deg2rad(arr["alpha_deg"])
h = rho * np.tan(alpha)
print("h = rho*tan(alpha) : min", h.min(), " max", h.max(), " h_blocked", model.h_blocked)
print("ecart max |h - h_blocked| =", np.max(np.abs(h - model.h_blocked)))
# Eq. 6 : rho/cos(alpha) doit croitre comme prod(1+dw)
l_fibre = rho / np.cos(alpha)
dw = arr["dw"]
l_rec = l_fibre[0] * np.cumprod(1.0 + dw[1:] )
print("Eq.6 (l_new=(1+dw)l) ecart relatif max:", np.max(np.abs(l_rec - l_fibre[1:]) / l_fibre[1:]))
