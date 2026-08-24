"""Audit numerique independant : solveur bloque + extremites compliantes."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

np.set_printoptions(precision=6, suppress=True)

print("=" * 70)
print("1. Run bloque rigide (uncoiled_length = 0)")
print("=" * 70)
model, arr = Base.run_blocked_actuation(
    eps=0.8, n_cycles=1, Pmax=1.3, dt=1.0, n_layers=4, n_phi=16, pre_steps=8
)
print("n steps:", len(arr["time"]))
print("max |residual| (N mm):", np.nanmax(np.abs(arr["residual"])))
print("force range (mN):", arr["force_mN"].min(), arr["force_mN"].max())
print("force[0] (baseline post-prestretch, mN):", arr["force_mN"][0])
print("torque_act range (uNm):", arr["torque_act_microNm"].min(), arr["torque_act_microNm"].max())
print("h varie ? min/max h_mm_per_rad:", arr["h_mm_per_rad"].min(), arr["h_mm_per_rad"].max())
print("h_blocked:", model.h_blocked, " h0:", model.h0)

print()
print("=" * 70)
print("2. Pre-tension : purement elastique ? (maintien P=0 pendant 600 s)")
print("=" * 70)
t_hold = np.linspace(0.0, 600.0, 61)
p_hold = np.zeros_like(t_hold)
m2, arr2 = Base.run_blocked_actuation(
    eps=0.8, n_layers=4, n_phi=16, pre_steps=8,
    pressure_time=t_hold, pressure_MPa=p_hold,
)
f = arr2["force_mN"]
print("force debut / fin (mN):", f[0], f[-1], " drift:", f[-1] - f[0])
print("=> si drift == 0 : aucune relaxation de la precontrainte (base elastique)")

print()
print("=" * 70)
print("3. Extremites compliantes serie (uncoiled_length = 10 mm, tangent_beam)")
print("=" * 70)
geom = Base.default_geometry_params(uncoiled_length=10.0)
m3, arr3 = Base.run_blocked_actuation(
    eps=0.8, n_cycles=1, Pmax=1.3, dt=1.0, n_layers=4, n_phi=16, pre_steps=8, geom=geom
)
print("compliance (mm/N):", m3.uncoiled_compliance_mm_per_N)
print("max |series_compatibility_residual| (mm):", np.max(np.abs(arr3["series_compatibility_residual_mm"])))
print("max |residual moment| (N mm):", np.nanmax(np.abs(arr3["residual"])))
# verification manuelle de la compatibilite : delta_active + delta_ext = 0
h_arr = arr3["h_mm_per_rad"]
dL_active = 2.0 * np.pi * m3.turns * (h_arr - m3.series_reference_h)
dL_ext = m3.uncoiled_compliance_mm_per_N * (arr3["force_N"] - m3.series_reference_force_N)
print("max |dL_active + dL_ext| (mm):", np.max(np.abs(dL_active + dL_ext)))
print("uncoiled_extension == compliance*dF ?",
      np.allclose(arr3["uncoiled_extension_mm"], dL_ext))
print("h bouge (respiration de la spire) ? min/max:", h_arr.min(), h_arr.max(), " ref:", m3.series_reference_h)
print("pic de force rigide vs serie (mN):", arr["force_mN"].max(), arr3["force_mN"].max())
print("force baseline rigide vs serie (mN):", arr["force_mN"][0], arr3["force_mN"][0])

print()
print("=" * 70)
print("4. Formule de compliance tangent_beam : re-derivation independante")
print("=" * 70)
g = m3.geom
mat = m3.mat
EA = mat.E_axial * np.pi * (g.Rout**2 - g.Rin**2) + mat.E_nylon * np.pi * g.r_nylon**2
EI = mat.E_axial * 0.25 * np.pi * (g.Rout**4 - g.Rin**4) + mat.E_nylon * 0.25 * np.pi * g.r_nylon**4
a0 = np.deg2rad(g.alpha0_deg)
Lh = 0.5 * g.uncoiled_length
c_hand = 2.0 * (Lh * np.sin(a0)**2 / EA + Lh**3 * np.cos(a0)**2 / (3.0 * EI))
print("EA (N):", EA, " EI (N mm^2):", EI)
print("compliance main:", c_hand, " modele:", m3.uncoiled_compliance_mm_per_N,
      " ecart rel:", abs(c_hand - m3.uncoiled_compliance_mm_per_N) / c_hand)
print("part flexion / part traction:",
      (Lh**3 * np.cos(a0)**2 / (3.0 * EI)) / (Lh * np.sin(a0)**2 / EA))
# mode axial_rod
geom_rod = Base.default_geometry_params(uncoiled_length=10.0, uncoiled_compliance_mode="axial_rod")
m_rod = Base.TCPAMaxwellBlockedModel(geom=geom_rod)
print("axial_rod: modele:", m_rod.uncoiled_compliance_mm_per_N, " main:", 10.0 / EA)

print()
print("=" * 70)
print("5. Multi-racines du residu de moment (balayage dw au 1er pas pressurise)")
print("=" * 70)
# reconstruire un modele frais et regarder le profil du residu
disc = Base.default_discretization(n_layers=4, n_phi=16, pre_steps=8, dw_bracket=(-0.05, 0.05))
m5 = Base.TCPAMaxwellBlockedModel(disc=disc)
m5.prestretch_to(0.8)
m5.step(0.0, 0.0, h_target=m5.h_blocked)
dP = 0.3
grid = np.linspace(-0.05, 0.05, 201)
vals = []
for x in grid:
    try:
        vals.append(float(m5._trial_state(x, dP, 1.0, m5.h_blocked)["residual"]))
    except Exception:
        vals.append(np.nan)
vals = np.array(vals)
finite = np.isfinite(vals)
sign_changes = np.sum(np.diff(np.sign(vals[finite])) != 0)
print("points finis:", finite.sum(), "/", len(grid))
print("changements de signe du residu:", sign_changes)
print("residu min/max:", np.nanmin(vals), np.nanmax(vals))
i_nan = np.where(~finite)[0]
print("zones NaN (indices):", i_nan[:5], "..." if len(i_nan) > 5 else "")

print()
print("=" * 70)
print("6. Historique de pression demarrant a P0 > 0 : verrou de reference serie")
print("=" * 70)
t6 = np.linspace(0.0, 60.0, 31)
p6 = np.full_like(t6, 0.5)
geom6 = Base.default_geometry_params(uncoiled_length=10.0)
m6, arr6 = Base.run_blocked_actuation(
    eps=0.8, n_layers=4, n_phi=16, pre_steps=8, geom=geom6,
    pressure_time=t6, pressure_MPa=p6,
)
print("force[0] au verrou (mN):", arr6["force_mN"][0])
print("uncoiled_extension[0] (mm):", arr6["uncoiled_extension_mm"][0])
print("=> la reference serie est verrouillee APRES le saut a P0=0.5 MPa")
print("   (l'allongement des extremites du au saut initial est compte comme nul)")

print()
print("=" * 70)
print("7. Prestretch et extremites compliantes : eps applique a la spire seule ?")
print("=" * 70)
print("h_end rigide:", model.h_blocked / model.h0, " serie:", m3.h_blocked / m3.h0,
      "(les deux = 1+eps => la precontrainte ignore la compliance des extremites)")
print("series_reference_force rigide vs serie (N):",
      model.series_reference_force_N, m3.series_reference_force_N)
