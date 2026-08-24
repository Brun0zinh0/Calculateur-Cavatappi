"""Audit numerique 2 : robustesse extremites tres compliantes + alpha0 vs alpha_tk."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

print("=" * 70)
print("A. Extremites tres compliantes (uncoiled_length = 40 mm), 2 cycles")
print("=" * 70)
geom = Base.default_geometry_params(uncoiled_length=40.0)
try:
    m, arr = Base.run_blocked_actuation(
        eps=0.8, n_cycles=2, Pmax=1.3, dt=1.0, n_layers=4, n_phi=16, pre_steps=8, geom=geom
    )
    print("OK. compliance (mm/N):", m.uncoiled_compliance_mm_per_N)
    print("max |series residual| mm:", np.max(np.abs(arr["series_compatibility_residual_mm"])))
    print("max |moment residual| Nmm:", np.nanmax(np.abs(arr["residual"])))
    print("force pic (mN):", arr["force_mN"].max(), " baseline:", arr["force_mN"][0])
    print("h min/max:", arr["h_mm_per_rad"].min(), arr["h_mm_per_rad"].max(), " ref:", m.series_reference_h)
    print("extension extremites max (mm):", arr["uncoiled_extension_mm"].max())
except Exception as e:
    print("ECHEC:", type(e).__name__, e)

print()
print("=" * 70)
print("B. alpha0 vs alpha_tk dans la compliance tangent_beam")
print("=" * 70)
g = Base.default_geometry_params(uncoiled_length=10.0)
m0 = Base.TCPAMaxwellBlockedModel(geom=g)
a0 = m0.alpha0
mat = m0.mat
EA = mat.E_axial * np.pi * (g.Rout**2 - g.Rin**2) + mat.E_nylon * np.pi * g.r_nylon**2
EI = mat.E_axial * 0.25 * np.pi * (g.Rout**4 - g.Rin**4) + mat.E_nylon * 0.25 * np.pi * g.r_nylon**4
Lh = 0.5 * g.uncoiled_length
def comp(alpha):
    return 2.0 * (Lh * np.sin(alpha)**2 / EA + Lh**3 * np.cos(alpha)**2 / (3.0 * EI))
# alpha_tk apres eps=0.8
m0.prestretch_to(0.8)
a_tk = m0.helix.alpha
print("alpha0 (deg):", np.rad2deg(a0), " alpha_tk (deg):", np.rad2deg(a_tk))
print("compliance(alpha0):", comp(a0), " compliance(alpha_tk):", comp(a_tk),
      " ecart rel:", (comp(a0) - comp(a_tk)) / comp(a0))

print()
print("=" * 70)
print("C. Chemin de repli de _find_dw (pas de changement de signe force)")
print("=" * 70)
# residu strictement positif sur tout le bracket ? on teste avec un bracket minuscule decale
disc = Base.default_discretization(n_layers=2, n_phi=8, pre_steps=4, dw_bracket=(0.02, 0.05))
m2 = Base.TCPAMaxwellBlockedModel(disc=disc)
m2.prestretch_to(0.5)
try:
    dw = m2._find_dw(0.3, 1.0, m2.h_blocked)
    print("dw trouve (repli):", dw)
except RuntimeError as e:
    print("RuntimeError attendue (solution en butee / residu trop grand):", e)

print()
print("=" * 70)
print("D. dt = 0 sur le pas initial : reponse elastique instantanee")
print("=" * 70)
disc = Base.default_discretization(n_layers=2, n_phi=8, pre_steps=4, dw_bracket=(-0.05, 0.05))
m3 = Base.TCPAMaxwellBlockedModel(disc=disc)
m3.prestretch_to(0.5)
r0 = m3.step(0.0, 0.0, h_target=m3.h_blocked)
print("residu pas dt=0 a P=0:", r0.residual, " dw:", r0.dw)
r1 = m3.step(0.5, 0.0, h_target=m3.h_blocked)
print("saut instantane a 0.5 MPa: dw:", r1.dw, " Ft (N):", r1.Ft, " residu:", r1.residual)
print("sigma_i nuls avant saut ? (base elastique) max|sigma_i| apres saut:",
      float(np.max(np.abs(m3.sigma_i))))
