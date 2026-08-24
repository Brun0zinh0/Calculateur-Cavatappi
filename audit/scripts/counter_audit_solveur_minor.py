# -*- coding: utf-8 -*-
"""Contre-expertise des 5 constats 'minor' d'audit_solveur.json (hors profil de biais)."""
import sys
import numpy as np

sys.path.insert(0, r"C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2")
import Base

print("MODEL_VERSION =", Base.MODEL_VERSION)

# ---------------------------------------------------------------- Test 1
# Constat: la precontrainte ignore la compliance des extremites.
print("\n=== Test 1 : precontrainte vs uncoiled_length ===")
for unc in (0.0, 10.0, 40.0):
    geom = Base.default_geometry_params(uncoiled_length=unc)
    disc = Base.default_discretization(n_layers=4, n_phi=16, pre_steps=60, dw_bracket=(-0.05, 0.05))
    m = Base.TCPAMaxwellBlockedModel(geom=geom, disc=disc)
    m.prestretch_to(0.8)
    print(f"unc={unc:5.1f} mm  h_blocked/h0={m.h_blocked/m.h0:.6f}  "
          f"F_fin_prestretch={m.history[-1].Ft*1000:.4f} mN  "
          f"alpha_fin={np.rad2deg(m.helix.alpha):.4f} deg  "
          f"compliance={m.uncoiled_compliance_mm_per_N:.6f} mm/N")

# ---------------------------------------------------------------- Test 2
# Constat: compliance tangent_beam evaluee a alpha0 et figee.
print("\n=== Test 2 : compliance a alpha0 vs alpha_tk ===")
geom = Base.default_geometry_params(uncoiled_length=10.0)
disc = Base.default_discretization(n_layers=4, n_phi=16, pre_steps=60, dw_bracket=(-0.05, 0.05))
m = Base.TCPAMaxwellBlockedModel(geom=geom, disc=disc)
alpha0 = m.alpha0
m.prestretch_to(0.8)
alpha_tk = m.helix.alpha
EA, EI = m._uncoiled_section_rigidities()
L = m.uncoiled_length
def comp(alpha):
    hl = 0.5 * L
    return 2.0 * (hl * np.sin(alpha)**2 / EA + hl**3 * np.cos(alpha)**2 / (3.0 * EI))
c0 = comp(alpha0)
ctk = comp(alpha_tk)
print(f"alpha0={np.rad2deg(alpha0):.4f} deg, alpha_tk={np.rad2deg(alpha_tk):.4f} deg")
print(f"compliance modele  = {m.uncoiled_compliance_mm_per_N:.6f} mm/N")
print(f"compliance(alpha0) = {c0:.6f} mm/N (doit egaler le modele)")
print(f"compliance(alpha_tk)= {ctk:.6f} mm/N  ecart = {100*(c0-ctk)/c0:.2f} %")
ax = 0.5*L*np.sin(alpha0)**2/EA*2
fl = 0.5*L**3/4*np.cos(alpha0)**2/(3*EI)  # not exact, just ratio check below
print(f"ratio flexion/traction a alpha0 = "
      f"{( (0.5*L)**3*np.cos(alpha0)**2/(3*EI) ) / ( 0.5*L*np.sin(alpha0)**2/EA ):.1f}")

# ---------------------------------------------------------------- Test 3
# Constat: longueur axiale ajoute la longueur desenroulee complete.
print("\n=== Test 3 : axial_length et longueur des extremites ===")
model, arr = Base.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3,
                                        geometry_overrides=None) if False else (None, None)
# run via config
cfg = Base.default_simulation_config()
cfg.geom = Base.default_geometry_params(uncoiled_length=10.0)
cfg.n_cycles = 1
model, arr = Base.run_blocked_actuation(cfg)
i = len(arr["time"]) // 3
print(f"active_axial[{i}]={arr['active_axial_length_mm'][i]:.4f} mm, "
      f"uncoiled_deformed[{i}]={arr['uncoiled_deformed_length_mm'][i]:.4f} mm, "
      f"axial_length[{i}]={arr['axial_length_mm'][i]:.4f} mm")
print("axial = active + uncoiled complet ?",
      np.allclose(arr['axial_length_mm'], arr['active_axial_length_mm'] + arr['uncoiled_deformed_length_mm']))
print(f"si projection axiale sin(alpha_tk): contribution = "
      f"{10.0*np.sin(arr['alpha_rad'][0]):.3f} mm au lieu de ~10 mm")
print(f"max |series_compatibility_residual| = {np.abs(arr['series_compatibility_residual_mm']).max():.3e} mm")

# ---------------------------------------------------------------- Test 4
# Constat: fragilites de _find_dw. (a) brentq avec NaN interieur, (b) granularite.
print("\n=== Test 4 : brentq et NaN interieur ===")
from scipy.optimize import brentq
def f_nan(x):
    if 0.3 < x < 0.42:
        return np.nan
    return x - 0.5
try:
    r = brentq(f_nan, 0.0, 1.0, xtol=1e-10, rtol=1e-9, maxiter=100)
    print(f"brentq avec NaN interieur -> renvoie {r} sans exception ; f(r) = {f_nan(r)}")
except Exception as e:
    print("brentq avec NaN interieur -> exception:", type(e).__name__, e)

def f_nan2(x):
    if 0.49 < x < 0.51:
        return np.nan
    return x - 0.5
try:
    r = brentq(f_nan2, 0.0, 1.0, xtol=1e-10, rtol=1e-9, maxiter=100)
    print(f"brentq NaN autour de la racine -> renvoie {r} ; f(r) = {f_nan2(r)}")
except Exception as e:
    print("brentq NaN autour de la racine -> exception:", type(e).__name__, e)

lo, hi = -0.05, 0.05
grid = np.linspace(lo, hi, 65)
print(f"pas du balayage = {grid[1]-grid[0]:.6e} (attendu 1.5625e-3)")

# verifier que step() n'a pas de garde de residu apres _find_dw :
import inspect
src = inspect.getsource(Base.TCPAMaxwellBlockedModel.step)
print("step() contient un raise sur residual ?", "residual" in src and "raise" in src)

# ---------------------------------------------------------------- Test 5
# Constat: verrou serie apres le pas initial a P(0) -> extension nulle si P0>0.
print("\n=== Test 5 : historique demarrant a P0 = 0.5 MPa ===")
cfg2 = Base.default_simulation_config()
cfg2.geom = Base.default_geometry_params(uncoiled_length=10.0)
t = np.linspace(0.0, 10.0, 21)
p = np.full_like(t, 0.5)
model2, arr2 = Base.run_blocked_actuation(cfg2, pressure_time=t, pressure_MPa=p)
print(f"force[0] = {arr2['force_mN'][0]:.2f} mN (verrou), "
      f"uncoiled_extension[0] = {arr2['uncoiled_extension_mm'][0]:.4e} mm")
print(f"force[-1] = {arr2['force_mN'][-1]:.2f} mN, "
      f"uncoiled_extension[-1] = {arr2['uncoiled_extension_mm'][-1]:.4e} mm")
# comparaison : meme historique mais precede d'un point a P=0
t3 = np.concatenate([[0.0], t + 1.0])
p3 = np.concatenate([[0.0], p])
model3, arr3 = Base.run_blocked_actuation(cfg2, pressure_time=t3, pressure_MPa=p3)
print(f"reference propre (P0=0) : force[0] = {arr3['force_mN'][0]:.2f} mN ; "
      f"a P=0.5 ensuite, extension = {arr3['uncoiled_extension_mm'][2]:.4e} mm")
# defaut : P(0) du profil cyclique
tc, pc = Base.cyclic_pressure_history(1, 1.3, 10.0, 1.5, 0.25, False)
print(f"profil cyclique par defaut : P(0) = {pc[0]} MPa (donc cas P0>0 non atteint par defaut)")
