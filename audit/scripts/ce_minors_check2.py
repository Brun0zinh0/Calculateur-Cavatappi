# -*- coding: utf-8 -*-
"""Contre-expertise : (a) longueur geometrique en suspendu ; (b) identite alpha V2 vs livrable ;
(c) compliance extremites figee a alpha0 ; (d) swap v12/v13 (ordre de grandeur d'impact)."""
import importlib.util
import sys
import numpy as np

ROOT = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

A = load("alpha_base2", ROOT + r"\alpha V2\Base.py")
L = load("livrable_base2", ROOT + r"\livrable\Base.py")
print(f"alpha V2 : {A.MODEL_VERSION} ; livrable : {L.MODEL_VERSION}")

print()
print("=== (a) Mode suspendu : axial_length_geometry_mm vs axial_length_mm ===")
model, arr = A.run_suspended_actuation(
    load_N=1.0, eps=0.8, n_cycles=1, Pmax=1.2, dt=1.0,
    n_layers=2, n_phi=8, pre_steps=6,
)
geomL = np.asarray(arr["axial_length_geometry_mm"], float)
actL = np.asarray(arr["axial_length_mm"], float)
d = np.abs(geomL - actL)
print(f"suspendu : max|L_geom - L| = {d.max():.4f} mm ({100*d.max()/actL.mean():.2f} %), "
      f"L moyen = {actL.mean():.2f} mm")
# controle du signe de la contraction
print(f"contraction max = {np.max(arr['free_contraction_mm']):.3f} mm")

print("--- meme comparaison en bloque ---")
model_b, arr_b = A.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.2, dt=1.0,
                                         n_layers=2, n_phi=8, pre_steps=6)
geomLb = np.asarray(arr_b["axial_length_geometry_mm"], float)
actLb = np.asarray(arr_b["axial_length_mm"], float)
print(f"bloque : max|L_geom - L| = {np.abs(geomLb-actLb).max():.2e} mm")

print()
print("=== (b) Identite numerique alpha V2 vs livrable (protocole bloque commun) ===")
common = dict(eps=0.8, n_cycles=1, Pmax=1.3, dt=1.0, n_layers=3, n_phi=8, pre_steps=6,
              nonlinear_pressure=True)
_, da = A.run_blocked_actuation(**common)
_, dl = L.run_blocked_actuation(**common)
fa, fl = np.asarray(da["force_mN"], float), np.asarray(dl["force_mN"], float)
ta, tl = np.asarray(da["torque_microNm"], float), np.asarray(dl["torque_microNm"], float)
n = min(len(fa), len(fl))
print(f"n_a={len(fa)} n_l={len(fl)}")
print(f"max|dF| = {np.abs(fa[:n]-fl[:n]).max():.6e} mN ; max|dT| = {np.abs(ta[:n]-tl[:n]).max():.6e} uNm")

print()
print("=== (c) Compliance des extremites : figee a alpha0 ===")
geom = A.default_geometry_params(uncoiled_length=4.0)
disc = A.default_discretization(n_layers=2, n_phi=8, pre_steps=6)
m2 = A.TCPAMaxwellBlockedModel(mat=A.default_material_params(), geom=geom, disc=disc)
c_init = m2.uncoiled_compliance_mm_per_N
m2.prestretch_to(0.8)
c_after = m2.uncoiled_compliance_mm_per_N
alpha_tk = np.rad2deg(m2.helix.alpha)
print(f"compliance init = {c_init:.6e} mm/N ; apres prestretch(0.8) = {c_after:.6e} (inchangee: {c_init==c_after})")
print(f"alpha0 = {np.rad2deg(m2.alpha0):.2f} deg ; alpha_tk = {alpha_tk:.2f} deg")
# compliance recalculee avec alpha courant
half = 0.5 * m2.uncoiled_length
axr, ber = m2._uncoiled_section_rigidities()
a_tk = m2.helix.alpha
c_tk = 2.0 * (half * np.sin(a_tk)**2 / axr + half**3 * np.cos(a_tk)**2 / (3.0 * ber))
print(f"compliance recalculee a alpha_tk = {c_tk:.6e} mm/N ; ecart = {100*(c_tk-c_init)/c_init:+.1f} %")
print(f"cos^2: {np.cos(m2.alpha0)**2:.4f} -> {np.cos(a_tk)**2:.4f}")

print()
print("=== (d) Swap v12/v13 : impact sur un run bloque (resolution moyenne) ===")
orig = A.effective_poissons
def swapped(Cbar):
    EL, v12, v13, v14 = orig(Cbar)
    return EL, v13, v12, v14
_, d_ref = A.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=0.5,
                                   n_layers=4, n_phi=16, pre_steps=12)
A.effective_poissons = swapped
try:
    _, d_sw = A.run_blocked_actuation(eps=0.8, n_cycles=1, Pmax=1.3, dt=0.5,
                                      n_layers=4, n_phi=16, pre_steps=12)
finally:
    A.effective_poissons = orig
f_ref = np.asarray(d_ref["force_mN"], float); f_sw = np.asarray(d_sw["force_mN"], float)
t_ref = np.asarray(d_ref["torque_act_microNm"], float); t_sw = np.asarray(d_sw["torque_act_microNm"], float)
print(f"force max : ref={f_ref.max():.2f} swap={f_sw.max():.2f} ecart={100*(f_sw.max()-f_ref.max())/f_ref.max():+.2f} %")
print(f"couple act max : ref={t_ref.max():.2f} swap={t_sw.max():.2f} ecart={100*(t_sw.max()-t_ref.max())/t_ref.max():+.2f} %")
