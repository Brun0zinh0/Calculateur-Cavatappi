# -*- coding: utf-8 -*-
"""Contre-expertise des constats minor : checks rapides sur alpha V2."""
import importlib.util
import sys
import numpy as np

ALPHA = r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2"

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

B = load("alpha_base", ALPHA + r"\Base.py")
P = load("alpha_pression", ALPHA + r"\pression.py")

print("=== 1. Poisson effectifs : ordre et valeurs ===")
C_local = B.ti_stiffness_from_paper(37.76, 8.82, 7.24, 0.205, 0.422)
theta = np.deg2rad(37.91)
Cbar = B.rotate_stiffness_bias(C_local, theta)
EL, v12, v13, v14 = B.effective_poissons(Cbar)
print(f"EL={EL:.4f}  v12(=nu s->phi)={v12:.4f}  v13(=nu s->r)={v13:.4f}  v14={v14:.4f}")
print(f"ratio v13/v12 = {v13/v12:.2f}")
# controle independant : compliance 6x6 complete
S6 = np.linalg.inv(Cbar)
eps = S6 @ np.array([1.0, 0, 0, 0, 0, 0])
print(f"controle 6x6 : nu(s->phi)={-eps[1]/eps[0]:.4f}  nu(s->r)={-eps[2]/eps[0]:.4f}")
# verification que la reduction 4x4 est exacte (decouplage epsilon4=epsilon5=0)
print(f"couplages residuels eps4={eps[3]:.2e} eps5={eps[4]:.2e}")

print()
print("=== 2. Degenerescence theta -> 0 ===")
C0r = B.rotate_stiffness_bias(C_local, 0.0)
num = C0r[0, 1] - C0r[0, 2]
den = C0r[2, 2] - C0r[1, 1]
print(f"theta=0 : C12-C13 = {num:.3e}  C33-C22 = {den:.3e}")
with np.errstate(invalid="ignore", divide="ignore"):
    print(f"B = (C12-C13)/(C33-C22) = {num/den!r}")
# validation accepte theta_f_deg=0 ?
geom = B.default_geometry_params(theta_f_deg=0.0)
disc = B.default_discretization(n_layers=2, n_phi=8, pre_steps=4)
mat = B.default_material_params()
try:
    with np.errstate(all="ignore"):
        model = B.TCPAMaxwellBlockedModel(mat=mat, geom=geom, disc=disc)
        print("construction theta_f=0 : OK (validation acceptee)")
        A_, B_, mu_ = model._layer_AB_mu(0)
        print(f"couche 0 : A={A_!r} B={B_!r} mu={mu_!r}")
        try:
            model.prestretch_to(0.2)
            print("prestretch: OK ??")
        except Exception as exc:
            print(f"prestretch echoue : {type(exc).__name__}: {exc}")
except Exception as exc:
    print(f"construction refusee : {type(exc).__name__}: {exc}")
# angle minimal atteint avec defauts
model_def = B.TCPAMaxwellBlockedModel(mat=B.default_material_params(),
                                      geom=B.default_geometry_params(),
                                      disc=B.default_discretization(n_layers=4, n_phi=8, pre_steps=4))
print(f"defauts n_layers=4 : angles de couche (deg) = {np.rad2deg(model_def.theta_layers)}")
den_layers = [model_def.C_total[j][2,2]-model_def.C_total[j][1,1] for j in range(4)]
print(f"denominateurs C33-C22 par couche = {den_layers}")
# petit angle : conditionnement
for th in (1.0, 0.1, 0.01):
    Ct = B.rotate_stiffness_bias(C_local, np.deg2rad(th))
    n_ = Ct[0,1]-Ct[0,2]; d_ = Ct[2,2]-Ct[1,1]
    print(f"theta={th:5.2f} deg : num={n_:.3e} den={d_:.3e} B={n_/d_:.4f} mu={np.sqrt(Ct[1,1]/Ct[2,2]):.6f}")

print()
print("=== 3. pression.py ===")
print(f"infer_force_unit('column1')      = {P.infer_force_unit('column1')}")
print(f"infer_force_unit('Force (N)')    = {P.infer_force_unit('Force (N)')}")
print(f"infer_force_unit('force_mn')     = {P.infer_force_unit('force_mn')}")
psi_exact = 4.4482216152605 / (0.0254**2) / 1e6  # lbf / in^2 en MPa
print(f"facteur psi code = 0.006894757293168361 ; exact = {psi_exact:.18f}")
cols = {"t": np.array([0.0, 1.0, 2.0]), "p": np.array([2.0, 5.0, 13.0])}
out = P.measured_pressure_payload(cols, "t", "p", "bar", subtract_initial=True)
print(f"bar avec offset 2.0 soustrait -> {out['pressure_MPa']}")
try:
    P.measured_pressure_payload({"t": np.array([0., 1., 2.]), "p": np.array([0., 8., 16.])}, "t", "p", "bar", False)
    print("16 bar accepte ??")
except ValueError as exc:
    print(f"16 bar rejete : {exc}")
# troncature a zero
cols2 = {"t": np.array([0.0, 1.0, 2.0]), "p": np.array([0.5, 0.2, 1.0])}
out2 = P.measured_pressure_payload(cols2, "t", "p", "MPa", subtract_initial=True)
print(f"zero puis troncature : {out2['pressure_MPa']} (0.2-0.5<0 tronque a 0)")

print()
print("=== 4. Historiques de pression : defauts et timings ===")
t, p = B.cyclic_pressure_history()  # appel direct, tous defauts
half = 60.0 * 1.5 / 10.0
print(f"demi-periode attendue = {half} s")
i_quarter = np.argmin(np.abs(t - 0.5 * half))
print(f"appel direct cyclic_pressure_history() : P(t=demi-charge)={p[i_quarter]:.4f} "
      f"(lineaire attendrait {1.3*0.5:.4f}, nonlin 1.3*0.5^3.5={1.3*0.5**3.5:.4f})")
peaks = p[np.isclose(t % (2*half), half)]
print(f"pics a demi-periodes : {peaks}")
cfg = B.default_simulation_config()
print(f"default_simulation_config().nonlinear_pressure = {cfg.nonlinear_pressure}")
print(f"defaut cyclic_pressure_history nonlinear = True (signature)")
Pm = load("alpha_parametres", ALPHA + r"\parametres.py")
print(f"parametres.DEFAULTS['nonlinear_pressure'] = {Pm.DEFAULTS.get('nonlinear_pressure') if hasattr(Pm,'DEFAULTS') else '???'}")
import dataclasses
sp_fields = {f.name: f.default for f in dataclasses.fields(Pm.SimulationParams)}
print(f"SimulationParams.nonlinear_pressure default = {sp_fields.get('nonlinear_pressure')}")
# run_hold_relaxation force la rampe non lineaire ?
import inspect
src = inspect.getsource(B.run_hold_relaxation)
print("run_hold_relaxation passe nonlinear_ramp ? ->", "nonlinear" in src)
t2, p2 = B.ramp_hold_pressure_history(P_hold=1.3, ramp_time=9.0, hold_time=10.0, dt=0.5)
i45 = np.argmin(np.abs(t2 - 4.5))
print(f"ramp_hold defaut : P(4.5 s)={p2[i45]:.4f} (nonlin 1.3*0.5^3.5={1.3*0.5**3.5:.4f}, lin=0.65)")
print(f"palier atteint a t=9 : P={p2[np.argmin(np.abs(t2-9.0))]:.4f}")
# timings du profil lineaire
t3, p3 = B.cyclic_pressure_history(n_cycles=2, Pmax=1.3, nonlinear=False)
print(f"lineaire : P(0)={p3[0]}, pics={p3[np.isclose(t3 % (2*half), half)]}, P(fin)={p3[-1]}, "
      f"pente={(p3[np.argmin(np.abs(t3-4.5))]/4.5):.5f} vs 1.3/9={1.3/9:.5f}")
