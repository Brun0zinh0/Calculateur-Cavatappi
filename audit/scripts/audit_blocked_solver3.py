"""Audit numerique 3 : pas dt=0 initial."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

disc = Base.default_discretization(n_layers=2, n_phi=8, pre_steps=4, dw_bracket=(-0.05, 0.05))
m3 = Base.TCPAMaxwellBlockedModel(disc=disc)
m3.prestretch_to(0.5)
print("max|sigma_i| apres prestretch (attendu 0):", float(np.max(np.abs(m3.sigma_i))))
print("max|sigma_reference| apres prestretch:", float(np.max(np.abs(m3.sigma_reference))))
r0 = m3.step(0.0, 0.0, h_target=m3.h_blocked)
print("pas dt=0 a P=0 : residu:", r0.residual, " dw:", r0.dw, " Ft:", r0.Ft)
r1 = m3.step(0.5, 0.0, h_target=m3.h_blocked)
print("saut instantane 0.5 MPa (dt=0) : dw:", r1.dw, " Ft (N):", r1.Ft, " residu:", r1.residual)
print("max|sigma_i| apres saut:", float(np.max(np.abs(m3.sigma_i))))
