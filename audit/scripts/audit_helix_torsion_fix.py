# -*- coding: utf-8 -*-
"""Torsion et courbure de l'helice par derivees analytiques exactes (Frenet)."""
import numpy as np

rho, alpha = 2.16, np.deg2rad(10.53)
h = rho * np.tan(alpha)
t0 = 0.7
r1 = np.array([-rho * np.sin(t0), rho * np.cos(t0), h])
r2 = np.array([-rho * np.cos(t0), -rho * np.sin(t0), 0.0])
r3 = np.array([rho * np.sin(t0), -rho * np.cos(t0), 0.0])
kappa = np.linalg.norm(np.cross(r1, r2)) / np.linalg.norm(r1) ** 3
tau = float(np.cross(r1, r2) @ r3) / float(np.cross(r1, r2) @ np.cross(r1, r2))
print(f"kappa Frenet = {kappa:.10f} ; cos^2a/rho   = {np.cos(alpha)**2/rho:.10f} ; ecart = {abs(kappa-np.cos(alpha)**2/rho):.2e}")
print(f"tau   Frenet = {tau:.10f} ; sin2a/(2rho) = {np.sin(2*alpha)/(2*rho):.10f} ; ecart = {abs(tau-np.sin(2*alpha)/(2*rho)):.2e}")
