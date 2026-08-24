"""Impact du profil d'angle de biais (lineaire en angle vs arctan Eq. 11) sur F/T bloques."""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail\alpha V2")
import Base

for profile in ("paper_linear", "uniform_twist"):
    geom = Base.default_geometry_params(bias_angle_profile=profile)
    m, arr = Base.run_blocked_actuation(
        eps=0.8, n_cycles=1, Pmax=1.3, dt=1.0, n_layers=6, n_phi=16, pre_steps=8, geom=geom
    )
    print(f"{profile:14s} angles couches (deg): {np.rad2deg(m.theta_layers).round(2)}")
    print(f"{'':14s} F baseline {arr['force_mN'][0]:8.1f} mN | F pic {arr['force_mN'].max():8.1f} mN"
          f" | T_act pic {arr['torque_act_microNm'].max():7.1f} uNm | T_act min {arr['torque_act_microNm'].min():8.1f} uNm")
