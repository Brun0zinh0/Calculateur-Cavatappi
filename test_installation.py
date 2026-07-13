from __future__ import annotations

import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import livrable_base_calcul as modele
import livrable_parametres as parametres


def main() -> None:
    base_path = Path(modele.BASE_PATH).resolve()
    expected_local_base = (APP_DIR / "Base.py").resolve()

    print("Dossier livrable :", APP_DIR)
    print("Base.py chargé   :", base_path)
    if base_path != expected_local_base:
        print("Attention : le moteur chargé n'est pas le Base.py local au dossier livrable.")
        print("Cela reste accepté si Base.py est volontairement dans le dossier parent.")

    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update(
        {
            "duration_s": 20.0,
            "n_cycles": 1,
            "dt": 5.0,
            "n_layers": 1,
            "n_phi": 4,
            "pre_steps": 2,
        }
    )

    config = parametres.build_config(settings)
    pressure_time, pressure_mpa = parametres.make_pressure_history(config)
    _, data = modele.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    resume = modele.summary(data)

    print("Mini-calcul OK")
    print(f"  points calculés : {len(data['time'])}")
    print(f"  pression max    : {float(data['pressure_MPa'].max()):.3f} MPa")
    print(f"  force max       : {resume['force_max_mN']:.3f} mN")


if __name__ == "__main__":
    main()
