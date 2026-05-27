import copy
import json
import time
from pathlib import Path

from config import CALIBRATION_PARAMS, SIMULATION_PARAMS
from simulated_annealing import recocido_simulado


def main():
    print("[short-annealing] preparando corrida...")
    params = copy.deepcopy(SIMULATION_PARAMS)
    params["duration_days"] = 140
    cfg = copy.deepcopy(CALIBRATION_PARAMS)

    cfg["n_replicas"] = 2
    cfg["annealing"]["max_iter"] = 5
    cfg["annealing"]["patience"] = 5
    cfg["annealing"]["temp_initial"] = 0.8
    cfg["annealing"]["temp_final"] = 0.1
    cfg["base_seed"] = None

    print("[short-annealing] configuración:")
    print(f"  duration_days={params['duration_days']}")
    print(f"  n_replicas={cfg['n_replicas']}")
    print(f"  base_seed={cfg['base_seed']}")
    print(f"  max_iter={cfg['annealing']['max_iter']}")
    print(f"  patience={cfg['annealing']['patience']}")
    print(f"  temp_initial={cfg['annealing']['temp_initial']}")
    print(f"  temp_final={cfg['annealing']['temp_final']}")
    print("[short-annealing] ejecutando recocido...")
    t0 = time.perf_counter()
    best_params, best_eval = recocido_simulado(params, objective_config=cfg)
    elapsed = time.perf_counter() - t0
    print(f"[short-annealing] terminado en {elapsed:.2f}s")

    print("Best energy:", best_eval["energy"])
    print("Best params:")
    for key in sorted(best_params.keys()):
        print(f"  {key}: {best_params[key]}")

    out_path = Path("sim/analysis/params_short_best.json")
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(best_params, f, indent=2)
    print(f"[short-annealing] best params guardados en: {out_path}")


if __name__ == "__main__":
    main()
