import copy
import argparse

from config import SIMULATION_PARAMS
from main import simulate_once


def parse_args():
    parser = argparse.ArgumentParser(description="Run one mock immune simulation.")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible runs. Omit for stochastic behavior.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    params = copy.deepcopy(SIMULATION_PARAMS)

    # Mock rápido: ajusta estos valores libremente
    params["duration_days"] = 140
    params["vaccination_schedule"] = [20, 40, 100]
    params["plasma_production_factor"] = {
        "1": 0.015,
        "5": 0.013,
        "6B": 0.012,
        "14": 0.009,
        "18C": 0.007,
        "19F": 0.011,
        "23F": 0.012,
    }
    params["memory_production_factor"] = {
        "1": 0.004,
        "5": 0.0035,
        "6B": 0.003,
        "14": 0.0025,
        "18C": 0.0018,
        "19F": 0.0032,
        "23F": 0.003,
    }
    params["mutation_p"] = 0.30
    params["temperature"] = 0.12
    params["THRESHOLD"] = 0.28

    out = simulate_once(params, seed=args.seed)

    print("Salida por serotipo (pre, post):")
    print("seed:", args.seed)
    for serotype in sorted(out.keys()):
        print(serotype, out[serotype])


if __name__ == "__main__":
    main()
