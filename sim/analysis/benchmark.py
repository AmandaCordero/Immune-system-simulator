import copy
import time

import numpy as np

from config import CALIBRATION_PARAMS, SIMULATION_PARAMS
from main import simulate_batch
from simulated_annealing import SEROTYPES, compute_objective, load_clinical_data


def _objective_for_weights(weights):
    clinical_df = load_clinical_data()
    t0 = time.perf_counter()
    summary = simulate_batch(
        copy.deepcopy(SIMULATION_PARAMS),
        n_replicas=CALIBRATION_PARAMS["n_replicas"],
        base_seed=CALIBRATION_PARAMS["base_seed"],
    )
    energy, detail = compute_objective(summary, clinical_df, weights)
    return energy, detail, time.perf_counter() - t0


def benchmark():
    one_weights = {k: 0.0 for k in SEROTYPES}
    one_weights["14"] = 1.0
    all_weights = {k: 1.0 for k in SEROTYPES}

    one_energy, _, one_time = _objective_for_weights(one_weights)
    all_energy, _, all_time = _objective_for_weights(all_weights)

    print("1-serotype objective")
    print("energy:", one_energy)
    print("seconds:", one_time)
    print("7-serotype objective")
    print("energy:", all_energy)
    print("seconds:", all_time)


def check_reproducibility():
    params = copy.deepcopy(SIMULATION_PARAMS)
    a = simulate_batch(params, n_replicas=4, base_seed=900)
    b = simulate_batch(params, n_replicas=4, base_seed=900)
    c = simulate_batch(params, n_replicas=4, base_seed=1200)
    same_seed_delta = 0.0
    diff_seed_delta = 0.0
    for st in a["mean"]:
        same_seed_delta += np.linalg.norm(np.array(a["mean"][st]) - np.array(b["mean"][st]))
        diff_seed_delta += np.linalg.norm(np.array(a["mean"][st]) - np.array(c["mean"][st]))
    print("same-seed delta:", float(same_seed_delta))
    print("diff-seed delta:", float(diff_seed_delta))


if __name__ == "__main__":
    benchmark()
    check_reproducibility()
