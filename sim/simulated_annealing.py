import copy
import csv
import hashlib
import json
import math
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from config import CALIBRATION_PARAMS, SIMULATION_PARAMS
from main import simulate_batch

SEROTYPES = {
    "1": ("pre PsC SP sp1", "post PsC SP sp1"),
    "14": ("pre PsC SP sp14", "post PsC SP sp14"),
    "18C": ("pre PsC SP sp18C", "post PsC SP sp18C"),
    "19F": ("pre PsC SP sp19F", "post PsC SP sp19F"),
    "23F": ("pre PsC SP sp23F", "post PsC SP sp23F"),
    "5": ("pre PsC SP sp5", "post PsC SP sp5"),
    "6B": ("pre PsC SP sp6B", "post PsC SP sp6B"),
}


def chamfer_distance(points1, points2):
    points1 = np.squeeze(points1)
    points2 = np.squeeze(points2)
    tree1 = cKDTree(points1)
    tree2 = cKDTree(points2)
    dist1, _ = tree2.query(points1)
    dist2, _ = tree1.query(points2)
    return float(np.mean(dist1) + np.mean(dist2))


def load_clinical_data(path="sim/data/train_VCN7-Tf_fit.csv"):
    return pd.read_csv(path)


def _score_serotype(sim_points, real_points, chamfer_w, landmark_w):
    chamfer = chamfer_distance(real_points, sim_points)
    sim_mean = np.mean(np.array(sim_points), axis=0)
    real_mean = np.mean(np.array(real_points), axis=0)
    landmark = float(np.linalg.norm(sim_mean - real_mean))
    return chamfer_w * chamfer + landmark_w * landmark, chamfer, landmark


def compute_objective(sim_summary, clinical_df, weights):
    chamfer_w = CALIBRATION_PARAMS["chamfer_weight"]
    landmark_w = CALIBRATION_PARAMS["landmark_weight"]

    serotype_scores = {}
    total_weight = 0.0
    total_score = 0.0
    for serotype, cols in SEROTYPES.items():
        pre_vals = clinical_df[cols[0]].values
        post_vals = clinical_df[cols[1]].values
        real_points = np.array([(pre_vals[i], post_vals[i]) for i in range(len(pre_vals))], dtype=float)

        run_points = sim_summary["runs"]
        sim_points = np.array([(run_points[i][serotype][0], run_points[i][serotype][1]) for i in range(len(run_points))], dtype=float)
        score, chamfer, landmark = _score_serotype(sim_points, real_points, chamfer_w, landmark_w)
        w = weights.get(serotype, 1.0)
        serotype_scores[serotype] = {
            "weighted_score": float(score * w),
            "score": float(score),
            "chamfer": float(chamfer),
            "landmark_l2": float(landmark),
            "weight": float(w),
        }
        total_weight += w
        total_score += score * w
    return float(total_score / max(total_weight, 1e-9)), serotype_scores


def _vector_keys():
    return list(CALIBRATION_PARAMS["param_bounds"].keys())


def _set_nested_param(params, dotted_key, value):
    parts = dotted_key.split(".")
    target = params
    for part in parts[:-1]:
        target = target[part]
    target[parts[-1]] = float(value)


def _get_nested_param(params, dotted_key):
    parts = dotted_key.split(".")
    value = params
    for part in parts:
        value = value[part]
    return float(value)


def vecino(parametros):
    nuevo = copy.deepcopy(parametros)
    bounds = CALIBRATION_PARAMS["param_bounds"]
    key = random.choice(_vector_keys())
    low, high = bounds[key]
    span = high - low
    delta = random.uniform(-0.12 * span, 0.12 * span)
    current_value = _get_nested_param(nuevo, key)
    mutated_value = float(min(high, max(low, current_value + delta)))
    _set_nested_param(nuevo, key, mutated_value)
    if nuevo["affinity_threshold_plasma"] > nuevo["affinity_threshold_memory"] - 0.02:
        nuevo["affinity_threshold_plasma"] = max(
            bounds["affinity_threshold_plasma"][0],
            nuevo["affinity_threshold_memory"] - 0.02,
        )
    return nuevo


def _param_hash(params):
    payload = {k: float(params[k]) if isinstance(params[k], (int, float, np.floating)) else params[k] for k in sorted(params.keys())}
    encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _evaluate(params, clinical_df, objective_config, cache):
    key = _param_hash(params)
    if key in cache:
        return cache[key], True

    t0 = time.perf_counter()
    summary = simulate_batch(
        params,
        n_replicas=objective_config["n_replicas"],
        base_seed=objective_config["base_seed"],
    )
    energy, details = compute_objective(summary, clinical_df, objective_config["serotype_weights"])
    elapsed = time.perf_counter() - t0
    result = {"energy": energy, "details": details, "summary": summary, "elapsed_sec": elapsed}
    cache[key] = result
    return result, False


def recocido_simulado(param_inicial, objective_config=None):
    if objective_config is None:
        objective_config = CALIBRATION_PARAMS
    annealing_cfg = objective_config["annealing"]
    clinical_df = load_clinical_data()

    current = copy.deepcopy(param_inicial)
    cache = {}
    current_eval, _ = _evaluate(current, clinical_df, objective_config, cache)
    best = copy.deepcopy(current)
    best_eval = current_eval

    T = annealing_cfg["temp_initial"]
    i = 0
    no_improvement = 0
    jsonl_path = Path("misc/annealing_log.jsonl")
    csv_path = Path("misc/annealing_log.csv")
    csv_exists = csv_path.exists()

    while T > annealing_cfg["temp_final"] and i < annealing_cfg["max_iter"]:
        candidate = vecino(current)
        cand_eval, from_cache = _evaluate(candidate, clinical_df, objective_config, cache)

        delta = cand_eval["energy"] - current_eval["energy"]
        accepted = False
        if delta < 0 or random.random() < math.exp(-delta / max(T, 1e-9)):
            current = candidate
            current_eval = cand_eval
            accepted = True

        if cand_eval["energy"] < best_eval["energy"]:
            best = copy.deepcopy(candidate)
            best_eval = cand_eval
            no_improvement = 0
        else:
            no_improvement += 1

        row = {
            "iter": i,
            "temp": T,
            "accepted": accepted,
            "from_cache": from_cache,
            "energy_current": current_eval["energy"],
            "energy_candidate": cand_eval["energy"],
            "best_energy": best_eval["energy"],
            "elapsed_candidate_sec": cand_eval["elapsed_sec"],
            "params_candidate": candidate,
            "serotype_scores": cand_eval["details"],
        }
        with jsonl_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(row) + "\n")

        flat_row = {
            "iter": i,
            "temp": T,
            "accepted": accepted,
            "from_cache": from_cache,
            "energy_current": current_eval["energy"],
            "energy_candidate": cand_eval["energy"],
            "best_energy": best_eval["energy"],
            "elapsed_candidate_sec": cand_eval["elapsed_sec"],
        }
        for k in _vector_keys():
            flat_row[k] = _get_nested_param(candidate, k)
        with csv_path.open("a", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(flat_row.keys()))
            if not csv_exists:
                writer.writeheader()
                csv_exists = True
            writer.writerow(flat_row)

        if no_improvement >= annealing_cfg["patience"]:
            break
        T *= annealing_cfg["alpha"]
        i += 1

    return best, best_eval


if __name__ == "__main__":
    base = copy.deepcopy(SIMULATION_PARAMS)
    best_params, best_eval = recocido_simulado(base)
    print("Best energy:", best_eval["energy"])
    print("Best params:", best_params)
