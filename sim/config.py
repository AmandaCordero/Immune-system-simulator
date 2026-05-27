# config.py
import numpy as np
SEROTYPES = ["1", "5", "6B", "14", "18C", "19F", "23F"]

SIMULATION_PARAMS = {
    "duration_days": 400,  # 2 años
    # "time_step": 1,        # días
    "vaccination_schedule": [60, 80, 120],  # días
    "decay_factor" : np.exp(-0.033),  # 0.033 ≈ ln(2)/21
    "plasma_production_factor": {s: 0.01 for s in SEROTYPES},
    "memory_production_factor": {s: 0.002 for s in SEROTYPES},
    "lf_decay": 0.15,
    "mutation_p": 0.2,
    "affinity_threshold_memory":0.6,
    "affinity_threshold_plasma":0.4,
    "mutation_rate":0.7,
    "mutation_strength":0.5,
    "temperature":0.1,
    "THRESHOLD" : 0.3,
    "naive_pool_size": 4000,
    "naive_recruitment_per_dose": 1200,
    "min_selection_fraction": 0.02,
    "max_selection_fraction": 0.35,
    "selection_growth_rate": 120.0,
    "max_gc_pool": 8000,
    "selection_schedule": {
        "base_temperature": 0.15,
        "min_temperature": 0.04,
        "decay_per_cycle": 0.005,
    },
    "mutation_schedule": {
        "base_strength": 0.45,
        "min_strength": 0.08,
        "decay_per_cycle": 0.004,
    },
}

CALIBRATION_PARAMS = {
    "param_bounds": {
        "lf_decay": (0.05, 0.35),
        "mutation_p": (0.05, 0.6),
        "mutation_rate": (0.2, 0.95),
        "mutation_strength": (0.05, 0.8),
        "temperature": (0.04, 0.35),
        "THRESHOLD": (0.15, 0.55),
        "affinity_threshold_memory": (0.5, 0.85),
        "affinity_threshold_plasma": (0.25, 0.6),
    },
    "n_replicas": 8,
    "base_seed": 12345,
    "serotype_weights": {
        "1": 1.0,
        "5": 1.0,
        "6B": 1.0,
        "14": 1.3,
        "18C": 1.0,
        "19F": 1.0,
        "23F": 1.0,
    },
    "chamfer_weight": 0.75,
    "landmark_weight": 0.25,
    "annealing": {
        "temp_initial": 1.2,
        "temp_final": 0.04,
        "alpha": 0.94,
        "max_iter": 120,
        "patience": 20,
    },
}

for ser in SEROTYPES:
    CALIBRATION_PARAMS["param_bounds"][f"plasma_production_factor.{ser}"] = (0.001, 0.06)
    CALIBRATION_PARAMS["param_bounds"][f"memory_production_factor.{ser}"] = (0.0, 0.03)
