# config.py
import numpy as np

GC_PARAMS = {
    "num_gcs": 5,
    "cycles_per_gc": 4,
    "lf_decay": 0.1,
    "THRESHOLD" : 0.3
    # ...
}

MUTATION_PARAMS = {
    "cdr_rate": 0.1,
    "fwr_rate": 0.02,
    "lethal_prob": 0.01,
    "block_prob": 0.05,
    # ...
}

SELECTION_PARAMS = {
    "boltzmann_temp": 0.1,
    "selection_pressure": 0.7,
    # ...
}

SIMULATION_PARAMS = {
    "duration_days": 730,  # 2 años
    "time_step": 1,        # días
    "vaccination_schedule": [30, 60, 300],  # días
    "decay_factor" : np.exp(-0.033),  # 0.033 ≈ ln(2)/21
    "plasma_production_factor": 0.01,
    "memory_production_factor": 0.002,
    "lf_decay": 0.1,
    "mutation_p": 0.1,
    # ...
}
