# config.py
import numpy as np

SIMULATION_PARAMS = {
    "duration_days": 730,  # 2 años
    # "time_step": 1,        # días
    "vaccination_schedule": [30, 150, 300],  # días
    "decay_factor" : np.exp(-0.033),  # 0.033 ≈ ln(2)/21
    "plasma_production_factor": 0.01,
    "memory_production_factor": 0.002,
    "lf_decay": 0.15,
    "mutation_p": 0.2,
    "affinity_threshold_memory":0.6,
    "affinity_threshold_plasma":0.4,
    "mutation_rate":0.7,
    "mutation_strength":0.5,
    "temperature":0.1,
    "THRESHOLD" : 0.3
}
