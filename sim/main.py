import copy
import random
from typing import Dict, List, Optional

import numpy as np

from agents.BCell import BCell
from agents.antigen import Antigen
from config import SIMULATION_PARAMS, SEROTYPES
from environment.immune_system import ImmuneSystem

def build_random_epitope_vectors(receptor_length: int = 5) -> Dict[str, List[float]]:
    return {serotype: [random.random() for _ in range(receptor_length)] for serotype in SEROTYPES}


def build_antigens(epitope_vectors: Dict[str, List[float]]) -> List[Antigen]:
    return [
        Antigen(id=1, serotype="1", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["1"]),
        Antigen(id=2, serotype="5", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["5"]),
        Antigen(id=3, serotype="6B", polysaccharide_ug=4.4, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["6B"]),
        Antigen(id=4, serotype="14", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["14"]),
        Antigen(id=5, serotype="18C", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["18C"]),
        Antigen(id=6, serotype="19F", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["19F"]),
        Antigen(id=7, serotype="23F", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["23F"]),
    ]


def initialize_naive(num_cells: int, receptor_length: int = 5) -> List[BCell]:
    return [BCell(id=i, receptors=[random.random() for _ in range(receptor_length)]) for i in range(num_cells)]


def _set_seed(seed: Optional[int]) -> None:
    if seed is None:
        return
    random.seed(seed)
    np.random.seed(seed)


def simulate_once(params: Dict, seed: Optional[int] = None) -> Dict[str, List[float]]:
    _set_seed(seed)
    sim_params = copy.deepcopy(params)
    epitope_vectors = build_random_epitope_vectors()
    antigens = build_antigens(epitope_vectors)
    initial_pool_size = sim_params.get("naive_pool_size", SIMULATION_PARAMS["naive_pool_size"])
    immune_system = ImmuneSystem(antigens, initialize_naive(initial_pool_size), sim_params)

    antibody_levels = {s.serotype: [] for s in antigens}
    vaccine_days = set(sim_params["vaccination_schedule"])
    schedule = sim_params["vaccination_schedule"]
    checkpoints = {dose_day + 30 for dose_day in schedule}

    for day in range(sim_params["duration_days"]):
        if day in vaccine_days:
            immune_system.vaccinate(antigens)
        if day in checkpoints:
            for serotype in antibody_levels:
                antibody_levels[serotype].append(immune_system.antibody_levels[serotype])
        immune_system.step()

    # Garantiza salida consistente de 2 puntos (pre/post) para calibración,
    # incluso si duration_days no alcanza todos los checkpoints.
    for serotype, values in antibody_levels.items():
        if len(values) >= 2:
            antibody_levels[serotype] = values[:2]
            continue
        final_level = immune_system.antibody_levels[serotype]
        if len(values) == 1:
            antibody_levels[serotype].append(final_level)
        else:
            antibody_levels[serotype] = [final_level, final_level]
    return antibody_levels


def simulate_batch(params: Dict, n_replicas: int = 8, base_seed: Optional[int] = 12345, seeds: Optional[List[int]] = None):
    if seeds is None:
        if base_seed is None:
            rng = random.SystemRandom()
            seeds = [rng.randrange(0, 2**32 - 1) for _ in range(n_replicas)]
        else:
            seeds = [base_seed + i for i in range(n_replicas)]
    outputs = [simulate_once(params, seed=s) for s in seeds]
    serotypes = outputs[0].keys() if outputs else []
    mean_output = {}
    std_output = {}
    for st in serotypes:
        points = np.array([run[st] for run in outputs], dtype=float)
        mean_output[st] = np.mean(points, axis=0).tolist()
        std_output[st] = np.std(points, axis=0).tolist()
    return {
        "runs": outputs,
        "mean": mean_output,
        "std": std_output,
        "seeds": [int(s) for s in seeds],
    }


# Compatibilidad con notebooks existentes
def la(params):
    return simulate_once(params)


def sim(params, times):
    return [simulate_once(params) for _ in range(times)]
