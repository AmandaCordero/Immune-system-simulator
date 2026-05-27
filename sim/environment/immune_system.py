# environment/immune_system.py
from typing import List
from environment.germinal_center import GerminalCenter
from agents.antigen import Antigen
from agents.BCell import BCell
from processes.affinity import compute_affinity
from config import SIMULATION_PARAMS
import numpy as np
import random

class ImmuneSystem:
    def __init__(self, antigens: List[Antigen], bcells_pool = List[BCell], params = SIMULATION_PARAMS):
        self.gcs = []
        self.bcells_pool = bcells_pool
        self.antigens = antigens
        self.memory_pool = {ag.serotype: [] for ag in antigens}
        self.plasma_pool = {ag.serotype: [] for ag in antigens}
        self.antibody_levels = {ag.serotype: 0.0 for ag in antigens}
        self.params = params
        self.receptor_length = len(antigens[0].epitope_vector) if antigens else 5

    def _factor_for_serotype(self, factor_name: str, serotype: str) -> float:
        value = self.params[factor_name]
        if isinstance(value, dict):
            return float(value.get(serotype, 0.0))
        if isinstance(value, (list, tuple, np.ndarray)):
            try:
                index = [ag.serotype for ag in self.antigens].index(serotype)
            except ValueError:
                return 0.0
            if index >= len(value):
                return 0.0
            return float(value[index])
        return float(value)

    def fillCGs(self, antigens: List[Antigen]):
        self.gcs = []
        self.gcs = [GerminalCenter(id=ag.id, antigen=ag, params=self.params) for ag in antigens]
        assigned_bcells = self._assign_bcells_vectorized(antigens)
        for gc in self.gcs:
            gc.seed_naive_cells(assigned_bcells[gc.id])

    def _assign_bcells_vectorized(self, antigens: List[Antigen]):
        assigned_bcells = {ag.id: [] for ag in antigens}
        if not self.bcells_pool:
            return assigned_bcells

        receptors = np.array([bc.receptors for bc in self.bcells_pool], dtype=float)
        epitopes = np.array([ag.epitope_vector for ag in antigens], dtype=float)

        diff = receptors[:, None, :] - epitopes[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        affinity_matrix = np.exp(-dist)

        max_aff = np.max(affinity_matrix, axis=1)
        threshold = self.params["THRESHOLD"]
        candidate_mask = max_aff >= threshold
        if not np.any(candidate_mask):
            return assigned_bcells

        candidate_rows = np.where(candidate_mask)[0]
        assigned_cols = np.argmax(affinity_matrix[candidate_mask], axis=1)

        for row_idx, col_idx in zip(candidate_rows, assigned_cols):
            bc = self.bcells_pool[row_idx]
            ag = antigens[int(col_idx)]
            bc.serotype = ag.serotype
            bc.affinity = float(affinity_matrix[row_idx, col_idx])
            assigned_bcells[ag.id].append(bc)
        return assigned_bcells
 
    def cycle(self):
        for ag in self.plasma_pool.keys():
            self.plasma_pool[ag] = [cell for cell in self.plasma_pool[ag] if random.random() >= self.params["decay_factor"]]

    def initialize_naive(self, num_cells=None, receptor_length=None) -> List['BCell']:
            if num_cells is None:
                num_cells = self.params["naive_pool_size"]
            if receptor_length is None:
                receptor_length = self.receptor_length
            """Genera células B naive con receptores aleatorios"""
            start_id = len(self.bcells_pool)
            return [
                BCell(
                    id=start_id + i,
                    receptors=[random.random() for _ in range(receptor_length)]
                ) for i in range(num_cells)
            ]


    def vaccinate(self, antigens: List[Antigen]):
        recruits = self.initialize_naive(num_cells=self.params["naive_recruitment_per_dose"])
        self.bcells_pool.extend(recruits)
        self.fillCGs(antigens)
            
    def step(self):
        for gc in self.gcs:
            memory, plasma = gc.run_cycle()
            for cell in memory:
                self.memory_pool[cell.serotype].append(cell)
            for cell in plasma:
                self.plasma_pool[cell.serotype].append(cell)
        
        for ag in self.antibody_levels.keys():
            plasma_factor = self._factor_for_serotype("plasma_production_factor", ag)
            memory_factor = self._factor_for_serotype("memory_production_factor", ag)
            plasma_production = len(self.plasma_pool[ag]) * plasma_factor
            memory_production = len(self.memory_pool[ag]) * memory_factor
            self.antibody_levels[ag] = (
                (self.antibody_levels[ag] + plasma_production + memory_production) * self.params["decay_factor"]
            )

        self.cycle()
