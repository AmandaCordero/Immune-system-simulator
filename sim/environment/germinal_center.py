# environment/germinal_center.py
from agents.BCell import BCell
from agents.antigen import Antigen
from typing import List
from processes.selection import boltzmann_selection
from processes.affinity import compute_affinity
from processes.differentiation import differentiate_bcell
from processes.mutation import mutate_bcell

class GerminalCenter:
    def __init__(self, id: int, antigen: Antigen):
        self.id = id
        self.antigen = antigen
        self.bcells: List[BCell] = []
        self.memory_cells: List[BCell] = []
        self.plasma_cells: List[BCell] = []
        self.factors: float = 1.0  # Recursos limitantes

    def seed_naive_cells(self, naive_pool: List[BCell]):
        self.bcells = naive_pool

    def run_cycle(self):
        if (len(self.bcells) > 0):
            selected = boltzmann_selection(self.bcells)
            for cell in selected:
                # diferenciar y sacar de gc (muerte por apoptosis)
                pass
            bc = self.bcells.copy()
            # division
            self.bcells = [mutate_bcell(cell) for cell in bc]
        return (self.memory_cells, self.plasma_cells)
