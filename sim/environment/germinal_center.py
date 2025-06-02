# environment/germinal_center.py
from agents.BCell import BCell
from agents.antigen import Antigen
from typing import List

class GerminalCenter:
    def __init__(self, id: int, antigen: Antigen):
        self.id = id
        self.antigen = antigen
        self.bcells: List[BCell] = []
        self.memory_cells: List[BCell] = []
        self.plasma_cells: List[BCell] = []
        self.factors: float = 1.0  # Recursos limitantes

    def seed_naive_cells(self, naive_pool: List[BCell]):
        # Selecciona y agrega células naïve al GC
        pass

    def run_cycle(self, mutation_engine, selection_process, differentiation_process):
        # 1. Mutación
        # 2. Selección (Boltzmann)
        # 3. Diferenciación (memoria/plasmática)
        pass
