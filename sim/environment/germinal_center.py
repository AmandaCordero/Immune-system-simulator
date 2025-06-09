# environment/germinal_center.py
from agents.BCell import BCell
from agents.antigen import Antigen
from typing import List
from processes.selection import boltzmann_selection
from processes.affinity import compute_affinity
from processes.differentiation import differentiate_bcell
from processes.mutation import mutate_bcell
import random
from config import SIMULATION_PARAMS

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

    def run_cycle(self, p):
        if len(self.bcells) > 0:
            selected = boltzmann_selection(self.bcells)
            
            for cell in selected:
                tipo = differentiate_bcell(cell)  # debe devolver 'memoria', 'plasma' o None
                if tipo == 'memoria':
                    self.memory_cells.append(cell)
                    if cell in self.bcells:
                        self.bcells.remove(cell)
                    continue
                if tipo == 'plasma':
                    self.plasma_cells.append(cell)
                    if cell in self.bcells:
                        self.bcells.remove(cell)
                    continue
                if cell in self.bcells and random.random() < SIMULATION_PARAMS["lf_decay"]:
                    self.bcells.remove(cell)
            
            bc = []
            for cell in self.bcells_pool:
                # Probabilidad de morir
                if random.random() < SIMULATION_PARAMS["lf_decay"]:
                    # La célula muere, no se añade a bc
                    continue
                # Sobrevive, con probabilidad q muta
                if random.random() < SIMULATION_PARAMS["mutation_p"]:
                    bc.append(cell)
                    cell = mutate_bcell(cell)
                    cell.affinity = compute_affinity(self.antigen.epitope_vector, cell.receptors)
                bc.append(cell)
            self.bcells_pool = bc
        
        return (self.memory_cells, self.plasma_cells)
