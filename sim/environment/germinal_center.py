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
        self.cycles : int = 0

    def seed_naive_cells(self, naive_pool: List[BCell]):
        self.bcells = naive_pool

    def run_cycle(self):
        self.cycles += 1
        if len(self.bcells) > 0:
            self.memory_cells = []
            self.plasma_cells = []
            selected = boltzmann_selection(self.bcells, max_survivors= int(self.cycles**2 / 1000 * len(self.bcells)))
            with open("log_simulacion.txt", "a") as f:
                f.write(f"selected {len(selected)}\n")

            for cell in selected:
                # with open("log_simulacion.txt", "a") as f:
                #     f.write(f"afinidad {cell.affinity}\n")
                tipo = differentiate_bcell(cell)  # debe devolver 'memoria', 'plasma' o None
                if tipo == 'memory':
                    self.memory_cells.append(cell)
                    if cell in self.bcells:
                        self.bcells.remove(cell)
                    continue
                if tipo == 'plasma':
                    self.plasma_cells.append(cell)
                    if cell in self.bcells:
                        self.bcells.remove(cell)
                    continue
            

            self.bcells = [cell for cell in self.bcells if random.random() > SIMULATION_PARAMS["lf_decay"]]

            

            new = []
            for cell in self.bcells:
                if random.random() < SIMULATION_PARAMS["mutation_p"]:
                    m = mutate_bcell(cell)
                    m.affinity = compute_affinity(self.antigen.epitope_vector, m.receptors)
                    
                    new.append(m)
            self.bcells += new
            with open("log_simulacion.txt", "a") as f:
                f.write(f"b desp mutar {len(self.bcells)}\n")
            # bc = []
            # for cell in self.bcells:
            #     # Probabilidad de morir
            #     if random.random() < SIMULATION_PARAMS["lf_decay"]:
            #         with open("log_simulacion.txt", "a") as f:
            #             f.write("muere\n")
            #         # La célula muere, no se añade a bc
            #         continue
            #     # Sobrevive, con probabilidad q muta
            #     if random.random() < SIMULATION_PARAMS["mutation_p"]:
            #         with open("log_simulacion.txt", "a") as f:
            #             f.write("muta\n")
            #         bc.append(cell)
            #         cell = mutate_bcell(cell)
            #         cell.affinity = compute_affinity(self.antigen.epitope_vector, cell.receptors)
            #     bc.append(cell)
            # with open("log_simulacion.txt", "a") as f:
            #     f.write(f"{len(bc)}\n")
            # self.bcells_pool = bc
        
        return (self.memory_cells, self.plasma_cells)
