# environment/germinal_center.py
from agents.BCell import BCell
from agents.antigen import Antigen
from typing import List
from processes.selection import boltzmann_selection
from processes.affinity import compute_affinity
from processes.differentiation import differentiate_bcell
from processes.mutation import mutate_bcell
import random
import math
from config import SIMULATION_PARAMS

class GerminalCenter:
    def __init__(self, id: int, antigen: Antigen, params = SIMULATION_PARAMS):
        self.id = id
        self.antigen = antigen
        self.bcells: List[BCell] = []
        self.memory_cells: List[BCell] = []
        self.plasma_cells: List[BCell] = []
        self.factors: float = 1.0  # Recursos limitantes
        self.cycles : int = 0
        self.params = params

    def seed_naive_cells(self, naive_pool: List[BCell]):
        self.bcells = naive_pool

    def run_cycle(self):
        self.cycles += 1
        if len(self.bcells) > 0:
            self.memory_cells = []
            self.plasma_cells = []
            selection_schedule = self.params["selection_schedule"]
            base_temp = selection_schedule["base_temperature"]
            min_temp = selection_schedule["min_temperature"]
            temp_decay = selection_schedule["decay_per_cycle"]
            adaptive_temperature = max(
                min_temp,
                base_temp * math.exp(-temp_decay * max(self.cycles - 1, 0)),
            )

            min_fraction = self.params["min_selection_fraction"]
            max_fraction = self.params["max_selection_fraction"]
            growth = self.params["selection_growth_rate"]
            frac = min_fraction + (max_fraction - min_fraction) * (
                1.0 - math.exp(-self.cycles / max(growth, 1e-6))
            )
            survivors = max(1, int(frac * len(self.bcells)))
            survivors = min(survivors, len(self.bcells))
            selected = boltzmann_selection(
                self.bcells,
                max_survivors=survivors,
                temperature=adaptive_temperature,
            )
            # with open("log_simulacion.txt", "a") as f:
            #     f.write(f"selected {len(selected)}\n")

            for cell in selected:
                # with open("log_simulacion.txt", "a") as f:
                #     f.write(f"afinidad {cell.affinity}\n")
                tipo = differentiate_bcell(cell, affinity_threshold_memory= self.params["affinity_threshold_memory"], affinity_threshold_plasma=self.params["affinity_threshold_plasma"])  # debe devolver 'memoria', 'plasma' o None
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
            

            self.bcells = [cell for cell in self.bcells if random.random() > self.params["lf_decay"]]

            

            mutation_schedule = self.params["mutation_schedule"]
            base_strength = mutation_schedule["base_strength"]
            min_strength = mutation_schedule["min_strength"]
            strength_decay = mutation_schedule["decay_per_cycle"]
            cycle_strength = max(
                min_strength,
                base_strength * math.exp(-strength_decay * max(self.cycles - 1, 0)),
            )

            new = []
            for cell in self.bcells:
                if random.random() < self.params["mutation_p"]:
                    affinity_factor = max(0.1, 1.0 - cell.affinity)
                    adaptive_strength = cycle_strength * affinity_factor
                    m = mutate_bcell(
                        cell,
                        mutation_rate=self.params["mutation_rate"],
                        mutation_strength=adaptive_strength,
                    )
                    m.affinity = compute_affinity(self.antigen.epitope_vector, m.receptors)

                    new.append(m)
            self.bcells += new

            max_gc_pool = self.params["max_gc_pool"]
            if len(self.bcells) > max_gc_pool:
                random.shuffle(self.bcells)
                self.bcells = self.bcells[:max_gc_pool]
            # with open("log_simulacion.txt", "a") as f:
            #     f.write(f"b desp mutar {len(self.bcells)}\n")
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
