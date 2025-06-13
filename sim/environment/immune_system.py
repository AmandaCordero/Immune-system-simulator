# environment/immune_system.py
from typing import List
from environment.germinal_center import GerminalCenter
from agents.antigen import Antigen
from agents.BCell import BCell
from processes.affinity import compute_affinity
from processes.mutation import mutate_bcell
from config import SIMULATION_PARAMS
import numpy as np
import simpy
import random

class ImmuneSystem:
    def __init__(self, antigens: List[Antigen], bcells_pool = List[BCell]):
        self.gcs = []
        self.bcells_pool = bcells_pool
        self.antigens = antigens
        self.memory_pool = {ag.serotype: [] for ag in antigens}
        self.plasma_pool = {ag.serotype: [] for ag in antigens}
        self.antibody_levels = {ag.serotype: 0.0 for ag in antigens}
        
    def fillCGs(self, antigens: List[Antigen]):
        self.gcs = []
        self.gcs = [GerminalCenter(id=ag.id, antigen=ag)for ag in antigens]
        
        env = simpy.Environment()

        bcell_store = simpy.FilterStore(env, capacity=len(self.bcells_pool))
        bcell_store.items = self.bcells_pool.copy()

        assigned_bcells = {ag.id: [] for ag in antigens}
        access_lock = simpy.Resource(env, capacity=1)

        def antigen_process(ag: Antigen):
            while True:
                with access_lock.request() as req:
                    yield req
                    if not bcell_store.items:
                        break

                    bc = random.choice(bcell_store.items)
                    affinity = compute_affinity(ag.epitope_vector, bc.receptors)
                    print(affinity)

                    if affinity >= SIMULATION_PARAMS["THRESHOLD"]:
                        assigned_bcells[ag.id].append(bc)
                        bc.serotype = ag.serotype
                        bc.affinity = affinity
                        yield bcell_store.get(lambda x: x == bc)
                        # Espera tras un match exitoso
                        print("a")
                        yield env.timeout(0.1)  
                    else:
                        # No hay unión, no espera y sigue intentando inmediatamente
                        print("b")
                        pass

        for ag in antigens:
            env.process(antigen_process(ag))

        env.run(until=200)

        for gc in self.gcs:
            gc.seed_naive_cells(assigned_bcells[gc.id])
            # print(gc.bcells)
        
        self.bcells_pool = list(bcell_store.items)
 
    def cycle(self):
        # bc = []
        # for cell in self.bcells_pool:
        #     # Probabilidad de morir
        #     if random.random() < SIMULATION_PARAMS["lf_decay"]:
        #         # La célula muere, no se añade a bc
        #         continue
        #     # Sobrevive, con probabilidad q muta
        #     if random.random() < SIMULATION_PARAMS["mutation_p"]:
        #         bc.append(cell)
        #         cell = mutate_bcell(cell)
        #     bc.append(cell)
        # self.bcells_pool = bc
        # pass
        # for ag in self.memory_pool.keys():
        #     self.memory_pool[ag] = [cell for cell in self.memory_pool[ag] if random.random() >= SIMULATION_PARAMS["decay_factor"]]
        for ag in self.plasma_pool.keys():
            self.plasma_pool[ag] = [cell for cell in self.plasma_pool[ag] if random.random() >= SIMULATION_PARAMS["decay_factor"]]

    def initialize_naive(self,num_cells: int = 10000, receptor_length: int = 5) -> List['BCell']:
            """Genera células B naive con receptores aleatorios"""
            return [
                BCell(
                    id=i,
                    receptors=[random.random() for _ in range(receptor_length)]
                ) for i in range(num_cells)
            ]


    def vaccinate(self, antigens: List[Antigen]):
        self.bcells_pool = self.initialize_naive()
        self.fillCGs(antigens)
        # for gc in self.gcs:
        #     with open("log_simulacion.txt", "a") as f:
        #         f.write(f"serotipo {gc.antigen.serotype}, b antes{len(gc.bcells)}\n")
        #     memory, plasma = gc.run_cycle()
        #     with open("log_simulacion.txt", "a") as f:
        #         f.write(f"serotipo {gc.antigen.serotype}, b despues{len(gc.bcells)}\n")
        #         f.write(f"memoria {len(memory)}\n")
        #         f.write(f"plasma {len(plasma)}\n")
        #     for cell in memory:
        #         self.memory_pool[cell.serotype].append(cell)
        #     for cell in plasma:
        #         self.plasma_pool[cell.serotype].append(cell)
            
    def step(self):
        for gc in self.gcs:
            with open("log_simulacion.txt", "a") as f:
                f.write(f"serotipo {gc.antigen.serotype}, b antes{len(gc.bcells)}\n")
            
            memory, plasma = gc.run_cycle()
            with open("log_simulacion.txt", "a") as f:
                f.write(f"serotipo {gc.antigen.serotype}, b despues{len(gc.bcells)}\n")
                f.write(f"memoria {len(memory)}\n")
                f.write(f"plasma {len(plasma)}\n")
            for cell in memory:
                self.memory_pool[cell.serotype].append(cell)
            for cell in plasma:
                self.plasma_pool[cell.serotype].append(cell)
        
        for ag in self.antibody_levels.keys():
            plasma_production = len(self.plasma_pool[ag]) * SIMULATION_PARAMS["plasma_production_factor"]
            memory_production = len(self.memory_pool[ag]) * SIMULATION_PARAMS["memory_production_factor"]
            self.antibody_levels[ag] = (
                (self.antibody_levels[ag] + plasma_production + memory_production) * SIMULATION_PARAMS["decay_factor"]
            )

        self.cycle()
