# environment/immune_system.py
from typing import List
from environment.germinal_center import GerminalCenter
from agents.antigen import Antigen
from agents.BCell import BCell
from processes.affinity import compute_affinity
from config import GC_PARAMS, SIMULATION_PARAMS
import numpy as np
import simpy
import random

class ImmuneSystem:
    def __init__(self, antigens: List[Antigen], bcells_pool = List[BCell]):
        self.gcs = []
        self.bcells_pool = bcells_pool
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

                    if affinity >= GC_PARAMS.THRESHOLD:
                        assigned_bcells[ag.id].append(bc)
                        yield bcell_store.get(lambda x: x == bc)
                        # Espera tras un match exitoso
                        yield env.timeout(0.1)  
                    else:
                        # No hay unión, no espera y sigue intentando inmediatamente
                        pass

        for ag in antigens:
            env.process(antigen_process(ag))

        env.run(until=100)

        for gc in self.gcs:
            gc.seed_naive_cells(assigned_bcells[gc.id])
        
        self.bcells_pool = list(bcell_store.items)
 
    def cycle(self):
        # division de celulas B
        pass

    def vaccinate(self, antigens: List[Antigen]):
        self.fillCGs(antigens)
        for gc in self.gcs:
            memory, plasma = gc.run_cycle()
            for ag in self.memory_pool.keys:
                self.memory_pool[ag].append(memory[ag])
            for ag in self.plasma_pool.keys:
                self.plasma_pool[ag].append(plasma[ag])

    def step(self):
        for gc in self.gcs:
            memory, plasma = gc.run_cycle()
            for ag in self.memory_pool.keys:
                self.memory_pool[ag].append(memory[ag])
            for ag in self.plasma_pool.keys:
                self.plasma_pool[ag].append(plasma[ag])
        
        for ag in self.antibody_levels.keys:
            plasma_production = self.plasma_pool[ag] * SIMULATION_PARAMS.plasma_production_factor
            memory_production = self.memory_pool[ag] * SIMULATION_PARAMS.memory_production_factor
            self.antibody_levels[ag] = (
                (self.antibody_levels[ag] + plasma_production + memory_production) * SIMULATION_PARAMS.decay_factor
            )

        self.cycle()
