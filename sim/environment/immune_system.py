# environment/immune_system.py
from typing import List
from environment.germinal_center import GerminalCenter
from agents.antigen import Antigen

class ImmuneSystem:
    def __init__(self, antigens: List[Antigen]):
        self.gcs = [GerminalCenter(i, ag) for i, ag in enumerate(antigens)]
        self.memory_pool = []
        self.plasma_pool = []
        self.antibody_levels = {ag.serotype: 0.0 for ag in antigens}
        # ...

    def vaccinate(self, antigen: Antigen, dose: int):
        # Distribuye células naïve a los GCs correspondientes
        pass

    def step(self, days: int):
        # Ejecuta ciclos en todos los GCs, actualiza pools y niveles de anticuerpos
        pass
