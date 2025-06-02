# main.py
from config import *
from agents.antigen import Antigen
from environment.immune_system import ImmuneSystem

def main():
    # 1. Inicializa antígenos (serotipos de PCV)
    antigens = [
        Antigen(serotype="4", charge=0.7, hydrophobicity=0.6, size=0.8),
        # ...
    ]
    # 2. Inicializa sistema inmune
    immune_system = ImmuneSystem(antigens)
    # 3. Simula esquema de vacunación y seguimiento
    for day in range(SIMULATION_PARAMS["duration_days"]):
        if day in SIMULATION_PARAMS["vaccination_schedule"]:
            for antigen in antigens:
                immune_system.vaccinate(antigen, dose=5)
        immune_system.step(days=1)
    # 4. Guarda y analiza resultados
    # ...

if __name__ == "__main__":
    main()
