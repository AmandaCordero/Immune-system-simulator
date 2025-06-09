# main.py
from matplotlib import pyplot as plt
from config import *
from agents.antigen import Antigen
from environment.immune_system import ImmuneSystem

def main():
    # 1. Inicializa antígenos (serotipos de PCV)
    from agents.antigen import Antigen

    epitope_vectors = {
        "1":  [0.72, 0.15, 0.40, 0.85, 0.50],
        "5":  [0.68, 0.20, 0.35, 0.80, 0.55],
        "6B": [0.75, 0.18, 0.42, 0.88, 0.48],
        "14": [0.70, 0.22, 0.38, 0.82, 0.52],
        "18C": [0.65, 0.25, 0.36, 0.78, 0.58],
        "19F": [0.73, 0.17, 0.41, 0.84, 0.49],
        "23F": [0.69, 0.19, 0.37, 0.81, 0.53],
    }

    
    # Inicialización de la vacuna con los antígenos indicados
    antigens = [
        Antigen(serotype="1", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["1"]),
        Antigen(serotype="5", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["5"]),
        Antigen(serotype="6B", polysaccharide_ug=4.4, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["6B"]),
        Antigen(serotype="14", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["14"]),
        Antigen(serotype="18C", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["18C"]),
        Antigen(serotype="19F", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["19F"]),
        Antigen(serotype="23F", polysaccharide_ug=2.2, carrier_protein="TT", immunogenicity_factor=1.0, epitope_vector=epitope_vectors["23F"]),
    ]

    # 2. Inicializa sistema inmune
    immune_system = ImmuneSystem(antigens)
    # 3. Simula esquema de vacunación y seguimiento
    
    antibody_levels = {s.serotype: [] for s in antigens}
    time_points = list(range(SIMULATION_PARAMS["duration_days"]))

    for day in range(SIMULATION_PARAMS["duration_days"]):
        if day in SIMULATION_PARAMS["vaccination_schedule"]:
            for antigen in antigens:
                immune_system.vaccinate(antigen, dose=5)
        # Registrar datos diariamente
        for serotype in antibody_levels:
            antibody_levels[serotype].append(immune_system.antibody_levels[serotype])
        immune_system.step()
    
    plt.figure(figsize=(12, 7))
    for serotype, levels in antibody_levels.items():
        plt.plot(time_points, levels, label=f"Serotipo {serotype}", alpha=0.8)
    
    plt.axhline(y=0.35, color='gray', linestyle='--', label="Umbral protector")
    plt.title("Respuesta Inmune Diaria a Vacuna Conjugada")
    plt.xlabel("Días post-vacunación")
    plt.ylabel("Nivel de IgG (µg/mL)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    main()
