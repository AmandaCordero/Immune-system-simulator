import numpy as np
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
from dataclasses import dataclass

# =====================
# ESTRUCTURAS DE DATOS
# =====================

@dataclass
class Antigen:
    """Estructura que representa un antígeno bacteriano"""
    serotype: str          # Ej: "6B", "23F"
    charge: float          # Carga eléctrica (afecta unión a anticuerpos)
    hydrophobicity: float  # Hidrofobicidad superficial
    size: float            # Tamaño molecular

@dataclass
class BCell:
    """Célula B con receptores específicos"""
    receptors: List[float] # Perfil de unión a antígenos
    affinity: float = 0.0   # Afinidad actual (0-1)
    is_memory: bool = False # ¿Es célula de memoria?
    clones: int = 0         # Número de clones generados

# =====================
# MÓDULO DE AFINIDAD
# =====================

class AffinityCalculator:
    """Calcula fuerza de unión antígeno-anticuerpo"""
    
    @staticmethod
    def compute(antigen: Antigen, bcells: List[BCell]) -> List[float]:
        """Calcula afinidad para todas las células B"""
        return [
            np.exp(-np.sqrt(
                (antigen.charge - b.receptors[0])**2 +
                (antigen.hydrophobicity - b.receptors[1])**2 +
                (antigen.size - b.receptors[2])**2
            ))
            for b in bcells
        ]

# =====================
# MÓDULO DE MUTACIÓN
# =====================

class MutationEngine:
    """Simula hipermutación somática en centros germinales"""
    
    def __init__(self, mutation_rate: float = 0.1):
        self.rate = mutation_rate
    
    def mutate(self, cell: BCell) -> BCell:
        """Aplica mutaciones aleatorias a los receptores"""
        new_receptors = [
            r + np.random.normal(0, 0.1) if np.random.random() < self.rate else r
            for r in cell.receptors
        ]
        return BCell(
            receptors=new_receptors,
            is_memory=cell.is_memory,
            clones=cell.clones + 1
        )

# =====================
# MÓDULO DE SELECCIÓN
# =====================

class SelectionProcess:
    """Selección darwiniana en centros germinales"""
    
    @staticmethod
    def select(antigen: Antigen, bcells: List[BCell], threshold: float = 0.7) -> Tuple[List[BCell], List[BCell]]:
        """Divide células en seleccionadas (alta afinidad) y descartadas"""
        affinities = AffinityCalculator.compute(antigen, bcells)
        selected = []
        discarded = []
        
        for b, aff in zip(bcells, affinities):
            b.affinity = aff
            (selected if aff > threshold else discarded).append(b)
        
        return selected, discarded

# =====================
# SISTEMA INMUNE
# =====================

class ImmuneSystem:
    """Modelo completo del sistema inmunológico"""
    
    def __init__(self):
        self.memory_pool: Dict[str, List[BCell]] = {}  # Células de memoria por serotipo
        self.plasma_cells: Dict[str, List[BCell]] = {} # Células productoras de anticuerpos
        self.antibody_levels: Dict[str, float] = {}    # Niveles actuales de IgG
        
        # Módulos
        self.mutator = MutationEngine()
        self.selector = SelectionProcess()
    
    def initialize_serotype(self, serotype: str):
        """Prepara estructuras para un nuevo serotipo"""
        self.memory_pool[serotype] = []
        self.plasma_cells[serotype] = []
        self.antibody_levels[serotype] = 0.0
    
    def generate_naive_cells(self, antigen: Antigen, n_cells: int = 10) -> List[BCell]:
        """Genera células B naive con receptores aleatorios"""
        return [
            BCell(receptors=[
                antigen.charge + np.random.normal(0, 0.2),
                antigen.hydrophobicity + np.random.normal(0, 0.2),
                antigen.size + np.random.normal(0, 0.1)
            ])
            for _ in range(n_cells)
        ]
    
    def germinal_center_reaction(self, antigen: Antigen, bcells: List[BCell], cycles: int = 3):
        """Simula múltiples ciclos de mutación y selección"""
        selected_cells = bcells
        
        for _ in range(cycles):
            # 1. Mutación
            mutated = [self.mutator.mutate(c) for c in selected_cells]
            
            # 2. Selección
            selected_cells, _ = self.selector.select(antigen, mutated)
            
            # 3. Proliferación (aumentar pool)
            selected_cells.extend([self.mutator.mutate(c) for c in selected_cells[:2]])
        
        return selected_cells
    
    def vaccinate(self, antigen: Antigen, dose: int = 1):
        """Administra una dosis de vacuna"""
        if antigen.serotype not in self.memory_pool:
            self.initialize_serotype(antigen.serotype)
        
        for _ in range(dose):
            # 1. Generar células naive
            naive_cells = self.generate_naive_cells(antigen)
            
            # 2. Reacción en centro germinal
            selected = self.germinal_center_reaction(antigen, naive_cells)
            
            # 3. Diferenciación:
            for cell in selected:
                if cell.affinity > 0.8:  # Alta afinidad -> memoria
                    cell.is_memory = True
                    self.memory_pool[antigen.serotype].append(cell)
                else:  # Media afinidad -> plasmablastos
                    self.plasma_cells[antigen.serotype].append(cell)
    
    def update_antibodies(self, days: int = 1):
        """Actualiza niveles de anticuerpos en el tiempo"""
        for serotype in self.plasma_cells:
            # Producción de anticuerpos
            plasma_production = len(self.plasma_cells[serotype]) * 0.1
            memory_production = len(self.memory_pool[serotype]) * 0.02
            
            # Decaimiento (vida media ~21 días)
            decay_factor = np.exp(-0.033 * days)  # 0.033 ≈ ln(2)/21
            
            self.antibody_levels[serotype] = (
                (self.antibody_levels[serotype] + plasma_production + memory_production) * decay_factor
            )
    
    def time_step(self, days: int = 1):
        """Avanza la simulación en el tiempo"""
        self.update_antibodies(days)
        
        # Eliminar células viejas (turnover)
        for serotype in self.plasma_cells:
            self.plasma_cells[serotype] = [
                c for c in self.plasma_cells[serotype]
                if np.random.random() > 0.01 * days  # 1% de muerte celular por día
            ]

# =====================
# SIMULACIÓN Y VISUALIZACIÓN
# =====================

def run_simulation():
    """Ejemplo completo de simulación"""
    # 1. Configuración inicial
    system = ImmuneSystem()
    pcv_serotypes = [
        Antigen(serotype="4", charge=0.7, hydrophobicity=0.6, size=0.8),
        Antigen(serotype="6B", charge=0.5, hydrophobicity=0.7, size=1.0),
        Antigen(serotype="23F", charge=0.65, hydrophobicity=0.55, size=0.75)
    ]
    
    # 2. Esquema de vacunación (0, 2, 6 meses)
    for antigen in pcv_serotypes:
        system.vaccinate(antigen)  # Dosis inicial
    
    # 3. Seguimiento por 2 años
    time_points = []
    antibody_levels = {s.serotype: [] for s in pcv_serotypes}
    
    for month in range(24):
        system.time_step(days=30)
        
        # Refuerzos (2da y 3ra dosis)
        if month == 2 or month == 6:
            for antigen in pcv_serotypes:
                system.vaccinate(antigen)
        
        # Registrar datos
        time_points.append(month)
        for serotype in antibody_levels:
            antibody_levels[serotype].append(system.antibody_levels[serotype])
    
    # 4. Visualización
    plt.figure(figsize=(10, 6))
    for serotype, levels in antibody_levels.items():
        plt.plot(time_points, levels, label=f"Serotipo {serotype}")
    
    plt.axhline(y=0.35, color='gray', linestyle='--', label="Umbral protector")
    plt.title("Respuesta Inmune a Vacuna Conjugada")
    plt.xlabel("Meses post-vacunación")
    plt.ylabel("Nivel de IgG (µg/mL)")
    plt.legend()
    plt.grid()
    plt.show()

    # =====================
# SIMULACIÓN Y VISUALIZACIÓN (VERSIÓN DIARIA)
# =====================

def run_daily_simulation():
    """Versión diaria de la simulación"""
    # 1. Configuración inicial (igual)
    system = ImmuneSystem()
    pcv_serotypes = [
        Antigen(serotype="4", charge=0.7, hydrophobicity=0.6, size=0.8),
        Antigen(serotype="6B", charge=0.5, hydrophobicity=0.7, size=1.0),
        Antigen(serotype="23F", charge=0.65, hydrophobicity=0.55, size=0.75)
    ]
    
    # 2. Esquema de vacunación ajustado a días
    days_schedule = [90, 180, 365]  
    # for day in days_schedule:
    #     for antigen in pcv_serotypes:
    #         system.vaccinate(antigen)
    
    # 3. Seguimiento por 720 días (2 años)
    total_days = 720
    time_points = list(range(total_days))
    antibody_levels = {s.serotype: [] for s in pcv_serotypes}
    
    for antigen in pcv_serotypes:
        system.vaccinate(antigen, dose = 0)

    for current_day in range(total_days):
        # Avanzar un día
        system.time_step(days=1)
        
        if current_day in days_schedule:
            for antigen in pcv_serotypes:
                system.vaccinate(antigen)

        # Registrar datos diariamente
        for serotype in antibody_levels:
            antibody_levels[serotype].append(system.antibody_levels[serotype])
    
    # 4. Visualización ajustada
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
    run_daily_simulation()