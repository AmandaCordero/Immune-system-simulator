# processes/mutation.py
import random
import numpy as np
from agents.BCell import BCell
from config import SIMULATION_PARAMS

def mutate_bcell(bcell: BCell, mutation_rate=SIMULATION_PARAMS["mutation_rate"], mutation_strength=SIMULATION_PARAMS["mutation_strength"]):
    """
    Aplica hipermutación somática a la célula B modificando su vector receptor.
    
    Args:
        bcell (BCell): célula B a mutar.
        mutation_rate (float): probabilidad de mutar cada componente del vector.
        mutation_strength (float): desviación estándar de la mutación gaussiana.
    
    Retorna:
        BCell mutada (nueva instancia o modificada in-place).
    """
    new_receptors = bcell.receptors.copy()
    
    for i in range(len(new_receptors)):
        if random.random() < mutation_rate:
            # Añade un cambio gaussiano pequeño
            delta = np.random.normal(0, mutation_strength)
            new_receptors[i] += delta
            # Opcional: limitar valores a rango válido (ej. 0 a 1)
            new_receptors[i] = min(max(new_receptors[i], 0.0), 1.0)
    
    # Crear una nueva célula B mutada (puedes modificar in-place si prefieres)
    mutated_bcell = BCell(
        id=bcell.id,
        receptors=new_receptors,
        affinity=0.0,
        serotype=bcell.serotype
    )
    return mutated_bcell
