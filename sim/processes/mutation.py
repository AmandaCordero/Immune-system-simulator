# processes/mutation.py
from random import random
import numpy as np
from agents.BCell import BCell

def mutate_bcell(bcell, mutation_rate=0.1, mutation_strength=0.05):
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
        affinity=bcell.affinity.copy(),
        is_memory=bcell.is_memory,
        clones=bcell.clones,
        family_id=bcell.family_id,
        mutations=bcell.mutations + 1,
        blocked=bcell.blocked,
        alive=bcell.alive
    )
    return mutated_bcell
