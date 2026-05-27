# processes/selection.py
import numpy as np
from config import SIMULATION_PARAMS

def boltzmann_selection(bcells, temperature=SIMULATION_PARAMS["temperature"], max_survivors=None):
    """
    Selección estocástica tipo Boltzmann basada en afinidad.
    
    Args:
        bcells (list): lista de objetos BCell con atributo 'affinity'.
        temperature (float): parámetro que regula la presión selectiva (menor temp = selección más estricta).
        max_survivors (int): número máximo de células que pueden sobrevivir (recurso limitado).
        
    Returns:
        list: células B seleccionadas para sobrevivir.
    """
    if not bcells:
        return []
    if max_survivors is not None and max_survivors <= 0:
        return []

    affinities = np.array([b.affinity for b in bcells], dtype=float)
    
    # Evitar overflow: normalizar afinidades restando el máximo
    affinities_norm = affinities - np.max(affinities)
    
    # Calcular probabilidades tipo Boltzmann
    temp = max(float(temperature), 1e-6)
    probs = np.exp(affinities_norm / temp)
    probs /= probs.sum()
    
    # Número de sobrevivientes
    if max_survivors is None or max_survivors > len(bcells):
        max_survivors = len(bcells)
    
    # Selección sin reemplazo según probabilidades
    selected_indices = np.random.choice(
        len(bcells),
        size=max_survivors,
        replace=False,
        p=probs,
    )
    
    selected = [bcells[i] for i in selected_indices]
    return selected
