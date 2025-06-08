# processes/selection.py
import numpy as np


def boltzmann_selection(bcells, temperature=0.1, max_survivors=None):
    """
    Selección estocástica tipo Boltzmann basada en afinidad.
    
    Args:
        bcells (list): lista de objetos BCell con atributo 'affinity'.
        temperature (float): parámetro que regula la presión selectiva (menor temp = selección más estricta).
        max_survivors (int): número máximo de células que pueden sobrevivir (recurso limitado).
        
    Returns:
        list: células B seleccionadas para sobrevivir.
    """
    affinities = np.array([b.affinity for b in bcells])
    
    # Evitar overflow: normalizar afinidades restando el máximo
    affinities_norm = affinities - np.max(affinities)
    
    # Calcular probabilidades tipo Boltzmann
    probs = np.exp(affinities_norm / temperature)
    probs /= probs.sum()
    
    # Número de sobrevivientes
    if max_survivors is None or max_survivors > len(bcells):
        max_survivors = len(bcells)
    
    # Selección sin reemplazo según probabilidades
    selected_indices = np.random.choice(len(bcells), size=max_survivors, replace=False, p=probs)
    
    selected = [bcells[i] for i in selected_indices]
    return selected

