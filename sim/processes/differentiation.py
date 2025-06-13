# processes/differentiation.py
from config import SIMULATION_PARAMS

def differentiate_bcell(bcell, affinity_threshold_memory = SIMULATION_PARAMS["affinity_threshold_memory"], affinity_threshold_plasma=SIMULATION_PARAMS["affinity_threshold_plasma"]):
    """
    Decide el destino de una célula B tras selección, basado en su afinidad.
    
    Args:
        bcell: instancia de BCell con atributo 'affinity' (float entre 0 y 1).
        affinity_threshold_memory: afinidad mínima para diferenciarse en célula de memoria.
        affinity_threshold_plasma: afinidad mínima para diferenciarse en célula plasmática.
        
    Returns:
        destino: str, uno de ['memory', 'plasma', 'NO']
    """
    if bcell.affinity >= affinity_threshold_memory:
        return 'memory'
    elif bcell.affinity >= affinity_threshold_plasma:
        return 'plasma'
    else:
        return 'NO'
