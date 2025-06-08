# processes/selection.py
import numpy as np

def agrupar_por_serotipo(bcells):
    grupos = {}
    for b in bcells:
        if not b.affinity:
            continue  # Ignorar células sin afinidad
        serotipo_dominante = max(b.affinity, key=b.affinity.get)
        if serotipo_dominante not in grupos:
            grupos[serotipo_dominante] = []
        grupos[serotipo_dominante].append(b)
    return grupos

def seleccion_por_serotipo(bcells, antigenos_activos):
    grupos = agrupar_por_serotipo(bcells)
    seleccionadas = []
    for serotipo, grupo in grupos.items():
        if serotipo not in antigenos_activos:
            continue  # Solo considerar serotipos presentes en la vacuna
        afinidades = [b.affinity[serotipo] for b in grupo]
        seleccion_serotipo = boltzmann_selection(grupo, temperature=0.1, affinities=afinidades)
        seleccionadas.extend(seleccion_serotipo)
    return seleccionadas

def es_cross_reactiva(bcell, umbral=0.7):
    return sum(afinidad > umbral for afinidad in bcell.affinity.values()) >= 2

def seleccion_con_cross_reactividad(bcells, antigenos_activos):
    seleccionadas = seleccion_por_serotipo(bcells, antigenos_activos)
    cross_reactivas = [b for b in bcells if es_cross_reactiva(b)]
    # Añadir cross-reactivas con prioridad (ej. 20% adicional de probabilidad)
    seleccionadas.extend(cross_reactivas)
    return list(set(seleccionadas))  # Eliminar duplicados


def boltzmann_selection(affinities, temperature=0.1, max_survivors=None):
    """
    Selección estocástica tipo Boltzmann basada en afinidad.
    
    Args:
        bcells (list): lista de objetos BCell con atributo 'affinity'.
        temperature (float): parámetro que regula la presión selectiva (menor temp = selección más estricta).
        max_survivors (int): número máximo de células que pueden sobrevivir (recurso limitado).
        
    Returns:
        list: células B seleccionadas para sobrevivir.
    """
    
    # Evitar overflow: normalizar afinidades restando el máximo
    affinities_norm = affinities - np.max(affinities)
    
    # Calcular probabilidades tipo Boltzmann
    probs = np.exp(affinities_norm / temperature)
    probs /= probs.sum()
    
    # Número de sobrevivientes
    if max_survivors is None or max_survivors > len(affinities):
        max_survivors = len(affinities)
    
    # Selección sin reemplazo según probabilidades
    selected_indices = np.random.choice(len(affinities), size=max_survivors, replace=False, p=probs)
    
    return selected_indices











# En GerminalCenter.run_cycle()
def run_cycle(self, antigenos_activos):
    # ... (mutación, cálculo de afinidad)
    # Selección
    bcells_seleccionadas = seleccion_con_cross_reactividad(self.bcells, antigenos_activos)
    self.bcells = bcells_seleccionadas
