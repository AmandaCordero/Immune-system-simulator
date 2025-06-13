import math
import random
import numpy as np
from scipy.spatial import cKDTree
from config import SIMULATION_PARAMS
from main import sim
import pandas as pd

serotipos = {
    'sp1': ('pre PsC SP sp1', 'post PsC SP sp1'),
    'sp14': ('pre PsC SP sp14', 'post PsC SP sp14'),
    'sp18C': ('pre PsC SP sp18C', 'post PsC SP sp18C'),
    'sp19F': ('pre PsC SP sp19F', 'post PsC SP sp19F'),
    'sp23F': ('pre PsC SP sp23F', 'post PsC SP sp23F'),
    'sp5': ('pre PsC SP sp5', 'post PsC SP sp5'),
    'sp6B': ('pre PsC SP sp6B', 'post PsC SP sp6B')
}


def vecino(parametros):
    nuevo = parametros.copy()
    # Modificar un parámetro aleatoriamente dentro de un rango
    clave = random.choice(["plasma_production_factor", "memory_production_factor"])
    delta = random.uniform(-0.005, 0.005)  # ajuste pequeño
    nuevo[clave] = max(0, min(0.1, nuevo[clave] + delta))  # asegurar rango [0,1]
    return nuevo

def recocido_simulado(param_inicial, temp_inicial=1.0, temp_final=1e-3, alpha=0.9, max_iter=100):
    estado_actual = param_inicial
    energia_actual = metrica(sim(estado_actual, 48))
    T = temp_inicial
    mejor_estado = estado_actual
    mejor_energia = energia_actual

    while T > temp_final:
        for _ in range(max_iter):
            estado_nuevo = vecino(estado_actual)
            energia_nueva = metrica(sim(estado_nuevo))
            delta = energia_nueva - energia_actual
            if delta < 0 or random.random() < math.exp(-delta / T):
                estado_actual = estado_nuevo
                energia_actual = energia_nueva
                if energia_nueva < mejor_energia:
                    mejor_estado = estado_nuevo
                    mejor_energia = energia_nueva
        T *= alpha  # enfriamiento
    return mejor_estado, mejor_energia

def metrica(antibodies):
    df = pd.read_csv('data/train_VCN7-T.csv')
    distancias_chamfer = 0

    for serotipo, (col_pre, col_post) in serotipos.items():
        # Extraer valores pre y post para ese serotipo
        pre_vals = df[col_pre].values.reshape(-1, 1)
        post_vals = df[col_post].values.reshape(-1, 1)
        
        points1 = [(pre_vals[i], post_vals[i]) for i in range(len(pre_vals))]
        points2 = [(antibodies[i][serotipo][0], antibodies[i][serotipo][0]) for i in range(len(antibodies))]

        distancias_chamfer += chamfer_distance(points1, points2)

    return distancias_chamfer/7





def chamfer_distance(points1, points2):
    """
    Calcula la distancia de Chamfer entre dos nubes de puntos 2D.
    points1, points2: arrays numpy de forma (N,2) y (M,2)
    Devuelve un valor escalar que indica la similitud (menor es más parecido).
    """
    tree1 = cKDTree(points1)
    tree2 = cKDTree(points2)

    # Para cada punto en points1, distancia al punto más cercano en points2
    dist1, _ = tree2.query(points1)
    # Para cada punto en points2, distancia al punto más cercano en points1
    dist2, _ = tree1.query(points2)

    # Promedio de ambas direcciones
    chamfer_dist = np.mean(dist1) + np.mean(dist2)
    return 1/(1 +chamfer_dist)


if __name__ == "__main__":
    a = recocido_simulado(SIMULATION_PARAMS)
    print(a)