# processes/affinity.py
import numpy as np

def compute_affinity(antigen_epitopo, bcell_receptor):
    # bcell_receptor y antigen_epitopo son vectores numéricos
    # distancia euclídea y función exponencial
    distancia = np.linalg.norm(np.array(bcell_receptor) - np.array(antigen_epitopo))
    afinidad = np.exp(-distancia)  # Afinidad entre 0 y 1, mayor si distancia es pequeña
    return afinidad
