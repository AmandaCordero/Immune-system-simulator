# agents/antigen.py
from dataclasses import dataclass
from typing import List

@dataclass
class Antigen:
    serotype: str
    polysaccharide_ug: float          # Cantidad de polisacárido en µg en la vacuna
    carrier_protein: str              # Ej: 'TT', 'CRM197'
    immunogenicity_factor: float      # Parámetro relativo (ej. 1.0 normal, >1 más inmunogénico)
    protective_threshold: float       # Nivel mínimo de anticuerpos (µg/ml) para protección
    virulence_index: float = 1.0      # Opcional para modelar competencia o prevalencia
    epitope_vector: List[float]
