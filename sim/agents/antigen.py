# agents/antigen.py
from dataclasses import dataclass
from typing import List

@dataclass
class Antigen:
    id:int
    serotype: str
    polysaccharide_ug: float          # Cantidad de polisacárido en µg en la vacuna
    carrier_protein: str              # Ej: 'TT', 'CRM197'
    immunogenicity_factor: float      # Parámetro relativo (ej. 1.0 normal, >1 más inmunogénico)
    epitope_vector: List[float]
