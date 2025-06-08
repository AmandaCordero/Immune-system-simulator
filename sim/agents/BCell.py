# agents/bcell.py
from dataclasses import dataclass, field
from typing import Dict, List, Optional

@dataclass
class BCell:
    id: int
    receptors: List[float]
    affinity: List[float]
    is_memory: bool = False
    clones: int = 0
    family_id: Optional[int] = None
    mutations: int = 0
    blocked: bool = False
    alive: bool = True
    # Puedes agregar historial de mutaciones, origen, etc.
