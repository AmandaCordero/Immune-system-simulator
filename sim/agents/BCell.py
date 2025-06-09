# agents/bcell.py
from dataclasses import dataclass, field
from typing import Dict, List, Optional

@dataclass
class BCell:
    id: int
    receptors: List[float]
    affinity: float = 0.0
    serotype: str = ""
    