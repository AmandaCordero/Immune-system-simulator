from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class TCell:
    id: int
    specificity: List[float]  # Perfil de reconocimiento de epítopos (puedes usar la misma estructura que el antígeno)
    activation_state: str = "naive"  # naive, activated, exhausted
    help_capacity: float = 1.0       # Cuánta "ayuda" puede dar en un ciclo GC
    memory: bool = False
    interactions: int = 0            # Número de células B ayudadas en el ciclo actual
