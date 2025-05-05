import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parámetros mejorados según el documento
params = {
    # Sistema inmunitario materno
    "IgG_maternal_init": 1.2,       # Nivel inicial IgG1 (principal transferida)
    "decay_maternal": 0.033,        # Vida media 21 días (ln(2)/21 ≈ 0.033)
    "IgA_lactation": 0.8,           # Nivel constante durante lactancia
    
    # Componentes inmunológicos bebé
    "k_dendritic": 0.4,             # Activación células dendríticas
    "k_T_activation": 0.25,         # Activación linfocitos T CD4+
    "k_B_proliferation": 0.15,      # Proliferación linfocitos B
    "decay_IgG": 0.015,             # Vida media IgG bebé (~46 días)
    "maturation_factor": 0.002,     # Maduración gradual del sistema
    
    # Vacunación
    "doses": [60, 120, 360],        # Esquema 2+1 (2, 4, 12 meses)
    "adjuvant_boost": 1.5,          # Potenciación por adyuvantes
    
    # Correlatos de protección por serotipo (µg/ml)
    "thresholds": {
        "1": 0.35, "5": 0.23, 
        "6B": 0.10, "19A": 0.12
    }
}

def immune_response(y, t, params):
    IgG_maternal, dendritic, T_cells, B_cells, IgG_baby_1, IgG_baby_6B = y
    
    # Efecto blunting (supresión por IgG materna)
    blunting = 1 - (IgG_maternal/params["IgG_maternal_init"])
    
    # Maduración progresiva del sistema inmunitario
    maturation = 1 + params["maturation_factor"] * t
    
    # Señal de vacunación (con adyuvantes)
    vaccine_signal = 0
    for dose in params["doses"]:
        if 0 < t - dose < 7:
            vaccine_signal = params["adjuvant_boost"]
    
    # Ecuaciones diferenciales
    dIgG_maternal = -params["decay_maternal"] * IgG_maternal
    
    ddendritic = (params["k_dendritic"] * vaccine_signal 
                 - 0.1 * dendritic)  # Decaimiento células dendríticas
    
    dT = (params["k_T_activation"] * dendritic * maturation 
         - 0.05 * T_cells)
    
    dB = (blunting * params["k_B_proliferation"] * T_cells * vaccine_signal 
         - 0.07 * B_cells)
    
    # Producción de anticuerpos para diferentes serotipos
    dIgG1 = 0.3 * B_cells * (IgG_baby_1 < params["thresholds"]["1"]) - params["decay_IgG"] * IgG_baby_1
    
    dIgG6B = 0.15 * B_cells * (IgG_baby_6B < params["thresholds"]["6B"]) - params["decay_IgG"] * IgG_baby_6B
    
    return [dIgG_maternal, ddendritic, dT, dB, dIgG1, dIgG6B]

# Condiciones iniciales
y0 = [
    params["IgG_maternal_init"],  # IgG maternal
    0.0,   # Células dendríticas
    0.0,   # Linfocitos T
    0.0,   # Linfocitos B
    0.0,   # IgG serotipo 1
    0.0    # IgG serotipo 6B
]

t = np.linspace(0, 365, 365)

# Resolver el sistema
solution = odeint(immune_response, y0, t, args=(params,))
IgG_maternal = solution[:,0]
IgG_1 = solution[:,4]
IgG_6B = solution[:,5]

# Visualización
plt.figure(figsize=(12, 6))

# Curvas de anticuerpos
plt.plot(t, IgG_maternal, '--', label='IgG Materna')
plt.plot(t, IgG_1, label='IgG Bebé (Serotipo 1)')
plt.plot(t, IgG_6B, label='IgG Bebé (Serotipo 6B)')

# Umbrales de protección
for sero, thresh in params["thresholds"].items():
    plt.axhline(thresh, linestyle=':', alpha=0.7, 
                label=f'Umbral {sero} ({thresh} µg/ml)')

# Eventos de vacunación
for dose in params["doses"]:
    plt.axvline(dose, color='gray', linestyle='--', alpha=0.5)

plt.title('Respuesta Inmunológica Mejorada a PCV\n(Incluye Múltiples Serotipos y Componentes Inmunológicos)')
plt.xlabel('Días desde nacimiento')
plt.ylabel('Nivel de Anticuerpos (Unidades Arbitrarias)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()

# Verificación de protección
for sero, thresh in params["thresholds"].items():
    idx = np.argmax(solution[:,4 if sero=="1" else 5] >= thresh)
    if solution[:,4 if sero=="1" else 5][idx] >= thresh:
        print(f'Protección alcanzada para {sero} a los {t[idx]:.0f} días')
    else:
        print(f'Protección NO alcanzada para {sero}')