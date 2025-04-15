import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parámetros del modelo (ajustables según literatura)
params = {
    # Anticuerpos maternos (IgG)
    "IgG_maternal_init": 1.0,  # Nivel inicial de IgG materna (arbitrario)
    "decay_maternal": 0.03,    # Tasa de degradación de IgG materna (vida media ~21 días)

    # Vacuna PCV (respuesta del bebé)
    "k_activation": 0.5,       # Tasa de activación de linfocitos B por antígeno
    "k_proliferation": 0.2,    # Tasa de proliferación de células plasmáticas
    "decay_IgG": 0.01,         # Tasa de degradación de IgG del bebé
    "threshold_protection": 0.35,  # Correlato de protección (ej. 0.35 µg/mL para PCV20)

    # Tiempo
    "days": 365,               # Duración de la simulación (1 año)
    "doses": [30, 60, 180]     # Días de administración de dosis (ej. esquema 3+1)
}

def immune_response(y, t, params):
    IgG_maternal, IgG_baby = y  # Variables de estado

    # Degradación de IgG materna
    dIgG_maternal = -params["decay_maternal"] * IgG_maternal

    # Respuesta a la vacuna (activación linfocitos B -> producción de IgG)
    vaccine_stimulus = 0
    for dose_day in params["doses"]:
        if abs(t - dose_day) < 7:  # Estímulo de 7 días tras cada dosis
            vaccine_stimulus = params["k_activation"]

    dIgG_baby = (params["k_proliferation"] * vaccine_stimulus) - params["decay_IgG"] * IgG_baby

    return [dIgG_maternal, dIgG_baby]

# Condiciones iniciales
y0 = [params["IgG_maternal_init"], 0]  # Nivel inicial de IgG materna y 0 en el bebé
t = np.linspace(0, params["days"], params["days"])

# Resolver ecuaciones diferenciales
solution = odeint(immune_response, y0, t, args=(params,))
IgG_maternal = solution[:, 0]
IgG_baby = solution[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(t, IgG_maternal, label="IgG materna", linestyle="--")
plt.plot(t, IgG_baby, label="IgG del bebé (vacuna)", color="red")
plt.axhline(y=params["threshold_protection"], color="green", linestyle=":", label="Umbral de protección")

# Marcar dosis
for dose in params["doses"]:
    plt.axvline(x=dose, color="gray", linestyle=":", alpha=0.5)

plt.title("Respuesta Inmunológica del Bebé a la Vacuna PCV")
plt.xlabel("Días desde el nacimiento")
plt.ylabel("Nivel de anticuerpos (IgG)")
plt.legend()
plt.grid(True)
plt.show()

# Verificar si se alcanza el correlato de protección
protection_day = np.where(IgG_baby >= params["threshold_protection"])[0]
if len(protection_day) > 0:
    print(f"Protección alcanzada a los {protection_day[0]} días")
else:
    print("No se alcanzó el umbral de protección")