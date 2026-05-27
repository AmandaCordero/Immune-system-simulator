# Resumen Punto a Punto de Cambios en la Simulación

## 1. Reestructuración del motor de simulación

- Se reemplazó el flujo antiguo centrado en `la()`/`sim()` por funciones explícitas:
  - `simulate_once(params, seed=None)` para una cohorte.
  - `simulate_batch(params, n_replicas, base_seed, seeds=None)` para múltiples réplicas.
- Se mantuvo compatibilidad retroactiva con notebooks/código previo:
  - `la(params)` ahora delega a `simulate_once`.
  - `sim(params, times)` ahora ejecuta varias corridas con `simulate_once`.

## 2. Reproducibilidad y estocasticidad controlada

- Se unificó el manejo de aleatoriedad:
  - Se fijan `random.seed` y `numpy.random.seed` cuando se pasa `seed`.
- Efecto práctico:
  - Misma `seed` => mismos resultados.
  - `seed` distinta o `seed=None` => variación estocástica real.

## 3. Eliminación del cuello de botella con `simpy`

- Se eliminó la asignación de células B por procesos/locks de `simpy`.
- Se implementó asignación vectorizada de afinidad en `ImmuneSystem`:
  - Se construye matriz de afinidad `N_bcells x N_antigens`.
  - Se asigna por `argmax` con umbral (`THRESHOLD`).
- Resultado:
  - Menor sobrecarga de concurrencia artificial.
  - Mejor rendimiento en la fase de sembrado de centros germinales.

## 4. Cambio en vacunación y estado inmune

- Antes: se reinicializaba por completo el pool naive en cada vacunación.
- Ahora:
  - Se mantiene estado inmune entre dosis.
  - Se reclutan naive incrementales por evento (`naive_recruitment_per_dose`).
- Esto evita reinicios agresivos y conserva dinámica acumulativa (memoria/plasma).

## 5. Rediseño de dinámica en centros germinales (GC)

- Se introdujo selección adaptativa por ciclo:
  - Temperatura efectiva decreciente (`selection_schedule`).
  - Fracción de supervivientes con crecimiento gradual.
  - Supervivencia mínima garantizada (evita `max_survivors=0` temprano).
- Se rediseñó mutación:
  - Fuerza base decreciente por ciclo (`mutation_schedule`).
  - Fuerza modulada por afinidad de cada célula (más exploración en baja afinidad).
- Se añadió límite de tamaño de pool GC (`max_gc_pool`) para evitar crecimiento descontrolado.

## 6. Robustez de selección estocástica

- `boltzmann_selection` ahora maneja mejor casos borde:
  - Lista vacía de células.
  - `max_survivors <= 0`.
  - Temperatura muy pequeña (protección numérica).

## 7. Configuración extendida del modelo

- Se amplió `SIMULATION_PARAMS` con nuevos controles:
  - `naive_pool_size`
  - `naive_recruitment_per_dose`
  - `min_selection_fraction`, `max_selection_fraction`
  - `selection_growth_rate`
  - `max_gc_pool`
  - `selection_schedule`
  - `mutation_schedule`
- Se añadió bloque `CALIBRATION_PARAMS`:
  - Rangos de parámetros (`param_bounds`) para optimización.
  - Réplicas y semilla base.
  - Pesos por serotipo.
  - Pesos de métrica compuesta.
  - Configuración del recocido (temperaturas, iteraciones, paciencia).

## 8. Rediseño completo de calibración (recocido simulado)

- Se reescribió `sim/simulated_annealing.py` para:
  - Optimizar un vector amplio de parámetros (no solo 2 factores).
  - Evaluar cada candidato con múltiples réplicas (`simulate_batch`).
  - Aplicar caché de evaluaciones por hash de parámetros.
  - Usar early stopping por `patience`.
- Se añadieron logs estructurados:
  - `misc/annealing_log.jsonl` (detalle completo por iteración).
  - `misc/annealing_log.csv` (resumen tabular por iteración).

## 9. Nueva métrica objetivo de optimización

- Antes: se usaba solo Chamfer por serotipo aislado.
- Ahora: objetivo compuesto por serotipo:
  - `score = 0.75 * chamfer + 0.25 * landmark_l2`
  - Promedio ponderado entre serotipos (`serotype_weights`), con mayor peso en `14`.
- Se conservó cálculo de Chamfer y se añadió término de distancia entre centroides.

## 10. Corrección para horizontes cortos (`duration_days` bajos)

- Se corrigió `simulate_once` para garantizar siempre 2 puntos de salida por serotipo.
- Esto evita errores (`IndexError`) cuando `duration_days` no alcanza todos los checkpoints.
- Impacto: el recocido funciona también con configuraciones rápidas (ej. `duration_days=140`).

## 11. Scripts nuevos de ejecución y soporte

- `sim/run_short_annealing.py`
  - Ejecuta recocido corto configurable para pruebas rápidas.
  - Incluye logs de progreso (inicio, configuración, tiempo total, resultados).
- `sim/run_single_mock.py`
  - Corre una sola simulación mock.
  - Acepta `--seed` opcional para modo reproducible o estocástico.
- `sim/show_serotype_scores.py`
  - Extrae y muestra scores por serotipo desde `annealing_log.jsonl`.
  - Permite consultar última iteración o una iteración específica.
- `sim/analysis/benchmark.py`
  - Pruebas de benchmark y chequeo de reproducibilidad.
- `sim/analysis/run_new_visuals.py`
  - Genera visuales automáticos comparando real vs simulado.
  - Exporta imágenes por serotipo y grid general en `images/new_model/`.
  - Exporta métricas por serotipo en `images/new_model/metrics.json`.

## 12. Salidas visuales nuevas

- Se generan automáticamente:
  - `images/new_model/serotypes_grid.png`
  - `images/new_model/sp1.png`, `sp14.png`, etc.
  - `images/new_model/metrics.json`
- Cada plot incluye:
  - Nube de puntos real (clínica).
  - Nube simulada.
  - Centroides real/simulado.
  - Chamfer y distancia entre centroides (C-L2).

## 13. Resultado técnico global de esta iteración

- Se modernizó la arquitectura para separar simulación simple, batch y calibración.
- Se eliminó una fuente principal de lentitud (`simpy` en asignación).
- Se mejoró control de variabilidad y reproducibilidad.
- Se amplió y estabilizó la optimización con recocido.
- Se añadió trazabilidad y herramientas de análisis visual/numérico para comparar resultados.
