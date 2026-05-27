# Immune System Simulator (PCV)

Simulador computacional de respuesta inmune para vacunas antineumocócicas conjugadas (PCV), con enfoque en:

- dinámica de células B en centros germinales,
- producción de anticuerpos por serotipo,
- calibración automática contra datos clínicos,
- análisis visual y métricas de ajuste.

Este repositorio implementa una versión optimizada del modelo, con mejoras de rendimiento, trazabilidad y control de estocasticidad.

## Estado Actual del Proyecto

La versión actual incluye:

- motor de simulación separado en corrida simple (`simulate_once`) y corrida por lotes (`simulate_batch`),
- eliminación de `simpy` en el camino crítico de asignación de B-cells (ahora vectorizado con NumPy),
- recocido simulado rediseñado con:
  - optimización multi-parámetro,
  - evaluación por múltiples semillas,
  - caché de evaluaciones,
  - early stopping,
  - logs estructurados (`jsonl` y `csv`),
- métricas de ajuste compuestas:
  - Chamfer distance
  - error de centroides (landmark L2),
- factores de producción de anticuerpos por serotipo:
  - `plasma_production_factor`
  - `memory_production_factor`
- scripts utilitarios para:
  - corridas mock,
  - corridas de calibración corta,
  - benchmarks,
  - visualizaciones comparativas por serotipo,
  - inspección de métricas/logs.

## Estructura Principal

- `sim/main.py`
  - API principal del simulador.
  - `simulate_once(params, seed=None)`: una cohorte.
  - `simulate_batch(params, n_replicas, base_seed, seeds=None)`: varias cohortes.
- `sim/environment/immune_system.py`
  - Modelo del sistema inmune global.
  - Asignación vectorizada de células B a antígenos.
  - Actualización diaria de niveles de anticuerpos.
- `sim/environment/germinal_center.py`
  - Dinámica por centro germinal:
  - selección, diferenciación, mutación y expansión.
- `sim/processes/`
  - `selection.py`: selección Boltzmann.
  - `mutation.py`: hipermutación somática.
  - `affinity.py`: cálculo de afinidad.
  - `differentiation.py`: destino memoria/plasma.
- `sim/simulated_annealing.py`
  - Calibración del modelo con recocido simulado.
- `sim/config.py`
  - Parámetros del simulador (`SIMULATION_PARAMS`).
  - Parámetros de calibración (`CALIBRATION_PARAMS`).
- `sim/analysis/`
  - `benchmark.py`: pruebas de rendimiento/reproducibilidad.
  - `run_new_visuals.py`: generación de figuras comparativas.
- `images/new_model/`
  - Salidas gráficas y métricas exportadas.

## Modelo y Flujo de Simulación

### 1. Inicialización

Se crean:

- antígenos (serotipos PCV),
- pool inicial naive de células B,
- estructuras de memoria y plasma por serotipo.

### 2. Vacunación

En cada día del `vaccination_schedule`:

- se recluta un lote incremental de células naive,
- se asignan a antígenos usando matriz de afinidad vectorizada,
- se siembran centros germinales por antígeno.

### 3. Dinámica diaria

Por cada centro germinal:

- selección estocástica tipo Boltzmann,
- diferenciación de células B (memoria/plasma),
- muerte estocástica,
- mutación adaptativa (fuerza decreciente por ciclo y modulada por afinidad).

Luego:

- se actualizan niveles de anticuerpos por serotipo,
- decae fracción de células plasmáticas por `decay_factor`.

### 4. Salida

Para cada serotipo se devuelve un par `[pre, post]` por corrida.

Nota: para evitar fallos con ventanas cortas (`duration_days` bajos), la salida se fuerza a longitud consistente (2 puntos por serotipo).

## Parámetros por Serotipo

`plasma_production_factor` y `memory_production_factor` ahora pueden ser:

- `dict` por serotipo (recomendado),
- `list/array` (indexado por orden de serotipos internos),
- escalar (modo legacy).

Ejemplo recomendado:

```python
params["plasma_production_factor"] = {
    "1": 0.015, "5": 0.013, "6B": 0.012, "14": 0.009,
    "18C": 0.007, "19F": 0.011, "23F": 0.012,
}
```

## Calibración (Recocido Simulado)

Archivo: `sim/simulated_annealing.py`

Características:

- vector de optimización ampliado (incluye parámetros biológicos y factores por serotipo),
- evaluación robusta por múltiples réplicas (`n_replicas`),
- función objetivo compuesta:
  - `0.75 * Chamfer + 0.25 * centroid_l2`,
- ponderación por serotipo (`serotype_weights`),
- caché de evaluaciones por hash de parámetros,
- detención temprana (`patience`).

Logs generados:

- `misc/annealing_log.jsonl` (detalle por iteración),
- `misc/annealing_log.csv` (resumen tabular).

## Scripts de Ejecución

### 1. Corrida mock simple

Archivo: `sim/run_single_mock.py`

- estocástico (seed libre):

```bash
python3 sim/run_single_mock.py
```

- reproducible:

```bash
python3 sim/run_single_mock.py --seed 123
```

### 2. Recocido corto (debug rápido)

Archivo: `sim/run_short_annealing.py`

```bash
python3 sim/run_short_annealing.py
```

### 3. Benchmark/reproducibilidad

Archivo: `sim/analysis/benchmark.py`

```bash
python3 sim/analysis/benchmark.py
```

### 4. Visualización comparativa (nuevo modelo vs datos)

Archivo: `sim/analysis/run_new_visuals.py`

```bash
python3 sim/analysis/run_new_visuals.py --replicas 24 --duration-days 140
```

Salida:

- `images/new_model/serotypes_grid.png`
- `images/new_model/sp*.png`
- `images/new_model/metrics.json`

### 5. Inspección de scores por serotipo desde log

```bash
python3 sim/show_serotype_scores.py
python3 sim/show_serotype_scores.py --iter 3
```

### 6. Inspección de factores plasma/memoria por serotipo

```bash
python3 sim/show_serotype_factors.py
python3 sim/show_serotype_factors.py --iter 3
```

## Datos y Métricas

Datasets relevantes:

- `sim/data/train_VCN7-Tf_fit.csv`
- `sim/data/test_VCN7-Tf_fit.csv`

Métricas usadas:

- **Chamfer distance** entre nube real y simulada por serotipo.
- **Landmark L2** entre centroides real/simulado.
- **Objetivo global**: promedio ponderado entre serotipos.

## Interpretación de Resultados

- `best_energy` menor implica mejor ajuste global según la métrica compuesta.
- Revisar siempre:
  - score por serotipo,
  - Chamfer por serotipo,
  - centroid_l2 por serotipo,
  - dispersión entre seeds.

Si un serotipo (p. ej. `18C`) queda muy desviado:

- revisar horizonte (`duration_days`) vs esquema de dosis,
- revisar pesos por serotipo en calibración,
- revisar factores plasma/memoria específicos de ese serotipo.

## Requisitos y Ejecución

Entorno esperado:

- Python 3
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

Ejecutar desde la raíz del repo:

```bash
cd /home/amanda/simulaciones/Immune-system-simulator
```

