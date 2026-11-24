# Parcial 2 – Optimización del corredor NQS de TransMilenio con lluvia

Parcial 2 aplicado de Optimización (Universidad Nacional de Colombia, 2026), Escenario 1. ¿Qué servicios (corriente y expresos) operar en la troncal NQS (líneas E y G), con qué paradas y frecuencias, para minimizar el costo de operación más el tiempo de los usuarios en la hora pico? ¿Y cómo cambia esa decisión con lluvia o con una estación cerrada?

## Modelos
| Bloque | Modelo | Herramienta |
|---|---|---|
| A | MILP de selección de servicios, paradas y frecuencias; ruta más corta con Dijkstra y A\* | PuLP + HiGHS |
| B | Frecuencia óptima no lineal con espera y hacinamiento; condiciones KKT | KKT a mano y `scipy.optimize.minimize` (SLSQP) |
| C | Simulated Annealing contra el MILP exacto; frontera de Pareto por ε-restricción | Python + PuLP |
| Sensibilidad | Lluvia +20 % y +40 %, estación cerrada (Ricaurte, Paloquemao), Pareto con lluvia | – |

## Estructura
- `preparar_datos.py`: arma la red NQS (27 estaciones) y la demanda a partir de los datos abiertos.
- `nqs_modelo.py`: todos los modelos (matriz OD, MILP, Dijkstra/A\*, KKT, SA, Pareto, sensibilidad).
- `parcial2_nqs.ipynb`: notebook que corre de principio a fin, con la formulación de cada bloque.
- `supuestos.md`: supuestos con su justificación y fuente. Los marcados PROVISIONAL deben validarse.
- `datos/`: datos abiertos descargados (estaciones, trazados, rutas, salidas por estación, flota, informe de oferta y demanda), más la red procesada.
- `documento/`: avance del Bloque A y resumen del proyecto (Word).

## Cómo ejecutar
**En Google Colab:** abrir [el notebook en Colab](https://colab.research.google.com/github/YamidGT/Parcial-2-optimizaci-n---Optimizaci-n-de-tranmilenio/blob/main/parcial2_nqs.ipynb) y usar *Entorno de ejecución → Ejecutar todo*. La primera celda de código clona este repositorio, instala las dependencias y prepara los datos. Para solo ver los resultados no hace falta ejecutar: el notebook está guardado con sus salidas.

**En local:**
```bash
pip install "pulp>=4.1" highspy numpy scipy pandas matplotlib networkx jupyter
python preparar_datos.py
jupyter nbconvert --to notebook --execute --inplace parcial2_nqs.ipynb
```
Probado con Python 3.13, PuLP 4.1 (solver HiGHS incluido), NumPy 2.5, SciPy 1.18, pandas 3.0, matplotlib 3.10 y networkx 3.7. La ejecución completa tarda unos 12 minutos, sobre todo por el MILP de la instancia grande del Bloque C (límite de 180 s).

## Fuentes de datos
- Portal de Datos Abiertos de Bogotá: salidas del sistema troncal por franja horaria, flota vinculada.
- Portal de datos abiertos de TransMilenio: estaciones troncales, trazados troncales, servicios.
- TRANSMILENIO S.A., Estadísticas de oferta y demanda del SITP, trimestre II de 2026 (Informe #99).
- Decreto 007 de 2026 (tarifa).
