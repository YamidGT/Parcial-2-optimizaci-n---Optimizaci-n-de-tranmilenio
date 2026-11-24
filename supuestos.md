# Supuestos del modelo (Escenario 1: TransMilenio, corredor NQS con lluvia)

Cada fila dice de dónde sale el valor o por qué se asume. Las filas marcadas **POR CONFIRMAR** son parámetros sin dato: el grupo debe fijar el valor y escribir la justificación antes de entregar (un dato sin explicación cuenta como error de modelado).

## Alcance y red
| # | Supuesto | Justificación / fuente |
|---|---|---|
| S1 | La red es un corredor lineal de 27 estaciones: La Castellana → Ricaurte (zona E, 11) y Comuneros → San Mateo (zona G, 16). Longitud del eje: 21,6 km. | PDF del parcial pide las líneas E y G. Orden norte→sur por coordenadas de `estaciones_troncales.geojson`. |
| S2 | Tygua–San José y Guatoque–Veraguas quedan fuera del eje. | Están en el ramal de la Carrera 10 (trazado TZ007), no sobre la NQS. Verificar en el mapa oficial vigente. |
| S3 | SENA se excluye de la red base. | Figura cerrada en el Informe #99 (trimestre II de 2026) y no aparece en las salidas de noviembre de 2024. |
| S4 | "Ampliación San Mateo" (57503) se suma a San Mateo (07503). | Es la misma terminal ampliada. Se suman sus salidas. |
| S5 | Distancia entre estaciones = distancia haversine entre coordenadas. | Aproxima el trazado; el error es pequeño en tramos casi rectos. |

## Demanda
| # | Supuesto | Justificación / fuente |
|---|---|---|
| D1 | Demanda por estación = salidas promedio de día hábil por intervalo de 15 min (noviembre de 2024). | `salidas_troncal_nov2024.csv`. Es el dato oficial más cercano a afluencia. Es de hace dos años; se declara. |
| D2 | Día hábil: se excluyen sábados, domingos y los festivos del 4 y 11 de noviembre de 2024 (19 días). | Calendario de festivos de Colombia. |
| D3 | Matriz origen-destino de la hora pico de la mañana (07:15–08:15). Destinos = salidas en ese pico. Orígenes = salidas en el pico de la tarde (18:30–19:30), por simetría del viaje de regreso. Distribución gravitacional doblemente acotada (Furness). Resultado: 33.373 viajes/h. | Las salidas muestran dirección: pico de la mañana al norte y centro (07:00 a 07:45), pico de la tarde al sur (17:45 a 18:45). |
| D4 | Hora pico de cada estación = mayor suma de 4 intervalos consecutivos. | Definido en `preparar_datos.py`. |

## Operación y costos (sin combustible, según el PDF)
| # | Supuesto | Justificación / fuente |
|---|---|---|
| O1 | Velocidad comercial base 25 km/h. Rango de prueba 21 a 27 km/h. | 25,2 km/h en 2024 (22,2 en pico); 21 km/h en 2023. Fuentes con metodologías distintas. |
| O2 | Capacidad del bus: articulado 160, biarticulado 250 pasajeros. En el modelo, Q = 220. | Cifras oficiales de TransMilenio. Q es el promedio ponderado de la flota troncal: 1.320 biarticulados y 610 articulados (Informe #99). |
| O3 | **PROVISIONAL.** Tiempo de parada τ = 0,5 min; regulación en terminal = 5 min. | Valor típico de BRT entre 20 y 40 s; falta una referencia. |
| O4 | **PROVISIONAL.** Costo por bus-hora c = 550.000 COP, sin combustible. | Cadena: costo técnico ≈ 2 × tarifa (el Distrito subsidia ~50 %) ≈ 7.100 COP/pax × IPK troncal ≈ 4,5 pax/km (Informe #99) × 25 km/h ≈ 800.000 COP/h; menos ~30 % de combustible (supuesto). Probar un rango de 400.000 a 800.000. |
| O5 | **PROVISIONAL.** Valor del tiempo v = 9.600 COP/h. | Salario mínimo 2026 de $1.750.905 (Decreto 1469 de 2025) / ~182 h al mes con jornada de 42 h. Verificar si el Decreto 0159 de 2026 cambió el monto. |
| O6 | Flota disponible N = 414 buses en la NQS. | 2.203 buses troncales × 18,8 % (participación de la NQS en las salidas troncales de día hábil, calculada en `nqs_modelo.participacion_nqs`). El modelo base usa 168, así que el supuesto no es crítico. |
| O8 | **PROVISIONAL.** Disuasión del modelo gravitacional β = 0,10 por km. | Da una distancia media de viaje de 9,2 km en un eje de 21,6 km. Probar 0,05 y 0,20. |
| O7 | Tarifa $3.550 (no se optimiza). | Decreto 007 de 2026. Fijada por decreto. |
| O9 | **PROVISIONAL.** Hacinamiento (Bloque B): el tiempo en vehículo se pondera por 1 + κ·ρ², con ρ = ocupación del tramo crítico y κ = 1 (a bus lleno el tiempo vale el doble). | La literatura de valoración del hacinamiento reporta multiplicadores de 1,5 a 2 para viajar de pie en alta densidad (ej. Wardman y Whelan, 2011). Falta verificar la cita y probar κ = 0,5 y 1,5. Sin este término (κ = 0), la capacidad fija todas las frecuencias y el Bloque B no decide nada. |

## Estrés (lluvia)
| # | Supuesto | Justificación / fuente |
|---|---|---|
| E1 | Lluvia: tiempos de recorrido +20 % y +40 %. | Niveles indicados en el PDF del parcial. |
| E2 | Estación cerrada: una estación de alta demanda (Ricaurte o Paloquemao). | SENA ya es el cierre real de la base; el estrés cierra otra para ver cambios en la decisión. |
| E3 | **PROVISIONAL.** Con una estación cerrada, sus viajes se reparten mitad y mitad entre la estación anterior y la siguiente. La caminata se cuenta a 4,5 km/h sobre la distancia media a esas vecinas. Los buses pasan sin parar. | Velocidad típica de caminata urbana. No se modela el transbordo con otras troncales (en Ricaurte, con Américas), así que el impacto del cierre de Ricaurte queda subestimado. |

## Fuera de alcance (según el PDF)
Combustible y densidad poblacional no se modelan.
