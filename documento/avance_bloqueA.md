# Avance – 24 de noviembre de 2026

**Escenario elegido:** 1. TransMilenio y SITP en el corredor de la NQS con lluvia.
**Pregunta de planeación:** ¿qué servicios (corrientes y expresos) operar en la NQS, con qué paradas, para minimizar el costo de operación más el tiempo de los usuarios, y cómo cambia esa decisión cuando llueve?

## 1. Alcance
- **Red:** 27 estaciones de la troncal NQS. Zona E: La Castellana → Ricaurte (11). Zona G: Comuneros → San Mateo (16). Eje de 21,6 km.
- **Excluidas:** SENA (cerrada según el Informe #99) y el ramal de la Carrera 10 (Tygua y Guatoque). Detalle en `supuestos.md`.
- **Plus (opcional, por decidir):** troncal Américas por Ricaurte o Calle 26 por CAD / Av. El Dorado, como estaciones de intercambio.
- **Fuera de alcance:** combustible y densidad poblacional. Costo por bus-hora sin combustible.

## 2. Datos
| Dato | Fuente |
|---|---|
| Estaciones, coordenadas, capacidad | Portal de datos abiertos de TransMilenio (estaciones troncales, oct 2026) |
| Demanda por estación cada 15 min | Salidas del sistema troncal por franja horaria, nov 2024 (datos abiertos Bogotá) |
| Velocidad, capacidad de buses, flota, tarifa | Informe #99, fuentes oficiales y prensa (ver `supuestos.md`) |

**Lo que muestran las salidas (día hábil promedio):**
- Portal Sur (65.622), San Mateo (49.135) y Terreros (31.237) concentran la demanda del sur.
- En la zona E, Av. Chile (18.715), U. Nacional (17.132) y Paloquemao (16.002) lideran.
- **El pico tiene dirección:** en la zona E, las salidas tienen pico de 07:00 a 07:45 (gente que llega a trabajar o estudiar). En la zona G, el pico es de 17:45 a 18:45 (regreso a casa). Esto justifica servicios expresos en sentido del pico.

## 3. Bloque A – Modelo núcleo (programación entera)

### 3.1 Selección de servicios y paradas (MILP)
**Conjuntos**
- \(N = \{0, \dots, 26\}\): estaciones en orden norte → sur.
- \(S\): servicios candidatos. Cada \(s\) tiene origen \(a_s\), destino \(b_s\), patrón de paradas \(P_s \subseteq N\) y sentido. Incluye el servicio corriente (para en todas) y una familia de expresos generados de antemano.
- \(K\): pares origen-destino \((o, d)\) con demanda positiva.
- \(L\): tramos entre estaciones consecutivas \((i, i+1)\).

**Parámetros**
- \(D_{od}\): pasajeros por hora en la hora pico de \(o\) a \(d\) (supuesto D3).
- \(t_{i,i+1} = \frac{\text{km}_{i,i+1}}{v_c}(1 + r)\): tiempo de recorrido del tramo, donde \(r\) es el factor de lluvia (0, 0,2 o 0,4).
- \(\tau\): tiempo de parada. \(T_s\): tiempo de ciclo del servicio \(s\) (ida y vuelta, recorrido más paradas).
- \(f_s\): frecuencia fija en el Bloque A (buses/h); se optimiza en el Bloque B.
- \(Q\): capacidad por bus. \(c\): costo por bus-hora. \(v\): valor del tiempo. \(N^{max}\): flota.
- \(\theta_{od,s}\): tiempo en vehículo de \(o\) a \(d\) en \(s\). Espera \(w_s = 1/(2 f_s)\).
- \(a_{od,s} = 1\) si \(s\) para en \(o\) y en \(d\) y va en ese sentido.

**Variables**
- \(y_s \in \{0, 1\}\): operar el servicio \(s\).
- \(x_{od,s} \ge 0\): pasajeros por hora de \((o, d)\) asignados a \(s\).

**Función objetivo**
\[
\min \; \underbrace{\sum_{s} c\, f_s\, T_s\, y_s}_{\text{costo de operación (buses en circulación)}} \;+\; v \underbrace{\sum_{(o,d)} \sum_{s} \left(\theta_{od,s} + w_s\right) x_{od,s}}_{\text{tiempo total de usuarios}}
\]

**Restricciones**
1. Atender la demanda: \(\sum_s x_{od,s} = D_{od}\) para todo \((o, d)\).
2. Solo servicios que paran en ambos extremos y operan: \(x_{od,s} \le a_{od,s}\, D_{od}\, y_s\).
3. Capacidad por tramo y servicio: \(\sum_{(o,d)\,\text{que cruzan}\,\ell} x_{od,s} \le Q\, f_s\, y_s\) para todo \(\ell \in L\), \(s \in S\).
4. Flota: \(\sum_s f_s\, T_s\, y_s \le N^{max}\).
5. El servicio corriente siempre opera: \(y_{\text{corriente}} = 1\). Garantiza cobertura de todas las estaciones.

**Interpretación pedida**
- **Precio sombra de la flota (restricción 4)** en la relajación LP: cuánto baja el costo total con un bus más.
- **Efecto de fijar una variable:** fijar un expreso en \(y_s = 1\) (o retirarlo) y comparar el objetivo.

**Solver:** PuLP con CBC (u OR-Tools).

### 3.2 Ruta más corta (Dijkstra y A\*), apoyo
- Grafo dirigido: nodos = estaciones. Arcos entre consecutivas con peso \(t_{i,i+1}(1+r)+\tau\).
- **Heurística admisible de A\*:** \(h(i) = \frac{|\text{km}_i - \text{km}_{\text{destino}}|}{v_{\max}}\). Nunca sobreestima porque usa la distancia sobre el eje y la velocidad máxima.
- **Comparación:** nodos expandidos, Dijkstra vs A\*, con y sin lluvia, por ejemplo de Portal Sur a U. Nacional.

## 4. Conexión con los otros bloques (adelanto)
- **Bloque B (U7):** con los servicios elegidos, frecuencia óptima \(\min \sum_s c f_s T_s + \sum_s v D_s/(2 f_s)\) sujeto a \(\sum_s f_s T_s \le N^{max}\). Se resuelve con KKT y con `scipy.optimize.minimize`. \(\lambda\) es el valor de un bus adicional.
- **Bloque C (U8):** algoritmo genético o simulated annealing sobre los patrones de paradas de los expresos, contra el MILP exacto en la misma instancia. Frontera de Pareto entre costo de operación y tiempo de usuarios.
- **Sensibilidad:** lluvia +20 %, +40 % y cierre de una estación (Ricaurte o Paloquemao).

## 5. Pendientes antes de resolver
Fijar y justificar \(c\), \(v\), \(\tau\), \(N^{max}\) y la matriz origen-destino \(D_{od}\). Están marcados como POR CONFIRMAR en `supuestos.md`.
