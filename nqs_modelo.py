"""Modelos del Parcial 2 (Escenario 1: corredor NQS con lluvia).

Bloque A:
  - construir_od(): matriz origen-destino de la hora pico (modelo gravitacional).
  - generar_servicios(): servicios candidatos (corriente + expresos) con frecuencias.
  - resolver_milp(): seleccion de servicios, paradas y frecuencias (PuLP / HiGHS).
  - dijkstra() y a_estrella(): ruta mas corta en la red de servicios.

Todos los parametros estan en PARAMS. Los marcados PROVISIONAL deben quedar
justificados en supuestos.md antes de la entrega.
Requiere haber corrido preparar_datos.py (genera datos/nqs_red.csv).
"""
from __future__ import annotations

import csv
import heapq
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pulp

BASE = Path(__file__).resolve().parent / "datos"

PARAMS = {
    "v_c": 25.0,            # km/h velocidad comercial (O1)
    "tau_min": 0.5,         # min de parada por estacion (O3, PROVISIONAL)
    "layover_min": 5.0,     # min de regulacion en cada terminal (PROVISIONAL)
    "Q": 220,               # pax/bus: promedio ponderado articulados y biarticulados (O2)
    "c_bus_h": 550_000,     # COP por bus-hora sin combustible (O4, PROVISIONAL)
    "v_tiempo_h": 9_600,    # COP por hora de usuario (O5, PROVISIONAL)
    "beta": 0.10,           # 1/km, disuasion del modelo gravitacional (D3, PROVISIONAL)
    "flota_troncal": 2203,  # buses troncales (Informe #99)
    "frecuencias": (6, 10, 15, 20, 30, 40),  # buses/h por servicio
}

DIAS_NO_HABILES = {2, 3, 4, 9, 10, 11, 16, 17, 23, 24, 30}


# ---------------------------------------------------------------- datos
@dataclass
class Red:
    codigos: list[str]
    nombres: list[str]
    km: np.ndarray                 # km acumulado sobre el eje
    salidas_dia: np.ndarray
    demanda_15: dict[str, np.ndarray]   # intervalo -> salidas por estacion

    @property
    def n(self) -> int:
        return len(self.codigos)


def cargar_red() -> Red:
    with open(BASE / "nqs_red.csv", encoding="utf-8") as fh:
        filas = list(csv.DictReader(fh))
    codigos = [f["codigo"] for f in filas]
    with open(BASE / "nqs_demanda_15min.csv", encoding="utf-8") as fh:
        lector = csv.reader(fh)
        cab = next(lector)
        assert cab[1:] == codigos, "el orden de estaciones no coincide"
        dem = {r[0]: np.array([float(x) for x in r[1:]]) for r in lector}
    return Red(
        codigos=codigos,
        nombres=[f["nombre"] for f in filas],
        km=np.array([float(f["km_acum"]) for f in filas]),
        salidas_dia=np.array([float(f["salidas_dia_habil"]) for f in filas]),
        demanda_15=dem,
    )


def participacion_nqs(red: Red) -> float:
    """Fraccion de las salidas troncales (dia habil) que ocurren en la NQS."""
    with open(BASE / "salidas" / "salidas_troncal_nov2024.csv", encoding="cp1252", newline="") as fh:
        filas = list(csv.reader(fh, delimiter=";"))
    cols = [i for i, d in zip(range(5, 36), range(1, 32)) if d not in DIAS_NO_HABILES]
    total = 0.0
    for r in filas[1:]:
        if len(r) < 36 or not r[1]:
            continue
        for i in cols:
            v = r[i].strip()
            if v.replace(".", "", 1).isdigit():
                total += float(v)
    total /= len(cols)
    return float(red.salidas_dia.sum() / total)


def ventana_pico(red: Red, manana: bool = True) -> list[str]:
    """Cuatro intervalos consecutivos con mas salidas en el corredor."""
    horas = sorted(red.demanda_15)
    tot = [red.demanda_15[h].sum() for h in horas]
    mejor, k_mejor = -1.0, 0
    for k in range(len(horas) - 3):
        if (horas[k] < "12:00") != manana:
            continue
        s = sum(tot[k:k + 4])
        if s > mejor:
            mejor, k_mejor = s, k
    return horas[k_mejor:k_mejor + 4]


def construir_od(red: Red, beta: float | None = None, iteraciones: int = 200) -> np.ndarray:
    """Matriz OD de la hora pico de la manana (pax/h).

    Destinos D_j: salidas en la hora pico de la manana.
    Origenes O_i: salidas en la hora pico de la tarde (simetria del viaje de regreso),
    escaladas para que sum(O) = sum(D). Distribucion doblemente acotada (Furness)
    con disuasion exp(-beta * distancia).
    """
    beta = PARAMS["beta"] if beta is None else beta
    am, pm = ventana_pico(red, True), ventana_pico(red, False)
    D = sum(red.demanda_15[h] for h in am)
    O = sum(red.demanda_15[h] for h in pm)
    O = O * D.sum() / O.sum()
    dist = np.abs(red.km[:, None] - red.km[None, :])
    F = np.exp(-beta * dist)
    np.fill_diagonal(F, 0.0)
    T = F * O[:, None]
    for _ in range(iteraciones):
        T *= (D / np.maximum(T.sum(axis=0), 1e-12))[None, :]
        T *= (O / np.maximum(T.sum(axis=1), 1e-12))[:, None]
    return T


# ------------------------------------------------------------ servicios
@dataclass
class Servicio:
    nombre: str
    paradas: tuple[int, ...]        # indices de estaciones, en orden
    f: float                        # buses/h
    patron: str                     # identificador del patron (sin frecuencia)
    corriente: bool = False
    T_ciclo_h: float = 0.0
    buses: int = 0
    pos: dict[int, int] = field(default_factory=dict)


def tiempo_tramo_h(red: Red, i: int, j: int, lluvia: float) -> float:
    return abs(red.km[j] - red.km[i]) / PARAMS["v_c"] * (1.0 + lluvia)


def tiempo_en_vehiculo_h(red: Red, s: Servicio, o: int, d: int, lluvia: float) -> float:
    """Recorrido de o a d mas las paradas intermedias del servicio."""
    intermedias = abs(s.pos[d] - s.pos[o]) - 1
    return tiempo_tramo_h(red, o, d, lluvia) + intermedias * PARAMS["tau_min"] / 60.0


def _completar(red: Red, s: Servicio, lluvia: float) -> Servicio:
    a, b = s.paradas[0], s.paradas[-1]
    ida = tiempo_tramo_h(red, a, b, lluvia) + (len(s.paradas) - 1) * PARAMS["tau_min"] / 60.0
    s.T_ciclo_h = 2 * ida + 2 * PARAMS["layover_min"] / 60.0
    s.buses = math.ceil(s.f * s.T_ciclo_h)
    s.pos = {e: k for k, e in enumerate(s.paradas)}
    return s


TERMINALES_NORTE = (0, 6, 10)    # La Castellana, U. Nacional, Ricaurte
TERMINALES_SUR = (21, 26)         # Portal Sur, San Mateo
PARADAS_K = (2, 4, 6, 8)


def patrones_candidatos(red: Red, od: np.ndarray, cerradas: frozenset[int] = frozenset(),
                        norte=TERMINALES_NORTE, sur=TERMINALES_SUR, ks=PARADAS_K):
    """Corriente + expresos entre terminales con las k estaciones de mas demanda."""
    abiertas = [i for i in range(red.n) if i not in cerradas]
    pats = {"Corriente": tuple(abiertas)}
    flujo = od.sum(axis=0) + od.sum(axis=1)          # pax que suben o bajan
    for a in norte:
        for b in sur:
            if a in cerradas or b in cerradas:
                continue
            medio = [i for i in range(a + 1, b) if i not in cerradas]
            orden = sorted(medio, key=lambda i: -flujo[i])
            for k in ks:
                paradas = tuple(sorted([a, b] + orden[:k]))
                nombre = f"E{red.codigos[a][-3:]}-{red.codigos[b][-3:]}-k{k}"
                pats[nombre] = paradas
    return pats


def generar_servicios(red: Red, od: np.ndarray, lluvia: float = 0.0,
                      cerradas: frozenset[int] = frozenset(), **kw) -> list[Servicio]:
    """kw: norte, sur, ks (para instancias mas grandes en el Bloque C)."""
    servicios = []
    for nombre, paradas in patrones_candidatos(red, od, cerradas, **kw).items():
        for f in PARAMS["frecuencias"]:
            s = Servicio(f"{nombre}@{f}", paradas, float(f), nombre, nombre == "Corriente")
            servicios.append(_completar(red, s, lluvia))
    return servicios


# ----------------------------------------------------------------- MILP
@dataclass
class ResultadoMILP:
    estado: str
    objetivo: float
    costo_operacion: float
    costo_usuarios: float
    elegidos: list[Servicio]
    buses: int
    flota: int
    segundos: float
    horas_usuario: float
    precio_sombra_flota: float | None = None
    detalle_carga: dict = field(default_factory=dict)
    duales_demanda: dict = field(default_factory=dict)   # (o,d) -> COP por pax adicional
    flujos: dict = field(default_factory=dict)   # servicio -> pasajeros, carga_max, horas_vehiculo


def resolver_milp(red: Red, od: np.ndarray, flota: int, lluvia: float = 0.0,
                  cerradas: frozenset[int] = frozenset(), relajar: bool = False,
                  fijar: dict[str, int] | None = None,
                  servicios: list[Servicio] | None = None,
                  peso_operacion: float = 1.0, peso_usuarios: float = 1.0,
                  limite_operacion: float | None = None,
                  limite_seg: int = 120) -> ResultadoMILP:
    """Seleccion de servicios, paradas y frecuencia (formulacion en el notebook).

    relajar=True resuelve la relajacion LP (y en [0,1]) para obtener el precio
    sombra de la restriccion de flota. fijar={nombre_servicio: 0/1} fija variables.
    """
    t0 = time.perf_counter()
    if servicios is None:
        servicios = generar_servicios(red, od, lluvia, cerradas)
    c, v, Q = PARAMS["c_bus_h"], PARAMS["v_tiempo_h"], PARAMS["Q"]

    # pares OD con demanda (las estaciones cerradas no generan ni atraen viajes)
    pares = [(o, d) for o in range(red.n) for d in range(red.n)
             if o != d and od[o, d] > 1e-6 and o not in cerradas and d not in cerradas]

    m = pulp.LpProblem("seleccion_servicios_NQS", pulp.LpMinimize)
    cat = pulp.LpContinuous if relajar else pulp.LpBinary
    # PuLP >= 4 crea variables desde el problema (m.add_variable)
    y = {s.nombre: m.add_variable(f"y_{k}", lowBound=0, upBound=1, cat=cat)
         for k, s in enumerate(servicios)}
    sirve = {}   # (o,d) -> lista de servicios que paran en ambas
    for s in servicios:
        for (o, d) in pares:
            if o in s.pos and d in s.pos:
                sirve.setdefault((o, d), []).append(s)
    idx = {s.nombre: k for k, s in enumerate(servicios)}
    x = {(o, d, s.nombre): m.add_variable(f"x_{o}_{d}_{idx[s.nombre]}", lowBound=0)
         for (o, d), lst in sirve.items() for s in lst}

    costo_op = pulp.lpSum(c * s.buses * y[s.nombre] for s in servicios)
    costo_us = pulp.lpSum(
        v * (tiempo_en_vehiculo_h(red, s, o, d, lluvia) + 1.0 / (2.0 * s.f)) * x[o, d, s.nombre]
        for (o, d), lst in sirve.items() for s in lst)
    m += peso_operacion * costo_op + peso_usuarios * costo_us

    # 1. demanda atendida
    for (o, d) in pares:
        m += pulp.lpSum(x[o, d, s.nombre] for s in sirve.get((o, d), [])) == od[o, d], f"dem_{o}_{d}"
    # 2. solo servicios operados
    for (o, d), lst in sirve.items():
        for s in lst:
            m += x[o, d, s.nombre] <= od[o, d] * y[s.nombre]
    # 3. capacidad por tramo, sentido y servicio
    for s in servicios:
        for sentido in (1, -1):
            for k in range(len(s.paradas) - 1):
                i, j = s.paradas[k], s.paradas[k + 1]
                cruzan = [x[o, d, s.nombre] for (o, d) in pares
                          if o in s.pos and d in s.pos and (d - o) * sentido > 0
                          and min(o, d) <= min(i, j) and max(o, d) >= max(i, j)]
                if cruzan:
                    m += pulp.lpSum(cruzan) <= Q * s.f * y[s.nombre]
    # 4. flota
    m += pulp.lpSum(s.buses * y[s.nombre] for s in servicios) <= flota, "flota"
    # 5. el corriente opera con exactamente una frecuencia; cada patron a lo sumo una
    patrones = {}
    for s in servicios:
        patrones.setdefault(s.patron, []).append(y[s.nombre])
    for p, ys in patrones.items():
        if p == "Corriente":
            m += pulp.lpSum(ys) == 1, "corriente"
        else:
            m += pulp.lpSum(ys) <= 1
    for nombre, valor in (fijar or {}).items():
        m += y[nombre] == valor, f"fijar_{idx[nombre]}"
    if limite_operacion is not None:          # epsilon-restriccion (Pareto, Bloque C)
        m += costo_op <= limite_operacion, "epsilon_operacion"

    stats = m.solve(pulp.HiGHS(msg=False, timeLimit=limite_seg))   # PuLP 4 trae HiGHS
    estado = stats.status.name          # Optimal, Infeasible, TimeLimit, ...
    seg = time.perf_counter() - t0
    if not stats.has_solution or estado not in ("Optimal", "TimeLimit"):
        return ResultadoMILP(estado, math.inf, math.inf, math.inf, [], 0, flota, seg, math.inf)

    elegidos = [s for s in servicios if (y[s.nombre].value() or 0) > 1e-6]
    op, us = pulp.value(costo_op), pulp.value(costo_us)
    sombra, duales = None, {}
    if relajar:
        for cons in m.constraints():
            if cons.name == "flota":
                sombra = cons.pi
            elif cons.name.startswith("dem_"):
                _, o, d = cons.name.split("_")
                duales[(int(o), int(d))] = cons.pi
    carga = {s.nombre: float(y[s.nombre].value()) for s in elegidos}
    # flujos de cada servicio elegido: pasajeros totales y carga del tramo critico
    flujos = {}
    for s in elegidos:
        xs = {(o, d): x[o, d, s.nombre].value() or 0.0 for (o, d) in pares if (o, d, s.nombre) in x}
        tramo_max = 0.0
        for sentido in (1, -1):
            for k in range(len(s.paradas) - 1):
                i, j = s.paradas[k], s.paradas[k + 1]
                tramo_max = max(tramo_max, sum(val for (o, d), val in xs.items()
                                               if (d - o) * sentido > 0 and min(o, d) <= min(i, j)
                                               and max(o, d) >= max(i, j)))
        flujos[s.nombre] = {"pasajeros": sum(xs.values()), "carga_max": tramo_max,
                            "horas_vehiculo": sum(val * tiempo_en_vehiculo_h(red, s, o, d, lluvia)
                                                  for (o, d), val in xs.items())}
    return ResultadoMILP(estado, op + us, op, us, elegidos,
                         int(round(sum(s.buses * (y[s.nombre].value() or 0) for s in servicios))),
                         flota, seg, us / v, sombra, carga, duales, flujos)


# ------------------------------------------------- ruta mas corta (A2)
def grafo_servicios(red: Red, servicios: list[Servicio], lluvia: float = 0.0):
    """Grafo dirigido: nodo estacion 'i' y nodos de anden (servicio, i).

    Abordar: i -> (s,i) con costo espera 1/(2f). Viajar: (s,i) -> (s,j) siguiente
    parada en cada sentido. Bajar: (s,i) -> i con costo 0. Costos en horas.
    """
    G: dict = {}

    def arco(u, w, c):
        G.setdefault(u, []).append((w, c))
        G.setdefault(w, [])

    for s in servicios:
        for e in s.paradas:
            arco(("est", e), ("srv", s.nombre, e), 1.0 / (2.0 * s.f))
            arco(("srv", s.nombre, e), ("est", e), 0.0)
        for k in range(len(s.paradas) - 1):
            i, j = s.paradas[k], s.paradas[k + 1]
            c = tiempo_tramo_h(red, i, j, lluvia) + PARAMS["tau_min"] / 60.0
            arco(("srv", s.nombre, i), ("srv", s.nombre, j), c)
            arco(("srv", s.nombre, j), ("srv", s.nombre, i), c)
    return G


def _estacion(nodo) -> int:
    return nodo[1] if nodo[0] == "est" else nodo[2]


def a_estrella(G, red: Red, origen: int, destino: int, heuristica: bool = True):
    """A*; con heuristica=False es Dijkstra. Devuelve (costo_h, ruta, expandidos)."""
    v_max = PARAMS["v_c"]           # ningun arco es mas rapido que v_c (admisible)

    def h(nodo):
        return abs(red.km[_estacion(nodo)] - red.km[destino]) / v_max if heuristica else 0.0

    ini, fin = ("est", origen), ("est", destino)
    g = {ini: 0.0}
    padre = {ini: None}
    cola = [(h(ini), 0, ini)]
    cerrados, contador, expandidos = set(), 1, 0
    while cola:
        _, _, u = heapq.heappop(cola)
        if u in cerrados:
            continue
        cerrados.add(u)
        expandidos += 1
        if u == fin:
            break
        for w, c in G[u]:
            nuevo = g[u] + c
            if nuevo < g.get(w, math.inf):
                g[w], padre[w] = nuevo, u
                heapq.heappush(cola, (nuevo + h(w), contador, w))
                contador += 1
    ruta, u = [], fin
    while u is not None:
        ruta.append(u)
        u = padre.get(u)
    return g.get(fin, math.inf), ruta[::-1], expandidos


def dijkstra(G, red: Red, origen: int, destino: int):
    return a_estrella(G, red, origen, destino, heuristica=False)


def describir_ruta(red: Red, ruta) -> str:
    tramos, actual = [], None
    for nodo in ruta:
        if nodo[0] == "srv":
            srv, e = nodo[1], nodo[2]
            if actual is None or actual[0] != srv:
                actual = [srv, e, e]
                tramos.append(actual)
            else:
                actual[2] = e
        else:
            actual = None
    return " | ".join(f"{srv.split('@')[0]}: {red.nombres[a]} -> {red.nombres[b]}"
                      for srv, a, b in tramos)


# ------------------------------------------------- Bloque B (no lineal)
# Frecuencia optima con espera y hacinamiento:
#   min Z(f) = sum_s  c T_s f_s  +  v D_s/(2 f_s)  +  v H_s (1 + kappa (L_s/(Q f_s))^2)
#   s.a.  sum_s T_s f_s <= N        (lambda >= 0)   flota
#         f_s >= L_s / Q            (mu_s  >= 0)    capacidad del tramo critico
# D_s pax/h del servicio, H_s horas-pasajero en vehiculo, L_s carga del tramo critico.
# Z es convexa en f > 0, asi que KKT es necesario y suficiente.
PARAMS.setdefault("kappa", 1.0)   # PROVISIONAL: a bus lleno el tiempo en vehiculo vale el doble


@dataclass
class ResultadoFrecuencias:
    metodo: str
    f: np.ndarray                # buses/h por servicio
    lam: float                   # multiplicador de la flota (COP/h por bus)
    mu: np.ndarray               # multiplicadores de capacidad minima
    costo_operacion: float
    costo_espera: float
    costo_viaje: float           # tiempo en vehiculo con hacinamiento
    buses: float
    factible: bool = True
    segundos: float = 0.0

    @property
    def objetivo(self) -> float:
        return self.costo_operacion + self.costo_espera + self.costo_viaje


def datos_frecuencias(res: ResultadoMILP, red: Red | None = None, lluvia: float = 0.0):
    """Insumos de B para los servicios elegidos en A: T, D, H, L y nombres.

    Con lluvia se recalcula el ciclo T_s y las horas en vehiculo H_s escalan con el
    tiempo de viaje en un sentido (sin la regulacion en terminal). D_s y L_s se mantienen.
    """
    T, D, H, L, nombres = [], [], [], [], []
    lay = 2 * PARAMS["layover_min"] / 60.0
    for s in res.elegidos:
        fl = res.flujos[s.nombre]
        h = fl["horas_vehiculo"]
        if red is not None and lluvia:
            s2 = _completar(red, Servicio(s.nombre, s.paradas, s.f, s.patron, s.corriente), lluvia)
            h *= (s2.T_ciclo_h - lay) / (s.T_ciclo_h - lay)
            s = s2
        T.append(s.T_ciclo_h)
        D.append(fl["pasajeros"])
        H.append(h)
        L.append(fl["carga_max"])
        nombres.append(s.patron)
    return np.array(T), np.array(D), np.array(H), np.array(L), nombres


def _costos(f, T, D, H, L, kappa=None):
    c, v, Q = PARAMS["c_bus_h"], PARAMS["v_tiempo_h"], PARAMS["Q"]
    kappa = PARAMS["kappa"] if kappa is None else kappa
    return (float(c * (f * T).sum()), float(v * (D / (2 * f)).sum()),
            float(v * (H * (1 + kappa * (L / (Q * f)) ** 2)).sum()))


def _derivada(f, T, D, H, L, kappa, lam):
    """dZ/df_s + lambda T_s  (estacionariedad sin el termino mu_s)."""
    c, v, Q = PARAMS["c_bus_h"], PARAMS["v_tiempo_h"], PARAMS["Q"]
    return (c + lam) * T - v * D / (2 * f ** 2) - 2 * v * kappa * H * L ** 2 / (Q ** 2 * f ** 3)


def _infactible(metodo, T, t0):
    nan = np.full_like(T, np.nan)
    return ResultadoFrecuencias(metodo, nan, math.nan, nan, math.inf, math.inf, math.inf,
                                math.inf, False, time.perf_counter() - t0)


def frecuencias_kkt(T, D, H, L, flota, kappa=None) -> ResultadoFrecuencias:
    """Resuelve el sistema KKT.

    Estacionariedad: (c+lambda) T_s - v D_s/(2 f^2) - 2 v kappa H_s L_s^2/(Q^2 f^3) - mu_s = 0.
    Para un lambda dado, cada f_s es la raiz (unica: la derivada crece con f) o fmin_s
    si la raiz cae por debajo (entonces mu_s > 0). lambda se ajusta por biseccion
    para que la flota se cumpla con igualdad cuando esta activa (complementariedad).
    Con kappa = 0 la raiz es cerrada: f_s = sqrt(v D_s / (2 T_s (c + lambda))) (regla de Mohring).
    """
    from scipy.optimize import brentq

    t0 = time.perf_counter()
    kappa = PARAMS["kappa"] if kappa is None else kappa
    fmin = L / PARAMS["Q"]
    if (fmin * T).sum() > flota + 1e-9:
        return _infactible("KKT", T, t0)

    def f_de(lam):
        f = np.empty_like(T)
        for i in range(len(T)):
            g = lambda x: _derivada(x, T[i], D[i], H[i], L[i], kappa, lam)
            lo = max(fmin[i], 1e-6)
            if g(lo) >= 0:          # la raiz esta por debajo de fmin: capacidad activa
                f[i] = lo
                continue
            hi = lo * 2
            while g(hi) < 0:
                hi *= 2
            f[i] = brentq(g, lo, hi, xtol=1e-12)
        return f

    lam = 0.0
    if (f_de(0.0) * T).sum() > flota:
        lo, hi = 0.0, PARAMS["c_bus_h"]
        while (f_de(hi) * T).sum() > flota:
            hi *= 2
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if (f_de(mid) * T).sum() > flota else (lo, mid)
        lam = hi
    f = f_de(lam)
    activa = f <= fmin * (1 + 1e-9)
    mu = np.where(activa, _derivada(f, T, D, H, L, kappa, lam), 0.0)
    op, es, vi = _costos(f, T, D, H, L, kappa)
    return ResultadoFrecuencias("KKT", f, lam, np.maximum(mu, 0.0), op, es, vi,
                                float((f * T).sum()), True, time.perf_counter() - t0)


def frecuencias_scipy(T, D, H, L, flota, kappa=None) -> ResultadoFrecuencias:
    """Mismo problema con scipy.optimize.minimize (SLSQP).

    lambda se estima por diferencias finitas: -(Z(N+d) - Z(N-d)) / (2d).
    """
    from scipy.optimize import minimize

    t0 = time.perf_counter()
    kappa = PARAMS["kappa"] if kappa is None else kappa
    fmin = L / PARAMS["Q"]
    escala = 1e6
    if (fmin * T).sum() > flota + 1e-9:
        return _infactible("scipy SLSQP", T, t0)

    def resolver(N):
        obj = lambda f: sum(_costos(f, T, D, H, L, kappa)) / escala
        grad = lambda f: _derivada(f, T, D, H, L, kappa, 0.0) / escala
        cons = [{"type": "ineq", "fun": lambda f: N - (f * T).sum(), "jac": lambda f: -T}]
        return minimize(obj, fmin * 1.01, jac=grad, method="SLSQP", constraints=cons,
                        bounds=[(m, None) for m in fmin], options={"ftol": 1e-14, "maxiter": 2000})

    r = resolver(flota)
    d = 0.5
    lam = 0.0
    if (fmin * T).sum() <= flota - d:
        lam = max(-(resolver(flota + d).fun - resolver(flota - d).fun) * escala / (2 * d), 0.0)
    activa = r.x <= fmin * (1 + 1e-6)
    mu = np.where(activa, np.maximum(_derivada(r.x, T, D, H, L, kappa, lam), 0.0), 0.0)
    op, es, vi = _costos(r.x, T, D, H, L, kappa)
    return ResultadoFrecuencias("scipy SLSQP", r.x, lam, mu, op, es, vi, float((r.x * T).sum()),
                                bool(r.success), time.perf_counter() - t0)


# ------------------------------------------------- Bloque C (metaheuristica)
# Espacio de busqueda: para cada patron (corriente + 24 expresos) un nivel de
# frecuencia en {apagado, 6, 10, 15, 20, 30, 40}. El corriente no se puede apagar.
# Es la MISMA instancia que el MILP del Bloque A (mismos candidatos).
class EvaluadorRapido:
    """Evalua un diseno en milisegundos con una asignacion voraz con capacidad.

    Cada par OD se asigna a sus servicios activos del mas barato al mas caro, hasta
    agotar la capacidad del tramo en ese sentido. Los pares con menos servicios
    activos van primero (los que solo puede llevar el corriente no quedan sin cupo).
    La demanda que no cabe y el exceso de flota se penalizan. Es una cota superior
    del costo exacto: la solucion final se reevalua con el MILP (y fijas).
    """

    def __init__(self, red: Red, od: np.ndarray, servicios: list[Servicio], flota: int,
                 lluvia: float = 0.0, penal_pax: float = 1e6, penal_bus: float = 5e6):
        self.red, self.od, self.flota = red, od, flota
        self.servicios = servicios
        self.penal_pax, self.penal_bus = penal_pax, penal_bus
        v = PARAMS["v_tiempo_h"]
        self.patrones = list(dict.fromkeys(s.patron for s in servicios))
        self.niveles = {p: sorted([k for k, s in enumerate(servicios) if s.patron == p],
                                  key=lambda k: servicios[k].f) for p in self.patrones}
        self.pares = sorted(((o, d) for o in range(red.n) for d in range(red.n)
                             if o != d and od[o, d] > 1e-6), key=lambda p: -od[p])
        # opciones de cada par: (indice servicio, costo por pasajero)
        self.opciones = {p: [] for p in self.pares}
        for k, s in enumerate(servicios):
            for (o, d) in self.pares:
                if o in s.pos and d in s.pos:
                    costo = v * (tiempo_en_vehiculo_h(red, s, o, d, lluvia) + 1.0 / (2.0 * s.f))
                    self.opciones[(o, d)].append((k, costo))
        for p in self.opciones:
            self.opciones[p].sort(key=lambda t: t[1])
        self.evaluaciones = 0

    def activos(self, estado: dict[str, int]) -> set[int]:
        return {self.niveles[p][n] for p, n in estado.items() if n >= 0}

    def evaluar(self, estado: dict[str, int]):
        """Devuelve (costo_penalizado, costo_operacion, costo_usuarios, pax_sin_atender, buses)."""
        self.evaluaciones += 1
        act = self.activos(estado)
        Q = PARAMS["Q"]
        cap = {}
        us, sin = 0.0, 0.0
        # primero los pares con menos servicios activos (los mas restringidos)
        orden = sorted(self.pares, key=lambda p: (sum(k in act for k, _ in self.opciones[p]), -self.od[p]))
        for (o, d) in orden:
            resto = self.od[o, d]
            lo, hi = min(o, d), max(o, d)
            sentido = 1 if d > o else -1
            for k, costo in self.opciones[(o, d)]:
                if k not in act:
                    continue
                arr = cap.get((k, sentido))
                if arr is None:
                    arr = cap[(k, sentido)] = np.full(self.red.n - 1, Q * self.servicios[k].f)
                libre = arr[lo:hi].min()
                a = min(resto, libre)
                if a > 1e-9:
                    arr[lo:hi] -= a
                    us += a * costo
                    resto -= a
                    if resto <= 1e-9:
                        break
            sin += max(resto, 0.0)
        buses = sum(self.servicios[k].buses for k in act)
        op = PARAMS["c_bus_h"] * buses
        pen = self.penal_pax * sin + self.penal_bus * max(0, buses - self.flota)
        return op + us + pen, op, us, sin, buses

    def a_fijar(self, estado: dict[str, int]) -> dict[str, int]:
        """Convierte un diseno en variables fijas para reevaluarlo con el MILP."""
        act = self.activos(estado)
        return {s.nombre: int(k in act) for k, s in enumerate(self.servicios)}


@dataclass
class ResultadoSA:
    estado: dict
    costo: float
    historia: list
    segundos: float
    evaluaciones: int


def recocido_simulado(ev: EvaluadorRapido, iteraciones: int = 3000, semilla: int = 0,
                      t0_rel: float = 0.01, tf_rel: float = 1e-5) -> ResultadoSA:
    """Simulated Annealing con enfriamiento geometrico.

    Vecindario: cambiar el nivel de frecuencia de un patron al azar (o apagarlo);
    con probabilidad 0,3 ademas se apaga otro expreso activo (intercambio).
    Aceptacion de Metropolis: exp(-delta / T).
    """
    import random

    rng = random.Random(semilla)
    t_ini = time.perf_counter()
    ev.evaluaciones = 0
    estado = {p: -1 for p in ev.patrones}
    estado["Corriente"] = len(ev.niveles["Corriente"]) - 1
    costo = ev.evaluar(estado)[0]
    mejor, mejor_costo = dict(estado), costo
    T = t0_rel * costo
    alfa = (tf_rel / t0_rel) ** (1.0 / iteraciones)
    historia = [mejor_costo]
    expresos = [p for p in ev.patrones if p != "Corriente"]
    for _ in range(iteraciones):
        vec = dict(estado)
        p = rng.choice(ev.patrones)
        tope = len(ev.niveles[p]) - 1
        opciones = [n for n in range(0 if p == "Corriente" else -1, tope + 1) if n != vec[p]]
        vec[p] = rng.choice(opciones)
        if rng.random() < 0.3:
            encendidos = [q for q in expresos if vec[q] >= 0 and q != p]
            if encendidos:
                vec[rng.choice(encendidos)] = -1
        c_vec = ev.evaluar(vec)[0]
        delta = c_vec - costo
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            estado, costo = vec, c_vec
            if costo < mejor_costo:
                mejor, mejor_costo = dict(estado), costo
        T *= alfa
        historia.append(mejor_costo)
    return ResultadoSA(mejor, mejor_costo, historia, time.perf_counter() - t_ini, ev.evaluaciones)


def frontera_pareto(red: Red, od: np.ndarray, flota: int, servicios: list[Servicio],
                    puntos: int = 8, lluvia: float = 0.0):
    """Frontera costo de operacion vs horas-usuario con epsilon-restriccion.

    Extremos: (1) minimo costo de operacion factible (peso usuarios ~0),
    (2) minimo tiempo de usuarios (peso operacion ~0). Entre ambos se minimiza el
    tiempo de usuarios con costo de operacion <= epsilon.
    """
    r_op = resolver_milp(red, od, flota, lluvia, servicios=servicios,
                         peso_operacion=1.0, peso_usuarios=1e-6)
    r_us = resolver_milp(red, od, flota, lluvia, servicios=servicios,
                         peso_operacion=1e-6, peso_usuarios=1.0)
    eps = np.linspace(r_op.costo_operacion, r_us.costo_operacion, puntos)
    salida = []
    for e in eps:
        r = resolver_milp(red, od, flota, lluvia, servicios=servicios, peso_operacion=1e-6,
                          peso_usuarios=1.0, limite_operacion=float(e) + 1.0)
        if r.estado in ("Optimal", "TimeLimit"):
            salida.append({"epsilon": float(e), "costo_operacion": r.costo_operacion,
                           "horas_usuario": r.horas_usuario, "buses": r.buses,
                           "servicios": len(r.elegidos), "resultado": r})
    return salida


# ------------------------------------------------- Sensibilidad: estacion cerrada
def redistribuir_od(od: np.ndarray, cerradas) -> np.ndarray:
    """Los viajes de una estacion cerrada pasan a sus vecinas abiertas (mitad y mitad).

    Supuesto E3: quien usaba la estacion cerrada camina a la anterior o a la siguiente.
    Los viajes que quedan con origen = destino (vecinos entre si) se descartan.
    """
    od = od.copy()
    n = od.shape[0]
    for c in sorted(cerradas):
        vecinas = [i for i in (c - 1, c + 1) if 0 <= i < n and i not in cerradas]
        for v in vecinas:
            od[v, :] += od[c, :] / len(vecinas)
            od[:, v] += od[:, c] / len(vecinas)
        od[c, :] = 0.0
        od[:, c] = 0.0
    np.fill_diagonal(od, 0.0)
    return od
