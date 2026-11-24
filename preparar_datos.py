"""Prepara la red NQS (lineas E y G) y la demanda para el parcial.

Entradas (carpeta datos/):
  geo/estaciones_troncales.geojson, salidas/salidas_troncal_nov2024.csv
Salidas:
  datos/nqs_red.csv          estaciones en orden de norte a sur, coordenadas,
                             distancia acumulada y demanda (salidas) en dia habil
  datos/nqs_demanda_15min.csv  salidas promedio por intervalo de 15 min y estacion
Supuestos que implementa este script (ver supuestos.md):
  S1 Corredor lineal Castellana -> Ricaurte -> Comuneros -> San Mateo (28 nodos).
     Tygua y Guatoque (ramal Carrera 10) quedan fuera del eje.
  S2 SENA se excluye: figura cerrada en el Informe #99 y no aparece en las salidas.
  S3 "Ampliacion San Mateo" (57503) se suma a San Mateo (07503).
  S4 Dia habil de noviembre 2024: se excluyen sabados, domingos y los festivos
     del 4 y 11 de noviembre.
"""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent / "datos"

# Orden norte -> sur sobre el eje NQS (codigo num_est). Ricaurte conecta E con G.
ORDEN = [
    "07101", "07102", "07103", "07104", "07105", "07106", "07107", "07108",
    "07109", "07110", "07111",                      # zona E (hasta Ricaurte)
    "07112", "07113", "07008", "07007", "07006", "07005", "07004", "07003",
    "07002", "07001", "07000", "07010", "07506", "07505", "07504", "07503",  # zona G
]
FUSIONES = {"57503": "07503"}  # S3
DIAS_NO_HABILES = {2, 3, 4, 9, 10, 11, 16, 17, 23, 24, 30}  # S4


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def cargar_estaciones():
    gj = json.loads((BASE / "geo" / "estaciones_troncales.geojson").read_text(encoding="utf8"))
    est = {}
    for f in gj["features"]:
        p = f["properties"]
        est[str(p["num_est"])] = p
    return est


def cargar_demanda():
    """Devuelve {codigo: {intervalo: promedio_salidas_dia_habil}}."""
    with open(BASE / "salidas" / "salidas_troncal_nov2024.csv", encoding="cp1252", newline="") as fh:
        filas = list(csv.reader(fh, delimiter=";"))
    cols_dia = list(range(5, 36))  # DIA 01 .. DIA 31
    habiles = [i for i, d in zip(cols_dia, range(1, 32)) if d not in DIAS_NO_HABILES]
    acum = defaultdict(lambda: defaultdict(float))
    for r in filas[1:]:
        if len(r) < 36 or not r[1]:
            continue
        codigo = r[1].split(")")[0].strip("( ")
        codigo = FUSIONES.get(codigo, codigo)
        if codigo not in ORDEN:
            continue
        total = 0.0
        for i in habiles:
            v = r[i].strip()
            if v.replace(".", "", 1).isdigit():
                total += float(v)
        acum[codigo][r[4]] += total / len(habiles)
    return acum


def hora_pico(serie):
    """Mayor suma de 4 intervalos consecutivos de 15 min."""
    horas = sorted(serie)
    vals = [serie[h] for h in horas]
    mejor, inicio = 0.0, horas[0]
    for k in range(len(vals) - 3):
        s = sum(vals[k:k + 4])
        if s > mejor:
            mejor, inicio = s, horas[k]
    return mejor, inicio


def main():
    est = cargar_estaciones()
    dem = cargar_demanda()
    faltan = [c for c in ORDEN if c not in est or c not in dem]
    if faltan:
        raise SystemExit(f"Estaciones sin datos: {faltan}")

    filas, acum_km, previo = [], 0.0, None
    for k, c in enumerate(ORDEN):
        p = est[c]
        lat, lon = p["latitud"], p["longitud"]
        if previo is not None:
            acum_km += haversine_km(previo[0], previo[1], lat, lon)
        previo = (lat, lon)
        total = sum(dem[c].values())
        pico, inicio = hora_pico(dem[c])
        filas.append({
            "orden": k, "codigo": c, "nombre": p["nom_est"], "zona": "E" if k <= 10 else "G",
            "lat": lat, "lon": lon, "km_acum": round(acum_km, 3),
            "cap_biart": p["cap_biart"], "cap_art": p["cap_art"], "num_vagones": p["num_vag"],
            "salidas_dia_habil": round(total), "salidas_hora_pico": round(pico),
            "inicio_hora_pico": inicio,
        })
    campos = list(filas[0])
    with open(BASE / "nqs_red.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)

    intervalos = sorted({h for c in ORDEN for h in dem[c]})
    with open(BASE / "nqs_demanda_15min.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["intervalo"] + ORDEN)
        for h in intervalos:
            w.writerow([h] + [round(dem[c].get(h, 0.0), 2) for c in ORDEN])

    print(f"{len(filas)} estaciones, longitud del eje {acum_km:.1f} km")
    for f in filas:
        print(f"{f['orden']:2d} {f['codigo']} {f['nombre'][:34]:34s} km={f['km_acum']:5.1f} "
              f"dia={f['salidas_dia_habil']:7d} pico={f['salidas_hora_pico']:6d} ({f['inicio_hora_pico']})")


if __name__ == "__main__":
    main()
