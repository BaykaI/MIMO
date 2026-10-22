"""Анализ экспорта ДН из CST (этапы 2–4): фазовые центры элементов, подрешёток и произведений.

  python cst/analyze_stage.py cst/stage2            # все частоты, сектор 30°
  python cst/analyze_stage.py cst/stage2 --sector 20

Выход: таблица в консоль + <папка>/phase_centers.md (Markdown) + phase_centers.json.
Ожидание (геометрия): Rx k -> (x_k, 0, 0); Tx n -> (x_n, 30, 0); сумма подрешётки -> среднее по
элементам (синфазное сложение одинаковых элементов даёт фазовый центр в центре подрешётки);
произведение F_tx·F_rx -> сумма фазовых центров (гипотеза Г1).
"""
import sys, os, json, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from phase_center_fit import read_cst_farfield, fit_phase_center, C0

RX_X = [-105, -75, -45, -15, 15, 45, 75, 105]
TX_X = [-127.5, -112.5, 112.5, 127.5]
Y_TX = 30.0
GEOM = {f"rx{i+1}": np.array([x, 0.0, 0.0]) for i, x in enumerate(RX_X)}
GEOM.update({f"tx{i+1}": np.array([x, Y_TX, 0.0]) for i, x in enumerate(TX_X)})

ap = argparse.ArgumentParser()
ap.add_argument("folder")
ap.add_argument("--sector", type=float, default=30.0, help="сектор для элементов и пар, °")
ap.add_argument("--sector-sub", type=float, default=3.0,
                help="сектор для подрешёток и их произведений, ° (внутри главного лепестка: "
                     "8 элементов с шагом λ — первый нуль 7,2°)")
ap.add_argument("--freqs", type=float, nargs="*", default=[9.8, 10.0, 10.2])
ap.add_argument("--component", default="copol",
                help="theta | phi | copol (Ludwig-3 для диполя вдоль y: E_th·sinφ + E_ph·cosφ — без скачков знака)")
a = ap.parse_args()

rows = []          # (freq, name, kind, expected[3], fitted[3], rms, dev)
patterns = {}
phases = {}        # (freq, name) -> φ0 элемента, °


def load(name, f):
    key = (name, f)
    if key not in patterns:
        path = os.path.join(a.folder, f"ff_{name}_{f:.1f}GHz.txt")
        if a.component == "copol":
            th, ph, Et = read_cst_farfield(path, "theta")
            _, _, Ep = read_cst_farfield(path, "phi")
            phr = np.radians(ph)
            patterns[key] = (th, ph, Et * np.sin(phr) + Ep * np.cos(phr))
        else:
            patterns[key] = read_cst_farfield(path, a.component)
    return patterns[key]


def fit(E, th, ph, f, sector, fix_z):
    """Если остаток > 30° (ложный минимум на изрезанной встроенной ДН крайнего элемента),
    повторить с секторами 20° и 15° и взять лучший."""
    k = 2 * np.pi * f * 1e9 / C0
    best = None
    for s in [sector] + ([20.0, 15.0] if not fix_z else []):
        xyz, rms, ph0 = fit_phase_center(th, ph, E, k, s, (0.0, 0.0), fix_z=fix_z, with_phase=True)
        if best is None or rms < best[1]:
            best = (xyz, rms, ph0, s)
        if rms < 30:
            break
    xyz, rms, ph0, s = best
    if s != sector:
        print(f"  [{f} ГГц] сектор уменьшен до {s:g}° (rms {rms:.1f}°)")
    return xyz * 1e3, rms, ph0


def add(f, name, kind, exp, E, th, ph, sector=None):
    """Элементы и пары: полный фит (x, y, z) в широком секторе; подрешётки и их произведения —
    узкий сектор внутри главного лепестка, z фиксирован (не обусловлен)."""
    xyz, rms, ph0 = fit(E, th, ph, f, sector or a.sector, fix_z=sector is not None)
    dev_xy = float(np.abs(xyz[:2] - exp[:2]).max())
    rows.append((f, name, kind, exp, xyz, rms, dev_xy))
    phases[(f, name)] = ph0


for f in a.freqs:
    lam = C0 / (f * 1e9) * 1e3
    # элементы
    for name, g in GEOM.items():
        th, ph, E = load(name, f)
        add(f, name, "элемент", g, E, th, ph)
    # подрешётки (синфазная сумма) — ожидание: центр подрешётки (среднее координат)
    subs = {"rxsum": [f"rx{i}" for i in range(1, 9)], "txsum": [f"tx{i}" for i in range(1, 5)],
            "txL": ["tx1", "tx2"], "txR": ["tx3", "tx4"], "rxL4": ["rx1", "rx2", "rx3", "rx4"]}
    sumpat = {}
    for sname, members in subs.items():
        th, ph, E = load(members[0], f)
        E = E.copy()
        for mname in members[1:]:
            E += load(mname, f)[2]
        sumpat[sname] = (th, ph, E)
        exp = np.mean([GEOM[mname] for mname in members], axis=0)
        add(f, sname, "подрешётка", exp, E, th, ph, a.sector_sub)
    # произведения: пара, подрешётка×подрешётка — ожидание: сумма фазовых центров (Г1)
    fitted = {r[1]: r[4] for r in rows if r[0] == f}
    prods = [("tx1*rx8", "tx1", "rx8"), ("tx4*rx1", "tx4", "rx1"), ("tx2*rx4", "tx2", "rx4"),
             ("txsum*rxsum", "txsum", "rxsum"), ("txL*rxsum", "txL", "rxsum")]
    for pname, p1, p2 in prods:
        th1, ph1, E1 = sumpat[p1] if p1 in sumpat else load(p1, f)
        th2, ph2, E2 = sumpat[p2] if p2 in sumpat else load(p2, f)
        add(f, pname, "произведение", fitted[p1] + fitted[p2], E1 * E2, th1, ph1,
            a.sector_sub if (p1 in sumpat or p2 in sumpat) else None)

# ---------------------------------------------------------------- вывод
lines = ["| f, ГГц | объект | тип | ожидание x,y,z, мм | фит x,y,z, мм | откл. по x,y, мм | rms, ° |",
         "|---|---|---|---|---|---|---|"]
for f, name, kind, exp, xyz, rms, dev in rows:
    e = ", ".join(f"{v:.1f}" for v in exp)
    x = ", ".join(f"{v:.2f}" for v in xyz)
    lines.append(f"| {f:.1f} | {name} | {kind} | {e} | {x} | {dev:.2f} | {rms:.1f} |")
worst = {}
for f, name, kind, exp, xyz, rms, dev in rows:
    worst[kind] = max(worst.get(kind, 0), dev)
lam0 = 30.0
summ = ["", "Максимальное отклонение по x, y по типам (λ/20 = 1,5 мм, λ/50 = 0,6 мм):"]
for kind, d in worst.items():
    summ.append(f"- {kind}: {d:.2f} мм = λ/{lam0/d:.0f}" if d > 0 else f"- {kind}: 0")
ph_lines = ["", "Фаза элементов φ0 относительно rx4 (°; разнофазность — калибровочная поправка ВАР):",
            "| f, ГГц | " + " | ".join(GEOM) + " |", "|---|" + "---|" * len(GEOM)]
for f in a.freqs:
    ref = phases[(f, "rx4")]
    ph_lines.append(f"| {f:.1f} | " + " | ".join(
        f"{np.degrees(np.angle(np.exp(1j*np.radians(phases[(f, n)] - ref)))):.1f}" for n in GEOM) + " |")
txt = "\n".join(lines + summ + ph_lines)
print(txt.replace(".", ","))
with open(os.path.join(a.folder, "phase_centers.md"), "w", encoding="utf-8") as fh:
    fh.write(f"# Фазовые центры ({os.path.basename(os.path.abspath(a.folder))}: сектор ±{a.sector:g}° "
             f"для элементов и пар, ±{a.sector_sub:g}° для подрешёток; компонента {a.component})\n\n")
    fh.write(txt.replace(".", ",") + "\n")
with open(os.path.join(a.folder, "phase_centers.json"), "w", encoding="utf-8") as fh:
    json.dump([dict(f=f, name=n, kind=k, expected=list(map(float, e)), fit=list(map(float, x)),
                    rms=float(r), dev=d) for f, n, k, e, x, r, d in rows], fh, ensure_ascii=False, indent=1)
