"""Проверка статьи [1]: ДН виртуальной решётки «по формуле» против «по CST».

Формула ВАР ([1], ф. (4)/(10)): F_var(u) = Σ_n Σ_k exp(j·k·(r_tx,n + r_rx,k)·u) — изотропные элементы в
геометрических точках. CST: F_cst(u) = Σ_n Σ_k F_n(u)·F_k(u) — произведение встроенных ДН каждой пары
(это и есть отклик пары «передатчик–приёмник» MIMO-РЛС в дальней зоне). Калибровка каналов: каждая
встроенная ДН делится на свою комплексную величину по нормали F_n(θ=0) (так калибруют реальную MIMO-РЛС).
Формула умножается на квадрат среднего множителя элемента (средняя |F_n|; пара — произведение двух ДН),
чтобы сравнивать только решётку.

  python cst/var_compare.py cst/stage3 --f 10 --cut 0        # линейная 8r4t, сечение φ=0 (вдоль оси x)
  python cst/var_compare.py cst/ring_R5 --f 10 --cut 0 90    # кольцо, два сечения

Выход: <папка>/var_<f>GHz.png, var_<f>GHz.md (ширина луча по −3 дБ, УБЛ, положение максимума).
"""
import sys, os, argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from phase_center_fit import read_cst_farfield, C0

ap = argparse.ArgumentParser()
ap.add_argument("folder")
ap.add_argument("--f", type=float, default=10.0)
ap.add_argument("--cut", type=float, nargs="+", default=[0.0])
ap.add_argument("--span", type=float, default=80.0, help="полуширина сечения, °")
ap.add_argument("--zoom", type=float, default=8.0, help="полуширина увеличенного сечения, °")
a = ap.parse_args()
LAM = C0 / (a.f * 1e9) * 1e3
k = 2 * np.pi / LAM

# раскладка
lay = os.path.join(a.folder, "layout_mm.txt")
if os.path.exists(lay):
    elems = [(p[0], float(p[1]), float(p[2])) for p in (l.split() for l in open(lay, encoding="utf-8")) if p and p[0] != "name"]
else:
    elems = [(f"rx{i+1}", x, 0.0) for i, x in enumerate([-105, -75, -45, -15, 15, 45, 75, 105])] + \
            [(f"tx{i+1}", x, 30.0) for i, x in enumerate([-127.5, -112.5, 112.5, 127.5])]
rx = [(n, x, y) for n, x, y in elems if n.startswith("rx")]
tx = [(n, x, y) for n, x, y in elems if n.startswith("tx")]


def load(name):
    path = os.path.join(a.folder, f"ff_{name}_{a.f:.1f}GHz.txt")
    th, ph, Et = read_cst_farfield(path, "theta")
    _, _, Ep = read_cst_farfield(path, "phi")
    phr = np.radians(ph)
    return th, ph, Et * np.sin(phr) + Ep * np.cos(phr)


pats = {}
for n, x, y in elems:
    th, ph, E = load(n)
    pats[n] = E
u = np.column_stack([np.sin(np.radians(th)) * np.cos(np.radians(ph)),
                     np.sin(np.radians(th)) * np.sin(np.radians(ph)), np.cos(np.radians(th))])
i0 = np.argmin(th + ph)  # θ=0, φ=0 — нормаль
# калибровка: по комплексному значению на нормали (берём среднее по φ при θ=0 — там все φ дают одну точку)
norm = {n: np.mean(pats[n][th == 0]) for n in pats}

F_var = np.zeros(len(th), complex)
F_cst = np.zeros(len(th), complex)
F_cal = np.zeros(len(th), complex)
for tn, tx_x, tx_y in tx:
    for rn, rx_x, rx_y in rx:
        F_var += np.exp(1j * k * ((tx_x + rx_x) * u[:, 0] + (tx_y + rx_y) * u[:, 1]))
        F_cst += pats[tn] * pats[rn]
        F_cal += (pats[tn] / norm[tn]) * (pats[rn] / norm[rn])
elem = np.mean([np.abs(pats[n]) for n in pats], axis=0)
F_var_el = F_var * (elem / elem.max()) ** 2   # пара = произведение двух ДН элемента


def cut(F, phi0):
    """Сечение: θ от −span до +span при φ = phi0 (отрицательные θ = φ+180)."""
    sel1 = (ph == phi0) & (th <= a.span)
    sel2 = (ph == (phi0 + 180) % 360) & (th <= a.span) & (th > 0)
    t = np.concatenate([-th[sel2], th[sel1]]); v = np.concatenate([F[sel2], F[sel1]])
    o = np.argsort(t)
    return t[o], v[o]


def metrics(t, v):
    d = 20 * np.log10(np.abs(v) / np.abs(v).max() + 1e-12)
    im = int(np.argmax(d))
    above = np.where(d >= -3)[0]
    # ширина по −3 дБ: непрерывный участок вокруг максимума с интерполяцией по краям
    lo, hi = im, im
    while lo > 0 and d[lo - 1] >= -3: lo -= 1
    while hi < len(d) - 1 and d[hi + 1] >= -3: hi += 1
    def interp(i, j):
        return t[i] + (t[j] - t[i]) * (-3 - d[i]) / (d[j] - d[i]) if j != i and d[j] != d[i] else t[i]
    tl = interp(lo, lo - 1) if lo > 0 else t[lo]
    tr = interp(hi, hi + 1) if hi < len(d) - 1 else t[hi]
    # УБЛ: максимум вне главного лепестка (за первыми нулями)
    l, r = im, im
    while l > 0 and d[l - 1] < d[l]: l -= 1
    while r < len(d) - 1 and d[r + 1] < d[r]: r += 1
    outside = np.concatenate([d[:l], d[r + 1:]])
    sll = outside.max() if outside.size else float("nan")
    return dict(peak_deg=float(t[im]), hpbw=float(tr - tl), sll_db=float(sll))


rows = []
fig, axes = plt.subplots(len(a.cut), 2, figsize=(11, 3.4 * len(a.cut)), squeeze=False)
for ci, phi0 in enumerate(a.cut):
    for label, F, style in (("формула ВАР [1] × множитель элемента", F_var_el, dict(color="0.25", lw=1.4)),
                            ("CST: Σ F_tx·F_rx без калибровки", F_cst, dict(color="tab:red", lw=1.0, ls="--")),
                            ("CST: с калибровкой каналов по нормали", F_cal, dict(color="tab:blue", lw=1.0))):
        t, v = cut(F, phi0)
        d = 20 * np.log10(np.abs(v) / np.abs(v).max() + 1e-12)
        mt = metrics(t, v); mt.update(cut=phi0, model=label); rows.append(mt)
        for ax, zoom in ((axes[ci][0], a.span), (axes[ci][1], a.zoom)):
            ax.plot(t, d, label=label, **style)
            ax.set_xlim(-zoom, zoom)
    for ax, ttl in ((axes[ci][0], f"сечение φ = {phi0:g}°, ±{a.span:g}°"), (axes[ci][1], f"главный лепесток, ±{a.zoom:g}°")):
        ax.set_ylim(-45 if ax is axes[ci][0] else -30, 1); ax.axhline(-3, color="r", lw=0.6, ls=":")
        ax.set_xlabel("θ, °"); ax.set_ylabel("дБ"); ax.set_title(ttl, fontsize=9); ax.grid(alpha=0.3)
    axes[ci][0].legend(fontsize=7, loc="lower left")
fig.suptitle(f"{os.path.basename(os.path.abspath(a.folder))}: ДН виртуальной решётки, {a.f:g} ГГц "
             f"({len(tx)} Tx × {len(rx)} Rx = {len(tx)*len(rx)} вирт. элементов)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(a.folder, f"var_{a.f:g}GHz.png"), dpi=170)

lines = ["| сечение | модель | макс., ° | ширина −3 дБ, ° | УБЛ, дБ |", "|---|---|---|---|---|"]
for r in rows:
    lines.append(f"| φ={r['cut']:g}° | {r['model']} | {r['peak_deg']:.1f} | {r['hpbw']:.2f} | {r['sll_db']:.1f} |")
txt = "\n".join(lines).replace(".", ",")
print(txt)
open(os.path.join(a.folder, f"var_{a.f:g}GHz.md"), "w", encoding="utf-8").write(txt + "\n")
json.dump(rows, open(os.path.join(a.folder, f"var_{a.f:g}GHz.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
