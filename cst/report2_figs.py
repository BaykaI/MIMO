"""Рисунки для отчёта-2 (проверка статьи [1] на линейной решётке и кольце).

  python cst/report2_figs.py cst/ring_R5   -> fig2/<папка>_layout.png, _map.png, _phase.png
Карта ДН: |F(θ,φ)| в координатах u = sinθ cosφ, v = sinθ sinφ — формула ВАР, CST, разность.
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from phase_center_fit import read_cst_farfield, fit_phase_center, C0

folder = sys.argv[1]
F0 = 10.0
LAM = C0 / (F0 * 1e9) * 1e3
k = 2 * np.pi / LAM          # 1/мм -> координаты фита в мм
tag = os.path.basename(os.path.abspath(folder))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2")
os.makedirs(OUT, exist_ok=True)
NL = chr(10)

lay = os.path.join(folder, "layout_mm.txt")
if os.path.exists(lay):
    elems = [(p[0], float(p[1]), float(p[2])) for p in (l.split() for l in open(lay, encoding="utf-8")) if p and p[0] != "name"]
else:
    elems = [(f"rx{i+1}", x, 0.0) for i, x in enumerate([-105, -75, -45, -15, 15, 45, 75, 105])] + \
            [(f"tx{i+1}", x, 30.0) for i, x in enumerate([-127.5, -112.5, 112.5, 127.5])]
rx = [e for e in elems if e[0].startswith("rx")]; tx = [e for e in elems if e[0].startswith("tx")]


def load(name):
    path = os.path.join(folder, f"ff_{name}_{F0:.1f}GHz.txt")
    th, ph, Et = read_cst_farfield(path, "theta"); _, _, Ep = read_cst_farfield(path, "phi")
    phr = np.radians(ph)
    return th, ph, Et * np.sin(phr) + Ep * np.cos(phr)


pats = {}
for n, x, y in elems:
    th, ph, pats[n] = load(n)
u = np.sin(np.radians(th)) * np.cos(np.radians(ph)); v = np.sin(np.radians(th)) * np.sin(np.radians(ph))
F_var = np.zeros(len(th), complex); F_cst = np.zeros(len(th), complex)
for tn, tx_x, tx_y in tx:
    for rn, rx_x, rx_y in rx:
        F_var += np.exp(1j * k * ((tx_x + rx_x) * u + (tx_y + rx_y) * v))
        F_cst += pats[tn] * pats[rn]
elem = np.mean([np.abs(pats[n]) for n in pats], axis=0)
F_var *= (elem / elem.max()) ** 2

# --- 1. фазовые центры и φ0 элементов
fits = {}
for n, x, y in elems:
    best = None
    for s in (30, 20, 15):
        xyz, rms, p0 = fit_phase_center(th, ph, pats[n], k, s, (0, 0), with_phase=True)
        if best is None or rms < best[1]:
            best = (xyz, rms, p0)
        if rms < 30:
            break
    fits[n] = best
ref = fits[elems[0][0]][2]
wrap = lambda a: np.degrees(np.angle(np.exp(1j * np.radians(a - ref))))
fig, ax = plt.subplots(figsize=(7.2, 5.4 if os.path.exists(lay) else 3.4))
for lst, mk, lab in ((rx, "o", "Rx"), (tx, "s", "Tx")):
    ax.plot([e[1] for e in lst], [e[2] for e in lst], mk, mfc="none", mec="0.3", ms=11, mew=1.2, label=f"{lab}: геометрия")
    sc = ax.scatter([fits[e[0]][0][0] for e in lst], [fits[e[0]][0][1] for e in lst], c=[wrap(fits[e[0]][2]) for e in lst],
                    cmap="coolwarm", vmin=-180, vmax=180, s=28, marker=mk, zorder=3)
cb = fig.colorbar(sc, ax=ax, shrink=0.8); cb.set_label(f"φ0 элемента относительно {elems[0][0]}, °")
dev = max(max(abs(fits[n][0][0] - x), abs(fits[n][0][1] - y)) for n, x, y in elems)
ax.set_title(f"{tag}: фазовые центры (закрашено, цвет — фаза φ0) и геометрия (контур)" + NL +
             f"макс. расхождение центра и геометрии {dev:.2f} мм = λ/{LAM/dev:.0f}", fontsize=9)
ax.set_xlabel("x, мм"); ax.set_ylabel("y, мм"); ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{tag}_layout.png"), dpi=190); plt.close(fig)

# --- 2. карта ДН ВАР на полусфере (u, v)
sel = th <= 90
U, V = u[sel], v[sel]
N = 361
gu = np.linspace(-1, 1, N); GU, GV = np.meshgrid(gu, gu)
maps = []
for F in (F_var, F_cst):
    d = 20 * np.log10(np.abs(F[sel]) / np.abs(F).max() + 1e-12)
    g = griddata((U, V), d, (GU, GV), method="linear")
    g[GU ** 2 + GV ** 2 > 1] = np.nan
    maps.append(g)
fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3))
titles = ("формула ВАР [1]:" + NL + "геометрия × множитель элемента", "CST: Σ F_tx·F_rx" + NL + "полноволновые встроенные ДН")
for ax, g, ttl in zip(axes[:2], maps, titles):
    im = ax.imshow(g, extent=[-1, 1, -1, 1], origin="lower", vmin=-40, vmax=0, cmap="viridis")
    ax.set_title(ttl, fontsize=9); ax.set_xlabel("u = sinθ·cosφ"); ax.set_ylabel("v = sinθ·sinφ")
fig.colorbar(im, ax=axes[:2], shrink=0.8, label="дБ")
diff = maps[1] - maps[0]
m = (maps[0] > -20)
im2 = axes[2].imshow(np.where(m, diff, np.nan), extent=[-1, 1, -1, 1], origin="lower", vmin=-3, vmax=3, cmap="coolwarm")
axes[2].set_title("CST − формула, дБ" + NL + "(где ДН выше −20 дБ)", fontsize=9); axes[2].set_xlabel("u"); axes[2].set_ylabel("v")
fig.colorbar(im2, ax=axes[2], shrink=0.8, label="дБ")
rms = np.sqrt(np.nanmean(diff[m] ** 2))
fig.suptitle(f"{tag}: ДН виртуальной решётки на полусфере, {F0:g} ГГц. "
             f"Расхождение CST − формула по лепесткам выше −20 дБ: {rms:.2f} дБ (скз)", fontsize=10)
fig.savefig(os.path.join(OUT, f"{tag}_map.png"), dpi=170, bbox_inches="tight"); plt.close(fig)

# --- 3. пример: фаза встроенной ДН одного элемента и модель точечного источника в его центре
n0, x0, y0 = elems[0]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
for phi0, c in ((0, "tab:blue"), (90, "tab:red")):
    for side, sgn in ((phi0, 1), ((phi0 + 180) % 360, -1)):
        s = (ph == side) & (th <= 60)
        t = sgn * th[s]; o = np.argsort(t)
        phase = np.degrees(np.angle(pats[n0][s]))[o]
        model = np.degrees(np.angle(np.exp(1j * (k * (x0 * u[s] + y0 * v[s]) + np.radians(fits[n0][2])))))[o]
        ax[0].plot(t[o], phase, ".", color=c, ms=3, label=f"CST, φ = {phi0}°" if sgn > 0 else None)
        ax[0].plot(t[o], model, "-", color=c, lw=0.8, alpha=0.7,
                   label=f"точечный источник в ({x0:.0f}; {y0:.0f}) мм" if (sgn > 0 and phi0 == 0) else None)
        resid = np.degrees(np.angle(np.exp(1j * np.radians(phase - model))))
        ax[1].plot(t[o], resid, "-", color=c, lw=0.9, label=f"φ = {phi0}°" if sgn > 0 else None)
ax[0].set_title(f"{n0}: фаза встроенной ДН (от глобального нуля) и модель k·r·u", fontsize=9)
ax[0].set_xlabel("θ, °"); ax[0].set_ylabel("фаза, °"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
ax[1].set_title("остаток: фаза CST − модель точечного источника", fontsize=9); ax[1].set_xlabel("θ, °"); ax[1].set_ylabel("°")
ax[1].set_ylim(-40, 40); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{tag}_phase.png"), dpi=170); plt.close(fig)
json.dump({n: dict(xyz=list(map(float, f[0])), rms=float(f[1]), phi0=float(f[2])) for n, f in fits.items()},
          open(os.path.join(folder, "phase_centers_elems.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("ok", tag, f"dev={dev:.2f} mm, map rms={rms:.2f} dB")
