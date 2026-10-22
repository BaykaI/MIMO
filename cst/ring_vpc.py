"""Координаты виртуальных (эквивалентных) фазовых центров MIMO-кольца: расчёт двумя методами и сравнение.

Метод 1 (геометрический, концепция ВАР, [1] ф. (4), [3] с. 240): r_v(n,k) = r_tx,n + r_rx,k.
Метод 2 (электродинамический): по полноволновому отклику пары F_tx,n(u)·F_rx,k(u) из CST — МНК-подгонка фазы
отклика моделью точечного источника k·(r_v·u) в секторе ±30° вокруг нормали (tools/phase_center_fit.py).
Для всех N_tx × N_rx пар: таблица расхождений, рисунок «виртуальная решётка по формуле и по CST».

  python cst/ring_vpc.py cst/ring_R5 [--f 10]
Выход: <папка>/vpc_<f>GHz.json/.md, fig2/<папка>_vpc.png
"""
import sys, os, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from phase_center_fit import read_cst_farfield, fit_phase_center, C0

ap = argparse.ArgumentParser(); ap.add_argument("folder"); ap.add_argument("--f", type=float, default=10.0)
ap.add_argument("--sector", type=float, default=30.0)
a = ap.parse_args()
folder = a.folder; tag = os.path.basename(os.path.abspath(folder))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2"); os.makedirs(OUT, exist_ok=True)
LAM = C0 / (a.f * 1e9) * 1e3; k = 2 * np.pi / LAM
NL = chr(10)

lay = os.path.join(folder, "layout_mm.txt")
elems = [(p[0], float(p[1]), float(p[2])) for p in (l.split() for l in open(lay, encoding="utf-8")) if p and p[0] != "name"]
pos = {n: np.array([x, y]) for n, x, y in elems}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]


def load(name):
    path = os.path.join(folder, f"ff_{name}_{a.f:.1f}GHz.txt")
    th, ph, Et = read_cst_farfield(path, "theta"); _, _, Ep = read_cst_farfield(path, "phi")
    phr = np.radians(ph)
    return th, ph, Et * np.sin(phr) + Ep * np.cos(phr)


pats = {}
for n in pos:
    th, ph, pats[n] = load(n)

rows = []
for tn in tx:
    for rn in rx:
        geom = pos[tn] + pos[rn]
        xyz, rms = fit_phase_center(th, ph, pats[tn] * pats[rn], k, a.sector, (0.0, 0.0))
        rows.append(dict(tx=tn, rx=rn, geom=list(map(float, geom)), cst=[float(xyz[0]), float(xyz[1])], z=float(xyz[2]),
                         dev=float(np.hypot(xyz[0] - geom[0], xyz[1] - geom[1])), rms=float(rms)))
dev = np.array([r["dev"] for r in rows]); rmsv = np.array([r["rms"] for r in rows])
summary = dict(f=a.f, n_pairs=len(rows), dev_max=float(dev.max()), dev_mean=float(dev.mean()), dev_median=float(np.median(dev)),
               dev_p95=float(np.percentile(dev, 95)), rms_mean=float(rmsv.mean()), rms_max=float(rmsv.max()),
               n_over_lam20=int((dev > LAM / 20).sum()), n_over_lam50=int((dev > LAM / 50).sum()))
json.dump(dict(summary=summary, pairs=rows), open(os.path.join(folder, f"vpc_{a.f:g}GHz.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

worst = sorted(rows, key=lambda r: -r["dev"])[:5]
lines = [f"# Виртуальные фазовые центры {tag}, {a.f:g} ГГц, сектор ±{a.sector:g}°", "",
         f"Пар: {len(rows)}. Расхождение CST − формула: среднее {dev.mean():.2f} мм, медиана {np.median(dev):.2f}, "
         f"95 % ≤ {np.percentile(dev, 95):.2f}, максимум {dev.max():.2f} мм (λ/{LAM/dev.max():.0f}). "
         f"Пар с расхождением > λ/50: {summary['n_over_lam50']}, > λ/20: {summary['n_over_lam20']}. Остаток фазы {rmsv.mean():.1f}° в среднем.", "",
         "| пара | формула x; y, мм | CST x; y, мм | расхождение, мм | остаток, ° |", "|---|---|---|---|---|"]
for r in worst:
    lines.append(f"| {r['tx']}·{r['rx']} | {r['geom'][0]:.1f}; {r['geom'][1]:.1f} | {r['cst'][0]:.2f}; {r['cst'][1]:.2f} | {r['dev']:.2f} | {r['rms']:.1f} |")
txt = NL.join(lines).replace(".", ",")
print(txt)
open(os.path.join(folder, f"vpc_{a.f:g}GHz.md"), "w", encoding="utf-8").write(txt + NL)

# рисунок: виртуальная решётка по формуле (контуры) и по CST (точки, цвет — расхождение) + гистограмма
fig, ax = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={"width_ratios": [1.45, 1]})
g = np.array([r["geom"] for r in rows]); c = np.array([r["cst"] for r in rows])
ax[0].plot([pos[n][0] for n in rx], [pos[n][1] for n in rx], "o", mfc="none", mec="0.6", ms=6, label="Rx (физические)")
ax[0].plot([pos[n][0] for n in tx], [pos[n][1] for n in tx], "s", mfc="none", mec="0.6", ms=6, label="Tx (физические)")
ax[0].plot(g[:, 0], g[:, 1], "^", mfc="none", mec="0.25", ms=9, mew=1.0, label="виртуальные центры по формуле r_tx + r_rx")
sc = ax[0].scatter(c[:, 0], c[:, 1], c=dev, cmap="viridis", s=16, zorder=3, vmin=0, vmax=max(1.5, dev.max()), label="виртуальные центры по CST")
cb = fig.colorbar(sc, ax=ax[0], shrink=0.85); cb.set_label("расхождение CST − формула, мм")
ax[0].set_aspect("equal"); ax[0].grid(alpha=0.3); ax[0].set_xlabel("x, мм"); ax[0].set_ylabel("y, мм")
ax[0].set_title(f"{tag}: {len(rows)} виртуальных фазовых центров ({len(tx)} Tx × {len(rx)} Rx)", fontsize=10)
ax[0].legend(fontsize=7, loc="upper right")
ax[1].hist(dev, bins=25, color="0.5")
ax[1].axvline(LAM / 50, color="tab:blue", ls="--", lw=1); ax[1].text(LAM / 50, ax[1].get_ylim()[1] * 0.9, " λ/50", color="tab:blue", fontsize=8)
ax[1].axvline(LAM / 20, color="r", ls="--", lw=1); ax[1].text(LAM / 20, ax[1].get_ylim()[1] * 0.9, " λ/20", color="r", fontsize=8)
ax[1].set_xlabel("расхождение координат CST − формула, мм"); ax[1].set_ylabel("число пар")
ax[1].set_title(f"среднее {dev.mean():.2f} мм, максимум {dev.max():.2f} мм", fontsize=10); ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{tag}_vpc.png"), dpi=170); plt.close(fig)
print("fig", os.path.join(OUT, f"{tag}_vpc.png"))
