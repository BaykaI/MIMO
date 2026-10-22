"""Фазовый центр MIMO-решётки в целом (а не отдельных элементов).

Отклик MIMO-РЛС в дальней зоне на направление u — сумма откликов всех пар: F(u) = Σ_n Σ_k F_tx,n(u)·F_rx,k(u).
Фазовый центр решётки — точка r_c, относительно которой фаза F(u) в главном лепестке не зависит от u.
Концепция ВАР предсказывает: r_c = r_c(Tx-подрешётки) + r_c(Rx-подрешётки) (центр виртуальной решётки),
а НЕ центр тяжести физических элементов. Скрипт находит r_c по полноволновым ДН из CST (МНК по фазе в
главном лепестке) и сравнивает с обоими кандидатами — для всей решётки и для несимметричных наборов каналов
(часть Tx), где кандидаты различаются сильнее всего.

  python cst/mimo_pc.py cst/stage3            # линейная: вся решётка, левая пара Tx, один Tx
  python cst/mimo_pc.py cst/ring_R5           # кольцо: вся решётка, половина Tx, один Tx
Выход: <папка>/mimo_pc.md, mimo_pc.json, fig2/<папка>_mimo_pc.png (фаза отклика вдоль сечения и две модели).
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
from phase_center_fit import read_cst_farfield, fit_phase_center, C0

folder = sys.argv[1]
tag = os.path.basename(os.path.abspath(folder))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2"); os.makedirs(OUT, exist_ok=True)
FREQS = (9.8, 10.0, 10.2)
NL = chr(10)

lay = os.path.join(folder, "layout_mm.txt")
if os.path.exists(lay):
    elems = [(p[0], float(p[1]), float(p[2])) for p in (l.split() for l in open(lay, encoding="utf-8")) if p and p[0] != "name"]
else:
    elems = [(f"rx{i+1}", x, 0.0) for i, x in enumerate([-105, -75, -45, -15, 15, 45, 75, 105])] + \
            [(f"tx{i+1}", x, 30.0) for i, x in enumerate([-127.5, -112.5, 112.5, 127.5])]
pos = {n: np.array([x, y]) for n, x, y in elems}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]
ring = os.path.exists(lay)

# наборы каналов: имя, список Tx, список Rx, сектор фита (внутри главного лепестка), сечение для рисунка
if ring:
    sets = [("вся решётка: 10 Tx × 20 Rx", tx, rx, 2.0, 0.0),
            ("половина Tx (tx1…tx5, правая) × 20 Rx", tx[:5], rx, 2.0, 0.0),
            ("один Tx (tx1) × 20 Rx", tx[:1], rx, 3.0, 0.0),
            ("один Tx (tx1) × один Rx (rx11, напротив)", tx[:1], rx[10:11], 30.0, 0.0)]
else:
    sets = [("вся решётка: 4 Tx × 8 Rx", tx, rx, 3.0, 90.0),
            ("левая пара Tx (tx1, tx2) × 8 Rx", tx[:2], rx, 3.0, 0.0),
            ("один Tx (tx1) × 8 Rx", tx[:1], rx, 3.0, 0.0),
            ("один Tx (tx1) × один Rx (rx8)", tx[:1], rx[7:8], 30.0, 0.0)]


def load(name, f):
    path = os.path.join(folder, f"ff_{name}_{f:.1f}GHz.txt")
    th, ph, Et = read_cst_farfield(path, "theta"); _, _, Ep = read_cst_farfield(path, "phi")
    phr = np.radians(ph)
    return th, ph, Et * np.sin(phr) + Ep * np.cos(phr)


rows = []
figdata = []
for f in FREQS:
    lam = C0 / (f * 1e9) * 1e3; k = 2 * np.pi / lam
    pats = {}
    for n in pos:
        th, ph, pats[n] = load(n, f)
    for sname, tset, rset, sector, cutphi in sets:
        F = np.zeros(len(th), complex)
        for tn in tset:
            for rn in rset:
                F += pats[tn] * pats[rn]
        c_var = np.mean([pos[n] for n in tset], axis=0) + np.mean([pos[n] for n in rset], axis=0)   # ВАР
        c_phys = np.mean([pos[n] for n in tset + rset], axis=0)                                      # центр тяжести элементов
        xyz, rms, ph0 = fit_phase_center(th, ph, F, k, sector, (0.0, 0.0), fix_z=(sector < 10), with_phase=True)
        rows.append(dict(f=f, set=sname, sector=sector, fit=[float(xyz[0]), float(xyz[1])], var=list(map(float, c_var)),
                         phys=list(map(float, c_phys)), dev_var=float(np.abs(xyz[:2] - c_var).max()),
                         dev_phys=float(np.abs(xyz[:2] - c_phys).max()), rms=float(rms)))
        if f == 10.0:
            figdata.append((sname, F.copy(), c_var, c_phys, xyz, ph0, cutphi, sector, k, th, ph))

# ---------------------------------------------------------------- таблица
lines = ["| f, ГГц | набор каналов | сектор фита | центр по CST x; y, мм | предсказание ВАР x; y | центр тяжести элементов x; y | откл. от ВАР, мм | откл. от центра тяжести, мм |",
         "|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append(f"| {r['f']:.1f} | {r['set']} | ±{r['sector']:g}° | {r['fit'][0]:.2f}; {r['fit'][1]:.2f} | {r['var'][0]:.1f}; {r['var'][1]:.1f} | "
                 f"{r['phys'][0]:.1f}; {r['phys'][1]:.1f} | {r['dev_var']:.2f} | {r['dev_phys']:.2f} |")
txt = NL.join(lines).replace(".", ",")
print(txt)
open(os.path.join(folder, "mimo_pc.md"), "w", encoding="utf-8").write(txt + NL)
json.dump(rows, open(os.path.join(folder, "mimo_pc.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ---------------------------------------------------------------- рисунок: фаза отклика решётки вдоль сечения и две модели
fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), squeeze=False)
for ax, (sname, F, c_var, c_phys, xyz, ph0, cutphi, sector, k, th, ph) in zip(axes.ravel(), figdata):
    span = max(sector, 3.0)
    sel1 = (ph == cutphi) & (th <= span); sel2 = (ph == (cutphi + 180) % 360) & (th <= span) & (th > 0)
    t = np.concatenate([-th[sel2], th[sel1]]); o = np.argsort(t); t = t[o]
    Fc = np.concatenate([F[sel2], F[sel1]])[o]
    ux = np.sin(np.radians(t)) * np.cos(np.radians(cutphi)); uy = np.sin(np.radians(t)) * np.sin(np.radians(cutphi))
    phase = np.degrees(np.angle(Fc))
    ref = phase[np.argmin(np.abs(t))]
    amp = 20 * np.log10(np.abs(Fc) / np.abs(Fc).max() + 1e-12)
    for c, lab, st in ((c_var, f"ВАР: ({c_var[0]:.0f}; {c_var[1]:.0f}) мм", "-"), (c_phys, f"центр тяжести: ({c_phys[0]:.0f}; {c_phys[1]:.0f}) мм", "--")):
        model = np.degrees(k * (c[0] * ux + c[1] * uy)); model -= model[np.argmin(np.abs(t))]
        ax.plot(t, np.degrees(np.angle(np.exp(1j * np.radians(model)))), st, lw=1.2, label=lab)
    ax.plot(t, np.degrees(np.angle(np.exp(1j * np.radians(phase - ref)))), "k.", ms=4, label="CST: фаза Σ F_tx·F_rx")
    ax2 = ax.twinx(); ax2.fill_between(t, -60, amp, color="0.85", zorder=0); ax2.set_ylim(-60, 0); ax2.set_yticks([])
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    ax.set_xlim(-span, span); ax.set_ylim(-180, 180); ax.grid(alpha=0.3)
    ax.set_title(sname + NL + f"найдено по CST: ({xyz[0]:.1f}; {xyz[1]:.1f}) мм, сечение φ = {cutphi:g}°", fontsize=9)
    ax.set_xlabel("θ, °"); ax.set_ylabel("фаза отклика, °"); ax.legend(fontsize=7, loc="lower right")
fig.suptitle(f"{tag}: фаза суммарного отклика MIMO-решётки в главном лепестке (точки) и её предсказание из двух точек "
             f"(серым — амплитуда лепестка), 10 ГГц", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{tag}_mimo_pc.png"), dpi=170); plt.close(fig)
print("fig", os.path.join(OUT, f"{tag}_mimo_pc.png"))
