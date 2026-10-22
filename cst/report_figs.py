"""Рисунки для отчёта: схема решётки и смещения/фазы элементов (этап 3, 10 ГГц)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "fig")
os.makedirs(OUT, exist_ok=True)
RX_X = [-105, -75, -45, -15, 15, 45, 75, 105]
TX_X = [-127.5, -112.5, 112.5, 127.5]

# 1. схема решётки + виртуальная решётка
fig, ax = plt.subplots(1, 2, figsize=(7.5, 2.6), gridspec_kw={"width_ratios": [1.3, 1]})
a = ax[0]
a.plot(RX_X, [0] * 8, "o", color="0.35", ms=7, label="Rx (8), y = 0")
a.plot(TX_X, [30] * 4, "s", color="tab:blue", ms=7, label="Tx (4), y = 30 мм")
for x in RX_X:
    a.plot([x, x], [-6.25, 6.25], "-", color="0.35", lw=1.2)
for x in TX_X:
    a.plot([x, x], [23.75, 36.25], "-", color="tab:blue", lw=1.2)
a.set_xlabel("x, мм"); a.set_ylabel("y, мм"); a.set_title("Решётка 8 Rx + 4 Tx (диполи вдоль y)", fontsize=9)
a.set_ylim(-20, 50); a.set_xlim(-145, 145); a.legend(fontsize=7, loc="upper center", ncol=2); a.grid(alpha=0.3)
a.set_aspect("equal")
b = ax[1]
vx = sorted(x + t for x in RX_X for t in TX_X)
b.plot(vx, [0] * len(vx), "d", color="tab:blue", ms=5)
b.set_yticks([]); b.set_xlabel("x_tx + x_rx, мм"); b.set_title("Виртуальная решётка: 32 элемента, шаг 15 мм", fontsize=9)
b.grid(alpha=0.3, axis="x")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "layout.png"), dpi=200); plt.close(fig)

# 2. смещения по x и фазы φ0 (этап 3, 12,5 мм; 10 ГГц)
rows = json.load(open(os.path.join(ROOT, "stage3", "phase_centers.json"), encoding="utf-8"))
names = [f"rx{i}" for i in range(1, 9)] + [f"tx{i}" for i in range(1, 5)]
el = {(r["f"], r["name"]): r for r in rows if r["kind"] == "элемент"}
fig, ax = plt.subplots(1, 2, figsize=(7.5, 2.7))
for f, c in ((9.8, "0.6"), (10.0, "tab:blue"), (10.2, "0.3")):
    dx = [el[(f, n)]["fit"][0] - el[(f, n)]["expected"][0] for n in names]
    ax[0].plot(range(12), dx, "o-", color=c, ms=4, lw=1, label=f"{f:.1f} ГГц".replace(".", ","))
ax[0].axhline(1.5, color="r", ls="--", lw=0.8); ax[0].axhline(-1.5, color="r", ls="--", lw=0.8)
ax[0].text(0.2, 1.6, "λ/20", color="r", fontsize=7)
ax[0].set_xticks(range(12)); ax[0].set_xticklabels(names, rotation=60, fontsize=7)
ax[0].set_ylabel("Δx фазового центра, мм"); ax[0].set_title("Смещение по оси решётки", fontsize=9)
ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
# φ0 относительно rx4 из phase_centers.md (пересчитать здесь нельзя — берём из таблицы)
md = open(os.path.join(ROOT, "stage3", "phase_centers.md"), encoding="utf-8").read().split("\n")
ph = {}
for line in md:
    if line.startswith("| 9,8 |") or line.startswith("| 10,0 |") or line.startswith("| 10,2 |"):
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 13:
            ph[float(cells[0].replace(",", "."))] = [float(c.replace(",", ".")) for c in cells[1:]]
for f, c in ((9.8, "0.6"), (10.0, "tab:blue"), (10.2, "0.3")):
    ax[1].plot(range(12), ph[f], "s-", color=c, ms=4, lw=1, label=f"{f:.1f} ГГц".replace(".", ","))
ax[1].set_xticks(range(12)); ax[1].set_xticklabels(names, rotation=60, fontsize=7)
ax[1].set_ylabel("φ0 относительно rx4, °"); ax[1].set_title("Фаза встроенной ДН элемента", fontsize=9)
ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "stage3_dx_phi.png"), dpi=200); plt.close(fig)
print("ok", os.listdir(OUT))
