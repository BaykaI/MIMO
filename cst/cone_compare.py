"""Конус: расчёт Чапурского (ОФН по ф. (12), виртуальные центры r_tx + r_rx) против полноволнового расчёта CST.

CST даёт встроенные ДН F_n(u) всех 32 элементов (фаза от центра основания). Отклик пары n–k в дальней зоне —
F_n(u)·F_k(u). Обработка по Чапурскому (24)–(25) с опорной точкой на горизонте в направлении u0 (φ = 0):
    Ψ_cst(u) = |Σ_nk F_n(u)·F_k(u)·exp(−j·k·(r_n + r_k)·u0)| / max,
то есть опорные фазы — геометрические (как в РЛС), а сигналы — из электродинамики. Для сравнения: Ψ по ф. (12)
с изотропными элементами и та же формула × квадрат среднего множителя элемента.
Виртуальные центры пар: МНК-подгонка фазы F_n·F_k моделью k·(r·u) в полосе углов места ±30° на всех азимутах, все три координаты.

  python cst/cone_compare.py [--f 10]
Выход: cone/cone_compare.json/.md, fig2/cone_cuts.png, fig2/cone_vpc.png
"""
import sys, os, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools")); sys.path.insert(0, HERE)
from phase_center_fit import read_cst_farfield, fit_phase_center, C0
from chap_ofn import metrics_cut

ap = argparse.ArgumentParser(); ap.add_argument("--f", type=float, default=10.0); ap.add_argument("--sector", type=float, default=30.0)
a = ap.parse_args()
folder = os.path.join(HERE, "cone"); OUT = os.path.join(HERE, "fig2")
LAM = C0 / (a.f * 1e9) * 1e3; k = 2 * np.pi / LAM
NL = chr(10)

elems = [(p[0], float(p[1]), float(p[2]), float(p[3])) for p in (l.split() for l in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")) if p and p[0] != "name"]
pos = {n: np.array([x, y, z]) for n, x, y, z in elems}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]

pats = {}
for n in pos:
    th, ph, E = read_cst_farfield(os.path.join(folder, f"ff_{n}_{a.f:.1f}GHz.txt"), "theta")   # z-диполь: ко-поляризация E_θ
    pats[n] = E
thr, phr = np.radians(th), np.radians(ph)
u = np.column_stack([np.sin(thr) * np.cos(phr), np.sin(thr) * np.sin(phr), np.cos(thr)])
u0 = np.array([1.0, 0.0, 0.0])                                   # горизонт, φ = 0

# --- ОФН тремя способами
F_cst = np.zeros(len(th), complex); F_geo = np.zeros(len(th), complex)
for tn in tx:
    for rn in rx:
        rv = pos[tn] + pos[rn]
        F_cst += pats[tn] * pats[rn] * np.exp(-1j * k * rv @ u0)
        F_geo += np.exp(1j * k * ((u - u0) @ rv))
elem = np.mean([np.abs(pats[n]) for n in pos], axis=0)
F_geo_el = F_geo * (elem / elem.max()) ** 2


def cut_az(F):
    s = th == 90
    o = np.argsort(((ph[s] + 180) % 360) - 180)
    az = (((ph[s] + 180) % 360) - 180)[o]; return az, np.abs(F[s])[o] / np.abs(F).max()


def cut_el(F):
    s1 = (ph == 0)                                   # вертикальная плоскость φ = 0: угол места = 90° − θ
    el = 90 - th[s1]; v = F[s1]
    m = (el >= -89) & (el <= 89)
    o = np.argsort(el[m]); return el[m][o], np.abs(v[m])[o] / np.abs(F).max()


rows = {}
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for ax_, cutf, name, span in ((ax[0], cut_az, "азимут", 90), (ax[1], cut_el, "угол места", 60)):
    for F, lab, st in ((F_geo, "Чапурский, ф. (12): геометрия, изотропные", dict(color="0.6", lw=1.0)),
                       (F_geo_el, "ф. (12) × множитель элемента", dict(color="k", lw=1.2)),
                       (F_cst, "CST: Σ F_tx·F_rx с геометрическими опорными фазами", dict(color="tab:red", lw=1.0, ls="--"))):
        x_, v = cutf(F); m_ = (x_ >= -span) & (x_ <= span)
        ax_.plot(x_[m_], 20 * np.log10(v[m_] + 1e-9), label=lab, **st)
        rows[(name, lab)] = metrics_cut(x_[m_], v[m_])
    ax_.set_ylim(-40, 1); ax_.set_xlim(-span, span); ax_.grid(alpha=0.3); ax_.set_ylabel("Ψ, дБ"); ax_.set_xlabel(name + ", °")
    ax_.set_title(f"конус, {a.f:g} ГГц: ОФН, {name} (опорное направление — горизонт, φ = 0)", fontsize=9); ax_.legend(fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "cone_cuts.png"), dpi=170); plt.close(fig)

# --- виртуальные центры всех пар
vrows = []
band = np.abs(th - 90) <= a.sector                   # полоса углов места ±sector на всех азимутах: x, y — из азимутальной
for tn in tx:                                        # зависимости фазы, z — из угломестной (сектор-конус вокруг одного
    for rn in rx:                                    # направления плохо определяет координату вдоль луча)
        rv = pos[tn] + pos[rn]
        xyz, rms = fit_phase_center(th[band], ph[band], (pats[tn] * pats[rn])[band], k, 180.0, (90.0, 0.0))
        vrows.append(dict(tx=tn, rx=rn, geom=list(map(float, rv)), cst=list(map(float, xyz)), dev=float(np.linalg.norm(xyz - rv)), rms=float(rms)))
dev = np.array([r["dev"] for r in vrows]); dev_xyz = np.array([np.abs(np.array(r["cst"]) - np.array(r["geom"])) for r in vrows])
summ = dict(n_pairs=len(vrows), dev_mean=float(dev.mean()), dev_p95=float(np.percentile(dev, 95)), dev_max=float(dev.max()),
            dev_xyz_max=list(map(float, dev_xyz.max(axis=0))), rms_mean=float(np.mean([r["rms"] for r in vrows])),
            n_over_lam20=int((dev > LAM / 20).sum()))
g = np.array([r["geom"] for r in vrows]); c = np.array([r["cst"] for r in vrows])
fig, ax = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
for i, (p, q, lab) in enumerate(((0, 1, "вид сверху: x, y"), (0, 2, "вид сбоку: x, z"))):
    ax[i].plot(g[:, p], g[:, q], "^", mfc="none", mec="0.3", ms=8, label="r_tx + r_rx (Чапурский)")
    sc = ax[i].scatter(c[:, p], c[:, q], c=dev, cmap="viridis", s=12, zorder=3, vmin=0, vmax=max(1.5, dev.max()), label="по CST")
    ax[i].set_aspect("equal" if i == 0 else "auto"); ax[i].grid(alpha=0.3); ax[i].set_title(f"виртуальные центры {len(vrows)} пар, {lab}", fontsize=9)
    if i == 1:
        ax[i].set_ylim(-20, 240)
    ax[i].set_xlabel(["x", "y", "z"][p] + ", мм"); ax[i].set_ylabel(["x", "y", "z"][q] + ", мм"); ax[i].legend(fontsize=7)
fig.colorbar(sc, ax=ax[:2], shrink=0.8, label="расхождение, мм")
ax[2].hist(dev, bins=25, color="0.5"); ax[2].axvline(LAM / 20, color="r", ls="--"); ax[2].text(LAM / 20, ax[2].get_ylim()[1] * 0.9, " λ/20", color="r", fontsize=8)
ax[2].set_xlabel("расхождение CST − формула, мм"); ax[2].set_ylabel("число пар"); ax[2].grid(alpha=0.3)
ax[2].set_title(f"среднее {dev.mean():.2f} мм, максимум {dev.max():.2f} мм", fontsize=9)
fig.savefig(os.path.join(OUT, "cone_vpc.png"), dpi=170, bbox_inches="tight"); plt.close(fig)

lines = ["| сечение | модель | ширина −3 дБ, ° | УБЛ, дБ |", "|---|---|---|---|"]
for (name, lab), m_ in rows.items():
    lines.append(f"| {name} | {lab} | {m_['hpbw']:.2f} | {m_['sll_db']:.1f} |")
lines += ["", f"Виртуальные центры: пар {summ['n_pairs']}, расхождение среднее {summ['dev_mean']:.2f} мм, 95 % ≤ {summ['dev_p95']:.2f}, "
          f"макс. {summ['dev_max']:.2f} мм (по осям x/y/z: {summ['dev_xyz_max'][0]:.2f}/{summ['dev_xyz_max'][1]:.2f}/{summ['dev_xyz_max'][2]:.2f}); "
          f"пар > λ/20: {summ['n_over_lam20']}; остаток фазы {summ['rms_mean']:.1f}°"]
txt = NL.join(lines).replace(".", ","); print(txt)
open(os.path.join(folder, "cone_compare.md"), "w", encoding="utf-8").write(txt + NL)
json.dump(dict(cuts={f"{n}|{l}": m_ for (n, l), m_ in rows.items()}, vpc=summ, pairs=vrows),
          open(os.path.join(folder, "cone_compare.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
