"""Оформление графиков по шаблону отчётов (СТО МГТУ / ГОСТ 2.105):
чёрно-белая гамма (чёрный и серые, различие типом линии и маркером), Times New Roman 11, ширина 165 мм,
подписи осей с прописной буквы и с единицами через запятую, десятичная запятая, без заголовков на поле графика
(название идёт в подрисуночную подпись), кривые помечаются номерами 1, 2, 3 прямо на графике и расшифровываются
в подписи «Рисунок N – …: 1 – …; 2 – …». Подрисунки помечаются буквами а), б) под осью.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

MM = 1 / 25.4
W_MM = 165

STYLES = [dict(color="k", ls="-", lw=1.4), dict(color="k", ls="--", lw=1.3), dict(color="0.5", ls="-", lw=1.4),
          dict(color="k", ls=":", lw=1.6), dict(color="0.5", ls="--", lw=1.3), dict(color="0.35", ls="-.", lw=1.3)]
MARKERS = ["o", "s", "^", "v", "D", "x"]


def setup():
    plt.rcParams.update({
        "font.family": "Times New Roman", "font.size": 11, "mathtext.fontset": "stix",
        "axes.grid": True, "grid.color": "0.82", "grid.linewidth": 0.5, "axes.linewidth": 0.8,
        "axes.unicode_minus": True, "xtick.direction": "in", "ytick.direction": "in",
        "legend.frameon": False, "savefig.dpi": 200, "lines.markersize": 4.5})


def comma(ax, axis="both", nd=None):
    def fmt(v, _):
        if nd is None:
            s = f"{v:g}"
        else:
            s = f"{v:.{nd}f}"
        return s.replace(".", ",").replace("-", "−")
    if axis in ("x", "both"):
        ax.xaxis.set_major_formatter(FuncFormatter(fmt))
    if axis in ("y", "both"):
        ax.yaxis.set_major_formatter(FuncFormatter(fmt))


def figure(n_cols=1, n_rows=1, h_mm=70, w_mm=W_MM, **kw):
    setup()
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(w_mm * MM, h_mm * MM), squeeze=False, **kw)
    return fig, ax


def style(i, marker=False):
    s = dict(STYLES[i % len(STYLES)])
    if marker:
        s["marker"] = MARKERS[i % len(MARKERS)]; s["mfc"] = "white" if i % 2 else s["color"]
    return s


def num(ax, x, y, text, dx=0, dy=0):
    """номер кривой на графике (в белой подложке, чтобы читался поверх сетки)"""
    ax.annotate(str(text), (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=11, ha="center", va="center",
                bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="none"))


def num_on(ax, xs, ys, at_x, text, dx=0, dy=8):
    """номер у кривой в точке с абсциссой at_x"""
    xs = np.asarray(xs); ys = np.asarray(ys)
    i = int(np.nanargmin(np.abs(xs - at_x)))
    num(ax, xs[i], ys[i], text, dx, dy)


def letter(ax, s, y=-0.36):
    """подпись подрисунка а), б) под подписью оси абсцисс"""
    ax.text(0.5, y, s, transform=ax.transAxes, ha="center", va="top", fontsize=11)


def save(fig, path, letters=False):
    fig.tight_layout(h_pad=1.8 if letters else 1.0, w_pad=1.5)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def callout(ax, x, y, text, tx, ty):
    """номер кривой с тонкой выноской: точка (x, y) на кривой, текст в (tx, ty) в координатах данных"""
    ax.annotate(str(text), xy=(x, y), xytext=(tx, ty), fontsize=11, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="k", shrinkA=0, shrinkB=0),
                bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))


def callout_on(ax, xs, ys, at_x, text, tx, ty):
    xs = np.asarray(xs); ys = np.asarray(ys); i = int(np.nanargmin(np.abs(xs - at_x)))
    callout(ax, xs[i], ys[i], text, tx, ty)


def no_minor_labels(ax):
    from matplotlib.ticker import NullFormatter
    ax.yaxis.set_minor_formatter(NullFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())


def auto_label(ax, curves, labels=None, frac=0.14):
    """Номера кривых с выносками. curves: список (x, y). Для каждой кривой берётся точка, где она дальше всего
    от остальных (в долях осей), номер ставится выше или ниже на frac высоты оси, в сторону большего зазора."""
    fig = ax.figure; fig.canvas.draw()
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    logy = ax.get_yscale() == "log"
    def ny(v):
        v = np.asarray(v, float)
        return (np.log10(v) - np.log10(y0)) / (np.log10(y1) - np.log10(y0)) if logy else (v - y0) / (y1 - y0)
    def dy_back(t):
        return 10 ** (np.log10(y0) + t * (np.log10(y1) - np.log10(y0))) if logy else y0 + t * (y1 - y0)
    labels = labels or list(range(1, len(curves) + 1))
    grid = np.linspace(x0 + 0.06 * (x1 - x0), x1 - 0.06 * (x1 - x0), 60)
    used = []
    for ci, (x, y) in enumerate(curves):
        x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(y)
        if m.sum() < 1:
            continue
        best = None
        for gx in grid:
            if gx < x[m].min() or gx > x[m].max():
                continue
            yy = ny(np.interp(gx, x[m], y[m]))
            if not (0.04 < yy < 0.96):
                continue
            others = []
            for cj, (xo, yo) in enumerate(curves):
                if cj == ci:
                    continue
                xo = np.asarray(xo, float); yo = np.asarray(yo, float); mo = np.isfinite(yo)
                if mo.sum() and xo[mo].min() <= gx <= xo[mo].max():
                    others.append(ny(np.interp(gx, xo[mo], yo[mo])))
            up = min([o - yy for o in others if o > yy] + [0.96 - yy]); dn = min([yy - o for o in others if o < yy] + [yy - 0.04])
            gap = max(up, dn)
            if gap < frac + 0.05:          # текст не помещается между кривыми: точка хуже
                gap *= 0.5
            penalty = min([abs((gx - ux) / (x1 - x0)) for ux in used] + [1.0])
            score = gap + 0.3 * min(penalty, 0.3)
            if best is None or score > best[0]:
                best = (score, gx, yy, 1 if up >= dn else -1, gap)
        if best is None:
            continue
        _, gx, yy, sgn, gap = best; used.append(gx)
        ty = yy + sgn * min(frac, max(gap - 0.05, 0.06))
        callout(ax, gx, dy_back(yy), labels[ci], gx, dy_back(ty))


def key(ax, lines, labels=None):
    """Ключ над полем графика (вне данных): образец линии и номер; расшифровка номеров в подрисуночной подписи."""
    labels = labels or [str(i) for i in range(1, len(lines) + 1)]
    ax.legend(lines, labels, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=len(lines), frameon=False,
              handlelength=2.6, columnspacing=1.6, handletextpad=0.5, borderaxespad=0.1, fontsize=11)
