"""Рисунки деградации (оформление gost_plot) -> fig5/*.png
  python cst/degr_figs.py
"""
import os, json
import numpy as np
import gost_plot as G

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig5"); os.makedirs(OUT, exist_ok=True)
D = json.load(open(os.path.join(HERE, "degr.json"), encoding="utf-8"))
T = json.load(open(os.path.join(HERE, "degr_thr.json"), encoding="utf-8"))
Z = dict(np.load(os.path.join(HERE, "degr.npz")))
SP = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0]
db = lambda x: 20 * np.log10(x / x.max() + 1e-9)
lab = lambda d: f"{d:g}".replace(".", ",") + "λ"


def draw(a, curves, styles, marker=True):
    ls = []
    for (x, y), st in zip(curves, styles):
        ls += a.plot(x, y, **G.style(st, marker))
    return ls


# b_loss: потеря пика при узкополосной обработке от полосы, 6 шагов (формулы) + CST точки
fig, ax = G.figure(h_mm=80); a = ax[0, 0]
ls = draw(a, [([r["rel"] * 100 for r in T[f"{d:g}"]["rows"]], [r["loss_db"] for r in T[f"{d:g}"]["rows"]]) for d in SP], range(6), marker=False)
for d, mk in ((0.5, "o"), (1.5, "s")):
    rows = D["band"]["cst"][f"{d:g}"]["rows"]
    a.plot([r["rel"] * 100 for r in rows], [r["nb"]["loss_db"] for r in rows], ls="none", marker=mk, color="k", mfc="white", ms=4)
for y in (-1, -3):
    a.axhline(y, color="k", lw=0.6, ls=(0, (1, 3)))
a.set_xlabel("Относительная полоса, %"); a.set_ylabel("Потеря пика, дБ"); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "b_loss.png"))

# b_thr: пороги 1 и 3 дБ от шага и закон λ/ΔL
fig, ax = G.figure(h_mm=66); a = ax[0, 0]
x = SP; lam = 29.98
ls = draw(a, [(x, [T[f"{d:g}"]["thr1"] * 100 for d in x]), (x, [T[f"{d:g}"]["thr3"] * 100 for d in x])], (0, 1))
xx = np.linspace(0.3, 2, 100); dL = np.interp(xx, x, [T[f"{d:g}"]["dL_lam"] for d in x])
ls += a.plot(xx, 100 / dL, color="0.5", lw=0.9, ls="-")
a.set_xlabel("Шаг приёмников по дуге d/λ"); a.set_ylabel("Пороговая полоса, %"); a.set_ylim(bottom=0); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "b_thr.png"))

# b_sll: УБЛ по азимуту от полосы, многочастотная и узкополосная обработка, CST 0,5λ и 1,5λ
fig, ax = G.figure(2, 1, h_mm=72)
for col, d in enumerate((0.5, 1.5)):
    a = ax[0, col]; rows = D["band"]["cst"][f"{d:g}"]["rows"]; xr = [r["rel"] * 100 for r in rows]
    ls = draw(a, [(xr, [r["wb"]["sll_az"] for r in rows]), (xr, [r["nb"]["sll_az"] for r in rows])], (0, 1))
    a.set_xlabel("Относительная полоса, %"); a.set_ylabel("УБЛ по азимуту, дБ"); G.comma(a); G.key(a, ls); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "b_sll.png"), letters=True)

# b_res: разрешение по дальности от полосы, 1,5λ CST: wb, nb, оценка c/(2B)
fig, ax = G.figure(h_mm=64); a = ax[0, 0]; rows = [r for r in D["band"]["cst"]["1.5"]["rows"] if r["nf"] > 1]; xr = [r["rel"] * 100 for r in rows]
ls = draw(a, [(xr, [r["wb"]["r_res"] for r in rows]), (xr, [r["nb"]["r_res"] for r in rows]), (xr, [149.9 / (r["rel"] * 10) for r in rows])], (0, 1, 2))
a.set_yscale("log"); G.no_minor_labels(a); a.set_yticks([10, 20, 50, 100]); a.set_yticklabels(["10", "20", "50", "100"])
a.set_xlabel("Относительная полоса, %"); a.set_ylabel("Разрешение по дальности, мм"); G.comma(a, "x"); G.key(a, ls)
G.save(fig, os.path.join(OUT, "b_res.png"))

# b_az: сечения по азимуту, CST 1,5λ, полоса 0/20/60/100 %: а) многочастотная; б) узкополосная
az = Z["az"]
fig, ax = G.figure(2, 1, h_mm=72)
for col, mode in enumerate(("wb", "nb")):
    a = ax[0, col]; ref = Z[f"c1.5_r0_{mode}_az"].max()
    cur = []
    for r in (0, 20, 60, 100):
        y = Z[f"c1.5_r{r}_{mode}_az"]; ref_w = Z[f"c1.5_r{r}_wb_az"].max()
        cur.append((az, 20 * np.log10(y / ref_w + 1e-9)))
    ls = draw(a, cur, range(4), marker=False)
    a.set_xlim(0, 360); a.set_xticks(range(0, 361, 90)); a.set_ylim(-30, 1); a.set_xlabel("Азимут, град."); a.set_ylabel("Ψ, дБ")
    G.comma(a); G.key(a, ls); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "b_az.png"), letters=True)

# b_rng: профиль по дальности, CST 1,5λ, полоса 20 и 60 %, wb и nb
dr = Z["dr"]; m = np.abs(dr) <= 200
fig, ax = G.figure(2, 1, h_mm=70)
for col, r in enumerate((20, 60)):
    a = ax[0, col]; ref = Z[f"c1.5_r{r}_wb_rng"].max()
    ls = draw(a, [(dr[m], 20 * np.log10(Z[f"c1.5_r{r}_wb_rng"][m] / ref + 1e-9)), (dr[m], 20 * np.log10(Z[f"c1.5_r{r}_nb_rng"][m] / ref + 1e-9))], (0, 1), marker=False)
    a.set_ylim(-30, 1); a.set_xlabel("Дальность относительно цели, мм"); a.set_ylabel("Ψ, дБ"); G.comma(a); G.key(a, ls); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "b_rng.png"), letters=True)

# ---------------- шаг: сетки 3 × 2
LET = ["а)", "б)", "в)", "г)", "д)", "е)"]


def grid(key_fn, xs, xlabel, ylabel, ylim, fname, xlim=None, xt=None):
    fig, ax = G.figure(2, 3, h_mm=190)
    for n_, d in enumerate(SP):
        a = ax[n_ // 2, n_ % 2]
        ls = draw(a, key_fn(d), range(len(key_fn(d))), marker=False)
        if ylim: a.set_ylim(*ylim)
        if xlim: a.set_xlim(*xlim)
        if xt is not None: a.set_xticks(xt)
        a.set_xlabel(xlabel); a.set_ylabel(ylabel); G.comma(a)
        if n_ == 0: G.key(a, ls)
        G.letter(a, LET[n_] + " " + lab(d), y=-0.42)
    fig.subplots_adjust(hspace=0.75)
    G.save(fig, os.path.join(OUT, fname), letters=True)


ar = Z["az_resid"]
grid(lambda d: [(ar, Z[f"s{d:g}_dipole_resid"]), (ar, Z[f"s{d:g}_point_resid"])], ar, "Азимут, град.", "Остаток фазы, град.", (-180, 180), "s_resid.png",
     xlim=(0, 360), xt=range(0, 361, 90))
grid(lambda d: [(az, db(Z[f"s{d:g}_formula_img"])), (az, db(Z[f"s{d:g}_dipole_img"]))], az, "Азимут, град.", "Ψ, дБ", (-30, 1), "s_img.png",
     xlim=(0, 360), xt=range(0, 361, 90))
fi = Z["fine"]
grid(lambda d: [(fi, db(Z[f"s{d:g}_formula_df"])), (fi, db(Z[f"s{d:g}_dipole_df"]))], fi, "Азимут, град.", "Спектр, дБ", (-30, 1), "s_df.png",
     xlim=(0, 360), xt=range(0, 361, 90))
print("ok")
