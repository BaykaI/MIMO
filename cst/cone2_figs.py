"""Рисунки по конусу (оформление gost_plot) -> fig6/*.png
  python cst/cone2_figs.py
Широкополосный приём: cone_pt_wb, cone_dip_wb (шаг постоянный). Шаг: cone_pt_d*, cone_dip_d* (+ s = 1 из *_wb), полоса 9-11 ГГц.
"""
import os, json, glob
import numpy as np
import gost_plot as G

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig6"); os.makedirs(OUT, exist_ok=True)
db = lambda x: 20 * np.log10(x / x.max() + 1e-9)


def load(tag):
    p = os.path.join(HERE, tag, "c2_result.json")
    if not os.path.exists(p):
        return None
    r = json.load(open(p, encoding="utf-8")); r["cuts"] = dict(np.load(os.path.join(HERE, tag, "c2_cuts.npz"))); return r


def draw(a, curves, styles, marker=True):
    ls = []
    for (x, y), st in zip(curves, styles):
        ls += a.plot(x, y, **G.style(st, marker))
    return ls


def img(r, rel, src="cst"):
    return [x for x in r["image"] if abs(x["rel"] - rel) < 1e-6][0][src]


WP, WD = load("cone_pt_wb"), load("cone_dip_wb")
WB = [w for w in (WP, WD) if w]

# ---------- широкополосный приём
if WB:
    # w_vpc: смещение виртуальных центров от r_t + r_r по частоте (среднее и наибольшее по 252 парам)
    fig, ax = G.figure(h_mm=70); a = ax[0, 0]; cur = []
    for w in WB:
        f = [v["f"] for v in w["vpc"]]
        cur += [(f, [v["dev_mean"] for v in w["vpc"]]), (f, [v["dev_max"] for v in w["vpc"]])]
    ls = draw(a, cur, range(len(cur))); a.set_ylim(bottom=0)
    a.axhline(1.5, color="k", lw=0.6, ls=(0, (1, 3)))
    a.set_xlabel("Частота, ГГц"); a.set_ylabel("Смещение центра, мм"); G.comma(a); G.key(a, ls)
    G.save(fig, os.path.join(OUT, "w_vpc.png"))

    # w_xyz: координаты виртуального центра одной пары относительно r_t + r_r по частоте
    fig, ax = G.figure(len(WB), 1, h_mm=70)
    for col, w in enumerate(WB):
        a = ax[0, col]; pr = w["vpc_examples"]["pairs"][0]; dx = np.array(w["vpc_examples"]["dxyz"][pr]); f = [v["f"] for v in w["vpc"]]
        ls = draw(a, [(f, dx[:, 0]), (f, dx[:, 1]), (f, dx[:, 2])], (0, 1, 2))
        a.set_xlabel("Частота, ГГц"); a.set_ylabel("Отклонение координаты, мм"); G.comma(a); G.key(a, ls)
        if len(WB) > 1: G.letter(a, ("а)", "б)")[col])
    G.save(fig, os.path.join(OUT, "w_xyz.png"), letters=len(WB) > 1)

    # w_rng: профиль по дальности при разной полосе (точечные, CST)
    c = WP["cuts"]; dr = c["dr"]; m = np.abs(dr) <= 150
    fig, ax = G.figure(h_mm=66); a = ax[0, 0]
    ls = draw(a, [(dr[m], db(c[f"b{b}_cst_rng"])[m]) for b in (10, 20, 30, 50)], range(4), marker=False)
    a.set_ylim(-30, 1); a.set_xlabel("Дальность относительно цели, мм"); a.set_ylabel("Ψ, дБ"); G.comma(a); G.key(a, ls)
    G.save(fig, os.path.join(OUT, "w_rng.png"))

    # w_ang: сечения по азимуту и углу места при одной частоте и полосах 20 и 50 % (точечные, CST)
    fig, ax = G.figure(2, 1, h_mm=72)
    for col, (nm, x, lim) in enumerate((("az", c["az"], (-180, 180)), ("el", c["el"], (-60, 60)))):
        a = ax[0, col]
        ls = draw(a, [(x, db(c[f"b{b}_cst_{nm}"])) for b in (0, 20, 50)], range(3), marker=False)
        a.set_xlim(*lim); a.set_ylim(-30, 1); a.set_xlabel(("Азимут, град.", "Угол места, град.")[col]); a.set_ylabel("Ψ, дБ")
        if nm == "az": a.set_xticks(range(-180, 181, 90))
        G.comma(a); G.key(a, ls); G.letter(a, ("а)", "б)")[col])
    G.save(fig, os.path.join(OUT, "w_ang.png"), letters=True)

    # w_band: а) разрешение по дальности; б) УБЛ по азимуту и углу места; от полосы
    fig, ax = G.figure(2, 1, h_mm=72)
    a = ax[0, 0]; cur = []
    rows = [x for x in WP["image"] if x["nf"] > 1]; xb = [x["rel"] * 100 for x in rows]
    cur.append((xb, [149.9 / (x["rel"] * 10) for x in rows]))
    for w in WB:
        rr = [x for x in w["image"] if x["nf"] > 1]; cur.append(([x["rel"] * 100 for x in rr], [x["cst"]["r_res"] for x in rr]))
    ls = draw(a, cur, range(len(cur))); a.set_yscale("log"); G.no_minor_labels(a); a.set_yticks([20, 50, 100]); a.set_yticklabels(["20", "50", "100"])
    a.set_xlabel("Относительная полоса, %"); a.set_ylabel("Разрешение по дальности, мм"); G.comma(a, "x"); G.key(a, ls); G.letter(a, "а)")
    a = ax[0, 1]; cur = []
    for w in WB:
        xx = [x["rel"] * 100 for x in w["image"]]
        cur += [(xx, [x["cst"]["az_sll"] for x in w["image"]]), (xx, [x["cst"]["el_sll"] for x in w["image"]])]
    ls = draw(a, cur, range(len(cur)))
    a.set_xlabel("Относительная полоса, %"); a.set_ylabel("УБЛ, дБ"); G.comma(a); G.key(a, ls); G.letter(a, "б)")
    G.save(fig, os.path.join(OUT, "w_band.png"), letters=True)
    print("полоса: ok")

# ---------- шаг (полоса 9-11 ГГц)
SPD = [0.5, 0.75, 1, 1.5, 2]


def series(kind):
    out = []
    for d in SPD:
        r = load(f"cone_{kind}_d{d:g}")
        if r: out.append(r)
    w = WP if kind == "pt" else WD
    if w: out.append(w)
    return sorted(out, key=lambda r: r["d_rx_lam"])


SP, SDP = series("pt"), series("dip")
if SP:
    f10 = lambda r: [v for v in r["vpc"] if abs(v["f"] - 10) < 1e-6][0]
    BAND = 0.2
    # s_vpc: смещение центров от шага на 10 ГГц
    fig, ax = G.figure(h_mm=66); a = ax[0, 0]; cur = []
    for S_ in (SP, SDP):
        if S_:
            x = [r["d_rx_lam"] for r in S_]; cur += [(x, [f10(r)["dev_mean"] for r in S_]), (x, [f10(r)["dev_max"] for r in S_])]
    ls = draw(a, cur, range(len(cur))); a.set_yscale("log"); G.no_minor_labels(a); a.axhline(1.5, color="k", lw=0.6, ls=(0, (1, 3)))
    a.set_yticks([0.01, 0.1, 1, 10]); a.set_yticklabels(["0,01", "0,1", "1", "10"])
    a.set_xlabel("Шаг приёмников нижнего кольца d/λ"); a.set_ylabel("Смещение центра, мм"); G.comma(a, "x"); G.key(a, ls)
    G.save(fig, os.path.join(OUT, "s_vpc.png"))

    # s_beam: а) ширина луча; б) УБЛ; по азимуту и углу места, полоса 20 %
    fig, ax = G.figure(2, 1, h_mm=72)
    for col, key in enumerate(("bw", "sll")):
        a = ax[0, col]; cur = []
        x = [r["d_rx_lam"] for r in SP]
        cur += [(x, [img(r, BAND, "formula")[f"az_{key}"] for r in SP]), (x, [img(r, BAND, "formula")[f"el_{key}"] for r in SP])]
        if SDP:
            xd = [r["d_rx_lam"] for r in SDP]
            cur += [(xd, [img(r, BAND)[f"az_{key}"] for r in SDP]), (xd, [img(r, BAND)[f"el_{key}"] for r in SDP])]
        ls = draw(a, cur, range(len(cur)))
        if key == "bw":
            a.set_yscale("log"); G.no_minor_labels(a); a.set_yticks([2, 5, 10, 20, 50]); a.set_yticklabels(["2", "5", "10", "20", "50"])
        a.set_xlabel("Шаг d/λ"); a.set_ylabel(("Ширина луча, град.", "УБЛ, дБ")[col]); G.comma(a, "x" if key == "bw" else "both"); G.key(a, ls)
        G.letter(a, ("а)", "б)")[col])
    G.save(fig, os.path.join(OUT, "s_beam.png"), letters=True)

    # s_az / s_el: сечения для всех шагов, 3 × 2, формулы и CST (вибраторы, если есть, иначе точечные)
    LET = ["а)", "б)", "в)", "г)", "д)", "е)"]
    for nm, xlab, lim in (("az", "Азимут, град.", (-180, 180)), ("el", "Угол места, град.", (-90, 90))):
        S_ = SDP if len(SDP) == len(SP) else SP
        fig, ax = G.figure(2, 3, h_mm=190)
        for n_, r in enumerate(S_[:6]):
            a = ax[n_ // 2, n_ % 2]; c = r["cuts"]; x = c[nm]
            ls = draw(a, [(x, db(c[f"b20_formula_{nm}"])), (x, db(c[f"b20_cst_{nm}"]))], (0, 1), marker=False)
            a.set_xlim(*lim); a.set_ylim(-30, 1); a.set_xlabel(xlab); a.set_ylabel("Ψ, дБ"); G.comma(a)
            if nm == "az": a.set_xticks(range(-180, 181, 90))
            if n_ == 0: G.key(a, ls)
            G.letter(a, LET[n_] + " " + f"{r['d_rx_lam']:.2g}".replace(".", ",") + "λ", y=-0.42)
        fig.subplots_adjust(hspace=0.75)
        G.save(fig, os.path.join(OUT, f"s_{nm}.png"), letters=True)
    print("шаг: ok", [round(r["d_rx_lam"], 2) for r in SP], [round(r["d_rx_lam"], 2) for r in SDP])
