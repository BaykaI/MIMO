"""Рисунки для раздельного исследования двух факторов (оформление gost_plot):
  фактор шага, одна частота 10 ГГц:  fig4/s_*.png   (spacing_1f.json, ringz_result.json)
  фактор полосы, точечные излучатели: fig4/w_*.png   (ringz_d*_ptwb/wb_result.json, wb_cuts.npz)
  python cst/factor_figs.py
"""
import os, json, glob
import numpy as np
import gost_plot as G

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig4"); os.makedirs(OUT, exist_ok=True)
XL = "Шаг приёмников по дуге d/λ"
db = lambda x: 20 * np.log10(x / x.max() + 1e-9)


def draw(a, curves, styles, marker=True):
    ls = []
    for (x, y), st in zip(curves, styles):
        ls += a.plot(x, y, **G.style(st, marker))
    return ls


# ================= фактор шага
S = json.load(open(os.path.join(HERE, "spacing_1f.json"), encoding="utf-8"))
fo = S["formula"]; dip = S["cst"]["dipole"]; pt = S["cst"]["point"]
xd = [r["d"] for r in dip]; xp = [r["d"] for r in pt]; xf = [r["d"] for r in fo]

# s_beam: ширина луча (формулы) MIMO и кольца Rx
fig, ax = G.figure(h_mm=66); a = ax[0, 0]
ls = draw(a, [(xf, [r["rx_bw"] for r in fo]), (xf, [r["mimo_bw"] for r in fo])], (0, 1), marker=False)
a.set_yscale("log"); G.no_minor_labels(a); a.set_yticks([2, 5, 10, 20, 50]); a.set_yticklabels(["2", "5", "10", "20", "50"])
a.set_xlabel(XL); a.set_ylabel("Ширина луча по азимуту, град."); G.comma(a, "x"); G.key(a, ls)
G.save(fig, os.path.join(OUT, "s_beam.png"))

# s_sll: а) MIMO-изображение, б) спектр пеленгатора (кольцо Rx), одна частота
fig, ax = G.figure(2, 1, h_mm=72)
a = ax[0, 0]
ls = draw(a, [(xf, [r["mimo_sll"] for r in fo])], (0,), marker=False)
ls += draw(a, [(xp, [r["cst_geo"]["sll"] for r in pt]), (xd, [r["cst_geo"]["sll"] for r in dip]), (xd, [r["cst_full"]["sll"] for r in dip])], (1, 2, 3))
for l in ls[1:2]:
    l.set_linestyle("none")
a.set_xlabel(XL); a.set_ylabel("УБЛ изображения, дБ"); G.comma(a); G.key(a, ls); G.letter(a, "а)")
a = ax[0, 1]
ls = draw(a, [(xf, [r["rx_sll"] for r in fo])], (0,), marker=False)
ls += draw(a, [(xp, [r["df_cst"]["sll_db"] for r in pt]), (xd, [r["df_cst"]["sll_db"] for r in dip])], (1, 2))
ls[1].set_linestyle("none")
a.set_xlabel(XL); a.set_ylabel("Ложный максимум, дБ"); G.comma(a); G.key(a, ls); G.letter(a, "б)")
G.save(fig, os.path.join(OUT, "s_sll.png"), letters=True)

# s_phase: а) фазовая ошибка формулы на 10 ГГц; б) ошибка пеленга
fig, ax = G.figure(2, 1, h_mm=70)
a = ax[0, 0]
ls = draw(a, [(xd, [r["v10"]["phase_err_rms_deg"] for r in dip]), (xd, [r["v10"]["global_rms_deg"] for r in dip]),
              (xp, [r["v10"]["phase_err_rms_deg"] for r in pt])], (0, 1, 2))
a.set_ylim(bottom=0); a.set_xlabel(XL); a.set_ylabel("Фазовая ошибка, град."); G.comma(a); G.key(a, ls); G.letter(a, "а)")
a = ax[0, 1]
ls = draw(a, [(xd, [r["df_cst"]["max"] for r in dip]), (xp, [r["df_cst"]["max"] for r in pt])], (0, 2))
a.set_ylim(bottom=0); a.set_xlabel(XL); a.set_ylabel("Ошибка пеленга, град."); G.comma(a); G.key(a, ls); G.letter(a, "б)")
G.save(fig, os.path.join(OUT, "s_phase.png"), letters=True)
print("шаг: ok")

# ================= фактор полосы
WB = []
for p in sorted(glob.glob(os.path.join(HERE, "ringz_d*_ptwb", "wb_result.json"))):
    r = json.load(open(p, encoding="utf-8")); r["cuts"] = dict(np.load(os.path.join(os.path.dirname(p), "wb_cuts.npz"))); WB.append(r)
WB.sort(key=lambda r: r["d_lam"])
if not WB:
    raise SystemExit("нет данных широкой полосы")

# w_pair: фаза пары по частоте в 4 направлениях (CST точки, формула линии), шаг из первого прогона
r = WB[0]; c = r["cuts"]; f = c["freqs"]
fig, ax = G.figure(h_mm=76); a = ax[0, 0]; ls = []
for n, j in enumerate((0, 3, 6, 9)):
    st = G.style(n, True)
    ls += a.plot(f, c["pair_phase_geo"][j], color=st.get("color", "k"), ls=st.get("ls", "-"), lw=st.get("lw", 1.0))
    a.plot(f, c["pair_phase_cst"][j], ls="none", marker=st.get("marker", "o"), color=st.get("color", "k"), ms=4, mfc="white")
a.set_xlabel("Частота, ГГц"); a.set_ylabel("Фаза отклика пары, град."); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "w_pair.png"))

# w_perr: фазовая ошибка формулы по частоте (точечные, оба шага)
fig, ax = G.figure(h_mm=62); a = ax[0, 0]
ls = draw(a, [([v["f"] for v in r["by_f"]], [v["phase_err_rms_deg"] for v in r["by_f"]]) for r in WB], range(len(WB)))
a.set_ylim(bottom=0); a.set_xlabel("Частота, ГГц"); a.set_ylabel("Фазовая ошибка, град."); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "w_perr.png"))

# w_band: а) разрешение по дальности; б) УБЛ изображения по азимуту; от полосы
fig, ax = G.figure(2, 1, h_mm=72)
a = ax[0, 0]; rows = [x for x in WB[0]["image_vs_band"] if x["nf"] > 1]
xb = [x["rel"] * 100 for x in rows]
cur = [(xb, [x["formula"]["range_res_mm"] for x in rows])] + [([x["rel"] * 100 for x in R["image_vs_band"] if x["nf"] > 1],
                                                               [x["cst"]["range_res_mm"] for x in R["image_vs_band"] if x["nf"] > 1]) for R in WB]
ls = draw(a, cur, range(len(cur))); a.set_yscale("log"); G.no_minor_labels(a)
a.set_yticks([10, 20, 50, 100]); a.set_yticklabels(["10", "20", "50", "100"])
a.set_xlabel("Относительная полоса, %"); a.set_ylabel("Разрешение по дальности, мм"); G.comma(a, "x"); G.key(a, ls); G.letter(a, "а)")
a = ax[0, 1]; cur = []
for R in WB:
    xx = [x["rel"] * 100 for x in R["image_vs_band"]]
    cur += [(xx, [x["formula"]["sll_az_db"] for x in R["image_vs_band"]]), (xx, [x["cst"]["sll_az_db"] for x in R["image_vs_band"]])]
ls = draw(a, cur, range(len(cur)))
a.set_xlabel("Относительная полоса, %"); a.set_ylabel("УБЛ по азимуту, дБ"); G.comma(a); G.key(a, ls); G.letter(a, "б)")
G.save(fig, os.path.join(OUT, "w_band.png"), letters=True)

# w_az: сечение изображения по азимуту, большой шаг: одна частота и полоса 100 %
R = WB[-1]; c = R["cuts"]; az = c["az"]
fig, ax = G.figure(h_mm=66); a = ax[0, 0]
ls = draw(a, [(az, db(c["band0_cst_az"])), (az, db(c["band30_cst_az"])), (az, db(c["band100_cst_az"]))], (0, 1, 2), marker=False)
a.set_xlim(0, 360); a.set_ylim(-30, 1); a.set_xticks(range(0, 361, 60)); a.set_xlabel("Азимут, град."); a.set_ylabel("Ψ, дБ"); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "w_az.png"))

# w_range: профиль по дальности (CST), полосы 10, 30, 100 %
fig, ax = G.figure(h_mm=62); a = ax[0, 0]; dr = c["dr"]; m = np.abs(dr) <= 150
ls = draw(a, [(dr[m], db(c[f"band{b}_cst_range"])[m]) for b in (10, 30, 100)], (0, 1, 2), marker=False)
a.set_ylim(-30, 1); a.set_xlabel("Дальность относительно цели, мм"); a.set_ylabel("Ψ, дБ"); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "w_range.png"))

# w_df: ложный максимум спектра пеленгатора от полосы
fig, ax = G.figure(h_mm=62); a = ax[0, 0]
ls = draw(a, [([x["rel"] * 100 for x in R["df_vs_band"]], [x["sll_db"] for x in R["df_vs_band"]]) for R in WB], range(len(WB)))
a.set_xlabel("Относительная полоса, %"); a.set_ylabel("Ложный максимум, дБ"); G.comma(a); G.key(a, ls)
G.save(fig, os.path.join(OUT, "w_df.png"))
print("полоса: ok", [r["folder"] for r in WB])
