"""Рисунки по результатам CST (серия колец ringz_d*, эталон ringz_d*_pt) в оформлении gost_plot.
  python cst/cst_figs.py   -> fig3/c_*.png, ringz_summary.json
Номера кривых на рисунках расшифровываются в подписях отчёта.
"""
import os, json, glob
import numpy as np
import gost_plot as G

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig3"); os.makedirs(OUT, exist_ok=True)


def load(pattern):
    out = []
    for p in glob.glob(os.path.join(HERE, pattern, "ringz_result.json")):
        r = json.load(open(p, encoding="utf-8")); r["cuts"] = dict(np.load(os.path.join(os.path.dirname(p), "cuts.npz"))); out.append(r)
    return sorted(out, key=lambda r: r["d_lam"])


runs = [r for r in load("ringz_d*") if not r["folder"].endswith("_pt")]
pts = load("ringz_d*_pt")
D = np.array([round(r["d_lam"], 2) for r in runs]); DP = np.array([round(r["d_lam"], 2) for r in pts])
form = json.load(open(os.path.join(HERE, "formula_band_spacing.json"), encoding="utf-8"))["spacing"]
KF, KG, KC, KM = "формулы (идеал)", "CST, обработка по геометрии", "CST, геометрия + калибровка каналов", "CST, полная модель CST"
XL = "Шаг приёмников по дуге d/λ"


def v10(r):
    return [v for v in r["vpc_by_f"] if abs(v["f"] - 10.0) < 1e-6][0]


def pick(lst, d):
    c = [r for r in lst if abs(r["d_lam"] - d) < 0.05]
    return c[0] if c else None


def db(x):
    return 20 * np.log10(x / x.max() + 1e-9)


def plot_curves(a, curves, styles, marker=True):
    lines = []
    for (x, y), st in zip(curves, styles):
        lines += a.plot(x, y, **G.style(st, marker))
    a._my_lines = lines
    return lines


summ = dict(dipole=[], point=[])
for key, lst in (("dipole", runs), ("point", pts)):
    for r in lst:
        summ[key].append(dict(d=float(r["d_lam"]), d_mm=r["d_mm"], R=r["R_rx"], coupling=r["coupling"], v10=v10(r),
                              perr_band=[v["phase_err_rms_deg"] for v in r["vpc_by_f"]], delay=r["delay"], df=r["df"],
                              image=r["image"], image_vs_band=r["image_vs_band"]))
json.dump(summ, open(os.path.join(HERE, "ringz_summary.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# 1. связь
fig, ax = G.figure(h_mm=68); a = ax[0, 0]
cur = [(D, [r["coupling"]["rx_rx_db"] for r in runs]), (D, [r["coupling"]["rx_tx_db"] for r in runs])]
if len(pts):
    cur.append((DP, [r["coupling"]["rx_rx_db"] for r in pts]))
plot_curves(a, cur, (0, 1, 2)); a.set_xlabel(XL); a.set_ylabel("Модуль коэффициента связи, дБ"); G.comma(a); G.key(a, a._my_lines)
G.save(fig, os.path.join(OUT, "c_coupling.png"))

# 2. фазовая ошибка (а) и смещение единого центра и центра задержки (б)
fig, ax = G.figure(2, 1, h_mm=70)
a = ax[0, 0]
cur = [(D, [v10(r)["phase_err_rms_deg"] for r in runs]), (D, [max(v["phase_err_rms_deg"] for v in r["vpc_by_f"]) for r in runs])]
if len(pts):
    cur.append((DP, [v10(r)["phase_err_rms_deg"] for r in pts]))
plot_curves(a, cur, (0, 1, 2)); a.set_xlabel(XL); a.set_ylabel("Фазовая ошибка, град."); a.set_ylim(bottom=0); G.comma(a); G.key(a, a._my_lines); G.letter(a, "а)")
a = ax[0, 1]
cur = [(D, [v10(r)["global_dev_mean"] for r in runs]), (D, [r["delay"]["dr_mean"] for r in runs])]
if len(pts):
    cur += [(DP, [v10(r)["global_dev_mean"] for r in pts]), (DP, [r["delay"]["dr_mean"] for r in pts])]
plot_curves(a, cur, (0, 1, 2, 4)); a.axhline(1.5, color="k", lw=0.6, ls=(0, (1, 3)))
a.set_yscale("log"); G.no_minor_labels(a); a.set_xlabel(XL); a.set_ylabel("Смещение центра, мм"); G.comma(a, "x"); G.key(a, a._my_lines); G.letter(a, "б)")
G.save(fig, os.path.join(OUT, "c_phase_err.png"), letters=True)

# 3. фазовая ошибка по частоте
fig, ax = G.figure(h_mm=72); a = ax[0, 0]
cur = [([v["f"] for v in r["vpc_by_f"]], [v["phase_err_rms_deg"] for v in r["vpc_by_f"]]) for r in runs]
plot_curves(a, cur, range(len(cur))); a.set_xlabel("Частота, ГГц"); a.set_ylabel("Фазовая ошибка, град."); a.set_ylim(bottom=0); G.comma(a); G.key(a, a._my_lines)
G.save(fig, os.path.join(OUT, "c_phase_vs_f.png"))

# 4. остаточная фаза пары tx1·rx1 по частоте в четырёх направлениях
fig, ax = G.figure(2, 1, h_mm=70)
for col, d in enumerate((0.5, 2.0)):
    r = pick(runs, d); a = ax[0, col]
    if r is None:
        continue
    c = r["cuts"]; f = c["freqs"]; P = c["pair_resid_phase"]
    Pn = P - P[:, [6]]                       # фаза относительно 10 ГГц
    Pm = Pn.mean(axis=0)                      # общая для всех направлений часть (общая задержка пары, вычитается)
    cur = [(f, Pn[j] - Pm) for j in (0, 3, 6, 9)]
    plot_curves(a, cur, range(4), marker=False); a.set_xlabel("Частота, ГГц"); a.set_ylabel("Остаточная фаза, град."); G.comma(a); G.key(a, a._my_lines)
    G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "c_pair_phase.png"), letters=True)

# 5. профиль по дальности через цель (полоса 30 %)
fig, ax = G.figure(2, 1, h_mm=70)
for col, d in enumerate((0.3, 1.5)):
    r = pick(runs, d); a = ax[0, col]
    if r is None:
        continue
    c = r["cuts"]; m = np.abs(c["dr"]) <= 200
    cur = [(c["dr"][m], db(c[k])[m]) for k in ("img_geo_geo_range", "img_geo_cst_range", "img_cal_cst_range")]
    plot_curves(a, cur, (0, 1, 2), marker=False); a.set_ylim(-30, 1); a.set_xlabel("Дальность относительно цели, мм"); a.set_ylabel("Ψ, дБ")
    G.comma(a); G.key(a, a._my_lines); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "c_range_profile.png"), letters=True)

# 6. разрешение по дальности от ширины полосы
fig, ax = G.figure(h_mm=68); a = ax[0, 0]
rows = [b for b in runs[0]["image_vs_band"] if b["nf"] > 1]; x = [b["rel"] * 100 for b in rows]
cur = [(x, [b["formula"]["range_res_mm"] for b in rows])]
for d in (0.3, 2.0):
    r = pick(runs, d)
    if r:
        rr = [b for b in r["image_vs_band"] if b["nf"] > 1]; cur.append((x, [b["cst_geo"]["range_res_mm"] for b in rr]))
plot_curves(a, cur, (0, 1, 2)); a.set_yscale("log"); G.no_minor_labels(a); a.set_yticks([30, 50, 100, 200]); a.set_yticklabels(["30", "50", "100", "200"])
a.set_xlabel("Относительная полоса сигнала, %"); a.set_ylabel("Разрешение по дальности, мм"); G.comma(a, "x"); G.key(a, a._my_lines)
G.save(fig, os.path.join(OUT, "c_range_vs_band.png"))

# 7. ошибка пеленга
fig, ax = G.figure(h_mm=68); a = ax[0, 0]
cur = [(D, [r["df"]["1 частота 10 ГГц"]["max"] for r in runs]), (D, [r["df"]["полоса 8,5-11,5 ГГц"]["max"] for r in runs])]
if len(pts):
    cur.append((DP, [r["df"]["1 частота 10 ГГц"]["max"] for r in pts]))
plot_curves(a, cur, (0, 1, 2)); a.set_xlabel(XL); a.set_ylabel("Ошибка пеленга, град."); a.set_ylim(bottom=0); G.comma(a); G.key(a, a._my_lines)
G.save(fig, os.path.join(OUT, "c_df_err.png"))

# 8. УБЛ изображения от шага
fig, ax = G.figure(h_mm=78); a = ax[0, 0]
cur = [(D, [r["image"][k]["sll_az_db"] for r in runs]) for k in (KF, KG, KC, KM)]
if len(pts):
    cur.append((DP, [r["image"][KG]["sll_az_db"] for r in pts]))
plot_curves(a, cur, (0, 1, 2, 3, 4)); a.set_xlabel(XL); a.set_ylabel("УБЛ изображения по азимуту, дБ"); G.comma(a); G.key(a, a._my_lines)
G.save(fig, os.path.join(OUT, "c_img_sll.png"))

# 9. сечение изображения по азимуту
fig, ax = G.figure(2, 1, h_mm=70)
for col, d in enumerate((0.3, 1.5)):
    r = pick(runs, d); a = ax[0, col]
    if r is None:
        continue
    c = r["cuts"]; m = (c["az"] >= -60) & (c["az"] <= 95)
    cur = [(c["az"][m], db(c[k])[m]) for k in ("img_geo_geo_az", "img_geo_cst_az", "img_cst_cst_az")]
    plot_curves(a, cur, (0, 1, 2), marker=False); a.set_ylim(-30, 1); a.set_xlabel("Азимут, град."); a.set_ylabel("Ψ, дБ")
    G.comma(a); G.key(a, a._my_lines); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "c_img_az.png"), letters=True)
print("ok", D, DP)

# 0. схема колец: шаг 0,5λ и 2λ (Rx — кружки, Tx — квадраты), с виртуальными точками r_tx + r_rx
fig, ax = G.figure(2, 1, h_mm=78)
for col, d in enumerate((0.5, 2.0)):
    r = pick(runs, d); a = ax[0, col]
    if r is None:
        continue
    lay = [l.split() for l in open(os.path.join(HERE, r["folder"], "layout_mm.txt"), encoding="utf-8")][1:]
    rx = np.array([[float(p[1]), float(p[2])] for p in lay if p[0].startswith("rx")]); tx = np.array([[float(p[1]), float(p[2])] for p in lay if p[0].startswith("tx")])
    v = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2)
    l3 = a.plot(v[:, 0], v[:, 1], ".", color="0.6", ms=3)[0]
    l1 = a.plot(rx[:, 0], rx[:, 1], "o", color="k", mfc="white", ms=6)[0]; l2 = a.plot(tx[:, 0], tx[:, 1], "s", color="k", ms=5)[0]
    a.set_aspect("equal"); a.set_xlabel("x, мм"); a.set_ylabel("y, мм"); G.comma(a); G.key(a, [l1, l2, l3]); G.letter(a, ("а)", "б)")[col])
G.save(fig, os.path.join(OUT, "c_layout.png"), letters=True)
