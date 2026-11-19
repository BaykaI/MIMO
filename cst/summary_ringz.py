"""Сводка серии колец (build_ringz.py + analyze_ringz.py) по двум вопросам: полоса и шаг.
  python cst/summary_ringz.py   -> fig2/ringz_spacing.png, fig2/ringz_band.png, fig2/ringz_image.png, ringz_summary.json
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig2")
runs = []
for p in glob.glob(os.path.join(HERE, "ringz_d*", "ringz_result.json")):
    runs.append(json.load(open(p, encoding="utf-8")))
runs.sort(key=lambda r: r["d_lam"])
form = json.load(open(os.path.join(HERE, "formula_band_spacing.json"), encoding="utf-8"))
D = np.array([r["d_lam"] for r in runs])
print("шаги:", D)


def at(r, f):
    return [v for v in r["vpc_by_f"] if abs(v["f"] - f) < 1e-6][0]


summ = []
for r in runs:
    v10 = at(r, 10.0)
    im = r["image"]
    summ.append(dict(
        d=r["d_lam"], d_mm=r["d_mm"], R=r["R_rx"],
        rxrx=r["coupling"]["rx_rx_db"], rxtx=r["coupling"]["rx_tx_db"], s11w=r["coupling"]["s11_db_worst_band"],
        perr10=v10["phase_err_rms_deg"], perr_max_band=max(v["phase_err_rms_deg"] for v in r["vpc_by_f"]),
        perr_min_band=min(v["phase_err_rms_deg"] for v in r["vpc_by_f"]),
        loc10=v10.get("local_dev_mean", float("nan")), loc10max=v10.get("local_dev_max", float("nan")),
        glob10=v10.get("global_dev_mean", float("nan")), glob10rms=v10.get("global_rms_deg", float("nan")),
        delay_dr=r["delay"]["dr_mean"], delay_resid=r["delay"]["resid_mm"], tau_spread=r["delay"]["tau0_spread_mm"], tau_std=r["delay"]["tau0_std_mm"],
        df1=r["df"]["1 частота 10 ГГц"], dfb=r["df"]["полоса 8,5-11,5 ГГц"],
        img={k: v for k, v in im.items()}, img_band=r["image_vs_band"]))
json.dump(summ, open(os.path.join(HERE, "ringz_summary.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ---------------- рисунок 1: зависимость от шага
fs = form["spacing"]; fd = np.array([x["d"] for x in fs])
fig, ax = plt.subplots(2, 3, figsize=(15, 8))
a = ax[0, 0]
a.plot(D, [s["rxrx"] for s in summ], "ko-", label="соседние приёмники"); a.plot(D, [s["rxtx"] for s in summ], "rs-", label="приёмник и ближайший передатчик")
a.set_title("взаимная связь (|S| на 10 ГГц)", fontsize=10); a.set_ylabel("дБ"); a.legend(fontsize=8)
a = ax[0, 1]
a.plot(D, [s["perr10"] for s in summ], "ko-", label="10 ГГц")
a.fill_between(D, [s["perr_min_band"] for s in summ], [s["perr_max_band"] for s in summ], color="0.8", label="разброс по полосе 8,5-11,5 ГГц")
a.set_title("фазовая ошибка геометрической модели r_tx + r_rx (скз по углам)", fontsize=10); a.set_ylabel("градусы"); a.legend(fontsize=8)
a = ax[0, 2]
a.plot(D, [s["loc10"] for s in summ], "ko-", label="местный центр, поперёк луча (среднее)")
a.plot(D, [s["delay_dr"] for s in summ], "b^-", label="центр задержки (широкая полоса)")
a.axhline(1.5, color="r", ls=":", lw=0.8); a.text(D.max() * 0.7, 1.7, "λ/20", color="r", fontsize=8)
a.set_title("смещение фазового центра пары от r_tx + r_rx", fontsize=10); a.set_ylabel("мм"); a.legend(fontsize=8); a.set_yscale("log")
a = ax[1, 0]
a.plot(D, [s["df1"]["max"] for s in summ], "ko-", label="CST, модель геометрии, 1 частота")
a.plot(D, [s["dfb"]["max"] for s in summ], "k^--", label="CST, модель геометрии, полоса")
a.set_title("ошибка пеленга кольцом Rx (максимум по кругу)", fontsize=10); a.set_ylabel("градусы"); a.legend(fontsize=8)
a = ax[1, 1]
key_f, key_g, key_c, key_m = "формулы (идеал)", "CST, обработка по геометрии", "CST, геометрия + калибровка каналов", "CST, полная модель CST"
for key, st, lab in ((key_f, "k-", "формулы Чапурского"), (key_g, "r--", "CST, обработка по геометрии"), (key_c, "b:", "CST, + калибровка каналов"), (key_m, "g-.", "CST, полная модель")):
    a.plot(D, [s["img"][key]["sll_az_db"] for s in summ], st, marker="o", ms=4, label=lab)
a.plot(fd, [x["mimo_band_sll"] for x in fs], color="0.6", lw=0.8, label="формулы, сплошной перебор шага")
a.set_title("MIMO-изображение, полоса 30 %: УБЛ по азимуту", fontsize=10); a.set_ylabel("дБ"); a.legend(fontsize=7)
a = ax[1, 2]
for key, st, lab in ((key_f, "k-", "формулы Чапурского"), (key_g, "r--", "CST, обработка по геометрии"), (key_m, "g-.", "CST, полная модель")):
    a.plot(D, [s["img"][key]["az_res"] for s in summ], st, marker="o", ms=4, label=lab)
a.plot(fd, [x["mimo_band_bw"] for x in fs], color="0.6", lw=0.8, label="формулы, сплошной перебор шага")
a.set_yscale("log"); a.set_title("MIMO-изображение: ширина по азимуту (−3 дБ)", fontsize=10); a.set_ylabel("градусы"); a.legend(fontsize=7)
for a in ax.ravel():
    a.grid(alpha=0.3, which="both"); a.set_xlabel("шаг приёмников по дуге d/λ")
fig.suptitle("Кольцо 12 Rx + 6 Tx, вертикальные диполи: что меняется с шагом между излучателями (CST)", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ringz_spacing.png"), dpi=160); plt.close(fig)

# ---------------- рисунок 2: зависимость от частоты и полосы
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
cm = plt.cm.viridis(np.linspace(0, 0.9, len(runs)))
for r, c in zip(runs, cm):
    ax[0].plot([v["f"] for v in r["vpc_by_f"]], [v["phase_err_rms_deg"] for v in r["vpc_by_f"]], "o-", color=c, ms=3, label=f"d = {r['d_lam']:.2g}λ")
ax[0].set_xlabel("частота, ГГц"); ax[0].set_ylabel("скз фазовой ошибки, °"); ax[0].set_title("фазовая ошибка геометрии на каждой частоте", fontsize=10)
ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
ax[1].plot(D, [s["delay_dr"] for s in summ], "b^-", label="смещение центра задержки")
ax[1].plot(D, [s["delay_resid"] for s in summ], "kv-", label="остаток задержки по направлениям")
ax[1].plot(D, [s["tau_std"] for s in summ], "rs-", label="разброс собственной задержки пар (СКО)")
ax[1].set_xlabel("шаг d/λ"); ax[1].set_ylabel("путь c·τ, мм"); ax[1].set_title("широкополосный сигнал: задержки пар", fontsize=10)
ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3); ax[1].set_yscale("log")
for r, c in zip(runs, cm):
    rows = [b for b in r["image_vs_band"] if b["nf"] > 1]
    ax[2].plot([b["rel"] * 100 for b in rows], [b["cst_geo"]["range_res_mm"] for b in rows], "o-", color=c, ms=3, label=f"CST, d = {r['d_lam']:.2g}λ")
rows = [b for b in runs[0]["image_vs_band"] if b["nf"] > 1]
ax[2].plot([b["rel"] * 100 for b in rows], [b["formula"]["range_res_mm"] for b in rows], "k--", lw=1.5, label="формулы")
ax[2].set_xlabel("относительная полоса, %"); ax[2].set_ylabel("разрешение по дальности, мм"); ax[2].set_yscale("log")
ax[2].set_title("разрешение по дальности от ширины полосы", fontsize=10); ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ringz_band.png"), dpi=160); plt.close(fig)

# ---------------- рисунок 3: изображения дальность-азимут для двух шагов
pick = [r for r in runs if abs(r["d_lam"] - 0.5) < 0.05] + [r for r in runs if abs(r["d_lam"] - 1.5) < 0.05]
if pick:
    fig, ax = plt.subplots(len(pick), 2, figsize=(11, 4 * len(pick)), squeeze=False)
    for row, r in enumerate(pick):
        ax_ = r["image_axes"]; dr = np.linspace(*ax_["dr"]); az = np.linspace(*ax_["az"])
        for col, (fn, lab) in enumerate((("img_formula.npy", "формулы Чапурского"), ("img_cst_geo.npy", "CST, обработка по геометрии"))):
            p = os.path.join(HERE, r["folder"], fn)
            if not os.path.exists(p):
                continue
            I = np.load(p); im = ax[row, col].imshow(20 * np.log10(I + 1e-9), extent=[az[0], az[-1], dr[0], dr[-1]], origin="lower", aspect="auto", vmin=-30, vmax=0, cmap="viridis")
            ax[row, col].set_xlim(-40, 75); ax[row, col].set_ylim(-150, 150)
            ax[row, col].set_title(f"d = {r['d_lam']:.2g}λ, {lab}", fontsize=10); ax[row, col].set_xlabel("азимут, °"); ax[row, col].set_ylabel("дальность, мм")
    fig.colorbar(im, ax=ax, shrink=0.6, label="дБ")
    fig.savefig(os.path.join(OUT, "ringz_image.png"), dpi=150, bbox_inches="tight"); plt.close(fig)

for s in summ:
    print(f"d={s['d']:.2f}λ: связь {s['rxrx']:.1f}/{s['rxtx']:.1f} дБ, фаз. ошибка {s['perr10']:.1f}° ({s['perr_min_band']:.1f}-{s['perr_max_band']:.1f}), "
          f"местный центр {s['loc10']:.2f} мм, центр задержки {s['delay_dr']:.2f} мм, разброс задержки {s['tau_std']:.2f} мм, "
          f"пеленг {s['df1']['max']:.3f}°/{s['dfb']['max']:.3f}°, УБЛ изобр. формула {s['img'][key_f]['sll_az_db']:.1f} / CST {s['img'][key_g]['sll_az_db']:.1f} / полн. {s['img'][key_m]['sll_az_db']:.1f}")
