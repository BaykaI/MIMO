"""Ожидание по формулам Чапурского (точечные излучатели) для двух вопросов.

1) Широкополосный сигнал. Отклик пары на частоте f с направления u: exp(j·2πf/c·(r_tx + r_rx)·u).
   Фазовый центр пары (по наклону фазы по углу) и центр задержки (по наклону фазы по частоте) у точечных излучателей
   совпадают и не зависят от частоты. Широкая полоса меняет не центры, а изображение: появляется разрешение по
   дальности c/(2B), а разные частоты дают дифракционные лепестки в разных местах, и они размываются.
   ОФН по ф. (11) статьи [7]: Ψ(Δ, φ) = |Σ_m Σ_nk exp{j 2π f_m [2Δ/c − (r_v·(u(φ) − u0))/c]}| / (Nt Nr M).
2) Шаг. Кольцо 12 Rx на R = 12 d/(2π), 6 Tx на 1,35 R. Перебор d/λ: ширина луча и наибольший ложный
   максимум по азимуту для пассивного кольца Rx и для MIMO-эквивалента, на одной частоте и в полосе 8,5-11,5 ГГц.

  python cst/formula_band_spacing.py   -> fig2/formula_band.png, fig2/formula_spacing.png, formula_band_spacing.json
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C0 = 299_792_458.0
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig2"); os.makedirs(OUT, exist_ok=True)
F0 = 10e9; LAM = C0 / F0 * 1e3
BAND = np.round(np.arange(8.5, 11.51, 0.25), 2) * 1e9
NRX, NTX = 12, 6


def ring(d_lam):
    R = NRX * d_lam * LAM / (2 * np.pi)
    a = 2 * np.pi * np.arange(NRX) / NRX; b = 2 * np.pi * (np.arange(NTX) + 0.5) / NTX
    rx = np.column_stack([R * np.cos(a), R * np.sin(a)]); tx = np.column_stack([1.35 * R * np.cos(b), 1.35 * R * np.sin(b)])
    return rx, tx


def af_az(pts, freqs, az, az0=0.0):
    """|Σ_f Σ_p exp(j k_f p·(u − u0))| на горизонте, нормированный"""
    u = np.column_stack([np.cos(np.radians(az)), np.sin(np.radians(az))]); u0 = np.array([np.cos(np.radians(az0)), np.sin(np.radians(az0))])
    s = np.zeros(len(az), complex)
    for f in freqs:
        k = 2 * np.pi * f / C0 / 1e3
        s += np.exp(1j * k * ((u - u0) @ pts.T)).sum(axis=1)
    return np.abs(s) / (len(pts) * len(freqs))


def metrics(az, p):
    d = 20 * np.log10(p / p.max() + 1e-12); i = int(np.argmax(d)); lo = hi = i
    while d[lo - 1] >= -3: lo -= 1
    while d[(hi + 1) % len(d)] >= -3: hi += 1
    l = i; r = i
    while d[l - 1] < d[l]: l -= 1
    while d[(r + 1) % len(d)] < d[r % len(d)]: r += 1
    mask = np.ones(len(d), bool); mask[np.arange(l, r + 1) % len(d)] = False
    return float((hi - lo) * (az[1] - az[0])), float(d[mask].max())


az = np.arange(-180, 180, 0.05)
res = {"spacing": [], "band": {}}

# ---------------- 2) шаг
ds = np.round(np.arange(0.2, 3.01, 0.05), 3)
for dl in ds:
    rx, tx = ring(dl); v = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2)
    row = dict(d=float(dl))
    for lab, pts in (("rx", rx), ("mimo", v)):
        for flab, fr in (("1f", [F0]), ("band", BAND)):
            bw, sl = metrics(az, af_az(pts, fr, az)); row[f"{lab}_{flab}_bw"] = bw; row[f"{lab}_{flab}_sll"] = sl
    res["spacing"].append(row)
S = res["spacing"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for key, st, lab in (("rx_1f", "k-", "кольцо Rx, 1 частота"), ("rx_band", "k--", "кольцо Rx, полоса 8,5-11,5 ГГц"),
                     ("mimo_1f", "r-", "MIMO-эквивалент, 1 частота"), ("mimo_band", "r--", "MIMO-эквивалент, полоса")):
    ax[0].plot(ds, [r[key + "_bw"] for r in S], st, lw=1.2, label=lab); ax[1].plot(ds, [r[key + "_sll"] for r in S], st, lw=1.2, label=lab)
ax[0].set_yscale("log"); ax[0].set_xlabel("шаг приёмников по дуге d/λ"); ax[0].set_ylabel("ширина луча по −3 дБ, °"); ax[0].grid(alpha=0.3, which="both")
ax[0].set_title("ширина луча по азимуту", fontsize=10); ax[0].legend(fontsize=8)
ax[1].axhline(-6, color="0.5", lw=0.7, ls=":"); ax[1].text(2.4, -5.5, "−6 дБ", color="0.4", fontsize=8)
ax[1].set_xlabel("шаг приёмников по дуге d/λ"); ax[1].set_ylabel("наибольший ложный максимум, дБ"); ax[1].grid(alpha=0.3)
ax[1].set_title("неоднозначность пеленга (ложные максимумы)", fontsize=10); ax[1].legend(fontsize=8)
fig.suptitle("Кольцо 12 Rx + 6 Tx по формулам Чапурского: зависимость от шага", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "formula_spacing.png"), dpi=170); plt.close(fig)

# ---------------- 1) полоса: ОФН дальность-азимут для d = 0.5λ и 1.5λ, и кольцо ДЗ1 (R = 50λ) с 1/4/16 частотами
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for i, dl in enumerate((0.5, 1.5)):
    rx, tx = ring(dl); v = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2)
    dr = np.linspace(-150, 150, 301); azs = np.arange(-40, 40.01, 0.25)
    U = np.column_stack([np.cos(np.radians(azs)), np.sin(np.radians(azs))]); u0 = np.array([1.0, 0.0])
    img = np.zeros((len(dr), len(azs)), complex)
    for f in BAND:
        k = 2 * np.pi * f / C0 / 1e3
        geo = np.exp(1j * k * ((U - u0) @ v.T)).sum(axis=1)               # азимутальная часть
        img += np.exp(1j * 2 * k * dr)[:, None] * geo[None, :]               # 2Δ — путь туда и обратно
    I = 20 * np.log10(np.abs(img) / np.abs(img).max() + 1e-9)
    im = ax[i].imshow(I, extent=[azs[0], azs[-1], dr[0], dr[-1]], origin="lower", aspect="auto", vmin=-30, vmax=0, cmap="viridis")
    ax[i].set_xlabel("азимут, °"); ax[i].set_ylabel("дальность относительно цели, мм")
    prof = np.abs(img[:, np.argmin(np.abs(azs))]); prof /= prof.max()
    j = int(np.argmax(prof)); lo = hi = j
    while prof[lo - 1] >= 0.707: lo -= 1
    while prof[hi + 1] >= 0.707: hi += 1
    rres = float(dr[hi] - dr[lo]); res["band"][f"d{dl}"] = dict(range_res_mm=rres, range_res_theory_mm=C0 / (2 * (BAND[-1] - BAND[0])) * 1e3)
    ax[i].set_title(f"d = {dl}λ, полоса 8,5-11,5 ГГц" + chr(10) + f"разрешение по дальности {rres:.0f} мм", fontsize=9)
fig.colorbar(im, ax=ax[:2], shrink=0.8, label="дБ")
# кольцо ДЗ1: 20 Rx на 50λ, 10 Tx на 67,5λ, азимут, 1/4/16 частот
a = 2 * np.pi * np.arange(20) / 20; b = 2 * np.pi * np.arange(10) / 10
rx = np.column_stack([50 * LAM * np.cos(a), 50 * LAM * np.sin(a)]); tx = np.column_stack([67.5 * LAM * np.cos(b), 67.5 * LAM * np.sin(b)])
v = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2); az2 = np.arange(-180, 180, 0.01)
for M, st in ((1, "k-"), (4, "b-"), (16, "r-")):
    fr = np.linspace(10e9, 16e9, M) if M > 1 else [10e9]
    p = af_az(v, fr, az2); bw, sl = metrics(az2, p); res["band"][f"ring50_M{M}"] = dict(bw=bw, sll=sl)
    ax[2].plot(az2, 20 * np.log10(p / p.max() + 1e-9), st, lw=0.8, label=f"{M} част. (10-16 ГГц): УБЛ {sl:.1f} дБ")
ax[2].set_xlim(-10, 10); ax[2].set_ylim(-30, 1); ax[2].grid(alpha=0.3); ax[2].legend(fontsize=7)
ax[2].set_xlabel("азимут, °"); ax[2].set_ylabel("дБ"); ax[2].set_title("кольцо ДЗ1 (R = 50λ), азимут" + chr(10) + "многочастотный сигнал гасит ложные максимумы", fontsize=9)
fig.savefig(os.path.join(OUT, "formula_band.png"), dpi=170, bbox_inches="tight"); plt.close(fig)
json.dump(res, open(os.path.join(HERE, "formula_band_spacing.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
for r in S:
    if abs(r["d"] * 100 % 25) < 1e-6:
        print(f"d={r['d']:.2f}λ  Rx: {r['rx_1f_bw']:.2f}° {r['rx_1f_sll']:6.1f} дБ | полоса {r['rx_band_sll']:6.1f} дБ   MIMO: {r['mimo_1f_bw']:.2f}° {r['mimo_1f_sll']:6.1f} дБ | полоса {r['mimo_band_sll']:6.1f} дБ")
print(json.dumps(res["band"], ensure_ascii=False))
