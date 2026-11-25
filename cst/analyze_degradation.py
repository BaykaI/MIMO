"""Деградация: отдельно от ширины полосы и отдельно от шага.

  python cst/analyze_degradation.py  -> cst/degr.json, cst/degr.npz

ПОЛОСА (шаг постоянный, точечные излучатели = идеальные широкополосные; CST ringz_d*_ptwb 5-15 ГГц + формулы
для 6 шагов). Полоса 0-100 % вокруг 10 ГГц через 10 %. Две обработки изображения:
  wb  многочастотная (Чапурский): на каждой частоте фаза пары считается на своей частоте;
  nb  узкополосная: фазы пар берутся на центральной частоте 10 ГГц для всех составляющих (как фазовращатели).
Метрики через пик изображения: потеря пика относительно wb, ошибка азимута и дальности, ширина по азимуту и
дальности, УБЛ по азимуту (весь круг) и по дальности.

ШАГ (одна частота 10 ГГц, CST ringz_d* полуволновые и ringz_d*_pt точечные, 6 шагов): сечения изображения по
азимуту, спектр пеленгатора, остаточная фаза пары по азимуту (искажение фазового центра).
"""
import os, re, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 299_792_458.0
kf = lambda f: 2 * np.pi * f * 1e9 / C0 / 1e3
NRX, NTX = 12, 6
F0 = 10.0
PHS = 17.0
AZ = np.arange(0, 360, 0.25); U = np.column_stack([np.cos(np.radians(AZ)), np.sin(np.radians(AZ)), np.zeros(AZ.size)])
us = np.array([np.cos(np.radians(PHS)), np.sin(np.radians(PHS)), 0.0])
DR = np.arange(-250, 250.1, 1.0)          # период по дальности при шаге частот 0,5 ГГц 300 мм: повторы пика вне окна
RELS = [round(0.1 * i, 1) for i in range(11)]
FGRID = [round(5 + 0.5 * i, 2) for i in range(21)]
SPACINGS = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0]


def geometry(d):
    R = NRX * d * 30.0 / (2 * np.pi)
    rx = np.array([[R * np.cos(2 * np.pi * i / NRX), R * np.sin(2 * np.pi * i / NRX), 0] for i in range(NRX)])
    tx = np.array([[1.35 * R * np.cos(2 * np.pi * (j + .5) / NTX), 1.35 * R * np.sin(2 * np.pi * (j + .5) / NTX), 0] for j in range(NTX)])
    return rx, tx


def lobe(x, p, circ):
    d = 20 * np.log10(p / p.max() + 1e-12); n = len(d); i = int(np.argmax(d))
    idx = (lambda j: j % n) if circ else (lambda j: min(max(j, 0), n - 1))
    lo = hi = 0
    while lo > -(n // 2) and d[idx(i + lo - 1)] >= -3 and (circ or i + lo > 0): lo -= 1
    while hi < n // 2 and d[idx(i + hi + 1)] >= -3 and (circ or i + hi < n - 1): hi += 1
    l = r = 0
    while l > -(n // 2) and d[idx(i + l - 1)] < d[idx(i + l)] and (circ or i + l > 0): l -= 1
    while r < n // 2 and d[idx(i + r + 1)] < d[idx(i + r)] and (circ or i + r < n - 1): r += 1
    main = {idx(i + j) for j in range(l, r + 1)}
    out = [d[j] for j in range(n) if j not in main]
    return float((hi - lo) * (x[1] - x[0])), float(max(out)) if out else float("nan"), float(x[i])


def image(RV, xs, freqs, mode):
    """xs[f]: сигналы 72 пар от цели; возвращает комплексное изображение (DR, AZ)"""
    img = np.zeros((len(DR), len(AZ)), complex)
    A0 = np.exp(1j * kf(F0) * (RV @ U.T))
    for f in freqs:
        A = np.exp(1j * kf(f) * (RV @ U.T)) if mode == "wb" else A0
        b = (A.conj().T @ xs[f]) / len(RV)
        img += np.exp(1j * 2 * kf(f) * DR)[:, None] * b[None, :]
    return img / len(freqs)


def metrics(I, I_ref_peak):
    a = np.abs(I); i, j = np.unravel_index(np.argmax(a), a.shape)
    bw_az, sll_az, az_at = lobe(AZ, a[i, :], True)
    if a.shape[0] > 1:
        bw_r, sll_r, r_at = lobe(DR, a[:, j], False)
    out = dict(loss_db=float(20 * np.log10(a.max() / I_ref_peak)), az_err=float(((az_at - PHS + 180) % 360) - 180),
               az_bw=bw_az, sll_az=sll_az, r_err=float(r_at), r_res=bw_r, sll_r=sll_r)
    return out, i, j


sub = lambda freqs, rel: [f for f in freqs if abs(f - F0) <= rel * F0 / 2 + 1e-9]
OUT = dict(band=dict(formula={}, cst={}), spacing={})
NPZ = {}

# ======================= ПОЛОСА
def band_study(RV, xs, freqs, tag, keep_cuts=False):
    rows = []
    for rel in RELS:
        fr = sub(freqs, rel)
        Iw = image(RV, xs, fr, "wb"); pk = np.abs(Iw).max()
        mw, iw, jw = metrics(Iw, pk)
        In = image(RV, xs, fr, "nb"); mn, i_n, j_n = metrics(In, pk)
        if len(fr) == 1:
            for m_ in (mw, mn):
                m_.update(r_res=float("nan"), sll_r=float("nan"), r_err=float("nan"))
        rows.append(dict(rel=rel, nf=len(fr), wb=mw, nb=mn))
        if keep_cuts and rel in (0.0, 0.2, 0.4, 0.6, 1.0):
            r10 = int(rel * 100)
            NPZ[f"{tag}_r{r10}_wb_az"] = np.abs(Iw[iw, :]); NPZ[f"{tag}_r{r10}_nb_az"] = np.abs(In[i_n, :])
            NPZ[f"{tag}_r{r10}_wb_rng"] = np.abs(Iw[:, jw]); NPZ[f"{tag}_r{r10}_nb_rng"] = np.abs(In[:, j_n])
    return rows


for d in SPACINGS:
    rx, tx = geometry(d); RV = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 3)
    xs = {f: np.exp(1j * kf(f) * (RV @ us)) for f in FGRID}
    OUT["band"]["formula"][f"{d:g}"] = dict(d=d, dL_mm=float(np.ptp(RV @ us)), rows=band_study(RV, xs, FGRID, f"f{d:g}", keep_cuts=True))
    print("полоса, формулы, шаг", d, "ok", flush=True)


def read_ff(folder, name, f):
    dd = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, f"ff_{name}_{f:.1f}GHz.txt"), encoding="utf-8")
                   if re.match(r"^\s*[-+]?\d", l)])
    hz = dd[:, 0] == 90; o = np.argsort(dd[hz, 1])
    return (10 ** (dd[hz, 3] / 20) * np.exp(1j * np.radians(dd[hz, 4])))[o]


def layout(folder):
    lay = [p.split() for p in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")][1:]
    pos = {p[0]: np.array(list(map(float, p[1:4]))) for p in lay}
    return pos, [n for n in pos if n.startswith("rx")], [n for n in pos if n.startswith("tx")]


for folder in sorted(glob.glob(os.path.join(HERE, "ringz_d*_ptwb"))):
    pos, rxn, txn = layout(folder); pairs = [(t, r) for t in txn for r in rxn]
    RV = np.array([pos[t] + pos[r] for t, r in pairs])
    wb = json.load(open(os.path.join(folder, "wb_result.json"), encoding="utf-8")); tau = wb["delay"]["tau0_mean_mm"]
    freqs = [float(x) for x in open(os.path.join(folder, "freqs.txt")).read().split()]
    xs = {}
    for f in freqs:
        P = {n: read_ff(folder, n, f) for n in pos}
        xs[f] = np.array([P[t][int(PHS)] * P[r][int(PHS)] for t, r in pairs]) * np.exp(1j * kf(f) * tau)
    d = round(wb["d_lam"], 2)
    OUT["band"]["cst"][f"{d:g}"] = dict(d=d, dL_mm=float(np.ptp(RV @ us)), rows=band_study(RV, xs, freqs, f"c{d:g}", keep_cuts=True))
    print("полоса, CST, шаг", d, "ok", flush=True)

# ======================= ШАГ (10 ГГц)
fine = np.arange(0, 360, 0.1)
Uf = np.column_stack([np.cos(np.radians(fine)), np.sin(np.radians(fine)), np.zeros(fine.size)])
for d in SPACINGS:
    row = dict(d=d)
    for kind, suf in (("dipole", ""), ("point", "_pt")):
        folder = os.path.join(HERE, f"ringz_d{d:g}{suf}")
        pos, rxn, txn = layout(folder); pairs = [(t, r) for t in txn for r in rxn]
        RV = np.array([pos[t] + pos[r] for t, r in pairs]); R = np.array([pos[r] for r in rxn])
        P = {n: read_ff(folder, n, F0) for n in pos}
        x = np.array([P[t][int(PHS)] * P[r][int(PHS)] for t, r in pairs])
        A = np.exp(1j * kf(F0) * (RV @ U.T))
        pat = np.abs(A.conj().T @ x); pat_f = np.abs(A.conj().T @ np.exp(1j * kf(F0) * (RV @ us)))
        NPZ[f"s{d:g}_{kind}_img"] = pat / pat.max(); NPZ[f"s{d:g}_formula_img"] = pat_f / pat_f.max()
        bw, sll, at = lobe(AZ, pat, True)
        # спектр пеленгатора (кольцо Rx), источник в 17°
        Ar = np.exp(1j * kf(F0) * (R @ Uf.T)); xr = np.array([P[r][int(PHS)] for r in rxn])
        sp = np.abs(Ar.conj().T @ xr); spf = np.abs(Ar.conj().T @ np.exp(1j * kf(F0) * (R @ us)))
        NPZ[f"s{d:g}_{kind}_df"] = sp / sp.max(); NPZ[f"s{d:g}_formula_df"] = spf / spf.max()
        dbw, dsll, dat = lobe(fine, sp, True)
        # остаточная фаза пары tx1·rx1 по азимуту (за вычетом формулы и средней), град.
        rv = pos["tx1"] + pos["rx1"]
        z = P["tx1"] * P["rx1"] * np.exp(-1j * kf(F0) * (U[::4] @ rv))
        w = np.abs(z) ** 2; c = np.angle(np.sum(w * z / np.abs(z)))
        NPZ[f"s{d:g}_{kind}_resid"] = np.degrees(np.angle(z * np.exp(-1j * c)))
        NPZ[f"s{d:g}_{kind}_amp"] = 20 * np.log10(np.abs(z) / np.abs(z).max())
        row[kind] = dict(img_bw=bw, img_sll=sll, img_az_err=float(((at - PHS + 180) % 360) - 180),
                         df_bw=dbw, df_sll=dsll, df_err=float(((dat - PHS + 180) % 360) - 180),
                         resid_pp=float(np.ptp(NPZ[f"s{d:g}_{kind}_resid"])))
    bwf, sllf, _ = lobe(AZ, pat_f, True); dbwf, dsllf, _ = lobe(fine, spf, True)
    row["formula"] = dict(img_bw=bwf, img_sll=sllf, df_bw=dbwf, df_sll=dsllf)
    OUT["spacing"][f"{d:g}"] = row
    print(f"шаг {d}: изобр. УБЛ форм {sllf:.1f} / вибр {row['dipole']['img_sll']:.1f} / точ {row['point']['img_sll']:.1f}; "
          f"пеленг. ложн. форм {dsllf:.1f} / вибр {row['dipole']['df_sll']:.1f}; остаток фазы p-p вибр {row['dipole']['resid_pp']:.0f}°", flush=True)

NPZ["az"] = AZ; NPZ["dr"] = DR; NPZ["fine"] = fine; NPZ["az_resid"] = AZ[::4]
np.savez(os.path.join(HERE, "degr.npz"), **NPZ)
json.dump(OUT, open(os.path.join(HERE, "degr.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
for src in ("formula", "cst"):
    for k, v in OUT["band"][src].items():
        print(f"полоса {src} шаг {k} (ΔL {v['dL_mm']:.0f} мм): " + "; ".join(
            f"{r['rel']*100:.0f}%: нб потеря {r['nb']['loss_db']:.1f} дБ, УБЛ {r['nb']['sll_az']:.1f}/{r['wb']['sll_az']:.1f}, разр {r['nb']['r_res']:.0f}/{r['wb']['r_res']:.0f}"
            for r in v["rows"][::2]))
