"""Фактор шага отдельно от полосы: одна частота 10 ГГц, кольца ringz_d* (полуволновые вибраторы) и ringz_d*_pt
(точечные излучатели). Плюс расчёт по формулам на частой сетке шагов с теми же определениями.

  python cst/analyze_spacing_1f.py   -> spacing_1f.json
Изображение: MIMO (72 пары), цель на горизонте в азимуте 17°, азимут по всему кругу, УБЛ вне главного лепестка.
Пеленг: кольцо из 12 приёмников, источник через 3° по всему кругу.
"""
import os, re, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 299_792_458.0
F = 10.0; K = 2 * np.pi * F * 1e9 / C0 / 1e3
NRX, NTX = 12, 6
AZ = np.arange(0, 360, 0.25); U = np.column_stack([np.cos(np.radians(AZ)), np.sin(np.radians(AZ)), np.zeros(AZ.size)])
PHS = 17.0


def lobe(x, p):
    """ширина по −3 дБ и наибольший уровень вне главного лепестка (до первых минимумов), круговая ось"""
    d = 20 * np.log10(p / p.max() + 1e-12); n = len(d); i = int(np.argmax(d))
    lo = hi = 0
    while lo > -n // 2 and d[(i + lo - 1) % n] >= -3: lo -= 1
    while hi < n // 2 and d[(i + hi + 1) % n] >= -3: hi += 1
    l = r = 0
    while l > -n // 2 and d[(i + l - 1) % n] < d[(i + l) % n]: l -= 1
    while r < n // 2 and d[(i + r + 1) % n] < d[(i + r) % n]: r += 1
    main = {(i + j) % n for j in range(l, r + 1)}
    out = [d[j] for j in range(n) if j not in main]
    return float((hi - lo) * (x[1] - x[0])), float(max(out)), float(x[i])


def geometry(d_lam):
    R = NRX * d_lam * 30.0 / (2 * np.pi)
    rx = np.array([[R * np.cos(2 * np.pi * i / NRX), R * np.sin(2 * np.pi * i / NRX), 0] for i in range(NRX)])
    tx = np.array([[1.35 * R * np.cos(2 * np.pi * (j + .5) / NTX), 1.35 * R * np.sin(2 * np.pi * (j + .5) / NTX), 0] for j in range(NTX)])
    return rx, tx


def mimo_pattern(RV, x, A=None):
    A = np.exp(1j * K * (RV @ U.T)) if A is None else A
    return np.abs(A.conj().T @ x) / np.sqrt(np.sum(np.abs(A) ** 2, axis=0))


def df_metrics(R, Hsrc):
    """Hsrc(i) -> сигналы 12 приёмников от источника в азимуте i°"""
    fine = np.arange(0, 360, 0.05); A = np.exp(1j * K * (R @ np.column_stack([np.cos(np.radians(fine)), np.sin(np.radians(fine)), np.zeros(fine.size)]).T))
    errs, slls = [], []
    for s_ in np.arange(0, 360, 3.0):
        p = np.abs(A.conj().T @ Hsrc(int(s_))) ** 2
        j = int(np.argmax(p)); jm, jp = (j - 1) % len(p), (j + 1) % len(p)
        dd = 0.5 * (p[jm] - p[jp]) / (p[jm] - 2 * p[j] + p[jp] + 1e-30)
        errs.append(((fine[j] + dd * 0.05) - s_ + 180) % 360 - 180)
        slls.append(lobe(fine, np.sqrt(p))[1])
    e = np.array(errs); big = np.abs(e) > 5
    return dict(max=float(np.abs(e[~big]).max()) if (~big).any() else float("nan"), n_wrong=int(big.sum()), sll_db=float(np.max(slls)),
                bw=float(lobe(fine, np.sqrt(np.abs(A.conj().T @ Hsrc(0)) ** 2))[0]))


out = dict(formula=[], cst={})
# ---- формулы на частой сетке
us = np.array([np.cos(np.radians(PHS)), np.sin(np.radians(PHS)), 0])
for d in np.round(np.arange(0.2, 3.001, 0.05), 3):
    rx, tx = geometry(d); RV = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 3)
    bw, sll, at = lobe(AZ, mimo_pattern(RV, np.exp(1j * K * (RV @ us))))
    bwr, sllr, _ = lobe(AZ, mimo_pattern(rx, np.exp(1j * K * (rx @ us))))
    out["formula"].append(dict(d=float(d), mimo_bw=bw, mimo_sll=sll, rx_bw=bwr, rx_sll=sllr))

# ---- CST
for path in sorted(glob.glob(os.path.join(HERE, "ringz_d*", "ringz_result.json"))):
    folder = os.path.dirname(path); tag = os.path.basename(folder)
    if tag.endswith("_ptwb"):
        continue
    kind = "point" if tag.endswith("_pt") else "dipole"
    r0 = json.load(open(path, encoding="utf-8"))
    lay = [p.split() for p in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")][1:]
    pos = {p[0]: np.array(list(map(float, p[1:4]))) for p in lay}
    rxn = [n for n in pos if n.startswith("rx")]; txn = [n for n in pos if n.startswith("tx")]
    P = {}
    for n in pos:
        dd = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, f"ff_{n}_10.0GHz.txt"), encoding="utf-8")
                       if re.match(r"^\s*[-+]?\d", l)])
        hz = dd[:, 0] == 90; o = np.argsort(dd[hz, 1])
        P[n] = (10 ** (dd[hz, 3] / 20) * np.exp(1j * np.radians(dd[hz, 4])))[o]          # горизонт, азимут 0..359
    pairs = [(t, r) for t in txn for r in rxn]
    RV = np.array([pos[t] + pos[r] for t, r in pairs])
    H = np.array([P[t] * P[r] for t, r in pairs])                                          # (72, 360)
    x = H[:, int(PHS)]
    Hf = np.repeat(H, 4, axis=1); fr = np.tile(np.arange(4) / 4, 360)                   # встроенные ДН на сетке 0,25°
    Hint = (1 - fr) * Hf + fr * np.roll(Hf, -4, axis=1)
    row = dict(d=float(r0["d_lam"]), coupling=r0["coupling"]["rx_rx_db"],
               v10=[v for v in r0["vpc_by_f"] if abs(v["f"] - 10) < 1e-6][0])
    for key, A in (("cst_geo", None), ("cst_full", Hint)):
        bw, sll, at = lobe(AZ, mimo_pattern(RV, x, A))
        row[key] = dict(bw=bw, sll=sll, az_err=at - PHS)
    bw, sll, at = lobe(AZ, mimo_pattern(RV, np.exp(1j * K * (RV @ us))))
    row["formula"] = dict(bw=bw, sll=sll, az_err=at - PHS)
    R = np.array([pos[r] for r in rxn])
    row["df_cst"] = df_metrics(R, lambda i: np.array([P[r][i] for r in rxn]))
    row["df_formula"] = df_metrics(R, lambda i: np.exp(1j * K * (R @ np.array([np.cos(np.radians(i)), np.sin(np.radians(i)), 0]))))
    out["cst"].setdefault(kind, []).append(row)
    print(f"{tag}: MIMO УБЛ формула {row['formula']['sll']:.1f}, CST {row['cst_geo']['sll']:.1f}, CST полн. {row['cst_full']['sll']:.1f} дБ; "
          f"пеленг макс {row['df_cst']['max']:.3f}°, ложных {row['df_cst']['n_wrong']}, УБЛ спектра {row['df_cst']['sll_db']:.1f} (форм. {row['df_formula']['sll_db']:.1f})")
for k in out["cst"]:
    out["cst"][k].sort(key=lambda r: r["d"])
json.dump(out, open(os.path.join(HERE, "spacing_1f.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
for r in out["formula"][::5]:
    print(f"формулы d {r['d']:.2f}: MIMO луч {r['mimo_bw']:.1f}°, УБЛ {r['mimo_sll']:.1f}; Rx луч {r['rx_bw']:.1f}°, УБЛ {r['rx_sll']:.1f}")
