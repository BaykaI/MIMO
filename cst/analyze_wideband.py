"""Фактор полосы отдельно от шага: кольцо 12 Rx + 6 Tx с точечными (широкополосными) излучателями,
полоса 5-15 ГГц (шаг 0,5 ГГц), шаг кольца фиксирован (build_ringz.py ... _ptwb 5:15:0.5).

  python cst/analyze_wideband.py cst/ringz_d0.5_ptwb
Выход: <папка>/wb_result.json, <папка>/wb_cuts.npz.

 1. Фаза пары на каждой частоте: фазовая ошибка формулы r_tx + r_rx и единый фазовый центр (5, 10, 15 ГГц).
 2. Фаза пары по частоте в разных направлениях: у каждой составляющей своя фаза, наклон по частоте = время прихода.
 3. Центр задержки по всей полосе и по подполосам.
 4. Изображение (ОФН, MIMO, 72 пары) в зависимости от полосы 0-100 % вокруг 10 ГГц, азимут по всему кругу.
 5. Пеленг кольцом Rx в зависимости от полосы.
"""
import sys, os, re, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from phase_center_fit import fit_phase_center, C0

folder = sys.argv[1]
lay = [p.split() for p in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")][1:]
pos = {p[0]: np.array(list(map(float, p[1:4]))) for p in lay}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]
FREQS = [float(x) for x in open(os.path.join(folder, "freqs.txt")).read().split()]
F0 = 10.0
kf = lambda f: 2 * np.pi * f * 1e9 / C0 / 1e3            # рад/мм


def read(name, f):
    d = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, f"ff_{name}_{f:.1f}GHz.txt"), encoding="utf-8")
                  if re.match(r"^\s*[-+]?\d", l)])
    return d[:, 0], d[:, 1], 10 ** (d[:, 3] / 20) * np.exp(1j * np.radians(d[:, 4]))


P = {}
for f in FREQS:
    for n in pos:
        th, ph, P[(n, f)] = read(n, f)
thr, phr = np.radians(th), np.radians(ph)
u = np.column_stack([np.sin(thr) * np.cos(phr), np.sin(thr) * np.sin(phr), np.cos(thr)])
R_rx = float(np.linalg.norm(pos["rx1"][:2])); d_mm = 2 * np.pi * R_rx / len(rx)
res = dict(folder=os.path.basename(folder), freqs=FREQS, R_rx=R_rx, d_mm=d_mm, d_lam=d_mm / (C0 / 10e9 * 1e3))
pairs = [(t, r) for t in tx for r in rx]
RV = np.array([pos[t] + pos[r] for t, r in pairs])
band = np.abs(th - 90) <= 30

# ---- 1. фаза пары на каждой частоте
res["by_f"] = []
for f in FREQS:
    k = kf(f); perr = []; gdev = []; grms = []
    for (t, r), rv in zip(pairs, RV):
        E = P[(t, f)] * P[(r, f)]
        z = (E * np.exp(-1j * k * (u @ rv)))[band]; w = np.abs(z) ** 2
        c = np.angle(np.sum(w * np.exp(1j * np.angle(z)))); e = np.angle(z * np.exp(-1j * c))
        perr.append(np.degrees(np.sqrt(np.average(e ** 2, weights=w))))
        if abs(f - 5) < 1e-6 or abs(f - 10) < 1e-6 or abs(f - 15) < 1e-6:
            xyz, rms = fit_phase_center(th[band], ph[band], E[band], k, 180.0, (90.0, 0.0), init=rv)
            gdev.append(np.linalg.norm(xyz - rv)); grms.append(rms)
    row = dict(f=f, d_lam_f=d_mm / (C0 / (f * 1e9) * 1e3), phase_err_rms_deg=float(np.mean(perr)), phase_err_max_deg=float(np.max(perr)))
    if gdev:
        row.update(global_dev_mean=float(np.mean(gdev)), global_dev_max=float(np.max(gdev)), global_rms_deg=float(np.mean(grms)))
    res["by_f"].append(row)

# ---- 3. центр задержки (по всей полосе и по подполосам вокруг 10 ГГц)
hor = band; ut = u[hor]


def delay_fit(freqs):
    Fh = np.array(freqs) * 1e9; out = []
    for (t, r), rv in zip(pairs, RV):
        H = np.array([(P[(t, f)] * P[(r, f)])[hor] * np.exp(-1j * kf(f) * (ut @ rv)) for f in freqs])
        ps = np.unwrap(np.angle(H), axis=0)
        slope = np.linalg.lstsq(np.column_stack([np.ones(len(Fh)), Fh]), ps, rcond=None)[0][1]
        L = -slope / (2 * np.pi) * C0 * 1e3
        w = np.mean(np.abs(H) ** 2, axis=0)
        M = np.column_stack([np.ones(len(ut)), -ut[:, 0], -ut[:, 1], -ut[:, 2]])
        sol = np.linalg.lstsq(M * np.sqrt(w)[:, None], L * np.sqrt(w), rcond=None)[0]
        out.append((sol[0], np.linalg.norm(sol[1:3]), np.sqrt(np.average((L - M @ sol) ** 2, weights=w))))
    o = np.array(out)
    return dict(tau0_mean_mm=float(o[:, 0].mean()), tau0_spread_mm=float(np.ptp(o[:, 0])), dr_mean=float(o[:, 1].mean()),
                dr_max=float(o[:, 1].max()), resid_mm=float(o[:, 2].mean()))


RELS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
sub = lambda rel: [f for f in FREQS if abs(f - F0) <= rel * F0 / 2 + 1e-9]
res["delay"] = delay_fit(FREQS)
res["delay_vs_band"] = [dict(rel=rel, nf=len(sub(rel)), **delay_fit(sub(rel))) for rel in RELS if len(sub(rel)) >= 3]
TAU = res["delay"]["tau0_mean_mm"]                          # общая задержка (опорная добавка CST + своя у излучателя)

# ---- 2. фаза пары tx1·rx1 по частоте в направлениях через 30°: CST и формула, общая задержка вычтена
rv = pos["tx1"] + pos["rx1"]; DIRS = np.arange(0, 360, 30.0)
ph_cst, ph_geo = [], []
for dct in DIRS:
    ii = np.where((th == 90) & (ph == dct))[0][0]
    z = np.array([P[("tx1", f)][ii] * P[("rx1", f)][ii] * np.exp(1j * kf(f) * TAU) for f in FREQS])
    zg = np.array([np.exp(1j * kf(f) * (u[ii] @ rv)) for f in FREQS])
    a = np.unwrap(np.angle(z)); g = np.unwrap(np.angle(zg))
    a += 2 * np.pi * np.round((g[0] - a[0]) / (2 * np.pi)) + 0.0          # одна ветвь с формулой на нижней частоте
    off = np.mean(a - g); off = np.angle(np.exp(1j * off))                   # постоянная фаза излучателей
    ph_cst.append(np.degrees(a - off)); ph_geo.append(np.degrees(g))
CUTS = dict(freqs=np.array(FREQS), pair_dirs=DIRS, pair_phase_cst=np.array(ph_cst), pair_phase_geo=np.array(ph_geo),
            pair_rv=rv)

# ---- 4. изображение в зависимости от полосы
sel = th == 90
azg = ph[sel]; order = np.argsort(azg)
Hcst = {f: np.array([(P[(t, f)] * P[(r, f)])[sel][order] for t, r in pairs]) for f in FREQS}
phs = 17.0; i_s = 17
dr = np.arange(-250, 250.1, 1.0); az_im = np.arange(0, 360, 0.25)
U = np.column_stack([np.cos(np.radians(az_im)), np.sin(np.radians(az_im)), np.zeros(az_im.size)])
us = np.array([np.cos(np.radians(phs)), np.sin(np.radians(phs)), 0.0])


def image(signal, freqs):
    img = np.zeros((len(dr), len(az_im)), complex)
    for f in freqs:
        k = kf(f); A = np.exp(1j * k * (RV @ U.T))
        x = Hcst[f][:, i_s] * np.exp(1j * k * TAU) if signal == "cst" else np.exp(1j * k * (RV @ us))
        b = (A.conj().T @ x) / np.sqrt(len(RV))
        img += np.exp(1j * 2 * k * dr)[:, None] * b[None, :]
    return np.abs(img) / np.abs(img).max()


def lobe(x, p, circ=False):
    d = 20 * np.log10(p / p.max() + 1e-12); i = int(np.argmax(d)); n = len(d)
    idx = (lambda j: j % n) if circ else (lambda j: min(max(j, 0), n - 1))
    lo = hi = 0
    while lo > -n // 2 and d[idx(i + lo - 1)] >= -3: lo -= 1
    while hi < n // 2 and d[idx(i + hi + 1)] >= -3: hi += 1
    l = r = 0
    while l > -n // 2 and d[idx(i + l - 1)] < d[idx(i + l)] and (circ or i + l > 0): l -= 1
    while r < n // 2 and d[idx(i + r + 1)] < d[idx(i + r)] and (circ or i + r < n - 1): r += 1
    main = {idx(i + j) for j in range(l, r + 1)}
    out = [d[j] for j in range(n) if j not in main]
    step = x[1] - x[0]
    flat = (not circ) and (i + lo <= 0 and i + hi >= n - 1)
    return (float("nan") if flat else float((hi - lo) * step)), (float(max(out)) if out else float("nan")), float(x[i])


res["image_vs_band"] = []
for rel in RELS:
    fr = sub(rel); row = dict(rel=rel, nf=len(fr), f_lo=fr[0], f_hi=fr[-1])
    for key in ("formula", "cst"):
        I = image(key, fr); i, j = np.unravel_index(np.argmax(I), I.shape)
        rres, sll_r, r_at = lobe(dr, I[:, j]) if len(fr) > 1 else (float('nan'), float('nan'), 0.0)
        ares, sll_a, a_at = lobe(az_im, I[i, :], circ=True)
        row[key] = dict(range_res_mm=rres, sll_r_db=sll_r, range_err_mm=r_at, az_res=ares, sll_az_db=sll_a, az_err=a_at - phs)
        CUTS[f"band{int(rel*100)}_{key}_range"] = I[:, j]; CUTS[f"band{int(rel*100)}_{key}_az"] = I[i, :]
    res["image_vs_band"].append(row)
# одна частота на краях полосы: электрический шаг другой
res["image_1f"] = []
for f in [x for x in (5.0, 7.5, 10.0, 12.5, 15.0) if x in Hcst]:
    I = image("cst", [f]); Ig = image("formula", [f])
    row = dict(f=f, d_lam_f=d_mm / (C0 / (f * 1e9) * 1e3))
    for key, J in (("cst", I), ("formula", Ig)):
        ares, sll_a, a_at = lobe(az_im, J[len(dr) // 2, :], circ=True)
        row[key] = dict(az_res=ares, sll_az_db=sll_a, az_err=a_at - phs)
    res["image_1f"].append(row)
CUTS["dr"] = dr; CUTS["az"] = az_im

# ---- 5. пеленг кольцом Rx в зависимости от полосы
fine = np.arange(0, 360, 0.05); uf = np.column_stack([np.cos(np.radians(fine)), np.sin(np.radians(fine)), np.zeros(fine.size)])
Rr = np.array([pos[r] for r in rx])
Hr = {f: np.array([P[(r, f)][sel][order] for r in rx]) for f in FREQS}


def df(freqs):
    errs, slls = [], []
    for s_ in np.arange(0, 360, 3.0):
        p = np.zeros(len(fine)); i = int(s_)
        for f in freqs:
            A = np.exp(1j * kf(f) * (Rr @ uf.T)); p += np.abs(A.conj().T @ Hr[f][:, i]) ** 2
        j = int(np.argmax(p)); jm, jp = (j - 1) % len(p), (j + 1) % len(p)
        dd = 0.5 * (p[jm] - p[jp]) / (p[jm] - 2 * p[j] + p[jp] + 1e-30)
        errs.append(((fine[j] + dd * 0.05) - s_ + 180) % 360 - 180)
        slls.append(lobe(fine, np.sqrt(p), circ=True)[1])
    e = np.array(errs); big = np.abs(e) > 5
    return dict(rms=float(np.sqrt(np.mean(e[~big] ** 2))), max=float(np.abs(e[~big]).max()), n_wrong=int(big.sum()),
                sll_db=float(np.max(slls)))


res["df_vs_band"] = [dict(rel=rel, nf=len(sub(rel)), **df(sub(rel))) for rel in RELS]
np.savez(os.path.join(folder, "wb_cuts.npz"), **CUTS)
json.dump(res, open(os.path.join(folder, "wb_result.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

print(f"{res['folder']}: шаг {res['d_lam']:.2f}λ на 10 ГГц ({d_mm:.1f} мм)")
print("  фазовая ошибка по частоте: " + "; ".join(f"{v['f']:g}: {v['phase_err_rms_deg']:.2f}°" for v in res["by_f"]))
for v in res["by_f"]:
    if "global_dev_mean" in v:
        print(f"  {v['f']:g} ГГц (шаг {v['d_lam_f']:.2f}λ): единый центр откл. {v['global_dev_mean']:.3f} мм, остаток {v['global_rms_deg']:.2f}°")
print(f"  центр задержки: {res['delay']}")
for row in res["image_vs_band"]:
    print(f"  полоса {row['rel']*100:.0f} % ({row['nf']} ч.): " + "; ".join(
        f"{k} разр. {row[k]['range_res_mm']:.0f} мм, УБЛ дальн. {row[k]['sll_r_db']:.1f}, луч {row[k]['az_res']:.1f}°, УБЛ аз. {row[k]['sll_az_db']:.1f} дБ"
        for k in ("formula", "cst")))
for row in res["image_1f"]:
    print(f"  {row['f']:g} ГГц (шаг {row['d_lam_f']:.2f}λ): УБЛ аз. CST {row['cst']['sll_az_db']:.1f}, формула {row['formula']['sll_az_db']:.1f} дБ")
for row in res["df_vs_band"]:
    print(f"  пеленг, полоса {row['rel']*100:.0f} %: макс. ошибка {row['max']:.3f}°, ложных {row['n_wrong']}, УБЛ спектра {row['sll_db']:.1f} дБ")
