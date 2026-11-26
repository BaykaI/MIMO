"""Разбор одного расчёта конуса (build_cone2.py): виртуальные фазовые центры пар на каждой частоте и изображение
по Чапурскому при разной ширине полосы (полоса симметрично вокруг 10 ГГц, из частот расчёта).

  python cst/analyze_cone2.py cst/<папка> [полосы через запятую, %]      напр. ... cone_pt_wb 0,10,20,30,40,50
Выход: <папка>/c2_result.json, <папка>/c2_cuts.npz.

1. Виртуальный центр пары n-k на частоте f: МНК-подгонка фазы F_n·F_k моделью k·(r·u) + φ0 в полосе углов места ±30°
   по всем азимутам (как в отчёте по конусу), старт из r_tn + r_rk. Отклонение от r_tn + r_rk по 252 парам.
2. Центр задержки пары (широкополосный приём): наклон фазы F_n·F_k·exp(−jk(r_tn + r_rk)·u) по частоте в каждом
   направлении даёт добавочную задержку; модель c·τ = c·τ0 − δr·u, δr – смещение центра задержки.
3. Изображение по Чапурскому (ф. (11)-(12), дальняя зона): цель на горизонте в азимуте 0°, сигналы пар из CST,
   опорные фазы по r_tn + r_rk. Ψ(Δ, u) = |Σ_f Σ_nk x_nk(f)·exp(−jk(r_tn + r_rk)·u)·exp(j2kΔ)|. Сечения через цель:
   по азимуту (горизонт), по углу места (φ = 0), по дальности. То же по формулам (идеальные точечные излучатели).
"""
import sys, os, re, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from phase_center_fit import fit_phase_center, C0

folder = sys.argv[1]
BANDS = [float(x) / 100 for x in (sys.argv[2] if len(sys.argv) > 2 else "0").split(",")]
F0 = 10.0
kf = lambda f: 2 * np.pi * f * 1e9 / C0 / 1e3
lay = [p.split() for p in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")][1:]
pos = {p[0]: np.array(list(map(float, p[1:4]))) for p in lay}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]
pairs = [(t, r) for t in tx for r in rx]; RV = np.array([pos[t] + pos[r] for t, r in pairs])
FREQS = [float(x) for x in open(os.path.join(folder, "freqs.txt")).read().split()]


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
band = np.abs(th - 90) <= 30 + 1e-9
d_rx = float(2 * np.pi * np.linalg.norm(pos["rx1"][:2]) / 10)          # шаг по дуге нижнего кольца (10 приёмников)
res = dict(folder=os.path.basename(folder), freqs=FREQS, d_rx_mm=d_rx, d_rx_lam=d_rx / 29.979,
           R1_mm=float(np.linalg.norm(pos["rx1"][:2])), H_mm=float(max(p[2] for p in pos.values())))

# ---- 1. виртуальные центры на каждой частоте
EX = [("tx1", "rx1"), ("tx5", "rx6"), ("tx12", "rx15")]          # пары для показа координат
res["vpc"] = []; ex_xyz = {f"{t}·{r}": [] for t, r in EX}
for f in FREQS:
    k = kf(f); dev = []; rms = []
    for (t, r), rv in zip(pairs, RV):
        E = P[(t, f)] * P[(r, f)]
        xyz, e = fit_phase_center(th[band], ph[band], E[band], k, 180.0, (90.0, 0.0), init=rv)
        dev.append(np.linalg.norm(xyz - rv)); rms.append(e)
        if (t, r) in EX:
            ex_xyz[f"{t}·{r}"].append((xyz - rv).tolist())
    res["vpc"].append(dict(f=f, dev_mean=float(np.mean(dev)), dev_max=float(np.max(dev)), rms_mean=float(np.mean(rms))))
    print(f"{f:g} ГГц: центры откл. ср. {np.mean(dev):.3f} мм, макс {np.max(dev):.3f} мм, остаток {np.mean(rms):.2f}°", flush=True)
res["vpc_examples"] = dict(pairs=list(ex_xyz), rv={f"{t}·{r}": (pos[t] + pos[r]).tolist() for t, r in EX}, dxyz=ex_xyz)

# ---- 2. центр задержки
ut = u[band]


def delay_fit(freqs):
    Fh = np.array(freqs) * 1e9; out = []
    for (t, r), rv in zip(pairs, RV):
        H = np.array([(P[(t, f)] * P[(r, f)])[band] * np.exp(-1j * kf(f) * (ut @ rv)) for f in freqs])
        ps = np.unwrap(np.angle(H), axis=0)
        slope = np.linalg.lstsq(np.column_stack([np.ones(len(Fh)), Fh]), ps, rcond=None)[0][1]
        L = -slope / (2 * np.pi) * C0 * 1e3; w = np.mean(np.abs(H) ** 2, axis=0)
        M = np.column_stack([np.ones(len(ut)), -ut])
        sol = np.linalg.lstsq(M * np.sqrt(w)[:, None], L * np.sqrt(w), rcond=None)[0]
        out.append((sol[0], np.linalg.norm(sol[1:]), np.sqrt(np.average((L - M @ sol) ** 2, weights=w))))
    o = np.array(out)
    return dict(tau0_mm=float(o[:, 0].mean()), tau0_spread_mm=float(np.ptp(o[:, 0])), dr_mean=float(o[:, 1].mean()),
                dr_max=float(o[:, 1].max()), resid_mm=float(o[:, 2].mean()))


sub = lambda rel: [f for f in FREQS if abs(f - F0) <= rel * F0 / 2 + 1e-9]
res["delay"] = [dict(rel=rel, nf=len(sub(rel)), **delay_fit(sub(rel))) for rel in BANDS if len(sub(rel)) >= 3]
TAU = (res["delay"][-1] if res["delay"] else dict(tau0_mm=0.0))["tau0_mm"]
print("центр задержки:", res["delay"], flush=True)

# ---- 3. изображение
i0 = int(np.where((th == 90) & (ph == 0))[0][0])
AZ = np.arange(-180, 180, 0.25); EL = np.arange(-90, 90.01, 0.25); DR = np.arange(-200, 200.1, 1.0)
uaz = np.column_stack([np.cos(np.radians(AZ)), np.sin(np.radians(AZ)), np.zeros(AZ.size)])
uel = np.column_stack([np.cos(np.radians(EL)), np.zeros(EL.size), np.sin(np.radians(EL))])
u0 = np.array([1.0, 0, 0])


def image(freqs, src):
    a_az = np.zeros(len(AZ), complex); a_el = np.zeros(len(EL), complex); a_r = np.zeros(len(DR), complex)
    for f in freqs:
        k = kf(f)
        x = (np.array([P[(t, f)][i0] * P[(r, f)][i0] for t, r in pairs]) * np.exp(1j * k * TAU) if src == "cst"
             else np.exp(1j * k * (RV @ u0)))
        a_az += np.exp(-1j * k * (uaz @ RV.T)) @ x
        a_el += np.exp(-1j * k * (uel @ RV.T)) @ x
        a_r += np.exp(1j * 2 * k * DR) * (np.exp(-1j * k * (RV @ u0)) @ x)
    return np.abs(a_az), np.abs(a_el), np.abs(a_r)


def lobe(x, p, circ=False):
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
    flat = (not circ) and i + lo <= 0 and i + hi >= n - 1
    return (float("nan") if flat else float((hi - lo) * (x[1] - x[0]))), (float(max(out)) if out else float("nan")), float(x[i])


CUTS = dict(az=AZ, el=EL, dr=DR, freqs=np.array(FREQS))
res["image"] = []
for rel in BANDS:
    fr = sub(rel); row = dict(rel=rel, nf=len(fr))
    for src in ("cst", "formula"):
        a, e, r = image(fr, src)
        bw_a, sll_a, at_a = lobe(AZ, a, True); bw_e, sll_e, at_e = lobe(EL, e); bw_r, sll_r, at_r = lobe(DR, r)
        if len(fr) == 1:
            bw_r = sll_r = float("nan")
        row[src] = dict(az_bw=bw_a, az_sll=sll_a, az_err=at_a, el_bw=bw_e, el_sll=sll_e, el_err=at_e, r_res=bw_r, r_sll=sll_r)
        for nm, v in (("az", a), ("el", e), ("rng", r)):
            CUTS[f"b{int(round(rel * 100))}_{src}_{nm}"] = v / v.max()
    res["image"].append(row)
    c = row["cst"]; g = row["formula"]
    print(f"полоса {rel*100:.0f} % ({len(fr)} ч.): азимут {c['az_bw']:.2f}°/{c['az_sll']:.1f} дБ (форм. {g['az_bw']:.2f}°/{g['az_sll']:.1f}), "
          f"место {c['el_bw']:.2f}°/{c['el_sll']:.1f} (форм. {g['el_bw']:.2f}/{g['el_sll']:.1f}), дальность {c['r_res']:.0f} мм/{c['r_sll']:.1f} дБ "
          f"(форм. {g['r_res']:.0f}/{g['r_sll']:.1f})", flush=True)
np.savez(os.path.join(folder, "c2_cuts.npz"), **CUTS)
json.dump(res, open(os.path.join(folder, "c2_result.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("ok", folder)
