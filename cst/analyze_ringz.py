"""Разбор одного кольца из CST (build_ringz.py) по двум вопросам.

  python cst/analyze_ringz.py cst/ringz_d0.5
Выход: <папка>/ringz_result.json (все числа), печать сводки.

Что считается (12 Rx + 6 Tx, вертикальные диполи, горизонт, 13 частот 8,5-11,5 ГГц):
 A. Связь: |S| соседних приёмников и приёмника с ближайшим передатчиком, |S11| по полосе.
 B. Фазовые центры 72 пар на каждой частоте (подгонка фазы F_tx·F_rx в полосе углов места ±30° по всем азимутам)
    против r_tx + r_rx.
 C. Центр задержки (широкополосный сигнал). Фаза пары после вычитания геометрии k·(r_tx + r_rx)·u разворачивается
    по частоте; её наклон даёт добавочную групповую задержку τ(u). Модель τ(u) = τ0 − δr·u/c: δr смещение центра
    задержки от r_tx + r_rx, c·τ0 собственная задержка пары (дисперсия излучателей). Отдельно разброс c·τ0 по парам.
 D. Широкополосное изображение (ОФН дальность-азимут, MIMO, все 72 пары, 13 частот), цель в азимуте 17°:
    сигналы из CST, обработка тремя моделями: геометрия (Чапурский), геометрия + калибровка каналов по одному
    направлению на каждой частоте, полная модель CST. Плюс идеал: сигналы и обработка по формулам.
 E. Пеленг пассивным кольцом Rx (12 эл.) по всему кругу: геометрическая модель, 1 частота и вся полоса.
"""
import sys, os, re, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools")); sys.path.insert(0, HERE)
from phase_center_fit import fit_phase_center, C0

folder = sys.argv[1]
lay = [p.split() for p in open(os.path.join(folder, "layout_mm.txt"), encoding="utf-8")][1:]
pos = {p[0]: np.array(list(map(float, p[1:4]))) for p in lay}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]
# точные частоты мониторов build_ringz.py (8,5-11,5 ГГц шаг 0,25); в именах файлов они округлены до десятых
FREQS = [round(8.5 + 0.25 * i, 2) for i in range(13)]


def read(name, f):
    d = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, f"ff_{name}_{f:.1f}GHz.txt" if f"{f:.1f}" == f"{f}" else f"ff_{name}_{f}GHz.txt"), encoding="utf-8")
                  if re.match(r"^\s*[-+]?\d", l)])
    return d[:, 0], d[:, 1], 10 ** (d[:, 3] / 20) * np.exp(1j * np.radians(d[:, 4]))


def fname(name, f):
    for cand in (f"ff_{name}_{f:.1f}GHz.txt",):
        if os.path.exists(os.path.join(folder, cand)):
            return cand
    raise FileNotFoundError(name, f)


def read(name, f):
    d = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, fname(name, f)), encoding="utf-8")
                  if re.match(r"^\s*[-+]?\d", l)])
    return d[:, 0], d[:, 1], 10 ** (d[:, 3] / 20) * np.exp(1j * np.radians(d[:, 4]))


P = {}
for f in FREQS:
    for n in pos:
        th, ph, P[(n, f)] = read(n, f)
thr, phr = np.radians(th), np.radians(ph)
u = np.column_stack([np.sin(thr) * np.cos(phr), np.sin(thr) * np.sin(phr), np.cos(thr)])
F = np.array(FREQS) * 1e9
i_f0 = int(np.argmin(np.abs(F - 10e9)))
res = dict(folder=os.path.basename(folder), freqs=FREQS, R_rx=float(np.linalg.norm(pos["rx1"][:2])),
           d_mm=float(2 * np.pi * np.linalg.norm(pos["rx1"][:2]) / len(rx)))
res["d_lam"] = res["d_mm"] / (C0 / 10e9 * 1e3)

# ---------------- A. связь
sp = [f_ for f_ in os.listdir(folder) if f_.endswith(f".s{len(pos)}p")][0]
raw = [l for l in open(os.path.join(folder, sp), encoding="utf-8", errors="replace") if l.strip() and not l.startswith("!")]
opt = [l for l in raw if l.startswith("#")][0].upper().split()
nums = np.array([float(v) for l in raw if not l.startswith("#") for v in l.split()]); n_ = len(pos)
S = nums.reshape(-1, 1 + 2 * n_ * n_); fS = S[:, 0] * {"HZ": 1, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}[opt[1]]; S = S[:, 1:].reshape(-1, n_, n_, 2)
Sc = S[..., 0] * np.exp(1j * np.radians(S[..., 1])) if opt[3] == "MA" else 10 ** (S[..., 0] / 20) * np.exp(1j * np.radians(S[..., 1]))
names = rx + tx; ix = {n: i for i, n in enumerate(names)}
j0 = int(np.argmin(np.abs(fS - 10e9)))
dist_tx = {t: np.linalg.norm(pos[t] - pos["rx1"]) for t in tx}; tn = min(dist_tx, key=dist_tx.get)
inb = (fS >= 8.5e9) & (fS <= 11.5e9)
res["coupling"] = dict(rx_rx_db=float(20 * np.log10(abs(Sc[j0, ix["rx2"], ix["rx1"]]))),
                       rx_tx_db=float(20 * np.log10(abs(Sc[j0, ix[tn], ix["rx1"]]))),
                       s11_db_10=float(20 * np.log10(abs(Sc[j0, 0, 0]))),
                       s11_db_worst_band=float(20 * np.log10(np.abs(Sc[inb, 0, 0]).max())))

# ---------------- B. фазовые центры по частоте: три меры
# B1. остаточная фазовая ошибка геометрической модели: arg(F_tx·F_rx·exp(−jk r_v·u)) минус лучшая постоянная,
#     скз по полосе углов места ±30° на всех азимутах с весом |F|²; в мм: δ = Δψ/k.
# B2. «местный» фазовый центр: подгонка в секторе ±30° вокруг 12 направлений по горизонту (через 30°);
#     берётся только смещение ПОПЕРЁК направления взгляда (вдоль него узкий сектор координату не определяет),
#     именно оно даёт ошибку пеленга. Секторы, где сигнал пары ниже −10 дБ от её максимума, пропускаются.
# B3. единый центр по всему кругу (как для кольца 5λ и конуса) с остатком подгонки: показывает, есть ли у пары
#     один фазовый центр для всех направлений.
band = np.abs(th - 90) <= 30
res["vpc_by_f"] = []
DIRS = np.arange(0, 360, 30.0)
for fi, f in enumerate(FREQS):
    k = 2 * np.pi * f * 1e9 / C0 / 1e3
    perr, loc_dev, loc_dev_max, glob_dev, glob_rms = [], [], [], [], []
    for t in tx:
        for r in rx:
            rv = pos[t] + pos[r]; E = P[(t, f)] * P[(r, f)]
            z = (E * np.exp(-1j * k * (u @ rv)))[band]; w = np.abs(z) ** 2
            c = np.angle(np.sum(w * np.exp(1j * np.angle(z))))                       # лучшая постоянная фаза
            e = np.angle(z * np.exp(-1j * c)); perr.append(np.degrees(np.sqrt(np.average(e ** 2, weights=w))))
            if fi % 6 == 0:                                                         # 8,5 / 10 / 11,5 ГГц: местные центры
                dd = []
                hz = th == 90; amax = np.abs(E[hz]).max()
                for dct in DIRS:
                    a_c = np.abs(E[hz & (ph == dct)])[0]
                    if a_c < amax * 10 ** (-10 / 20):        # сектор в провале диаграммы пары (ниже −10 дБ): фаза не определена
                        continue
                    xyz, rms = fit_phase_center(th, ph, E, k, 30.0, (90.0, float(dct)))
                    tr = np.array([-np.sin(np.radians(dct)), np.cos(np.radians(dct)), 0.0])   # поперёк направления взгляда
                    dd.append(abs((xyz - rv) @ tr))
                if dd:
                    loc_dev.append(np.mean(dd)); loc_dev_max.append(np.max(dd))
                # единый центр: старт от r_tx + r_rx и от развёрнутой фазы, берётся лучший по остатку
                c1 = fit_phase_center(th[band], ph[band], E[band], k, 180.0, (90.0, 0.0), init=rv)
                c2 = fit_phase_center(th[band], ph[band], E[band], k, 180.0, (90.0, 0.0))
                xyz, rms = min((c1, c2), key=lambda c: c[1])
                glob_dev.append(np.linalg.norm(xyz - rv)); glob_rms.append(rms)
    perr = np.array(perr)
    row = dict(f=f, phase_err_rms_deg=float(perr.mean()), phase_err_max_deg=float(perr.max()),
               phase_err_mm=float(np.radians(perr.mean()) / k))
    if loc_dev:
        row.update(local_dev_mean=float(np.mean(loc_dev)), local_dev_max=float(np.max(loc_dev_max)),
                   global_dev_mean=float(np.mean(glob_dev)), global_rms_deg=float(np.mean(glob_rms)))
    res["vpc_by_f"].append(row)

# ---------------- C. центр задержки и собственная задержка пар
hor = np.abs(th - 90) <= 30
ut = u[hor]
gd = []
for t in tx:
    for r in rx:
        rv = pos[t] + pos[r]
        H = np.array([(P[(t, f)] * P[(r, f)])[hor] * np.exp(-1j * 2 * np.pi * f * 1e9 / C0 / 1e3 * (ut @ rv)) for f in FREQS])   # (Nf, Ndir)
        ps = np.unwrap(np.angle(H), axis=0)
        A = np.column_stack([np.ones(len(F)), F])
        slope = np.linalg.lstsq(A, ps, rcond=None)[0][1]           # dψ/df, рад/Гц, по каждому направлению
        L = -slope / (2 * np.pi) * C0 * 1e3                       # добавочная групповая задержка как путь c·τ, мм
        w = np.mean(np.abs(H) ** 2, axis=0)
        M = np.column_stack([np.ones(len(ut)), -ut[:, 0], -ut[:, 1], -ut[:, 2]])   # c·τ = c·τ0 − δr·u  (мм)
        sol = np.linalg.lstsq(M * np.sqrt(w)[:, None], L * np.sqrt(w), rcond=None)[0]
        resid = L - M @ sol
        gd.append(dict(pair=f"{t}·{r}", tau0_mm=float(sol[0]), dr=list(map(float, sol[1:])),
                       dr_norm=float(np.linalg.norm(sol[1:3])), resid_mm=float(np.sqrt(np.average(resid ** 2, weights=w)))))
t0 = np.array([g["tau0_mm"] for g in gd]); drn = np.array([g["dr_norm"] for g in gd])
res["delay"] = dict(tau0_mean_mm=float(t0.mean()), tau0_spread_mm=float(t0.max() - t0.min()), tau0_std_mm=float(t0.std()),
                    dr_mean=float(drn.mean()), dr_max=float(drn.max()), resid_mm=float(np.mean([g["resid_mm"] for g in gd])))

# ---------------- D. широкополосное изображение дальность-азимут (MIMO)
sel = th == 90
azg = ph[sel]; order = np.argsort(azg); azg = azg[order]
Hcst = {}                      # (Nf, Npairs, Naz) на горизонте
pairs = [(t, r) for t in tx for r in rx]
RV = np.array([pos[t] + pos[r] for t, r in pairs])
for f in FREQS:
    Hcst[f] = np.array([(P[(t, f)] * P[(r, f)])[sel][order] for t, r in pairs])
phs = 17.0; i_s = int(np.argmin(np.abs(azg - phs)))     # цель (узел сетки 1°)
dr = np.arange(-400, 400.1, 2.0); az_im = np.arange(-90, 124.01, 0.25)
U = np.column_stack([np.cos(np.radians(az_im)), np.sin(np.radians(az_im)), np.zeros(az_im.size)])
ucal = 0                                            # калибровка по азимуту 0°


def interp_manifold(f, az):
    i = np.floor(az).astype(int) % 360; tt = az - np.floor(az)
    return (1 - tt) * Hcst[f][:, i] + tt * Hcst[f][:, (i + 1) % 360]


TAU_COMMON = res["delay"]["tau0_mean_mm"]                # общая для всех пар задержка (путь), мм


def image(model, signal, freqs=None):
    img = np.zeros((len(dr), len(az_im)), complex)
    for f in (FREQS if freqs is None else freqs):
        k = 2 * np.pi * f * 1e9 / C0 / 1e3
        geo = np.exp(1j * k * (RV @ U.T))                          # (Np, Naz) модель отклика по геометрии
        if signal == "cst":
            x = Hcst[f][:, i_s]
            if model == "geo":                                     # одна общая поправка задержки на всю систему
                x = x * np.exp(1j * k * TAU_COMMON)       # ψ = −ω·τ0 убирается множителем exp(+jω·τ0)
        elif signal == "cst_raw":
            x = Hcst[f][:, i_s]
        else:
            x = np.exp(1j * k * (RV @ np.array([np.cos(np.radians(phs)), np.sin(np.radians(phs)), 0])))
        if model == "geo":
            A = geo
            if signal == "cst_raw":
                pass
        elif model == "cal":
            c = Hcst[f][:, ucal] / np.exp(1j * k * (RV @ np.array([1.0, 0, 0])))
            A = c[:, None] * geo
        else:
            A = interp_manifold(f, az_im)
        b = (A.conj().T @ x) / np.sqrt(np.sum(np.abs(A) ** 2, axis=0))   # нормировка на мощность модели (иначе максимум тянется к сильным направлениям)
        img += np.exp(1j * 2 * k * dr)[:, None] * b[None, :]        # дальность: путь туда и обратно 2Δ, цель на Δ = 0
    return np.abs(img) / np.abs(img).max()


def lobe_cut(x, p):
    """ширина по −3 дБ и наибольший боковой лепесток (за первыми минимумами) для одного сечения"""
    d = 20 * np.log10(p / p.max() + 1e-12); i = int(np.argmax(d)); lo = hi = i
    while lo > 0 and d[lo - 1] >= -3: lo -= 1
    while hi < len(d) - 1 and d[hi + 1] >= -3: hi += 1
    l = r = i
    while l > 0 and d[l - 1] < d[l]: l -= 1
    while r < len(d) - 1 and d[r + 1] < d[r]: r += 1
    out = np.concatenate([d[:l], d[r + 1:]])
    flat = (lo == 0 and hi == len(d) - 1)
    return (float("nan") if flat else float(x[hi] - x[lo])), (float(out.max()) if out.size else float("nan"))


def img_metrics(I):
    i, j = np.unravel_index(np.argmax(I), I.shape)
    rres, sll_r = lobe_cut(dr, I[:, j]); ares, sll_a = lobe_cut(az_im, I[i, :])
    return dict(az_err=float(az_im[j] - phs), range_err_mm=float(dr[i]), range_res_mm=rres, az_res=ares,
                sll_az_db=sll_a, sll_r_db=sll_r)


res["image"] = {}
CUTS = dict(dr=dr, az=az_im)
for key, model, signal in (("формулы (идеал)", "geo", "geo"), ("CST, геометрия без поправки задержки", "geo", "cst_raw"),
                           ("CST, обработка по геометрии", "geo", "cst"),
                           ("CST, геометрия + калибровка каналов", "cal", "cst"), ("CST, полная модель CST", "cst", "cst")):
    I = image(model, signal); res["image"][key] = img_metrics(I)
    i0 = int(np.argmin(np.abs(dr))); j0 = int(np.argmin(np.abs(az_im - phs)))
    CUTS[f"img_{model}_{signal}_range"] = I[:, j0]; CUTS[f"img_{model}_{signal}_az"] = I[i0, :]
# зависимость от ширины полосы (симметрично вокруг 10 ГГц)
res["image_vs_band"] = []
for lo_, hi_ in ((10.0, 10.0), (9.75, 10.25), (9.5, 10.5), (9.0, 11.0), (8.5, 11.5)):
    fr = [f for f in FREQS if lo_ - 1e-9 <= f <= hi_ + 1e-9]
    row = dict(band=f"{lo_:g}-{hi_:g}", rel=float((hi_ - lo_) / 10.0), nf=len(fr))
    for key, model, signal in (("formula", "geo", "geo"), ("cst_geo", "geo", "cst"), ("cst_cal", "cal", "cst")):
        I = image(model, signal, fr); row[key] = img_metrics(I)
        CUTS[f"band{len(fr)}_{key}_range"] = I[:, int(np.argmin(np.abs(az_im - phs)))]
    res["image_vs_band"].append(row)
res["image_axes"] = dict(dr=[float(dr[0]), float(dr[-1]), len(dr)], az=[float(az_im[0]), float(az_im[-1]), len(az_im)], target_az=float(azg[i_s]))

# пример: остаточная фаза пары tx1·rx1 (после вычитания k·r_v·u) по частоте в направлениях через 30° по горизонту
rv = pos["tx1"] + pos["rx1"]; ex = []
for dct in np.arange(0, 360, 30.0):
    ii = np.where((th == 90) & (ph == dct))[0][0]
    z = [P[("tx1", f)][ii] * P[("rx1", f)][ii] * np.exp(-1j * 2 * np.pi * f * 1e9 / C0 / 1e3 * (u[ii] @ rv)) for f in FREQS]
    ex.append(np.unwrap(np.angle(z)))
CUTS["pair_resid_phase"] = np.degrees(np.array(ex)); CUTS["pair_dirs"] = np.arange(0, 360, 30.0); CUTS["freqs"] = np.array(FREQS)
np.savez(os.path.join(folder, "cuts.npz"), **CUTS)

# ---------------- E. пеленг пассивным кольцом Rx
fine = np.arange(0, 360, 0.05); uf = np.column_stack([np.cos(np.radians(fine)), np.sin(np.radians(fine)), np.zeros(fine.size)])
R = np.array([pos[r] for r in rx])


def df_err(freqs):
    errs = []
    for s_ in np.arange(0, 360, 3.0):
        p = np.zeros(len(fine))
        for f in freqs:
            k = 2 * np.pi * f * 1e9 / C0 / 1e3
            Hr = np.array([P[(r, f)][sel][order] for r in rx]); i = int(s_) % 360
            x = Hr[:, i]; A = np.exp(1j * k * (R @ uf.T))
            p += np.abs(A.conj().T @ x) ** 2
        j = int(np.argmax(p)); jm, jp = (j - 1) % len(p), (j + 1) % len(p)
        dd = 0.5 * (p[jm] - p[jp]) / (p[jm] - 2 * p[j] + p[jp] + 1e-30)
        e = ((fine[j] + dd * 0.05) - s_ + 180) % 360 - 180; errs.append(e)
    e = np.array(errs); big = np.abs(e) > 5
    return dict(rms=float(np.sqrt(np.mean(e[~big] ** 2))) if (~big).any() else float("nan"), max=float(np.abs(e[~big]).max()) if (~big).any() else float("nan"),
                n_wrong=int(big.sum()), n=len(e))


res["df"] = {"1 частота 10 ГГц": df_err([FREQS[i_f0]]), "полоса 8,5-11,5 ГГц": df_err(FREQS)}

json.dump(res, open(os.path.join(folder, "ringz_result.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"{res['folder']}: d = {res['d_lam']:.2f}λ ({res['d_mm']:.1f} мм), R = {res['R_rx']:.1f} мм")
print(f"  связь: Rx-Rx {res['coupling']['rx_rx_db']:.1f} дБ, Rx-Tx {res['coupling']['rx_tx_db']:.1f} дБ, S11(10 ГГц) {res['coupling']['s11_db_10']:.1f} дБ, худший в полосе {res['coupling']['s11_db_worst_band']:.1f} дБ")
print("  фазовая ошибка геометрии по частоте (скз, °): " + "; ".join(f"{v['f']:g}: {v['phase_err_rms_deg']:.1f}" for v in res["vpc_by_f"]))
for v in res["vpc_by_f"]:
    if "local_dev_mean" in v:
        print(f"  {v['f']:g} ГГц: местный центр (сектор ±30°) откл. в среднем {v['local_dev_mean']:.2f} мм, макс {v['local_dev_max']:.2f} мм; единый центр по кругу откл. {v['global_dev_mean']:.2f} мм при остатке {v['global_rms_deg']:.1f}°")
dl = res["delay"]
print(f"  центр задержки: смещение от r_tx+r_rx среднее {dl['dr_mean']:.2f} мм, макс {dl['dr_max']:.2f} мм, остаток {dl['resid_mm']:.2f} мм; c·τ0 пар {dl['tau0_mean_mm']:.1f} мм (в т. ч. опорная добавка CST), разброс {dl['tau0_spread_mm']:.2f} мм (СКО {dl['tau0_std_mm']:.2f})")
for k_, v_ in res["image"].items():
    print(f"  изображение, {k_}: ошибка азимута {v_['az_err']:+.2f}°, дальности {v_['range_err_mm']:+.0f} мм, разрешение {v_['range_res_mm']:.0f} мм × {v_['az_res']:.1f}°, УБЛ по азимуту {v_['sll_az_db']:.1f} дБ, по дальности {v_['sll_r_db']:.1f} дБ")
for row in res["image_vs_band"]:
    print(f"  полоса {row['band']} ГГц ({row['rel']*100:.0f} %, {row['nf']} част.): " + "; ".join(f"{k_} разр. {row[k_]['range_res_mm']:.0f} мм, ошибка дальн. {row[k_]['range_err_mm']:+.0f} мм, УБЛ аз. {row[k_]['sll_az_db']:.1f}" for k_ in ("formula", "cst_geo", "cst_cal")))
for k_, v_ in res["df"].items():
    print(f"  пеленг кольцом Rx, {k_}: скз {v_['rms']:.3f}°, макс {v_['max']:.3f}°, ложных пеленгов {v_['n_wrong']} из {v_['n']}")
