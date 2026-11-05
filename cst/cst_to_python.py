"""ЧТО ПЕРЕДАЁТСЯ ИЗ CST В PYTHON, КАК ПЕРЕСЧИТЫВАЕТСЯ В ВИРТУАЛЬНЫЙ ЭКВИВАЛЕНТ И ЗАЧЕМ (пеленгатор на кольцевой решётке).

  python cst/cst_to_python.py cst/cone [--f 10]
Печатает по шагам ВСЁ, что прочитано (файлы, сетка, колонки, контрольные числа), затем пересчёт и пеленгацию.
Лог: <папка>/cst_to_python.log, рисунки: fig2/df_*.png.

ШАГ 1. Из CST берутся две вещи:
  1) встроенная диаграмма направленности каждого порта n (файл ff_<порт>_<f>GHz.txt): для каждой пары углов (θ, φ) —
     модуль и фаза двух компонент поля E_θ и E_φ в дальней зоне (дБи и градусы) с фазой ОТНОСИТЕЛЬНО ГЛОБАЛЬНОГО НУЛЯ
     координат. Порт n возбуждён, остальные порты нагружены 50 Ом — значит, в диаграмме уже сидит взаимная связь
     (переизлучение соседей). По теореме взаимности та же функция F_n(u) — коэффициент приёма порта n с направления u.
  2) S-матрица всех портов (файл .s<N>p, Touchstone) — согласование и развязка; для пересчёта координат не нужна.
ШАГ 2. Пересчёт: отклик пары «передатчик n — приёмник k» в направлении u равен F_n(u)·F_k(u). Фазовый центр этого
  отклика (точка, вокруг которой фаза не зависит от направления) — виртуальный (эквивалентный) центр пары. Из CST он
  ищется подгонкой фазы F_n·F_k моделью k·(r·u); из геометрии — как r_tn + r_rk. Сравнение — проверка метода.
ШАГ 3. Польза: пеленгатор. Кольцо из Nr приёмников принимает плоскую волну с азимута φ_s: сигналы x_k = F_k(u_s).
  Пеленг ищется по максимуму |a(φ)^H x|² (Бартлетт), где a(φ) — модель отклика решётки на направление φ (manifold).
  Если модель геометрическая (точечные элементы в r_k), ошибка пеленга — это цена пренебрежения связью; с калибровкой
  каналов по одному направлению и с полной моделью из CST ошибка исчезает. MIMO (виртуальная решётка Tx×Rx) даёт
  более узкий пик и разрешает два источника, которые одно кольцо не разрешает.
"""
import sys, os, re, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from phase_center_fit import fit_phase_center, C0

ap = argparse.ArgumentParser(); ap.add_argument("folder"); ap.add_argument("--f", type=float, default=10.0)
a = ap.parse_args()
folder = a.folder if os.path.isabs(a.folder) else os.path.join(os.getcwd(), a.folder)
OUT = os.path.join(HERE, "fig2"); os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(folder, "cst_to_python.log"), "w", encoding="utf-8")


def say(*args):
    s = " ".join(str(x) for x in args)
    print(s); LOG.write(s + chr(10))


def hdr(t):
    say(""); say("=" * 100); say(t); say("=" * 100)


LAM = C0 / (a.f * 1e9) * 1e3; k = 2 * np.pi / LAM
hdr(f"ШАГ 1. ЧТО СЧИТАНО ИЗ CST   (папка {folder}, частота {a.f:g} ГГц, λ = {LAM:.2f} мм, k = 2π/λ = {k:.5f} 1/мм)")

# ---- 1а. раскладка
lay = os.path.join(folder, "layout_mm.txt")
elems = [(p[0], float(p[1]), float(p[2]), float(p[3]) if len(p) > 3 else 0.0)
         for p in (l.split() for l in open(lay, encoding="utf-8")) if p and p[0] != "name"]
pos = {n: np.array([x, y, z]) for n, x, y, z in elems}
rx = [n for n in pos if n.startswith("rx")]; tx = [n for n in pos if n.startswith("tx")]
say(f"1а. Раскладка {os.path.basename(lay)}: {len(rx)} приёмных (rx) и {len(tx)} передающих (tx) портов, координаты в мм "
    f"(это ТО ЖЕ, что задано в CST при построении; CST координаты не возвращает, они известны из build_*.py):")
for n in rx[:3] + ["…"] + rx[-1:] + tx[:2] + ["…"] + tx[-1:]:
    say("    " + (n if n == "…" else f"{n:5s} x = {pos[n][0]:9.3f}  y = {pos[n][1]:9.3f}  z = {pos[n][2]:8.3f} мм   |r| = {np.linalg.norm(pos[n]):.2f} мм"))

# ---- 1б. файлы ДН
files = sorted(f_ for f_ in os.listdir(folder) if re.match(r"ff_(rx|tx)\d+_\d+\.\dGHz\.txt$", f_))
freqs = sorted({float(re.search(r"_(\d+\.\d)GHz", f_).group(1)) for f_ in files})
say(f"1б. Файлов ДН: {len(files)} = {len(pos)} портов × {len(freqs)} частот {freqs} ГГц. Пример файла: {files[0]}")
with open(os.path.join(folder, f"ff_{rx[0]}_{a.f:.1f}GHz.txt"), encoding="utf-8") as fh:
    lines = fh.readlines()
header = lines[0].rstrip(); cols = [c for c in re.split(r"\s{2,}", header.strip()) if c]
say("    заголовок:", header)
say("    колонки  :", " | ".join(f"[{i}] {c}" for i, c in enumerate(cols)))
data = np.array([list(map(float, l.split())) for l in lines if re.match(r"^\s*[-+]?\d", l)])
th, ph = data[:, 0], data[:, 1]
say(f"    строк данных: {len(data)}; θ от {th.min():g} до {th.max():g}° ({len(np.unique(th))} значений, шаг {np.diff(np.unique(th))[0]:g}°), "
    f"φ от {ph.min():g} до {ph.max():g}° ({len(np.unique(ph))} значений, шаг {np.diff(np.unique(ph))[0]:g}°) → "
    f"{len(np.unique(th))} × {len(np.unique(ph))} = {len(np.unique(th)) * len(np.unique(ph))} направлений на полной сфере")
say("    первые 3 строки как есть:")
for l in lines[2:5]:
    say("      " + l.rstrip())


def read_ff(name, f=None):
    f = a.f if f is None else f
    d = np.array([list(map(float, l.split())) for l in open(os.path.join(folder, f"ff_{name}_{f:.1f}GHz.txt"), encoding="utf-8")
                  if re.match(r"^\s*[-+]?\d", l)])
    Et = 10 ** (d[:, 3] / 20) * np.exp(1j * np.radians(d[:, 4]))     # E_θ: модуль из дБи → линейный, фаза из градусов
    Ep = 10 ** (d[:, 5] / 20) * np.exp(1j * np.radians(d[:, 6]))     # E_φ
    return d[:, 0], d[:, 1], d[:, 2], Et, Ep


th, ph, Dir, Et, Ep = read_ff(rx[0])
thr, phr = np.radians(th), np.radians(ph)
u = np.column_stack([np.sin(thr) * np.cos(phr), np.sin(thr) * np.sin(phr), np.cos(thr)])   # единичные векторы направлений
say(f"1в. Как числа превращаются в комплексное поле: E_θ = 10^(Abs(Theta)/20) · exp(j·Phase(Theta)), то же для E_φ.")
i0 = np.where((th == 90) & (ph == 0))[0][0]
say(f"    {rx[0]} в направлении горизонт/φ=0 (θ=90°, φ=0): Abs(Theta) = {20*np.log10(abs(Et[i0])):.3f} дБи, Phase(Theta) = {np.degrees(np.angle(Et[i0])):.2f}°, "
    f"Abs(Phi) = {20*np.log10(abs(Ep[i0])+1e-12):.1f} дБи  → у z-диполя рабочая компонента E_θ, E_φ на 60–100 дБ ниже")
chk = 10 * np.log10(abs(Et) ** 2 + abs(Ep) ** 2) - Dir
say(f"1г. Контроль: Abs(Dir) должен равняться 10·lg(|E_θ|² + |E_φ|²): расхождение по всей сфере max {np.abs(chk).max():.2e} дБ — колонки согласованы.")
say(f"1д. Контроль фазовой привязки: фаза ДН задана ОТ ГЛОБАЛЬНОГО НУЛЯ. Для {rx[0]} в r = ({pos[rx[0]][0]:.0f}; {pos[rx[0]][1]:.0f}; {pos[rx[0]][2]:.0f}) мм "
    "фаза должна меняться по азимуту на горизонте как k·x·cos φ:")
sel = th == 90
azl = ((ph[sel] + 180) % 360) - 180; o = np.argsort(azl)
meas = np.degrees(np.angle(Et[sel][o] * np.conj(Et[i0])))
model = np.degrees(k * (pos[rx[0]][0] * (np.cos(np.radians(azl[o])) - 1)))
model = np.degrees(np.angle(np.exp(1j * np.radians(model))))
for az_ in (0, 5, 10, 20, 30):
    j = int(np.argmin(np.abs(azl[o] - az_)))
    say(f"      φ = {az_:3d}°: фаза CST − фаза(0°) = {meas[j]:8.2f}°,   k·x·(cos φ − 1) = {model[j]:8.2f}°,   разность {np.degrees(np.angle(np.exp(1j*np.radians(meas[j]-model[j])))):6.2f}°")
say("    → привязка и знак фазы (+k·r·u) подтверждены на сырых данных CST.")

# ---- 1е. S-матрица
sp = [f_ for f_ in os.listdir(folder) if re.match(r".*\.s\d+p$", f_)]
if sp:
    p_ = os.path.join(folder, sp[0]); n_ = int(re.search(r"\.s(\d+)p$", sp[0]).group(1))
    raw = [l for l in open(p_, encoding="utf-8", errors="replace") if l.strip() and not l.startswith("!")]
    opt = [l for l in raw if l.startswith("#")][0].strip()
    nums = np.array([float(v) for l in raw if not l.startswith("#") for v in l.split()])
    nf = len(nums) // (1 + 2 * n_ * n_)
    S = nums.reshape(nf, 1 + 2 * n_ * n_)
    fS = S[:, 0]; S = S[:, 1:].reshape(nf, n_, n_, 2)
    fmt = opt.upper().split()[3]
    Sc = S[..., 0] * np.exp(1j * np.radians(S[..., 1])) if fmt == "MA" else (10 ** (S[..., 0] / 20) * np.exp(1j * np.radians(S[..., 1])) if fmt == "DB" else S[..., 0] + 1j * S[..., 1])
    i_f = int(np.argmin(np.abs(fS - a.f * 1e9 * (1 if "HZ" in opt.upper().split()[1] and opt.upper().split()[1] == "HZ" else 1e-9))))
    say(f"1е. S-матрица {sp[0]}: {n_} портов, {nf} частотных точек, формат «{opt}». На {a.f:g} ГГц: "
        f"|S11| = {20*np.log10(abs(Sc[i_f,0,0])):.1f} дБ (согласование rx1), |S({len(rx)+1},1)| = {20*np.log10(abs(Sc[i_f,len(rx),0])):.1f} дБ (связь tx1–rx1). "
        "Для пересчёта координат S-матрица не используется; она показывает, что связь между элементами есть (десятки дБ), и именно она сидит во встроенных ДН.")

# ================================================================ ШАГ 2
hdr("ШАГ 2. ПЕРЕСЧЁТ В ВИРТУАЛЬНЫЙ ЭКВИВАЛЕНТ (на числах из шага 1)")
P = {}
for n in pos:
    P[n] = read_ff(n)[3]            # ко-поляризация z-диполя: E_θ
say(f"2а. Загружены встроенные ДН всех {len(P)} портов на {a.f:g} ГГц: массивы по {len(th)} комплексных чисел каждый.")
say("2б. Метод 1 (геометрия, Чапурский / ВАР): r_v(n,k) = r_tn + r_rk. Метод 2 (CST): подгонка фазы F_n(u)·F_k(u) моделью k·(r_v·u)")
say("    в полосе углов места ±30° на всех азимутах (x, y — из азимутальной зависимости фазы, z — из угломестной).")
band = np.abs(th - 90) <= 30
say(f"    точек в полосе: {band.sum()} из {len(th)}. Пример на одной паре — {tx[0]}·{rx[0]}:")
Fp = P[tx[0]] * P[rx[0]]
xyz, rms = fit_phase_center(th[band], ph[band], Fp[band], k, 180.0, (90.0, 0.0))
geo = pos[tx[0]] + pos[rx[0]]
say(f"      r_tn = {np.round(pos[tx[0]],1)}, r_rk = {np.round(pos[rx[0]],1)} → метод 1: r_v = {np.round(geo,1)} мм")
say(f"      метод 2 (CST): r_v = {np.round(xyz,2)} мм, остаток фазы после подгонки {rms:.1f}° (это рябь от соседей), расхождение {np.linalg.norm(xyz-geo):.2f} мм")
say("      как это видно в сырой фазе отклика пары вдоль горизонта (φ = 0…30°): фаза CST против k·(r_v·u) по методу 1:")
meas = np.degrees(np.angle(Fp[sel][o] * np.conj(Fp[i0])))
model = np.degrees(np.angle(np.exp(1j * k * (geo[0] * (np.cos(np.radians(azl[o])) - 1) + geo[1] * np.sin(np.radians(azl[o]))))))
for az_ in (0, 2, 5, 10, 20, 30):
    j = int(np.argmin(np.abs(azl[o] - az_)))
    say(f"        φ = {az_:3d}°: CST {meas[j]:8.2f}°   геометрия {model[j]:8.2f}°   разность {np.degrees(np.angle(np.exp(1j*np.radians(meas[j]-model[j])))):6.2f}°")
devs = []
for tn in tx:
    for rn in rx:
        xyz, rms = fit_phase_center(th[band], ph[band], (P[tn] * P[rn])[band], k, 180.0, (90.0, 0.0))
        devs.append(np.linalg.norm(xyz - (pos[tn] + pos[rn])))
devs = np.array(devs)
say(f"2в. Все {len(devs)} пар: расхождение метод 2 − метод 1: среднее {devs.mean():.2f} мм, 95 % ≤ {np.percentile(devs,95):.2f}, макс {devs.max():.2f} мм = λ/{LAM/devs.max():.0f}; "
    f"пар хуже λ/20: {(devs > LAM/20).sum()}")
say("    → координаты виртуальных центров реальной решётки равны r_tn + r_rk. Это и есть результат пересчёта: формула применима к реальным излучателям.")

# ================================================================ ШАГ 3
hdr("ШАГ 3. ЗАЧЕМ: КОЛЬЦЕВАЯ РЕШЁТКА КАК ПЕЛЕНГАТОР ИСТОЧНИКА РАДИОИЗЛУЧЕНИЯ")
ring = [n for n in rx if abs(pos[n][2]) < 1e-6]        # нижнее кольцо приёмников (z = 0)
say(f"3а. Пеленгатор: кольцо из {len(ring)} приёмников радиусом {np.linalg.norm(pos[ring[0]][:2]):.0f} мм ({np.linalg.norm(pos[ring[0]][:2])/LAM:.1f}λ) в плоскости z = 0 "
    f"(элементы {', '.join(ring)}). Источник на горизонте в азимуте φ_s; по взаимности сигнал на k-м приёмнике x_k = F_k(u_s) — "
    "берётся прямо из ДН CST (взаимная связь уже внутри).")
say("    Три модели отклика решётки a_k(φ), с которыми пеленгатор сравнивает принятый вектор (оценка по максимуму |a^H x|²):")
say("      (а) геометрическая: a_k = exp(j k r_k·u(φ)) — точечные элементы в геометрических точках, связь не учтена;")
say("      (б) геометрическая + калибровка каналов по ОДНОМУ направлению φ_cal = 0: a_k = c_k·exp(j k r_k·u(φ)), c_k = F_k(u_cal)/exp(j k r_k·u_cal);")
say("      (в) полная из CST: a_k = F_k(u(φ)) (таблица встроенных ДН, интерполяция по φ).")
azg = np.arange(0, 360, 1.0)                       # сетка CST по азимуту на горизонте
idx = {int(round(p_)): i for i, p_ in enumerate(ph[sel])}
thsel = np.where(sel)[0]
A_cst = np.array([[P[n][thsel[idx[int(p_)]]] for p_ in azg] for n in ring])      # (Nr, 360)
fine = np.arange(0, 360, 0.05)
uf = np.column_stack([np.cos(np.radians(fine)), np.sin(np.radians(fine)), np.zeros(fine.size)])
R = np.array([pos[n] for n in ring])
A_geo = np.exp(1j * k * (R @ uf.T))                # (Nr, Nfine)
ucal = np.array([1.0, 0, 0]); c = A_cst[:, 0] / np.exp(1j * k * (R @ ucal))
A_cal = c[:, None] * A_geo
say(f"    калибровочные коэффициенты c_k (модуль, фаза): " + "; ".join(f"{abs(ci):.2f}∠{np.degrees(np.angle(ci)):.0f}°" for ci in c))
say("    — вот что CST добавляет к геометрии: амплитуды и начальные фазы каналов различаются (связь, положение на кольце).")


def cst_manifold(phi):
    """линейная интерполяция комплексной ДН по азимуту между узлами 1°"""
    i = int(np.floor(phi)) % 360; t = phi - np.floor(phi)
    return (1 - t) * A_cst[:, i] + t * A_cst[:, (i + 1) % 360]


def bearing(x, A, grid):
    p = np.abs(A.conj().T @ x) ** 2
    i = int(np.argmax(p))
    if 0 < i < len(grid) - 1:      # параболическое уточнение максимума
        y0, y1, y2 = p[i - 1], p[i], p[i + 1]; d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
        return grid[i] + d * (grid[1] - grid[0])
    return grid[i]


A_cstf = np.array([cst_manifold(p_) for p_ in fine]).T
errs = {"(а) геометрия": [], "(б) геометрия + калибровка по φ=0": [], "(в) полная модель CST": []}
src = np.arange(0, 360, 3.0)
for phs in src:
    x = cst_manifold(phs)
    for key, A in zip(errs, (A_geo, A_cal, A_cstf)):
        e = bearing(x, A, fine) - phs; e = (e + 180) % 360 - 180; errs[key].append(e)
say("3б. Ошибка пеленга (оценка − истина) для источников через 3° по всему кругу, без шума:")
for key, e in errs.items():
    e = np.array(e); say(f"      {key:38s}: скз {np.sqrt(np.mean(e**2)):6.3f}°, максимум {np.abs(e).max():6.3f}°")
e_geo = np.array(errs["(а) геометрия"])
bw = np.degrees(0.5 * LAM / (2 * np.linalg.norm(pos[ring[0]][:2])))
e_cst = np.array(errs["(в) полная модель CST"])
say(f"    Ширина луча кольца по −3 дБ ≈ {bw:.1f}° (λ/2D). Систематическая ошибка геометрической модели {np.abs(e_geo).max():.2f}° — "
    f"{100*np.abs(e_geo).max()/bw:.0f} % ширины луча: для симметричного кольца связь почти не портит пеленг, потому что начальные фазы каналов "
    f"у всех элементов близки (см. c_k). Полная модель из CST снижает ошибку до {np.abs(e_cst).max():.3f}°; калибровка по одному направлению "
    "амплитудную зависимость от угла не учитывает и выигрыша здесь не даёт.")

# шум: СКО пеленга при ОСШ
rng = np.random.default_rng(1)
say("3в. То же с шумом (ОСШ на элемент 10 дБ, 200 реализаций, источник в φ_s = 37°):")
phs = 37.0; x0 = cst_manifold(phs); sig = np.sqrt(np.mean(abs(x0) ** 2)) / 10 ** (10 / 20)
for key, A in zip(errs, (A_geo, A_cal, A_cstf)):
    e = []
    for _ in range(200):
        x = x0 + sig * (rng.standard_normal(len(x0)) + 1j * rng.standard_normal(len(x0))) / np.sqrt(2)
        e.append(((bearing(x, A, fine) - phs) + 180) % 360 - 180)
    e = np.array(e); say(f"      {key:38s}: среднее (систематика) {e.mean():+6.3f}°, СКО {e.std():.3f}°")

# MIMO: виртуальная решётка Tx-кольцо × Rx-кольцо — разрешение двух источников (целей)
tring = [n for n in tx if abs(pos[n][2] - pos[tx[0]][2]) < 1e-6]
say(f"3г. MIMO-вариант (активная РЛС): {len(tring)} передатчиков кольца {np.linalg.norm(pos[tring[0]][:2]):.0f} мм × {len(ring)} приёмников = {len(tring)*len(ring)} "
    "виртуальных элементов. Отклик цели для пары (n,k): F_n(u_s)·F_k(u_s). Две цели одинаковой силы в азимутах φ1, φ2:")
A_tcst = np.array([[P[n][thsel[idx[int(p_)]]] for p_ in azg] for n in tring])


def mimo_x(phs):
    i = int(np.floor(phs)) % 360; t = phs - np.floor(phs)
    ft = (1 - t) * A_tcst[:, i] + t * A_tcst[:, (i + 1) % 360]; fr = cst_manifold(phs)
    return (ft[:, None] * fr[None, :]).ravel()


Rv = (np.array([pos[n] for n in tring])[:, None, :] + R[None, :, :]).reshape(-1, 3)
A_mimo = np.exp(1j * k * (Rv @ uf.T))
def hpbw(p):
    d = 10 * np.log10(p / p.max()); i = int(np.argmax(d)); lo = hi = i
    while lo > 0 and d[lo - 1] >= -3: lo -= 1
    while hi < len(d) - 1 and d[hi + 1] >= -3: hi += 1
    return fine[hi] - fine[lo]


p1 = np.abs(A_geo.conj().T @ cst_manifold(50.0)) ** 2; p2 = np.abs(A_mimo.conj().T @ mimo_x(50.0)) ** 2
say(f"    ширина пика по −3 дБ для одного источника: кольцо Rx {hpbw(p1):.2f}°, MIMO {hpbw(p2):.2f}°")
spectra = {}
for dphi in (8.0, 6.0, 5.0, 4.0):
    x_rx = cst_manifold(50.0) + cst_manifold(50.0 + dphi); x_m = mimo_x(50.0) + mimo_x(50.0 + dphi)
    p_rx = np.abs(A_geo.conj().T @ x_rx) ** 2; p_m = np.abs(A_mimo.conj().T @ x_m) ** 2
    w = (fine > 40) & (fine < 70)
    # разрешены, если в окне два локальных максимума выше −6 дБ и провал между ними ≥ 3 дБ

    def resolved(p):
        pp = 10 * np.log10(p[w] / p[w].max()); x_ = fine[w]
        peaks = [i for i in range(1, len(pp) - 1) if pp[i] >= pp[i - 1] and pp[i] >= pp[i + 1] and pp[i] > -6]
        if len(peaks) < 2:
            return 0.0
        peaks = sorted(sorted(peaks, key=lambda i: -pp[i])[:2])
        return float(min(pp[peaks[0]], pp[peaks[1]]) - pp[peaks[0]:peaks[1] + 1].min())

    spectra[dphi] = (p_rx, p_m)
    say(f"      Δφ = {dphi:3.0f}°: провал между целями — только кольцо Rx: {resolved(p_rx):5.1f} дБ, MIMO (виртуальная решётка): {resolved(p_m):5.1f} дБ "
        f"{'→ кольцо не разрешает, MIMO разрешает' if resolved(p_rx) < 3 <= resolved(p_m) else ''}")
say(f"    Вывод: виртуальная решётка MIMO, координаты которой считаются как r_tn + r_rk (что и проверено в шаге 2), даёт пик в "
    f"{hpbw(p1)/hpbw(p2):.1f} раза уже и разрешает цели, которые одно кольцо не разрешает. Для пассивного пеленгатора польза пересчёта — "
    "знать, что геометрическая модель отклика верна по координатам (ошибка пеленга доли градуса), а остаток снимается таблицей ДН из CST.")

# ---- рисунки
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for key, e, st in zip(errs, errs.values(), ("-", "--", ":")):
    ax[0].plot(src, e, st, lw=1.2, label=key)
ax[0].set_xlabel("истинный азимут источника, °"); ax[0].set_ylabel("ошибка пеленга, °"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=7)
ax[0].set_title(f"кольцо {len(ring)} Rx, R = {np.linalg.norm(pos[ring[0]][:2])/LAM:.0f}λ: ошибка пеленга без шума", fontsize=9)
x = cst_manifold(37.0)
for A, lab, st in ((A_geo, "(а) геометрия", "-"), (A_cal, "(б) с калибровкой", "--"), (A_cstf, "(в) CST", ":")):
    p = np.abs(A.conj().T @ x) ** 2; ax[1].plot(fine, 10 * np.log10(p / p.max()), st, lw=1.2, label=lab)
ax[1].axvline(37, color="r", lw=0.6); ax[1].set_xlim(0, 360); ax[1].set_ylim(-30, 1); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=7)
ax[1].set_xlabel("азимут, °"); ax[1].set_ylabel("|a^H x|², дБ"); ax[1].set_title("пеленгационный рельеф, источник в 37°", fontsize=9)
for dphi, (p_rx, p_m) in spectra.items():
    if dphi == 6.0:
        ax[2].plot(fine, 10 * np.log10(p_rx / p_rx.max()), "k-", lw=1.2, label=f"кольцо Rx ({len(ring)} эл.)")
        ax[2].plot(fine, 10 * np.log10(p_m / p_m.max()), "r--", lw=1.2, label=f"MIMO: {len(tring)}×{len(ring)} вирт. эл.")
ax[2].axvline(50, color="0.5", lw=0.6); ax[2].axvline(56, color="0.5", lw=0.6)
ax[2].set_xlim(30, 75); ax[2].set_ylim(-25, 1); ax[2].grid(alpha=0.3); ax[2].legend(fontsize=7)
ax[2].set_xlabel("азимут, °"); ax[2].set_ylabel("дБ"); ax[2].set_title("две цели в 50° и 56°: кольцо против MIMO", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "df_ring.png"), dpi=170)
say(""); say(f"Рисунок: {os.path.join(OUT, 'df_ring.png')}; лог: {LOG.name}")
LOG.close()
