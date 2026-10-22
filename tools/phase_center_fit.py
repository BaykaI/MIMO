"""Оценка фазового центра по комплексной ДН (этапы 2–4 плана ДЗ2).

Модель: в дальней зоне фаза поля, отнесённая к началу координат,
    phase(θ, φ) ≈ φ0 + k·(x·sinθ·cosφ + y·sinθ·sinφ + z·cosθ),
где (x, y, z) — фазовый центр. Решается линейный МНК по (x, y, z, φ0)
в заданном угловом секторе, с весами |E|², после unwrap по θ.

Тот же фиттер применяется к:
  * ДН одного элемента            -> фазовый центр элемента (сверка с кнопкой CST);
  * ДН подрешётки (Combine Results) -> фазовый центр подрешётки;
  * произведению F_tx·F_rx (считается здесь) -> положение виртуального элемента
    пары / всей MIMO-решётки (гипотеза 1: r_pc = r_pc(Tx) + r_pc(Rx)).

Экспорт CST: Farfield -> Export -> Plot data (ASCII). Ожидаются колонки
  Theta[deg] Phi[deg] ... Abs(Theta)[dBi] Phase(Theta)[deg] Abs(Phi)[dBi] Phase(Phi)[deg] ...
Парсер ищет колонки по подстрокам в заголовке; проверьте вывод `--info`
на первом файле — порядок колонок зависит от версии CST.

Запуск:
  python tools/phase_center_fit.py --selftest
  python tools/phase_center_fit.py farfield_port1.txt --sector 30 --freq-ghz 10
  python tools/phase_center_fit.py tx_sum.txt --rx rx_sum.txt --sector 30 --freq-ghz 10   # произведение
"""
import argparse
import re
import numpy as np

C0 = 299_792_458.0


# ---------------------------------------------------------------- чтение CST ASCII
def read_cst_farfield(path, component="theta"):
    """Возвращает theta[deg], phi[deg], E (комплексная, линейная амплитуда)."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    header = next(l for l in lines if re.search(r"theta", l, re.I) and re.search(r"phi", l, re.I))
    cols = re.split(r"\s{2,}|\t", header.strip())
    cols = [c for c in cols if c]

    def find(*subs):
        for i, c in enumerate(cols):
            cl = c.lower().replace(" ", "")
            if all(s in cl for s in subs):
                return i
        raise KeyError(f"колонка {subs} не найдена в заголовке: {cols}")

    data = np.array([list(map(float, l.split())) for l in lines
                     if re.match(r"^\s*[-+]?\d", l)])
    i_th, i_ph = find("theta["), find("phi[")
    comp = component.lower()
    i_abs = find("abs(" + comp)
    i_pha = find("phase(" + comp)
    th, ph = data[:, i_th], data[:, i_ph]
    amp_db = data[:, i_abs]
    amp = 10 ** (amp_db / 20)
    E = amp * np.exp(1j * np.radians(data[:, i_pha]))
    return th, ph, E


def info(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        for i, l in enumerate(f):
            print(l.rstrip())
            if i > 4:
                break


# ---------------------------------------------------------------- фиттер
def direction_cosines(theta_deg, phi_deg):
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    return np.column_stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])


def initial_guess(th, ph, u, E, w, k):
    """Начальное приближение: фаза разворачивается (unwrap) вдоль θ на каждой φ-линии,
    затем линейный МНК по (x, y, z) с отдельной константой на каждую линию.
    Требование: шаг по θ такой, что k·|r|·Δθ < π (на 1° и λ = 3 см это |r| < 0,86 м)."""
    phase = np.angle(E)
    lines = []
    for p in np.unique(ph):
        idx = np.where(ph == p)[0]
        idx = idx[np.argsort(th[idx])]
        if len(idx) < 3:
            continue
        lines.append((idx, np.unwrap(phase[idx])))
    n = sum(len(i) for i, _ in lines)
    A = np.zeros((n, 3 + len(lines)))
    b = np.zeros(n)
    ww = np.zeros(n)
    r = 0
    for j, (idx, up) in enumerate(lines):
        s = slice(r, r + len(idx))
        A[s, :3] = k * u[idx]
        A[s, 3 + j] = 1.0
        b[s] = up
        ww[s] = w[idx]
        r += len(idx)
    sw = np.sqrt(ww)[:, None]
    x, *_ = np.linalg.lstsq(A * sw, b * sw[:, 0], rcond=None)
    return x[:3]


def fit_phase_center(theta_deg, phi_deg, E, k, sector_deg=30.0, center=(0.0, 0.0), fix_z=False,
                     with_phase=False):
    """МНК по фазе в секторе радиусом sector_deg вокруг center=(θc, φc).

    Начальное приближение — по развёрнутой вдоль θ фазе (initial_guess), затем
    итерационное уточнение линейным МНК по свёрнутому в (-π, π] остатку.
    fix_z=True — искать только x, y (z = 0): в узком секторе вокруг нормали z не обусловлен
    (чувствительность ∝ 1 − cosθ), а для плоской решётки нужны именно x, y.
    Возвращает (xyz в метрах, rms остатка фазы в градусах); with_phase=True — ещё и φ0 в градусах
    (фаза элемента в начале координат; разнофазность элементов смещает фазовый центр их суммы).
    """
    th, ph = np.asarray(theta_deg, float), np.asarray(phi_deg, float)
    E = np.asarray(E, complex)
    thc, phc = center
    u = direction_cosines(th, ph)
    uc = direction_cosines(np.array([thc]), np.array([phc]))[0]
    ang = np.degrees(np.arccos(np.clip(u @ uc, -1, 1)))
    m = (ang <= sector_deg) & (np.abs(E) > 1e-12)
    if m.sum() < 8:
        raise ValueError("слишком мало точек в секторе")
    phase = np.angle(E[m])
    A = np.column_stack([k * u[m], np.ones(m.sum())])
    w = np.abs(E[m]) ** 2
    x = np.zeros(4)
    x[:3] = initial_guess(th[m], ph[m], u[m], E[m], w, k)
    if fix_z:
        x[2] = 0.0
        A[:, 2] = 0.0
    for _ in range(50):
        resid = np.angle(np.exp(1j * (phase - A @ x)))
        dx, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], resid * np.sqrt(w), rcond=None)
        x += dx
        if np.abs(dx[:3]).max() * k < 1e-6:
            break
    resid = np.angle(np.exp(1j * (phase - A @ x)))
    rms = np.degrees(np.sqrt(np.average(resid ** 2, weights=w)))
    if with_phase:
        return x[:3], rms, float(np.degrees(np.angle(np.exp(1j * x[3]))))
    return x[:3], rms


def product_pattern(th1, ph1, E1, th2, ph2, E2):
    """Поэлементное произведение двух ДН на одинаковой сетке."""
    if not (np.allclose(th1, th2) and np.allclose(ph1, ph2)):
        raise ValueError("сетки (theta, phi) двух ДН не совпадают")
    return th1, ph1, E1 * E2


# ---------------------------------------------------------------- самопроверка
def synth_point_source(xyz, k, step=1.0):
    """Изотропный точечный источник в xyz: E = exp(+jk r·u) (фаза относительно начала координат)."""
    th = np.arange(0, 180 + 1e-9, step)
    ph = np.arange(0, 360, step)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    th, ph = TH.ravel(), PH.ravel()
    u = direction_cosines(th, ph)
    E = np.exp(1j * k * (u @ np.asarray(xyz)))
    return th, ph, E


def selftest():
    lam = 0.03
    k = 2 * np.pi / lam
    rng = np.random.default_rng(0)
    ok = True
    for r in ([0.004, -0.002, 0.0], [0.11, 0.0, 0.001], [-0.1275, 0.03, 0.0]):
        th, ph, E = synth_point_source(r, k, step=1.0)
        E *= np.exp(1j * rng.normal(0, np.radians(3), E.shape))  # 3° фазового шума
        est, rms = fit_phase_center(th, ph, E, k, sector_deg=30, center=(0, 0))
        err = np.abs(est - r).max()
        print(f"ист. {np.round(r,4)}  оценка {np.round(est,4)}  макс.ошибка {err*1e3:.3f} мм  rms {rms:.2f}°")
        ok &= err < lam / 200
    # гипотеза 1 на точечных источниках: произведение ДН -> сумма координат
    rt, rr = [-0.1275, 0.03, 0], [0.105, 0, 0]
    th, ph, Et = synth_point_source(rt, k)
    _, _, Er = synth_point_source(rr, k)
    est, rms = fit_phase_center(th, ph, Et * Er, k, sector_deg=30)
    print(f"произведение: ожидание {np.add(rt, rr)}  оценка {np.round(est,4)}  rms {rms:.2f}°")
    ok &= np.abs(est - np.add(rt, rr)).max() < lam / 200
    print("SELFTEST", "OK" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?")
    ap.add_argument("--rx", help="вторая ДН: будет взято произведение file·rx")
    ap.add_argument("--freq-ghz", type=float)
    ap.add_argument("--sector", type=float, default=30.0)
    ap.add_argument("--center", type=float, nargs=2, default=(0.0, 0.0), metavar=("THETA", "PHI"))
    ap.add_argument("--component", default="theta", help="theta | phi — поляризационная компонента")
    ap.add_argument("--info", action="store_true", help="показать заголовок файла")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        raise SystemExit(0 if selftest() else 1)
    if a.info:
        info(a.file); return
    if a.freq_ghz is None:
        ap.error("--freq-ghz обязателен")
    k = 2 * np.pi * a.freq_ghz * 1e9 / C0
    th, ph, E = read_cst_farfield(a.file, a.component)
    if a.rx:
        th2, ph2, E2 = read_cst_farfield(a.rx, a.component)
        th, ph, E = product_pattern(th, ph, E, th2, ph2, E2)
    xyz, rms = fit_phase_center(th, ph, E, k, a.sector, tuple(a.center))
    print(f"фазовый центр, мм: x={xyz[0]*1e3:.3f} y={xyz[1]*1e3:.3f} z={xyz[2]*1e3:.3f};  rms остатка {rms:.2f}°")


if __name__ == "__main__":
    main()
