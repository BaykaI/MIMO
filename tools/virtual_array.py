"""Виртуальная решётка MIMO и её двумерная ДН (этап 0 плана ДЗ2).

Исправленная и расширенная версия dz1/DZ_py.py:
  * координаты виртуальных элементов = попарная сумма координат Tx и Rx (как в [1], [3]);
  * ДН считается в полном 2D (theta, phi) по exp(jk(x sinθ cosφ + y sinθ sinφ)),
    без ошибки квадранта arctan(y/x), которая в DZ_py.py искажала кольцевую решётку;
  * численно оцениваются ширина луча по -3 дБ и УБЛ.

Формат файла данных — как в dz1/*.txt: колонки Rx_x Rx_y Tx_x Tx_y [0], строки с Nan
в колонках Tx пропускаются. Координаты в длинах волн.

Запуск:  python tools/virtual_array.py dz1/kolso.txt [--phi 0] [--plot]
"""
import argparse
import numpy as np


def load_layout(path):
    d = np.loadtxt(path)
    rx = d[:, :2]
    tx = d[:, 2:4]
    tx = tx[~np.isnan(tx[:, 0])]
    return rx, tx


def virtual_array(rx, tx, r0=(0.0, 0.0)):
    """Координаты виртуальных элементов: r_tx + r_rx + r0, форма (Ntx*Nrx, 2)."""
    return (tx[:, None, :] + rx[None, :, :] + np.asarray(r0)).reshape(-1, 2)


def pattern_cut(coords, theta_deg, phi_deg=0.0, weights=None):
    """Нормированная ДН |F|, дБ, в плоскости phi для набора изотропных элементов (λ = 1)."""
    k = 2 * np.pi
    th = np.radians(theta_deg)
    ph = np.radians(phi_deg)
    u = np.sin(th) * np.cos(ph)
    v = np.sin(th) * np.sin(ph)
    w = np.ones(len(coords)) if weights is None else np.asarray(weights)
    F = (w[None, :] * np.exp(1j * k * (np.outer(u, coords[:, 0]) + np.outer(v, coords[:, 1])))).sum(1)
    F = np.abs(F)
    return 20 * np.log10(F / F.max())


def pattern_2d(coords, theta_deg, phi_deg):
    """ДН на сетке (theta, phi), дБ; форма (len(theta), len(phi))."""
    out = np.empty((len(theta_deg), len(phi_deg)))
    for j, ph in enumerate(phi_deg):
        out[:, j] = pattern_cut(coords, theta_deg, ph)
    return out


def beam_metrics(theta_deg, FdB):
    """Ширина по -3 дБ (град) и УБЛ (дБ) для сечения с максимумом в theta=0."""
    i0 = int(np.argmin(np.abs(theta_deg)))
    lo = hi = i0
    while lo > 0 and FdB[lo - 1] >= -3:
        lo -= 1
    while hi < len(FdB) - 1 and FdB[hi + 1] >= -3:
        hi += 1
    bw = theta_deg[hi] - theta_deg[lo]
    # первый минимум справа от главного лепестка -> всё правее него считаем боковыми
    j = i0
    while j < len(FdB) - 1 and FdB[j + 1] <= FdB[j]:
        j += 1
    sll = FdB[j:].max() if j < len(FdB) - 1 else np.nan
    return bw, sll


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("layout")
    ap.add_argument("--phi", type=float, default=0.0, help="плоскость сечения, град")
    ap.add_argument("--n", type=int, default=16001)
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()

    rx, tx = load_layout(a.layout)
    va = virtual_array(rx, tx)
    theta = np.linspace(-80, 80, a.n)
    FdB = pattern_cut(va, theta, a.phi)
    bw, sll = beam_metrics(theta, FdB)
    print(f"Rx: {len(rx)}, Tx: {len(tx)}, виртуальных элементов: {len(va)}")
    print(f"phi = {a.phi:.1f}°: ширина -3 дБ = {bw:.3f}°, УБЛ = {sll:.2f} дБ")

    if a.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].scatter(rx[:, 0], rx[:, 1], marker="o", c="r", label="Rx")
        ax[0].scatter(tx[:, 0], tx[:, 1], marker="s", c="b", label="Tx")
        ax[0].scatter(va[:, 0], va[:, 1], marker="^", c="g", s=12, label="виртуальные")
        ax[0].set_aspect("equal"); ax[0].legend(); ax[0].set_xlabel("x, λ"); ax[0].set_ylabel("y, λ")
        ax[1].plot(theta, FdB); ax[1].axhline(-3, c="r", lw=0.8)
        ax[1].set_ylim(-60, 1); ax[1].set_xlabel("θ, град"); ax[1].set_ylabel("дБ"); ax[1].grid()
        plt.tight_layout(); plt.show()


if __name__ == "__main__":
    main()
