"""Расчёт по Чапурскому [7]: бистатические задержки (4)–(6), корреляционный интеграл (9), ОФН (10)–(12).

Ψ(r, r0) = (1/(Nt·Nr·M)) |Σ_n Σ_k Σ_m exp{j ξ_m ω0 [τ_nk(r) − τ_nk(r0)]}|,   τ_nk(r) = (|r − r_tn| + |r − r_rk|)/c.
Одна частота (M = 1) — формула (12). Всё считается честно по задержкам, без приближения дальней зоны.

Режимы:
  python cst/chap_ofn.py chap      — воспроизведение его примера (с. 83–84): 9×9 Tx шаг 0,25 м + 9×9 Rx шаг 0,288 м (81 + 81),
                                     f0 = 6 ГГц, точка на 1 м по нормали; сечение ОФН в плоскости (x, y) 1×1 м;
                                     у него: разрешение «единицы сантиметров», УБЛ ≤ 0,06 (−24 дБ) в области 2×2 м.
  python cst/chap_ofn.py cone      — конусная MIMO-решётка (геометрия из cone_geometry()), ОФН в угловых координатах
                                     вокруг опорного направления (горизонт, φ = 0) на дальности R0, и карта в плоскости.
Выход: fig2/chap_*.png, cst/cone_chap.json (числа для отчёта).
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C0 = 299_792_458.0
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig2"); os.makedirs(OUT, exist_ok=True)
NL = chr(10)


def tau(pts, tx, rx):
    """Бистатические задержки (4)–(6): pts (P,3), tx (Nt,3), rx (Nr,3) -> (P, Nt, Nr), с."""
    dt = np.linalg.norm(pts[:, None, :] - tx[None, :, :], axis=2)
    dr = np.linalg.norm(pts[:, None, :] - rx[None, :, :], axis=2)
    return (dt[:, :, None] + dr[:, None, :]) / C0


def ofn(pts, r0, tx, rx, freqs, weights=None, block=20000):
    """ОФН (10)–(11): нормированный модуль корреляционного интеграла (9) для точек pts относительно опорной r0.
    Считается блоками точек и по одному передатчику: память ≤ block × Nr комплексных чисел (≈ 25 МБ при Nr = 81)."""
    w = np.ones(len(freqs)) if weights is None else np.asarray(weights)
    pts = np.asarray(pts, float); r0 = np.asarray(r0, float)
    t0 = tau(r0[None, :], tx, rx)[0]                      # (Nt, Nr)
    q = np.zeros(len(pts), complex)
    for b in range(0, len(pts), block):
        p = pts[b:b + block]
        dr = np.linalg.norm(p[:, None, :] - rx[None, :, :], axis=2) / C0      # (P, Nr)
        for n in range(tx.shape[0]):
            dt = np.linalg.norm(p - tx[n], axis=1) / C0                       # (P,)
            dtau = dt[:, None] + dr - t0[n][None, :]                           # (P, Nr)
            for f, wm in zip(freqs, w):
                q[b:b + block] += wm * np.exp(1j * 2 * np.pi * f * dtau).sum(axis=1)
    return np.abs(q) / (tx.shape[0] * rx.shape[0] * w.sum())


def cone_geometry(lam=0.03):
    """Конусная MIMO-решётка (усечённый конус, ось z, основание в z = 0):
    четыре кольца на боковой поверхности (радиус 5λ…2,5λ, высота 0…4,2λ, полуугол 30°): Rx — кольца 1 и 3 (10 + 8), Tx — кольца 2 и 4 (8 + 6)."""
    rings = [(5.0 * lam, 0.0, 10, "rx"), (4.17 * lam, 1.4 * lam, 8, "tx"),
             (3.33 * lam, 2.8 * lam, 8, "rx"), (2.5 * lam, 4.2 * lam, 6, "tx")]
    tx, rx, names = [], [], []
    for R, z, n, kind in rings:
        for i in range(n):
            a = 2 * np.pi * i / n + (np.pi / n if kind == "tx" else 0.0)
            p = [R * np.cos(a), R * np.sin(a), z]
            (tx if kind == "tx" else rx).append(p)
    return np.array(tx), np.array(rx)


def metrics_cut(ang, psi):
    d = 20 * np.log10(psi / psi.max() + 1e-12)
    im = int(np.argmax(d)); lo = hi = im
    while lo > 0 and d[lo - 1] >= -3: lo -= 1
    while hi < len(d) - 1 and d[hi + 1] >= -3: hi += 1
    # границы −3 дБ с линейной интерполяцией между соседними отсчётами
    al = ang[lo] + (ang[lo - 1] - ang[lo]) * (-3 - d[lo]) / (d[lo - 1] - d[lo]) if lo > 0 and d[lo - 1] != d[lo] else ang[lo]
    ar = ang[hi] + (ang[hi + 1] - ang[hi]) * (-3 - d[hi]) / (d[hi + 1] - d[hi]) if hi < len(d) - 1 and d[hi + 1] != d[hi] else ang[hi]
    l = r = im
    while l > 0 and d[l - 1] < d[l]: l -= 1
    while r < len(d) - 1 and d[r + 1] < d[r]: r += 1
    outside = np.concatenate([d[:l], d[r + 1:]])
    return dict(hpbw=float(ar - al), sll_db=float(outside.max()) if outside.size else float("nan"))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "chap"
    if mode == "chap":
        # ---- пример Чапурского, с. 83–84 ----
        f0 = 6.0e9
        g = (np.arange(9) - 4.0)          # Nt + 1 = Nr + 1 = 81 элемент ([7], с. 83): 9×9
        xt, yt = np.meshgrid(g * 0.25, g * 0.25); xr, yr = np.meshgrid(g * 0.288, g * 0.288)
        tx = np.column_stack([xt.ravel(), yt.ravel(), np.zeros(81)]); rx = np.column_stack([xr.ravel(), yr.ravel(), np.zeros(81)])
        r0 = np.array([0.0, 0.0, 1.0])
        # сечение в плоскости (x, y) на z = 1 м: область 1×1 м (его рис. 3) и 2×2 м (его оценка УБЛ)
        res = {}
        for span, n in ((0.5, 401), (1.0, 601)):
            x = np.linspace(-span, span, n); X, Y = np.meshgrid(x, x)
            pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, 1.0)])
            psi = ofn(pts, r0, tx, rx, [f0]).reshape(X.shape)
            res[span] = (x, psi)
        x, psi = res[0.5]
        cut = psi[psi.shape[0] // 2]
        m = metrics_cut(x, cut)
        # УБЛ в области 2×2 м без главного лепестка (радиус 3 см)
        x2, psi2 = res[1.0]; X2, Y2 = np.meshgrid(x2, x2)
        sl = psi2[(X2 ** 2 + Y2 ** 2) > 0.03 ** 2].max()
        print(f"Чапурский, 9×9+9×9, 6 ГГц, точка на 1 м: ширина по −3 дБ {m['hpbw']*100:.1f} см; "
              f"макс. боковой лепесток в 2×2 м: {sl:.3f} = {20*np.log10(sl):.1f} дБ (у него ≤ 0,06 = −24 дБ)")
        fig, ax = plt.subplots(1, 3, figsize=(13, 4))
        ax[0].plot(tx[:, 0], tx[:, 1], "o", mfc="none", mec="0.3", ms=5, label="Tx 9×9, шаг 0,25 м")
        ax[0].plot(rx[:, 0], rx[:, 1], "s", mfc="none", mec="tab:blue", ms=5, label="Rx 9×9, шаг 0,288 м")
        ax[0].set_aspect("equal"); ax[0].legend(fontsize=7); ax[0].set_title("решётка из [7], рис. 2", fontsize=9); ax[0].set_xlabel("x, м"); ax[0].set_ylabel("y, м")
        im = ax[1].imshow(psi, extent=[-0.5, 0.5, -0.5, 0.5], origin="lower", cmap="hot", vmin=0, vmax=1)
        ax[1].set_title("ОФН Ψ(x, y) на z = 1 м по ф. (12) — как рис. 3, а в [7]", fontsize=9); ax[1].set_xlabel("x, м"); ax[1].set_ylabel("y, м")
        fig.colorbar(im, ax=ax[1], shrink=0.8)
        ax[2].plot(x * 100, 20 * np.log10(cut + 1e-6), "k", lw=1)
        ax[2].axhline(-24, color="r", ls="--", lw=0.8); ax[2].text(-48, -23, "−24 дБ: оценка УБЛ в [7]", color="r", fontsize=8)
        ax[2].set_ylim(-50, 1); ax[2].set_xlabel("x, см (y = 0)"); ax[2].set_ylabel("Ψ, дБ"); ax[2].grid(alpha=0.3)
        ax[2].set_title(f"сечение: ширина {m['hpbw']*100:.1f} см, УБЛ в 2×2 м {20*np.log10(sl):.1f} дБ", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "chap_example.png"), dpi=170)
        json.dump(dict(hpbw_cm=m["hpbw"] * 100, sll_db=float(20 * np.log10(sl))), open(os.path.join(os.path.dirname(OUT), "chap_example.json"), "w"))
    else:
        # ---- конус по тому же методу ----
        lam = 0.03; f0 = C0 / lam
        tx, rx = cone_geometry(lam)
        np.savetxt(os.path.join(os.path.dirname(OUT), "cone_layout_m.txt"),
                   np.vstack([np.column_stack([np.zeros(len(tx)), tx]), np.column_stack([np.ones(len(rx)), rx])]),
                   header="kind(0=tx,1=rx) x y z [м]", fmt="%.6f")
        out = {}
        for R0 in (3.0, 1000.0):
            r0 = np.array([R0, 0.0, 0.0])                     # опорная точка: горизонт, φ = 0, дальность R0
            az = np.radians(np.linspace(-90, 90, 1441)); el = np.radians(np.linspace(-60, 60, 961))
            pts_az = np.column_stack([R0 * np.cos(az), R0 * np.sin(az), np.zeros(az.size)])
            pts_el = np.column_stack([R0 * np.cos(el), np.zeros(el.size), R0 * np.sin(el)])
            psi_az = ofn(pts_az, r0, tx, rx, [f0]); psi_el = ofn(pts_el, r0, tx, rx, [f0])
            # 3 частоты как в CST (9,8/10/10,2 ГГц) — для сравнения формы
            psi_az3 = ofn(pts_az, r0, tx, rx, [f0 * 0.98, f0, f0 * 1.02])
            out[R0] = dict(az=np.degrees(az), el=np.degrees(el), psi_az=psi_az, psi_el=psi_el, psi_az3=psi_az3,
                           m_az=metrics_cut(np.degrees(az), psi_az), m_el=metrics_cut(np.degrees(el), psi_el))
            print(f"конус, R0 = {R0:g} м: азимут — ширина {out[R0]['m_az']['hpbw']:.2f}°, УБЛ {out[R0]['m_az']['sll_db']:.1f} дБ; "
                  f"угол места — ширина {out[R0]['m_el']['hpbw']:.2f}°, УБЛ {out[R0]['m_el']['sll_db']:.1f} дБ")
        fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
        a3 = fig.add_subplot(1, 3, 1, projection="3d"); ax[0].remove()
        a3.scatter(tx[:, 0] * 1e3, tx[:, 1] * 1e3, tx[:, 2] * 1e3, c="tab:red", marker="s", s=18, label=f"Tx ({len(tx)})")
        a3.scatter(rx[:, 0] * 1e3, rx[:, 1] * 1e3, rx[:, 2] * 1e3, c="tab:blue", marker="o", s=18, label=f"Rx ({len(rx)})")
        v = tx[:, None, :] + rx[None, :, :]; v = v.reshape(-1, 3) * 1e3
        a3.scatter(v[:, 0], v[:, 1], v[:, 2], c="0.6", marker="^", s=6, alpha=0.5, label=f"виртуальные центры r_tx + r_rx ({len(tx)*len(rx)})")
        a3.set_xlabel("x, мм"); a3.set_ylabel("y, мм"); a3.set_zlabel("z, мм"); a3.legend(fontsize=6, loc="upper left")
        a3.set_title("конусная MIMO-решётка и её виртуальный эквивалент", fontsize=9)
        for R0, c in ((1000.0, "k"), (3.0, "tab:orange")):
            ax[1].plot(out[R0]["az"], 20 * np.log10(out[R0]["psi_az"] + 1e-6), c, lw=1, label=f"R0 = {R0:g} м")
            ax[2].plot(out[R0]["el"], 20 * np.log10(out[R0]["psi_el"] + 1e-6), c, lw=1, label=f"R0 = {R0:g} м")
        ax[1].plot(out[1000.0]["az"], 20 * np.log10(out[1000.0]["psi_az3"] + 1e-6), "tab:blue", lw=0.8, ls="--", label="3 частоты 9,8/10/10,2 ГГц")
        for a_, t in ((ax[1], "азимут (горизонт), опорное направление φ = 0"), (ax[2], "угол места, плоскость φ = 0")):
            a_.set_ylim(-40, 1); a_.set_ylabel("Ψ, дБ"); a_.grid(alpha=0.3); a_.legend(fontsize=7); a_.set_title("ОФН по ф. (12): " + t, fontsize=9)
        ax[1].set_xlabel("азимут, °"); ax[2].set_xlabel("угол места, °")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "chap_cone.png"), dpi=170)
        json.dump({str(k): dict(m_az=v_["m_az"], m_el=v_["m_el"]) for k, v_ in out.items()},
                  open(os.path.join(os.path.dirname(OUT), "cone_chap.json"), "w", encoding="utf-8"), indent=1)
        np.savez(os.path.join(os.path.dirname(OUT), "cone_chap_cuts.npz"), **{f"{k}_{n}": v_[n] for k, v_ in out.items() for n in ("az", "el", "psi_az", "psi_el", "psi_az3")})
