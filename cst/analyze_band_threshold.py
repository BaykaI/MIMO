"""Порог деградации по ширине полосы при узкополосной обработке (фазы пар по 10 ГГц), по формулам
(точечные излучатели в CST совпали с формулами). Полоса 0-100 % через 2 %, частоты через 0,1 ГГц.
  python cst/analyze_band_threshold.py -> cst/degr_thr.json
Для каждого шага: потеря пика и УБЛ по азимуту от полосы, пороги потери 1 и 3 дБ, ΔL (разброс пути пар
вдоль направления на цель) и произведение порог·ΔL/λ.
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 299_792_458.0
kf = lambda f: 2 * np.pi * f * 1e9 / C0 / 1e3
F0, PHS = 10.0, 17.0
AZ = np.arange(0, 360, 0.5); U = np.column_stack([np.cos(np.radians(AZ)), np.sin(np.radians(AZ)), np.zeros(AZ.size)])
us = np.array([np.cos(np.radians(PHS)), np.sin(np.radians(PHS)), 0.0])
DR = np.arange(-300, 300.1, 2.0)
RELS = [round(0.02 * i, 2) for i in range(51)]
SPACINGS = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0]


def geometry(d):
    R = 12 * d * 30.0 / (2 * np.pi)
    rx = np.array([[R * np.cos(2 * np.pi * i / 12), R * np.sin(2 * np.pi * i / 12), 0] for i in range(12)])
    tx = np.array([[1.35 * R * np.cos(2 * np.pi * (j + .5) / 6), 1.35 * R * np.sin(2 * np.pi * (j + .5) / 6), 0] for j in range(6)])
    return (tx[:, None, :] + rx[None, :, :]).reshape(-1, 3)


def sll_circ(p):
    d = 20 * np.log10(p / p.max() + 1e-12); n = len(d); i = int(np.argmax(d)); l = r = 0
    while l > -(n // 2) and d[(i + l - 1) % n] < d[(i + l) % n]: l -= 1
    while r < n // 2 and d[(i + r + 1) % n] < d[(i + r) % n]: r += 1
    main = {(i + j) % n for j in range(l, r + 1)}
    return float(max(d[j] for j in range(n) if j not in main))


out = {}
for d in SPACINGS:
    RV = geometry(d); A0 = np.exp(1j * kf(F0) * (RV @ U.T)); rows = []
    for rel in RELS:
        nf = max(1, int(round(rel * F0 / 0.1)) + 1)
        fr = np.linspace(F0 * (1 - rel / 2), F0 * (1 + rel / 2), nf) if nf > 1 else np.array([F0])
        Iw = np.zeros((len(DR), len(AZ)), complex); In = np.zeros_like(Iw)
        for f in fr:
            x = np.exp(1j * kf(f) * (RV @ us)); rng = np.exp(1j * 2 * kf(f) * DR)[:, None]
            Iw += rng * ((np.exp(1j * kf(f) * (RV @ U.T)).conj().T @ x) / len(RV))[None, :]
            In += rng * ((A0.conj().T @ x) / len(RV))[None, :]
        aw, an = np.abs(Iw), np.abs(In)
        i, j = np.unravel_index(np.argmax(an), an.shape); iw, jw = np.unravel_index(np.argmax(aw), aw.shape)
        rows.append(dict(rel=rel, loss_db=float(20 * np.log10(an.max() / aw.max())), sll_nb=sll_circ(an[i, :]), sll_wb=sll_circ(aw[iw, :])))
    loss = np.array([r["loss_db"] for r in rows]); rr = np.array(RELS)
    thr = lambda L: float(np.interp(L, -loss, rr)) if (-loss).max() >= L else float("nan")
    dL = float(np.ptp(RV @ us)); lam = C0 / (F0 * 1e9) * 1e3
    out[f"{d:g}"] = dict(d=d, dL_mm=dL, dL_lam=dL / lam, thr1=thr(1.0), thr3=thr(3.0), rows=rows,
                         k1=thr(1.0) * dL / lam, k3=thr(3.0) * dL / lam)
    print(f"шаг {d}: ΔL {dL:.0f} мм = {dL/lam:.1f}λ; порог 1 дБ {thr(1.0)*100:.0f} %, 3 дБ {thr(3.0)*100:.0f} %; "
          f"порог·ΔL/λ: {thr(1.0)*dL/lam:.2f} и {thr(3.0)*dL/lam:.2f}", flush=True)
json.dump(out, open(os.path.join(HERE, "degr_thr.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
