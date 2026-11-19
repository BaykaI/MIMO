"""Рисунки по формулам Чапурского для отчёта о полосе и шаге (оформление gost_plot).
  python cst/formula_figs.py   (нужен formula_band_spacing.json от formula_band_spacing.py)
Выход: fig3/f_*.png
"""
import os, json
import numpy as np
import gost_plot as G

C0 = 299_792_458.0
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "fig3"); os.makedirs(OUT, exist_ok=True)
J = json.load(open(os.path.join(HERE, "formula_band_spacing.json"), encoding="utf-8"))
S = J["spacing"]; d = np.array([r["d"] for r in S])
LAM = C0 / 10e9 * 1e3

# 1. ширина луча от шага
fig, ax = G.figure(h_mm=70); a = ax[0, 0]
for i, (k, lab, at, tx, ty) in enumerate((("rx_1f_bw", 1, 0.8, 1.0, 25), ("mimo_1f_bw", 2, 1.2, 1.0, 3.5))):
    y = [r[k] for r in S]; a.plot(d, y, **G.style(i)); 
G.key(a, a.get_lines()); a.set_yscale("log"); G.no_minor_labels(a); a.set_xlabel("Шаг приёмников по дуге d/λ"); a.set_ylabel("Ширина луча по азимуту, град.")
a.set_yticks([2, 3, 5, 10, 20, 50]); a.set_yticklabels(["2", "3", "5", "10", "20", "50"]); G.comma(a, "x")
G.save(fig, os.path.join(OUT, "f_spacing_bw.png"))

# 2. ложные максимумы от шага
fig, ax = G.figure(h_mm=80); a = ax[0, 0]
spec = (("rx_1f_sll", 0, 1.15, 1.32, 1.5), ("rx_band_sll", 1, 2.1, 2.35, -12.6), ("mimo_1f_sll", 2, 0.95, 0.75, -20.0), ("mimo_band_sll", 3, 2.5, 2.6, -19.5))
for n, (k, st, at, tx, ty) in enumerate(spec, start=1):
    y = [r[k] for r in S]; a.plot(d, y, **G.style(st)); 
G.key(a, a.get_lines()[:4]); a.axhline(-6, color="k", lw=0.6, ls=(0, (1, 3)))
a.set_xlabel("Шаг приёмников по дуге d/λ"); a.set_ylabel("Наибольший ложный максимум, дБ"); G.comma(a)
G.save(fig, os.path.join(OUT, "f_spacing_sll.png"))

# 3. сечение по азимуту MIMO d = 1,5λ: одна частота и полоса; кольцо ДЗ1 1/4/16 частот
NRX, NTX = 12, 6


def ring(dl):
    R = NRX * dl * LAM / (2 * np.pi)
    a_ = 2 * np.pi * np.arange(NRX) / NRX; b_ = 2 * np.pi * (np.arange(NTX) + 0.5) / NTX
    rx = np.column_stack([R * np.cos(a_), R * np.sin(a_)]); tx = np.column_stack([1.35 * R * np.cos(b_), 1.35 * R * np.sin(b_)])
    return (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2)


def af(pts, freqs, az):
    u = np.column_stack([np.cos(np.radians(az)), np.sin(np.radians(az))]); u0 = np.array([1.0, 0.0]); s = 0
    for f in freqs:
        k = 2 * np.pi * f / C0 / 1e3; s = s + np.exp(1j * k * ((u - u0) @ pts.T)).sum(axis=1)
    s = np.abs(s); return 20 * np.log10(s / s.max() + 1e-9)


fig, ax = G.figure(2, 1, h_mm=62)
v = ring(1.5); az = np.arange(-180, 180, 0.05)
for n, (fr, st) in enumerate(((np.array([10e9]), 0), (np.arange(8.5, 11.51, 0.25) * 1e9, 2)), start=1):
    y = af(v, fr, az); ax[0, 0].plot(az, y, **G.style(st))
    
G.key(ax[0, 0], ax[0, 0].get_lines()); ax[0, 0].set_xlim(-180, 180); ax[0, 0].set_xticks([-180, -90, 0, 90, 180]); ax[0, 0].set_ylim(-30, 1); ax[0, 0].set_xlabel("Азимут, град."); ax[0, 0].set_ylabel("Ψ, дБ"); G.comma(ax[0, 0])
G.letter(ax[0, 0], "а)")
a_ = 2 * np.pi * np.arange(20) / 20; b_ = 2 * np.pi * np.arange(10) / 10
rx = np.column_stack([50 * LAM * np.cos(a_), 50 * LAM * np.sin(a_)]); tx = np.column_stack([67.5 * LAM * np.cos(b_), 67.5 * LAM * np.sin(b_)])
v50 = (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2); az2 = np.arange(-10, 10, 0.005)
win = {}
for n, (M, st) in enumerate(((1, 2), (16, 0)), start=1):
    fr = np.linspace(10e9, 16e9, M) if M > 1 else np.array([10e9])
    y = af(v50, fr, az2); ax[0, 1].plot(az2, y, color=G.STYLES[st]["color"], ls=G.STYLES[st]["ls"], lw=0.9)
    pk = np.abs(az2) > 0.25                     # вне главного лепестка (ширина 0,24°)
    win[f"M{M}"] = float(y[pk].max())
import json as _j
_j.dump(dict(ring50_window10=win), open(os.path.join(HERE, "f_band_az.json"), "w"))

G.key(ax[0, 1], ax[0, 1].get_lines()); ax[0, 1].set_xlim(-10, 10); ax[0, 1].set_ylim(-30, 1); ax[0, 1].set_xlabel("Азимут, град."); ax[0, 1].set_ylabel("Ψ, дБ"); G.comma(ax[0, 1])
G.letter(ax[0, 1], "б)")
G.save(fig, os.path.join(OUT, "f_band_az.png"), letters=True)
print("ok", sorted(os.listdir(OUT)))
