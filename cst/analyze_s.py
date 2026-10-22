"""Разбор S-матрицы 12×12 из Touchstone (этапы 3–4): согласование и развязка Tx–Rx.

  python cst/analyze_s.py cst/stage3/mimo_half.s12p [--f 10]
  python cst/analyze_s.py cst/stage4/sweep_ytx*/mimo_patch.s12p --f 10     # кривая «развязка vs сдвиг ряда Tx»

Порты 1–8 = Rx1…Rx8, 9–12 = Tx1…Tx4. Выход: таблица в консоль + <файл>.md рядом с каждым .s12p.
"""
import sys, os, re, glob, argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("files", nargs="+")
ap.add_argument("--f", type=float, default=10.0, help="частота отчёта, ГГц")
a = ap.parse_args()
NAMES = [f"Rx{i}" for i in range(1, 9)] + [f"Tx{i}" for i in range(1, 5)]


def read_touchstone(path):
    """Возвращает f[ГГц], S[nf, n, n] (комплексная)."""
    n = int(re.search(r"\.s(\d+)p$", path, re.I).group(1))
    fmt, fscale = "MA", 1e9
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.split("!")[0].strip()
            if not s:
                continue
            if s.startswith("#"):
                t = s.upper().split()
                fscale = {"HZ": 1, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}[t[1]]
                fmt = t[3]
                continue
            rows += [float(v) for v in s.split()]
    per = 1 + 2 * n * n
    data = np.array(rows).reshape(-1, per)
    f = data[:, 0] * fscale / 1e9
    v = data[:, 1:].reshape(-1, n * n, 2)
    if fmt == "MA":
        S = v[..., 0] * np.exp(1j * np.radians(v[..., 1]))
    elif fmt == "DB":
        S = 10 ** (v[..., 0] / 20) * np.exp(1j * np.radians(v[..., 1]))
    else:
        S = v[..., 0] + 1j * v[..., 1]
    return f, S.reshape(-1, n, n)


summary = []
for path in sorted(sum((glob.glob(p) for p in a.files), [])):
    f, S = read_touchstone(path)
    i = int(np.argmin(np.abs(f - a.f)))
    Sdb = 20 * np.log10(np.abs(S[i]) + 1e-12)
    n = S.shape[1]
    lines = [f"# S-матрица {os.path.relpath(path)} на {f[i]:.2f} ГГц, дБ", "",
             "| | " + " | ".join(NAMES[:n]) + " |", "|---|" + "---|" * n]
    for r in range(n):
        lines.append(f"| {NAMES[r]} | " + " | ".join(f"{Sdb[r, c]:.1f}" for c in range(n)) + " |")
    # согласование: минимум |Sii| по частоте
    lines += ["", "| порт | min\\|Sii\\|, дБ | f мин, ГГц | \\|Sii\\| на f отчёта, дБ |", "|---|---|---|---|"]
    for r in range(n):
        d = 20 * np.log10(np.abs(S[:, r, r]) + 1e-12)
        j = int(np.argmin(d))
        lines.append(f"| {NAMES[r]} | {d[j]:.1f} | {f[j]:.3f} | {Sdb[r, r]:.1f} |")
    # развязка Tx–Rx: худшая по каждому Tx и в целом
    if n == 12:
        worst = []
        for t in range(8, 12):
            c = int(np.argmax(Sdb[t, :8]))
            worst.append((NAMES[t], NAMES[c], Sdb[t, c]))
        lines += ["", "Худшая развязка каждого Tx с Rx (дБ): " +
                  "; ".join(f"{t}–{r} {v:.1f}" for t, r, v in worst)]
        wmax = max(v for _, _, v in worst)
        rxrx = max(Sdb[r, c] for r in range(8) for c in range(8) if r != c)
        lines.append(f"Худшая развязка Tx–Rx: {wmax:.1f} дБ; худшая Rx–Rx: {rxrx:.1f} дБ")
        m = re.search(r"ytx(\d+(?:\.\d+)?)", path)
        summary.append((float(m.group(1)) if m else float("nan"), wmax, rxrx))
    txt = "\n".join(lines).replace(".", ",")
    print(txt + "\n")
    with open(os.path.splitext(path)[0] + "_S.md", "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")

if len(summary) > 1:
    print("| сдвиг ряда Tx, мм | худшая развязка Tx–Rx, дБ | худшая Rx–Rx, дБ |\n|---|---|---|")
    for y, w, r in sorted(summary):
        print(f"| {y:g} | {w:.1f} | {r:.1f} |".replace(".", ","))
