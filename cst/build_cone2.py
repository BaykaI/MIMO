"""Конусная MIMO-решётка (геометрия chap_ofn.cone_geometry, как в отчёте по конусу) для двух исследований:
широкополосный приём при неизменном шаге и шаг при неизменной полосе.

  python build_cone2.py <излучатель pt|dip> <масштаб s> <f1:f2:шаг ГГц> <папка> [шаг ДН, град.] [n1,n2,n3,n4]
    n1..n4 – число излучателей на кольцах (по умолчанию 10, 8, 8, 6); окружение: CONE_TH="60,120" – пояс выгрузки ДН,
    PT_GAP – длина дискретного порта точечного источника, мм.
    pt  – точечный источник: дискретный порт длиной 1 мм без металла (элементарный электрический диполь вдоль z);
    dip – полуволновый вибратор 12,5 мм, провод 0,4 мм, порт в зазоре 0,5 мм.
    s   – масштаб всех размеров конуса (радиусы колец и высоты); s = 1 – конус из отчёта (R1 = 5λ, шаг Rx 3,14λ).
Конус: ось z, четыре кольца; Rx – кольца 1 и 3 (10 + 8), Tx – кольца 2 и 4 (8 + 6), Tx смещены на полшага.
Порты: rx1…rx18, tx1…tx14. Фаза ДН отсчитывается от начала координат (центр основания).
Выход: cst/<папка>/ (проект, ff_<порт>_<f>GHz.txt, layout_mm.txt, freqs.txt, build.log).
"""
import sys, os, time, atexit, math
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface
from export_ff import export_all

KIND = sys.argv[1]; S = float(sys.argv[2]); FSPEC = sys.argv[3]; TAG = sys.argv[4]
STEP = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
COUNTS = [int(x) for x in sys.argv[6].split(",")] if len(sys.argv) > 6 else None   # число излучателей на кольцах (постоянный размер)
TH = tuple(map(float, os.environ.get("CONE_TH", "0,120").split(",")))
GRID = float(os.environ.get("CONE_GRID", "0"))      # шаг округления координат излучателей, мм (0 – без округления)
LAM = 30.0
_f1, _f2, _fs = map(float, FSPEC.split(":"))
FREQS = tuple(round(_f1 + _fs * i, 3) for i in range(int(round((_f2 - _f1) / _fs)) + 1))
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, TAG); os.makedirs(OUT, exist_ok=True)
PRJ = os.path.join(OUT, "cone.cst")

rings = [(5.0 * LAM, 0.0, 10, "rx"), (4.17 * LAM, 1.4 * LAM, 8, "tx"), (3.33 * LAM, 2.8 * LAM, 8, "rx"), (2.5 * LAM, 4.2 * LAM, 6, "tx")]
tx, rx = [], []
if COUNTS:
    rings = [(R, z, c, kd) for (R, z, n, kd), c in zip(rings, COUNTS)]
for R, z, n, kind in rings:
    for i in range(n):
        a = 2 * math.pi * i / n + (math.pi / n if kind == "tx" else 0.0)
        x_, y_ = S * R * math.cos(a), S * R * math.sin(a)
        if GRID:      # координаты на сетке GRID мм: порты ложатся на общие линии сетки CST (нет тонких ячеек)
            x_, y_ = round(x_ / GRID) * GRID, round(y_ / GRID) * GRID
        (tx if kind == "tx" else rx).append((x_, y_, S * z))
elems = [(f"rx{i+1}",) + p for i, p in enumerate(rx)] + [(f"tx{i+1}",) + p for i, p in enumerate(tx)]
with open(os.path.join(OUT, "layout_mm.txt"), "w", encoding="utf-8") as fh:
    fh.write("name x_mm y_mm z_mm\n" + "".join(f"{n} {x:.4f} {y:.4f} {z:.4f}\n" for n, x, y, z in elems))
open(os.path.join(OUT, "freqs.txt"), "w").write(" ".join(f"{f:g}" for f in FREQS))
d_rx = 2 * math.pi * 5.0 * LAM * S / rings[0][2]
print(f"{KIND}, s = {S}: шаг Rx по нижнему кольцу {d_rx:.1f} мм = {d_rx / LAM:.2f}λ, {len(FREQS)} частот {FREQS[0]}-{FREQS[-1]} ГГц", flush=True)

for f in os.listdir(OUT):
    if f.endswith(".lok"):
        os.remove(os.path.join(OUT, f))
de = cst.interface.DesignEnvironment.new()
de.set_quiet_mode(True)
atexit.register(lambda: de.close())
prj = de.new_mws()
m = prj.model3d
cnt = [0]


def h(name, vba):
    cnt[0] += 1
    m.add_to_history(f"{cnt[0]:03d} {name}", vba)


GAP = float(os.environ.get("PT_GAP", "1.0")) if KIND == "pt" else 0.5
h("par", f'StoreParameter("gap", "{GAP}")')
h("units", 'With Units\n.Geometry "mm"\n.Frequency "GHz"\n.Time "ns"\nEnd With')
h("freq", f'Solver.FrequencyRange "{FREQS[0] - 0.5:g}", "{FREQS[-1] + 0.5:g}"')
h("background", 'With Background\n.Type "Normal"\n.Epsilon "1"\n.Mu "1"\n'
  '.XminSpace "15"\n.XmaxSpace "15"\n.YminSpace "15"\n.YmaxSpace "15"\n.ZminSpace "15"\n.ZmaxSpace "15"\nEnd With')
h("boundary", 'With Boundary\n.Xmin "expanded open"\n.Xmax "expanded open"\n.Ymin "expanded open"\n'
  '.Ymax "expanded open"\n.Zmin "expanded open"\n.Zmax "expanded open"\n.ApplyInAllDirections "False"\nEnd With')
if KIND == "pt":      # без металла CST не строит сетку: маленький вакуумный кубик задаёт расчётную область
    h("box", 'With Brick\n.Reset\n.Name "box"\n.Component "air"\n.Material "Vacuum"\n'
      f'.Xrange "{-S*5*LAM-2:.3f}", "{S*5*LAM+2:.3f}"\n.Yrange "{-S*5*LAM-2:.3f}", "{S*5*LAM+2:.3f}"\n'
      f'.Zrange "-2", "{S*4.2*LAM+2:.3f}"\n.Create\nEnd With')
for k, (name, x, y, z) in enumerate(elems, start=1):
    if KIND == "dip":
        for arm, sg in (("a", 1), ("b", -1)):
            z1, z2 = f"{z:.4f} + {sg}*gap/2", f"{z:.4f} + {sg}*6.25"
            h(f"arm {name}{arm}", f'''With Cylinder
.Reset
.Name "{name}{arm}"
.Component "dip"
.Material "PEC"
.OuterRadius "0.4"
.InnerRadius "0"
.Axis "z"
.Zrange "{z1 if sg > 0 else z2}", "{z2 if sg > 0 else z1}"
.Xcenter "{x:.4f}"
.Ycenter "{y:.4f}"
.Segments "0"
.Create
End With''')
    h(f"port {name}", f'''With DiscretePort
.Reset
.PortNumber "{k}"
.Type "SParameter"
.Label "{name}"
.Impedance "50.0"
.SetP1 "False", "{x:.4f}", "{y:.4f}", "{z:.4f} - gap/2"
.SetP2 "False", "{x:.4f}", "{y:.4f}", "{z:.4f} + gap/2"
.LocalCoordinates "False"
.Monitor "True"
.Radius "0"
.Create
End With''')
for f in FREQS:
    h(f"monitor {f}", f'''With Monitor
.Reset
.Name "ff_{f:g}"
.Domain "Frequency"
.FieldType "Farfield"
.Frequency "{f:g}"
.Create
End With''')
h("mesh", '''With Mesh
.MeshType "PBA"
.SetCreator "High Frequency"
.AutomeshRefineAtPecLines "True", "4"
.UseRatioLimit "True"
.RatioLimit "20"
.LinesPerWavelength "10"
.MinimumStepNumber "10"
.Automesh "True"
End With''')
h("solver", '''With Solver
.Method "Hexahedral"
.CalculationType "TD-S"
.StimulationPort "All"
.StimulationMode "All"
.SteadyStateLimit "-35"
.MeshAdaption "False"
.AutoNormImpedance "True"
.NormingImpedance "50"
.SParaSymmetry "False"
.StoreTDResultsInCache "False"
End With''')
prj.save(PRJ, allow_overwrite=True)
print("built", PRJ, len(elems), "elements", flush=True)
t0 = time.time()
m.add_to_history("run solver", "Solver.Start")
print(f"solved in {time.time()-t0:.0f} s", flush=True)
prj.save()
t1 = time.time()
export_all(prj, OUT, FREQS, STEP, names=[e[0] for e in elems], th_range=TH)
print(f"exported in {time.time()-t1:.0f} s", flush=True)
prj.save()
print("DONE", flush=True)
