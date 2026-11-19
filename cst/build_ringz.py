"""Кольцевая MIMO-система (пеленгатор по азимуту) в CST 2024 для двух проверок:
  1) широкополосный сигнал: 13 частот 8,5-11,5 ГГц (шаг 0,25), фазовые центры и групповые задержки по частоте;
  2) шаг между излучателями: радиус задаётся шагом приёмников по дуге d (в длинах волн на 10 ГГц).

Геометрия: плоскость z = 0, вертикальные (вдоль z) полуволновые диполи 12,5 мм, порт 50 Ом в разрыве.
Nrx = 12 приёмников на радиусе R = Nrx*d/(2*pi), Ntx = 6 передатчиков на радиусе 1,35*R, смещены на половину шага.
Обзор: горизонт (theta = 90°), азимут 0-359°. Выгрузка ДН только в полосе theta 60-120° (углы места ±30°).

  python build_ringz.py <d в λ> [линий/λ] [r_wire мм] [L_dip мм] [зазор мм] [суффикс]
  напр. python build_ringz.py 0.5 ; эталон с точечными источниками: python build_ringz.py 0.5 10 0.1 1.5 0.3 _pt
Результаты: cst\ringz_d<d>\ (проект, ff_<порт>_<f>GHz.txt, .s18p, layout_mm.txt, build.log).
Порты: 1-12 rx1…rx12, 13-18 tx1…tx6.
"""
import sys, os, time, atexit, math
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface
from export_ff import export_all, export_touchstone

D_LAM = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
LPW = sys.argv[2] if len(sys.argv) > 2 else "10"
R_WIRE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.4
L_DIP = float(sys.argv[4]) if len(sys.argv) > 4 else 12.5      # 1.5 = короткий вибратор λ/20 (точечный источник)
GAP = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
TAG = sys.argv[6] if len(sys.argv) > 6 else ""                  # суффикс папки, напр. _pt
LAM = 30.0
NRX, NTX = 12, 6
R_RX = NRX * D_LAM * LAM / (2 * math.pi)
R_TX = 1.35 * R_RX
FREQS = tuple(round(8.5 + 0.25 * i, 2) for i in range(13))
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, f"ringz_d{D_LAM:g}{TAG}"); os.makedirs(OUT, exist_ok=True)
PRJ = os.path.join(OUT, "ringz.cst")

elems = [(f"rx{i+1}", R_RX * math.cos(2 * math.pi * i / NRX), R_RX * math.sin(2 * math.pi * i / NRX), 0.0) for i in range(NRX)]
elems += [(f"tx{j+1}", R_TX * math.cos(2 * math.pi * (j + 0.5) / NTX), R_TX * math.sin(2 * math.pi * (j + 0.5) / NTX), 0.0) for j in range(NTX)]
with open(os.path.join(OUT, "layout_mm.txt"), "w", encoding="utf-8") as fh:
    fh.write("name x_mm y_mm z_mm\n" + "".join(f"{n} {x:.4f} {y:.4f} {z:.4f}\n" for n, x, y, z in elems))
print(f"d = {D_LAM}λ: R_rx = {R_RX:.2f} мм, R_tx = {R_TX:.2f} мм, шаг Rx по дуге {2*math.pi*R_RX/NRX:.2f} мм", flush=True)

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


for name, val in (("l_dip", L_DIP), ("r_wire", R_WIRE), ("gap", GAP)):
    h(f"par {name}", f'StoreParameter("{name}", "{val}")')
h("units", 'With Units\n.Geometry "mm"\n.Frequency "GHz"\n.Time "ns"\nEnd With')
h("freq", 'Solver.FrequencyRange "8", "12"')
h("background", 'With Background\n.Type "Normal"\n.Epsilon "1"\n.Mu "1"\n'
  '.XminSpace "15"\n.XmaxSpace "15"\n.YminSpace "15"\n.YmaxSpace "15"\n.ZminSpace "15"\n.ZmaxSpace "15"\nEnd With')
h("boundary", 'With Boundary\n.Xmin "expanded open"\n.Xmax "expanded open"\n.Ymin "expanded open"\n'
  '.Ymax "expanded open"\n.Zmin "expanded open"\n.Zmax "expanded open"\n.ApplyInAllDirections "False"\nEnd With')
for k, (name, x, y, z) in enumerate(elems, start=1):
    for arm, s in (("a", 1), ("b", -1)):
        z1, z2 = f"{z:.4f} + {s}*gap/2", f"{z:.4f} + {s}*l_dip/2"
        h(f"arm {name}{arm}", f'''With Cylinder
.Reset
.Name "{name}{arm}"
.Component "dip"
.Material "PEC"
.OuterRadius "r_wire"
.InnerRadius "0"
.Axis "z"
.Zrange "{z1 if s > 0 else z2}", "{z2 if s > 0 else z1}"
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
h("mesh", f'''With Mesh
.MeshType "PBA"
.SetCreator "High Frequency"
.AutomeshRefineAtPecLines "True", "4"
.UseRatioLimit "True"
.RatioLimit "20"
.LinesPerWavelength "{LPW}"
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
export_all(prj, OUT, FREQS, 1.0, names=[e[0] for e in elems], th_range=(60, 120))
print(f"exported in {time.time()-t1:.0f} s", flush=True)
export_touchstone(prj, os.path.splitext(PRJ)[0])
print("touchstone exported", flush=True)
prj.save()
print("DONE", flush=True)
