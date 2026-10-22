"""Плоское кольцо MIMO (ДЗ1, kolso.txt) в CST 2024: 20 Rx на радиусе R, 10 Tx на 1,35·R, все в плоскости z = 0,
диполи вдоль y, ДН по нормали z. Масштаб задаётся радиусом Rx в длинах волн (в статье 50λ — не считается).

  python build_ring.py <R_rx в λ> [L_dip мм] [линий/λ] [r_wire мм]   напр. python build_ring.py 5 ; python build_ring.py 10 12.5 8 0.5
Результаты: cst\ring_R<R>\ (проект, ff_<port>_<f>GHz.txt, .s30p, build.log). Порты 1–20 = rx1…rx20 (угол 18°·(i−1)),
21–30 = tx1…tx10 (угол 36°·(j−1)). λ = 30 мм (10 ГГц).
"""
import sys, os, time, atexit, math
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface
from export_ff import export_all, export_touchstone

R_LAM = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
L_DIP = float(sys.argv[2]) if len(sys.argv) > 2 else 12.5
LPW = sys.argv[3] if len(sys.argv) > 3 else "12"       # линий сетки на длину волны
R_WIRE_ARG = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
LAM = 30.0
R_RX, R_TX = R_LAM * LAM, 1.35 * R_LAM * LAM
N_RX, N_TX = 20, 10
R_WIRE, GAP = R_WIRE_ARG, 0.5
FREQS = (9.8, 10.0, 10.2)
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, f"ring_R{R_LAM:g}")
os.makedirs(OUT, exist_ok=True)
PRJ = os.path.join(OUT, "mimo_ring.cst")

elems = [(f"rx{i+1}", R_RX * math.cos(2 * math.pi * i / N_RX), R_RX * math.sin(2 * math.pi * i / N_RX)) for i in range(N_RX)]
elems += [(f"tx{j+1}", R_TX * math.cos(2 * math.pi * j / N_TX), R_TX * math.sin(2 * math.pi * j / N_TX)) for j in range(N_TX)]
with open(os.path.join(OUT, "layout_mm.txt"), "w", encoding="utf-8") as fh:
    fh.write("name x_mm y_mm\n" + "".join(f"{n} {x:.4f} {y:.4f}\n" for n, x, y in elems))

for f in os.listdir(OUT):
    if f.endswith(".lok"):
        os.remove(os.path.join(OUT, f))
de = cst.interface.DesignEnvironment.new()
de.set_quiet_mode(True)
atexit.register(lambda: de.close())
prj = de.new_mws()
m = prj.model3d
n = [0]


def h(name, vba):
    n[0] += 1
    m.add_to_history(f"{n[0]:03d} {name}", vba)


for name, val in (("l_dip", L_DIP), ("r_wire", R_WIRE), ("gap", GAP), ("r_rx", R_RX), ("r_tx", R_TX)):
    h(f"par {name}", f'StoreParameter("{name}", "{val}")')
h("units", 'With Units\n.Geometry "mm"\n.Frequency "GHz"\n.Time "ns"\nEnd With')
h("freq", 'Solver.FrequencyRange "9", "11"')
h("background", 'With Background\n.Type "Normal"\n.Epsilon "1"\n.Mu "1"\n'
  '.XminSpace "15"\n.XmaxSpace "15"\n.YminSpace "15"\n.YmaxSpace "15"\n.ZminSpace "20"\n.ZmaxSpace "20"\nEnd With')
h("boundary", 'With Boundary\n.Xmin "expanded open"\n.Xmax "expanded open"\n.Ymin "expanded open"\n'
  '.Ymax "expanded open"\n.Zmin "expanded open"\n.Zmax "expanded open"\n.ApplyInAllDirections "False"\nEnd With')
for k, (name, x, y) in enumerate(elems, start=1):
    for arm, s in (("a", 1), ("b", -1)):
        y1, y2 = f"{y:.4f} + {s}*gap/2", f"{y:.4f} + {s}*l_dip/2"
        h(f"arm {name}{arm}", f'''With Cylinder
.Reset
.Name "{name}{arm}"
.Component "dip"
.Material "PEC"
.OuterRadius "r_wire"
.InnerRadius "0"
.Axis "y"
.Yrange "{y1 if s > 0 else y2}", "{y2 if s > 0 else y1}"
.Xcenter "{x:.4f}"
.Zcenter "0"
.Segments "0"
.Create
End With''')
    h(f"port {name}", f'''With DiscretePort
.Reset
.PortNumber "{k}"
.Type "SParameter"
.Label "{name}"
.Impedance "50.0"
.SetP1 "False", "{x:.4f}", "{y:.4f} - gap/2", "0"
.SetP2 "False", "{x:.4f}", "{y:.4f} + gap/2", "0"
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
.AutomeshRefineAtPecLines "True", "6"
.UseRatioLimit "True"
.RatioLimit "30"
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
export_all(prj, OUT, FREQS, 1.0, names=[e[0] for e in elems])
export_touchstone(prj, os.path.splitext(PRJ)[0])
print("touchstone exported", flush=True)
prj.save()
print("DONE", flush=True)
