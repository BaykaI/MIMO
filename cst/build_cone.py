"""Конусная MIMO-решётка в CST 2024 (геометрия как в chap_ofn.cone_geometry): усечённый конус, ось z, основание в z = 0,
четыре кольца диполей на боковой поверхности; Rx — кольца 1 и 3 (10 + 8), Tx — кольца 2 и 4 (8 + 6). Диполи вдоль z
(вертикальная поляризация, в азимуте всенаправленные). Фаза ДН — от центра основания (глобальный ноль).

  python build_cone.py [L_dip мм] [линий/λ] [r_wire мм]        напр. python build_cone.py 12.5 10 0.4
Результаты: cst\cone\ (проект, ff_<port>_<f>GHz.txt, .s32p, layout_mm.txt, build.log). Порты: rx1…rx18, tx1…tx14.
"""
import sys, os, time, atexit, math
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface
from export_ff import export_all, export_touchstone

L_DIP = float(sys.argv[1]) if len(sys.argv) > 1 else 12.5
LPW = sys.argv[2] if len(sys.argv) > 2 else "10"
R_WIRE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.4
GAP = 0.5
LAM = 30.0
FREQS = (9.8, 10.0, 10.2)
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "cone"); os.makedirs(OUT, exist_ok=True)
PRJ = os.path.join(OUT, "mimo_cone.cst")

rings = [(5.0 * LAM, 0.0, 10, "rx"), (4.17 * LAM, 1.4 * LAM, 8, "tx"), (3.33 * LAM, 2.8 * LAM, 8, "rx"), (2.5 * LAM, 4.2 * LAM, 6, "tx")]
tx, rx = [], []
for R, z, n, kind in rings:
    for i in range(n):
        a = 2 * math.pi * i / n + (math.pi / n if kind == "tx" else 0.0)
        (tx if kind == "tx" else rx).append((R * math.cos(a), R * math.sin(a), z))
elems = [(f"rx{i+1}",) + p for i, p in enumerate(rx)] + [(f"tx{i+1}",) + p for i, p in enumerate(tx)]
with open(os.path.join(OUT, "layout_mm.txt"), "w", encoding="utf-8") as fh:
    fh.write("name x_mm y_mm z_mm\n" + "".join(f"{n} {x:.4f} {y:.4f} {z:.4f}\n" for n, x, y, z in elems))

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
h("freq", 'Solver.FrequencyRange "9", "11"')
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
export_all(prj, OUT, FREQS, 1.0, names=[e[0] for e in elems])
export_touchstone(prj, os.path.splitext(PRJ)[0])
print("touchstone exported", flush=True)
prj.save()
print("DONE", flush=True)
