"""Построение MIMO-решётки 8 Rx + 4 Tx в CST Studio Suite 2024 (ДЗ2, этапы 2–3).

Запуск (Python 3.9 из Visual Studio, см. README):
  python build_array.py short   -> этап 2: короткие диполи 1,5 мм (точечные источники)
  python build_array.py half    -> этап 3: полуволновые диполи (длина L_DIP)

Создаёт проект cst\stage<N>\mimo_<stage>.cst, запускает TD-решатель по всем портам,
затем экспортирует ДН каждого порта на 9,8 / 10 / 10,2 ГГц в ASCII (θ 0…180°, φ 0…360°, шаг 1°)
с началом фазы в глобальном нуле (FarfieldPlot.Origin "zero").
Порты: 1–8 = Rx1…Rx8 (x = −105…105 мм, y = 0), 9–12 = Tx1…Tx4 (x = ∓127,5, ∓112,5, +112,5, +127,5; y = Y_TX).
"""
import sys, os, time, atexit
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface

STAGE = sys.argv[1] if len(sys.argv) > 1 else "short"
ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "stage2" if STAGE == "short" else "stage3")
os.makedirs(OUT, exist_ok=True)
PRJ = os.path.join(OUT, f"mimo_{STAGE}.cst")

RX_X = [-105, -75, -45, -15, 15, 45, 75, 105]
TX_X = [-127.5, -112.5, 112.5, 127.5]
Y_TX = 30.0
if STAGE == "short":
    L_DIP, R_WIRE, GAP = 1.5, 0.1, 0.3     # полная длина, радиус провода, зазор порта
else:
    L_DIP, R_WIRE, GAP = float(sys.argv[2]) if len(sys.argv) > 2 else 13.6, 0.25, 0.5
FREQS = (9.8, 10.0, 10.2)

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


def par(name, val):
    h(f"par {name}", f'StoreParameter("{name}", "{val}")')


par("l_dip", L_DIP); par("r_wire", R_WIRE); par("gap", GAP); par("y_tx", Y_TX)

h("units", 'With Units\n.Geometry "mm"\n.Frequency "GHz"\n.Time "ns"\nEnd With')
h("freq", 'Solver.FrequencyRange "9", "11"')
h("background", 'With Background\n.Type "Normal"\n.Epsilon "1"\n.Mu "1"\n'
  '.XminSpace "15"\n.XmaxSpace "15"\n.YminSpace "15"\n.YmaxSpace "15"\n.ZminSpace "20"\n.ZmaxSpace "20"\nEnd With')
h("boundary", 'With Boundary\n.Xmin "expanded open"\n.Xmax "expanded open"\n.Ymin "expanded open"\n'
  '.Ymax "expanded open"\n.Zmin "expanded open"\n.Zmax "expanded open"\n.Xsymmetry "none"\n'
  '.Ysymmetry "none"\n.Zsymmetry "none"\n.ApplyInAllDirections "False"\nEnd With')

# диполи вдоль y: два плеча (цилиндры PEC) + дискретный порт в зазоре
elems = [(f"rx{i+1}", x, 0.0) for i, x in enumerate(RX_X)] + [(f"tx{i+1}", x, Y_TX) for i, x in enumerate(TX_X)]
for k, (name, x, y) in enumerate(elems, start=1):
    for arm, s in (("a", 1), ("b", -1)):
        y1 = f"{y} + {s}*gap/2"
        y2 = f"{y} + {s}*l_dip/2"
        h(f"arm {name}{arm}", f'''With Cylinder
.Reset
.Name "{name}{arm}"
.Component "dip"
.Material "PEC"
.OuterRadius "r_wire"
.InnerRadius "0"
.Axis "y"
.Yrange "{y1 if s > 0 else y2}", "{y2 if s > 0 else y1}"
.Xcenter "{x}"
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
.SetP1 "False", "{x}", "{y} - gap/2", "0"
.SetP2 "False", "{x}", "{y} + gap/2", "0"
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
.AutomeshRefineAtPecLines "True", "6"
.UseRatioLimit "True"
.RatioLimit "30"
.LinesPerWavelength "15"
.MinimumStepNumber "10"
.Automesh "True"
End With''')
h("solver", '''With Solver
.Method "Hexahedral"
.CalculationType "TD-S"
.StimulationPort "All"
.StimulationMode "All"
.SteadyStateLimit "-40"
.MeshAdaption "False"
.AutoNormImpedance "True"
.NormingImpedance "50"
.CalculateModesOnly "False"
.SParaSymmetry "False"
.StoreTDResultsInCache "False"
.FullDeembedding "False"
.SuperimposePLWExcitation "False"
.UseSensitivityAnalysis "False"
End With''')

prj.save(PRJ, allow_overwrite=True)
print("built", PRJ, flush=True)

t0 = time.time()
m.add_to_history("run solver", "Solver.Start")
print(f"solved in {time.time()-t0:.0f} s", flush=True)
prj.save()

# экспорт ДН (в той же сессии — нужна гекс-сетка в памяти, см. export_ff.py)
from export_ff import export_all, export_touchstone
export_all(prj, OUT, FREQS, 1.0)
export_touchstone(prj, os.path.splitext(PRJ)[0])
print("touchstone exported", flush=True)
prj.save()
print("DONE", flush=True)
