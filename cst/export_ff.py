"""Экспорт ДН всех портов через FarfieldPlot.GetList (θ 0…180°, φ 0…359°, шаг 1°,
фаза относительно глобального нуля, Directivity).

ВАЖНО (CST 2024): ДН считается только в той сессии, где в памяти есть гекс-сетка, т.е. сразу после
решателя. После переоткрытия проекта CalculateList падает с «Farfield: No HEX mesh found».
Поэтому export_all() вызывается из build_array.py; автономный запуск
  python export_ff.py <проект.cst> [частоты] [шаг]
перед экспортом перезапускает решатель (Solver.Start).

Ключи GetList в CST 2024: "Spherical linear theta abs/phase", "Spherical linear phi abs/phase"
(без слова linear — «Invalid field component»).

Формат файла (читает tools/phase_center_fit.py):
  Theta [deg.]  Phi [deg.]  Abs(Dir.)[dBi]  Abs(Theta)[dBi]  Phase(Theta)[deg.]  Abs(Phi)[dBi]  Phase(Phi)[deg.]
"""
import sys, os, atexit, time
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface

NAMES = [f"rx{i}" for i in range(1, 9)] + [f"tx{i}" for i in range(1, 5)]
HEADER = ("Theta [deg.]  Phi [deg.]  Abs(Dir.)[dBi]  Abs(Theta)[dBi]  Phase(Theta)[deg.]  "
          "Abs(Phi)[dBi]  Phase(Phi)[deg.]" + chr(10) + "-" * 110 + chr(10))


def export_all(prj, out_dir, freqs=(9.8, 10.0, 10.2), step=1.0, names=NAMES):
    m = prj.model3d
    tmp = os.path.join(out_dir, "_ff_tmp.txt")
    tmp_vba = tmp.replace("\\", "\\\\")
    for k, name in enumerate(names, start=1):
        for f in freqs:
            ff = f"ff_{f:g} [{k}]"
            t0 = time.time()
            vba = f'''Sub Main
Dim th As Double, ph As Double, i As Long
Dim a As Variant, ta As Variant, tp As Variant, pa As Variant, pp As Variant, vt As Variant, vp As Variant
SelectTreeItem "Farfields\\{ff}"
With FarfieldPlot
.Reset
.SetPlotMode "Directivity"
.SetScaleLinear "False"
.Origin "zero"
For th = 0 To 180 Step {step:g}
  For ph = 0 To 360 - {step:g} Step {step:g}
    .AddListEvaluationPoint th, ph, 0, "spherical", "frequency", "-1"
  Next ph
Next th
.CalculateList "{ff}"
vt = .GetList("Point_T")
vp = .GetList("Point_P")
a = .GetList("Spherical abs")
ta = .GetList("Spherical linear theta abs")
tp = .GetList("Spherical linear theta phase")
pa = .GetList("Spherical linear phi abs")
pp = .GetList("Spherical linear phi phase")
End With
Open "{tmp_vba}" For Output As #1
For i = LBound(a) To UBound(a)
  Print #1, CStr(vt(i)) & " " & CStr(vp(i)) & " " & CStr(a(i)) & " " & CStr(ta(i)) & " " & CStr(tp(i)) & " " & CStr(pa(i)) & " " & CStr(pp(i))
Next i
Close #1
End Sub'''
            m._execute_vba_code(vba, timeout=None)
            out = os.path.join(out_dir, f"ff_{name}_{f:.1f}GHz.txt")
            with open(tmp, encoding="utf-8", errors="replace") as src, open(out, "w", encoding="utf-8") as dst:
                dst.write(HEADER)
                for line in src:
                    dst.write(line.replace(",", "."))
            print(f"exported {os.path.basename(out)} ({time.time()-t0:.0f} s)", flush=True)
    os.remove(tmp)


def export_touchstone(prj, path_noext):
    """S-матрица всех портов в Touchstone (<path>.s<N>p), нормировка 50 Ом."""
    p = path_noext.replace("\\", "\\\\")
    prj.model3d._execute_vba_code(f'''Sub Main
With TOUCHSTONE
.Reset
.FileName "{p}"
.Impedance "50"
.FrequencyRange "Full"
.Renormalize "True"
.UseARResults "False"
.Write
End With
End Sub''', timeout=None)


if __name__ == "__main__":
    PRJ = os.path.abspath(sys.argv[1])
    FREQS = [float(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else [9.8, 10.0, 10.2]
    STEP = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    OUT = os.path.dirname(PRJ)
    for f in os.listdir(OUT):
        if f.endswith(".lok"):
            os.remove(os.path.join(OUT, f))
    de = cst.interface.DesignEnvironment.new()
    de.set_quiet_mode(True)
    atexit.register(lambda: de.close())
    prj = de.open_project(PRJ)
    t0 = time.time()
    prj.model3d.add_to_history("run solver", "Solver.Start")
    print(f"solved in {time.time()-t0:.0f} s", flush=True)
    export_all(prj, OUT, FREQS, STEP)
    export_touchstone(prj, os.path.splitext(PRJ)[0])
    prj.save()
    print("DONE", flush=True)
