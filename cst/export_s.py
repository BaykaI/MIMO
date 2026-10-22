"""Выгрузка S-матрицы из уже посчитанного проекта CST в Touchstone (<проект>.s<N>p).
Сетка в памяти не нужна — 1D-результаты хранятся в проекте.

  python export_s.py <проект.cst>
"""
import sys, os, atexit
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
import cst.interface
from export_ff import export_touchstone

PRJ = os.path.abspath(sys.argv[1])
OUT = os.path.dirname(PRJ)
for f in os.listdir(OUT):
    if f.endswith(".lok"):
        os.remove(os.path.join(OUT, f))
de = cst.interface.DesignEnvironment.new()
de.set_quiet_mode(True)
atexit.register(lambda: de.close())
prj = de.open_project(PRJ)
export_touchstone(prj, os.path.splitext(PRJ)[0])
print("DONE", flush=True)
