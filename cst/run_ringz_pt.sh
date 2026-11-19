#!/bin/bash
# Эталонная серия: те же кольца, короткие вибраторы λ/20 (точечные источники).
cd "$(dirname "$0")"
PY="/c/Program Files (x86)/Microsoft Visual Studio/Shared/Python39_64/python.exe"
for d in "$@"; do
  out="ringz_d${d}_pt"
  if grep -q "^DONE" "$out/build.log" 2>/dev/null; then echo "skip $d"; continue; fi
  mkdir -p "$out"; echo "START $d $(date +%H:%M)"
  PYTHONIOENCODING=utf-8 "$PY" build_ringz.py "$d" 10 0.1 1.5 0.3 _pt > "$out/build.log" 2>&1
  echo "END $d rc=$? $(date +%H:%M) $(grep -cE '&H|Traceback' "$out/build.log") ошибок"
  PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=2 python analyze_ringz.py "$out" > "$out/analyze.log" 2>&1
  echo "ANALYZED $d rc=$?"
done
echo "QUEUE DONE $(date +%H:%M)"
