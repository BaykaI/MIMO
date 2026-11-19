#!/bin/bash
# Очередь колец с разным шагом: по одному расчёту CST за раз, готовые пропускаются.
cd "$(dirname "$0")"
PY="/c/Program Files (x86)/Microsoft Visual Studio/Shared/Python39_64/python.exe"
for d in "$@"; do
  out="ringz_d$d"
  if grep -q "^DONE" "$out/build.log" 2>/dev/null; then echo "skip $d"; continue; fi
  mkdir -p "$out"
  echo "START $d $(date +%H:%M)"
  PYTHONIOENCODING=utf-8 "$PY" build_ringz.py "$d" > "$out/build.log" 2>&1
  echo "END $d rc=$? $(date +%H:%M) $(grep -cE '&H|Traceback' "$out/build.log") ошибок"
  PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=2 python analyze_ringz.py "$out" > "$out/analyze.log" 2>&1
  echo "ANALYZED $d rc=$?"
done
echo "QUEUE DONE $(date +%H:%M)"
