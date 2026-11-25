#!/bin/bash
# Фактор полосы отдельно: точечные излучатели (вибраторы λ/20 на 10 ГГц), полоса 5-15 ГГц шаг 0,5, шаг кольца фиксирован.
cd "$(dirname "$0")"
PY="/c/Program Files (x86)/Microsoft Visual Studio/Shared/Python39_64/python.exe"
for d in "$@"; do
  out="ringz_d${d}_ptwb"
  if grep -q "^DONE" "$out/build.log" 2>/dev/null; then echo "skip $d"; continue; fi
  mkdir -p "$out"; echo "START $d $(date +%H:%M)"
  PYTHONIOENCODING=utf-8 "$PY" build_ringz.py "$d" 10 0.1 1.5 0.3 _ptwb 5:15:0.5 > "$out/build.log" 2>&1
  echo "END $d rc=$? $(date +%H:%M) $(grep -cE '&H|Traceback' "$out/build.log") ошибок"
done
echo "QUEUE DONE $(date +%H:%M)"
