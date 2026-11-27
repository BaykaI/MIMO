#!/bin/bash
# Шаг при неизменном размере конуса: число точечных излучателей на кольцах меняется, полоса 9-11 ГГц (5 частот).
cd "$(dirname "$0")"
PY="/c/Program Files (x86)/Microsoft Visual Studio/Shared/Python39_64/python.exe"
export CONE_TH="60,120" PT_GAP="2.0"
run() {  # d counts
  tag="cone_fa_d$1"
  if grep -q "^DONE" "$tag/build.log" 2>/dev/null; then echo "skip $tag"; return; fi
  for t in 1 2; do
    mkdir -p "$tag"; echo "START $tag $(date +%H:%M) попытка $t"
    PYTHONIOENCODING=utf-8 "$PY" build_cone2.py pt 1 9:11:0.5 "$tag" 2 "$2" > "$tag/build.log" 2>&1
    echo "END $tag rc=$? $(date +%H:%M)"
    if grep -q "^DONE" "$tag/build.log"; then
      PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=3 python analyze_cone2.py "$tag" 0,20 > "$tag/analyze.log" 2>&1
      echo "ANALYZED $tag rc=$? $(date +%H:%M)"; return
    fi
    mv "$tag" "${tag}_failed$t"
  done
}
run 3.14 10,8,8,6
run 2 16,13,10,8
run 1.5 21,17,14,10
run 1 31,26,21,16
run 0.75 42,35,28,21
run 0.5 63,52,42,31
echo "QUEUE DONE $(date +%H:%M)"
