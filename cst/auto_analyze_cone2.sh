#!/bin/bash
# Разбирает каждый готовый расчёт конуса, как только в build.log появится DONE.
cd "$(dirname "$0")"
while true; do
  for d in cone_pt_d* cone_dip_d* cone_dip_wb; do
    [ -d "$d" ] || continue; case "$d" in *failed*) continue;; esac
    grep -q "^DONE" "$d/build.log" 2>/dev/null || continue
    [ -f "$d/c2_result.json" ] && continue
    bands="0,20"; [ "$d" = "cone_dip_wb" ] && bands="0,10,20,30,40,50"
    PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=2 python analyze_cone2.py "$d" "$bands" > "$d/analyze.log" 2>&1
    echo "ANALYZED $d rc=$? $(date +%H:%M)"
  done
  grep -q "QUEUE DONE" cone2_queue2.log && ! ls -d cone_*/ | grep -v failed | while read x; do [ -f "${x}c2_result.json" ] || echo m; done | grep -q m && { echo "ALL ANALYZED"; exit 0; }
  sleep 60
done
