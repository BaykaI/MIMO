#!/bin/bash
# Очередь конуса (вторая): ждёт текущий расчёт cone_pt_d0.5, затем перезапуск cone_pt_wb (прерван на порту 25) и остальное.
cd "$(dirname "$0")"
PY="/c/Program Files (x86)/Microsoft Visual Studio/Shared/Python39_64/python.exe"
while kill -0 4048 2>/dev/null; do sleep 30; done
echo "END cone_pt_d0.5 $(date +%H:%M) $(grep -c '^DONE' cone_pt_d0.5/build.log) DONE"
run() {  # kind s fspec tag
  if grep -q "^DONE" "$4/build.log" 2>/dev/null; then echo "skip $4"; return; fi
  for t in 1 2; do
    mkdir -p "$4"; echo "START $4 $(date +%H:%M) попытка $t"
    PYTHONIOENCODING=utf-8 "$PY" build_cone2.py "$1" "$2" "$3" "$4" > "$4/build.log" 2>&1
    echo "END $4 rc=$? $(date +%H:%M) $(grep -cE '&H|Traceback' "$4/build.log") ошибок"
    grep -q "^DONE" "$4/build.log" && return
    mv "$4" "${4}_failed$t"
  done
}
run pt 1 7.5:12.5:0.25 cone_pt_wb
for d in 0.5 0.75 1 1.5 2; do
  s=$(python -c "print(round($d/3.14159265,5))"); run pt "$s" 9:11:0.25 "cone_pt_d$d"
done
run dip 1 7.5:12.5:0.25 cone_dip_wb
for d in 0.5 0.75 1 1.5 2; do
  s=$(python -c "print(round($d/3.14159265,5))"); run dip "$s" 9:11:0.25 "cone_dip_d$d"
done
echo "QUEUE DONE $(date +%H:%M)"
