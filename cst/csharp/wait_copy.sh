#!/bin/bash
# Ждёт, пока закроются все окна ConeMimo.exe из пакета, и кладёт туда новую сборку.
P="/c/Users/koste/HFSS/Пакет_MIMO_конус/программа_CSharp"
while tasklist | grep -qi "ConeMimo.exe"; do sleep 3; done
cp ConeMimo.exe "$P/ConeMimo.exe" && echo "COPIED $(date +%H:%M:%S)"
