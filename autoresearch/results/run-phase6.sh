#!/bin/bash
set -u
JAR=$(ls lib/build/libs/*-jmh.jar)
OUT=autoresearch/results/phase6-floatbuild-tess3.txt
: > "$OUT"
echo "Phase 6: V24 (float build + double scan, C1) vs DISTINCT_PACKED_V3 @ tess 3, consume=AREA" | tee -a "$OUT"
echo "jar=$JAR  start=$(date -Is)" | tee -a "$OUT"
for T in 1 16; do
  echo "" | tee -a "$OUT"
  echo "===== threads=$T  load-before=$(cat /proc/loadavg) =====" | tee -a "$OUT"
  java --add-modules jdk.incubator.vector -Xmx4g -jar "$JAR" SurfaceBench \
    -f 3 -wi 4 -w 2s -i 6 -r 2s -t $T \
    -p variantId=V24,DISTINCT_PACKED_V3 -p tess=3 -p consume=AREA 2>&1 \
    | grep -E "Benchmark|V24|DISTINCT_PACKED_V3|Run complete" | grep -vE "^#" | tee -a "$OUT"
  echo "----- threads=$T  load-after=$(cat /proc/loadavg) -----" | tee -a "$OUT"
done
echo "" | tee -a "$OUT"; echo "DONE end=$(date -Is)" | tee -a "$OUT"
