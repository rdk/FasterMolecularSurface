#!/bin/bash
set -u
JAR=$(ls lib/build/libs/*-jmh.jar)
OUT=autoresearch/results/phase8-points-path.txt
: > "$OUT"
echo "Phase 8: V3 AREA vs POINTS @ tess 2&3, -t 1/16, -prof gc" | tee -a "$OUT"
echo "jar=$JAR start=$(date -Is)" | tee -a "$OUT"
for T in 1 16; do
  echo "" | tee -a "$OUT"; echo "===== threads=$T load-before=$(cat /proc/loadavg) =====" | tee -a "$OUT"
  java --add-modules jdk.incubator.vector -Xmx4g -jar "$JAR" SurfaceBench \
    -f 2 -wi 3 -w 2s -i 6 -r 2s -t $T \
    -p variantId=DISTINCT_PACKED_V3 -p tess=2,3 -p consume=AREA,POINTS \
    -prof gc 2>&1 | grep -E "Benchmark|DISTINCT_PACKED_V3|gc.alloc.rate.norm|Run complete" | grep -vE "^#" | tee -a "$OUT"
  echo "----- threads=$T load-after=$(cat /proc/loadavg) -----" | tee -a "$OUT"
done
echo "" | tee -a "$OUT"; echo "DONE end=$(date -Is)" | tee -a "$OUT"
