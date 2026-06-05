#!/bin/bash
set -u
JAR=$(ls lib/build/libs/*-jmh.jar)
OUT=autoresearch/results/phase2-float-tess3.txt
: > "$OUT"
echo "Phase 2: FLOAT_V2 vs DISTINCT_PACKED_V3 @ tess 3, -prof gc, consume=AREA" | tee -a "$OUT"
echo "jar=$JAR  start=$(date -Is)" | tee -a "$OUT"
for T in 1 4 8 16; do
  echo "" | tee -a "$OUT"
  echo "===== threads=$T  load-before=$(cat /proc/loadavg) =====" | tee -a "$OUT"
  java --add-modules jdk.incubator.vector -Xmx4g -jar "$JAR" SurfaceBench \
    -f 2 -wi 3 -w 2s -i 6 -r 2s -t $T \
    -p variantId=FLOAT_V2,DISTINCT_PACKED_V3 -p tess=3 -p consume=AREA \
    -prof gc 2>&1 | grep -E "Benchmark|FLOAT_V2|DISTINCT_PACKED_V3|gc.alloc.rate|gc.alloc.rate.norm|Iteration|Result|Run complete" | tee -a "$OUT"
  echo "----- threads=$T  load-after=$(cat /proc/loadavg) -----" | tee -a "$OUT"
done
echo "" | tee -a "$OUT"
echo "DONE end=$(date -Is)" | tee -a "$OUT"
