# autoresearch META — durable lessons across rounds (READ FIRST)

Accumulated wisdom that outlives any single round. Update this at the end of every round: add closed
leads, new meta-lessons, and tooling gotchas. A new agent reads this **before** arming a round so it
never re-treads closed ground or re-learns the operating rules the hard way.

---

## 1. Permanently-closed leads — do NOT re-try (with the one-line reason)

| Lead | Round | Why it's dead |
|---|---|---|
| Neighbor-major **cap→direction bitmask** / fully-buried-atom culling (§1b) | 2026-06-05 | Floor is only ~12–18% of scalar tests, but the baseline is 4-wide SIMD; a sound cheap LUT is 85–92% false positives → gather-bound, loses to SIMD. |
| **DCLM** second dot-lattice (A5) | 2026-06-05 | Same neighbor-major spatial-narrowing family; survivor saving already inside the bitmask total, loses to SIMD. |
| **Region-binned hint** (C3, earlier) | pre-2026-06-05 | Single last-occluder hint was already enough; per-octant was ~2–3% slower. |
| **Float SCAN at tess ≥ 3, multi-thread** (C9) | 2026-06-05 | Vector-API float (8-lane) intrinsics deopt to boxing under concurrency (13.6 MB/op → 1.08 GB/op at ≥2 threads, 7–23× slower). NOT source-fixable. Float surfaces are **tess-2-only**. |
| **Analytic Gauss-Bonnet** per-atom SASA (B1) | 2026-06-05 | Area-only (no point set); ~21% per-atom off the sampled surface; O(neighbors²) trig, slower. |
| **Power-diagram / weighted-alpha-complex** exact SASA | 2026-06-05 | Area-only (no point set); robust regular/power Delaunay is native-C++-only (CGAL GPLv3 / Geogram BSD-3), JNI from JVM. |
| **Tighter neighbor grid** (A2) | pre-2026-06-05 | Smaller cells → bigger stencil + bookkeeping; lost. |
| **Morton/Hilbert atom relabeling** | 2026-06-05 | Modest; V18 already sequentializes candidate-coordinate reads (the main locality win); output inverse-permutation may eat the gain; gains "likely small" for a cache-resident already-cell-ordered build. |
| **Float BUILD (C1/V24) at tess 3** | 2026-06-05 | Bit-identical on the corpus and wins at tess 2, but build-only wins dilute at tess 3 where the double scan dominates. Tess-2-only. |
| **Point-materialization (`consume=POINTS`) path** | 2026-06-05 | Zero-copy `surfacePointsXYZ()`/`PackedSurfaceAccess` is free — POINTS ≈ AREA in time and allocation at every tess/thread cell. No headroom. |
| A7 padded-tail (double V22 + float V25) | pre-2026-06-05 | Padded loop defeats Vector intrinsics. |

**Net:** the bit-exact SASA **scan is a dead well** and the **build is tapped** for tess-2/3 single-protein
point-set throughput. As of 2026-06-05 the surface is at a **strong local optimum**; new wins need a
**scope change** (batch/GPU C4, single-protein latency C5, an area-only consumer, or a new operating point).

## 2. Meta-lessons (how to evaluate ideas — these generalize)

- **Counting headroom ≠ wall-clock headroom when the baseline is vectorized.** Decisive lens: a
  scalar-test-count win must beat a **4-wide-SIMD sequential** baseline → ~4 scalar tests ≈ 1 unit of
  wall-clock, and **gather is strictly worse than sequential**. Every scan alternative that wins on test
  count but turns sequential-SIMD into gather has lost.
- **A Vector-API speedup that holds single-threaded can INVERT under concurrency** if the intrinsics
  deopt to the boxing fallback (C9). Always validate vectorized variants at the **deployment thread count**
  and watch `gc.alloc.rate.norm`, not just single-thread wall-clock.
- **An optimization's value = (its phase's share of total) × (its speedup).** Build-only wins **dilute at
  higher tessellation** because the per-direction scan grows (42→162→642 dirs) and dominates. Always ask
  "what fraction of total is this phase, at the *target* operating point?" before investing.
- **De-risk cheaply FIRST.** Most negatives here were killable by a load-immune *count*, not a built
  variant. Prototype the cap-mask / instrument the cost / `-prof gc` probe before scaffolding + benchmarking.
- **The fully-buried / survivor structure is real but unexploitable:** fully-buried atoms are cheap
  *because* their directions bury fast (early-exit); survivors are expensive but you can't skip a survivor's
  full scan without first paying a spatial pre-filter that costs more than it saves.

## 3. Machine & tooling gotchas (this box)

- **Shared box** (`lab`, AMD Ryzen 9 9950X, 16C/32T). 1-min load swings **1.0–25** through the day.
  **Poll `/proc/loadavg`; benchmark only if 1-min < 1.5.** Counting/instrumentation/correctness/code-reading
  are load-immune — do them anytime; defer timing to idle windows.
- **`perf` hardware counters are BLOCKED** (`perf_event_paranoid=4`) and **async-profiler is not installed** —
  `-prof async`/`perfnorm` fail. **`-prof gc` works** and `alloc.rate.norm` (bytes/op) is the load-immune
  proxy that diagnosed the C9 boxing. Cache-miss deltas can't be measured directly here.
- **JMH:** `./gradlew jmhJar` builds `lib/build/libs/*-jmh.jar`. Test-only additions don't change the jar,
  so an existing jar stays valid for production-variant benchmarks. Add a variant with one line in
  `SurfaceCatalog` (test source). Invoke: `java --add-modules jdk.incubator.vector -Xmx4g -jar …-jmh.jar SurfaceBench -f 3 -wi 4 -w 2s -i 8 -r 2s [-t 16] -p variantId=… -p tess=2,3 -p consume=AREA|POINTS`.
  Bracket every run with load-before/after; discard cells where load rose.
- **Scorecard test stdout is captured by gradle** — read it from
  `lib/build/test-results/scorecard/*<TestName>*.xml` (`<system-out>` CDATA), not the console.
- **Background benchmarks notify on completion** (harness-tracked). Don't poll them with short wakeups;
  set a long fallback heartbeat. Load-immune work can fill the wait.
- **Frozen-class policy** (root `CLAUDE.md`): never edit a frozen surface/strategy class. **Instrument via
  test-scope CLONES**; the sanctioned mutation point is new seam implementations on `DevSurfaceV1Soa`.

## 4. Reusable instrumentation (test scope, from the 2026-06-05 round — extend, don't rebuild)

- `BitmaskFeasibilityScan` + `BitmaskFeasibilityTest` + `CapDirectionLut` — neighbor-major occlusion
  cost model + a sound cube-face/cosθ cap→direction LUT with a soundness sweep.
- `DclmFeasibilityScan` + `DclmFeasibilityTest` — survivor-direction cost share + spatial-bin candidate count.
- `AnalyticGateTest` — sampled tess-N area vs a converged (≈analytic) reference, total + per-atom error.
- `BuildProfileTest` — phase-decomposed clone of the V3 SIMD build (grid / expandedR / SoA repack /
  distance pass / edges→CSR) + candidate-to-survivor ratio.
- `FloatBuildToleranceGateTest` — tolerance envelope of a float-build variant vs the double oracle.
- Pre-existing: `ScanInstrumentationTest`/`InstrumentingOcclusionScan`/`ScanStats`, `SurfaceAccuracy`,
  `VariantEquivalence`, `SurfaceQuality`.

All are `@Tag("scorecard")`: run with `./gradlew scorecard --tests '*<Name>'`.

## 5. Archive index (completed rounds)

| Round | Scope | Outcome | Report |
|---|---|---|---|
| `2026-06-05_tess2-3-throughput` | tess 2 & 3, single + 16t, point-set, bit-exact/tolerance | 8 phases, 0 promotable wins; surface at a strong local optimum; C9 confirmed; §1b prize closed | `archive/2026-06-05_tess2-3-throughput/REPORT.md` |
