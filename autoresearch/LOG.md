# Autoresearch log — FasterMolecularSurface (tess 2 & 3)

Append-only phase journal. Each phase: hypothesis · cheap kill-experiment · what was done · measurement
(with load-before) · verdict · meta-lesson · commit. Newest at the bottom.

See `KICKOFF.md` for the mission and operating rules (it is rewritten each phase to reflect the latest
state and lead ranking).

---

## Phase 0 — kickoff (setup)

- Created the `autoresearch` branch and this folder. Authored `KICKOFF.md` (self-updating mission +
  operating rules + prioritized leads, scoped to tess 2 & 3).
- Inherited state: bit-exact default `DistinctPackedNumericalSurfaceV3`; tess-2 float
  `FloatNumericalSurfaceV2`. Wins so far: A6 (SIMD build), C1 (float build). Negatives: A2, A7 (double +
  float), C3, naive A1, A4.
- Biggest open lead: the LUT-bitmask / fully-buried-atom prize (61% of atoms / 57% of scan time at tess 2).
- Next: Lead #1 — prototype the cheap cap→direction bitmask + instrument its cost before building a variant.
