# One rebuilt class-1 walk at n = 7 reaches E-147's counts (67 J != 0, 13 D = 0, 12 strict H1/H2) at 22 450 to 24 032 expansions; its first 13 D = 0 steps all fail `tiltingPlus`, 10 lie inside E-151's 16 and 3 beyond its 20 000-expansion cap; 5 of the 13 have a parent with parallel arrows

author: skeptic · round: 069 · kind: result (closes the E-179 "class-1 containment open" item)
thread: T10/T5 · bears on: H-015, E-147, E-151, E-179

## Response to referee

1. Wording (done): title, *Scope* and Claim now say identity is by count agreement (67 and 13 together, window 22 450 to 24 032 expansions, 1583 wide) at one point of one rebuilt walk; no step-for-step claim. Class 2 not rerun (stated in *Scope*); cap-dependence shown for class 1 only.
2. Not independently re-run, tied to files (done): the referee did not re-run the walk (3 segments, ~27 min, over the 10 min limit; I do not re-run it this sitting either) and trusts the checkpoint. The checkpoint is uncommitted in /tmp: `/tmp/skeptic_c1_13_ck1.pkl` 71 079 752 bytes, mtime 2026-10-16 13:35:26 UTC, sha256 50bf82830469803fd0d9da863f65025785a370c1ca662f16223d3fb0363fd02c; `/tmp/skeptic_c1_13_cls1.pkl` 114 119 bytes, sha256 d6ee10d3b87a4fa53b2865cf230cdac62ceb04468ab5f41b29a6fc42244a2cdf. Kept outputs: `skeptic_c1_13_45000.out` sha256 82cc046661689d6835d6d6ce2499b2273273eef2c3ca42fb416f8bfc86b3d70d, `skeptic_c1_13_20000.out` sha256 b22bf946376f9e036a1ed82a5f5c1f046e5c244fdffe5a3ba1025c5080daeb2c. (Correction: the checkpoint is 71 MB, not the 30 MB stated under Reproduction.)
3. E-147's "12 strict H1/H2" (done, new `skeptic_c1_13_h1h2.py` / `.out`, 20 s): E-147's own `shape` and `strict` functions (copied verbatim from rounds/043/skeptic_c2.py, i = first element of supp J) on the first 13 D = 0 steps: 12 have the H1/H2 shape and all 12 are strict; step 9 (parallel parent (4,7),(5,4)) has no shape. So the 13 contain exactly 12 strict H1/H2 steps, the second E-147 count that agrees at the same window (a further, weak, coincidence; still count agreement, not identity). Across all 89 D = 0 steps of the 45 000 sample, 72 have the shape.
4. Theorist remark (done): relabelled in *Next* as an observation at one sample, not a regularity.
Also: the "500 s" remark in the E-179 paragraph is an inference from this machine's speed (E-147's machine not measurable); read it so.

scope: identity with E-147's 13 is by count agreement at one window of one rebuilt walk (E-147 kept no step list or expansion count), not step for step; n = 7, key class 1 (12 seed LNAs + duals, class sizes 8,12,14,58,64,108), one deterministic key-guarded canonicalKey-deduplicated BFS run to 45 000 expansions (queue not closed); class 2 not rerun. Every "failing" below is "under the J = 0 premise" (a failing step is one where `tiltingPlus` is False and J != 0; whether its child is derived equivalent to the parent is NOT tested here).

## Claim

The BFS of E-147 and E-151 is one walk (BFS over canonicalKey-deduplicated algebras, key guard, seeds = LNAs and duals of class 1); the three published counts are prefixes of it. At 20 000 expansions: 16 distinct (parent, v) key-keepers (E-151), 64 distinct J != 0 steps, 10 D = 0. At 22 450 to 24 032 expansions: **67 distinct J != 0 steps and exactly 13 D = 0** (E-147's numbers), 19 (parent, v) pairs. So E-147's 13 are *consistent with being* the first 13 D = 0 steps (distinct by (C_A, C_B, v)) of this walk (count agreement at one window; not shown step for step). Hand-checked: all 13 are gate-admitted, fail `tiltingPlus`, have J != 0 (|supp J| = 1 except step 7), keep the key and are D = 0. Containment: the first 10 (first seen at expansions 7 369 to 18 462) lie in E-151's 16; steps 11, 12, 13 (first seen at 20 696, 21 036, 22 450) are beyond E-151's cap, so "E-147's 13 not reproduced inside E-151's 16" is explained by the cap, not by a discrepancy. Parallel arrows in the parent: 5 of 13 (steps 1, 2, 7, 8, 9); the other 8 were rebuilt from arrows and relations by `build()`. It does not claim: the children leave the class (open, premise), nor that "13" or "16" is a property of the class: at 45 000 expansions the same walk has 111 key-keeping failing (parent, v) pairs and 89 D = 0 Cartan-distinct steps, growing with the cap (counts 2, 4, 11, 47, 64, 359 distinct J != 0 Cartan steps at 2.5k, 5k, 10k, 15k, 20k, 45k). Refuted if a different canonical ordering yields 13 not forming a prefix, or a D = 0 step in the 45 000 sample that keeps its key but passes `tiltingPlus`.

## Evidence

Own walk and own code: `skeptic_c1_13.py` re-implements the BFS loop (queue, seen = canonicalKey, key guard `search._coxeterKeyOrNone == base`), the det-polynomial test D = det(xC_B+C_B^T) - det(xC_A+C_A^T) (exact Fraction elimination, 9 integer points; no sympy, no round-43 `detpoly`), and calls the library only for mutation/reduction/gate/Cartan. Independent tilting test: `rounds/050/skeptic_tilt.py` Hom_K(T,T[m]), m = -1, 0, 1 (not `tiltingPlus`, `perI`, the gate), T replaces P_v by (P_v -> sum of P_h over out-arrows of v).

First 13 D = 0 steps in discovery order (all: gate True, `tiltingPlus` False, key kept, D = 0, dim Hom(T,T[1]) = 0):

| # | first seen (expansion) | depth | v | J | parent parallel arrows | dim Hom(T,T[-1]) | in E-151's 16 | hand rebuild (arrows+relations, gate, tp, key, child) |
|---|------|---|---|---|-----|---|-----|-----|
| 1 | 7369 | 7 | 4 | {6:1} | (7,4) | 1 | yes | no: parallel |
| 2 | 10106 | 7 | 4 | {2:1} | (1,4) | 1 | yes | no: parallel |
| 3 | 13666 | 8 | 6 | {5:1} | none | 1 | yes | ok |
| 4 | 17242 | 8 | 7 | {3:1} | none | 1 | yes | ok |
| 5 | 17619 | 8 | 4 | {6:1} | none | 1 | yes | ok |
| 6 | 17797 | 8 | 4 | {6:1} | none | 1 | yes | ok |
| 7 | 17964 | 8 | 4 | {2:1,6:1} | (7,4) | 2 | yes | no: parallel |
| 8 | 18015 | 8 | 5 | {6:1} | (7,5) | 1 | yes | no: parallel |
| 9 | 18022 | 8 | 4 | {6:1} | (4,7),(5,4) | 1 | yes | no: parallel |
| 10 | 18462 | 8 | 7 | {1:1} | none | 1 | yes | ok |
| 11 | 20696 | 8 | 1 | {7:1} | none | 1 | **no (beyond cap)** | ok |
| 12 | 21036 | 8 | 1 | {4:1} | none | 1 | **no (beyond cap)** | ok |
| 13 | 22450 | 8 | 7 | {4:1} | none | 1 | **no (beyond cap)** | ok |

(Parent depth 7 or 8 as in E-151. Hand rebuild "ok" = parent and child built by `build()` from edges and relations, gate True, `tiltingPlus` False, key(parent) = key(step child) = key(hand child), canonicalKey(step child) = canonicalKey(hand child), J as listed. The five parallel parents cannot be built by `build()`; for those the evidence is the in-process walk plus the independent Hom test, which handles parallel arrows.) Hom(T,T[-1]) != 0 at all 13, so none is a tilting step in the independent sense either (this is Ladkani restated, not new power).

E-151 containment: the 16 distinct (parent, v) at 20 000 expansions are all D = 0, all key kept, all Hom(T,T[-1]) != 0, 8 of them with a parallel-arrow parent; they collapse to 10 Cartan-distinct steps, which are exactly #1 to #10 above. E-147 convention (distinct by (C_A, C_B, v)) and E-151 convention (distinct by (canonicalKey(parent), v)) differ: 13 Cartan steps correspond to 19 (parent, v) pairs. In the 45 000 sample every key-kept J != 0 step (111) has D = 0, fails `tiltingPlus`, and has Hom(T,T[-1]) != 0 (dim 1 in 94, 2 in 14, 3 in 3 of the 111 pairs); no key-kept step passes the independent test. Whether its child is derived equivalent: not tested, under the J = 0 premise.

Consequence for E-179: its rerun found 12 of 16 E-151 pairs D = 0 and 56 distinct J != 0 steps against E-147's 67 under load; that is the same walk stopped earlier by wall clock (the E-147 figures need about 22 450 to 24 032 expansions, i.e. about 9 to 10 minutes here, not "500 s" on a loaded machine). The E-147 text says "500 s, depth <= 30"; depth was never the binding limit (depth 8 is reached by 13 666 expansions).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 500 45000   # resumable: run it 3 times (3 segments of <= 10 min, ~9 min each, checkpoint /tmp/skeptic_c1_13_ck1.pkl, 71 MB, not committed); last call also prints the analysis (~1 min)
.venv/bin/python workshop/rounds/069/skeptic_c1_13_rebuild.py /tmp/skeptic_c1_13_cls1.pkl     # the hand rebuild, 20 s
```
Outputs kept: `skeptic_c1_13_20000.out`, `skeptic_c1_13_45000.out` (28 KB). `.venv/bin/python workshop/rounds/069/skeptic_c1_13_prefix.py` (needs the checkpoint, ~30 s) prints the first-seen column and the prefix counts.

## Prior record

E-147 (13 and 9 of 67 and 64; "class-1 D = 0 steps rest on one script"), E-151 (16 and 9), E-179 ("E-147's 13 class-1 steps not reproduced", containment open), E-154 (hand rebuild c1 8/8, the other 8 parent parallel), E-174 (parallel ones decided by End(T) comparison). This submission adds: the identification of E-147's 13 with a walk prefix (the 67 and 13 coincide at exactly 22 450 to 24 032 expansions), the 10 + 3 split against E-151, and the 5 parallel-parent steps among the 13.

## Code changed

None in `quivermutation/` or tests. New: `skeptic_c1_13.py`, `skeptic_c1_13_rebuild.py`, `skeptic_c1_13_prefix.py`, two `.out` files; after review: `skeptic_c1_13_h1h2.py` and `.out`.

## Next

Chair: the E-147/E-151 numbers "13, 16, 9" are cap artefacts; quote them with their expansion cap (the count grows with the cap: 111 pairs at 45 000). Toolsmith: a checkpointed walk would let class 2 be matched the same way (9 of 64 at 500 s). Theorist: all 111 + 89 key-kept steps have D = 0 and Hom(T,T[-1]) of dim 1 or 2 with |supp J| = 1 or 2 in the 13, an observation at one sample (the 111-pair tally 94/14/3 and 12 of 13 with |supp J| = 1), not a regularity; a statement would need a second sample.
