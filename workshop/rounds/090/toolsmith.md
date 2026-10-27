# The 7 + 7 F-010 NO at n = 9 (join97) dropped 9 and 11 children for a None key, not 26 353 / 15 401, none can hide a meet; the NO still has no power beyond the 3304000 move component

author: toolsmith · round: 090 · kind: result (correction of a README figure, plus a power statement)
thread: T10/T5 · bears on: F-010, E-219, E-226, E-233, E-234, E-238, E-244
scope: n = 9, one pair (3060000, 3304000), J = 0 key-kept F/R moves, balls to depth 7 (join97 run). None-key children found by re-expanding the complete depth-6 frontier of both balls (19 999 + 11 355 parents, not a sample). No new meet attempted; J != 0 steps still not walked; 2400230 not rerun; join97 depth-7 checkpoints not available.

## Claim

1. **README misreading.** `workshop/overnight/2026-10-26-join97/README.md` gives "children dropped for a `None` key: 26 353 / 15 401". The log prints `... expansions N None-key children dropped M secs ...` (`report` in `skeptic_f010join.py`: `'expansions', nexp, 'None-key children dropped', none`). So 26 353 and 15 401 are the **expansion counts** (equal to the ball sizes through depth 6, hence the look of a coincidence); the drops are **M = 9 (A) and 11 (B)**, out of 284 665 and 151 733 kept children. The README's "removes a large number of children" is false for this run.
2. **The drops are depth-7 leaves.** Counter 0 after all depth <= 5 parents (my rebuilt depth-6 balls, both roots, log line `dropped 0`; consistent with the join97 logs). `canonicalKey` returns None when the product of factorials of the parallel-bundle sizes exceeds the cap 720 (`fingerprint.py` l. 316-322), which depends only on bundle sizes, so a None algebra is None in both balls. All other vertices of depth <= 6 have keys, so a J = 0 path of total <= 14 is missed only if its meet vertex is a depth-7 vertex of both balls with a None key. (A path of total <= 13, or a different meet vertex, is seen.)
3. **That blind spot is empty.** Re-expanding all depth-6 parents gives exactly 9 (A) and 11 (B) None-key children, matching the logs, in 8 and 11 isomorphism classes (bundle shapes such as [2,2,3,3,4], [2,2,6], [3,4,6]); none of the 88 A-B class pairs is isomorphic under the meet's own test (WL hash then `iso`). The `canonicalKey` gap of E-160 is not what limits the 7 + 7 NO at n = 9.
4. **Power.** The limits are the move set and the J = 0 restriction, not the key. (a) E-226/E-238: 2400230 ~ 3304000 is certified (quipu theorem, Hom-verified `lnaMoves`) but E-226 found no meet at 6 + 8, total 14. A depth-d ball is the J = 0 key-kept geodesic ball, so no meet at 6 + 8 and none at 7 + 7 are the same statement about total 14: there is a **known false NO at total 14 inside this ball family**. The inequivalence of 3060000 and 3304000 rests on F-010 (E-057, Ladkani Cor 3.15) and E-238, not on the join. (b) YES controls reach total 6 (E-226, cross-class) and 7-8 (E-233, E-234, E-244: a child against its own root); nothing was ever joined at 9-14, so there is no positive control at the tested radius. (c) Gate-admitted moves refused for J != 0: 4 802 of 289 467 (1.66%) and 2 253 of 153 986 (1.46%), pre-key counts (E-229). The key-kept ones are the real perI J != 0 steps of E-234/E-239/E-244 (A: 10 of 4 770 depth-5 parents; B: 1 of 2 958 at depth 5, 6 of 11 355 at depth 6); every one of the 17 observed children rejoined its own root by J = 0 moves at total 6-8. For A at depth 6 I did not count them; by those rates "tens" would be skipped (estimate only).

**For T10 under the J = 0 premise.** The NO is the expected answer (independent certificate) and is **not evidence about the premise**: only J = 0 edges were walked, each taken as a derived equivalence; an empty meet is what inequivalence predicts, and a meet would have been the only informative outcome. Proposed E-entry wording: "no key-kept J = 0 meet at 7 + 7 (balls 117 839 / 62 230); None-key children 9 / 11 (README's 26 353 / 15 401 are expansion counts), all depth-7, 8 + 11 iso classes, no A-B match, so no hidden meet; a false NO at total 14 is known (2400230 ~ 3304000, E-226 6 + 8); the NO says nothing on the J = 0 premise."

## Evidence

| quantity | A = 3060000 | B = 3304000 |
|---|---|---|
| expansions at the end (log) = ball size at depth 6 | 26 353 | 15 401 |
| None-key children dropped (log counter) | 9 | 11 |
| None-key children found by re-expanding all depth-6 parents | 9 (shards 2, 2, 5, 0) | 11 (4, 2, 0, 5) |
| depth-6 parents re-expanded / children generated | 19 999 / 217 602 | 11 355 / 112 440 |
| drops while building depth 0-6 (rebuilt) | 0 | 0 |
| isomorphism classes of the dropped | 8 | 11 |
| A-class isomorphic to a B-class | 0 of 88 pairs | |

Parent index, move and bundle shapes per dropped child: `workshop/rounds/090/toolsmith_nonekey_out.txt`. The rebuilt depth-6 balls have the logged sizes (26 353, 15 401), which is the check that the rebuild is the join97 ball. Weak points: (i) the iso is the meet's own (relation nodes, signs not tracked), not the sign-aware one; (ii) "None depends only on bundle sizes" is a reading of `canonicalKey`, not tested by mutating a child; (iii) depth 7 was not rebuilt, only the children of depth-6 parents were scanned, which suffices because the drops are exactly those.

## Reproduction

```
D=<scratch dir>
for t in A:3060000 B:3304000; do timeout 10m .venv/bin/python -u -W ignore workshop/rounds/079/skeptic_f010join.py side x ${t#*:} 6 540 $D/ck_${t%%:*}.pkl & done; wait   # A needs the call repeated once; about 10 + 10 min wall, in parallel
for s in 0 1 2 3; do timeout 10m .venv/bin/python -u -W ignore workshop/rounds/090/toolsmith_nonekey.py scan $D/ck_A.pkl $s 4 470 $D/nA_$s.pkl & done; wait   # repeat once with budget 200; about 10 + 3 min; B same with ck_B.pkl and nB_, one call, 7 min
timeout 5m .venv/bin/python -u -W ignore workshop/rounds/090/toolsmith_nonekey.py report $D/nA_?.pkl -- $D/nB_?.pkl   # about 1 min
```
4 cores, about 10 parents per second per process.

## Prior record

E-219/E-226 (no meet 5 + 5, 6 + 6), E-233/E-234/E-239/E-244 (J != 0 steps rejoin), E-238 (NO has power only inside the 3304000 component), E-160 (None-key blind side, n = 7). `grep` of `research/EXPERIMENTS.md` finds no entry for the 7 + 7 run yet, and the 9 / 11 correction is not recorded. Claim 4(a) restates E-226/E-238; that 6 + 8 and 7 + 7 coincide at total 14 follows from the balls being BFS balls of one graph (the None-key and hash details are closed by claims 2-3).

## Code changed

New `workshop/rounds/090/toolsmith_nonekey.py` (scan and report modes, resumable) and its output `toolsmith_nonekey_out.txt` (small). No library code touched; no tests apply.

## Next

- Chair: correct the join97 README row and sentence; record the E-entry with the wording above.
- Toolsmith/theorist: a NO with power needs a calibration: an equivalent key-equal pair at n = 8 whose J = 0 distance exceeds 8, or the key-kept perI J != 0 steps added as edges to test whether 2400230 then meets 3304000.
- Experimentalist: count the key-kept perI J != 0 steps among the 19 999 depth-6 parents of 3060000 (the "tens" above is unmeasured).
