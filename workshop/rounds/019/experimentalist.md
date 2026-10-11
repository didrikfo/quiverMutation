# On every gate-admitted step tested the Cartan discrepancy is exactly -dim ker: nonzero iff non-tilting, with dim ker g_i never above 1 on the guarded walks

author: experimentalist · round: 019 · kind: result
thread: T5 · bears on: H-015, E-095, E-087, E-086

## Claim

On the guarded-BFS steps of n = 5, 6, 7 (class 0 at n = 6, 7, time-capped) and the E-080 family, dim ker (the map p -> (p beta)_beta, summed over i != k) is 0 on every tilting step and >= 1 on every non-tilting step: there is no non-tilting step with dim ker 0 and no tilting step with dim ker > 0 (0 of 1 056 non-tilting, 0 of 167 000+ tilting). So "Cartan discrepancy = -dim ker, nonzero iff non-tilting" holds on this data. On the three guarded walks the non-tilting dim ker is always exactly 1 (a single entry -1), and each rejecting parent has exactly one non-tilting vertex; totals 2 and 3 occur only in the hand-built E-080 family. Not claimed: any theorem; independence (the walks are not independent, E-086); coverage of other key classes at n = 6, 7; per-i values above 1 (never seen).

## Evidence

Script = scholar's `scholar_cartan_vs_tilt.py` plus a dim ker tally on every step (tilting too) and a parent tally (parent = `fingerprint.canonicalKey`). Disagreements tilt xor cong: 0 in all runs.

| run | steps tilting | steps non-tilting | distinct parents | parents with a non-tilting step | non-tilting with dim ker 0 | tilting with dim ker > 0 |
|---|---|---|---|---|---|---|
| E-080 family (18 algebras) | 75 | 6 | 18 | 6 | 0 | 0 |
| n = 5, both classes, closed | 30 300 | 0 | 11 700 | 0 | - | 0 |
| n = 6 class 0 (480 s cap, 39 921 algebras) | 83 591 | 907 | 22 869 | 907 | 0 | 0 |
| n = 7 class 0 (480 s cap, 28 921 algebras) | 53 502 | 143 | 13 564 | 143 | 0 | 0 |

Histogram of total dim ker (sum over i != k) on non-tilting steps:

| total dim ker | E-080 | n = 6 | n = 7 |
|---|---|---|---|
| 1 | 3 | 907 | 143 |
| 2 | 2 | 0 | 0 |
| 3 | 1 | 0 | 0 |
| >= 4 | 0 | 0 | 0 |

Per-i values g_i: on the walks, 3 628 zero and 907 equal to 1 (n = 6); 715 and 143 (n = 7); E-080: 22 zero, 10 equal to 1. On tilting steps every one of the 417 955 / 321 012 / 121 200 / 413 per-i values is 0. Each non-tilting step also satisfies E-095's identity (diff is row k off-diagonal, entry (k,i) = -dim ker g_i), counted in the raw output.

Reading. Since cong False means the diagonal-free discrepancy is nonzero and the identity says it equals -dim ker, "no non-tilting step with dim ker 0" is a corollary of E-095 on its data, and the converse (tilting => dim ker 0) is the new content: it is the same as "congruent". The step counts are not E-095's (696/111) because the caps were 480 s, not 420 s on a differently loaded machine; the n = 6 non-tilting number 907 versus E-095's 696 is a time effect, same pattern.

Observation worth the theorist's time: 1 is the only value seen on 1 050 walk steps. Either dim ker <= 1 holds for the A5-shaped rejecting parents that dominate (E-086 says all rejecting parents were of that shape; I did not check shape per parent here), or the guarded walk never reaches a parent with two bad predecessors at i. Parents are all distinct, so the 1 050 are not one repeated case.

## Reproduction

```
.venv/bin/python workshop/rounds/019/experimentalist_kerhist.py --e078     # 3 s
.venv/bin/python workshop/rounds/019/experimentalist_kerhist.py 5 --all    # 82 s
timeout 10m .venv/bin/python workshop/rounds/019/experimentalist_kerhist.py 6 --class 0 --budget-sec 480   # 8 min
timeout 10m .venv/bin/python workshop/rounds/019/experimentalist_kerhist.py 7 --class 0 --budget-sec 480   # 8 min
```
Outputs: `workshop/rounds/019/experimentalist_kerhist_{e078,n5,n6_c0,n7_c0}.txt`. The wall-clock cap makes n = 6, 7 counts non-deterministic (a time budget is not a verdict).

## Prior record

E-095 reports the identity and says the dim ker distribution and the number of distinct parents were not reported; this answers that. "dim ker never 0 on a non-tilting step" is implied by E-095 for its data and I find no other record (grep of `research/` for "dim ker" gives E-095 only). The A5-shape of rejecting parents is E-086. Nothing in RETRACTIONS bears on it.

## Code changed

New file `workshop/rounds/019/experimentalist_kerhist.py` (copy of the scholar script, tallies added). No library change, no tests needed.

## Next

- Theorist: prove "dim ker g_i in {0,1}" or find a parent where it is 2 (needs a non-A5 rejecting parent; the walks do not produce one).
- Toolsmith: a sampled search for gate-admitted parents with total dim ker >= 2 that are not hand-built; also n = 6 classes 1-3 and n = 7 classes 1+ (not run here; 480 s was class 0 only).
- Experimentalist next: classify the 1 050 rejecting parents by shape (A5 or not) to test the reading above.
