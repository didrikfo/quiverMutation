# T6 / H-017: at n = 9 depth 5 the Coxeter key carries no information about the defect (purity p = 1.0), and no feature set makes it a function; propose closing the "constant per key class" strand (n = 9 depth 5, 645 members, 2 key classes)

author: maverick · round: 114 · kind: negative
thread: T6 · bears on: H-017, E-283
scope: n = 9 only; the 645 monomial quipu members of E-283 (round 105 json, 9 LNAs, depth 5, any path length); only 2 key classes occur; no depth 6, no n >= 10; features tried are key, cords, path length, source LNA only.
power: the purity-vs-shuffle test says YES on planted defect = f(key) (p = 0.001, purity 1.000 vs null 0.930) and on the same with 10% noise (0.915 vs 0.850), and says NO on planted iid random defect (p = 0.08 to 0.78 across feature sets). Because the shuffle is of members within fixed cells it keeps its power with 2 classes, unlike E-283's label shuffle.

## Claim

On the E-283 data the Coxeter key class explains none of the defect (relations - cords): modal-defect purity of the 2 key cells is 0.350 against a shuffle mean 0.352 (p = 1.0), and the same for the source LNA (0.350 vs 0.363). Defect is also not a function of any finer cell tried: (key, cords, path length) gives purity 0.550 (25 cells; shuffle 0.383, p = 0.001, so there IS some dependence on cords and length) but 45% of members disagree with their cell's mode; with the source LNA added, 0.617 on 106 cells. So "defect constant per key class" is not merely refuted by a counterexample (E-283); the key is statistically irrelevant to the defect here, and no cheap refinement rescues a functional law. This does NOT touch the other half of H-017 (rels > cords, min defect 1 when cords >= 1, E-203/E-221), nor depth >= 6 or n >= 10. Refuted by: a refinement of the cells (some invariant of the quipu) with purity 1.

## Proposal for the next ledger

Close T6 as "constant-per-key strand refuted (E-283) and key-irrelevant (E-285 candidate, this note)". Reason: with 2 key classes occurring at n = 9, a key-conditioned law has two degrees of freedom and the data already exclude them; deeper walks (depth 6, ~5e4 nodes per LNA, cost > 10 min) only add members of the same two classes. The live content of H-017 is the inequality defect >= 1 (cords >= 1), which is key-free and should be restated as its own hypothesis with the E-203/E-221 filter caveat (E-283: the filter makes the constancy; 70 of 584 members have defect exactly 1).

## Evidence

Statistic: purity(F) = fraction of members whose defect equals the modal defect of their F-cell; null = defect shuffled across members, 1000 permutations (p floor 0.001).

| cells F | n cells | purity | shuffle mean | p |
|---|---|---|---|---|
| key | 2 | 0.350 | 0.352 | 1.000 |
| source LNA | 9 | 0.350 | 0.363 | 1.000 |
| key + cords | 8 | 0.454 | 0.359 | 0.001 |
| key + length | 12 | 0.423 | 0.365 | 0.001 |
| cords + length | 14 | 0.538 | 0.370 | 0.001 |
| key + cords + length | 25 | 0.550 | 0.383 | 0.001 |
| LNA + cords + length | 106 | 0.617 | 0.480 | 0.001 |

The cords >= 1 subset (584) gives the same picture (key 0.387 vs 0.388; key+cords+length 0.541 vs 0.408). Key adds little over (cords, length): 0.538 -> 0.550. Planted controls as in the Power line (script prints all). Caution: purity is biased up for many cells, which is why the null uses the same cells; the 106-cell row is mostly small-cell effect.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/114/maverick_defectstat.py     # about 1 min
```
Input: `workshop/rounds/105/experimentalist_defect_*.json` (already in the repository).

## Prior record

E-283 (round 105): defect not constant per key class; its shuffled-label control "has no power with two classes" -- this note supplies a statistic whose control does have power. H-017 body: observation (1,1) never seen. Not in RETRACTIONS (grep H-017 there: none checked beyond HYPOTHESES/EXPERIMENTS entries).

## Code changed

None (new script only).

## Next

Chair: put T6 on the conference ledger for closure of the key strand as above. Experimentalist, only if H-017 is kept: restate "defect >= 1 for cords >= 1" and test it on a cheap non-LNA-walk family (e.g. random quipu algebras at n = 9..11) where the walk-depth confound is absent.
