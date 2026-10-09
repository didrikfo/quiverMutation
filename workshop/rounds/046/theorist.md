# Bystander observation for an interior overlapping pair: lowering needs a run of three (k <= 2 sweep, k = 3 partial, margin 2); the one-step locality L1 holds for all LNAs n <= 10; H-010 stays unproved

author: theorist · round: 046 · kind: result (observations, no proof)
thread: T7 · bears on: H-010, F-022, F-024, F-025, E-025, E-029, E-066, E-122, H-011
scope: pair `(8:3)(9:3)` (plus `(8:4)(9:4)`, `(8:3)(9:4)` for a few placements), >= 6 empty vertices on the left of the pattern and >= 6 on the right; ONE bystander `(t:m)`, m = 2..4, left and right; mutations only at vertices within 2 of the pattern (margin 2), left and right mutations; k <= 2 full sweep, k = 3 only the cells listed. L1: every LNA of length 4..10, every single mutation. Nothing for two bystanders, k > 3, or margin > 2. Not a proof of anything about H-010.

## Response to referee

1. *22 of 22.* Wrong; it is **21 of 21** at k = 2 (23 rows of the table minus two runs of three; 7 + 7 + 7 for m = 2, 3, 4). Placements per cell at k = 3 are now stated (table below): 10 inert (m = 2: 7; m = 3, g = -1, 0, 1: 3), m = 3, g >= 2 and m = 4 not run at k = 3 (cap), 2 runs of three on the right (4 with the left).
2. *Title/names.* "Lemma L" is dropped; it is an observation. Title and Scope now carry the scope: k <= 2 full, k = 3 partial, margin 2, one bystander, both sides, pair `(8:3)(9:3)` with three further shapes sampled.
3. *F-022.* The prior statement is F-022's table: `(1:3)(2:3)(3:3)` reaches overlap 0; `(1:3)(2:3)(4:2)` and `(1:3)(2:3)(5:2)` stay at 2 (bystander sharing 1 or 0 arrows, inert); `(1:4)(2:4)(3:4)` only reaches 2 from 3 (run of three, not always sufficient); `(1:2)(2:3)(3:3)` stuck (first relation shares one arrow). E-025 qualifies it for long runs. What the sweep adds: the same statement for every gap g = -2..5 and m = 2..4 rather than three fixed placements, both sides, and k = 3 for some cells. Nothing deeper than F-022.
4. *`(8:3)(9:3)(10:4)` at k = 3.* It **lowers**: 4 reached, 1 lowered, min overlap 0 (23 s). At k = 2 it did not (3 reached, none lowered). So a run of three is not "inert up to k = 2 and hence inert": the lowering for `(10:4)` needs 3 mutations. The accurate statement is: in the data, no bystander sharing <= 1 arrow lowers, and every run of three tried lowers by k = 3 (5 of 5 tried at k = 3: right `(10:3)`, `(10:4)`, left `(7:3)`, `(6:4)`, and `(8:4)(9:4)(11:3)`). It is a necessary condition in every cell where it was tested, and sufficient for k = 3 in all 5 runs of three tried. Not proved, and E-025 says it fails for long runs.
5. *Locality.* Marked as **conjecture** in the claim below (the k-step window). What I ran is L1 (the one-step half): `theorist_locality.py`. For every LNA of length n = 4..10 (5, 14, 42, 132, 429, 1430, 4862 algebras), every vertex v and direction with `mutationIsPossibleAtVertex` true (n per LNA, 48 620 events at n = 10), the result is an LNA in exactly 2 of the n cases per LNA (9 724 at n = 10) and in **all** of those the changed relation starts lie in {v-2, v-1, v} (offset histogram identical in shape for every n; max distance 2, never 3). So single-mutation locality holds there. Caveats: right and left directions are pooled in the histogram (not split); the mutation is applied to LNAs only, whereas a k-step sequence passes through non-LNA algebras, so L1 does **not** give the k-step window; that remains conjecture. Also the pooled counts (2 LNA results per LNA at every n) are regular enough to deserve an explanation; I did not look for one (`theorist_whichv.py` shows which vertices: LNA results occur at v = 1 and n, 2 and n-1, ..., with counts 42,14,10,10,14,42 at n = 7).
6. *Left bystanders.* Run (`theorist_leftsweep.py 2 2 2 7`, t = 2..7, m = 2..4; the right side is the 045 table). Left is the mirror image: no lowering for any bystander sharing <= 1 arrow (14 placements); the left runs of three `(7:3)` (4 reached, 2 lowered, min 1) and `(6:4)` (3 reached, none lowered) at k = 2 behave exactly like `(10:3)` and `(10:4)`; at k = 3 `(7:3)` gives 6 reached / 4 lowered / min 0, `(6:4)` gives 4 / 1 / 0 -- identical to the right. This confirms the referee's variation and removes "not run".
7. *Status.* I no longer recommend closing T7. Recommendation to the Chair: "H-010 unproved; no counterexample; one-step locality L1 verified n <= 10; a proof needs an invariant on non-LNA intermediates, none found in two sittings." My earlier "no proof in reach" is withdrawn as too strong: it meant "I did not find one". On E-066: I conflated two things. E-066 is a rejection at step 7 of the *E-032 walk* (a commutative square, `c = [8,6,4]+[8,10,4]`, `c*(4->9) = 0`), and it is tagged H-010 in the record; H-010 itself says "a proof from the procedure's step 7", which is the procedure's step, not E-066. I withdraw the sentence that E-066 "is not an instance of H-010"; whether it bears on H-010 is not something I established. The F-025 line (cheaper line in H-010) is still not engaged.
8. *Another pair shape and the k = 3 cap rows.* Done in part: other shapes at k = 3 (below). The cap-cut rows (m = 3, g >= 2; m = 4) are not rerun.

## Claim

For an isolated pair `(8:3)(9:3)` planted with >= 6 empty vertices on each side, one added bystander relation, and at most 2 mutations at vertices within 2 of the pattern, the maximum overlap is lowered only by a bystander sharing >= 2 arrows with the pair (a run of three): 21 of 21 right-side placements sharing <= 1 arrow stay at 2, as do 14 left-side ones; at k = 3 the 10 placements that finished are inert, and all 5 runs of three tried (right `(10:3)`, `(10:4)`; left `(7:3)`, `(6:4)`; and `(8:4)(9:4)` + `(11:3)`) reach overlap 0 or lower than the start. The same holds for the other pair shapes sampled. One-step locality L1 (a mutation of an LNA that gives an LNA changes relation starts only within {v-2, v-1, v}) holds for all LNAs with n <= 10. Does not claim: anything about two bystanders, k > 3, margin > 2, long runs (E-025 refutes the naive version), or H-010 for all k. Refuted by: a bystander sharing <= 1 arrow that lowers the overlap in the same window, or an LNA mutation whose result is an LNA and changes a relation start at distance >= 3.

## Evidence

Pair `(8:3)(9:3)`, bystander `(t:m)`; "shared" = arrows common with the pair's arrows 8..11. Right-side rows are the 045 table (reproduced by the referee).

| side | m | placements | shared | k = 2 reached / lowered | k = 3 reached / lowered |
|---|---|---|---|---|---|
| right | 2 | g = -1..5 (7) | 1 or 0 | 4,5,6,6,6,6,6 / 0 | 9,11,11,12,12,12,12 / 0 |
| right | 3 | g = -1, 0, 1 (3) | 1, 0, 0 | 1, 2, 2 / 0 | 1, 2, 2 / 0 |
| right | 3 | g = 2..5 (4) | 0 | 2 / 0 | not run |
| right | 4 | g = -1..5 (7) | 1 or 0 | 1,2,2,2,2,2,2 / 0 | not run |
| right | 3 | `(10:3)` | 2 | 4 / 2 (min 1) | 6 / 4 (min 0) |
| right | 4 | `(10:4)` | 2 | 3 / 0 | **4 / 1 (min 0)** |
| left | 2 | t = 2..7 (6) | 0 or 1 | 6,6,6,6,5,4 / 0 | not run |
| left | 3 | t = 2..6 (5) | 0 or 1 | 2,2,2,2,1 / 0 | `(6:3)`: 1 / 0 |
| left | 3 | `(7:3)` | 2 | 4 / 2 (min 1) | 6 / 4 (min 0) |
| left | 4 | t = 2..5 (4) | 0 or 1 | 2,2,2,1 / 0 | pair shifted, `(5:4)`: 1 / 0 |
| left | 4 | `(6:4)` | 2 | 3 / 0 | 4 / 1 (min 0) |

Other pair shapes at k = 3 (`theorist_pairshapes_k3.txt`; the pair alone is the first row of each group):

| pattern | base | reached / lowered / min |
|---|---|---|
| `(8:4)(9:4)` | 3 | 2 / 0 / 3 |
| `(8:4)(9:4)(11:3)` (shares 2 with `(9:4)`) | 3 | 4 / 2 / 0 |
| `(8:4)(9:4)(12:3)`, `(14:3)` (share 1, 0) | 3 | 1, 2 / 0 / 3 |
| `(8:3)(9:4)` | 2 | 3 / 0 / 2 |
| `(8:3)(9:4)(13:3)` (gap) | 2 | 3 / 0 / 2 |

L1 (`theorist_locality_4_10.txt`): n = 4..10, events 20, 70, 252, 924, 3432, 12870, 48620; LNA results 10, 28, 84, 264, 858, 2860, 9724; offset (changed start - v) takes values -2, -1, 0 only, in ratio 1:2:1 at every n.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/045/theorist_bystander.py 2 2          # 2 min (referee reproduced)
timeout 10m .venv/bin/python workshop/rounds/046/theorist_single.py 3 2 10:4        # 23 s
timeout 10m .venv/bin/python workshop/rounds/046/theorist_leftsweep.py 2 2 2 7      # 1 min
timeout 10m .venv/bin/python workshop/rounds/046/theorist_single.py 3 2 7:3         # 14 s (also 6:4, 6:3, 5:4)
PAIR=8:4,9:4 timeout 10m .venv/bin/python workshop/rounds/046/theorist_single.py 3 2 11:3   # 24 s
timeout 10m .venv/bin/python workshop/rounds/046/theorist_locality.py 4 10          # 3.5 min
timeout 10m .venv/bin/python workshop/rounds/046/theorist_whichv.py                 # seconds
```

## Prior record

F-022 table (prior statement), E-025 (long runs), E-029 item 8, F-024, H-010 text (the isolation hypothesis; F-025 "cheaper line" untouched), E-122 (gate J_i), E-066 (tagged H-010; see item 7). New and modest: the gap/length/side/shape sweep and the L1 one-step check. Not found in `research/`: a statement of one-step locality across all LNAs (grep "locality", "distance 1", "changed start": nothing), but the 2-per-LNA regularity is unexplained and might be a known fact I did not find.

## Code changed

None to `quivermutation/`; no tests run. New scripts in `workshop/rounds/046/`: `theorist_single.py`, `theorist_leftsweep.py`, `theorist_locality.py`, `theorist_whichv.py`, and their outputs.

## Next

- Theorist: explain why exactly 2 mutations per LNA give an LNA, and prove L1 from the gate (E-122): it is a statement about one step on a line.
- Theorist/skeptic: k-step window as a theorem needs control of non-LNA intermediates (e.g. compute, for all k <= 3 sequences from an LNA, whether support of the change is within 2k+1 of the vertices used, over all algebras on the way).
- Skeptic: the cap-cut k = 3 rows (m = 3 g >= 2; m = 4) with margin 1; two bystanders.
- Chair: keep T7 open with H-010 SUPPORTED, unproved.
