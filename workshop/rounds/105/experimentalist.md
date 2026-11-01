# The defect is not constant per Coxeter-key class on n = 9 depth-5 data (it ranges 1..4 among members with cords >= 1, 1..5 counting the 61 cord-free members, inside one key class); constant only on the defect <= 1 filter of E-203/E-221, where constancy is made by the filter (key classes: only 2 occur)

author: experimentalist · round: 105 · kind: negative
thread: T6 · bears on: H-017
scope: n = 9 only · defect range 1..4 for cords >= 1 members (645 - 61 = 584), 1..5 only with the 61 cords = 0 members · the 9 LNAs outside every quipu class · forward walk depth 5 (both orientations, Coxeter guard on) · all monomial quipu members reached, any path length, cords >= 0 · 645 members, 2 Coxeter-key classes · the E-203/E-221 member files hold only cords >= 1, L = 5, defect <= 1 (36 members), so the walk was re-collected, not read from them · no depth 6, no n >= 10
power: control D (planted error): in the filtered set (36 members, defect all 1, test says YES), one member's defect flipped to 2 -> test says NO (1 of 1 class non-constant). Control B (planted defect = f(key)) -> YES (0 of 2 non-constant). Shuffled key labels do NOT discriminate here (2 of 2 classes non-constant with real and with shuffled labels, 1000 shuffles): with 2 classes and defect 1..5 in both, the shuffle has no power on this data.
backlog: n/a

## Response to referee

1. cords = 0 count and cords >= 1 histogram: done. 61 of 645 members have cords = 0; histogram for cords >= 1 now in Claim (class 1 {1:70 2:200 3:152 4:118}, class 2 {2:26 3:6 4:12}); title and Scope say 1..4 for cords >= 1, 1..5 only with the 61. Verified from the saved JSON (recomputed defect = rels - cords).
2. Reproduction time: done. Re-ran all nine collects: 3033030 39 s, the others 74-126 s each (parallel, 4 cores, 126 s total wall). The "4-5 min" is replaced. The re-collected JSON gave the same `analyse` output.
3. Referee's "constant only on the filter is trivial": agreed, kept as the explanation, not as a finding.

## Claim

Defect = relations - cords of a monomial quipu member. Over all 645 quipu members reached by the depth-5 walk from the nine outside LNAs, the defect is not a function of the Coxeter key. 61 of the 645 members have cords = 0 (defect = relations; class 1: 14 at 3, 38 at 4, 8 at 5; class 2: 1 at 4); all 9 members with defect 5 are among them. With cords >= 1 (584 members, the population of H-017's table): class (1,1,0,0,...) {1: 70, 2: 200, 3: 152, 4: 118}, class (1,1,-2,-2,...) {2: 26, 3: 6, 4: 12}; non-constant in both, so the verdict is unchanged. All 645: key class (1,1,0,0,...) (8 LNAs, 600 members) has defects {1: 70, 2: 200, 3: 166, 4: 156, 5: 8}; key class (1,1,-2,-2,...) (3033030 alone, 45 members) has {2: 26, 3: 6, 4: 13}. The only sense in which the saved E-203/E-221 data show "constant defect" is the selection defect <= 1 (with cords >= 1, L = 5): all 36 members have defect exactly 1, but that is the filter, which kept the minimum. What stays plausible, and is what H-017 needs, is that the minimum defect per class is >= 1: observed minimum 1 (class 1) and 2 (class 2, a single LNA); 4444400 alone has minimum 2 at depth 5 (1 at depth 6, E-065). Not claimed: any statement about depth >= 6, n != 9, non-monomial members, or that min defect >= 1 is proved (it is not).

Every member's key equals its source LNA's key (guard working: 0 mismatches of 645), so "key class" here equals "set of LNAs sharing a key": only two classes exist among the nine, so the census has two degrees of freedom and cannot say more about key-dependence than this.

## Evidence

| key class | LNAs | members | defect histogram | min |
|---|---|---|---|---|
| (1,1,-2,-2,...) | 3033030 | 45 | 2:26 3:6 4:13 | 2 |
| (1,1,0,0,...) | 3345000 3455000 3505000 4444400 4550400 5040330 5050030 5504030 | 600 | 1:70 2:200 3:166 4:156 5:8 | 1 |

Per LNA minimum (all members): 1 for 3345000, 3455000, 3505000, 4550400, 5040330, 5050030, 5504030; 2 for 4444400 and 3033030. So the minimum is not even constant inside the big key class at depth 5. Opposite pairs have identical histograms (3345000/5040330, 3455000/5504030, 3505000/5050030), a consistency check on the dual walk.

Cross-check against the saved data: my filter (cords >= 1, L = 5, defect <= 1) gives per-LNA counts 3345000: 3, 3455000: 3, 3505000: 9, 4550400: 6, 5040330: 3, 5050030: 9, 5504030: 3, 4444400: 0; this equals the member lists in `rounds/069/maverick_defect1_*.txt` and `rounds/077/maverick_defect1_rest_*.txt` (36 members; all defect 1; cells (1,2), (2,3), (3,4)).

Controls (from `analyse`): planted defect = f(key): 0 of 2 non-constant (YES, correct). Planted single error in the filtered set: 1 of 1 non-constant (NO, correct). Classes split by cords: 6 of 8 sub-classes non-constant (defect varies with cords at fixed key; expected, since the cell (c, c+1) is what the filter picked). Shuffle: no power (see power line). Honest reading: the statistic "number of key classes with non-constant defect" is saturated (2 of 2); the verdict NO is carried by the histograms, which show defect 1..5 in one class, not by a comparison with a null.

## Reproduction

```
for n in 3033030 3345000 3455000 3505000 4444400 4550400 5040330 5050030 5504030; do
  timeout 10m .venv/bin/python workshop/rounds/105/experimentalist_defect_census.py collect $n; done   # in parallel, measured: 3033030 39 s; the other eight 74-126 s each; all nine in parallel on 4 cores 126 s wall
timeout 10m .venv/bin/python workshop/rounds/105/experimentalist_defect_census.py analyse   # seconds
```
JSON outputs `workshop/rounds/105/experimentalist_defect_<LNA>.json`, 68 KB in total.

## Prior record

HYPOTHESES.md H-017 already tabulates the (cords, relations) cells per LNA at depth 5 (e.g. 3345000: (1,2) ... (3,5)), so defect 1..4 within an LNA is recorded; E-065 says the Coxeter polynomial cannot see (cords, relations) and the minimum depends on depth. What is new: the per-member census tied to the key class, the explicit two-class structure, the 36-member match with E-203/E-221, and that the minimum differs between 4444400 and its key-mates. E-221 itself says its controls test walk invertibility only; unchanged. This does not test H-017's claim (min defect >= 1 over all depths).

## Code changed

New `workshop/rounds/105/experimentalist_defect_census.py` only. No library change, no tests run.

## Next

Skeptic/theorist: "constant defect per key class" is false as stated and should be dropped from the T6 line; the live question is whether the minimum is >= 1 per class, which depth-5 data cannot decide (4444400 min 2 at depth 5, 1 at depth 6). Experimentalist, only with an overnight slot: depth-6 minima for all 9 LNAs (E-083: ~5e4 nodes each, shards of ~10 min).
