# H1 (floating rules are length independent) holds for all 330 floating rules of window width 6..8 at one length each beyond the tests (0 failures in 43 116 filed applications; lengths 10-12 partly unchecked); the H6 ablation (no anchored rules) changes no verdict in 45 placements at n = 12, 13 in the reduced walk, where the anchored rules are expected to be redundant, but does in a rules-only walk

author: theorist · round: 051 · kind: result (null; two ablations of my own hypotheses)
thread: T4 · bears on: H-020, F-051, F-053, E-157, E-046/E-051
scope: H1: 0 failures for all 330 floating rules of width 6..8 of `lnaMoves.VERIFIED_MOVES` at one length each (12 for w = 6, 7; 13 for w = 8; w = 6, 7 rules also at 13 for 82 of 216), 43 116 filed applications; the tests (`lengthsToCheck`) cover w+1..w+4 for w <= 6 and w+1..w+2 for w >= 7, so lengths 11 (w = 6), 10-11 (w = 7) and 11-12 (w = 8) are unchecked, as are all lengths above 13 and widths 9..11 (54 rules) beyond w+2. H6: with the reduced walk (free, end edges, doubles) removing the anchored rules changes no orbit at 45 placements of 7 cores at n = 12, 13 (20 000-row cap, never hit); with a rules-only walk it does change them (core `45`, `skeptic_sensitivity.py`). Not n >= 14, not cores of three relations.

## Response to referee

1. Title, claim, scope fixed. "First length beyond the tests" was wrong: `tests/test_lna_moves.py:lengthsToCheck` gives w+1..w+4 only for w <= 6 and w+1..w+2 for w = 7, 8. Each rule was checked at one length (12 for w = 6, 7; 13 for w = 8; 13 additionally for 82 of the 216 w = 6, 7 rules). Lengths in the gaps, not checked: 11 for w = 6 (tests 7..10, checked 12); 10, 11 for w = 7 (tests 8, 9, checked 12); 11, 12 for w = 8 (tests 9, 10, checked 13). This is a spot check at one larger length, not length independence. The 82 are the first 20+22+20+20 rules of slices 0..3 of 8 of the w = 6, 7 list (last completed indices 19, 21, 19, 19 in the `13b` files, killed runs, no DONE line).
2. Prior record now cites F-032 (FINDINGS.md line 1388: the rule table on top of the double mutation and free move adds nothing at n = 10, same 262 rows in 16 orbits) and the finding "Which move does it matters, and the answer is not the rule table" (FINDINGS.md line 607: table alone / table + edges: no join at n = 10..16; table + edges + doubles: join). The "verdict is mostly free/edge/double moves" conclusion is theirs, not mine. Added here: only the anchored rules are dropped (the floating 414 kept), at 45 placements of 7 cores at n = 12, 13, orbit identity by row-set hash and verdict compared; no change. That narrower statement is the whole addition.
3. Stated: the ablation ran only in the reduced walk (free, end edges, doubles), where anchored rules are expected to be redundant, so it cannot refute H6 as a statement about the rule table alone. Rules-only sensitivity, re-run (`skeptic_sensitivity.py N 4 5`: `ALL_MOVES` vs `VERIFIED_MOVES` = floating only, with free=False, edges=False, doubles=False, limit 20000, so exactly the anchored rules are dropped). Core `45`, rows in orbit all vs floating: n = 12 DIFF at offsets 0 (1801 vs 2; almost-separate row reached vs not), 4 (91 vs 2), 5 (132 vs 1; verdict differs), SAME at 1, 2, 3; n = 13 DIFF at offsets 0 (3529 vs 2; verdict differs), 5 (177 vs 2), 6 (177 vs 1; verdict differs), SAME at 1..4. So in a rules-only walk the anchored rules do matter, at the head and tail and also at one interior offset (4 at n = 12, 5 at n = 13), and the redundancy found in the ablation comes from free/edge/double moves. "H6 is not even needed" and "the verdict is not coming from the rule table" are withdrawn as general statements; they hold for the reduced walk at these cores. I did not test H6 interior placements separately.
4. The width-8 length-14 probe was an unfiled inline command (132 applications, 20 s). It is excluded from the total: filed applications are 16 812 + 4 788 + 21 516 = 43 116 (the 43 248 earlier included the 132). Scope and claim are changed accordingly.

## Claim

(1) H1 survives. For each of the 330 floating rules of width 6..8, `verifyMove` (result equals the predicted LNA, every mutation of the sequence admissible, Coxeter polynomial unchanged) found 0 failures at one length each (not the first beyond the tests; see response): 16 812 applications (widths 6, 7 at length 12; this is 216 rules), 4 788 (width 8 at 13; 114 rules), 21 516 more (82 w = 6,7 rules at 13, killed by the 10-minute limit; the other 134 not done at 13), (one unfiled w = 8 rule at 14, 132 applications, is not counted): 43 116 in all. With round 049 (w <= 5 to length 11, 11 170 applications) every width except 9..11 has now been checked at one length >= w+5, with the gaps listed in Scope. This does not say the rules are length independent at all lengths or that widths 9..11 are; a failure would be a rule whose sequence is illegal or changes the Coxeter polynomial at some length >= 14 (or >= 13 for the 134 unfinished ones).

(2) H6 (anchored rules matter only at the head/tail), reduced walk only: dropping the anchored rules (floating table only) changes neither orbit size, orbit identity (hash of the row set) nor verdict (reaches an almost separate row or not) in any of the 45 placements, including the head and tail placements. So the anchored rules are redundant for these cores under the reduced walk that already has free moves, end edges and doubles, where they are expected to be; this cannot refute H6 as a statement about the table alone. In a rules-only walk dropping them changes the orbit and the verdict (see Response, item 3). Dropping the whole table (rules = []) changes the orbit at 2 of 45 placements (core `46`, n = 13, offsets 0 and 4: 842 rows instead of 2386, same verdict `inside`) and no verdict at any. So in the reduced walk at these cores the verdict is not coming from the rule table (the finding of F-032 and FINDINGS.md line 607, not new); F-051's "outside = no almost separate row in the closed orbit" is mostly a statement about free/edge/double moves. Falsifier for this reading: a core whose verdict changes between `ALL_MOVES` and `[]` (none in 45).

## Evidence

H1 (counts from the `DONE` lines in `theorist_rulelen_12_*.txt`, `_13_*.txt`, `_13b_*.txt`):

| width | rules | length | applications | failures |
|---|---|---|---|---|
| 6, 7 | 216 | 12 | 16 812 | 0 |
| 8 | 114 | 13 | 4 788 | 0 |
| 6, 7 | 82 of 216 | 13 | 21 516 | 0 |
| 8 | 1 (unfiled probe, excluded) | 14 | (132) | 0 |

Times: about 400 s per quarter at length 12 and 105-135 s per eighth of the 114 width-8 rules at 13 (4 processes); the w = 6,7 rules at 13 cost 20-40 s each, i.e. about 15 min on 4 cores for all 216 (would not fit one command).
H6 (`theorist_ablate_out.txt`, one line per offset: all / floating / no rules, then SAME or DIFF between the first two): cores (a,b) with relation lengths a,b adjacent: 45, 34, 44 at n = 12; 55, 35, 45, 46 at n = 13 (also 26 at 12); offsets where admissible. Placements 45; SAME 45/45 for all vs floating; all vs no-rules differ only for `46` at n = 13 offsets 0, 4. Orbit sizes for `45` at n = 13 reproduce round 049 (2386, 1127, 4217, 4217, 1127, 2386, 447). Many cores are in one orbit (e.g. 34, 44, 55, 35 all give 1410 or 447 rows), so these cores are weak tests: the only informative cores are `45`, `46` (outside placements exist).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 12 6,7,8 K 4   # K = 0..3, run 4 at once, ~400 s each
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 13 8 K 8       # K = 0..7, 4 at once, ~130 s
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 13 6,7 K 8     # K = 0..3: only ~20 rules each finish in 10 min
timeout 9m .venv/bin/python workshop/rounds/051/theorist_ablate.py 13 4 5            # ~3 min; other cores: `N A B [limit]`, I used limit 20000
```
(The saved `_12_*` outputs print progress each 10 rules; the script was later changed to print every rule. The saved 13b files are the killed runs: their last line is the last completed rule.) A width-8 rule at length 14 was a one-off timing probe in an inline command (132 applications, 20 s), not filed and excluded from the counts.

## Prior record

Round 049 (E-157): H1-H6 stated, H1 tested w <= 5. `verifyMove` is run by the tests for `lengthsToCheck(w)`: w+1..w+4 for w <= 6, w+1..w+2 for w = 7, 8. Grep of `research/` for rule length independence beyond that found nothing, so H1 for w = 6..8 beyond the tests is new but a null. Corrected after review: F-032 (FINDINGS.md line 1388) records that the table on top of the double mutation and free move adds nothing at n = 10 (262 rows, 16 orbits), and FINDINGS.md line 607 ("Which move does it matters, and the answer is not the rule table") records table alone / table + edges: no join at n = 10..16, table + edges + doubles: join. Added here: only the anchored rules dropped, orbit identity and verdict at 45 placements, n = 12, 13, plus the rules-only contrast.

## Code changed

None under `quivermutation/`. New: `workshop/rounds/051/theorist_rulelen.py`, `theorist_ablate.py`, outputs.

## Next

- OVERNIGHT: H1 for the 134 unfinished w = 6,7 rules at 13 (slices 4..7 of 8 and the tails of 0..3; ~10 min per slice alone) and length 14 for w = 8 (about 20 s per rule per process, 114 rules: ~10 min on 4 cores); widths 9..11 at their w+5 = 14..16 (needs ~1e6 LNAs at 15: probably infeasible by enumeration; propose random sampling of LNAs with the pattern embedded instead).
- Skeptic: re-run H6 on more informative cores (those with outside placements at n = 13, 14, from F-053's list), and check whether the rule table is ever needed for a verdict: the question H3 asks in reverse.
- Theorist (me): H4 (a map of orbits c@o and c@(o+1)); why the table is redundant for the verdict.
