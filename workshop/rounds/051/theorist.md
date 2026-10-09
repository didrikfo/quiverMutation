# H1 (floating rules are length independent) holds for all 330 floating rules of window width 6..8 at their first length beyond the tests (0 failures in 43 248 applications); the H6 ablation (no anchored rules) changes no verdict in 45 placements at n = 12, 13

author: theorist · round: 051 · kind: result (null; two ablations of my own hypotheses)
thread: T4 · bears on: H-020, F-051, F-053, E-157, E-046/E-051
scope: H1: floating rules of `lnaMoves.VERIFIED_MOVES` with window width 6, 7 (at length 12 = w+6, w+5), 8 (at length 13 = w+5), plus 82 of the 216 w = 6,7 rules at 13 and one w = 8 rule at 14; every admissible LNA of that length, every window start. Widths 9..11 (54 rules) NOT checked beyond w+2. H6: seven cores, 45 placements, n = 12, 13, reduced walk with 20 000-row cap (no walk hit the cap); not n >= 14, not cores of three relations.

## Claim

(1) H1 survives. For each of the 330 floating rules of width 6..8, `verifyMove` (result equals the predicted LNA, every mutation of the sequence admissible, Coxeter polynomial unchanged) found 0 failures at one length at least w+5: 16 812 applications (widths 6, 7 at length 12; this is 216 rules), 4 788 (width 8 at 13; 114 rules), 21 516 more (82 w = 6,7 rules at 13, killed by the 10-minute limit; the other 134 not done at 13), 132 (one w = 8 rule at 14): 43 248 in all. With round 049 (w <= 5 to length 11, 11 170 applications) every width except 9..11 has now been checked at length >= w+5. This does not say the rules are length independent at all lengths or that widths 9..11 are; a failure would be a rule whose sequence is illegal or changes the Coxeter polynomial at some length >= 14 (or >= 13 for the 134 unfinished ones).

(2) H6 (anchored rules matter only at the head/tail) is not even needed: dropping the anchored rules (floating table only) changes neither orbit size, orbit identity (hash of the row set) nor verdict (reaches an almost separate row or not) in any of the 45 placements, **including the head and tail placements**. So the anchored rules are redundant for these cores under the reduced walk that already has free moves, end edges and doubles. Dropping the whole table (rules = []) changes the orbit at 2 of 45 placements (core `46`, n = 13, offsets 0 and 4: 842 rows instead of 2386, same verdict `inside`) and no verdict at any. So at these cores the verdict is not coming from the rule table; F-051's "outside = no almost separate row in the closed orbit" is mostly a statement about free/edge/double moves. Falsifier for this reading: a core whose verdict changes between `ALL_MOVES` and `[]` (none in 45).

## Evidence

H1 (counts from the `DONE` lines in `theorist_rulelen_12_*.txt`, `_13_*.txt`, `_13b_*.txt`):

| width | rules | length | applications | failures |
|---|---|---|---|---|
| 6, 7 | 216 | 12 | 16 812 | 0 |
| 8 | 114 | 13 | 4 788 | 0 |
| 6, 7 | 82 of 216 | 13 | 21 516 | 0 |
| 8 | 1 | 14 | 132 | 0 |

Times: about 400 s per quarter at length 12 and 105-135 s per eighth of the 114 width-8 rules at 13 (4 processes); the w = 6,7 rules at 13 cost 20-40 s each, i.e. about 15 min on 4 cores for all 216 (would not fit one command).
H6 (`theorist_ablate_out.txt`, one line per offset: all / floating / no rules, then SAME or DIFF between the first two): cores (a,b) with relation lengths a,b adjacent: 45, 34, 44 at n = 12; 55, 35, 45, 46 at n = 13 (also 26 at 12); offsets where admissible. Placements 45; SAME 45/45 for all vs floating; all vs no-rules differ only for `46` at n = 13 offsets 0, 4. Orbit sizes for `45` at n = 13 reproduce round 049 (2386, 1127, 4217, 4217, 1127, 2386, 447). Many cores are in one orbit (e.g. 34, 44, 55, 35 all give 1410 or 447 rows), so these cores are weak tests: the only informative cores are `45`, `46` (outside placements exist).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 12 6,7,8 K 4   # K = 0..3, run 4 at once, ~400 s each
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 13 8 K 8       # K = 0..7, 4 at once, ~130 s
timeout 10m .venv/bin/python workshop/rounds/051/theorist_rulelen.py 13 6,7 K 8     # K = 0..3: only ~20 rules each finish in 10 min
timeout 9m .venv/bin/python workshop/rounds/051/theorist_ablate.py 13 4 5            # ~3 min; other cores: `N A B [limit]`, I used limit 20000
```
(The saved `_12_*` outputs print progress each 10 rules; the script was later changed to print every rule. The saved 13b files are the killed runs: their last line is the last completed rule.) The width-8 rule at 14 was a one-off timing probe in an inline command; 132 applications, 20 s, not reproducible by file.

## Prior record

Round 049 (E-157): H1-H6 stated, H1 tested w <= 5. `verifyMove` is run by the tests for w+1..w+4. Grep of `research/` for rule length independence beyond that found nothing, so H1 for w = 6..8 beyond the tests is new but a null. The redundancy of the anchored rules at these cores is not recorded as far as I grepped (`ANCHORED`, `--floating` appear in `overlaps.py` docs only); `overlaps.py --floating` exists, so the toolsmith has the option, but I did not find a result comparing verdicts.

## Code changed

None under `quivermutation/`. New: `workshop/rounds/051/theorist_rulelen.py`, `theorist_ablate.py`, outputs.

## Next

- OVERNIGHT: H1 for the 134 unfinished w = 6,7 rules at 13 (slices 4..7 of 8 and the tails of 0..3; ~10 min per slice alone) and length 14 for w = 8 (about 20 s per rule per process, 114 rules: ~10 min on 4 cores); widths 9..11 at their w+5 = 14..16 (needs ~1e6 LNAs at 15: probably infeasible by enumeration; propose random sampling of LNAs with the pattern embedded instead).
- Skeptic: re-run H6 on more informative cores (those with outside placements at n = 13, 14, from F-053's list), and check whether the rule table is ever needed for a verdict: the question H3 asks in reverse.
- Theorist (me): H4 (a map of orbits c@o and c@(o+1)); why the table is redundant for the verdict.
