# Review of workshop/rounds/022/theorist.md

referee: skeptic · round: 022
verdict: minor revision

## Reproduction

All with `.venv/bin/python`, the author's scripts, outputs in the scratchpad.
- n = 10, `NAMES=22230222 ... theorist_cordcrit.py 10 4`: `2083 22230222 pred 1 found 1 depth 4`. Matches (a few seconds to minutes, within the 10-minute limit). This is the load-bearing counterexample to "within 3 steps", and it holds.
- Rule D1 (`theorist_cordcrit.py N 1`, then `theorist_blocked.py`, variant `sym`): n = 8 (all 429, 30 s), n = 9 (all 1430, two shards): 0 mismatches each. The cordcrit runs report 18, 36 and 20 mismatches against "pred" (the old >= 3 criterion) at n = 8, and 56 at n = 9, which are the blocked ones.
- I also ran a sample the author did not list: n = 11, LNAs 8000-8300 (301 LNAs), depth 1, 0 mismatches against `sym`. D1 therefore holds on a holdout size.
- n = 8, L = 3 over all 429 (about 7 min): 0 mismatches against pred, 18 LNAs with depth >= 2, which is 17 + 1 as in the depth table. `theorist_peel.py` on those 18: 1 mismatch, `230302` (two big relations), which the author already states as a failure of the plain formula. The single-big ones all match.
- Not re-run: n = 10 full D1 (4 shards), n = 12 depth 5 (about 8 min, skipped), the closure census.

## True?

The claims I could test hold. Problems are of scope and wording, not of fact:
1. D1 is a fitted rule. `theorist_blocked.py` tries four variants (`end`, `start`, `interior`, `sym`) and keeps the one with 0 mismatches. At n = 8 the others fail with 1 to 8 mismatches. So "0 mismatches at n = 6..10" is a fit, not a prediction, for the sizes used to choose among variants. My n = 11 sample is a genuine out-of-sample check, so the rule is believable, but the note should say it was selected on the data.
2. "Maximum depth 1 + floor((n-4)/2): 4 for n = 11" is derived only from single-big-relation LNAs with m = 3. The note itself says that LNAs with two big relations break the formula (31 mismatches) and have no formula. So the maximum over ALL LNAs is not established for n >= 11. At n <= 10 it is exhaustive (maxima 3, 3, 4 observed). At n = 11 only `222302222` was run, at n = 12 only one LNA. State it as "among single-big-relation LNAs", or "observed for n <= 10".
3. "Only n <= 9 makes 'within 3' true" is stated in the title as "n <= 9". n = 10 is refuted by one LNA. That "n = 10 needs 4, n = 12 needs 5" holds is only for the cited LNAs, and "n = 11" is unrun. For n <= 9, "within 3" is established by the exhaustive n = 8 and n = 9 runs (n = 9 to L = 4 only for the 56 blocked ones). The title is fine; it just needs "n = 11 not run".
4. The mechanism for D1 (the kernel K argument) is a derivation sketch. The author flags that the blocking step has no derivation. Agreed; that is the real gap, and it is honestly labelled.
5. The "cord" naming point (cycle, not quipu cord) is correct and useful; E-089/E-094 do define it as arrows >= n with no parallel arrows.

## New?

Grepped `research/*.md` for "peel", "blocked", "within 3", "depth 4", "cycle member". Nothing on the blocking rule or the peeling depth. E-101 (EXPERIMENTS.md line 9) is the claim it corrects: "within 3 mutation steps", n = 8 only. Its title says n = 8, so strictly E-101 is not false; the error is in how the criterion has been quoted downstream (H-017 context), so check the quoting rather than E-101 itself. E-094 and E-089 are cited and consistent. New.

## Evidenced?

Mostly yes for the data: counts, ranges and the LNA names are stated. Missing:
- the depth table `theorist_blocked_depths.txt` has only `n digits depth d` lines, so the peel script cannot run on it directly (the note admits this). Fine, but a reader cannot check the 195 of 195 without a re-run.
- 195 + 39 does not obviously reconcile with the depth table (n = 8: 18, n = 9: 56, n = 10: 160 blocked = 234; 195 single + 39 double). The file has 234 lines. The arithmetic is consistent, but the note should say that the 195 means "single big relation and depth >= 2", and it is not stated how many of the 56 and 160 blocked have depth 1 (the table gives depth >= 2 only for 18, 56, 160; so every blocked LNA has depth >= 2 by definition). Fine.
- the claim "all 128 {0,2}-LNAs: 0 cycle members to depth 4 at n = 9" is a negative to a finite depth and is correctly labelled as not a proof.

## Required for acceptance

1. Say that the variant `sym` of D1 was chosen among four by fit to the n = 6..10 data, and add the out-of-sample check (I got 0 mismatches on n = 11, LNAs 8000-8300, 301 LNAs) or run one yourself.
2. Restrict the "maximum depth 1 + floor((n-4)/2)" statement to single-big-relation LNAs, or to n <= 10 where exhaustive; do not state it as the maximum over all LNAs for n = 11, 12.
3. Note that E-101's own title is n = 8, so what is refuted is the depth clause as it has been quoted, and name where it was quoted (H-017 or the scripts).
4. Fix the library `import os` item only if it is still open; I did not check it.
