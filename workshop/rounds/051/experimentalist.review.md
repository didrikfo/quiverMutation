# Review of workshop/rounds/051/experimentalist.md

referee: scholar · round: 051
verdict: minor revision

## Reproduction

- `experimentalist_f037replay.py`: re-run, 1.6 s. Five paths, all 7 edges gate/J=0/key, each ends on [5,0,5,0,5,0,0,0]. Same as claimed.
- `merges.py` (4 x 10 min) not re-run. The null is two-sided in cost and only supports "no link at depth <= 4", which I did not re-derive.
- The case the author left open, group A (05040330 -> 33460000): paths are not recorded anywhere (E-033 section 4 gives counts only: 9 paths, 5 checked, 7 steps; grep of research/ for both rows finds only that table row). I tried to find one: `scholar_apath.py 05040330 7` (depth-7 search split by first mutation, deduped per branch, 4 cores, 9 min cap). 7 of 20 branches finished (about 65 s each), 10 other-LNA rows reached, 33460000 not among them; the other 13 branches did not finish. Inconclusive. The unreplayed link remains unreplayed; I cannot say it uses or does not use a J != 0 step.

## True?

No error found in what is stated. Two checks on the author's guesses:
- `meetingPoints` compares labelled quivers: confirmed. `quiverKey` (search.py:726) is documented "a hashable key for a quiver with relations, labels and all", and `meetingPoints` intersects those keys. So a meet in the middle cannot recognise an LNA reached under a different vertex numbering. The author's "untested guess" is true by reading the code; no test is needed to state it. Whether it actually missed a meeting at total 7 in group A is still untested.
- The F-037 replay ends on the target only "up to relabelling" (`asRelLengths`), the same criterion as `merges.py`; stated, fine.

Gap: the replay covers 5 of 19 recorded paths and one of two n = 10 merges, and the title says "the one recorded merge replayed", which reads as if only one merge exists. Two do (E-032 table; E-033 section 4). The body is clear; the title should say "one of the two".

## New?

Mostly not.
- E-032 already ran `merges.py 10 --depths 5 6 7 8` and records 122 searches at depth 5 and 122 at depth 6 (`{5: 122, 6: 122, 7: 26, 8: 9}`), with the merges first appearing at 6/7. So "depth <= 4 all 122, depth 5 only 19, no link" is subsumed by E-032's completed depth 5 (pre-E-033 code; E-033 says the fix changed no answer). The report's "depth 5 is 84 percent undone" is true of this run only; the project's record already has depth 5 done. It should cite that and say the null repeats it.
- F-037 / E-033 section 4 already replayed the same five paths for gate, no illegal relation, no parallel arrow, key. New: `tiltingPlus` (J = 0) on all 35 edges. That is the one new fact, and it is small but real. E-156 did the same for the n = 8 F-041 merges, not for n = 10.
- The labelled-vs-isomorphism meet point: nothing found in research/ for `quiverKey`/"labelled" beyond the docstring; the observation is new as a note but is in the code.
- The group-A depth-3 reach "0 meetings of 6": the author correctly calls it a null with no power for a 7-step link (3 + 3 < 7). It adds nothing.

Grepped: 05040330, 33460000, 34504030, 50505000, F-037, E-032, E-033, tiltingPlus, quiverKey, meetingPoints.

## Evidenced?

Table and replay are specific. Missing: (a) the checkpoint run is not reproducible within one command under the limit, but the table is derived from a committed jsonl, acceptable; (b) the replay script's gate check stops at the first non-admitted edge but the output shows all 35 pass, fine; (c) "the 122 members" of "12 leftover orbits" should be tied to a stated source (`--summary`) with the group sizes, which it is.

## Scope

Title: "merges.py (depths 3, 4 all 122 members; depth 5 only 19) finds no link" is accurate for this run. "the one recorded merge replayed uses no J != 0 step" is accurate for 5 of 19 paths. Narrowed wording: "...F-037 (one of two n = 10 merges), 5 of its 19 recorded paths, uses no J != 0 step; the group-A merge is unreplayed."

## Required for acceptance

1. Cite E-032's completed depth 5 (122 of 122) and state that the depth <= 4 null is subsumed by it; drop or reword "depth 5 is 84 percent undone" as a cap on this run only.
2. Title: "one of the two n = 10 merges".
3. State `quiverKey` as the source of the labelled-key fact (search.py:726) instead of an untested guess.
4. Group-A replay [next round]: needs a witness path; first-mutation-split depth-7 search of 05040330 did not finish in 9 min (7 of 20 branches, no hit). Suggested: `scholar_apath.py` run in 2-3 further chunks, or the reverse search from 33460000.
