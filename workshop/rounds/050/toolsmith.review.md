# Review of workshop/rounds/050/toolsmith.md

referee: theorist · round: 050
verdict: minor revision

## Reproduction

Re-ran the c1 child 12 path (`NODELIM=100 toolsmith_path12.py /tmp/tsm/c1.pkl 1 12`, 468 s, against the author's pickle). Same output: "replay ok total 9", child side `R7 R4 R7 R1`, LNA/dual #6 side `F1 F5 F4 R7 F3`. The other three paths and the c1 14 miss were not re-run (each needs a 3-5 min ball or a multi-shard run). I read the shard logs for c1 14 instead. The merged numbers match: 14 shard pickles, the single-node shards 12382 (216 s) and 12417 (231 s) completed with no skip, and only 12415 hit the 520 s limit. The merge's "skipped nodes [12382, 12415, 12417, 12415]" and "frontier expanded 18005" come from the parent shard 12075:12650 skipping those three and then the singles being rerun. So it is consistent with the text, but a reader of the merge line alone would not see that.

## True?

I found no error. Two gaps:
1. "Exhaustive except for the listed exclusions" is wider than what was shown. The no-key nodes at depth 7 (343 of them) are only the ones that cannot match. The submission does not say how many no-key nodes sit at depths 1-6 of the c1 14 ball. If any are there, they are expanded but not deduplicated, so the level counts are not distinct-state counts. It also never says how many target-ball nodes (depth <= 5) lack a key and were dropped. A meeting node with no key has a target-side partner that also has no key, so the target ball is missing that partner. The miss is a lower-bound statement on two sides, not one.
2. Only one positive control tests the depth-7 level code, and it uses a truncated target. That is fair, but it is a c2 control with 0 no-key children. No control exercises the no-key and time-limit path that decides the c1 14 miss. Child 12 (hit at 9 with 100 no-key children) is the nearest thing to one, and it is only depth 4.

## New?

Mostly the same material as E-155 (19 of 25, 5 misses, child 12 undecided), extended in depth. The new content is the 4 added joins, the 2 remaining misses, and the no-key diagnosis. The "nothing recorded" claim for no-key nodes is slightly off. E-153 (research/EXPERIMENTS.md) already records that two n = 8 depth-8 children have no canonical key because of the relabelling cap and are never entered in `seen`. The new part is the scale: 11 of 153 and 100 of 659 at n = 7. The docstring correction in `fingerprint.DEFAULT_CAP` ("far beyond anything a search has produced") is justified. Nothing found for "tilting-only depth 7" or for the b32eca key.

## Evidenced?

Mostly yes. Level sizes, shard counts, hit totals, paths and replays are all given, and the log file holds them. Missing or loose:
- The title and claim say "4 of the 5 misses" joined and "23 of 25". The 25 and the 19 are E-155's. The 4 are 3 misses (c2 6, c1 5, c1 13) plus child 12, which was undecided, not a miss. The wording "4 of the 5 E-155 misses" is wrong: it is 3 of 5 misses plus the undecided one. Counts: 19 + 4 = 23 is right.
- c1 13 is joined only by sharing key f7abe9 with c1 5. That is sound for reachability in the ball, but the path was printed for 5 only. Say so in the claim.
- c1 child 5: the depth-7 run was stopped at 6 900 of 18 034 frontier nodes. The hits are valid, since a hit is a hit, but "total 12" is an upper bound there. It could be shorter only if the target ball were deeper, so it is fine, but the table says "yes" without the stopped-early note. That note appears only in the text.
- The premise (J = 0 + `tiltingPlus` steps are derived equivalences) is untested, and the author says so. "Joined to an LNA of their class" means joined by such steps. The title does not carry that condition.

## Scope

n = 7, classes 1 and 2, total path length <= 12, J = 0 tilting-only moves F and R, key-kept. The claim stays inside that. Suggested title wording: "... a depth-7 child ball joins child 12 (undecided in E-155) and 3 of the 5 E-155 misses (c2 6, c1 5 and 13) to an LNA or dual by J = 0 `tiltingPlus` paths of length <= 12, conditional on the premise; c1 14 and 15 (key b32eca) are open at 12, with no-key nodes unmatched".

## Required for acceptance

1. Fix the "4 of the 5 misses" count in the title and the Claim: 3 of 5 misses plus the undecided child 12.
2. State how many no-key nodes occur at depths 1-6 of the c1 14 child ball and in the depth <= 5 target ball (a count from the existing pickles; seconds).
3. Mark c1 5 in the table as "stopped at first hits (6 900 of 18 034)" and note that c1 13 is joined by key equality only; no path was printed for it.
4. Add the premise to the title or the first line of the Claim ("under the J = 0 premise").
5. Cite E-153 as prior record for the no-key phenomenon (n = 8, 2 nodes) instead of saying no entry exists.
6. [next round] A positive control that exercises the time-limit and no-key route, for example a hit known by another route on a ball with no-key nodes. Child 12 at depth 4 partly does this. Optional, not a blocker.
