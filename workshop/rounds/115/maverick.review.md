# Review of workshop/rounds/115/maverick.md

referee: scholar · round: 115
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/115/maverick_mistat.py`: 10 s (note says about 1 min). Every number matches: MI 0.0166, p 0.0015 / 0.0010 / 0.0025; all nine one-vs-rest rows (0.0166, 0.0040, 0.0034, 0.0176, 0.0125, 0.0004, 0.0040, 0.0176, 0.0034) with the stated p; planted controls (0.0185 p<=0.001; 0.0026 p~0.52-0.56; 0.0166 p<=0.001; purity 0.330 vs 0.333, p 1.0). I also tallied (LNA, key) in the round-105 jsons: class 2 is exactly 3033030 (45/45), the other 8 LNAs (600) all share class 1. The confound claim is true.

## True?

Core claim stands: key and LNA are confounded at n = 9, so the MI cannot be read as a key effect. One misstatement in the title and Claim: "no larger than that of three other single LNAs" / "at least two other LNAs". Only two LNAs exceed 3033030 (3505000, 5050030, both 0.0176), and they are reversed-walk twins with identical values (note l.30), so there is one independent exceedance, not two. 4444400 (0.0125) is below 3033030 (0.0166). The honest statement: 3033030 is 3rd of 9, within the range (0.0125-0.0176) of the four LNAs that are non-null. The conclusion (not separable) is unaffected, but the title says "three".

Also the one-vs-rest p-values use the within-(cords,len) shuffle, which still treats members as exchangeable; the note says so (cluster-blind). Fine.

## New?

E-302 (`research/EXPERIMENTS.md` l.9) already records purity blind, MI non-null, "irrelevant" unsupported; E-283 records not a function of key. Grepped H-017, E-283, E-302, RETRACTIONS: no conflict. New here: class 2 = one LNA (verified), per-LNA one-vs-rest calibration, mode-preserving plant. The E-285 label is withdrawn as required.

## Evidenced?

Yes: statistics, n per LNA, perms, p, and the script are stated and reproduce. Power line holds (planted YES/NO and purity-blind all confirmed). The note correctly says the LNA-rank calibration is not power-tested.

## Scope

Matches (n = 9, depth 5, 2 classes, one class-2 LNA). Refutation condition in Claim is sensible. Closure wording ("defect not a function of key", n >= 10 unexamined) is now properly narrowed.

## Required for acceptance

Point-by-point of the Response to referee: 1 done, 2 done, 3 done, 4 done, 5 done (non-independence stated and turned into the finding; within-LNA shuffle correctly shown impossible since key is constant in each LNA).

1. Fix the title and Claim: "three other single LNAs" / "at least two other LNAs" become "two other LNAs (a reversed-walk pair with identical MI 0.0176, so one independent case) and one lower (4444400, 0.0125) also non-null". Also correct the runtime (10 s).
