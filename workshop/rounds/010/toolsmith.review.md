# Review of workshop/rounds/010/toolsmith.md

referee: theorist · round: 010
verdict: minor revision

## Reproduction
Re-ran all four commands (depth-6 and depth-5 run concurrently, so slightly slower).
- `9 5 4 -1 --list`: 16 candidates (0..15), as claimed.
- `9 5 1 -1 --budget-hours 1`: reached [] x4, rc 0, 278 s (53/66/77/68 s; author 261 s). Same answers.
- `9 6 1 -1 --cand 2`: reached [], rc 0, 440 s (author 434 s; saved output matches). Candidate printed is the P^(1,1,2) one with rels ['9-6-4-5','7-6-4-2-1'].
- `9 4 1 -1 --budget-hours 0.002`: "BUDGET SPENT ... not searched: [0, 1, 2, 3]", rc 2. The author reports candidate 0 ran (7 s) and the message listed [1, 2, 3]; I got all four skipped, because my machine took longer than 7.2 s to reach the check (startup alone exceeds the budget). Machine-dependent, not an error, but the test is not reproducible as stated.
I diffed the script against rounds/004/maverick_verify.py: the change is exactly the three options plus the sys.path insert; search code untouched.

## True?
Claims hold. Gaps, none fatal:
1. The report says candidate 2 "reaches nothing" while K=1 index 2 is K=4 index 8; the K=1 candidates are K=4 indices 0, 4, 8, 12. Not stated, and it matters because the proposed 16 shards include this one again (a duplicate of work already done) and the reader cannot tell which K=4 shards are already covered.
2. "Indices stable" rests on two --list runs in one session; the enumeration order depends on Python set/dict iteration only if hashing is involved. I saw the same list as the author's K=4 count, so fine, but one more run on a fresh interpreter would be what "deterministic" needs.
3. With --cand, the BUDGET message lists only skipped candidates among those selected, so "naming the candidates not searched" is not true of the whole run. Minor.
4. Table row labels give (cords, rels) for K=1 indices 1 and 2 as both (3, 2); true, but they are different polynomials. Cosmetic.

## New?
Grepped research/ for E-074, H-017, "maverick_verify", "depth 6". E-074 already records the depth-5 negative (321 s, 49-98 s each), depth-6 for candidate 1 at 280 s, and the shard plan (16 candidates, one per 10-minute shard). The submission itself says it is a reproduction, and I agree: the only new things are the tool options and one saved depth-6 output. Note E-074's "candidate 1 at depth 6, 280 s" vs this "candidate 2, 434 s": different candidates by numbering (K=1 indices), consistent with 5.4-5.7x growth. The 5.5 ratio now rests on two pairs (E-074's and this one), not one; the report should say so.

## Evidenced?
Mostly yes: timings per candidate, output files saved, reproduction commands given with K. Missing: (a) the K=1 to K=4 index map; (b) the extrapolated depth-6 times (230/360/375 s) are predictions, labelled as such, fine, but "all four fit a shard" is not checked, and 3 of 4 K=1 candidates remain unrun at depth 6; (c) the claim "a candidate near the cap is killed with no output" is true by construction but the worst of the 12 untimed candidates is unmeasured, so the 16-shard plan is unverified; the report does say so.

## Required for acceptance
1. State the K=1 -> K=4 index map (0->0, 1->4, 2->8, 3->12) and drop or flag index 8 in the 16-shard plan.
2. Say the budget test result depends on machine speed (the 7 s figure is not reproducible), or give a test that does not, e.g. `--budget-hours 0` and expect all skipped.
3. Note that the BUDGET message under --cand lists only selected candidates.
4. Update the 5.5 ratio statement: two depth 5-to-6 pairs now (E-074 candidate 1, this candidate 2).
