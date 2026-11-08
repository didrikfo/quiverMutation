# The n = 12 key group of 47 is not merged across its SNF split by the move census (9 orbits, 0 mixed), and key-guarded balls to depth 5 from three 14-roots (and depth 4 from two 33-roots) never leave their SNF(S) class; the walk has no natural crossing control (no power at walk level: the guard never fires in these balls)

author: maverick · round: 121 · kind: negative
thread: T10 · bears on: H-015, H-032, E-314, H-024
scope: n = 12, ONE key group (47 LNAs, key (1,1,-2,-4,-2,3,6,3,-2,-4,-2,1,1)). Census: `freeMoves.derivedOrbits(12, free, edges, doubles)` over the whole length (plus the restricted in-group version, also with free ADD of a 2-relation). Balls: the library's `mutationSearchDepthFirst` (gate + key guard, `fingerprint.Visited`), from the root and from its opposite; roots 3040500030 (depth 5 both sides; depth 6 right side only), 0040033030 (depth 5), 0500033030 (depth 4) in the 14-class; 3040400030, 0400033030 (depth 4) in the 33-class. Not covered: depth >= 6 on both sides, roots 3033050000/3050040030 (mirrors only), J != 0 vs J = 0 classification of the steps, other key groups.
power: detector: planted crossing nodes (the 33 LNAs of the other class fed in as ball nodes) flagged 33 of 33, the 14 same-class LNAs 0 of 14 (NO where NO is true, YES where YES is true). Census: SNF(S) is preserved on all 136 in-group census edges (certified derived equivalences), so the SNF labels are consistent with known YES steps. WALK LEVEL: none. The guard-off control (same root, coxeterGuard=False, depth 3 and 4, two roots) gives the IDENTICAL ball (438 / 1431 / 1552 visits, all key-kept): no gate-admitted step in reach moves the key, so there is no step known to leave the class within these balls, and nothing here shows the walk could have reached the 33-class if a key-kept crossing existed.

## Claim

(1) Census, answered NO, no bug: the move census (rule table + free move + edge moves + double mutation, `freeMoves.derivedOrbits`) partitions the 47 into 9 orbits of sizes 13, 13, 7 (all 33-class) and 4, 4, 2, 2, 1, 1 (all 14-class); no orbit is mixed, no census move from a group member leaves the key group, and the whole-length partition agrees with the restricted one. So neither "H-015 refuted" nor "census bug" applies: the census respects the SNF split. (2) Balls: in key-guarded balls (both orientations) of 1590-4658 node visits per root, every node reached (not only LNAs) has the root's SNF(C + C^T); no step lands in the other SNF class, let alone on one of its LNAs. The balls reach no LNA other than the root, its opposite, and their free-move variants, so "lands among the 33" at LNA level was never likely at this depth; the node-level SNF test is the stronger one. (3) Not claimed: support for H-015 beyond "no crossing in these balls". The guard is inert here (guard-off ball = guard-on ball), so this is a test of the gate, and without a natural crossing control it cannot distinguish "no key-kept crossing exists" from "none within depth 5". Refuted by: a guarded-walk node from a 14-root with SNF tail 12 (or the converse), or an orbit of `derivedOrbits(12, ...)` containing both classes.

Speculation level: tested on small cases (one group, depth <= 5).

## Evidence

SNF(S), S = C + C^T, computed twice per LNA (python-flint `fmpz_mat.snf` and pari `matsnf`, asserted equal): (0, 1^10, 20) on 14 LNAs, (0, 1^10, 12) on 33, reproducing E-314. The group is mirror-closed and the mirror preserves the label (14 -> 14, 33 -> 33).

Census (`maverick_census47.py`, CHECK=1): move applications from the 47 members: rule 2, edge 16, double 64, strip 24, add2 30; 136 in-group edges, 0 targets outside the group, 0 edges across the split. Whole-length `derivedOrbits`: 9 orbits meet the group, sizes [1, 1, 2, 2, 4, 4, 7, 13, 13] (sum 47, so no orbit holds a non-member), 0 mixed.

| orbit size | class | members (first) |
|---|---|---|
| 13 | 33 | 0400033030 0400233030 0400330030 ... |
| 13 | 33 | 3003304000 3003304002 3023304000 ... |
| 7 | 33 | 3004040030 3024040030 3034403030 ... |
| 4 | 14 | 0040033030 0240033030 2040033030 2240033030 |
| 4 | 14 | 3033040000 3033040002 3033040020 3033040022 |
| 2 | 14 | 0500033030 2500033030 |
| 2 | 14 | 3033050000 3033050002 |
| 1 | 14 | 3040500030 |
| 1 | 14 | 3050040030 |

The 14-class splits into 3 mirror pairs of orbits; the 33-class into one mirror pair (13 + 13) and one self-mirror orbit (7). Within each class the census does not join all orbits (lower bound: 6 and 3 classes under known moves).

Balls (`maverick_ball47.py`; "visits" counts routes, "distinct" the deduped nodes expanded, per side):

| root | class | depth | guard | visits | distinct (right/opp) | SNF tail of every node | LNAs reached |
|---|---|---|---|---|---|---|---|
| 3040500030 | 14 | 2 | on | 106 | 35/35 | 20 | root, opposite |
| 3040500030 | 14 | 4 | on | 1431 | 343/342 | 20 | root, opposite |
| 3040500030 | 14 | 5 | on | 4157 | 953/957 | 20 | root, opposite |
| 3040500030 | 14 | 6 (right only) | on | 5873 | 2659/- | 20 | root |
| 0040033030 | 14 | 5 | on | 4658 | 1140/1010 | 20 | root, opposite + its 3 free variants |
| 0500033030 | 14 | 4 | on | 1590 | 411/381 | 20 | root, opposite, 3033050002 |
| 3040400030 | 33 | 4 | on | 1597 | 375/397 | 12 | 4 LNAs, all 33-class |
| 0400033030 | 33 | 4 | on | 1735 | 435/420 | 12 | 6 LNAs, all 33-class |
| 3040500030 | 14 | 3 | OFF | 438 | 117/117 | 20, key kept at all | root, opposite |
| 3040500030 | 14 | 4 | OFF | 1431 | 343/342 | 20, key kept at all | root, opposite |
| 0040033030 | 14 | 4 | OFF | 1552 | 394/368 | 20, key kept at all | as guarded |

Key-kept nodes with a different SNF: 0 in every run. No cyclic (key None) node occurred. Growth about x3 per level (per side: 1.5 s, 45 s, 145 s at depths 2, 4, 5 after the 84 s `lnaKeys(12)` table).

Depth 6, right side only, root 3040500030: 5873 visits, 2659 distinct nodes expanded, every node SNF tail 20, no LNA but the root; 520 s wall-clock including the 84 s table (close to the 10-minute cap, so the opposite side at depth 6 was not run).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_census47.py              # 84 s
CHECK=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_census47.py      # 111 s (whole-length derivedOrbits)
CTRL=planted timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_ball47.py   # 85 s, planted detector control
ROOT=3040500030 DEPTH=5 timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_ball47.py   # 373 s
ROOT=0040033030 DEPTH=5 timeout 10m .venv/bin/python -W ignore workshop/rounds/121/maverick_ball47.py   # 390 s
ROOT=0500033030 | 3040400030 | 0400033030 DEPTH=4 ...                                          # about 180 s each
ROOT=3040500030 DEPTH=4 CTRL=off ...; ROOT=0040033030 DEPTH=4 CTRL=off ...                   # about 150 s each
ROOT=3040500030 DEPTH=6 SIDES=right ...                                                      # 520 s
```

## Prior record

E-314 (the split itself; the SNF certificate; group index 3 of rounds/107); E-291 (amended by E-314); E-306 (the n = 9 F-010 instrument: J = 0 steps cannot say NO; no crossing step found); E-147/E-151 (the guard is not a tilting test); E-308 (Ladkani necessity chain behind "SNF(S) differs => not derived equivalent"); H-032 (my proposal; this is its first case, YES-side consistency only). grep of research/ for the member names and "group of 47" finds only E-314: the census partition of the group and the balls are new. Not a rediscovery.

## Code changed

None in the library. New scripts `workshop/rounds/121/maverick_census47.py`, `workshop/rounds/121/maverick_ball47.py`. No tests touched.

## Next

- The finding that matters for T10: in this group the guard never fires (guard-off ball = guard-on ball to depth 4), as at n = 9 (E-306), so a key-kept crossing, if one exists, must be a gate-admitted step that the guard could not have stopped anyway. The question "can the guard be fooled" needs a ball where the guard does fire; this group is not one at depth <= 4.
- Walk-level power control (toolsmith): a planted crossing edge (graft one step to a 33-class algebra into the walk and check the detector and downstream bookkeeping), and a ball from a root whose unguarded walk does move the key, so that "SNF flags key-moving steps" is shown on real steps.
- H-032 (maverick, experimentalist): within each SNF class the census leaves 6 (14-class) and 3 (33-class) orbits; are they joined by deeper walks (meetingPoints between orbits of one class), or does the Cartan level fail to separate walk components (H-032's second falsifier)? Pairwise E-314 singular search between the orbits is cheap.
- Overnight (if wanted): depth 6-7, both sides, from 3040500030 and 0040033030 (about 8 min per side at depth 6, x3 per level), with SIDES split per command.
