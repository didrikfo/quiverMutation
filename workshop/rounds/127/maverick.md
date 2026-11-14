# No pure-F route to a J != 0 crossing within the reach tested: 20 F-only depth-5 balls at n = 12 (one root per SNF class of the 10 split groups) hold 0 key-kept J != 0 steps in 129 813, and no LNA reaches any of the 13 known crossing parents (E-318, E-323) by <= 6 (n = 11; 5 at three parents incl. n = 12) J = 0 F steps

author: maverick · round: 127 · kind: negative
thread: T10 · bears on: H-032, H-015, R-017, E-318, E-323, E-317, E-151, E-204
scope: n = 12, 20 roots (first-listed LNA of each SNF class of the 10 (SNF, signature)-split key groups, `toolsmith_snf12_out.txt`), F-only (library) key-guarded balls from root and dual, depth 5; backward: the 12 n = 11 and 1 n = 12 F-type crossing parents of round 125 (`toolsmith_jxclass_out.jsonl`, deduplicated up to iso; R-type parents dualised), R-only J = 0 balls, depth 6 (10 targets) or 5 (n = 11 #5, #10 and n = 12, 10-min cap). Not covered: depth >= 6 forward, other roots, crossing parents not yet found, pure-F routes containing an earlier J != 0 (non-tilting) F step, multi-vertex moves. Labelled dedup (ISO = 0).
power: backward instrument: the inverse control finds the start LNA from 11 of 11 sampled F-depth-K nodes (n = 11 K = 4, n = 12 K = 3), each first at depth exactly K, never earlier (truth: YES at K, NO below K, given the F-distance is K); with the J = 0 filter off, it finds 59-class LNAs at depth 2 from 3 of the 12 n = 11 parents through the known mirror crossing step, which the filter correctly excludes (they are not F-routes in reverse). Forward detector: at all 13 parents the same step code flags the known F crossing (13 of 13, truth YES). There is no control on which "no route exists" is known true beyond the depth-minimality above.

## Prediction (written to this file before any run)

P1: F-only depth-5 balls at n = 12, one root per SNF class: 0 key-kept J != 0 steps. P2: no LNA reaches a known F-type crossing parent by a pure-F walk of length <= 5. P3: the inverse control says YES on every sampled node. Reason given then: in F/R balls the first key-kept J != 0 step appears at depth 4 (E-318, E-323), in F-only balls at n = 7 only at parent depth 7-8 (E-151, E-204).

Outcome: P1, P2, P3 hold (P2 to depth 6 at 10 of 13 parents). Speculation level of the reading below: tested on small cases.

## Claim

(1) In 20 F-only key-guarded balls of depth 5 at n = 12 (58 214 nodes, 129 813 gate-admitted key-kept steps) every step is J = 0: no key-kept J != 0 step, so no crossing. With E-317 (the 47-group, depth 5) this covers one root per SNF class of every split n = 12 group, both sides.
(2) None of the 13 known F-type crossing parents (12 at n = 11, P12 at n = 12) is the end of a pure-F walk of J = 0 steps from any LNA or dual of length <= 6 (n = 11 except #5, #10: <= 5; n = 12: <= 5). The R-only J = 0 balls from these parents (1 230 to 7 807 nodes) contain no LNA at all.
(3) So the answer to T10's question at this reach is NO: every known crossing needs an R move (E-323 (4)), and no crossing of the library walk was found from new roots. This does not say pure-F crossings cannot exist: E-151/E-204 show pure-F walks do reach key-kept J != 0 steps at n = 7, but only at parent depth 7-8; the analogous depth at n = 12 is beyond a 10-minute command.

Not claimed: anything at depth >= 6 forward; that a pure-F route through an earlier harmless J != 0 F step is excluded (the backward inverse is only valid for J = 0 steps); anything about crossing parents not in round 125's list. Refuted by: a pure-F path, from an LNA or its dual, that ends in a key-kept J != 0 SNF-changing step, at the depths stated.

**Reframing (idea, not tested).** A walk of J = 0 F steps from L is a descending chain L = T_0 > T_1 > ... in the silting poset of L (each right mutation lowers, a derived equivalence preserves the order), so pure-F walks see only tilting objects below L; F/R walks reach objects incomparable with L. Read this way, (2) says: no known crossing parent is End(T) for a T at most 5-6 irreducible mutations below some LNA. A crossing on a pure-F walk needs a tilting T <= L whose next irreducible mutation is silting, not tilting, with key kept. If this is the right picture, the deep-first J != 0 of E-151 would be a statement about how far below L the first non-tilting irreducible mutation sits. I have no proof the repo's F step is exactly the left/right silting mutation of End(T); the control only shows it is inverted by R steps.

## Evidence

| set | count | depth | key-kept steps | key-kept J != 0 | crossings |
|---|---|---|---|---|---|
| fwd n = 12, 20 roots | 58 214 nodes | 5 | 129 813 | 0 | 0 |
| back n = 11, 12 parents, J0 = 1 | 1 230-7 807 nodes each | 6 (10), 5 (#5, #10) | 2 448-12 956 each | (not followed) | LNAs reached: 0 |
| back n = 12, P12, J0 = 1 | 1 691 nodes | 5 | 3 340 | (not followed) | LNAs reached: 0 |
| back n = 11, J0 = 0 | 12 parents | 2 | - | - | 3 parents reach a 59-class LNA (300330300, 300330302, 330340030) via the mirror crossing R step |
| ctrl | 5 + 6 nodes | K = 4 (n = 11), 3 (n = 12) | - | - | start found 11 / 11, first at depth K |

Per-root rows and per-target lines: `workshop/rounds/127/maverick_pureF_out.txt` (9 KB). Forward cost 134-418 s per root (4 jobs on 4 cores). The three J0 = 0 hits are no news: P's dual has the R form of the same E-318/E-323 crossing step, which then lands next to the 59-class LNAs (cf. E-323: X two J = 0 steps from 030330030); F-only depth-3 balls of these three LNAs have no J != 0 step (796 steps).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py fwd 12 5 0040033030      # 134-418 s per root; the 20 roots are listed in maverick_pureF_out.txt
J0=1 TGT=0 timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py back 6 11      # TGT = 0..11 at n = 11: 89-587 s; #5, #10 cap at depth 6
J0=1 TGT=0 timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py back 5 12      # 257 s (depth 6 caps)
J0=0 timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py back 2 11            # 31 s, filter-off check
J0=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py ctrl 003033030 4 5   # 101 s
J0=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py ctrl 0040033030 3 6  # 56 s
.venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py fwd 11 3 300330300,330340030,300330302  # 20 s
```

Needs python-flint (`uv pip install --python .venv/bin/python python-flint` if missing).

## Prior record

E-323 (4): every one of the 50 crossings has an R move; E-317's F-only depth-5 balls from 0040033030 and its dual are all J = 0 (consistent here: 4 656 steps from 0040033030 and dual, 0 J != 0). E-318: the F-only depth-6 ball of 003033030 is clean. E-151 / E-204: pure-F walks at n = 7 meet key-kept J != 0 steps at parent depth 7-8 (only n = 7, where no crossing can exist). E-321: recorded F-only merges n <= 10 have no J != 0 step. New here: the other 9 split n = 12 groups (both classes) at F-only depth 5, and the backward test from the crossing parents. grep `research/` for "silting order", "silting poset", "pure-F": the reframing is not recorded.

## Code changed

None in the library. New script `workshop/rounds/127/maverick_pureF.py` (reuses round-122/125 toolsmith code by exec), output summary `maverick_pureF_out.txt`. No tests touched.

## Next

- Toolsmith / overnight: F-only depth 6-7 at n = 12 from the 20 roots (depth 6 ~ 15-20 min per root; overnight-sized), and the backward J = 0 balls at depth 7-8 (meet in the middle: forward 4 + backward 4 would reach 8 for a fraction of the cost).
- Theorist: is the repo's single-vertex F step on End(T) the irreducible silting mutation of T in D(L) (sign and direction)? If so, pure-F reachability is an interval question in the silting poset and E-151's depth 7-8 may have an explanation.
- Chair: for the refusal-rule decision, the library's F-only walk has no observed crossing; the 50 crossings are all on F/R walks.
