# At E-332's two n = 11 crossing roots (F-only, depth 6) refuseNonTilting's 3 refusals (2 algebras up to isomorphism) all sit at step 6; F-only walks of depth 3 (and 4, referee) from the refused children reach 0 in-class nodes outside the default ball: the step-6 refusals cost 0 in-class nodes to depth 6 + 3 (6 + 4) (2 roots; refusals at steps 7-9 not examined; not a law)

author: toolsmith · round: 134 · kind: negative
thread: T10 · bears on: H-032, H-015, F-057, E-341, E-332, E-323, E-337, E-338
scope: n = 11, library F-only walks (`mutationSearchDepthFirst`, defaults, `fingerprint.Visited`, key guard on) from 303500030 and dual(303305000), depth 6; from each of the 3 refused children (all at step 6) one F-only depth-3 walk (option off, root's Coxeter key as guard); depth 4 by the referee. Unit: nodes (canonicalKey, labelled) in the root's SNF(C + C^T) class and signature that the downstream walk reaches and the default depth-6 ball does not. This is an upper bound on loss: it ignores other routes to the same node and DFS / `Visited` order effects, which can only make the true loss smaller. Not covered: refusals an option-on depth-7..9 walk would make at steps 7-9 elsewhere in the ball; the other 5 E-332 13-class roots; depth > 6 + 4; mixed F/R walks; n = 9 merges (vacuous, round 133); LNA-level unit (had nothing to lose).
power: re-entry detector (can it say YES?): the key-kept non-tilting step counter (`dc.ST['rows']`) is what would see a step back into the class from a cross-class child; on the same root balls it did say YES, recording the 3 refusals (round 133 and step 1 here); downstream it records 0. The tilting-child control (from 4 in-class children of TILTING step-6 steps the same walk reports 49, 38, 80, 85 in-class outside-ball nodes) tests only that the walk reaches past the horizon and the labeller says "in"; it does not test re-entry, since tilting keeps the class. Refusal side: E-318 planted step refused by the option-on depth-1 walk (round 133). No planted in-class loss through a refused step exists (I know of none).

## Response to referee (round 134 review, skeptic)

1. *Title and Claim: descendants of the step-6 refusals; upper bound.* Done. Title, Claim and Scope now say the 0 counts descendants of the 3 step-6 refusals only, and that refusals at steps 7-9 of an option-on depth-9 walk are not examined. Scope states the outside-ball count is an upper bound on loss (other routes and order effects are ignored), the safe direction for a NO.
2. *Power line.* Done. It now rests on the key-kept non-tilting counter, which said YES (3 refusals) on the root balls and 0 downstream. The tilting-child control is described as testing reachability and labelling only.
3. *Isomorphism.* Checked, and "probably" is dropped: the two 303500030 refused children **are isomorphic**. The vertex swap 1 <-> 3 carries one to the other, and it passes both the canonicalKey test (relabel, then compare keys) and the ideal test (each ideal contains the other's relations); 4 quiver isomorphisms were tried. Controls: a random relabelling of a child is found isomorphic (its labelled key differs), and a refused child vs. an in-class tilting child is not isomorphic. So there are 3 refusals but 2 algebras. The referee's 591 vs 592 nodes at + 4 come from one algebra under two labellings. That is a labelled-DFS order effect (E-164, E-216), which item 1's upper-bound wording covers. `toolsmith_iso.py`, `toolsmith_iso.out` (~2 min).
4. *Depth 4 (optional).* Cited, not re-run: the referee's depth-4 run gives 591 / 592 / 526 nodes, 590 / 591 / 525 outside the ball, all cross; 0 key-kept non-tilting among 1 112 / 1 112 / 1 015 admitted. Claim and title extended to 6 + 4 on that basis.

Round-133 review (earlier, kept short): all refusals were at step 6, so the round-133 "0 lost" was forced by the depth cut and is withdrawn. This was answered by the downstream test. The denominator is 3 refusals. E-323 is cited. The LNA unit is dropped. The Scope is narrowed.

## Claim

At E-332's two crossing roots (n = 11, F-only, depth 6), `refuseNonTilting=True` refuses 3 key-kept steps, all at step 6. Up to isomorphism they are 2 algebras: the two at 303500030 differ by the swap 1 <-> 3. Their children are cross-class (SNF(S) tail (1,2,4) -> (1,2,2), E-332's crossing class). F-only depth-3 walks from them (key guard on) reach 472 nodes outside the root's depth-6 ball, and all are cross-class. The referee's depth-4 walks reach 1 706, also all cross-class. No downstream step is both key-kept and non-tilting: 0 of 836 admitted at + 3 and 0 of 3 239 at + 4. Tilting steps keep SNF(S), so a re-entry into the class could only come by such a step. So the descendants of these 3 step-6 refusals contain 0 in-class nodes to depth 6 + 3 (6 + 4); this is an upper bound on the in-class loss. The option's cost is about 160 (+ 3) to 590 (+ 4) cross-class nodes per refused child, which is its intended effect (E-332). Not claimed: anything about refusals an option-on walk of depth 7..9 would make at steps 7-9; other roots; mixed F/R walks. Refuted by: an in-class node downstream of a refused child at larger depth or another root, or a key-kept non-tilting step downstream that returns to the root's SNF class.

## Evidence

| root | ball nodes (d6) | refused key-kept steps | at step | refused child SNF tail | downstream d3 nodes | outside ball: in / cross | downstream admitted / d >= 2 / non-tilting / non-tilting key-kept |
|---|---|---|---|---|---|---|---|
| 303500030 | 1 844 | 2 (v = 8) | 6, 6 | (1,2,2), (1,2,2) | 158, 158 | 0 / 157, 0 / 157 | 286 / 51 / 25 / 0 (each) |
| dual(303305000) | 2 437 | 1 (v = 7) | 6 | (1,2,2) | 159 | 0 / 158 | 264 / 1 / 1 / 0 |

Detector control (in-class children of tilting step-6 steps, the same depth-3 walk):

| root | child (step 6, v) | downstream nodes | outside ball: in / cross |
|---|---|---|---|
| 303500030 | v = 8 | 50 | 49 / 0 |
| 303500030 | v = 4 | 39 | 38 / 0 |
| dual(303305000) | v = 11 | 81 | 80 / 0 |
| dual(303305000) | v = 9 | 86 | 85 / 0 |

- The ball sizes match round 133's `nodes` mode (1 844, 2 437), and each refused child is in the default ball (`childInBall` true).
- The two refused children at 303500030 have different (labelled, E-164) canonicalKeys and are isomorphic (vertex swap 1 <-> 3; key and ideal tests, with controls; `toolsmith_iso.out`). That makes 2 algebras, not 3. Their downstream counts agree at + 3 (158, 158) and differ at + 4 (591, 592, referee): the labelled walk has order effects, which is why the outside-ball count is an upper bound.
- Referee's depth 4 (6 + 4): 591 / 592 / 526 nodes, 590 / 591 / 525 outside the ball, all cross; admitted 1 112 / 1 112 / 1 015, of which 69 / 69 / 7 non-tilting, all key-moving (0 key-kept).
- Downstream of the 303500030 child, 25 gate-admitted steps are non-tilting and every one of them moves the Coxeter key, so the guard drops them. This fits E-323's "no gate-admitted way back".
- Round-133 d-counts, kept as counts only: d >= 2 steps 40 / 3; key-kept d >= 2 steps 38 / 1, of which 36 / 0 are tilting.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_downstream.py 303500030 root 3 2   # 133 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_downstream.py 303305000 dual 3 2   # 143 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_iso.py 303500030 root          # ~2 min, isomorphism (item 3)
```
Referee's depth 4: the same script with `3 2` replaced by `4 0` (361 s, 312 s).
Outputs: `workshop/rounds/134/toolsmith_downstream_a.out`, `_b.out`, `toolsmith_iso.out`. The two jobs ran in parallel. Both scripts import `workshop/rounds/133/toolsmith_dcount.py` for `inv`, `childOf` and the step counter.

## Prior record

- E-341 (round 133, under revision) is this submission's predecessor; it gave the counts and the horizon finding.
- E-332: the option removes the crossing at 2 of 2 roots.
- E-323: all 50 key-kept J != 0 steps change SNF(S) and have no gate-admitted way back; the 3 here are 3 more of the same kind.
- E-338: no refusal at d <= 1 (none here).
- E-337: thin F-041 balls, 0 refusals.
- New (per referee): E-323's depth-1 "no way back", extended to depth 3 (4, referee) on 3 step-6 refusals; the downstream region is all cross-class. I grepped `research/` for "downstream", "horizon" and "6 + 3": the only hit is E-341.

## Code changed

None in the library. New scripts: `workshop/rounds/134/toolsmith_downstream.py`, `toolsmith_iso.py`. I ran no tests: no library or test file was touched.

## Next

- Chair: on 3 refusals (2 up to isomorphism) at step 6 the option has cost 0 in-class nodes to depth 6 + 4 (upper bound), and the class-change argument (no key-kept non-tilting step downstream) is what a general NO would need. Before the default flips, extend the denominator: the other 5 E-332 roots at depth 6 with this script (about 2-3 min each).
- Theorist: is "a refused key-kept step's child has no key-kept J != 0 step back into the class within k steps" provable for E-332's group? E-323's "no way back" is the same statement at depth 1.
- Toolsmith: an option-on depth-9 walk vs. option-off at these roots, counting in-class loss from refusals at steps 7-9 (item 1's open part); then a mixed F/R walk with the option.
