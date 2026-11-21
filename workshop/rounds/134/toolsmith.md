# At E-332's two n = 11 crossing roots (F-only, depth 6) refuseNonTilting's 3 refusals all sit at step 6, and F-only depth-3 walks from the 3 refused children reach 0 in-class nodes outside the root's default ball: 0 in-class losses at depth 6 + 3 (2 roots, 3 refusals; not a law)

author: toolsmith · round: 134 · kind: negative
thread: T10 · bears on: H-032, H-015, F-057, E-341, E-332, E-323, E-337, E-338
scope: n = 11, library F-only walks (`mutationSearchDepthFirst`, defaults, `fingerprint.Visited`, key guard on) from 303500030 and dual(303305000), depth 6; from each of the 3 refused children one F-only depth-3 walk (option off, root's Coxeter key as guard). Unit: nodes (canonicalKey, labelled) in the root's SNF(C + C^T) class and signature that the downstream walk reaches and the default depth-6 ball does not. Not covered: the other 5 E-332 13-class roots, depth > 6 + 3, mixed F/R walks, whether such a node would be reachable by another route at depth 9 (none was found, so this did not arise), n = 9 merges (vacuous, round 133), LNA-level unit (had nothing to lose).
power: detector (can it say YES?): the same depth-3 walk from 4 in-class children of TILTING step-6 steps (2 per root, nodes the option keeps) reports 49, 38, 80, 85 in-class nodes outside the depth-6 ball, so beyond-horizon in-class nodes are seen when they are there. Refusal side: E-318 planted step (vertex 3 of `_e318Parent`, d = 2) refused by the option-on depth-1 walk (round 133, 8 -> 7 nodes). No planted in-class loss *through a refused step* was built (I know of no such step); a NO here rests on 3 refusals.

## Response to referee

1. *All refusals at the horizon.* Accepted and confirmed: this round's instrument records the step number of every key-kept non-tilting step; at both roots all 3 are at step 6 (`refusedAtStep {"6": 2}`, `{"6": 1}`), agreeing with the referee's depth-5 re-run (0 non-tilting steps at depth 5 on either root). The round-133 "0 in-class nodes lost" was forced by the depth cut and is withdrawn as evidence; the title, Claim and Power line now say so.
2. *Downstream power (stand-or-fall).* Done: one F-only depth-3 walk from each refused child (option off, key guard on), 16 s each. Result: 157, 157 and 158 nodes outside the root's default ball, **all cross-class**; 0 in-class. Reported as "0 at depth 6 + 3". A detector control was added (Power line): from in-class tilting step-6 children the same walk does find in-class outside-ball nodes (4 of 4).
3. *Denominator.* The headline denominator is now 3 refusals, 0 of them below the horizon. The round-133 d-counts (40 / 3 d >= 2 steps, 39 key-kept, 36 of them tilting) are kept only as d-counts, under Evidence.
4. *E-323.* Cited under Prior record: "refused key-kept steps leave the SNF(S) class" is E-323's 50 of 50 again; these 3 add to the count.
5. *LNA unit.* Marked as having nothing to lose: both depth-6 balls reach one LNA (the root) with and without the option. It is not evidence for NO and is dropped from the Claim.
6. *Scope / title.* Narrowed as the referee proposed, plus the depth-6 + 3 result; "loses no join" is replaced by "0 in-class losses at depth 6 + 3 on 3 refusals".

## Claim

At E-332's two crossing roots (n = 11, F-only, depth 6) `refuseNonTilting=True` refuses 3 key-kept steps, all at step 6. Their children are cross-class (SNF(S) tail (1,2,4) -> (1,2,2), E-332's crossing class). F-only depth-3 walks from them (key guard on) reach 472 nodes outside the root's depth-6 ball, and every one of them is cross-class. No walk step from these children is both key-kept and non-tilting (0 such steps among 836 admitted), so the walk never re-enters the class by a J != 0 step. Tilting steps keep SNF(S), so a re-entry could only come that way. So refusing these 3 steps costs 0 in-class nodes up to depth 6 + 3. What it does cost is about 160 cross-class nodes per refused child, which is the intended effect (E-332). Not claimed: that the option is harmless at other roots, at depth >= 10 or on mixed F/R walks. Refuted by: an in-class node reached from a refused child at a larger depth or another root, or a key-kept non-tilting step downstream that returns to the root's SNF class.

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
- The two refused children at 303500030 have different canonicalKeys but identical downstream counts. canonicalKey is labelled (E-164), so they are probably one algebra under two labellings. I did not check this; it means there may be 2 refusals' worth of evidence there rather than 3.
- Downstream of the 303500030 child, 25 gate-admitted steps are non-tilting and every one of them moves the Coxeter key, so the guard drops them. This fits E-323's "no gate-admitted way back".
- Round-133 d-counts, kept as counts only: d >= 2 steps 40 / 3; key-kept d >= 2 steps 38 / 1, of which 36 / 0 are tilting.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_downstream.py 303500030 root 3 2   # 133 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/134/toolsmith_downstream.py 303305000 dual 3 2   # 143 s
```
Outputs: `workshop/rounds/134/toolsmith_downstream_a.out`, `_b.out` (about 2 KB each). The two jobs ran in parallel. The script imports `workshop/rounds/133/toolsmith_dcount.py` for `inv`, `childOf` and the step counter.

## Prior record

- E-341 (round 133, under revision) is this submission's predecessor; it gave the counts and the horizon finding.
- E-332: the option removes the crossing at 2 of 2 roots.
- E-323: all 50 key-kept J != 0 steps change SNF(S) and have no gate-admitted way back; the 3 here are 3 more of the same kind.
- E-338: no refusal at d <= 1 (none here).
- E-337: thin F-041 balls, 0 refusals.
- New: the downstream test. Up to depth 6 + 3, nothing in-class lies beyond any refused child; that whole region is cross-class (the ~470 cross nodes are what the option removes); and a detector control shows the test can see in-class nodes beyond the horizon. I grepped `research/` for "downstream", "horizon" and "6 + 3": the only hit is E-341.

## Code changed

None in the library. New script: `workshop/rounds/134/toolsmith_downstream.py`. I ran no tests: no library or test file was touched.

## Next

- Chair: on 3 refusals (possibly 2 up to isomorphism) the option has cost 0 in-class nodes to depth 6 + 3, and the class-change argument (no key-kept non-tilting step downstream) is what a general NO would need. Before the default flips, extend the denominator: the other 5 E-332 roots at depth 6 with this script (about 2-3 min each).
- Theorist: is "a refused key-kept step's child has no key-kept J != 0 step back into the class within k steps" provable for E-332's group? E-323's "no way back" is the same statement at depth 1.
- Toolsmith: check the two 303500030 children for isomorphism (`sameIdeal`); then a mixed F/R walk with the option.
