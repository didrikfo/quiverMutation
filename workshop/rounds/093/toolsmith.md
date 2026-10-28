# 2400230(8) vs 3304000(6) at n = 9: no meet with one lnaMoves jump layer (jumps from LNA nodes only) on top of J = 0 F/R moves, implied by the component census and E-238; E-244 index 8427 is an isomorphic child of 8063 (MISS)

## Response to referee

Verdict minor revision; both required items done, no re-run needed (referee reproduced every line).
1. Required 1 (union MISS follows from the census): done. Claim now says the union MISS is implied by "0 common lnaMoves components" together with E-238 (a jump can only be taken from an LNA node, and the LNA components of the two balls are disjoint), not an independent test. The only new search content is the F/R meet at 8 + 6 (0 WL-bucket pairs), which repeats E-226/E-254 at larger radius. Title leads with this now.
2. Required 2 (A8 rerun): done. Reproduction now says the first A8 build exits 2 at the budget and leaves a depth-7 pickle (15 177 nodes); the build must be rerun to reach depth 8 (20 210 nodes). Using the depth-7 pickle, `lna` gives depths 7 6 and different counts (59 LNA nodes at depth 8, not 75).
3. Scope narrowing (referee): "no meet on the union" reads "no meet with one jump layer, jumps from LNA nodes only"; the 8427 point reads "isomorphic child of 8063, whose join at total 6 is E-244". Applied in title, Scope and Claim.
4. True? point 2 (2400230 isolated is about computed moves, `movesByRule` as complete undirected edge set, not about equivalence): accepted, stated in Claim.


author: toolsmith · round: 093 · kind: negative
thread: T10/T5 · bears on: F-010, F-047, E-226, E-238, E-244, E-254
scope: n = 9 · no meet means: one jump layer only, jumps taken from LNA nodes only; the 8427 point holds only as 'isomorphic child of 8063, whose join at total 6 is E-244' · E-212 `moves` balls (gate, J = 0 `tiltingPlus`, key kept): 2400230 to depth 8 (20 210), 3304000 to depth 6 (15 401) · one lnaMoves jump layer, over-approximated (any LNA node may jump to any member of its component at cost 1) · J != 0 steps only through the 8427 child · iso meets up to vertex relabelling (WL hash 3 iterations, then networkx iso)

## Stated before running

Sizing (`--plan`): both builds known from E-244 / E-254 (5 + 11 min, two commands each); census and meets under 1 min each. It fits in one sitting; no OVERNIGHT proposal needed.
- HIT would mean: a path 2400230 -F/R^a- LNA -lnaMoves- LNA -F/R^b- 3304000 with a <= 8, b <= 6 (or through the 8427 child), every edge Hom-verified (E-238) or J = 0 tilting. It would close E-226's 6 + 8 gap inside the union move set and show that the gap is a missing-edge-type effect, not a J = 0 premise effect. Consistent with the certified equivalence (E-238); no information on the premise.
- MISS means only: no such path of that shape. It cannot rule out paths with two or more jumps interleaved with F/R steps, F/R depth beyond 8 / 6, J != 0 steps other than index 8427 (E-254 scanned 2400230 only; 3304000 depth 6 is E-244, depth 7+ unscanned), non-LNA bridges the iso meet does not see, or anything about the J = 0 premise (2400230 ~ 3304000 is certified, so a NO is impossible as an inequivalence claim: a MISS is a statement about move sets).

## Claim

The union MISS follows from the component census (0 common lnaMoves components among the LNA nodes of the two balls) together with E-238 (the certificate connects them only through naming): a jump can only be taken from an LNA node. It is not an independent test; the only new search content is the F/R meet at 8 + 6, which repeats E-226/E-254 at the larger radius. 2400230 isolated is an observation about computed moves (`movesByRule`, undirected), not about equivalence.

At n = 9, with 2400230 and 3304000 in the union graph (J = 0 key-kept F/R edges plus lnaMoves jumps, jump cost 1 whatever the length), the radius-8 ball around 2400230 and the radius-6 ball around 3304000 (plus the child of E-244 index 8427 and its 13 J = 0 neighbours) have no common node up to isomorphism, no WL-hash collision at all, and no common lnaMoves component among their LNA nodes. Also: 2400230 is an isolated vertex of the lnaMoves graph (component size 1, all 1 430 rows scanned), so the jump layer cannot start at the root; and the child of index 8427 is isomorphic to the child of index 8063, which E-244 joined to 3304000 at total 6, so "8427 first" is not a new bridge candidate. It does not claim anything about the J = 0 premise or inequivalence (the pair is equivalent, E-238), and not that no longer path exists.

## Evidence

Builds reproduce the recorded sizes exactly (20 210, 15 401). Table, LNA nodes (rows standardise to a line; lnaMoves components over the full undirected graph on all 1 430 rows, edges = `movesByRule`):

| ball | LNA nodes by depth | LNA classes | components touched (full size in rows) |
|---|---|---|---|
| A = 2400230, depth 8 | 1,1,3,4,2,19,28,27,75 (160) | 9 | 3: {2400230} 1, {0400230 ...} 7, {2400330} 1 |
| B = 3304000, depth 6 | 1,1,4,4,4,22,32 (68) | 9 | 3: {3304000, 3004000, 3024000 ...} 7, {3304002} 1, {3024002} 1 |

Common components: 0. Control, same code, B against B: 3 common components, meet total 0 (95 449 WL pairs). F/R meet, A(8) against B(6): 0 WL-bucket pairs, 0 iso tests, so not even a hash collision (E-226 / E-254 form; baseline repeated). With the 8427 child (depth 7) and its radius-1 ball (depth 8, 13 nodes) added to B: still 0. Child of 8427: key kept, not an LNA, isomorphic (labelled quiver and relations, relabelling allowed) to the child of 8063, not labelled-equal (so the labelled `canonicalKey` would call them distinct: the E-164 effect, again). Both parents lie in ball B; each carries exactly one key-kept J != 0 candidate (R at 7, R at 9).
Reading of the 160 A-side LNA nodes: A's LNAs sit in its own 0400230 component (7 rows) and singletons; none of these rows is in 3304000's component, which the quipu certificate (E-238) says is only reached through naming, not through moves. Power: positive controls for this union form exist only inside 3304000's component (E-226) and B against B here; there is no control that crosses components, because no known pair does. The test is therefore as strong as the movable edge set, and that is the point it measures.

## Reproduction

All from the repository root; builds write to a scratch directory (pickles 20 MB and 17 MB, not committed):
```
DEPTH=8 timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py build 2400230 $S/A8.pkl --budget-hours 0.14   # exit 2 at 8 min (leaves a depth-7 pickle, 15 177 nodes; `lna` on it gives depths 7 6 and 59 LNA nodes at depth 8, not 75): the rerun of this same command is REQUIRED, 169-190 s, to reach 20 210 nodes
DEPTH=6 timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py build 3304000 $S/B6.pkl --budget-hours 0.14   # exit 2 at 8 min, rerun: 14 s (two builds in parallel: 4 cores)
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/093/toolsmith_unionjoin.py lna $S/A8.pkl $S/B6.pkl     # 6 s, toolsmith_lna_out.txt
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/093/toolsmith_unionjoin.py extra $S/B6.pkl               # < 1 min, toolsmith_extra_out.txt
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/093/toolsmith_unionjoin.py meet $S/A8.pkl $S/B6.pkl --extra   # 17 s, toolsmith_meet_out.txt (also without --extra, and B B control)
```

## Prior record

E-226 (6 + 8 no meet, J = 0 only), E-238 (lnaMoves components, certificate), E-244 (index 8427 listed "not run"; 087 `toolsmith_d6distinct_out.txt` already groups 8427 with 8063 and 11123 as isomorphic children; the 8063 and 11123 joins at total 6 were done), E-254 (2400230 depth 8 scan, proposed 8427 and "2400230(8) vs 3304000(6)"). New here: the union run itself, 2400230 being isolated in the lnaMoves graph, the component census of both balls, and the observation that the 8427 proposal is vacuous by isomorphism (the referee of E-254 already hinted this, review item 4). Not a rediscovery of anything in RETRACTIONS.md as far as grep shows.

## Code changed

New `workshop/rounds/093/toolsmith_unionjoin.py` only (execs the 087 and 082 scripts' definitions; no library change, no tests touched). Trap fixed on the way: that script must not use the names `A_` or `P` (087's prelude overwrites `A_`; `P` is the prime used by `rank`): a silent exit and a pow() TypeError respectively.

## Next

- Skeptic / theorist: since 2400230 has no lnaMoves edge at all, the gap in the move set is at 2400230 itself: which local move would connect 2400230 to its quipu (it names to P^(1,1,2)_(1,2) directly, E-238)? A `discoverMoves` run seeded with 2400230 would test it.
- Experimentalist: 3304000 depth 7 scan for J != 0 steps (about 4 x depth-6 cost) is the only unscanned exit on that side; 2400230 depth 9 (about 55 000 nodes) for A-side exits. Neither is a join, both only widen a MISS.
- Closing line for the thread: the F/R + jump join is a MISS at 8 + 6; with no cross-component control, it does not get stronger by being repeated at larger radius.
