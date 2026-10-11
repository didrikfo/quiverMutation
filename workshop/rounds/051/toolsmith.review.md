# Review of workshop/rounds/051/toolsmith.md

referee: skeptic · round: 051
verdict: accept (with the wording narrowing below)

## Reproduction

Not re-run in full: the ball builds take hours of wall clock (tpar, cball and 8 shards, 5+ min each). The claim is a hit, and a hit is checked by replaying the path, so I replayed it instead.
- `workshop/rounds/051/skeptic_replay13.py` (new) takes child 14 of `/tmp/tsm/c1.pkl` and applies `F1 F3 F1 F5 F4 R7 R1` to it. It applies `F7 R2 R1 F2 F7 R3` to `classes[base][9]`. Each step is tested with the E-161 Hom(T,T[m]) test from `rounds/050/skeptic_tilt.py`. That test does not use `tiltingPlus`, `perI` or the gate for the verdict.
- Result: all 13 edges have Hom(T,T[-1]) = Hom(T,T[1]) = 0, and the Cartan matrix of End(T) equals the child's (labelled). J = 0 and `tiltingPlus` hold on every edge. The two ends have the same non-None `canonicalKey`. This matches the printed "replay ok total 13".
- `skeptic_k1415.py` (new): children 14 and 15 are of the 16 c1 failures, both at parent depth 8, v = 7, and their `canonicalKey`s are equal and not None. So "15 joined by key equality" holds. `canonicalKey` equality means equal presentations up to gauge (docstring), not just similar invariants.

## True?

I found no error. The replay is independent of the search on the premise. It also uses fresh move generation: the author's own `follow` selects by (kind, v, key) from `moves()`, and mine builds each step directly.
- Nothing is wrong with the "not shortest" wording. 13 is the first total at which a hit appears with these two depths. 12 was missed at 6+6 and, by E-160, at 7+5.
- The 12 keyless nodes of cost >= 622 080 and the unfinished shard 11500:13800 can only hide further hits. The author is right that this does not matter for a hit.

## New?

Grepped `research/EXPERIMENTS.md` for b32eca, "children 14", "child 15". Only E-160 (open at 12) and E-157/E-161 come up. The total-13 path, the depth-6 target ball (62 297 keys) and the cap-5040 comparison are new. I have not checked FINDINGS or HYPOTHESES beyond what the title cites (H-015, T10).

## Evidenced?

Yes for the hit: the path is printed and I replayed it. The cap table and the no-key counts I did not re-derive. They are specific enough to believe, and they do not bear on the verdict.
- Weak controls: control 2 is a non-geodesic walk, and the "minimum exactly 13" control was not built. The author states both. The hit does not need a control.
- The "21 hit keys, all total 13" figure is consistent with the log, but only one path was printed and replayed. Child 15 and the other 20 meet keys rely on key equality or on the same code.

## Scope

The claim matches what was checked: n = 7, class 1, J = 0 `tiltingPlus` F/R steps. The "all 25 joined" in the title combines three things:
- the E-157 19, whose paths the skeptic replayed with the Hom test (E-161);
- the E-160 23 including its 3 new paths, which are still not Hom-tested;
- this new path, which I Hom-tested above.
The 3 E-160 paths (c2 6, c1 5, c1 12) have not been replayed with the independent test by anyone. Narrowed wording: "all 25 are joined by J = 0 `tiltingPlus` paths; the independent Hom test covers the E-157 paths and this one, not yet the three E-160 paths." The equivalence rests on the Okuyama-Rickard generation assumption (E-161 caveat).

## Required for acceptance

1. Put the narrowed wording above into the title or the Scope line (one sentence).
2. [next round, optional] Hom-test the three E-160 paths with `rounds/050/skeptic_replay.py`-style code. That closes T10 (i) at Cartan level for all 25.
