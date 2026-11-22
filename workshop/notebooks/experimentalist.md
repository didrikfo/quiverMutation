# Experimentalist notebook (rewritten round 135)

## What I now believe

S-1 / H-026, round 135 (`rounds/135/experimentalist_bridge.py`, `experimentalist_boundary.py`):
- E-325's transit claim held on 5 new lone-3 classes: n = 13 h = 1, 4 and n = 15 h = 2, 4, 5. Under rules, mut and der, 0 key-changing edges have both ends at K >= 1.
  - Power: the tail-deletion key gives 68-4972 such edges per run, and a random label gives thousands.
  - This is a code check of the definitional restriction property, not new mathematics.
- Boundary key changes (one end K = 0) are all `rules` moves directed K = 0 -> K = 1, onto a row starting 03, with m[4:] unchanged. The reverse direction gives 0 (the move generator is not symmetric; checked by `experimentalist_boundary_dir.py`, mut and der).
  - Counts: n = 13 h = 1..4 give 392 / 16 / 8 / 0; n = 15 h = 2..5 give 64 / 32 / 0 / 0.
  - At n = 13 h = 4 and n = 15 h = 4, 5 (the only h >= 4 classes checked), every key change sits on a K=0 to K=0 edge.
  - E-325's Power line already had the 8 / 16 / 32 counts. New in round 135: the shapes, the direction, h = 1, n = 15 h = 2 and the h >= 4 zeros.
  - Head shapes: h = 2 gives 3404 -> 035; h = 3 gives 3405 -> 036; h = 1 has 4 shapes.
- Bridge depth: minority keys reach at most K0 in every mut/der bridge. The claim that the majority key reaches the class max K is true for der bridges only. It FAILS for two small mut bridges: n = 13 h = 4, 23 rows, {0: 2, 1: 4}; n = 15 h = 4, 123 rows, {0: 4, 1: 4}. The referee caught this. Lesson: print every bridge's map, not only the largest one's. E-301 / E-304 are whole-class statements and do not imply anything bridge by bridge.
- `rules` bridges: none at h = 4, 5; several at h = 1, 2 (minority reaches K0).
- Cost: the n = 15 key scan takes about 170 s and is 90% of the run. n = 17 is about 8x that, so overnight or a faster key.
- Central classes n = 13 h = 4 and n = 15 h = 5:
  - 3 `mut` bridges each, plus one small bridge: 23 rows at n = 13 (starting 2027...) and 25 rows at n = 15 (2029...).

T7 / H-010, overlap 2, round 133 (still true, details in rounds/133):
- No clean interior (8,4)(10,3) cell lowers at k = 4: 15 decided, 9 capped at 560 s.
- Band move-closure: 0 of 5390 / 17738 clean cells lower at W = 13 / 14.
- A row lowers only if its longest overlap>=2 run has 3 relations and contains a gap-1 (3,3) pair. That is necessary, not sufficient.
- Double mutation = [-t,-t].

Practical:
- About 2 effective cores.
- Use Monitor or until-loops to wait.
- Write the submission early.
- Print vacuous counts beside every denominator.
- Power needs a control on the same side of the question. A planted-label control (tail key) works well for "edge-constancy" tests.

## What I tried
- 135: bridge and boundary censuses (above).
- 133: band / pick / reach / cox.
- 130: screen.
- 127: gap2fix.
- 125: transit.
- 123: classlevel.

## What I would do next
1. S-1: test whether "no boundary key change at h >= 4" holds at n = 17 (h = 4, 6, 7). Needs an overnight or faster key; propose for OVERNIGHT.md.
2. S-1: test whether the boundary counts scale by 2 per extra free tail vertex. 8 -> 32 and 16 -> 64 from n = 13 to n = 15 suggest it. Check at n = 14 for odd-h classes; note that the even-n self-mirror class fails H-026.
3. T7: the 9 capped two-bystander cells at a 1500 s cap; band closure at W = 15.
4. A faster Coxeter key (integer charpoly) would unlock n = 17 for S-1 and S_15 / S_17 for T7.
