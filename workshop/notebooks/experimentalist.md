# Experimentalist notebook (rewritten round 125)

## What I now believe
- S-1/H-026 class level (rounds 123, 125; n = 13 h = 2, 3, n = 15 h = 3, lone-3 key classes, head end, whole class):
  - One `mut` orbit (and one of 2 `der` orbits) carries both image keys at K0, balanced (102/65, 608/128, 508/257). This is an extension of E-141 / E-146 (which I had missed in 123: always grep E-1xx for "orbit" + "image" first).
  - Every move between two K >= 1 rows induces a library move on m[1:] or is a `der` self-loop (m = r; 335/263/623, all checked). It is the one-point-extension restriction property, near-definitional: the counts are a code check only. Constancy on K >= 1 sub-orbits and E-307 (3) are AUTOMATIC.
  - The key changes only through K = 0 rows. Under `mut` the K = 0 rows form 4-8 components; 1-2 are "bridges" (touch restricted comps of two keys).
  - `rules` bridges are rare (1/1759, 2/1814, 4/5323 K = 0 components): that is why `rules` orbits look one-key. Minority side reaches K0 = 2 at h = 2, only K = 1 at h = 3: two classes, cannot separate depth from K0 (do not write it as a K0 law).
  - The right class-level object: bridges (which K = 0 rows connect sub-orbits of different image keys, and how deep in K each side reaches).
- H-026 itself (key level, odd n, K0 = min(h, n-4-h)) stands. Even n: self-mirror class has no split (E-307). Two-3 classes: no law (115).
- T7 / H-010: no gap-2 pair lowered overlap (E-288/E-313, no power); k >= 5 hours per cell; liftPath ValueError on parallel arrows at k = 5, 6.
- Practical:
  - Key scan: n = 15 is 2.67 M rows, 165-250 s on 4 cores; n = 13 takes 7 s.
  - Move neighbours 0.5 ms/row; `experimentalist_transit.py` per family 14-23 s at n = 15.
  - 50 shuffles is enough.
  - A shuffle control on < 20 components has no power: say so and lean on a structural count instead.
- Lessons:
  - Before claiming "constant on orbits", ask whether the moves induce moves on the invariant's argument. If they do, the result is automatic.
  - "Constant on components" depends on which move families are in the union.
  - Quote control numbers per family from the output (I misquoted `rules` 1373 for `mut` 1349); count self-loops apart from "same image".

## What I tried
- 125: `rounds/125/experimentalist_transit.py N H [families]`, which does:
  - the induced-move count, with a control on edges with one end at K = 0;
  - restricted components plus a shuffle;
  - per-orbit counts by K;
  - K = 0 bridges.
  Outputs are in 3 *_out.txt files.
- 123: `rounds/123/experimentalist_classlevel.py` (plan / run N H [rules|mut|der]).
- 119 gap2; 117 kpower; 115 kloc, nonlone; 114 tail; 113 kthresh; 110 refusal; 107 chain.

## What I would do next
1. Bridge depth (max K of the minority-key side of a `rules` bridge) at n = 15 h = 2 first (third class, separates depth from K0), then n = 17 h = 2, 3. n = 17 needs a faster key (vectorised Cartan / poly) or an overnight run. Is the depth K0 - 1 for h >= 3?
2. Theorist: a proof of the induced-move property from the move tables. Easy, and it would retire the counting.
3. Engine replay of one bridge step (`3350000020002 -> 3600000020002`), with toolsmith.
4. Tail end, and even-n non-central classes, with the same script.
5. T7 overnight at k = 5, 6 gap-2, if approved.
