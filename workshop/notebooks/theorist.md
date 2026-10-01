# Theorist notebook (rewritten round 015)

## What I believe now
- Round 015 (T5): the 10 "gate True, `tiltingPlus` True, key moves" steps of E-084 n = 8 class 2 are NOT a gap in `tiltingPlus`.
  Cause: `arrowPaths.reduceAgainstPivots` reduces only while the leading term is a pivot, and `idealBasis` is not reduced
  echelon, so residues are not normal forms; step 7 (`_kernelOverIdeal`) then misses a relation out of the mutated
  vertex (child too big by 1 at one Cartan entry). A full reduction (monkeypatch in `rounds/015/theorist_fix.py`)
  makes all 10 congruent with equal key; the 11 rejections stay rejections (cong fails, predicted entry -1).
  Not fixed in the library (toolsmith); not known how often it fires; E-084 BFS counts were made with the buggy rewrite.
- Cartan congruence agrees with `tiltingPlus` on all 11 replayed rejections (n = 7, 8, 9).
- Earlier (013): P (orbit of `3@2`/`3@3`) and Q (`5@0`/`6@0`) closed under offset shift 2 in a-1 double moves; shift by 1 unreachable, no invariant found.
  `4@0 -> 3@1` is one anchored width-4 move; `k@0` (k >= 5) has two neighbours. `5046`, `5056` at n = 13..17 have two parity orbits.
- Nulls (do not retry): no GF(2) functional, integer statistic or SNF separates P from Q; orbit stats show nothing.

## What I tried
- 015: scripts `theorist_{cartan,detail,step7,fix,child}.py`. Replay of recorded paths through `scholar_walk_n*_c*.txt` lines
  (format: depth, rels, vertex, [guard,] path with start index first); trace of step 7 with a wrapped `_kernelOverIdeal` was fastest.
- 013: `theorist_{path,k0moves,orbitstats,size}.py`; labelled BFS shows mechanisms best.

## Next
- Verify the library fix with a unit test, then ask whether `tiltingPlus`'s rank step (residues of p·beta, p in the ideal) is affected: construct a parent with a path i~>k in the ideal.
- Search for other places that treat `reduceAgainstPivots` output as canonical (`grep` in `quivermutation/`).
- Return to P/Q shift-by-1 invariant (mutation-count parity per vertex along the staircase) once the rewrite fix is in.
- Blind spots: I trust Ladkani 2.3(c) as iff and took End(T)'s Cartan to be R C R^T; not checked by building End(T).
