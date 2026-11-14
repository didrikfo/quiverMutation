# Experimentalist notebook (rewritten round 127)

## What I now believe
- T7 / H-010 (round 127):
  - E-313's ValueError was a bug: `lnaMoves._copy` dropped `arrowRels`/`relCombinations`. Fixed, with a regression test in `tests/test_lna_moves.py`. Neutral where the old code did not raise (k = 4 end run identical).
  - Interior, margin 1, (8,l)(10,l'): still only 3 non-vacuous cells (k = 4: (8,4)(10,3), (8,5)(10,4); k = 5: (8,4)(10,3)), none lowers. (8,3)(10,3) is vacuous at k = 4, 5, 6.
  - Capped (undecided): (8,5)(10,5) k = 4 (280 s); (8,4)(10,4), (8,5)(10,4), (8,5)(10,5) at k = 5 (585 s).
  - Power control that works: plant the same shape at start 1 (`--end`). Lowering in 4 of 5 shapes at k = 4 and at k = 5. End runs are much faster (L = 13-15 against 20-22).
  - Lowering sequences at the end use vertex 1 or vertices 5-6 (4-5 from the end). So "end in reach" is broader than "mutate at 1".
- S-1/H-026 (rounds 123, 125). The class-level object is the K = 0 bridge.
  - One `mut` orbit carries both image keys at K0. The key changes only through K = 0 rows.
  - Constancy on K >= 1 is automatic (induced moves).
  - `rules` bridges are rare.
  - H-026 itself holds at key level for odd n. Even n fails at the self-mirror class (E-307).
- Practical:
  - 4 parallel jobs give about 1.7 cores of throughput on this box for `localMutationSequences`. Fewer jobs per cell is better for capped cells.
  - k = 4 -> 5 grew the cost 3-5x.
  - Key scan: n = 15 takes 165-250 s; n = 13 takes 7 s.
  - 50 shuffles is enough.
- Lessons:
  - "Documented design" in an E-entry can be a bug. Read the traceback before accepting it: here the bug was 3 frames up from the raise.
  - A placement control (the end) is NOT power for an interior "never lowers" claim (127 referee): the end changes what is reachable 4-5 vertices in, and the end mechanism is already known to work. Power needs a control on the same side of the question: an interior cell known to lower.
  - Dedup keys that conflate states bias toward NO. Check the key whenever a fix lets the search into new kinds of states (parallel arrows, `_stateKey`, F-039).
  - Ask whether moves induce moves on the argument of an invariant. If they do, constancy is automatic.

## What I tried
- 127: `experimentalist_nokey.py` (faithful dedup key on parallel states). `rounds/127/experimentalist_gap2fix.py STEPS MARGIN CAP cells [--end] [--jobs J]` (per-cell cap, errors recorded). Fixed `_copy`.
- 125: `experimentalist_transit.py` (induced moves, bridges).
- 123: `experimentalist_classlevel.py`.
- 119: gap2 (superseded by gap2fix).
- 117: kpower. 115: kloc, nonlone. 114: tail. 113: kthresh. 110: refusal. 107: chain.

## What I would do next
1. T7: overnight on the capped interior cells withdrawn (no power). First build an interior control: sweep of how far the start can be from vertex 1 and still lower (start 2, 3, 4 with margin 1): this measures the reach of the end.
2. S-1 bridge depth at n = 15 h = 2, then n = 17 (needs a faster key or an overnight run).
3. Look for other hand-built `PathAlgebra` copies that drop `arrowRels` (`quipuRelations.py` l. 591?).
