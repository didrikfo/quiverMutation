# Scholar's notebook (after round 014)

## Believe
- Ladkani 1001.4765 2.3(c) = Aihara-Iyama 1009.3370 2.32(b) = `tiltingPlus`: one map (p |-> (p beta)_beta). Authors'
  statement, not independent of the code.
- r011: A5 (square a>b>d, a>c>d, d>e, relation abde = acde) is gate-admitted and fails tiltingPlus at d (n = 5..7).
- r014 (T5): guarded BFS from LNAs reaches gate-admitted, tiltingPlus-False parents at n = 6 (distance 8) and in 10/10
  smallest classes at n = 7, 8, 9 (distance 5..7). Every rejection is guard-refused (key moves); about 1.29e6
  guard-admitted steps, 0 fail tiltingPlus. n = 5 closes (11.7k algebras) with none. E-055/E-057 missed it by depth only.
- Reached parents are A5-shaped but not E-078 up to padding (extra prefix relation).
- Open oddity: n = 8 class 2, 10 steps gate True, tiltingPlus True, key moves, parallel arrows (not studied).
- CHZ Cor 3.6 "monomial?" still UNVERIFIED; arxiv.org blocked (403) in r006, r011. Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.

## Did
- R001/R002/R006: scholar_h015*.py, nonmono, step7, sides. R011: scholar_square.py.
- R014: rounds/014/scholar_walk.py (per key class BFS, --plan, --stop-on-reject), scholar_sweep.sh, scholar_replay.py,
  scholar_a5key.py; tests/test_gate_without_tilting.py.

## Next
- Chair: isTilting decision (redundant under guard); theorist: the 10 parallel-arrow mismatches and the one-map identity.
- Cheap check: Cartan congruence for those 10; and rejections' guard status at the remaining 15 n = 9 classes.
- Someone with the CHZ PDF: Cor 3.6 hypothesis.
- Lessons: "none found at depth d" is only a bound; split a walk by invariant class (Coxeter key) to size and
  reach deeper; classes with 2 starts are cheapest for first-reach tests. Check the depth of earlier negatives first.
