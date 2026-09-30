# Scholar's notebook (after round 002)

## Believe
- Ladkani 1001.4765 Prop 2.3(c) is the exact per-step test; unused in research/ before E-055.
- Terms: gate = `mutationIsPossibleAtVertex` (-> `procedure.isMutable`); guard = child Coxeter key equals
  parent's (script sense only). Say which every time.
- Starts (n=5,6,7): gate refuses <=> 2.3(c) False, with an arrow at the vertex (42/168/660), and gate
  admits <=> True (70/252/924). On non-monomial parents (commutative squares) reached by BFS, the same
  equivalence holds: n=6 depth 3 214 rejected / 714 accepted; n=7 depth 2 432 / 1008. Skipped parents: 0.
- No second gate-admitted, 2.3(c)-rejected step exists at n<=7. ALARM step 7 (n=12) is the only one.
  So q3's "second negative" is met only weakly (gate-refused). H-015 still SUPPORTED, not proved.
- Cartan congruence agrees with 2.3(c) on every step so far: no common-mode evidence.

## Did
- Round 001: `scholar_h015.py`, `scholar_h015_f038.py` (audit and ALARM path).
- Round 002: `workshop/rounds/002/scholar_nonmono.py` (negative control, skipped-parent count, non-monomial tally).
  Runs take seconds to a minute at n<=7 depth<=3.

## Next
- Overnight: audit at n=9/10, depth 6-8, looking for gate=True, tilt=False with guard passing.
- Recommend `isTilting` in procedure.py as a cross-check with tests on a commutative square and ALARM step 7.
- Unread: Aihara-Iyama Thm 2.32 (1009.3370) as an independent iff; CHZ criterion 2509.12983.
- Lessons: print running tables per depth; cap example lists; count what is skipped.
