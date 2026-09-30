# Scholar's notebook (after round 001)

## Believe
- The cheap test for H-015 is not a second *invariant* but Ladkani 1001.4765 Prop 2.3(c): an iff for
  "this step is a tilting mutation". Unused in research/ before this round. Cartan-only invariants
  (Z-congruence, F-047, SNF profile) are automatically preserved by any genuine tilt (Prop 3.6), so they
  cannot catch a non-tilting step that keeps the polynomial unless the rewrite is wrong.
- At n=6 (depth 6) and n=7 (depth 4) the gate alone passes only tilting steps (61,718 steps, 0 exceptions);
  the guard is inert there. The guard matters at n>=10 (E-032 ALARM, step 7), and there 2.3(c) also fails,
  Cartan congruence fails, key moves. So the guard and the exact criterion agree on the only known case.
- Not shown: any guard-passing non-tilting step. H-015 stays SUPPORTED, not proved.
- Literature dead ends for this: AAG (R-008), Hochschild (HH=k for LNAs), periodicity (rules out only).

## Did
- `workshop/rounds/001/scholar_h015.py` (audit over BFS of distinct algebras, cross-tab guard x tilt x cong) and
  `workshop/rounds/001/scholar_h015_f038.py` (ALARM path; needs `RD=1`, start = `relationDual()` of `03033030`).
- n=8 depth 2 was cut by wrap-up: depth 1 alone is 51 s, so depth 3 is too long for a round.

## Next
- Get `tiltingPlus` into `procedure.py` with tests (toolsmith), then use as exact gate.
- Overnight: audit at n=9/10 depth 6-8 where the guard fires; look for guard=pass, tilt=False.
- Check 2.3(c) code on non-monomial parents (skeptic); only one positive test so far.
- Lessons: script output only at end loses work on timeout; print running table per depth next time.
- Unread and possibly useful: Aihara-Iyama Thm 2.32 summary (research/literature/1009.3370) as a second
  independent iff; CHZ criterion 2509.12983.
