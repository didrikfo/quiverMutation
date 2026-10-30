# Plan for cycle 001: Experimentalist

author: experimentalist · cycle: 001 (rounds 101-120) · not binding

## What I want to find out this cycle, and why

The workshop's empirical weak spots are controls and powers, not more
positive cells. T7 (H-010) has 8 of 8 two-bystander cells at overlap 2 and
0 lowering controls (E-271), so those cells say nothing yet. The Cartan-level
NOs of T3/T8 (E-265, E-269) have no power check on certified-equivalent pairs.
H-017's defect-per-key-class is a census on saved n = 9 data that nobody has
run. I would spend the cycle making each of these able to fail, then reading
the result. S-1's K-threshold law (K0 = 3, 4, 5 at n = 11, 13, 15) rests on
three lengths and needs a fourth or a reason.

## Phase 1, literature

I would not lead here; I would ask the scholar for sources that pin a number.
- Necessary condition "derived equivalent => Phi-equivariantly Z-congruent"
  (Nikulin 1.13.3 / 1.14.2 as cited in E-265, E-267, UNVERIFIED): checks the
  necessity step behind every Cartan-level NO. Search "Coxeter transformation
  derived equivalence Z-congruence Cartan form", Bondal-Polishchuk, Ladkani.
- Nakayama / gentle derived-equivalence papers with tables (Happel-Seidel,
  Bobinski-Geiss-Skowronski, Ladkani on Nakayama) giving explicit equivalent
  pairs at n = 9..12: independent positives for the power control, and counts
  to compare with H-020's rule table (E-225).
- Anything counting orbit sizes of mutation on A_n-type quivers with
  relations: checks |orbit B_n| = n 2^(n-5) (E-247, unproved).

## Phases 2-4, research

1. **T3/T8 power control (H-hypothesis to be written).** Run the E-269
   Z-congruence test on the 1 002 same-name n = 9 pairs of E-238 and the 6
   within-quipu n = 12 pairs. Step: a script over saved pairs, sized by
   `--plan`, a few minutes. One failure voids the 9 NO answers; zero failures
   gives them power for the first time. Either way a table.
2. **T7 lowering control (H-010).** Stepped chain starts 8, 10, 12 and other
   lengths, k = 3 margin 2, to find any case where bystander sharing lowers
   overlap. If none lowers across lengths, record H-010's cells as unfalsified
   at this shape and stop adding two-bystander cells; if one lowers, run the
   matching two-bystander cell. Margin 2 at k >= 4 goes to OVERNIGHT.
3. **H-017 defect per key class, n = 9 (saved data, E-203, E-221).** Group
   rows by key class, tabulate relations minus cords, report whether it is
   constant per class and how many classes were seen. Cheap; answers T6 yes
   or no. If time remains: S-1 at n = 17 only if a plan estimate is under the
   command limit, else the K-threshold law stays at three points.

## Phase 5, backlog

I expect left open: the 14 width-5 class-B rules of T4 (one-step rule), the
move-orbit count of the 18 328 lone-3 rows (S-1), and any overnight proposals
(E-204, E-253, E-260). I want those either sized and queued in `OVERNIGHT.md`
or retired, and the E-entries for the cycle's tables closed with commands.

## Progress

Progress: the power control run with a verdict (voided or surviving NOs); H-017
census done with a number; T7 either has a lowering control or a recorded
statement that none was found at three start gaps; one F-entry or status
change from these. Failure: more cells at n = 7 or two-bystander shapes with
no control; a timeout read as a verdict; a rerun of something in EXPERIMENTS.md.

## Disagreements

I disagree with spending more slots on join searches for T10 (E-250, E-255,
E-268 cannot test the J = 0 premise); I would hand that to the theorist and
skeptic for a non-move-set certificate and use my slots on the controls above.
