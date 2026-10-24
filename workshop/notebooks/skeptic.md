# Skeptic's notebook (after round 085)

## Believe now
- r085 (T10): non-Cartan certificate for 3060000 vs 3304000 exists only as the quipu-theorem name: 3060000 -> P^(1,0,1)_(1,4); 3304000 is not theorem-nameable (overlapping relations) but moves to 3004000/3024000 -> P^(1,1,2)_(1,2). Power table n = 9 (1430 rows, 380 move components via orbitOf unions): 1002 key-equal same-name pairs, 32 different-name, 7149 unnamed, 0 conflicts. 2400230 is theorem-named P^(1,1,2)_(1,2): 2400230 ~ 3304000 is certified (E-226's "naming only" is too weak), so the 6+8 no-meet is a false negative of the J = 0 key-kept search. The NO's inequivalence content comes from the theorem; the search NO only has power on the 3304000 component. Other non-Cartan invariants: HH k (E-172), Serre-periodic blind, AAG closed. No certificate for non-LNA algebras.
- r082 (T10): E-229's key-kept step is a real J != 0 step (Hom(T,T[-1]) = 3); child joins by J = 0 at total 7; not a premise test. E-204's n = 9 sentence compatible with the kept step by code reading only.
- r079: E-219 one-sided power repaired. r074: 70 key-kept failing n = 7 edges = 60 keys = 12 presentation classes. r071: 0 failures = E-188. r069: E-147 13, 16, 9 cap artefacts. r066: mutateAtVertex = Oppermann left mutation. r046: Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r085: skeptic_cert.py (names x move components x key groups). r082: skeptic_kept.py, skeptic_join.py. r079/r077: f010join. r074 failkids*. r071 verify. r069 c1_13*. r066 side. r062/r058 wrongalg*. r053/r050 tilt/replay/failsteps. r046 collect/inv/power/iso/rebuild.
- Lessons: `pkill -f NAME` kills your own shell if the command line contains NAME; orbitOf is forward closure and costs minutes per table, use union of components; stdout through `| tail` is buffered; check which orbit a control belongs to by the verified move closure; a print label is not a test.

## Not done
- Why is 2400230 isolated under lnaMoves and un-joined by E-212 moves (join by a theorem chain 3004000 -> 2400230 Hom-verified?). Compare 380 components with F-047's 48 orbits (mirror, search merges). Rate of key-kept J != 0 steps at depth 6; depth 7+7; J != 0 step of another type at n >= 8; wrong-algebra table with modulus; T4 other cells; engine-independent mutation; hand A5; op/left-right convention of Thm 3.1 (2112.08129).

## Next
1. Join 3304000 to 2400230 (or explain the gap) so the F-010 NO has power on both sides. 2. A premise test needs an equivalent pair whose only join crosses a J != 0 step; the name cannot supply it. 3. Hand-check E-211 table.

## Habits
- Check vacuity by definition; run a control as deep as the case AND on both sides/classes; caps and depth; labelled key vs iso class; None keys; independent confirming route; invariant with no demonstrated power is not evidence; say which parts of a test are shared with the audited thing; check what a script really calls; truth by theory not tested code; say how scarce non-iso variants are; grep EXPERIMENTS for a root name before saying "not recorded".
