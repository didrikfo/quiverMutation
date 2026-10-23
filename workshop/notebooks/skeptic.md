# Skeptic's notebook (after round 082)

## Believe now
- r082 (T10): E-229's key-kept `perI` step (3060000 ball, parent depth 4 via F1 R9 R8 R8, R at v = 8) is a REAL J != 0 step: `tiltingPlus` False (E-229's "also passes tiltingPlus-J" was a print label, never called), stepTest gives Hom(T,T[-1]) = 3 (S = -1) / Hom(T,T[1]) = 1 (S = 0), End dims = Cartan(child)^T anyway (Cartan data blind, as E-204). Child C is joined to 3060000 by Hom-verified J = 0 moves at total 7 (C: R1 R3 || root: F1 F7 F6 F5 F4) and to its parent at 6; no meet with the 3304000 ball at 5 + 5. So the premise survives this one event; it is not a premise test (join is J = 0 anyway). E-204's n = 9 sentence is from round-057 `walk` (pure F and pure R walks, quiverKey reach sets, depth 7), not the E-212 mixed-move balls: compatible with the kept step (code reading only; 057 walk not rerun).
- r079: E-219's one-sided power repaired; 5 + 5 NO has power vs equivalents at distance <= 6 (3304000 side) / 9 (6000030 side). No J != 0 step is walked in the F-010 balls; a premise test needs an equivalent pair whose only known join crosses a J != 0 key-kept step; none known at n = 9.
- r074: 70 key-kept failing n = 7 edges = 60 canonicalKeys = 12 presentation classes. r071: skeptic_verify 0 failures = E-188. r069: E-147 13, 16, 9 are cap artefacts. r066: mutateAtVertex = Oppermann LEFT mutation. r053/r050: Hom(T,T[-1]) != 0 at the 25 failing steps. r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r082: skeptic_kept.py (find/test), skeptic_join.py (copy of f010join with pkl:PATH:child|parent roots); balls C 13 982, R 6 354, P (depth 4) 5 242, Q 4 046. r079: f010join walks/controls. r077: f010join. r074: failkids*. r071: verify. r069: c1_13*. r066: side. r062/r058: wrongalg*. r053/r050: tilt/replay/failsteps. r046: collect/inv/power/iso/rebuild.
- Lessons: `pgrep -f X` in a wait loop matches its own shell (use ps|grep "[X]" or poll a log file); three balls in parallel make each ~2x slower; balls at depth 5 take 140-310 s alone, meets < 15 s; check which orbit a control belongs to by the verified move closure, not by naming; a "print label" in a script is not a test (read what is called).

## Not done
- Rate of key-kept J != 0 steps: depth 6 ball scan (about 20 min per ball, more with perI tests). Depth 7 + 7 of the F-010 pair. Rerun of round-057 walk to confirm the E-204 reconciliation; check whether E-212 mixed balls restricted to pure F/R reproduce 10 526. J != 0 step at n >= 8 of another type (here the long-arm vertices 4,5,6 carry the Hom); wrong-algebra table with modulus; T4 other cells; engine-independent mutation; hand A5; op/left-right convention of Thm 3.1 (2112.08129).

## Next
1. Premise test proper: an equivalent pair whose only joins go through a J != 0 step (search the E-212 balls for J != 0 kept children and ask whether C's J = 0 ball, minus the parent route, still meets; here it does, at 7 by a route not through P). 2. Hand-check E-211 table. 3. Convention note for Thm 3.1.

## Habits
- Check vacuity by definition; run a control as deep as the case AND on both sides/classes; check caps and depth; labelled key vs iso class; None keys; independent confirming route; an invariant with no demonstrated power is not evidence; say which parts of a test are shared with the audited thing; check the input of a test and what a script really calls; truth by theory not tested code; say how scarce non-iso variants are; grep EXPERIMENTS for a root name before saying "not recorded".
