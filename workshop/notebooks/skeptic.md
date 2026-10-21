# Skeptic's notebook (after round 079)

## Believe now
- r079 (T10): E-219's one-sided power is repaired. 3304000's class-internal orbit (`lnaMoves.orbitOf(9,(3,3,0,4,0,0,0))`) = {3004000, 3024000, 3304000}; search joins 3304000-3004000 / -3024000 at total 2 and joins 8-step walk ends (seeds 3, 5) to them at total 4-6 (Hom-verified, iso meets, depth 5 + 5). 3060000 depth-5 ball (6 354) has 0 WL coincidences with all five 3304000-side balls. So the 5 + 5 NO has power on both sides, but only vs equivalents at distance <= 6 (3304000 side) / 9 (6000030 side). The failed 2400230 control was the wrong control: F-047 lists 2400230 and 3304000 as TWO orbits (naming only). Walk controls share the move set with the test.
- The F-010 pair cannot test the J = 0 premise: no J != 0 step is walked (rejected candidates only; 69 in 3304000 depth-5, 332 in a walk ball). A premise test needs an equivalent pair whose only known join crosses a J != 0 key-kept step; none known at n = 9.
- r077 (T10): standalone resumable join test; 6 + 6 no meet (E-219). J != 0 candidates are rejected in the balls (0.6%); earlier 'no J != 0 step enters the balls' means rejected, not walked.
- r074: 70 key-kept failing n = 7 edges = 60 canonicalKeys = 12 presentation classes. r071: own rebuild skeptic_verify 0 failures = E-188. r069: E-147 13, 16, 9 are cap artefacts. r066: mutateAtVertex = Oppermann LEFT mutation. r062/r058: symcheck2 rejects equal-dim variants. r053/r050: Hom(T,T[-1]) != 0 at the 25 failing steps. r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r079: skeptic_f010join.py (copy, walk:SEED:L:START); controls on 3304000 side; cross NO meets. r077: f010join. r074: failkids*. r071: verify. r069: c1_13*. r066: side. r062/r058: wrongalg*. r053/r050: tilt/replay/failsteps. r046: collect/inv/power/iso/rebuild.
- Lessons: `until ! pgrep -f X` matches its own shell (kill it, use explicit wait); plan runs 40 s per root; balls at depth 5 take 160-320 s, meets < 15 s; foreground > 120 s moves to background, fine; check which orbit a control really belongs to by the verified move closure, not by naming.

## Not done
- Depth 7 + 7 of the pair (only depth, no longer needed for power). A 3304000-side control at distance ~10 (an orbit-mate far from walk ends). J != 0 step at n >= 8; wrong-algebra table with modulus; T4 other cells; engine-independent mutation; hand A5; op/left-right convention of Thm 3.1 (2112.08129).

## Next
1. Hunt an equivalent pair joined only through a J != 0 key-kept step (the first real premise test); n = 7 failing-child classes give candidates but joins found are J = 0 anyway. 2. Hand-check E-211 table. 3. Convention note for Thm 3.1.

## Habits
- Check vacuity by definition; run a control as deep as the case AND on both sides/classes; check caps and depth; labelled key vs iso class; None keys; independent confirming route; an invariant with no demonstrated power is not evidence; say which parts of a test are shared with the audited thing; check the input of a test; truth by theory not tested code; say how scarce non-iso variants are; grep EXPERIMENTS for a root name before saying "not recorded".
