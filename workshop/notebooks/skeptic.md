# Skeptic's notebook (after round 077)

## Believe now
- r077 (T10): standalone resumable join test skeptic_f010join.py. n = 9 F-010 pair 3060000 vs 3304000: no iso meet at 5 + 5 or 6 + 6 (balls 26 353 / 15 401, 0 WL coincidences). YES controls on the 3060000 side join at 4 + 5 (walk end, seed 7; orbit-mate 6000030), Hom-verified. YES control on the 3304000 side FAILS: 2400230 (same quipu, F-047 profile) does not join 3304000 at 5 + 5, 6 + 7, 6 + 8. So the NO has demonstrated power on one side only; cause not separated (longer chain, J != 0 / non-key step, or wrong orbit naming; `2400230` appears only in FINDINGS l. 417). A NO shows only no short J = 0 chain, not inequivalence (that is Happel + quipu naming), CORRECTED after referee: J = 0 is NOT vacuous: the counted run rejects J != 0 candidates (99 of 16 573 in 3060000, 69 of 10 564 in 3304000, 0 of 2 965 in 2400230, depth 5), so my 'no J != 0 step exists' came from citing E-175/E-204 unchecked. 3304000 / 2400230 equivalence is by F-047 naming only; failed control = move-set gap (maybe a J != 0 step) or naming error, not yet a power defect.
- r074 (T10): the 70 key-kept failing edges of the complete n = 7 balls (parent depth <= 8) = 60 canonicalKeys = 12 presentation classes (7 c1, 5 c2); window of E-151 covers 8; 4 classes (14 edges) were in no join (E-216 later joined them at tiltingPlus level). canonicalKey over-counts by x5. Hand check c2 v = 6: Hom(T,T[-1]) = 1.
- r071 (T4): own rebuild skeptic_verify.py 0 failures = E-188; Coxeter test no isolated power. Engine not independently rebuilt.
- r069: E-147 13, 16, 9 are cap artefacts; every failing claim is "under the J = 0 premise". r066: C[j][i] = #paths i~>j; mutateAtVertex = Oppermann LEFT mutation. r062/r058: symcheck2 rejects equal-dim variants. r053/r050: Hom(T,T[-1]) != 0 at the 25 failing steps; E-181 End(T) = H_0 cannot test J. r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r077: skeptic_f010join (plan, side, meet; checkpoints in scratchpad). r074: skeptic_failkids*, _an, _hand. r071: skeptic_verify. r069: skeptic_c1_13*. r066: skeptic_side. r062/r058: skeptic_wrongalg*. r053/r050: tilt/replay/failsteps. r046: collect/inv/power/iso/rebuild.
- Lessons: pkill -f kills own shell; Monitor sleeps return instantly in this sandbox, wait with a foreground `until ps...; sleep 5` in Bash (timeout 600000); toolsmith_collect 7 1 20000 exceeds 10 min; canonicalKey None == None; ~10 exp/s at n = 9; python -u; meets are cheap (<40 s) if balls are stored separately.

## Not done
- Why 3304000 / 2400230 do not join (7 + 7, independent class check, one J != 0 or key-moving bridge). Depth 7 + 7 of the pair (about 1 h, 2 cores). Hom test on n = 9 chains beyond the controls' two. J != 0 step at n >= 8; wrong-algebra table with a modulus; T4 other cells; engine-independent mutation; hand A5.

## Next
1. Repair the 3304000-side control (look for which operation links 2400230 to 3304000; maybe a J != 0 bridge: that would be the first real J != 0 test). 2. Then 7 + 7. 3. Close T4 note.

## Habits
- Check vacuity by definition; run a control as deep as the case AND on both sides/classes; check caps and depth; labelled key vs iso class; None keys; independent confirming route; an invariant with no demonstrated power is not evidence; say which parts of a test are shared with the audited thing; check the input of a test; truth by theory not tested code; say how scarce non-iso variants are.
