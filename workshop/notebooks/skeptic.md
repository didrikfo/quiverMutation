# Skeptic's notebook (after round 071)

## Believe now
- r071 (T4): own rebuild of verifyMove (skeptic_verify.py: enumeration, matching, application, Cartan/Coxeter polynomial; package engine only for mutation) gives w7 L10: 130 rules, 650 applications, 0 failures = E-188. Controls: bad sequence 560/650 flagged, bad result 650/650. Coxeter test never first to fire in a control (no isolated power demonstrated; separates 279/300 random pairs). Engine not independently rebuilt. Recommended T4 for closing (widths 6..8 done at every length; rest is filler).
- r069 (T10/T5): E-147 c1 numbers (13 D=0 of 67 distinct J!=0 steps) = a PREFIX of the E-151 BFS (expansions 22450..24032); E-151's 20000 cap gives 64 / 10 D=0 Cartan steps = 16 (parent,v) pairs. The 13 = 10 inside E-151's 16 + 3 beyond cap; all 13 fail tiltingPlus, key kept, D=0. Counts grow with cap (111 pairs / 89 steps at 45000, all failing). "13","16","9" are cap artefacts: quote the cap. Every failing claim is "under the J = 0 premise".
- r066: repo C = cartanMatrix has C[j][i] = #paths i~>j, C = H^T at T = A; mutateAtVertex = Oppermann's LEFT mutation (naming trap). Not done: hand A5 of Def 3.2, n >= 8.
- r062/r058: equal-dim wrong-algebra controls: symcheck2 does the work, dimcheck cannot reject equal-dim variants. r053/r050: Hom(T,T[-1]) != 0 at all 25 failing steps; E-181: End(T)=H_0 cannot test J. r046: Smith/signature/q mod m invariants have no power. r043, r037, r034 (gate blind at out-degree 3,4), r015 unit of evidence = orbit.

## Tried
- r071: skeptic_verify (+ _10_7.out). r069: skeptic_c1_13 (+ _rebuild, _prefix). r066: skeptic_side. r062: skeptic_wrongalg2. r058: skeptic_wrongalg. r053: replay3, key13. r050: tilt/selftest/agree/replay/failsteps. r046: collect/inv/power/randpower/iso/rebuild. Earlier as before.
- Bug lessons: pkill -f kills own shell (use kill PID); canonicalKey None == None; never `until pgrep -f name`; sympy det of many 7x7 too slow (Fraction elimination); checkpointed BFS segments of 520 s; ~40 expansions/s; fresh sessions may need `uv pip install -e '.[test]'`; python -u; LNA row has n-2 entries (off-by-one hit me r071).

## Not done
- T4: other four E-188 cells with my script; Coxeter-only sabotage control; engine-independent mutation. Class 2 prefix match; 45000 sample's 111 failing pairs not joined to LNAs; hand A5; wrong-algebra control on 370 path edges / 35 n=10 edges; n = 8; non-tilting key-keepers between inequivalent algebras.

## Next
1. Referee key-keepers-leave-class claims (cap quoted). 2. J != 0 step at n >= 8. 3. Class-2 prefix match. 4. Wrong-algebra table with a modulus. 5. If T4 survives conference, run skeptic_verify on w6 L11 etc.

## Habits
- Check vacuity by definition; run a control as deep as the case, and one that isolates each sub-check; check caps and depth; presentation dependence; pairs vs rows; None keys; confirm by an independent route; an invariant with no demonstrated power is not evidence; a test agreeing with the thing it audits adds only code independence (say which parts are shared); check what a test takes as INPUT; truth by theory not tested code; say how scarce non-iso variants are.
