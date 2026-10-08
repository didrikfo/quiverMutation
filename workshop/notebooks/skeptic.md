# Skeptic's notebook (after round 043)

## Believe now
- r043 (T5): independent build `skeptic_c2.py` (own det-poly D(x), own H1/H2, own orbit search). Class 0 at n = 6, 7: 316 / 300 distinct J != 0 steps, all H1/H2-ish, D low term x^2 coeff 1, s = 1 or -2 mod 8 (E-143 reproduces; no c_2 != 0 there). BUT n = 7 key-guarded walks (500 s) in class 1 (67 steps) and class 2 (64) contain D = 0 (key-preserving) J != 0 steps: 13 and 9, also under a tiltingPlus filter. So E-140 "none key-preserving" was a depth artefact; the law is class-0 only. Orbit relation fails (c1: 44 of 66, c2: 52 of 61 H1/H2 steps; random parents mostly); c2 D = 0 steps have s = 10. Random parents give c_2 != 0 matrices (x^2 coeff 2, 3, -1, -2) with key moving.
- Every J != 0 step fails tiltingPlus (by definition), so key-preservation there says nothing on derived equivalence; hand rebuild of 4 c2 steps: gate True, tiltingPlus False, Cartan incongruent.
- r039: E-134's meeting is not a derived-equivalence path (first step of hit -> M non-tilting, key-preserving). r037: E-127 rows = 3 algebras; J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r026/r025 shape necessary not sufficient. r023 redundant relations fooled hasLongSquare. r015 unit of evidence = orbit.

## Tried
- r043: skeptic_c2.py (walks, guard/tguard), skeptic_rand.py (random parents), skeptic_zero.py (hand rebuild of D = 0 steps; c1 pickle fails on parallel arrows). r039 skeptic_n6path/n6tilt.py. r037 dump/orbits/mult/kernel. r034 outdeg3. r031 loose. r029 dprime. r026 x. r025 n8. r023 offwalk/reach. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell; python -u for files; `| tail` hides progress; machine has 4 cores and 3+ jobs slow each other; first-N-found print only at N (print early); tool cuts waits at 120 s so use timeout 118 until-loops; hard-coded /tmp paths in scripts.

## Not done
- n = 8 c1, c2; whether the c1/c2 D = 0 parents are reached by a path whose steps all pass Cartan congruence (tguard says tiltingPlus yes, not replayed); c1 hand rebuild with parallel arrows; whether D = 0 children are derived equivalent (unknown, Cartan incongruent); F order > 8 at c2 (s = 10); m = 2 hand-built algebra; mirror of 7822/7831; n = 9 past cap.

## Next
1. Replay one c1 and one c2 D = 0 path from an LNA with tiltingPlus + Cartan congruence at every step, then test the child against the LNA side (E-137 style). 2. n = 8 c1, c2 sizing. 3. Referee any "key guard excludes J != 0" claim: demand class labels. 4. Old items above.

## Habits
- Check vacuity by definition; run a control; check caps AND depth (a sample of 2 steps is vacuous); presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test, not the gate/guard that generated it.
