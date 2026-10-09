# Skeptic's notebook (after round 046)

## Believe now
- r046 (T10 i): rebuilt the E-149 walk (`skeptic_collect.py`): n = 7 c1 16, c2 9 distinct key-keeping tiltingPlus-failers (all J != 0, parent depth 7-8; c1 8 have parallel arrows). Congruence invariants (Smith of xC+C^T x=-3..3, Smith of f(Phi), signature, q mod m = 2..7) are equal child/parent in all 25 (0 differ); explicit P in GL_7(Z) with P C_B P^T = C_A found for all 25 (box 1; one needs box 2), though not the step's R. Invariants are sound (random congruent pairs) and have power on random forms (23 of 100 equal-charpoly pairs separated), none on these. So: children are NOT shown outside the class and NOT shown inside. Bounded J = 0 reach to an LNA from children (6000 exp): 0 of 9 hits, but 0 of 4 for controls too (no power). Hand rebuild: c2 9/9, c1 8/8 buildable, all gate True, tiltingPlus False, J != 0, key kept.
- r043: independent D(x) build; class 0 n = 6, 7 steps all H1/H2 with x^2 coefficient 1; c1, c2 have D = 0 J != 0 steps (13, 9), E-140 "none" was a depth artefact. Orbit relation fails in c1, c2. Every J != 0 step fails tiltingPlus by definition so key-preservation there says nothing on derived equivalence.
- r039: E-134 meeting is not a derived-equivalence path. r037: J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r026/r025 shape necessary not sufficient. r023 redundant relations fooled hasLongSquare. r015 unit of evidence = orbit.

## Tried
- r046: skeptic_collect/inv/power/randpower/iso/iso12/rebuild/back.py. r043 skeptic_c2/rand/zero.py. r039 n6path/n6tilt. r037 dump/orbits/mult/kernel. r034 outdeg3. r031 loose. r029 dprime. r026 x. r025 n8. r023 offwalk/reach. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell; python -u for files; `| tail` hides progress; 4 cores, 3+ jobs slow each other; tool blocks `sleep N; cat`, use `timeout 118 bash -c 'until ...'`; hard-coded /tmp paths; control samples taken early are shallow (depth <= 4) vs failers at 7-8: match depth or say so; hand `build` rejects parallel arrows.

## Not done
- Whether any failing child is derived equivalent to its parent (needs HH^1 or a long tilting path or a proof of non-tilting-reachability); reverse tilting search with a deep positive control; matching the 13 + 9 E-145 steps one by one with E-149's 16 + 9; n = 8 c0, c1 (no failures at 1200 exp); c1 parallel-arrow hand rebuild; F order > 8 at c2 (s = 10); n = 9 past cap.

## Next
1. Referee any claim "failing key-keepers leave the class" or "key guard is safe": both unsupported. 2. Ask toolsmith for a reverse search with control at depth >= 7, or theorist whether the P found is induced by a derived autoequivalence. 3. Old items.

## Habits
- Check vacuity by definition; run a control and make it as deep as the case; check caps AND depth; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test, not the gate that generated it; an invariant with no demonstrated power (control also misses) is not evidence.
