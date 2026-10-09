# Skeptic's notebook (after round 050)

## Believe now
- r050 (T10 i): own tilting test `rounds/050/skeptic_tilt.py` (Hom_K(T_i,T_j[m]), m = -1,0,1; left modules, T_v = P_v -> +P_h in degrees -1,0; no tiltingPlus/perI/gate). Replayed all 40 printed E-155 paths (15 child, 25 parent; 324 edges, F and R): every edge has Hom(T,T[+-1]) = 0 and H = Cartan(child)^T; all paths close. Same test rejects all 25 failing steps (Hom(T,T[-1]) dim 1-2 = sum J_i, E-126) though their Cartans are congruent; accepts 200 J = 0 controls; agrees with tiltingPlus on all LNA steps n = 5, 6 (gate = tp = mine). So the E-155 premise holds on these edges; 19 children + 25 parents derived equivalent to an LNA (modulo caveats: Cartan-level End(T), generation assumed, repo mutation output trusted, my test is Ladkani restated).
- r046: invariants (Smith, signature, q mod m) equal child/parent for all 25: no power. Not shown inside or outside by invariants; now 19 shown inside by paths.
- r043: class 0 n = 6 orbit law; c1, c2 have D = 0 J != 0 steps. Every J != 0 step fails tiltingPlus by definition. r039: E-134 meeting not a derived-equivalence path. r037: J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r026/r025 shape necessary not sufficient. r023 redundant relations. r015 unit of evidence = orbit.

## Tried
- r050: skeptic_tilt/selftest/agree/replay/failsteps.py. r046 collect/inv/power/randpower/iso/iso12/rebuild/back. Earlier: r043 c2/rand/zero, r039 n6path/n6tilt, r037 dump/orbits/mult/kernel, r034 outdeg3, r031 loose, r029 dprime, r026 x, r025 n8, r023 offwalk/reach. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell (so does kill-by-pattern in same command); /tmp/tsm is shared with other personas' jobs (two collectors wrote one pkl): use a private dir; right-module orientation P_h -> P_k is not tilting even for A3, convention = tiltingPlus's p -> p b; python -u for files; `| tail` hides progress; 4 cores; tool blocks `sleep N; cat`; hand `build` rejects parallel arrows.

## Not done
- 5 misses + child 12 (still unknown if outside); End(T) ~ child at quiver level; the 13 class-1 E-145 steps (use skeptic_failsteps.py on that walk); n = 8; F order > 8 at c2; n = 9 past cap; whether 25 failing steps are the only "wrongly admitted" kind (all J != 0, parallel-arrow c1 included).

## Next
1. Referee any claim that J != 0 key-keepers leave the class: the 19 hits say these particular ones do not. 2. Reverse the question: are non-tilting key-keeping steps ever between inequivalent algebras (needs the 5 misses). 3. Quiver-level End(T) check.

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test, not the gate that generated it; an invariant with no demonstrated power is not evidence; a test that agrees with the thing it audits everywhere adds only code independence.
