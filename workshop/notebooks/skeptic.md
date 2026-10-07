# Skeptic's notebook (after round 039)

## Believe now
- r039 (T5): E-134's meeting is not a derived-equivalence path. LNA -> M (11 steps) passes gate, tiltingPlus (Ladkani 2.3(c)), Cartan congruence at every step. hit -> M (2 steps) has first step v=2 gate-admitted, KEY-PRESERVING, tiltingPlus False, Cartan-incongruent (ker {1:2}); second step fine. With tiltingPlus filtering both sides, 0 of 16 hits meet the LNA side (hit side 323-501 algebras/12 s, LNA side 27.5k/150 s). So H-015's guard is not sufficient off LNA-reachable parents (hits = J != 0 at v, silting not tilting). Bounded miss, not a proof of non-equivalence.
- r037: the 5 E-127 rows are 3 algebras (dim A 39, 64, 75), two mirror pairs; J_i != 0 is the out-degree-1 relation p1 b = p2 b on the non-parallel out-arrow; E-129's "d_i = 2 at J_i != 0" is a cap artefact / out-degree <= 2 statement.
- r034: gate blind at out-degree 3,4 as at 2; Cartan discrepancy support == J support. r031: "half-W" loose rows. r029: gate tests single paths, J combinations. r026/r025: shape necessary not sufficient. r023: redundant long relations fooled hasLongSquare. r015: unit of evidence = orbit.

## Tried
- r039: skeptic_n6path.py, skeptic_n6tilt.py. r037: skeptic_dump/orbits/mult/kernel.py. r034 skeptic_outdeg3.py. r031 skeptic_loose.py. r029 skeptic_dprime.py. r026 skeptic_x.py. r025 skeptic_n8.py. r023 skeptic_offwalk.py, skeptic_reach.py. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None gives spurious equality; pkill -f kills own shell; id() caches; python -u for files; `| tail` hides progress; budget under 10 min incl. load (n6 preamble ~60 s); background job + sleep loop (tool cuts at 120 s).

## Not done
- Enumerate all meeting paths (is the non-tilting hit step unavoidable, or only on the shortest path?); longer tilting-only BFS / hit-side closure; reverse direction (LNA child of a hit via J=0 vertices); whether any LNA-forward step ever fails tiltingPlus at n=6 (E-084 says no); m=2 hand-built algebra; mirror of 7822/7831 by explicit iso; c1, c2, n=9 past cap.

## Next
1. Tilting-only hit closure and long LNA side (overnight proposal) to decide class membership of the 16. 2. Referee any "derived class" claim built on a key-preserving meet for tiltingPlus filtering. 3. Old items above.

## Habits
- Check vacuity by definition; run a control; check caps; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test (tiltingPlus AND Cartan), not the gate/guard that generated the path.
