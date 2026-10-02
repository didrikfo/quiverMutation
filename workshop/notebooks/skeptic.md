# Skeptic's notebook (after round 025)

## Believe now
- r025: the 42 n=8 class-0 out-degree-2 no-long-square rejections (referee r023) reproduce (200 s: 42 with 7877 and 8004 algebras; cap-dependent). Parents are walk members, so Coxeter key == class base (42/42, trivial by BFS). alg.rels irredundant in all 42. Every parent: 2 out-arrows at v, a relation through v into each = E-103 two-out kind, not a new mechanism. New only that the walk reaches it at n=8 (none n<=7). Sample 1 has `4513=4573` beside `451=0` (monomial in disguise): presentation not minimal in the split sense. I did not extract the kernel element x.
- r023 (T5): off walks "genuine long relation => reject" holds; redundant long relations fooled hasLongSquare (362 steps). "reject => long square" false: 2-out vertex, long relation on each out-arrow (12/4463 n=6, 11/2388 n=7). Two n=6 examples have Coxeter polynomial in no LNA/dual class: unreachable by the key-preserving walk. Conjecture for theorist: fails iff x!=0 in e_aAe_v with x*beta in I for all out-arrows beta.
- r021 null: "exactly one placement in S" not a 4-property. r018: split counts 6/9/12/15/18/18 (n=12..17). r015: E-075 "20 of 25" not reproducible (11 of 25). Unit of evidence = orbit.

## Tried
- r025: skeptic_n8.py (n=8 class 0 walk, collect hits, class/shape/redundancy; 5 min each, two runs).
- r023: skeptic_offwalk.py, skeptic_reach.py (n=6 class 0, 400 s). r021 nulls; r018 rowset16; r015 rowset; r013 orbscan; r010 probe; r007 nulls.
- Bug lessons: cache by id() reuses ids; pkill/pgrep -f match own shell; python output to file is buffered (use -u); exec-ing another script clobbers my _argv (save under a new name); `| tail` hides progress.

## Not done
- Kernel element x for the 42; classes 1, 2 at n=8; n=7 class 1/2; exhaustive small-quiver enumeration; parallel arrows / 2-arrow sums off walks; whether off-class (r023 n=6) examples have any walk-reachable analogue.
- From r021: n=15..17 null; where other placements go.

## Next
1. Extract x for the 42, check whether it is a path-difference of fixed length.
2. Referee theorist proof of the x*beta in I statement.
3. Minimal-presentation test (zero relation vs commutativity) before any shape test.
## Habits
- Check vacuity by definition (class membership of walk parents is trivial); run a control; check timeouts/caps; check presentation dependence.
