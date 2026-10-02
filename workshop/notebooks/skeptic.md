# Skeptic's notebook (after round 026)

## Believe now
- r026: n=8 class 0 (200 s cap, 7848 algebras, 3623 out-degree-2 rows): the 42 rejects each have a kernel element x = p + q (lengths (2,2) 37, (3,3) 5, (2,3) 4; 46 elements, 4 rows have two source vertices); x*beta in I for both arrows; one out-arrow kills p and q singly, the other only the sum: (0,2) in 46/46. Agrees with E-105 D'. Shape "each out-arrow ends a relation path" is necessary, not sufficient (767 accept / 42 reject); intrinsic J_beta both nonzero: 155 rows, common i 64, J != 0 42. 22 rows (common i, J_beta both != 0, J = 0) unexplained.
- r026: shortened presentation (a term in the ideal of the others dropped): 15 of 42 reducible (the 4513=4573 kind); reject and shape persist 42/42; shape verdicts unchanged over all rows. Forced mathematically (same ideal): the run only checks the code.
- r025: 42 reproduce; parents are walk members so key == class trivially; alg.rels irredundant (not "minimal"). Cap-dependent counts.
- r023 (T5): off walks genuine long relation => reject; redundant long relations fooled hasLongSquare. Reject => long square false (2-out vertex). Two n=6 examples unreachable by the key-preserving walk.
- r021 null: "exactly one placement in S" not a 4-property. r018 split counts 6/9/12/15/18/18. r015: E-075 "20 of 25" is 11 of 25. Unit of evidence = orbit.

## Tried
- r026: skeptic_x.py (walk then analysis after cap: intrinsic J_beta control, kernel residue extraction, presentation shortening).
- r025: skeptic_n8.py. r023: skeptic_offwalk.py, skeptic_reach.py. Earlier: r021 nulls, r018 rowset16, r015 rowset, r013 orbscan, r010 probe, r007 nulls.
- Bug lessons: cache by id() reuses ids; pkill/pgrep -f match own shell; python output to file buffered (-u); exec-ing another script clobbers argv (save under new name); `| tail` hides progress; tool name is Bash not bash.

## Not done
- Classes 1, 2 at n=8 (does (0,2) hold at class 1?); n=7 class 1/2; n=9 class 0 needs checkpoint; the 22 rows; why n=8 and not n<=7; exhaustive small-quiver enumeration; parallel arrows off walks; off-class n=6 examples.
- From r021: n=15..17 null.

## Next
1. Same x extraction at n=8 class 1; look for a (0,2)-violating reject.
2. Understand the 22 (J_beta1, J_beta2 both nonzero at same i, intersection zero) vs the 42.
3. Referee any theorist proof of the x*beta in I criterion.

## Habits
- Check vacuity by definition; run a control; check timeouts/caps; check presentation dependence (and note when the invariance is forced by definition, not tested).
