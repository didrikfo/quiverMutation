# Skeptic's notebook (after round 029)

## Believe now
- r029 (T5): the 61 n=8 c0 D' rejects are admitted because the gate tests single paths and J != 0 is the same test on combinations; near tautology. Out-1 analogues (long square) at n=6,7,8 c0 (15, 26, 4). Out-2 rejects only in n=8 c0 in my capped sweep (c1..c10 n=8, all n=6,7 classes: 0). Loose shape L (commute into b1 + monomial into b2): n=8 c0 84 rows, 61 reject, 23 accept; n=7 c0 20 and n=6 c0 5, all accepted. W clause separates; L alone is not sufficient (the null). Total 53 162 algebras.
- r026: n=8 c0 42 rejects: x = p + q, lengths (2,2) 37, (3,3) 5, (2,3) 4; one out-arrow kills p, q singly, other only the sum. Shape "each out-arrow ends a relation path" necessary, not sufficient (767/42). 22 rows (both J_beta nonzero, J = 0) unexplained. 15 of 42 reducible; verdict unchanged (forced by definition).
- r025: 42 reproduce; key == class trivially; alg.rels irredundant not minimal; cap-dependent counts.
- r023: off walks genuine long relation => reject; redundant long relations fooled hasLongSquare. r021: "exactly one placement in S" not a 4-property. r015: E-075 "20 of 25" is 11 of 25. Unit of evidence = orbit.

## Tried
- r029: skeptic_dprime.py (L/W/K tally, witness patterns, class sweeps n=6,7,8). r026: skeptic_x.py. r025: skeptic_n8.py. r023: skeptic_offwalk.py, skeptic_reach.py. Earlier: nulls r021, r018 rowset16, r015, r013, r010, r007.
- Bug lessons: id() caches; pkill/pgrep -f match own shell; python -u for files; `| tail` hides progress; tool is Bash; six parallel jobs slow c0 to 462 s: run c0 first/solo; always budget under 10 min incl. load.

## Not done
- Why L-true rows at n=6,7 c0 are accepted (W false): inspect; the 22 rows; n=9 c0 (needs checkpoint); parallel-arrow positive control; exhaustive small-quiver enumeration; classes 4..10 beyond prefixes.

## Next
1. Look at the n=6/7 c0 L-true accepts (which W clause fails).
2. The 22 rows. 3. Referee any theorist proof of J != 0 <=> circuit in Gamma_i.

## Habits
- Check vacuity by definition (the gate/J relation is definitional); run a control; check caps; check presentation dependence.
