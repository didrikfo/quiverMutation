# Skeptic's notebook (after round 031)

## Believe now
- r031 (T5): the 25 loose-shape W-false rows at n=6,7 c0 (5+20) are all "half-W": the monomial kills one term of R=p1 b1=p2 b1, the other term survives b2, so x*b2 != 0 and J=0. Same pattern as the 23 L-true accepts at n=8 c0. Gate and W identical across n; K == Both-die in all 109 pairs. n=6,7 rows are all arrow-vs-path squares (lengths (2,1),(1,2)); n=8 half rows also (2,2),(2,3). No both-die row at n=6,7 in capped walks: walk reach or algebra, undecided.
- r029 (T5): 61 n=8 c0 D' rejects admitted because the gate tests single paths and J != 0 tests combinations; near tautology. Out-1 analogues at n=6,7,8 c0 (15,26,4). Out-2 rejects only n=8 c0 in capped sweeps. W's both-terms-die clause is what separates L.
- r026: 42 rejects x=p+q lengths (2,2) 37,(3,3) 5,(2,3) 4; shape necessary, not sufficient (767/42). 22 rows (both J_beta nonzero, J=0) unexplained (E-113: half-W 20, loose pendants 6).
- r025: 42 reproduce; key==class trivially; alg.rels irredundant not minimal; cap-dependent counts.
- r023: redundant long relations fooled hasLongSquare. r021: "exactly one placement in S" not a 4-property. r015: E-075 "20 of 25" is 11 of 25. Unit of evidence = orbit.

## Tried
- r031: skeptic_loose.py (per R,M pair pattern table; n=6,7,8 c0). r029: skeptic_dprime.py. r026: skeptic_x.py. r025: skeptic_n8.py. r023: skeptic_offwalk.py, skeptic_reach.py. Earlier nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: id() caches; pkill/pgrep -f match own shell; python -u for files; `| tail` hides progress; parallel jobs slow c0; budget under 10 min incl. load; foreground tool call is cut at 120 s unless timeout is passed.

## Not done
- Second-kill mechanism in Both rows (my one-M flag misses it); hand-built both-die n=6 control (walk reach test); n=9 c0 (checkpoint); the 22 rows; exhaustive small-quiver enumeration; classes 4..10 beyond prefixes.

## Next
1. Hand-build an n<=7 both-die square, check gate and reachability. 2. Classify Both rows' second kill. 3. Referee any proof of J != 0 <=> circuit in Gamma_i.

## Habits
- Check vacuity by definition; run a control; check caps; presentation dependence; count pairs vs rows before quoting equal numbers.
