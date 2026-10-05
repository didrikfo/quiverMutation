# Skeptic's notebook (after round 034)

## Believe now
- r034 (T5): on the n=8 c0 walk (8456 expansions) all 5 out-degree >=3 rows with J != 0 (4 at out-3, 1 at out-4, all parallel arrows into 8) FAIL checkCartan; 199 J=0 out>=3 rows congruent. Cartan discrepancy support == {(v,i): J_i != 0} in all 5: socle reading (E-122) not refuted, but only 5 rows, 1-2 orbits possibly (mirror pairs unchecked). So rejects exist at out-degree 3,4; gate blind there as at out-2. E-111's "0 rejects" is a prefix statement. Relations on these rows are messy (repeated paths, 3 terms): E-110 circuit lemma not applicable.
- r031 (T5): 25 loose-shape W-false rows at n=6,7 c0 are "half-W" (one term killed, J=0); same as 23 n=8 accepts; K == Both-die in 109 pairs. n=6,7 capped walks have no both-die row: reach or algebra, undecided.
- r029: 61 n=8 c0 D' rejects admitted because gate tests single paths, J != 0 tests combinations. r026: 42 rejects shape necessary not sufficient. r025: 42 reproduce; alg.rels irredundant not minimal; cap-dependent counts.
- r023: redundant long relations fooled hasLongSquare. r021: "exactly one placement in S" not a 4-property. r015: E-075 "20 of 25" is 11 of 25. Unit of evidence = orbit.

## Tried
- r034: skeptic_outdeg3.py. r031: skeptic_loose.py. r029: skeptic_dprime.py. r026: skeptic_x.py. r025: skeptic_n8.py. r023: skeptic_offwalk.py, skeptic_reach.py. Earlier nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: id() caches; pkill/pgrep -f match own shell; python -u for files; `| tail` hides progress; parallel jobs slow c0; budget under 10 min incl. load; background job + sleep with timeout param, foreground cut at 120 s otherwise.

## Not done
- Distinct orbits among the 5 rows (canonicalKey of mirror/dual); W/Gamma_i classification of them; is #J_i != 0 == parallel multiplicity (n=5 only); a J-independent socle test (coker, Hom(N,N[-1])); both-die n=6 hand control reach; n=9 c0; the 22 rows; exhaustive small enumeration.

## Next
1. Hand-build an out-degree 3 parallel row with 2 parallel arrows and test #J_i vs multiplicity (cheap, would kill or support the count). 2. Orbit-dedupe the 5. 3. Referee any proof of J != 0 <=> circuit in Gamma_i; check it covers 3-term relations and parallel arrows.

## Habits
- Check vacuity by definition (J vs socle is tautological; use Cartan as the independent side); run a control; check caps; presentation dependence; count pairs vs rows before quoting equal numbers.
