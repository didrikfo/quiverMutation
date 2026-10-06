# Skeptic's notebook (after round 037)

## Believe now
- r037 (T5): the 5 E-127 rows (n=8 c0, expansions 6820/7424/7822/7831/7836) are 3 algebras (dim A 39, 64, 75): two mirror pairs under the vertex swap 1<->3 and a singleton. Pairing of 7822/7831 is by invariant (product-rank table T), canonicalKey differs (non-minimal presentations). Reproduced E-127 (rerun 514 s, rows pickled in rounds/037/skeptic_rows.pkl).
- r037: J_i != 0 in all 8 pairs is the out-degree-1 relation p1 b = p2 b on the single NON-parallel out-arrow; e_iAe_8 = 0 exactly at those i (parallel arrows ground). #J_i != 0 = parallel multiplicity (2,2,3) is unexplained, may be chance. 3 of 5 rows have d_i = 4,4,5 with J_i != 0: E-129's "d_i = 2 at J_i != 0" is a cap artefact (E-129 stopped at 5958 expansions; these rows are at 7822+). L1 dimJ <= d-1 still holds.
- r034: all 5 fail checkCartan; 199 J=0 out>=3 control rows congruent; Cartan discrepancy support == J support. Gate blind at out-degree 3,4 as at 2.
- r031: 25 loose-shape W-false rows n=6,7 are "half-W". r029: 61 n=8 c0 D' rejects: gate tests single paths, J tests combinations. r026/r025: shape necessary not sufficient; irredundant != minimal; cap-dependent counts. r023: redundant long relations fooled hasLongSquare. r015: unit of evidence = orbit.

## Tried
- r037: skeptic_dump.py, skeptic_orbits.py, skeptic_mult.py, skeptic_kernel.py. r034: skeptic_outdeg3.py. r031 skeptic_loose.py. r029 skeptic_dprime.py. r026 skeptic_x.py. r025 skeptic_n8.py. r023 skeptic_offwalk.py, skeptic_reach.py. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None gives spurious equality (r037); pkill -f kills own shell, use pgrep for pid or timeout; id() caches; python -u for files; `| tail` hides progress; budget under 10 min incl. load; background job + `until grep -q done` loop.

## Not done
- Hand-built m=2 algebra with three i having e_iAe_8=0 (tests the m coincidence); mirror of 7822/7831 proved (explicit algebra iso, not invariants); other classes (c1, c2) and n=9 past the cap; the 22 rows; J-independent socle test; both-die n=6 reach; n=9 c0; E-129 histogram without cap.

## Next
1. Hand-build the m-coincidence test (cheap). 2. Re-run E-129's histogram on c0 to 8500+ expansions and see d >= 3 J != 0 counts per distinct algebra. 3. Referee any claim "d_i <= 2" or "no nn" for scope: these hold only below the cap / at out-degree 2.

## Habits
- Check vacuity by definition; run a control; check caps (E-129 cap vs row depth); presentation dependence; count pairs vs rows; check None keys before comparing keys; confirm a derived number (d_i) by an independent route (Cartan).
