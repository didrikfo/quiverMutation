# Skeptic's notebook (after round 023)

## Believe now
- r023 (T5): off the guarded walks "long square => reject" holds if the square is genuine (truncation x = sum c p[:-1] nonzero mod I; it is the kernel element of the one map). A long relation in `alg.rels` can be redundant (killed by e.g. abd=acd): then hasLongSquare is True at a tilting step (362 of ~65k gate-admitted tilting steps in hand-built parents). So E-100's test needs a minimal presentation off the walks.
- r023: "reject => hasLongSquare" is false literally: 2-out vertex with a long relation on each out-arrow (kernel killed by both), shared longer suffix abdef=acdef at v=e, sum relation reduced by a monomial. 12/4463 n=6, 11/2388 n=7 gate-admitted rejections. Two examples have Coxeter polynomials in no LNA/dual class at n=6 => unreachable by the key-preserving walk. Control: n=6 padded E-078 IS reached (class 0).
- Right statement (conjecture for theorist): step-7 map fails iff there is x != 0 in e_aAe_v with x*beta in I for all out-arrows beta.
- r021 null: "exactly one placement in S" for 444-words is no 4-property (also words without 4); gap pattern is a property of S. r018: split counts 6/9/12/15/18/18 (n=12..17), 0 partial vacuous. r015: E-075 "20 of 25" not reproducible (11 of 25). r013, r010, r007 as before. Unit of evidence = orbit.

## Tried
- r023: skeptic_offwalk.py (random acyclic quivers n=5..7, 0-3 relations, modes A seeded square / C none; 5 runs ~23k parents), skeptic_reach.py (guarded BFS membership by canonicalKey, n=6 class 0, 400 s).
- r021 skeptic_null*.py; r018 rowset16, partial_where; r015 rowset; r013 orbscan; r010 probe; r007 nulls.
- Bug lessons: cache by id() reuses ids; `pkill -f`/`pgrep -f` match own shell (killed it; use pid + kill -0); python output to file is buffered (use -u).

## Not done
- Exhaustive enumeration (not random) of small quivers; n=8; length-2 sum relations and parallel arrows not generated; reachability of kind c and of n=7 examples; whether kinds a-c have non-tilting child confirmed independent of tiltingPlus (End(T) not computed).
- From r021: n=15..17 null; where the other placements go.

## Next
1. Exhaustive n=6 enumeration of kinds a-c; add parallel arrows / 2-arrow sums.
2. Ask experimentalist to rerun E-100 with genuine-truncation test and count 2-out vertices on walks.
3. Referee any theorist proof of the x*beta in I statement at the multi-out case.
## Habits
- Check whether a pass is vacuous by definition; run a control stratum; check timeouts; check a test's presentation dependence (alg.rels vs relationsFrom, redundant relations).
