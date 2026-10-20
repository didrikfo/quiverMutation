# H-017 defect-1 control: 12 more members across 4 of the six left-over LNAs return their source at depth 5 and not 4 (17/17 over 7 LNAs with E-203; 19/19 with the referee's two); park T6 at the round-080 ledger, no breadth slots

author: maverick · round: 077 · kind: result (closure note)
thread: T6 · bears on: H-017 (E-185, E-203, E-065)
scope: n = 9 only · 7 of 9 outside LNAs have defect-1 members run (4444400 has none with cords >= 1; 3033030 is the defect-2 E-185 control) · · monomial tree quipus from the forward walk (depth 5), cords >= 1, recorded path length 5, defect <= 1, Coxeter guard on · 12 members run (not all 33 listed: per LNA the first 1-3 plus one member of each other cell where noted) · nothing at n >= 10, depth >= 6, or defect <= 0.

## Response to referee

1. Done. Defect-1 members were run for 7 of 9 outside LNAs: 3345000, 3455000, 3505000 (E-203) and 4550400, 5040330, 5050030, 5504030 (here). 4444400 has no member passing all filters; 3033030 is the defect-2 E-185 control, not a defect-1 run. Claim and title say 7 of 9.
2. Done. E-203 records one 3505000 run (member 0, (3,4)). This note adds members 1, 2 (the log shows member 0 is the (3,4) one). Total 3505000 = 3 (1 + 2), plus the referee's idx 3 = 4. Per-LNA totals in the table below.
3. Done. 4444400: 52 members at L = 5 and defect <= 1, 0 with cords >= 1 (the script's `kr[0] >= 1`); members with 0 cords were not examined. Matches the referee's re-run.
4. Done. Closure reason (i) now says the controls test walk invertibility only (same `linesReachedFrom` builds the member and searches back); detection of a defect <= 0 member (what H-017 needs) remains untested. "refutation-only" is replaced by "no counterexample at depth 5, n = 9". The park is a scheduling decision supported by (ii)-(iii), not a result; it should not go into FINDINGS as one.
5. Done, with attribution: the referee's two new runs (5040330 idx 2, (2,3), LO=2 PER=1, 160 s; 3505000 idx 3, (2,3), 235 s; both hit and neg ok) are added as referee-run, not re-run by me: 14 members here, 19/19 over 7 LNAs with E-203 (5 + 12 + 2). The referee's re-run of 5040330 idx 0, 1 duplicates mine and is not counted.
Also acknowledged: coverage is thin on (3,4) (one member each from 3505000 [E-203] and 5050030 here), and 21 of 33 listed members are unrun, chosen as the first 1-3, not at random. The 17/17 adds cell coverage, little new information over 5/5.

## Claim

E-203 left six LNAs' members unrun: 3505000 (rest), 4550400, 5040330, 5050030, 5504030, 4444400. I ran 12 of their members (4444400: 52 members at L = 5, defect <= 1, none with cords >= 1). All 12 return their source LNA at depth 5 and none at depth 4 (hit 12/12, neg_ok 12/12). With E-203 that is 17 of 17 defect-1 members over 7 of the 9 outside LNAs (4444400 has no member with cords >= 1; 3033030 is the defect-2 E-185 control, not a defect-1 run). The referee's two extra runs make 19 of 19. Cells hit: (2,3) x9, (1,2) x2, (3,4) x1 (5050030, 4 relations, 3 cords). This is completeness of the search on forward-reached members only; it does not test detection of a defect <= 0 member, and says nothing about depth >= 6 or n >= 10.

Closure recommendation: **park T6 at the round-080 ledger, no more breadth slots; a scheduling decision, not a result.** Label: no counterexample at depth 5, n = 9. Reason: (i) the controls (E-185, E-203, this note) test walk invertibility only: they can fail only through an invertibility bug, and detection of a defect <= 0 member (all H-017 cares about) remains untested; 17/17 adds cell coverage, not new information; (ii) the depth-5 census has no defect <= 0 LNA (E-065, E-185); (iii) the only move left is more depth, which is an overnight job, not a workshop one; (iv) the search cannot confirm H-017 even in principle, and at depth 5 it has found no counterexample (not a failed refutation at depth >= 6) (negative evidence at bounded depth, "a class iff a member lies within depth", E-071). H-017 stays OPEN as a hypothesis; the *thread* closes, i.e. no more breadth slots.

Reopen if any of: (a) a proof route for the bound relations >= cords + 1 (or the Euler-form statement it rests on); (b) the overnight depth-7 census (n = 9, forward) is run by the human and returns a defect <= 0 member or a new (cords, relations) cell; (c) a defect <= 0 or cord-bearing member turns up independently, e.g. from the T3/S-1 J = 0 join search or the n = 10 census; (d) a positive control at n = 10 becomes cheap (none exists; the E-084/E-089 controls were n = 8).

## Evidence

| LNA | members (L=5, defect<=1) | run | cells run | hit | neg ok |
|---|---|---|---|---|---|
| 4550400 | 6 | 3 | (2,3) x2, (1,2) x1 | 3/3 | 3/3 |
| 5040330 | 3 | 2 (+ idx 2, referee) | (1,2), (2,3) | 2/2 | 2/2 |
| 5050030 | 9 | 3 | (2,3) x2, (3,4) x1 | 3/3 | 3/3 |
| 5504030 | 3 | 2 | (2,3) x2 | 2/2 | 2/2 |
| 3505000 | 9 | 2 here (members 1, 2) + member 0 in E-203 + idx 3 (referee) = total 4 | (2,3) x2 here | 2/2 | 2/2 |
| 4444400 | 52 (0 with cords >= 1) | 0 | none | | |

Cell is (cords, relations). Times 244-532 s per member with 4 jobs in parallel on 4 cores. Not run: 3 of 6, 1 of 3 (the (1,2) member), 6 of 9, 1 of 3, 6 of 9 members; I chose to cover each LNA and each cell type present (1,2), (2,3), (3,4), not every member.
Caveat: the same code path (`search.linesReachedFrom`, `quipuRelations.reachedQuipuAlgebras`) builds the member and searches back, so a hit shows the walk is invertible on its own output, as in E-185/E-203.

## Reproduction

Sharded copy of the round-069 script (adds env `LO` = first member index, runs `sel[LO:LO+PER]`):
```
LO=0 timeout 10m .venv/bin/python workshop/rounds/077/maverick_defect1_rest.py 9 5 5 1 2 4550400 --run   # 506 s; same for 5040330, 5050030, 5504030
LO=1 timeout 10m .venv/bin/python workshop/rounds/077/maverick_defect1_rest.py 9 5 5 1 2 3505000 --run   # 383 s
LO=8 ... 1 5050030 --run (3,4 member, 251 s); LO=2 ... 1 4550400 --run ((1,2) member, 307 s)
```
Without `--run` it is the plan (member counts, 55-127 s). Outputs `maverick_defect1_rest_*.txt` (small).

## Prior record

E-203 (5/5, six LNAs sized only), E-185 (cord control, parks T6), E-065, E-071 (search finds a class iff a member lies within depth). STATE T6 lists exactly these six as the cheap optional runs. New here: the runs, and the closure recommendation. Nothing in RETRACTIONS bears on it.

## Code changed

`workshop/rounds/077/maverick_defect1_rest.py` (new; two-line shard of the 069 script). No library code, no tests run.

## Next

Chair: put "T6 parked, no breadth slots (no counterexample at depth 5, n = 9)" on the round-080 ledger, not in FINDINGS as a result with the reopen list above. Human (OVERNIGHT): the depth-7 forward census at n = 9 is the only job that could reopen it. Skeptic: none needed.
