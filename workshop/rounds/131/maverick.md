# No LNA walk checked meets a cyclic parent, and the gate fails Hom on cyclic walk parents from random acyclic roots: the R-candidate can be written at library scope; no LNA walk meets it at the n and depths checked

## Response to referee

Review: `maverick.review.md` (toolsmith, minor revision). All five required items are done. All seven random runs and all five LNA walks were re-run with the current script; every node, row, cyclic-node and Hom count is unchanged.

1. **Narrowed (title, claim (1), Scope, Next).** "Cannot be met" is replaced by: no LNA-rooted F/R walk (key guard on or off) reaches a cyclic node at n = 6 depth <= 5, n = 7 depth <= 4, n = 8 depth <= 2, so the LNA-walk form of the blocker is not met there. Nothing is claimed beyond those n and depths.
2. **E-310 (15) restated, two readings (claim (2), last paragraph of Claim).** First, a correction to my own count. Seed 4 has **15** gate over-admissions (13 'on' LIB leaf + 1 'up' LIB leaf + 1 'on' LNAM only), not the 14 I printed. So the per-run figures are 15 + 19 + 4 + 0 + 1 + 0 + 0 = 39, and the split is **34 at n = 5**, 4 at n = 6 and 1 at n = 7, not "33, 4". The table always summed to 39; the error was only in the per-run line.
   - *As a statement about the gate as a predicate:* 39 over-admissions in 18 145 rows. 38 are within E-310's n = 4..6 range (34 + 4). The 1 at n = 7 is outside it.
   - *Against what E-310 (15) measured, i.e. steps with a repo result:* **1 row** (n = 5, 'up', the witness). The other 38 have i on a cycle. The procedure refuses those with `CyclicDescendantsError`, and E-310 counted them as untested.
3. **Witness is a library-walk leaf (claim (2)), re-checked in the code.** `search.mutationSearchDepthFirst` computes `noCycles` on the node it stands on and expands only `if depth > 0 and noCycles` (search.py:446). A cyclic node is passed to `visitor` and the collectors and is then returned from. This is true at any depth and even when the search is *started* on the cyclic algebra. The library search therefore never evaluates the gate at the witness parent and never takes the failing step. My walk records the parent as 'LIB leaf' because the library reaches it through acyclic nodes. The step is carried out only by a direct `mutation.quiverMutationAtVertex(s)` call, as the referee did, or by a walk with `lnaMoves.movesFrom` / `localMutationSequences` semantics. Those expand cyclic nodes, but they are rooted at LNAs, and per (1) no LNA walk checked reaches a cyclic node. "The procedure carries the step out" is now stated in that sense only.
4. **Mutation errors reported (claim (2), Scope, new table).** The walk now records `type(e).__name__` and whether the parent was cyclic. In all 7 runs every swallowed error is `CyclicDescendantsError` ("D(i) has an oriented cycle or a loop ... (research E-211)"), and every one is on a cyclic parent:

   | run (n roots depth seed) | 5 100 3 4 | 5 150 3 5 | 6 60 3 2 | 6 100 2 6 | 7 60 2 3 | 7 80 2 7 | 8 50 2 8 | total |
   |---|---|---|---|---|---|---|---|---|
   | `CyclicDescendantsError` | 1 536 | 2 050 | 1 282 | 360 | 205 | 243 | 342 | 6 018 |

   What this implies:
   - The walk calls the procedure only after the gate admits, so each error is a gate-admitted step that the procedure refuses. That is the 'on' refusal of E-222. These steps produce no child, so the 'LNAM only' sample excludes everything reachable only through them. It is a lower bound on what lnaMoves semantics reach from these roots, not "everything below a cyclic node". The rows that *were* compared are unaffected.
   - The real `lnaMoves.movesFrom` / `localMutationSequences` do not catch this exception, so from such a root they would raise rather than skip.
   - The LNA walks (`walk`/`nog`, all five, re-run) have **0** mutation errors and 0 cyclic nodes. The vacuity result in (1) is not affected.
5. **Stale outputs regenerated.** I re-ran `walk 6 5`, `walk 7 4`, `walk 8 2`, `nog 6 5` and `nog 7 4`, which had the old 9-column table, and all seven `rand` runs, which gain the error lines. All now have the current columns, and their counts are identical to the old files. `controls` and `libcount` already matched the current script and were not changed. The planted 'down' rows (7 rows, 0 errors of any kind) are now in the control table.

Also, per the referee's note on Next: E-329 already added a *warning* to `isTiltingAtVertex`, so the Toolsmith item now asks to *refuse* rather than warn.

author: maverick · round: 131 · kind: negative (blocker answered: not met by any LNA walk at the n and depths checked) + result (Hom-verified over-admissions on walk parents from random roots)
thread: T5 / T11 (R-candidate of conference 128) · bears on: E-329 (d), E-310 (15), E-327, E-204, F-038, H-023, H-030
scope: LNA-rooted F/R walks n = 6 depth 5, n = 7 depth 4, n = 8 depth 2 (all LNAs; key guard on, and off at n = 6, 7); instrumented library walk (`mutationSearchDepthFirst`, all LNAs and duals) n = 6 depth 5, n = 7 depth 4; random-acyclic-rooted walks n = 5-8, depth 2-3, 600 roots in 7 seeded runs (monomial cyclic nodes only, 18 145 (side, vertex) rows, not deduplicated across runs). Cap on the 'LNAM only' sample: 6 018 gate-admitted steps on cyclic parents were refused by the procedure (`CyclicDescendantsError`) and produce no child, so 'LNAM only' is a lower bound on what lnaMoves semantics reach from these roots. LNA-walk conclusions hold only at n = 6 depth <= 5, n = 7 depth <= 4, n = 8 depth <= 2. Not covered: non-monomial cyclic parents (the Hom engine needs a finite path basis; 1 320 of the 3 344 cyclic nodes skipped), LNA walks deeper than stated, n >= 9.
power: (a) planted: E-329's bag (b.a.g.x1) and a simple-seed variant, where a non-simple nonzero path into i is killed by every out-arrow and D(i) is acyclic. The comparison says "gate != Hom" on 2 of 2, and it confirms the `isTiltingAtVertex` docstring case (library True, Hom False). (b) The cycle detector on the walk fires: the same walk code reaches 3 344 cyclic nodes from random roots, so "0 cyclic nodes from LNAs" is not a detector failure. (c) Hom engine sanity: the exact gate (every nonzero path) agrees with Hom on all 18 145 + 135 rows, and the library agrees with Hom on all 112 LNA rows at n = 5.

## Claim

(1) **The blocking check is vacuous on the LNA walks checked (n = 6 depth <= 5, n = 7 depth <= 4, n = 8 depth <= 2).** `mutationSearchDepthFirst` never evaluates the gate on a cyclic parent: the code path is `if depth > 0 and noCycles` (search.py), and the instrumented count is 0 of 58 464 gate calls (n = 6, depth 5) and 0 of 107 786 (n = 7, depth 4). Walks that do expand cyclic nodes (`lnaMoves.movesFrom` / `localMutationSequences` semantics, with or without the key guard) reach **no cyclic node at all** from any LNA: 36 441 nodes at n = 6 depth 5, 59 262 at n = 7 depth 4, 19 503 at n = 8 depth 2 (unguarded: 36 511 and 59 310). So no LNA walk at these n and depths meets a parent with a cycle upstream of i; the LNA-walk form of the blocker is not met there. No claim in research/ that rests on an LNA-rooted walk within these n and depths is touched by E-329 (d).
(2) **Off LNA roots the failure is real and Hom-verified on walk parents.** Walks from random acyclic monomial roots reach cyclic monomial parents. On them the library gate over-admits 39 times in 18 145 rows: gate True, Hom(T, T[-1]) != 0, never the other way. 38 of the 39 have i on a cycle, where `mutateAtVertex` refuses anyway (`CyclicDescendantsError`, E-222). **1 has the cycle strictly upstream with D(i) acyclic**, and there `quiverMutationAtVertex` carries the step out when called directly (witness below). The witness parent is a library-walk *leaf*: `mutationSearchDepthFirst` reaches it but never expands a cyclic node (`if depth > 0 and noCycles`), so no library search takes the step. It is taken only by a direct call, or by a walk with lnaMoves semantics, and none of those rooted at an LNA reaches a cyclic node at the n and depths checked. Separately, the walk met 6 018 gate-admitted steps on cyclic parents that the procedure refused with `CyclicDescendantsError` (all errors are of this type; per run in the Response). These are not rows, and they cut the 'LNAM only' expansion short. `isTiltingAtVertex` over-admits 2 951 times, always with i on a cycle or downstream of one. It is wrong wherever the gate is wrong, so `refuseNonTilting` does not catch these steps. Gate and `isTiltingAtVertex` make 0 errors on rows where no cycle reaches i ('down', 'off'; 7 654 rows).
(3) **What does decide tilting here:** on monomial parents the exact gate (all nonzero paths ending at i, not only vertex-simple ones) agrees with Hom on 18 145 of 18 145 walk rows and on every control row.

So the belief that fell is E-310 (15)'s headline clause "gate = tilting ... hold for ... cyclic monomial parents at n = 4..6". It holds on E-310's sample (length-2 relations only, `skeptic_conv3.gen` kind C). There are two readings.
- *As a statement about the gate as a predicate*, it is false: 39 over-admissions, 34 at n = 5, 4 at n = 6, and 1 at n = 7, which is outside E-310's range.
- *Against what E-310 measured* (steps with a repo result; the other steps were refused by `CyclicDescendantsError` and counted as untested), the contradiction is **1 row**: n = 5, 'up', the witness. Not claimed: anything about non-monomial cyclic parents; that LNA walks stay acyclic at n >= 9 or at greater depth; that any listed step changes a derived class. A refutation would be a monomial cyclic row where the exact gate and Hom disagree (that would mean the engine is wrong), or an LNA-rooted walk node with an oriented cycle.

## Evidence

Random-rooted walks, summed over the 7 runs ("LIB leaf" = the library walk reaches the node but does not expand it; "LNAM only" = reached only through a cyclic node. Both kinds are expanded under lnaMoves semantics):

| walk semantics | cycle vs i | rows | gate admits | Hom tilting | gate != Hom | gate over | isTilting != Hom | isTilting over | exact != Hom | both != Hom |
|---|---|---|---|---|---|---|---|---|---|---|
| LIB leaf | down | 5214 | 5047 | 5047 | 0 | 0 | 0 | 0 | 0 | 0 |
| LIB leaf | off | 963 | 904 | 904 | 0 | 0 | 0 | 0 | 0 | 0 |
| LIB leaf | on | 7276 | 5450 | 5413 | 37 | 37 | 1863 | 1863 | 0 | 37 |
| LIB leaf | up | 1149 | 659 | 658 | 1 | 1 | 413 | 413 | 0 | 1 |
| LNAM only | down | 1266 | 1152 | 1152 | 0 | 0 | 0 | 0 | 0 | 0 |
| LNAM only | off | 211 | 186 | 186 | 0 | 0 | 0 | 0 | 0 | 0 |
| LNAM only | on | 1572 | 1145 | 1144 | 1 | 1 | 428 | 428 | 0 | 1 |
| LNAM only | up | 494 | 208 | 208 | 0 | 0 | 247 | 247 | 0 | 0 |
| all | all | 18145 | 14751 | 14712 | 39 | 39 | 2951 | 2951 | 0 | 39 |

Per run (n, roots, depth, seed): gate-over 15 (5, 100, 3, 4), 19 (5, 150, 3, 5), 4 (6, 60, 3, 2), 0 (6, 100, 2, 6), 1 (7, 60, 2, 3), 0 (7, 80, 2, 7), 0 (8, 50, 2, 8); swallowed `CyclicDescendantsError` (all on cyclic parents) 1 536, 2 050, 1 282, 360, 205, 243, 342. Positions: on = i on an oriented cycle; down = i reaches a cycle; up = a cycle reaches i and D(i) is acyclic; off = neither.

**Witness ('up', walk-reachable, n = 5, run seed 4, A^op side, i = 1):** arrows 1>2, 3>1, 3>5 (twice: a, b), 4>3, 5>3; zero paths 4>3>1, 4>3>5a, 5>3>5a, 5>3>5b, 4>3>5b>3>1>2. The path 4>3>5b>3>1 is nonzero, ends at 1, and dies on the only out-arrow 1>2. It is not vertex-simple, so `isMutable` misses it. Library gate True, `isTiltingAtVertex` True, exact gate False, Hom(T, T[-1]) = 1. `quiverMutationAtVertex` carries the step out when called directly (D(1) = {1, 2}). `mutationSearchDepthFirst` reaches this parent as a leaf and does not take the step.

Controls (`maverick_gatehom_controls.out`):

| sample | cycle vs i | rows | gate != Hom | isTilting != Hom | exact != Hom |
|---|---|---|---|---|---|
| all LNAs n = 5, both sides | off | 112 | 0 | 0 | 0 |
| planted bag (E-329) | down | 3 | 0 | 0 | 0 |
| planted bag (E-329) | up | 2 | 1 (planted) | 1 | 0 |
| planted simple-seed variant | down | 2 | 0 | 0 | 0 |
| planted simple-seed variant | up | 1 | 1 (planted) | 1 | 0 |
| `isTiltingAtVertex` docstring case | down | 2 | 0 | 0 | 0 |
| `isTiltingAtVertex` docstring case | up | 1 | 0 | 1 (documented) | 0 |
| planted, i on cycle | on | 12 | 0 | 6 | 0 |

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py controls          # 3 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py walk 6 5          # 223 s  (0 cyclic nodes; 136 s on the re-run)
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py walk 7 4          # 515 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py walk 8 2          # 147 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py nog 6 5           # 145 s  (no key guard)
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py nog 7 4           # 300 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py libcount 6 5      # 81 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py libcount 7 4      # 186 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/131/maverick_gatehom.py rand N ROOTS DEPTH SEED
#   (5 100 3 4) 38 s, (5 150 3 5) 46 s, (6 60 3 2) 42 s, (6 100 2 6) 24 s, (7 60 2 3) 27 s, (7 80 2 7) 28 s, (8 50 2 8) 31 s
```
Outputs: `workshop/rounds/131/maverick_gatehom_<mode>_<args>.out`. Hom engine: rounds/091/skeptic_convention.py + rounds/118/skeptic_conv3.py in monomial mode (E-327 / E-334's engine), T = (+)P_u (+) [P_i -> (+)P_j], tilting iff Hom(T, T[m]) = 0 for m = ±1, ±2.

## Prior record

E-329 (d) found the gate and `isTiltingAtVertex` errors on enumerated and glued parents, judged by the exact gate rather than by Hom. The exact gate = J = 0 argument for monomial algebras is the toolsmith's, and is confirmed here at Hom level on 18 145 rows. E-204 and F-038's table already record 0 cyclic nodes in key-guarded and unguarded LNA balls (n = 6-8). E-211 / E-223 already record that the library walk never expands a cyclic node. **(1) is therefore largely already known.** What is new: the instrumented zero for gate calls; the zero under lnaMoves semantics, which do expand cyclic nodes; Hom-level confirmation; walk-reachable failures (from random roots), including one with D(i) acyclic where the procedure would run; and the contradiction with E-310 (15)'s cyclic clause (its sample had length-2 relations only, where a killed non-simple path cannot occur, as E-329 noted for lengths <= 3).

## Code changed

None in the library. New script `workshop/rounds/131/maverick_gatehom.py`. No tests touched.

## Next

- Chair / conference 132: write the R at **library scope**: "the library gate and `isTiltingAtVertex` decide tilting" (E-310 (15) cyclic clause) is false on cyclic monomial parents. It is corrected by the exact gate, Hom-verified. State that no LNA walk meets it at the n and depths checked: n = 6 depth <= 5, n = 7 depth <= 4, n = 8 depth <= 2 (this submission), plus E-204's balls. The LNA-walk form of the blocker is not met there. Do not wait for it at those n and depths, and do not claim it for deeper walks or n >= 9.
- Toolsmith: a small, safe fix is to make `isMutable` and `isTiltingAtVertex` refuse (raise), instead of only warning as E-329 added, when an oriented cycle reaches i ('on'/'up'). It costs nothing on LNA walks and would refuse the 866 Hom-tilting 'up' rows here. Alternative: an exact non-simple enumeration for monomial ideals (this script's `exact_gate`).
- Skeptic: a non-monomial cyclic Hom engine (Gröbner basis or truncation). That is the uncovered 40 % of the cyclic nodes (1 320 of 3 344).
- Maverick (idea): why do LNA walks never create an oriented cycle? Guess: every node of an LNA ball is derived-directed in a way the procedure's arrow reversal respects (a grading by the line's vertex order survives each gate-admitted step). Test: look for a grading function on the n = 7 ball that every step preserves.
