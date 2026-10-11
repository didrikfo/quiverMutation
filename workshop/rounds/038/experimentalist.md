# On the n = 8 c0 walk to 8 465 expansions, J_i != 0 with d_i >= 3 occurs in 9 rows / 5 algebras, all with out-degree(i) = 3; out-degree(i) = 3 is necessary in the sample, not sufficient

author: experimentalist · round: 038 · kind: result
thread: T5 · bears on: E-131, E-133, E-134, E-128 (L1)

## Claim

Re-running E-133's walk (n = 8, class c0, key-preserving BFS, acyclic algebras, gate-admitted v; 560 s, stopped by the cap at 8 465 expansions, 8 443 distinct keys, BFS level 9 of 10 partial) and saving every (alg, v, i) with a path i -> v and (J_i != 0 or d_i >= 3) gives 378 rows of 56 261. (d, dim J) counts: (2,1) 234, (3,1) 2, (4,1) 4, (5,1) 3, (3,0) 119, (4,0) 16. So J_i != 0 occurs at d = 2 (234 rows, 201 algebras), and at d >= 3 in 9 rows (5 algebras: dim A 47, 47, 64, 64, 75), always dim J_i = 1 (L1 holds, no row with j > d-1). Among J_i != 0 rows, the only invariant in my table that separates d = 2 from d >= 3 without being a restatement of d is out-degree(i): d >= 3 rows have out-degree(i) = 3 (9 of 9), d = 2 rows have 1 (35), 2 (191) or 3 (8). Not sufficient (8 d=2 rows have out-degree 3), and the d >= 3 rows are all in the last 670 expansions (exp >= 7798), so "necessary" is a statement about a thin tail. Not claimed: any statement past the cap; that the 5 algebras are pairwise non-isomorphic (ids are canonicalKey where not None, else a hash of quiver + Cartan matrix; two of the 5 have keys, three hashes).

New relative to E-133: two (3,1) rows (exp 7798, 7810; v = 7, i = 5, dim A 47, out(i) targets {1,3,4} all distinct, no parallel arrow, out-degree(v) = 2). E-133 filtered on out-degree(v) >= 3 and so missed them; the E-133 rows (7822, 7831, 7836) reproduce with (d, out(i), out(v)) = (4,3,3), (4,3,3), (5,3,4). So "d >= 3 at J_i != 0 needs parallel arrows" is false; and d = 3 with J != 0 occurs on a walk (E-134's hand (3,1) example was excluded from walks by key).

## Evidence

Table (rows with J_i != 0, split by d; "algs" = distinct algebra ids, "Cartan row" = dim e_iAe_j over all j, i.e. row of the Cartan matrix of i):

| d_i | rows | algs | out(i) | out(v) | dim A | Cartan row sum | max entry |
|---|---|---|---|---|---|---|---|
| 2 | 234 | 201 | 1:35, 2:191, 3:8 (of the 8, 4 have 2 distinct targets, 4 have 3) | 1:23, 2:207, 3:4 | 26..39 | 5..8 | 2 |
| 3 | 2 | 2 | 3 (targets 1,3,4 distinct) | 2 | 47 | 11 | 3 |
| 4 | 4 | 2 | 3 (targets {x,6,6}) | 3 | 64 | 13 | 4 |
| 5 | 3 | 1 | 3 (targets {4,4,4}) | 4 | 75 | 13 | 5 |

Rows with J_i = 0 and d >= 3 (135 rows, 59 algebras, 66 (alg, v)): out(i) from 1 to 7 (3:11, 4:6, 5:11, 6:4, 7:1, 2:93, 1:9); dim A 31..75; Cartan row sum 9..13. So d >= 3 alone does not force out(i) = 3; the J_i != 0 case is the one that picks out 3.
By BFS level: J != 0 rows at level 6: 4, 7: 12, 8: 81, 9: 146; d >= 3 and J != 0 at level 9 only; d >= 3 and J = 0 at levels 8 (18) and 9 (117). Everything d >= 3 sits in the last two levels, so the sample for d >= 3 is the tail of a capped walk.
Not separating: dim J (always 1), dim A and the Cartan row sum (both are dominated by d and by walk depth; d = 2 rows have dim A up to 39, d >= 3 rows 47..75, but the J = 0 d >= 3 rows go down to 31, so size alone is not it), out(v) (2 for the two new rows, 3 or 4 for the E-133 rows, mostly 2 for d = 2).
The 2 new and the 7 E-133 rows: Cartan rows have 2, 3 or 4 zero entries as the d = 2 rows do (1..4). The E-133 observation e_iAe_8 = 0 at the i with J_i != 0 is about other i and is not retested here.
Full per-row table (d, dim J, out(i), out(v), targets, Cartan row, algebra id, level, expansion) is in `experimentalist_d3table_n8c0.pkl` (378 dicts); the printed tables are `experimentalist_d3tab.txt`.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/038/experimentalist_d3table.py 8 0 560 12000 workshop/rounds/038/experimentalist_d3table_n8c0.pkl   # 560 s, cap-bound
.venv/bin/python workshop/rounds/038/experimentalist_d3tab.py workshop/rounds/038/experimentalist_d3table_n8c0.pkl > workshop/rounds/038/experimentalist_d3tab.txt
```

The cap is time-based: expansion counts are load-dependent (E-133 needed 7 840 to reach its rows; I reached 8 465 and the rows at 7798/7810 are earlier than those). A rerun under more load may stop before them.

## Prior record

E-131 (d = 2 at all 285 capped rows, d >= 3 only at J = 0), E-133 (3 d >= 3 rows past the cap, filtered on out-degree(v) >= 3, tables unsaved), E-134 (hand (3,1) example, n = 6 enumeration), E-128 (L1). Grep of EXPERIMENTS.md for "out-degree(i)" / "(3, 1)" on a walk finds nothing: the two (3,1) walk rows with distinct targets and the out(i) = 3 observation are new; d = 2 at 234 rows reproduces E-131's count order (n = 8 c0 had 109 there on 5.9k expansions).

## Code changed

None (new scripts only: `experimentalist_d3table.py`, `experimentalist_d3tab.py`). No tests run.

## Next

- theorist: why out(i) = 3 (three out-arrows of i with J_i != 0 and d >= 3); the two new rows have distinct targets, so a Gamma_i / parallel explanation is not needed there. What is the kernel (skeptic_kernel.py) on 7798/7810?
- skeptic: isomorphism of the 5 d >= 3 algebras beyond ids; 7798 vs 7810 (same dim A 47; both have keys) and 7822/7831.
- toolsmith/overnight: closed or deeper n = 8 c0 walk (levels 9, 10) to see whether out(i) = 3 survives or d >= 3 with out(i) = 2 appears; c1, c2 with this script (arguments `8 1`, `8 2`), and n = 9 c0 (`9 0`), each in its own command.
