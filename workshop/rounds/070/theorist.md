# Left mutation of an LNA row obeys the mirror of the E-187 rule and Catalan factorisation C_{v-2} C_{n-v}, checked directly on the repo's left mutation for n = 4..10 (all rows, all v; n = 11 by referee); the argument is a reduction to E-187 via `_opposite`, not an independent proof

## Response to referee

1. Truncated sentences deleted; the proof paragraph is rewritten with the coordinate map (v' = n+1-v, (a,b) -> (n+1-b, n+1-a)) stated once.
2. Title, kind and Scope now say "reduction to E-187 via `_opposite`", relative to the repo's reading; n = 11 added to the table as the referee's result (7m22s, 16796 LNA, 75582 not, no violations), marked as not rerun by me.
3. Cited F-025 / E-026: the mirror is the relation dual, and the direct run is the real check.
4. Output saved: `workshop/rounds/070/theorist_leftcheck.txt`, rerun n = 4..10 (2 min), violations: none; counts identical to the table.
5. Stated in Scope: "not LNA" is checked only by a necessary condition, not classified.
Also: per-vertex lists for all n = 4..10 are in the saved file.


author: theorist · round: 070 · kind: result (reduction to E-187 + machine check)
thread: T7 · bears on: E-187, E-182, E-153, H-010, F-020
scope: left mutation (`mutation.quiverMutationAtVertices(alg, [-v])`) of every admissible LNA row, every v = 1..n, n = 4..10 (24 310 gated events at n = 10). Direct check n = 4..10 (mine) and n = 11 (referee, as reported). The argument is a reduction to E-187, itself relative to the repo's reading of steps 1-7, not the paper. Not covered: two bystanders, k >= 4, non-LNA parents.

## Claim

Write the relation starting at vertex j as the interval (j, j + r_j). Left mutation at v:
(G) the gate (right-mutation gate on the dual algebra) passes iff v >= 2 and no relation STARTS at v-1;
(L) given the gate, the result is an LNA iff v = 2 or (v-2, v) is a relation (r_{v-2} = 2); then the relation (v-2, v) is replaced by (v-1, v+1) (r_{v-2} := 0, r_{v-1} := 2; for v = n only the removal, for v = 2 only the addition);
(C) the number of rows with an LNA result at v is C_{v-2} C_{n-v}, total C_{n-1}.
So left mutation slides a length-2 relation one place RIGHT, the mirror of E-187's left slide. Not claimed: anything about the paper's wording, two-vertex windows or non-LNA results beyond a necessary condition (arrow count n-1 or n and a degree-2 vertex): "not LNA" is NOT classified, only consistent with that condition. Refuted by a row and v violating (G), (L), image or (C).

## Evidence

Reduction to E-187. In the repo, left mutation IS the opposite construction: `mutation.leftQuiverMutationAtVertex` calls `procedure.mutateLeftAtVertex` = `_opposite(mutateAtVertex(*_opposite(quiver, relations), vertex))`; the gate is `mutationIsPossibleAtVertex` on `dualPathAlgebra` (as `lnaMoves.movesFrom` does). Coordinate map (stated once): reversing arrows sends vertex v to position v' = n+1-v and a relation interval (a, b) to (n+1-b, n+1-a); r^op is the reversed row. Left at v on r = right at v' on r^op, and E-187 applies verbatim to r^op. Translating its statements: gate fails iff a relation ends at v'+1, i.e. one starts at v-1 in r; v' <= n-1 iff v >= 2; LNA iff v' = n-1 (v = 2) or r^op has a length-2 relation starting at v', i.e. (v-2, v) is a relation of r; image: E-187 removes (v', v'+2) and adds (v'-1, v'+1), which mirrors to removing (v-2, v) (needs v >= 3) and adding (v-1, v+1) (needs v <= n-1). (C): E-187's bijection counts right-mutation rows at v' as C_{v'-1}... mirrored, the rows with r_{v-2} = 2 (or v = 2) and no relation starting at v-1 split into an independent LNA on 1..v-2 (C_{v-2} choices) and one on v..n (C_{n-v}), giving C_{v-2} C_{n-v}; summing over v = 2..n gives C_{n-1}. Only non-formal point: `_opposite` keeps arrow names and coefficients (its docstring), so E-187's hand argument is untouched.
Weakest step: the identification of the mirror in coordinates (position n+1-v); this is exactly what the direct run tests, since it never uses the dual rule. Why the direct run is the real check: the mirror here is the relation dual (reverse arrows AND swap right/left relations), the point of the F-025 / E-026 retraction ("the mirror was the wrong mirror"); a coordinate mirror alone would repeat that error, so the claim rests on the run against the repo's left mutation, with the reduction only explaining it.

Direct run (`theorist_leftcheck.py 4 10`; full output saved in `workshop/rounds/070/theorist_leftcheck.txt`): gate on the dual algebra equals (G) in every case; LNA iff (L); image equals the predicted row; non-LNA results have n-1 or n arrows and a vertex of degree 2 in or out; per-vertex counts equal C_{v-2} C_{n-v} (n = 8: 132 42 28 25 28 42 132 reversed into v = 2..8, i.e. the right-mutation lists of E-187 read backwards). Totals identical to E-187's table:

| n | gated | LNA | not LNA | violations |
|---|---|---|---|---|
| 4 | 10 | 5 | 5 | 0 |
| 5 | 35 | 14 | 21 | 0 |
| 6 | 126 | 42 | 84 | 0 |
| 7 | 462 | 132 | 330 | 0 |
| 8 | 1716 | 429 | 1287 | 0 |
| 9 | 6435 | 1430 | 5005 | 0 |
| 10 | 24310 | 4862 | 19448 | 0 |
| 11 (referee) | 92378 | 16796 | 75582 | 0 |

(n = 11 reported by the referee: `11 11`, 7m22s, per-vertex [4862 1430 858 660 588 588 660 858 1430 4862], "violations: none"; I did not rerun it. Gated = LNA + not, from the saved output; they equal the E-187 table. Consequence: the average of 2 LNA neighbours per LNA (E-182) is now right + left both proved, 2 C_{n-1}.)

Paper comparison (arXiv:2112.08129, cheap part only). The arXiv API shows it is the human author's paper "A combinatorial procedure for tilting mutation" (Fosse, v2); research/FINDINGS.md l. 2599 already records that `procedure.py` is steps 1-7 of it and EXPERIMENTS.md l. 3734 that it "is the procedure itself". I did not re-read its text against steps 3-7 as used in E-187; that remains the scholar's item (LaTeX is not in `sources/`; arXiv e-print available under the standing permission). So the proof is still "relative to the repo's reading".

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/070/theorist_leftcheck.py 4 10   # 2 m 03 s; last line "violations: none"
```

## Prior record

E-187: right mutation proved, "left mutation by duality, not re-proved"; E-182 verified the left mirror on all LNAs n = 5..10 (0 mismatches) but only as an observation. grep of `research/` for "slides ... right" / left rule: none stated as a theorem. New here: the left rule in coordinates, its proof by reduction to E-187 via `_opposite`, and the direct machine check of the repo's left mutation including (C). It closes the "left mutation" item of T7; it does not touch two bystanders, k >= 4, the n = 6 non-monomial remark or the paper comparison.

## Code changed

None. New script `workshop/rounds/070/theorist_leftcheck.py`; no tests touched.

## Next

Scholar: fetch the e-print of 2112.08129 and compare steps 3-7 (beta deleted, alpha flipped, rbar = r/alpha, step 7 kernel) with E-187's use. Theorist: two-bystander case via the same step bookkeeping on a non-LNA intermediate (the result of one mutation is the input of the next: a composite right-then-left at v, v+1 is the natural first window).
