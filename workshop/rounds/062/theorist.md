# Proof, relative to the repo's `procedure` steps 1-7 (not the paper), of the single-mutation LNA rule (E-182), with both "only if" cases; Catalan factorisation (C) proved and machine-checked per vertex for n = 4..10

author: theorist · round: 062 · kind: proof
thread: T7 · bears on: E-182, E-153, F-020, H-010
scope: right mutation of an LNA (monomial line algebra, n vertices, admissible row r) at one vertex v, as defined by `procedure` steps 1-7 followed by reduction. Proof is by hand, for all n, against the steps as the repo's `procedure` docstring states them (not against the paper); case analysis and the per-vertex count (C) machine-checked for n = 4..10 (all rows, all v). Left mutation by the dual (not re-proved). Not covered: LNAs with the procedure read from the paper itself (I did not re-read arXiv:2112.08129).

## Response to referee

Verdict minor revision; all four required items answered.

1. Per-vertex (C) check: DONE. `theorist_t7check.py` now counts, for each n and v = 1..n-1, the LNA results at v and compares with C_{v-1} C_{n-1-v}; mismatch is a violation. Re-run `timeout 10m .venv/bin/python workshop/rounds/062/theorist_t7check.py 4 10`: "violations: none" for n = 4..10 (totals unchanged, table below). Per-vertex counts printed, e.g. n = 8: 132 42 28 25 28 42 132; n = 10: 1430 429 264 210 196 210 264 429 1430. So (C) is now machine-checked per vertex, not only in total.
2. Title and scope now say "relative to the repo's `procedure` steps 1-7, not the paper". Claim narrowed accordingly: proof for all n relative to that reading; machine-checked n = 4..10.
3. Generation of the result ideal (step 7 / weakest step): argued, see the new paragraph "No other kernel" below. It is still relative to the reading that the result algebra is the line quiver of step 5 modulo the monomial ideal of vanishing paths. It does not rest on the machine run, but the machine run (image = r') remains the only check of that reading.
4. Noted in (F): E-182 observes that for n = 6 the non-LNA results (r_v >= 3) can carry non-monomial relations. I did not re-check that observation. (F) uses only the arrow count and the out-degree of v+1, so it does not depend on whether the relations are monomial; "not an LNA" is all that (F) asserts.

## Claim

Write r_j for the length of the relation starting at j (0 = none), e(j) = j + r_j for its end vertex. For v in 1..n:

(G) the gate `isMutable` passes iff v <= n-1 and no relation ends at v+1.
(L) If the gate passes, the reduced result is an LNA iff v = n-1 or r_v = 2. Then its row is r with (r_{v-1}, r_v) := (2, 0) (for v = 1 only r_1 := 0).
(F) If the gate passes and v <= n-2 with r_v = 0, the result has n-1 arrows and a vertex of out-degree 2; if r_v >= 3 it has n arrows (and out-degree 2): neither is a line.
(C) The number of rows r with an LNA result at v is C_{v-1} C_{n-1-v} (sum over v = 1..n-1 is C_{n-1}).

E-182's "no earlier relation covers arrow v" is two different things. A covering relation j < v with e(j) = v+1 is exactly the gate failing (no mutation, so no LNA). A covering relation with e(j) >= v+2 forces r_v != 2 by increasing ends, so it is already excluded by "r_v = 2". So the rule is: **gate passes and (r_v = 2 or v = n-1)**; "covered" adds nothing beyond that. Refuted by: a row, a v with gate passing and r_v = 2 whose result is not the stated row.

## Evidence (proof)

Notation: beta = a_{v-1}, alpha = a_v, gamma = a_{v+1} (arrows a_j : j -> j+1). Relations are paths (monomial).

(G) Gate: the only out-arrow of v is alpha (needs v <= n-1). It fails iff some nonzero path ending at v has zero composite with alpha. A nonzero path p: j -> v with p alpha zero contains, as it is nonzero but p alpha is zero, a relation as a suffix of p alpha, i.e. a relation j' ~> v+1 starting at j' <= v-1 (every relation has >= 2 arrows) -- conversely a relation ending at v+1 gives p = its proper prefix. Hence: fails iff a relation ends at v+1. (Checked by machine.)

(L), "if", case r_v = 2, v in 2..n-2, gate passing. No relation ends at v+1, so every relation starting before v ends at <= v (a relation j < v with e(j) >= v+2 would have e(j) >= v+2 = e(v), contradicting strictly increasing ends; e(j) = v+1 is excluded by the gate). In particular r_{v-1} = 0. Steps:
 - 1: composite c = alpha beta : v-1 -> v+1. 2: alpha* : v+1 -> v*. 3: the relation v ~> v+2 (alpha gamma) gives rbar : v* -> v+2. (beta is deleted, alpha is flipped.)
 - 4: the single beta gives c alpha* = 0 (path v-1 -> v+1 -> v*).
 - 5: rbar alpha* = r/alpha = gamma. Reduction eliminates gamma. Quiver: ..., v-1 -c-> v+1 -alpha*-> v* -rbar-> v+2, ... : a line; positions: old v+1 sits at position v, v* at position v+1, the rest unmoved.
 - 6: a relation j ~> v (the only ones that meet v, by the first paragraph) becomes j ... v-1 -c-> v+1: same length, same start, same end position v. Relations ending <= v-1 are carried unchanged.
 - Relations starting at s >= v+1 are carried; the one at s = v+1 (if any) starts with gamma = rbar alpha*, so it is alpha* times a path out of v*.
 - 7: the paths out of v* to a target t are unique: rbar a_{v+2} ... a_{t-1}. It is in the kernel iff gamma a_{v+2} ... a_{t-1} lies in I iff the path v+1 -> t contains a relation, iff t >= e* := the smallest end among relations starting at s >= v+1. The kernel is a monomial: one relation v* ~> e*. If the relation with end e* starts at s = v+1, this is the relation of position v+1 with r'_{v+1} = e* - (v+1) = r_{v+1}, and the carried relation at position v (alpha* times it) is redundant, so r'_v = 0. If s* >= v+2, the step-7 relation contains the carried relation s* ~> e* and is redundant, and r'_{v+1} = 0 = r_{v+1}; r'_v = 0 again.
 Result generators: left relations unchanged, c alpha* (start v-1, length 2), the right relations unchanged. Starts and ends are strictly increasing (left ends <= v < v+1 = end of the new one < every right end >= v+3), so the result is admissible, with row r' = r except (r'_{v-1}, r'_v) = (2, 0). (v = 1: no beta, no c, no step 4, so r'_1 = 0 only. v = n-1: no relation starts at v, so no rbar and no step 5/7; v* is a sink; step 4 gives r'_{n-2} = 2, step 6 as above.)
 **No other kernel (answer to referee 3).** In the line of step 5 there is at most one path between two vertices, so the ideal is monomial and is generated by the minimal vanishing paths; it suffices to decide which paths vanish. Classify a path P by where it sits (positions as above, c: v-1 -> v+1, alpha*: v+1 -> v*): (a) P inside 1..v-1 or inside v+2..n (old arrows): zero iff it contains an old relation, since the steps leave those arrows and relations unchanged (relations meeting v are step 6, below); (b) P ends at v+1 through c, start j <= v-2: P = (j..v-1) c, zero iff (j..v-1) beta alpha contains a relation; as none ends at v+1, this holds iff j..v contains one, i.e. iff the carried step-6 relation (or an old one) lies in it; (c) P contains c alpha*: zero by step 4, and is the generator of length 2, so every longer path through it is redundant; (d) P starts at v+1 or v* and goes out through alpha* rbar a_{v+2}...: unique per target, handled by step 7 (zero iff t >= e*); (e) P starts at v+1 via alpha* only (target v*): alpha is not in a relation of I (r_v = 2 starts at v, so alpha alone is nonzero), so nonzero. Every path is in one class, and in each the vanishing set is the one stated; the generators are exactly the minimal ones listed. **Weakest step:** the reading that the reduced algebra's ideal is the set of vanishing line paths as determined by steps 4-7 (a statement about the procedure's reduction, not re-derived from the paper); the machine run (image = r' in every case) agrees with it but is not independent of it.

(F), "only if". Reduction can only remove an arrow lying in a relation that is a single arrow modulo others. 
 - r_v = 0 (v <= n-2): no rbar, no step 5, no step 7. The relations are monomial of length >= 2 (a carried relation through beta alpha has length >= 3 since the length-2 case beta alpha is a relation ending at v+1, gated out; c alpha* has length 2). Nothing eliminates gamma, so v+1 has out-arrows alpha* and gamma: out-degree 2. (v = 1: alpha* : 2 -> 1*, a_2 : 2 -> 3, same.) Arrows: (n-1) - 2 + 2 = n-1 (v >= 2), a tree with a branch vertex, not a line.
 - r_v >= 3: step 5 gives rbar alpha* = a_{v+1} ... a_{v+r-1}, a path of >= 2 arrows, so no arrow is eliminated. Arrows (n-1) - 2 + 3 = n for v >= 2 (n-1 - 1 + 2 = n for v = 1) > n-1. Not a line.
 Both agree with the machine counts below. E-182 notes that for n = 6 the non-LNA results at r_v >= 3 can have non-monomial relations; (F) does not use this (only arrow count and out-degree).

(C) The map r -> r' is injective with inverse "(r'_{v-1}, r'_v) = (2, 0) -> (0, 2)". Its image is the set of admissible rows with r'_{v-1} = 2 and r'_v = 0 (the preimage has r_v = 2, no relation ending at v+1: left relations end <= v, right ones end >= v+3, r_v ends v+2). Such a row is a free pair: an LNA on vertices 1..v (relations starting <= v-2, ends <= v; C_{v-1} of them) and an LNA on vertices v+1..n (C_{n-v-1}); they do not interact since the middle relation has end v+1 and every left end <= v < v+1 < every right end. v = n-1: C_{n-2} C_0; v = 1: C_0 C_{n-2}.

Machine check (`theorist_t7check.py`, n = 4..10, all rows, all v): gate = (G) in every case; LNA iff (v = n-1 or r_v = 2); image = r'; (F) edge counts and out-degree; (C) per-vertex count of LNA results = C_{v-1} C_{n-1-v}. Gated cases and LNA results (count LNA = C_{n-1} as (C) says):

| n | gated events | LNA | not LNA | violations |
|---|---|---|---|---|
| 4 | 10 | 5 | 5 | 0 |
| 5 | 35 | 14 | 21 | 0 |
| 6 | 126 | 42 | 84 | 0 |
| 7 | 462 | 132 | 330 | 0 |
| 8 | 1716 | 429 | 1287 | 0 |
| 9 | 6435 | 1430 | 5005 | 0 |
| 10 | 24310 | 4862 | 19448 | 0 |

Total LNA results = C_{n-1} per direction, 2 C_{n-1} with the left, as E-182. Hence the average number of LNA results per LNA is 2 (right plus left), as E-182 found, and not constant: the count per row is #{v : r_v = 2 or v = n-1, gate passing}.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/062/theorist_t7check.py 4 10   # 1 m 43 s; prints the table, "violations: none"
```

## Prior record

E-182 states the rule and the factorisation as checked, not proved ("Open: proof of the rule and of the factorisation"); E-153 the weaker locality. I grepped `research/` for a proof of the rule and for "ends at v+1"; none found. New here: the proof, the observation that "covered" reduces to the gate plus r_v = 2 (so the rule has no independent covering clause), and (C). The reading of the procedure is the repo's `procedure` docstring; not compared with the paper. Mirror image for left mutation is by `dualRule` symmetry (E-182 verified it, I did not re-prove it).

## Code changed

None. New script `workshop/rounds/062/theorist_t7check.py` only; no tests touched.

## Next

Scholar: compare steps 3-7 as used above (beta deleted, alpha flipped, r/alpha) with the paper's statement, since the proof is exactly as good as that reading. Then this proof is a template for two mutations (the E-152 / H-010 question: a bystander that makes the window leave the LNA class): the same step-by-step bookkeeping on a non-LNA intermediate.
