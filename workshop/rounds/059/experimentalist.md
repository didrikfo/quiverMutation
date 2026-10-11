# For all LNAs of length n <= 10 (checked, not proved) a right LNA mutation stays an LNA iff it slides an uncovered two-arrow relation one place left; this explains E-153's "2 per LNA" (an average, not a per-LNA fact); the listed k = 3, one-bystander cap-cut cells of E-152 are inert

author: experimentalist · round: 059 · kind: result
thread: T7 · bears on: H-010, E-152, E-153, F-020, F-022
scope: Part A: every LNA of length n = 5..10, all single mutations, rule checked on all of them (n = 10: 4862 LNAs, 48 620 events); exact on this finite range, proof only sketched. Part B (listed cells only): pair (8:3)(9:3) + one bystander, k = 3, margin 2, 17 rows, capped at ~25 min wall each (rows beyond the listed ones were killed), nothing for k > 3, two bystanders, margin > 2.

## Response to referee

1. Title: the iff-claim now carries "for all LNAs of length n <= 10 (checked, not proved)". Done.
2. Title/claim/scope: Part B is now said to cover only the listed k = 3, one-bystander cells (title "listed ... cells"; claim B-paragraph unchanged, scope line adds "listed cells only"). Done.
3. Part C: relabelled a conjecture supported by per-v counts at n = 7, 9, 10 (not derived). Done.
4. Independent n = 6 failure-mode check (new `experimentalist_fail.py 6`, 3 s). Of 42 LNAs x 6 vertices = 252 right mutations: 42 are LNA results (exactly the rule's), 126 are not possible at all (42 at v = n, 84 at a left-covered arrow), and 84 are possible but not an LNA: r_v >= 3: 20, every one a quiver with 6 arrows on 6 vertices (one extra arrow, no parallel pair, no cycle) and a non-monomial relation set; left-covered arrow: 36, of which 28 give a 5-arrow line quiver whose relation set is not an LNA and 8 an extra arrow with non-monomial relations; r_v = 0 (v < n-1): 28, all 5-arrow acyclic non-LNA. So the earlier "parallel or extra arrow" guess was wrong in detail: no parallel arrows occur. This confirms the failure side at n = 6 (no non-LNA is mislabelled); it does not add an independent check of "only v-1, v change" beyond the n <= 10 image check. Done at n = 6 only.
5. Not done: re-running n = 10 and Part B (7 min and 25 min per row; reviewer did not require).

## Claim

**A. E-153's headline "of the n single mutations of an LNA exactly 2 give an LNA" is wrong per LNA; it is an average.** Number of single mutations (right + left) giving an LNA, per LNA: n = 5: 0..4 (counts of (#R,#L): 7 different pairs); n = 8: 0..7; n = 10: 0..9 (25 different pairs). The identity that is true is *total* LNA results = 2 * Catalan(n-1) (right and left each C_{n-1}), which is what E-153's 20/70/252/.../48 620 vs 10/28/84/.../9 724 shows. (The claim in the title needs rewording; the data in E-153 are fine.)

**B. Exact local rule, n <= 10.** Write r_j for the number of arrows of the zero relation starting at arrow j (0 = none; relation (j, m) covers arrows j..j+m-1, arrow v = v -> v+1). The right mutation at vertex v (1 <= v <= n-1) of an LNA r is an LNA iff
(i) no relation starting before v covers arrow v (in particular r_{v-1} = 0), and
(ii) r_v = 2, or v = n-1 (a "virtual" length-2 relation sticking out of the line).
The result is r with (r_{v-1}, r_v) := (2, 0): the length-2 relation is slid one place left (at v = 1 it is deleted, at v = n-1 it is created). Left mutations are the mirror image: L(v) of r = R(n+1-v) of the mirrored vector. No other vertex, no relation of length >= 3, ever gives an LNA by one mutation. So the changed starts lie in {v-1, v}, sharper than E-153's {v-2, v-1, v}; this is the single-mutation form of F-020's "lone two-arrow slide".

**C. The count (conjecture, supported by per-v counts at n = 7, 9, 10; not derived).** R(v) is conjectured to hold for C_{v-1} * C_{n-1-v} LNAs (left part: any LNA on arrows 1..v-1, relation (v,2), right part any LNA on the later vertices), and sum_v C_{v-1}C_{n-1-v} = C_{n-1} is the Catalan recurrence. Checked against the computed per-v counts at n = 9 (429, 132, 84, 70, 70, 84, 132, 429) and n = 10 (1430, 429, 264, 210, 196, 210, 264, 429, 1430); E-153's pooled n = 7 list 42, 14, 10, 10, 14, 42 is the same convolution. The "tail is an arbitrary LNA" factorisation is the unproved step; what was checked is the per-v counts and totals.

**Refutes if:** an LNA of any length where a mutation at v with r_v >= 3, or with arrow v covered from the left, gives an LNA; none at n <= 10.

## Evidence

| n | LNAs | R events | rule mismatches (R) | image mismatches | left-mirror mismatches | distinct (#R,#L) |
|---|---|---|---|---|---|---|
| 5 | 14 | 14 | 0 | 0 | 0 | 7 |
| 6 | 42 | 42 | 0 | 0 | 0 | 11 |
| 7 | 132 | 132 | 0 | 0 | 0 | 12 |
| 8 | 429 | 429 | 0 | 0 | 0 | 16 |
| 9 | 1430 | 1430 | 0 | 0 | 0 | 17 |
| 10 | 4862 | 4862 | 0 | 0 | 0 | 23 |

(R events = number of (LNA, v) with LNA result, from `mutationIsPossibleAtVertex` + `quiverMutationAtVertices`; left mirror: offset c = n+1 is the only one of 5 tried with 0 mismatches. The distinct-pair column is from the printed (#R,#L) dictionaries.) Per-LNA example, n = 5: [2,0,0] has R at v = 1, 4 and L at 3; [3,3,0] and [4,0,0] have none.

**Part B, E-152 cap-cut rows (k = 3, margin 2, base overlap 2).** Right, bystander (t:m): m = 3 at g = 2, 3 (nothing shared); m = 4 at g = -1 (1 shared), 0, 1. Left, bystander ending g arrows before the pair: m = 2, g = -1..3; m = 3, g = -1..2; m = 4, g = -1..1. All rows: reached LNAs 1..12, min overlap 2, lowered 0 (rows 124-387 s each under 7-way CPU sharing). So no placement sharing <= 1 arrow lowers the overlap at k = 3 in the rows now run, in agreement with E-152's 10 right + 2 left earlier; with this the formerly capped m = 4 rows (g >= -1) and m = 3, g >= 2 are partly filled (not right m = 3 g >= 4, m = 4 g >= 2, left m = 3 g >= 3, m = 4 g >= 2). It is a null result for those cells only, and inert is the expected outcome. I did not run the two-bystander or k >= 4 cells.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/059/experimentalist_which.py 8     # 7 s; n = 9: ~1.5 min, n = 10: ~7 min alone (needs a background job above 10 min on a loaded machine)
.venv/bin/python workshop/rounds/059/experimentalist_rule.py 8                   # needs the JSON from the line above; n = 9 JSON is committed (132 KB), n = 10 JSON (495 KB) deleted, rerun to regenerate
timeout 25m .venv/bin/python workshop/rounds/059/experimentalist_capcut.py right 3 4 -1 1    # --plan lists rows without running
```
Row outputs: `workshop/rounds/059/experimentalist_capcut_*.txt`.

## Prior record

E-153 records the counts and calls the "exactly 2" unexplained; its own totals (events 20..48 620, results 10..9 724) are exactly 2 * C_{n-1}, so the "exactly 2 per LNA" wording in its title is an average. The Catalan product by vertex was not noticed there (grep "Catalan" in research/ has no hit on this). F-020 (lone two-arrow slide) and F-022 contain the multi-step slide; the single-mutation form and its exact side condition (arrow not covered from the left) are, as far as grep shows, not recorded. E-152 rows: see above; no recorded result contradicted.

## Code changed

None in `quivermutation/`. New scripts in `workshop/rounds/059/`: `experimentalist_which.py`, `_rule.py`, `_capcut.py`.

## Next

- theorist: prove the rule (a mutation at v of a Nakayama algebra with r_v >= 3 or a left-covered arrow v gives a non-LNA: at r_v >= 3 an extra arrow with a non-monomial relation, otherwise a line quiver with a non-LNA relation set (n = 6 table in the response)) and the Catalan factorisation; then L1 for k steps becomes "iterated slides of length-2 relations", which is the natural route to H-010 (a longer relation never moves by a single mutation, so overlap changes need a multi-step window with non-LNA intermediates).
- chair/scribe: reword E-153's title ("exactly 2 on average; total 2 C_{n-1}"); not done by me (research/ is not mine).
- experimentalist (later): two-bystander cells and k = 4 need an OVERNIGHT proposal; cost of one k = 3 row is 2-6 min at L = 21-27, k = 4 is estimated one order more.
