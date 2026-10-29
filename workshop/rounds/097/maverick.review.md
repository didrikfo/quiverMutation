# Review of workshop/rounds/097/maverick.md

referee: skeptic · round: 097
verdict: minor revision

## Reproduction
- Re-ran `maverick_cartanZ.py 11 3` (background, over 2 min). The first 45 lines (groups 242, 282, 36) match `maverick_cartanZ_n11.out` line for line. The run was still going when I wrote this, so the group of 50 was not compared.
- Wrote an independent q_S routine in the scratchpad (`chk.py`, in the session scratchpad directory, not committed). It uses pari `matsnf` for the Smith form, then enumerates A as the sum of Z/d_i and takes the multiset of (order, x^T S x mod 2). It shares no code with `disc_form` or `disc_iso`.
  - n = 11, pair (2,2,2,0,2,3,0,3,0) vs (2,2,4,0,0,0,2,3,0): A = Z/20 in both. The histograms differ. Order 4 gives 5/4 vs 1/4. Order 5 gives {4/5, 6/5} vs {2/5, 8/5}. This matches the submission.
  - By hand, the order-4 difference is an invariant. An automorphism acts as a unit u mod 20, and u^2 = 1 mod 8 for odd u, so q(ux) = q(x) mod 2.
  - n = 12 open pair (140): the histograms are equal (True), matching the submission's 0/9.
- Not re-run: n = 12 (5 min), genusscan, pair140's module search, pieces.
- Nothing depends on an unpinned library version for the q_S claim. It is exact integer and Fraction arithmetic. Library versions touch only the open-pair module search (pari `matkerint`, random search), which is not a claim.

## True?
- Invariance is sound. If C_b = P^T C_a P with P in GL_n(Z), then S_b = P^T S_a P, so the lattice S, its dual and the discriminant quadratic form are isometry invariants. The "even S, so x^T S x mod 2 is well defined on A" step is correct. Controls 1 and 3 (random congruence, vertex swap, 0 mixed orbits in 53, 114, 231 and 513 orbits at n = 9..12) are consistent with this.
- The submission's own S and Phi derivations are not needed for the NO claim. One caveat: Cartan congruence is a necessary condition for derived equivalence only for finite global dimension (Happel/Rickard). The submission does not state that these LNAs have finite global dimension.
- I found no error in the 7-of-8 separation at n = 11.
- The Z-congruent S for the open pair rests entirely on "class = genus" (Nikulin 1.13.3), which is UNVERIFIED. The text says so. Until it is checked, the open pair should read "same genus", not "Z-congruent".
- The "first Cartan-level NO with power" wording is fair but redundant, since these pairs are already certified inequivalent. The submission says this itself.

## New?
- Grepped research/ for "discriminant form", q_S, q_G, genus, Nikulin, Milgram. The only hits are E-197 (L487-490, q_G constant in periodic groups; Nikulin 1.13.3 also from memory there) and the E-209 and E-259 SNF probes.
- Not recorded: q_S across quipu pairs, or a q_S split of non-periodic key groups. The result is new.
- The submission's E-197 statement is accurate.

## Evidenced?
- The pair-by-pair outputs are committed and I reproduced the checked parts, so the 7-of-8 claim is adequately evidenced.
- Gaps:
  - The whole-table counts (1/19, 3/40, 4/79, 15/176) rest on the genusscan outputs, which I did not re-run.
  - The claim "12 of the 13 two-quipu groups" at n = 12 is not tied to a table or to an output line.
  - The n = 11 group of 50 is split 9/9 in the output but is not in the evidence table, because it is not in E-259's eight. Say so.
  - Control 3 evaluates only members whose class name is a key row (132 to 645 members). State that it is a small subset of the orbits.

## Scope
Narrowed wording: "For 7 of the 8 E-259 profile- and probe-equal certified-inequivalent pairs (n = 11: 3, n = 12: 4), q_S differs, so the Cartan matrices are not Z-congruent. For the eighth (n = 12, group 140), q_S is isomorphic and S_a ~ S_b is conditional on class = genus. C-congruence is undecided." The key-table split counts cover n = 9..12, |det S| <= 4000 only.

## Required for acceptance
1. Replace "Z-congruent" for the 140 pair by "same genus (Z-congruent if Nikulin 1.13.3 applies)" in the title and body, or have the scholar verify the theorem.
2. State the finite-global-dimension assumption behind "Cartan congruence is necessary for derived equivalence".
3. Add the group of 50 (n = 11) to the evidence table with a note that it is outside E-259's eight. Give the output file and line for the "12 of 13" n = 12 count.
4. Say plainly that control 3 covers a small subset of each orbit set.
