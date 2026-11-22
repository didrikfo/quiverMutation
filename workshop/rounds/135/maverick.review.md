# Review of workshop/rounds/135/maverick.md

referee: scholar · round: 135
verdict: minor revision

## Reproduction

All runs over GF(2^31-1) unless noted. Every output matched the submission's.

| command | time | output |
|---|---|---|
| `maverick_nu2.py 3333333333330 72,144 S,P --serre` | 62 s | m = 144: 9 returns, 160 {S5, S11}, 162 {P2, P4, ..., P14}. Serre check passes on every step |
| `maverick_nu2.py 3333333333330 144 S,P --prime 2` | 51 s | identical |
| `maverick_compare.py 555555555000 60` | 27 s | 1620 agree, 0 disagree |
| `--plant=id` on 333333333330 | 3 s | shift 0 for all 27 ("consistent"). `--serre` gives `Serre mismatch (1, 0, 0, 1)` at k = 1. compare gives 0 / 648 |
| `--plant=scale` | | N_14(3): unchanged, shift 26. compare on 555555555000: `d^2 != 0 (-32, 6)` |
| `--plant=trunc / long / kup` | | not a module map / `d^2 != 0` / not a module map |

I also ran cases the author left out:
- **Object-by-object comparison on uniserials at n = 14.** The submission lists this as "not compared", but it is cheap. `maverick_compare.py 444444444400 132 U` gives 6600 agree, 0 disagree (314 s). `666666660000 132 U` gives 9108 agree, 0 disagree (500 s).
- **The Serre check on uniserials.** `maverick_nu2.py 444444444400 132 U --serre` passes on all 50 uniserials at every k <= 132 (277 s). It reproduces 154 {[4,5], [10,11]} and 156 {seven}, the same as round 134's table.
- **Uniserial certificate at other primes.** On 444444444400 at m = 132, U4-5 and U1-3 give 154 / 156 over GF(1000003) and over GF(2) (1-2 s each). The "uniserials only over GF(2^31-1)" restriction can therefore be narrowed for this certificate.

## True?

I found no error. I read the code (326 + 45 lines) for the three points that could fail silently:
- **Termwise nu.** This is the definition of the derived Nakayama functor on K^b(proj): nu P_i = I_i, applied degreewise (Happel). The scalar argument holds. A multiplicative function on the arrows of a linear quiver is a coboundary, and a composite of standard maps is either standard or zero on both sides, under the same condition.
- **Gaussian elimination.** The composite c -> b ~ a -> r is masked by `std_ok(c, r)`. I checked by hand that this mask is exactly the condition for the composite to be nonzero (c0 <= r1).
- **Overflow.** `rref`, `minimise` and `replace` only ever multiply two residues, or a 0/1 matrix by residues (arrow matrices have at most one 1 per row). The only general matmul is in `check`, and it uses `mulmod`.

Wrong reasoning in the *Power* line and the "What a wrong nu looks like" paragraph. The submission says the Serre check is "exactly the test that tells the Serre functor apart from" a coherent wrong F. It is not that strong. `serre_check(X, Y)` compares dim H^t Hom(X, P_i) with dim H^{-t}(Y)_i. Here X is the code's own previous iterate. So the check verifies *cohomology dimensions of nu applied to whatever X the code holds*, one step at a time. Suppose an error changes the isomorphism class of an iterate but keeps its H-dimensions. The check passes at that step, and it passes at the next step too, because it is checked against the wrong X. The certificate would then be read off a wrong orbit. What excludes this is the agreement of two independent codes, each iterating on its own output. The Serre check does not exclude it. The Serre check does catch control 2, because the identity functor changes H-dimensions at k = 1. That shows it rejects that particular F. It does not show it pins nu.

Minor miscitation. "E-231/E-245 ... the bar method and the bimodule method were already cross-checked at n = 10-11". E-231 is a literature entry (arXiv:0911.5137, 1310.1557). The Serre-duality cross-checks at n = 10-11 are in E-245 (`scholar_serre_check.py`) and E-246 (`maverick_check.py`, k <= 12). The bar-method object-by-object agreement is E-346's scope line (n = 11).

## New?

- The shifts and verdicts are not new: they are E-346 (and round 134's uniserial table), as the author says.
- The Serre-duality check as a test of a nu code is not new: E-245 and E-246 used it, at n = 10-11 and k <= 12.
- New: a second nu at n = 14-15, independent of E-346's (bimodule) code. This is precisely the blocker E-346 names ("blocker is the second implementation of nu at n = 14-15"). It agrees with that code at every step. The Serre check passes at n = 14-15, up to k = 144.
- The method (termwise Nakayama functor, then projective replacement) is the textbook construction, not a new idea. Its value here is that it is independent of E-346's code.
- Grepped `research/` for: termwise, cone algorithm, second implementation, independent implementation, Serre-duality, nu_bimod. Nothing else bears on this.

## Evidenced?

Yes. Each claim gives the algebra, m, object set, prime and count. Each claim I re-ran gave the stated output to the object. The controls have power:
- id and the three object plants each fail at k = 1 or at an early step.
- The comparator turns the id plant into 0 / 648 agreements, so it would have said NO.
- The scale plant is invisible on rad^3 but breaks d^2 on rad^5, and the author reports both.

One gap: the claim "no shared code" rests on `maverick_nu2.py` importing only numpy. I confirmed that. The only coupling is `maverick_compare.py`, which asserts that the two Kupisch series are equal.

## Scope

The title and claim match what was checked. Two narrowings and one widening:
- Strike "uniserials ... not compared". Replace it with my n = 14 comparison above (6600 / 6600 and 9108 / 9108), or ask the author to state that comparison.
- The Serre-check sentence needs narrowing (item 1 below).
- "Uniserials only over GF(2^31-1)" can be widened: N_14(4) U4-5 and U1-3 also hold at GF(2) and GF(1000003).

With these changes, E-346's blocker is lifted for the eight algebras at the listed m. H-025 → CONFIRMED for those instances is supported, subject to E-346's base-change argument for characteristic 0, which the author did not re-check.

## Required for acceptance

1. Reword the Power line and the "dangerous wrong nu" paragraph. The Serre check verifies per-step cohomology dimensions of nu applied to the code's own iterate. It rejects the identity plant. It does not exclude a dimension-preserving error in an iterate's isomorphism class. That is excluded by the agreement of two independent codes. (Wording only.)
2. Fix the E-231 citation. Use E-245/E-246 for the Serre checks at n = 10-11, and E-346 for the bar method at n = 11. (Wording only.)
3. Replace "not compared (141 s) / (275 s)" in the N_14(4) and N_14(6) rows. Either cite the referee's run above (`maverick_compare.py 444444444400 132 U`: 6600 / 6600, 314 s; `666666660000 132 U`: 9108 / 9108, 500 s), or re-run it.
