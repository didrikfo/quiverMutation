# Review of workshop/rounds/101/maverick.md

referee: skeptic · round: 101
verdict: minor revision

## Reproduction

Both scripts re-run unchanged. `maverick_litiso.py` took 18 s and `maverick_litpairs.py` 2.5 s, both far under the 10-minute limit.
- litiso: 3 pairs x 2 slacks, all `congruent: True`. nsol is 5/22 for N(10,8)~N(10,7), 4/17 for N(11,9)~N(11,8) and 4/15 for N(12,10)~N(12,9). All match the author's numbers.
- litpairs: all 10 pairs have charpoly equal, Phi-modules isomorphic, and the same sig/det as in the table. q_S is True on the three non-vacuous pairs with |det S| > 1 that were run (11,9/8; 12,5/4; 12,10/9), as claimed.
- The table's "q_S YES" count (3) matches, and it is vacuous on the other 7.

## True?

Literature (re-fetched LaTeX, own scratchpad dirs):
- **Ladkani math/0610685.** Counting theorem-type environments in section 3 gives Prop 3.11 (Euler = Moebius), Def 3.12, Cor 3.13 (X ~ Y => 1_X, 1_Y Z-congruent), Lemma 3.14 and Cor 3.15. The numbers are right. The proof of 3.13 shows [F]^t 1_Y^{-1} [F] = 1_X^{-1}, and the paper notes that inverses of congruent matrices are congruent. The extension to finite-global-dimension algebras is correct and is correctly flagged as ours. Small wording point: the statement of 3.13 is about 1_X (zeta matrix), and the Euler matrix 1_X^{-1} appears in the proof. The write-up's "Euler-form matrices ... Z-congruent" is equivalent and fine.
- **Ueda 2302.02880.**
  - Cor 1.3 reads: n = p(p+1)q + p(p-1)r, l = (p+1)q + pr, A(n,l+1) ~ A(n,l), with r in (1/2)Z for p = 2.
  - I recomputed every pair: (10,8/7) is p=2, q=1, r=2; (11,9/8) is r=2.5; (12,10/9) is r=3; (9,7/6) is r=1.5; (12,7/6) is p=2, q=2, r=0; (12,5/4) is p=3, q=1, r=0. All are valid.
  - Cor 1.6 gives A(pq+1,q+1) ~ A(pq+1,p+1), which yields (9,5/3) with p=2, q=4 and (11,6/3) with p=2, q=5. Valid.
  - Ladkani 1.2 gives (10,3/6) and (12,3/7). Valid, as cited by Ueda.
- **Bessenrodt-Holm math/0504547.** Theorem `unimod` is the first theorem of section 2, so "2.2" is right. It states derived invariance of the unimodular equivalence class of C.
- **Radical-power convention.** I checked that `cartan(n, "l"*(n-l)+"0"*(l-2))` is exactly the radical-power Cartan matrix, with entries 1 for 0 <= i-j < l and 0 otherwise (lower triangular). So the pairs really are the cited algebras.
- **The patched Salem threshold (1.5 -> 1.05).**
  - In `pq` the threshold only selects which real eigenvalue is lam. I checked the real eigenvalues > 1.0001 of Phi: exactly one per algebra (1.17628, 1.23039, 1.26123), and equal across each pair.
  - The bound in `test` depends on lam only through the sampled t-window [lam^-1/2, lam^1/2]. That window is one fundamental domain of Phi^Z acting on t, and it is valid for any lam > 1. So the patch is harmless.
  - A negative t gives the same quadratic form, so the sign is not an issue.
  - The 099 code was just over-restrictive at 1.5 (it hit IndexError).
  - The bound is not stress-tested on a NO at small lam, but all YES results have nvec = nsol at slack 0.5, so the bound is tight and not a loose slack.
  - The 3 of 3 YES result is unaffected either way.
- I found no counterexample.

Overreach, not error: the headline "no independent second invariant exists for LNAs" is not established. What was checked is:
- HH^* (trivial on LNAs, F-019),
- the q-Cartan determinant (gentle only, trivial here),
- and AAG, explicitly not addressed.
It does not rule out other invariants, for example the full derived Picard/Serre data or non-gentle LNAs for AAG. Narrow to: "the invariants reviewed (HH^*, q-Cartan det, unimodular class) give no second NO on LNAs".

## New?

- The literature-note half is only partly new. `research/literature/math-0610685-sheaves-over-finite-posets.md` already has "Proposition 3.11 and Corollary 3.13" and "Lemma 3.14 and Corollary 3.15", and says these transfer verbatim. FINDINGS F-047 and the EXPERIMENTS.md sections at lines ~3110 and ~3200 already use Cor 3.13/3.15. What is new is only the observation that H-024's necessity step is exactly this, plus the finite-gldim extension.
- H-024 (open) says "no source matched" and its falsifier is E-238's 1 002 pairs or the 6 n = 12 within-quipu pairs. Both are already YES (E-269 and the 099 table).
- Nothing found for "Cor 1.3 / 1.6 radical-power pairs as Cartan controls". The 3 extra E-269 YES runs on Ueda pairs are new but low information, since they are the same method on equivalent pairs.
- The Hochschild/AAG route is F-019 / idea 22, rediscovered. The author says so.

## Evidenced?

Yes for what is claimed: specific pairs, n, tests, counts, and an honest "YES-only" power line. The weak points:
- Power. The 10 pairs cannot have returned NO unless the method was broken. They confirm no false NO on these pairs, nothing more.
- Pairs. All 10 are radical powers, which share one Phi-class. They do not exercise the mixed-digit strings behind the E-259 NOs. The author says this.
- E-269 coverage. Only 3 of 10 pairs got the exhaustive search, and 4 of the 7 not run (two Salem pairs) are untouched.

## Scope

- Title: "Necessity of Z-congruence is a published one-line argument" should read "... for incidence algebras of posets (Ladkani Cor 3.13); the same proof gives it for LNAs".
- Title: "10 pairs pass the Cartan tests" should read "pass char-poly and Phi-module tests; q_S is non-vacuous on only 3; E-269 exhaustive on 3".
- Claim 4: narrow as above.

## Required for acceptance

1. Narrow the "no independent second invariant exists" claim in the title and Claim 4 to the invariants actually reviewed.
2. State in Claim 1 that Cor 3.13 itself is about 1_X, with the Euler matrix 1_X^{-1} in the proof, and that LNA applicability is the author's extension. The one-sentence proof for any finite-gldim algebra should be given, not just asserted.
3. Add to the Prior record the existing matches (literature/math-0610685 note, F-047, E-267) so the "new" part is only the Ueda/Ladkani control pairs.
4. Add one line saying the 1.5 -> 1.05 patch only changes which eigenvalue is selected, because the bound depends on lam only through the t-window. This makes the "needs the skeptic's eye" note resolvable.
5. [next round] Add a NO-power control for the small-lam bound: a cross-quipu n = 10 or 11 pair with a Salem number below 1.5. Without it the patched search is untested for NO at those lam.
