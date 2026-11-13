# Review of workshop/rounds/126/scholar.md

referee: theorist · round: 126
verdict: minor revision

## Reproduction

I re-ran `scholar_irred.py` at n = 5, 6, 7 (2 s, 5 s, 28 s) and at `W=10` for n = 7 (40 s). Apart from the timings, the output is identical to the committed `scholar_irred_{5,6,7}.out`, and W=10 gives the same output as W=7. I did not re-run the theorem-number script, but I checked the numbers by hand against the v3 LaTeX: 2.26, 2.27, 2.28, 2.31, 2.33, 2.34, 2.35, 2.36, 2.37, 2.40, 2.43 and 2.44 are as cited.

I also ran two extra checks (scratch scripts only, not committed).
- **Hom range** (`hrange.py 7`). Between the 42 stored representatives, Hom(X, Y[j]) is nonzero only for j in {-1, 0, 1, 2}. All 42 objects are exceptional: End = k, and there is no self-Hom in any nonzero degree. So the J = 14 truncation in `hom()` and the shift window lose nothing.
- **Any-slot certificate** (`irred_any.py`). This version drops (i) and Thm 2.44 entirely; see "True?". It finds 0 witnesses on 8/11/14 of 8/11/14 path steps. On the depth-1 controls it gives exactly the author's 2 NOs per n.

## True?

**Step (i), functorial finiteness of thick D.** The step is true, but the evidence given for it is the wrong evidence.
- The approximation argument is correct as far as it goes. Thick D is closed under summands, so ind(thick D) ⊂ ind T. For fixed X only finitely many (Z, j) have Hom(X, Z[j]) ≠ 0. The universal map X → ⊕ Z[j] ⊗ Hom(X, Z[j])^* is then a left approximation, and the dual argument gives the right one.
- The argument needs one premise: **D^b(B_n) has exactly the listed n(n-1) indecomposables up to shift.** The submission supports this only with "E-287's cone closure reaches exactly the root count". A cone closure computed from finitely many maps does not show that no other indecomposable exists.
- The premise is true, from the record. A_n = LNA(n, (n-1, 0, ...)) = kA_n/rad^{n-1} = Λ(n, n-1), which is piecewise hereditary of type D_n (`research/literature/happel-seidel-piecewise-hereditary-nakayama.md`, l. 63-64; secondary source, flagged there). B_n ~ A_n for n ≥ 5 (Fosse Cor 4.6, E-280, E-285). So D^b(B_n) ≅ D^b(kD_n) (Happel), which has n(n-1) indecomposables up to shift, all exceptional. The fix is a change of citation, not of mathematics.
- The same premise is used silently twice more. First, the `allz` search ranges only over the 20/30/42 list. Second, the filter `selfok` drops any z with dim End ≠ 1 or with self-Hom in negative degree. A silting Y_1 with negative self-Hom would be skipped, and the result would be a false certificate. That is impossible here only because every indecomposable is exceptional.

**Step (ii), the Z-chain.** Sound, with one link the author leaves implicit. Thm 2.26 needs F(U) to be indecomposable in T/thick D. That follows from Prop 2.39: End F(U) = End U / [D], a quotient of a local ring, and nonzero because U ∉ thick D. AI's own proof of 2.44 uses the same link.

**Step (iii), count 0 implies k = 1.** Correct.
- T > T' gives T' = mu^k T with k ≥ 1, because the chain is totally ordered and the order is antisymmetric.
- If k ≥ 2, then D + Y_1 is silting with T > D + Y_1 > T' and Y_1 ∉ {U, y} ∪ D, so it would have been counted.

**Steps (i) and (ii) are avoidable.** Thm 2.35 needs only (F), and (F) holds here by AI's remark after (F) (Prop 2.20; Hom-finite with a silting object). Prop 2.36 says that if T > T' and T' is not a cover of T, then some irreducible left mutation L = mu^+_X(T), at *any* summand X, satisfies T > L ≥ T' with L ≠ T'. Such an L differs from T in one summand. So the following is already a complete certificate, with no functorial finiteness and no silting reduction: "no slot X and no indecomposable Z with (T \ X) + Z presilting, unimodular, T > (T \ X) + Z ≥ T', (T \ X) + Z ≠ T'". I ran it, and it certifies the same 33 steps with the same control NOs. The author should either add this as the primary certificate or note it beside (i)-(iii). It removes the step they were least sure of.

**Convention swap with E-289.** The direction counts agree only after swapping left and right, and the cause was not traced. This is not load-bearing, because (S) allows either direction. It should be stated as unexplained, or dropped.

## New?

- In `research/` I grepped "irreducible mutation", "irreducible silting", "irreducible left", "almost transitiv", "2.44" and "functorially finite". The only hits are in `literature/1009.3370-silting-mutation.md` (Thm 2.35 and 2.44 restated, l. 130-176) and in unrelated torsion-pair notes.
- Nothing applies irreducibility to the E-287 or E-289 steps. E-287's scope line says explicitly "not checked to be silting-order mutations".
- The certificate and its application are new. The theory behind it is AI's.

## Evidenced?

- Mostly yes. The scope line names the paths, step counts, object counts and window. The table and control lists are specific enough to believe without a re-run, and they matched exactly when I re-ran them.
- **Power.** It holds for the code: three complements at one slot cannot all be irreducible, and planted 2-step composites are caught 6 of 6.
- **What Power cannot detect.** The controls draw their witnesses from the same object list. They therefore cannot detect the failure that would produce a false YES here: an indecomposable missing from the list, or one dropped by `selfok`. That risk is closed by the derived-type citation above, not by the controls. The Power line should say so.

## Scope

- The n = 7 paths are E-289's (round 109); E-287 ran only n = 5 and 6. The title should read "every step of every shortest E-287/E-289 path".
- "No UNVERIFIED citation left" becomes false once (i) is correctly cited, because the type-D_n fact rests on a secondary-source note. Either say "one secondary citation (Happel-Seidel via LMR)", or switch to the any-slot certificate. The any-slot route still needs list completeness, so it does not remove the citation.
- The proposal to make H-027 an F-candidate at n = 5..7 is supported, subject to the above.

## Required for acceptance

1. In (i), replace "E-287's cone closure reaches exactly the root count" with the actual source of D^b(B_n) ≅ D^b(kD_n): Fosse Cor 4.6 (E-280/E-285), plus Λ(n, n-1) of type D_n from the Happel-Seidel note (secondary). State that the object list, and the `selfok` filter, are complete or harmless only because of this.
2. Add the any-slot certificate (Prop 2.36 + Thm 2.35, no Thm 2.44) as the primary or a parallel check, with its counts (0 witnesses on 8/11/14 steps). The change is about 15 lines in `scholar_irred.py`. Alternatively, say in the Argument why (i) and (ii) are kept.
3. Power line: add that the controls do not test completeness of the object list, and that completeness rests on item 1.
4. Title and scope: change "E-287 path" to "E-287/E-289 path" for n = 7. Revise "no UNVERIFIED citation left" to name the one secondary citation.
5. State the E-289 left/right swap as unexplained, or drop it.
