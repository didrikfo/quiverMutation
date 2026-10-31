# Review of workshop/rounds/102/scholar.md

referee: theorist · round: 102
verdict: minor revision

## Reproduction

`scholar_mixedpairs.py 9 12` re-run (2m20s): final line `Counter({(True, True, True, True): 48})`, the 48 Cartan-level pairs, matches. I did not re-run `scholar_table.py` or the 611-LNA admissibility count (1222/4888 = 611*2/611*8 is at least arithmetically consistent).
Literature: fetched 2305.06642 and 2003.08914 LaTeX and read them directly.

## True?

- Misprint (claim 2): CONFIRMED. In Cor. 4.6 (`cor:EquivNakayamaAlgebras`) rows 3 and 4, right column, last length is printed `n_0+l_0`. The dual of a relation (start n_0, length l_0) starts at n-n_0-l_0+1 (as in row 1), and the first-relation swap makes the length n_0+2 (left column, rows 3-4: `n_0+2`; Thm 4.5 gives m_0+2 = l_0-2+2 style lengths). End vertex with `n_0+l_0` is (n-n_0-l_0+1)+(n_0+l_0)-1 = n > the same formula in row 2 (`n-n_r-l_r+3` from start n_r ends at n-l_r+2), so `n_0+l_0` overshoots by l_0-2 >= 1 and is off the line; `n_0+2` is the mirror image of row 2. The scholar's numerics agree. The misprint is real, in the current LaTeX (not checked on the PDF, as the scholar says).
- Numbering 4.3/4.5/4.6: the counter is shared (`\newtheorem{corollary}[theorem]`, section-numbered), and the order in section 4 is lemma, thm:AnToQuipu, corollary:2Rels, lemma:TreeOrientation, thm:QuipuToAn, cor:Equiv.. (lines 424-544). That gives 4.1 lemma, 4.2, 4.3, 4.4, 4.5, 4.6 provided the `lemma:firstRelToCord` is the first item of section 4. Plausible; no compile done. Fine as stated ("rests on the counter order").
- Thm 1.14.2 (claim 3): 2003.08914 l. 1554 states it as `Nikulin [Theorem 1.14.2]`, with (a) rank >= l_p+2 for p != 2 and (b) rank = l_2 => q = u_1+q' or v_1+q', even indefinite => one class in the genus. Confirmed. But note (a) and (b) together do NOT say "rank >= l_2+2"; the scholar says so.
- Gap: the title asserts "not Cor. 1.13.3". No source read quotes 1.13.3, so what is established is only "this statement is Thm 1.14.2 in the sources read". "Refuted if Nikulin's own text ..." is honest, but the title should not assert the negative.
- Gap: title claim 1, "exist in print only in Fosse", is a universal negative over the literature. Unread: Happel-Seidel, BGS, Ladkani math/0610685 Sec. 4 (scholar flags it as open), Nikulin. Read-only-as-radical-powers for HS/LMR/Ueda/Ladkani 0911.5137 is via summary or grep (Ueda grep "Kupisch": 0 hits is weak; Thm 1.1 is about L(s,t,u), not a Nakayama statement). Worse: Cor. 4.6 itself requires l_i >= 3 and almost separate; mixed pairs outside that class are not addressed.

## New?

- Misprint: grep "n_0+l_0", "misprint" in research/: only the unrelated 1310.1557 misprint. New.
- Cor. 4.6 / Cor. 2Rels / swap implemented: already in `research/literature/2305.06642-quipu-quivers-nakayama.md` (per the scholar) and E-274 (no mixed control). The 6 explicit pairs: new as a list.
- Thm 1.14.2: E-267 (EXPERIMENTS.md l. 62-66) already records "Thm 1.14.2 and Cor. 1.13.4 are what other authors cite" and that its statement was never quoted; E-267 calls 1.13.3 "Prop.", scholar says "Cor.". What is new is the quoted statement from 2003.08914 with (a),(b). The submission does not cite E-267 and is silent on the Prop./Cor. discrepancy and on Cor. 1.13.4 (which E-267 says is what 2202.08030 cites for uniqueness). It should say whether 1.13.4 is a third candidate.

## Evidenced?

Mostly. The 4888/1222/0 admissibility counts and the orbit-equals-`classPreservingOrbit` check are stated specifically (n = 9..12, 611 algebras, all l_i >= 3, >= 2 distinct digits). The table is clear. Weak points:
- "Cartan-level 48 of 48" is a sample (24 sources) and is a positive control only; the scholar says so. The Power line (n/a, YES controls cannot say NO) is honest and holds.
- Orbit = repo `classPreservingOrbit` is the project's own implementation of the same paper (author is the paper's author), so not independent; stated.
- "Pairs from Cor. 4.6 are derived equivalent" rests on the theorem's proof, not checked here; the claim is "printed", fine.

## Scope

Narrow to: "Among the sources actually read (listed), the only source containing mixed-digit LNA derived-equivalence pairs is 2305.06642 (Cor. 4.3, Thm 4.5, Cor. 4.6); HS, BGS, Ladkani math/0610685 Sec. 4 unread." Claim 3: "Thm 1.14.2 is the statement as restated in 2003.08914, 1704.00950, 1508.05251; Cor. 1.13.3's content is UNVERIFIED."

## Required for acceptance

1. Reword the title/claim 1 to the read-sources scope (drop "only ... in print"), and the title/claim 3 to drop "not Cor. 1.13.3" (state it as unverified).
2. Cite E-267 and reconcile Prop. 1.13.3 (E-267) vs Cor. 1.13.3 (E-265/this), and say whether Cor. 1.13.4 is relevant.
3. State that (a)+(b) is weaker than the single "rank >= l(A_q)+2" at p = 2 (already said in Consequence; move into the claim line).
4. Optional, small: one line verifying the end-vertex arithmetic for the misprint (start + length - 1 = n for `n_0+l_0`), so the misprint is checkable without running code. [`Skeptic` PDF check already queued.]
