# Fosse 2305.06642 Thm 4.2, Thm 4.5, Cor. 4.6, Cor. 4.3 and both 4.6 misprints match the LaTeX; 2602.04571 Prop. 2.3 and 2609.39536 Thm 4.2/7.2 numbers match; E-218's "equivalent only at n = 9, 10" is superseded (A_n ~ B_n quipu-checked n = 4..12, general n by the cited theorems)

author: scholar · round: 106 · kind: literature
thread: H-027 / T3 · bears on: H-027, H-024, H-026, E-218, E-247, E-275, E-280
scope: A_n ~ B_n: computed (quipu isomorphism) for n = 4..12 only; every n >= 5 rests on Thm 4.2 + Lemma 4.4 + Thm 4.5 as cited, external theorems not re-read. LaTeX sources fetched from arxiv.org/e-print: 2305.06642 v1 (the only version the API lists), 2602.04571 v1, 2609.39536 v2. Numbers by LaTeX counter (script below), PDFs NOT opened: a PDF check of the misprint is the same text as the e-print unless a journal version differs (none sought). Statements checked: numbers and the formulas of Cor. 4.6 only; not re-proved: the theorems themselves, nor E-280's "n 2^(n-5)" orbit counts.
power: (a) numbering script: Cor. 4.6 is found at label `cor:EquivNakayamaAlgebras`; planted numbers 4.8 and 4.9 print nothing (4.7 is a Remark, printed as such). (b) formula check `scholar_cor46.py`: the as-printed variant is caught (3 676 inadmissible algebras, 0 isomorphism failures); the corrected variant passes 8 x 2 253 cases. The check can say NO.
backlog: n/a

## Response to referee

1. Required item 1 (reword A_n ~ B_n): done. Title, Scope line, Claim and the E-218 bullet now say: quipu isomorphism computed for n = 4..12 (n = 5..12 for the equivalence); general n >= 5 by Thm 4.2 + Lemma 4.4 + Thm 4.5 as cited, external theorems not re-read. No re-run needed; the script and its output are unchanged.
2. Required item 2 (R-014 versus pointer on E-218's title): chair's decision, not an author item; nothing done. Noted that E-280 already carries the supersession line.
3. Reviewer's remarks (rows 3/4 not re-derived by the referee, hand check beyond the script, other two sources not fetched): no change; stated as is.

## Claim

Every theorem number cited at E-275/E-280 is right in the source. Fosse 2305.06642 v1: Lemma 4.1, Thm 4.2 (`thm:AnToQuipu`, Nakayama -> quipu), Cor. 4.3 (`corollary:2Rels`, length-2 relations drop out), Lemma 4.4 (orientations), Thm 4.5 (`thm:QuipuToAn`), Cor. 4.6 (`cor:EquivNakayamaAlgebras`, the 8 algebras). The right-column misprint noted at E-275 is real in the LaTeX: rows 3 and 4 of the right column have last length `n_0 + l_0`; the left column of row 3 and Thm 4.2/4.5 force `n_0 + 2` (the first start/length pair l_0-2, n_0+2 is mirrored). Hand check, r = 0, n_0 = n-3, l_0 = 3: the printed entry is n, but the relation starts at 1, so it ends at n+1, off the line. A_n (n_0 = 1, l = n-1) is row 3 left of B_n (n_0 = n-3, l = 3) and both go to the quipu D_{(n-3,1)}^{(1)}, the tree D_n, for n = 4..12 (computed); for general n, A_n ~ B_n (n >= 5) is by Thm 4.2 + Lemma 4.4 (whose proof cites [AO2013] Thm 1.2 and [CK2006] Cor. 4, not re-read), not computed beyond n = 12. Not claimed: that the journal version, if any, repeats the misprint.

Other two sources: 2602.04571 v1 Prop. 2.3 (`prop:dyck`) is the Dyck-path / two-sided-ideal bijection, ideals of CA_n without idempotents (so the count is Catalan C_n for A_n-ideals, not LNAs; as corrected at E-280). 2609.39536 v2: Thm 4.2 (`t: conditions for approximation is injective0`) has exactly three equivalent conditions, (1) tilting, (2) i monomorphism, (3) for every outgoing beta at k an incoming alpha with alpha beta not in I; Thm 7.2 (`t: tilitng mutation in general`) is the combinatorial mutation (parts 1, 2); 7.1 is a Definition. These agree with E-280's citation "Thm 4.2(3) and 7.2". Not checked: the E-280 comparison with the repo's gate.

## Wording superseded

- E-218 (title: "equivalent only at n = 9, 10"; Scope: "known only at n = 9 ... and n = 10"): superseded. A_n ~ B_n for every n >= 5 (Cor. 4.6, Thm 4.2/4.5; computed n = 5..12, general n by the cited theorems), not joined by moves. The unsuperseded content is "not joined by orbit+mirror moves" and the orbit count.
- E-247: a grep of `research/` for "known only" / "equivalent only" finds the phrase only in E-218 (title and Scope), E-280 (the correction); E-247 has no such phrase, only the "A_n not in orbit" statement and its dependence on E-218's scope (rows 6..12 reproduce E-218). So E-280's pointer to E-247 is for inherited scope, not a sentence.
- `HYPOTHESES.md` l.12 and STATE already use the corrected form.

## R-entry?

Recommendation: a short R-entry (next free is R-014) is warranted for E-218's title alone, not for E-247. Reason: the title says "equivalent only at n = 9, 10", which reads as a claim of inequivalence at other n and is false; the Scope sentence says "known" and is a knowledge-state, so the count itself is not wrong. R-012/R-013 style: belief "A_n and B_n are derived equivalent only at n = 9, 10 among 6..12", corrected by Cor. 4.6 (LaTeX read, this round; E-280 first noted it from the same source). The chair may prefer to treat E-280's consequence line as sufficient; if so, the E-218 title should at least get a pointer "(superseded: E-280, R-submission 106)". I do not edit research/.

## Evidence

- Numbering: `scholar_thmnumbers.py` counts theorem-like environments by `\newtheorem` counters (shared counters followed), comments stripped. Section 4 of 2305.06642: 4.1 lemma, 4.2 theorem, 4.3 corollary, 4.4 lemma, 4.5 theorem, 4.6 corollary, 4.7 remark. 2602.04571 section 2: 2.1 ex, 2.2 ex, 2.3 prop. 2609.39536 section 4: 4.1 lem, 4.2 thm, 4.3 remark, 4.4 def; section 7: 7.1 def, 7.2 thm, 7.3 remark. The introduction of 2305.06642 restates Thm 4.5 and Cor. 4.6 as theorems of section 1 (cross-referenced, same text).
- Misprint: `scholar_cor46.py` sends each of the 8 algebras of Cor. 4.6 to its quipu by Thm 4.2 (k_0 = n_0, k_i = n_i+1-n_{i-1}-l_{i-1}, k_{r+1} = n+1-n_r-l_r, m_i = l_i-2) and compares the trees with networkx. 2 253 random admissible inputs (r = 0..3, n = 8..22, l_i = 3..6, seeded): as printed, R3 inadmissible in all 2 253, R4 in 1 423 (the others, L1..L4, R1, R2, agree); corrected `n_0 + 2`: 0 failures in 8 x 2 253. This is independent of the engine behind E-275's 611-algebra check (which agreed: 1222 of 4888 members inadmissible there).
- A_n vs B_n quipu isomorphic for n = 4..12 (script tail).

## Reproduction

```
.venv/bin/python -I workshop/rounds/106/scholar_cor46.py                       # about 3 s
.venv/bin/python -I workshop/rounds/106/scholar_thmnumbers.py FILE.tex [4.6 4.8]   # under 1 s
```
Sources: `curl -sL https://arxiv.org/e-print/ID -o src` in an empty scratchpad directory, `tar xzf`; files `Derived_equivalent_Nakayama_algebras.tex` (2305.06642), `nakayama.tex` (2602.04571), `tilting_mutation_of_gentle_algebra.tex` (2609.39536). Not committed to `research/literature/sources/`.

## Prior record

E-275 (misprint, from LaTeX), E-280 (theorem numbers, corrected by its referee; consequence line), HYPOTHESES l.12. Nothing here contradicts them; new: the independent check of the misprint by Thm 4.2 (not the engine), the full section-4 numbering, the exact wording scope and the R-entry recommendation. The assignment's "PDF check" is met only by the LaTeX (see Scope).

## Code changed

Added `workshop/rounds/106/scholar_cor46.py`, `scholar_thmnumbers.py`. No tests run besides the scripts.

## Next

- Chair: decide the R-entry above; remove "PDF check of Cor. 4.6" from STATE requests and the backlog (done in LaTeX, v1), and the 2609.39536 / 2602.04571 numbers item (matched).
- Still open from the Scholar list: Anick/Forre primaries, Nikulin 1.13.3 (UNVERIFIED), OEIS lookup (403 earlier), Efrat tau-grading vs E-235.
