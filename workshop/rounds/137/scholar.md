# The 2112.15587 summary (l. 127, 235) matches the v2 LaTeX; of the four UNVERIFIED primaries only Anick 1982 is load-bearing (for H-023, OPEN): Nikulin, Happel-Seidel and BGS can be retired, Happel-Seidel with a re-derivation note (Λ(n,n-1) ≃ kD_n now proved for all n ≥ 4)

author: scholar · round: 137 · kind: backlog
thread: T3/T8, T5, H-027 literature · bears on: F-044, F-046, F-055, F-056, E-328, E-265, E-267, E-235, E-240, E-326, H-023
scope: LMR 2112.15587v2 LaTeX (Thm 6.3, its proof, Prop 3.5, 4.1, Thm 4.4, Conj 5.7, 6.2); KLM 1203.5505v1 LaTeX (Thm 6.11, Table 1); record scan by keyword over FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS only (not literature/ or workshop/); D_n re-derivation proved for all n ≥ 4 (hand proof above), computer-checked at n = 4..14. Primaries themselves still not read.
power: (1),(2) literature: n/a. The D_n computation says NO where NO is true: N_9(3) (affine, semidefinite) and N_12(3) (not PH) are refused as not positive definite, and N_7(3), N_8(3) are PD but matched to E, not D (4 of 4 controls); the script also asserts the signed-cycle identity and det 4 for each Λ(n,n-1), n = 4..14. The keyword scan's planted pattern returns 0 entries (weak: it shows only that the scan does not invent hits).
backlog: (1) LMR summary correction (STATE l. 35, 75) settled; (2) four UNVERIFIED primaries (cycle-1 backlog 18, STATE l. 38) Nikulin retired, BGS retired, Happel-Seidel retired with re-derivation note, Anick narrowed only (carry); (3) owed literature files: KLM 1203.5505 drafted below, others not done.

## Response to referee (maverick, minor revision)

1. **Brüstle citation: done.** Now "Thm 1.1 + Thm 1.2 (Table I, corank 0)" in the claim table, Next and the
   `scholar_dn.py` docstring (which also cites the record file l. 32-41).
2. **Scan hits: done.** E-314, E-275 (Nikulin) and E-308 (Anick) added to the table. E-314 cites Nikulin 1.9.4
   only as the name of the discriminant form q_S, to build planted pairs with a different q_S; its NO uses only
   that an integral isometry P^T S P = S' carries S^{-1}Z^n/Z^n and its form onto those of S' (elementary, the
   same fact as in my row). Not load-bearing; retirement stands. Whether "1.9.4" is the right number is not
   checked (primary unread). E-275 and E-308 are "not read" / "UNVERIFIED, no new source" scope lines. Counts
   now 12/12 Nikulin, 11/11 Anick. Nikulin's E-308 was already in the table.
3. **LMR Thm 4.4 wording: done.** Now "conjecture-free; Happel-Seidel enters only through Cor 3.3(3), which
   F-044's rows do not use", with Barot-Lenzing 2003 named as the replacement for that HS citation. The
   same narrowing applies to "unconditional" in my parenthetical in the 2112.15587 file (it means
   conjecture-free there, which is true); I leave that file alone this sitting, as instructed. Chair: read it so.
4. **arXiv command: done.** 2112.15587 is now `| tar xz` (the e-print is a tar); 1203.5505 stays `gunzip -c`.
5. **Side finding: done.** The record's own file `happel-seidel-piecewise-hereditary-nakayama.md` l. 41 has
   the correct row (`N_9(3), N_9(5) | [2,3,6]`); F-044 l. 678 is a transcription slip from it. Timing fixed
   (1.4 s on my re-run; referee 1.3 s).
6. **Optional general-n proof: adopted.** I re-derived the Euler form from projective resolutions
   (pd S_1 = 2 via 0 → P_n → P_2 → P_1; pd S_i = 1 otherwise), giving ⟨x,y⟩ = Σx_iy_i − Σx_iy_{i+1} + x_1y_n,
   so C^{-1}+C^{-T} is the referee's signed n-cycle and 2q = (x_1+x_n)² + Σ(x_i−x_{i+1})² is PD for all n.
   det = Π_{ω^n=−1}(2−ω−ω^{-1}) = 4 for all n, so by Brüstle 1.1 + 1.2 the type is D_n (only Dynkin type with
   det 4 at n ≥ 4). `scholar_dn.py` now also asserts the signed-cycle identity and det 4 at n = 4..14.
   Claim, Scope and Next changed accordingly; "n > 14 not claimed" and the theorist's Next item dropped.

## Claim

**(1)** The summary `research/literature/2112.15587-nakayama-fuchsian-singularities.md` l. 127 and l. 235 is
**correct** against the v2 LaTeX (`FuchsianNakayamaCategory.tex`, 2022-02-07): the lead-in at l. 1400 puts
Thm 6.3 (l. 1402) under Conj 6.2 and Conj 5.7; the proof of part (3) for `⟨2,3,10]` (l. 1516) uses Prop 3.5
(`N_{r+6}(r)`, r = 9, l. 670) and Prop 4.1(2) with (a,b) = (3,9), 4.1(3) with (3,8) (l. 783), and no
conjecture; the `⟨2,5,7]`, `⟨2,6,6]`, `⟨2,3,11]` bullets each say "under Conjecture 5.7" (Conj 5.7 is
exactly those three identifications, l. 1240); completeness uses (F2) (l. 1536). Theorem numbers 3.5, 4.1,
5.5, 5.7, 6.2, 6.3 re-matched by counter order. One adjacent imprecision, l. 129-131, I corrected (edit below):
the "non-PH certificate by Prop 4.3" is conditional for N_14(4), N_14(6), N_14(5), N_16(8), but those are
non-PH anyway by Thm 4.4 (conjecture-free: the proof at l. 893-950 cites no conjecture; Happel-Seidel enters only through Cor 3.3(3), which F-044's rows do not use), and N_14(4) ≃ N_14(6)
is unconditional (Prop 4.1(2), (a,b) = (4,6)). Only the three Fuchsian *types* are conditional.

**(2)** Which entries rest on each primary, and my recommendation:

| primary | entries that cite it | load-bearing? | recommendation |
|---|---|---|---|
| Nikulin 1979 (1.13.3 "class = genus"; 1.9.4) | E-197, E-265, E-267, E-272, E-308, E-320, E-334 (cite/scope); E-314 (1.9.4, as the name of q_S for planted pairs); E-275 ("not read" scope line); F-055, F-057, F-058 ("does not enter" or "UNVERIFIED") | **No.** Every use is the YES direction (same genus => Z-congruent), stated as "would follow". Every NO uses only that an isometry P^T S P = S' preserves the discriminant form q_S (immediate) or the Ladkani/Rickard/AI chain (E-308, E-272). E-314's planted NO uses only that same invariance. F-057 and F-058 already say Nikulin does not enter. | **Retire as "not load-bearing".** F-055's "the lattice step behind q_S (Nikulin Cor 1.13.3)" should read "q_S is an isometry invariant (elementary); Nikulin 1.13.3 would only give the converse, not used". |
| Happel-Seidel 2010 | F-044 (table), F-046 (the symmetry), F-056 + E-328 (Λ(n,n-1) of type D_n), E-040, E-231, E-274, E-275 | **Was**, for F-056 and F-044; **now replaceable** by sources read in LaTeX or PDF (next rows) | **Retire as "load-bearing but unverified", with re-derivation note:** |
| - the classification table (F-044) | | | LMR Thm 4.4 (v2 l. 860; conjecture-free; Happel-Seidel enters only through Cor 3.3(3), which F-044's rows do not use, and Barot-Lenzing 2003 can replace that citation): every one of the record's 11 module and 12 sheaf rows is in it (N_9(5) under [2,3,6], see side finding). |
| - the symmetry D^b(N_n(a)) ≃ D^b(N_n(b)), n = (a-1)(b-1) (F-046) | | | KLM 1203.5505v1 Thm 6.11 (l. 1683, via Prop 6.9), and LMR Prop 4.1(1) which cites it. |
| - Λ(n,n-1) ≃ kD_n (F-056, E-328) | | | Brüstle 2001 Thm 1.1 + Thm 1.2 (Table I, corank 0) (read from PDF; PD Euler form ⇒ derived-tame ⇒ derived equivalent to kΔ, Δ Dynkin, fixed by n and Dynkin type) + all-n proof: 2q = (x_1+x_n)² + Σ(x_i−x_{i+1})² is PD, det 4 ⇒ D_n (Response item 6). Checked at n = 4..14, 0 failures (`scholar_dn.py`). Covers F-056's n = 5, 6, 7 and all n ≥ 4. |
| BGS (Bobiński-Geiss-Skowroński) | E-274, E-275, only in "not read" scope lines | **No.** Nothing cites it for a statement. | **Retire** (drop from the PDF wish list). |
| Anick 1982 (Thm 3.1/3.2, strong freeness) | E-235 (Lemma A), E-237, E-240, E-273, E-278, E-290, E-322, E-326 ((H3)(ii) "two Anick facts marked weakest"), E-227, E-264; E-308 ("UNVERIFIED, no new source" scope line) | **Yes, for H-023 (OPEN) only**; no F-entry rests on it. The monomial case of Lemma A is re-derived in-house (E-240, refereed line by line); open joints: E-235 step 4 landing in E-240's hypotheses, loops of degree >= 1 at i, E-326 (H3)(ii). | **Carry**, re-labelled "Anick (monomial case) re-derived in E-240; general case UNVERIFIED, load-bearing for H-023 only". Retire if H-023 is retired. |

Not claimed: that the primaries say what is attributed to them (still unread); that the keyword scan finds
an entry that relies on a primary without naming it.

## Evidence

- LMR v2: lines quoted above read directly; `rounds/110/scholar_thmnumbers.py` gives 3.3 = Cor (familis),
  3.5 = Prop N-r+6-r, 4.1 = HS-symmetry, 5.5 = tiltingobjects, 5.7, 6.2, 6.3. Prop 4.1(1) is proved by
  citing KLM Thm 6.11; (2),(3) by perpendicular categories and Barot-Lenzing 2003. Cor 3.3(3) still cites
  Happel-Seidel for N[I] ~ N[P] (one-point extension), the same fact Prop 4.1's proof attributes to
  Barot-Lenzing; Thm 6.3(3) at (2,3,10) does not go through Cor 3.3.
- KLM v1: Thm 6.11 "Let a, b >= 2. Then for n = (a-1)(b-1) the algebras A_n(a), A_n(b) and the stable
  category vect-X of type (2,a,b) are derived-equivalent"; k algebraically closed (l. 214). Its Table 1
  "corresponds to Table 1 from [HS]" and states A_n(n-1) is of type D_n as an observation, not a proof.
- Record scan: `scholar_citedeps.py` (12 entries Nikulin, 8 Happel-Seidel, 2 BGS, 11 Anick, 0 planted); each
  hit's line read in context (E-197 l. 1143, E-265 l. 635-637, E-267 l. 615-619, E-308 l. 312, E-334 l. 118,
  F-055 l. 54, F-056, E-328, E-240 l. 839, E-326 l. 151, E-273 l. 572-574).
- D_n: Λ(n, n-1), n = 4..14: 11 of 11 signed-cycle form with det 4, PD, Coxeter polynomial = treeCoxeterKey(star(2,2,n-2)), not A_n's. General n: Response item 6.
- Side finding (record erratum, not fixed): F-044's module table row "`N_6(3)`, `N_9(5)` | `[2,3,r]`" cannot be
  right for N_9(5): [2,3,5] has 8 vertices. Computed: N_9(5) matches [2,3,6], not [2,3,5]; N_8(5) matches
  [2,3,5]. LMR Thm 4.4 lists N_9(5) under [2,3,6], and so does the record's own `happel-seidel-piecewise-hereditary-nakayama.md`
  l. 41 (`N_9(3), N_9(5) | [2,3,6]`): F-044 l. 678 is a transcription slip from it. The row should read `N_6(3)` alone (or `N_9(5)` moved to
  the `[2,3,6]` line); the 11/11 count is unaffected if the code compared [2,3,6].

## Reproduction

```
.venv/bin/python -I workshop/rounds/137/scholar_citedeps.py --plant      # < 1 s
timeout 5m .venv/bin/python -I workshop/rounds/137/scholar_dn.py 14      # about 1.4 s
curl -sL https://arxiv.org/e-print/2112.15587 | tar xz      # in a scratch dir; tar of .tex + .bbl
curl -sL https://arxiv.org/e-print/1203.5505  | gunzip -c > klm.tex
```

## Prior record

E-308 already removed Nikulin from the necessity chain; F-057/F-058 already say it does not enter. My R134
notebook proposed the l. 127/235 correction, applied since; this round only re-checks it. E-240 is the
in-house Anick re-derivation. The cycle-2 plan (plan/scholar.md l. 22) proposed retiring all four as
"not load-bearing"; that is right for Nikulin and BGS, wrong for Anick (H-023 needs it) and, for
Happel-Seidel, true only after the re-derivation above.

## Code changed

New scripts `workshop/rounds/137/scholar_citedeps.py`, `workshop/rounds/137/scholar_dn.py`. No library
change, no tests needed. **Literature edit (permitted under (1)):** `research/literature/2112.15587-nakayama-fuchsian-singularities.md`,
after l. 131, a five-line parenthetical (conditional Prop 4.3 route for the four CI algebras; Thm 4.4 and
Prop 4.1(2) conjecture-free). L. 127 and l. 235 unchanged.

## Draft literature file (for the chair): 1203.5505 Kussin-Lenzing-Meltzer

*Triangle singularities, ADE-chains, and weighted projective lines*, Adv. Math. 237 (2013); arXiv v1 LaTeX read
R137. k algebraically closed. For X of weight type (a,b,c): vect-X is Frobenius with the line bundles as
projective-injectives, the stable category is triangulated, fractionally CY (Sec. 7, Prop 7.1, 7.3, 7.5), has a
tilting "cuboid" (Sec. 6). **Prop 6.9 / Thm 6.11 (Happel-Seidel symmetry):** for n = (a-1)(b-1), A_n(a), A_n(b)
and vect-X(2,a,b) are derived equivalent (A_n(r) = our N_n(r)). Table 1 = HS Table 1 extended to n <= 12; it
says A_n(n-1) is of type D_n (observation). Spontaneous series A_{(r-1)(rn-2)}(rn), A_{(r-1)rn}(rn) are
conjectural (checked to n <= 50). Use: the sourced form of the HS symmetry (F-046); its fractional CY
statements are the background of F-058. Caveat: radical powers only.

## Next

- Chair: in F-055 replace the Nikulin clause as in the table; F-056 / E-328 "Happel-Seidel, UNVERIFIED" ->
  "Brüstle Thm 1.1 + 1.2 (Table I, corank 0) + signed-cycle PD proof, all n ≥ 4 (R137)"; fix F-044's N_9(5) row
  (slip from the HS file l. 41); drop BGS and Nikulin from the PDF wish.
- Scholar: the remaining owed files (LdlP math/0611532, Lenzing 1909.10362, 1310.1907, math/0101145, 0906.3422).
