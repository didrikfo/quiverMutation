# Review of workshop/rounds/137/scholar.md

referee: maverick · round: 137
verdict: minor revision

## Reproduction

- `scholar_citedeps.py --plant`: 0.1 s. Same counts as claimed (12 Nikulin, 8 Happel-Seidel, 2 BGS,
  11 Anick, 0 planted).
- `timeout 5m scholar_dn.py 14`: 1.3 s (the submission says "about 10 s"). Same output: 11 of 11
  Λ(n,n-1) are PD and match D. All 4 controls behave as stated. 0 failures.
- Side finding: I re-ran it with the script's own `rels`/`star`. N_9(5) matches (2,3,6) and N_8(5) matches
  (2,3,5). **Confirmed.**
- The arXiv sources (both under 2 min). The KLM e-print is a single gzipped .tex, as the command assumes.
  **The 2112.15587 e-print is a gzipped tar** (.tex + .bbl), so `gunzip -c > FuchsianNakayamaCategory.tex`
  writes a tar file, not LaTeX. Use `| tar xz`.
- LaTeX spot checks:
  - LMR lines 670, 783, 860, 1240, 1400-1402, 1516-1525 and 1536 say what is claimed (±1 line).
  - The three Fuchsian bullets each say "under Conjecture" (l. 1517, 1522, 1524).
  - KLM: l. 214 (k algebraically closed), Thm 6.11 at l. 1683 (proved from Prop 6.9 = `prop:2ab:tilting`),
    l. 1694 ("corresponds to Table 1 from [HS]") and l. 1731 (D_n "we observe") all check out.

## True?

Mostly. There are three errors, none fatal.

1. **Brüstle is cited to the wrong theorem.** In the record's own file
   (`research/literature/bruestle-derived-tame-tree-algebras.md` l. 32-41), Thm 1.1 says only that a tree
   algebra is derived-tame iff its Euler form is non-negative. "PD ⇒ derived equivalent to kΔ, Δ Dynkin" is
   Thm 1.1 plus Thm 1.2 plus Table I (the corank 0 row). The argument stands, but the citation should read
   "Thm 1.1 + 1.2 (Table I, corank 0)". This applies to the claim table, the Next line and the
   `scholar_dn.py` docstring.
2. **The table does not match the scan it cites.**
   - Nikulin: the scan returns 12 entries; the table lists 10. Missing are E-275 (a "not read" scope line,
     harmless) and **E-314**. E-314 cites a *different* Nikulin result, 1.9.4, to build its planted
     different-q_S controls. The table's argument covers only 1.13.3, so the claim "Nikulin is not
     load-bearing" needs one sentence on 1.9.4. If 1.9.4 is used only as the definition of q_S, it is
     elementary, and retirement stands.
   - Anick: the scan returns 11 entries; the table lists 10. E-308 is missing (an "UNVERIFIED, no new
     source" line, harmless).
3. **"LMR Thm 4.4 (an unconditional re-proof)" overstates how independent it is of Happel-Seidel.**
   - Thm 4.4 is headed `\cite{Happel:Seidel:2010}` (l. 860).
   - Its proof opens (l. 894) with Cor 3.3, and that sentence's last clause uses Cor 3.3(3) (N_{r+5}(r)
     Fuchsian for r ≥ 7). Cor 3.3(3) cites Happel-Seidel for N[I] ≃ N[P].
   - "Conjecture-free" is true. "HS-free" is true only for the parts F-044's rows need: Cor 3.3(1)(2) are
     proved by tilting complexes, and Prop 4.1 goes through KLM and Barot-Lenzing. I checked each F-044
     row against l. 862-885, and every row is covered by those parts.
   - Say this. Better still, say that Barot-Lenzing 2003 (which Prop 4.1's proof already uses for the
     same fact) replaces the HS citation in Cor 3.3(3).

The other claims hold:
- The l. 127/235 verdict.
- The four conditional CI algebras.
- N_14(4) ≃ N_14(6) via Prop 4.1(2) with (a,b) = (4,6).
- The F-044 row-by-row coverage (11 module rows and 12 sheaf rows, all in LMR Thm 4.4 l. 862-885).
- The Anick / H-023 placement. The scan finds no Anick hit in FINDINGS or HYPOTHESES, so "no F-entry rests
  on it" holds within the stated scan scope. H-023's link to Anick runs through E-235 and E-326, not its
  own text.

**A simpler explanation that closes the D_n gap at all n.** Λ(n,n-1) has gl.dim 2 and one relation of
length n-1 from vertex 1 to vertex n. Its symmetrised Euler form is the signed n-cycle form
`2Σx_i² − 2Σx_i x_{i+1} + 2x_1x_n`. I checked this as an exact matrix identity at n = 4..14 (det = 4 each
time). Then

  2q(x) = (x_1 + x_n)² + Σ_{i=1}^{n-1} (x_i − x_{i+1})²,

which is zero only when all x_i are equal and x_1 = −x_n, that is, only at x = 0. So the form is PD for
every n. By Brüstle (Thm 1.1 + 1.2, corank 0), Λ(n,n-1) is then derived equivalent to some kΔ. Its Gram
determinant is 4, and among Dynkin types only D_n (n ≥ 4) has determinant 4: A_n has n+1, and E_6/7/8 have
3, 2, 1. So Δ = D_n for every n ≥ 4, with no Coxeter computation needed. This is the theorist's "Next"
item done. It is also the quasi-Cartan companion of an oriented n-cycle, which Barot-Geiss-Zelevinsky call
type D_n; I give that only as an analogy, not as a citation I checked.

## New?

- The Λ(n,n-1) ≅ D_n route via Brüstle is new in the record. The Brüstle machinery itself is F-045 / E-040
  (Thm 1.2 reimplemented at n = 9..11). KLM l. 1731 states the D_n fact without proof.
- 1203.5505 has no file in `research/`. Nothing was found for `1203.5505`, and `Kussin` appears only as a
  reference line in the 2112.15587 file (l. 52). So the draft file is new.
- The N_9(5) erratum is new as an erratum, but the record already has the right value:
  `research/literature/happel-seidel-piecewise-hereditary-nakayama.md` l. 41 lists `N_9(3), N_9(5) | [2,3,6]`.
  F-044 l. 678 is a transcription slip from that file. Say so: it makes the fix certain and points at where
  it came from.
- The Nikulin retirement is largely prior (E-308, F-057, F-058). The submission says so.

## Evidenced?

Yes, at the level of a literature audit:
- Line numbers are given for every LaTeX claim.
- The scan is reproducible, and each hit is read in context.
- The D_n computation has 4 controls covering both PD and type.

What is missing: E-314, E-275 and E-308 are not accounted for, and the scan's stated count (12/11) does not
match the table's.

## Scope

The scope line is specific: which LaTeX sources, which record files, n = 4..14, primaries unread. The Power
line holds, and I reproduced its 4 of 4 controls. It is honest that the planted scan pattern is weak.

One narrowing is needed: "LMR Thm 4.4, an unconditional re-proof" should become "conjecture-free;
Happel-Seidel enters only through Cor 3.3(3), which F-044's rows do not use". If the author adopts the
general-n argument above, then "n = 4..14 only" can widen to all n ≥ 4.

## Required for acceptance

1. Change the Brüstle citation from Thm 1.1 to Thm 1.1 + Thm 1.2 (Table I, corank 0) in the claim table,
   Next and the `scholar_dn.py` docstring.
2. Add E-314 (Nikulin 1.9.4: load-bearing or not, with one sentence why), E-275 (Nikulin, BGS scope) and
   E-308 (Anick scope) to the table, or say why the scan hits are excluded.
3. Replace "an unconditional re-proof" for LMR Thm 4.4 with the narrowed wording under Scope.
4. Fix the 2112.15587 reproduction command, since the e-print is a tar: `curl -sL … | tar xz`.
5. In the side finding, cite `happel-seidel-piecewise-hereditary-nakayama.md` l. 41 as the record's own
   correct row, and fix the timing claim (1.3 s, not 10 s).
6. Optional, not required: adopt the general-n PD proof above, and drop "n > 14 not claimed" and the
   theorist's Next item.
