# Ready-to-apply citation patches (E-171), round 058, scholar

Each patch: file, exact OLD text (one or more whole lines, copied from the file as of HEAD 1eb46ea), exact NEW text. Source checked: `research/literature/sources/1009.3370v3-aihara-iyama-silting-mutation.tex` (AI) and `2509.12983v2-pavon-chz-criterion.tex` (CHZ). Nothing was fetched from arXiv: the LaTeX in the repo was enough. Numbers are LaTeX's (one counter per section). Patches 1-4 are in `research/literature/`; patch 5 (provenance only) and 6 (syntheses, a dated checkpoint, annotated rather than rewritten) are optional. Patches 2 and 4 carry the "no monomial" correction.

Checks behind the wording (all by grep in the source): CHZ Cor 3.6 (label `cor:path-algebra`, line 1410) reads "Let Lambda = kQ/I be a path algebra with relations" -- the words "monomial" and "kA_n" do not occur in the file; "admissible ideal of relations" is in Example 3.4 (line 1299, label `ex:path-algebra`); Example 3.4 asserts that elements of soc P_i "are represented by tail-maximal paths" (line 1313-1315). AI Thm 2.32 (label `when silting is tilting`, line 982): M tilting; (a) D covariantly finite, mu^+ tilting iff each M has a left D-approximation f with Hom(D,f) injective; (b) dual. AI Prop 2.3 (line 375-380) is "Let T be an algebraic triangulated category. If T has a tilting object M, then T is triangle equivalent to K^b(proj End(M))", cited to Keller.

## Patch 1 -- `research/literature/README.md`, line 26 (the 1009.3370 row)

OLD:
```
| [1009.3370](1009.3370-silting-mutation.md) — Aihara–Iyama, Silting mutation in triangulated categories | Theorem 2.32: an **iff** for "this mutation is a derived equivalence", where our gate is necessary only; and transitivity for piecewise hereditary algebras |
```
NEW:
```
| [1009.3370](1009.3370-silting-mutation.md) — Aihara–Iyama, Silting mutation in triangulated categories | Theorem 2.32(b): an **iff** for "the mutated object is *tilting*" (injectivity of `Hom(g, D)`; tilting then gives a derived equivalence `End(T) ~ A` by Prop. 2.3 / Keller, the converse is not claimed), where our gate is necessary only; and transitivity for piecewise hereditary algebras (Thm 3.1, hereditary or canonical algebras only) |
```
Reason: 2.32 is an iff for tilting, not for "derived equivalence" (audit row 10); Thm 3.1 scope is row 18 (the summary already says so).

## Patch 2 -- `research/literature/README.md`, line 37 (the 2509.12983 row)

OLD:
```
| [2509.12983](2509.12983-chz-criterion-derived-equivalences.md) — Pavon, Detecting derived equivalences with the CHZ criterion | Cor. 3.6: an iff for an HRS-tilt to be a derived equivalence, stated for `kA_n/I`; a set-mutation our engine lacks |
```
NEW:
```
| [2509.12983](2509.12983-chz-criterion-derived-equivalences.md) — Pavon, Detecting derived equivalences with the CHZ criterion | Prop. 3.5 (any artin algebra) and Cor. 3.6 (stated for `kQ/I`, no "monomial" hypothesis printed; its Example 3.4 silently needs one): an iff for an HRS-tilt to be a derived equivalence. Cor. 3.6 is an iff for monomial `I`, hence for every `kA_n/I` (our derivation, E-171); a set-mutation our engine lacks |
```

## Patch 3 -- `research/literature/1504.02617-quivers-for-silting-mutation.md`, lines 132-133

OLD:
```
  surviving negative-degree arrow is precisely the failure of the tilting
  condition `Hom(T, T[<0]) = 0` of Aihara–Iyama Theorem 2.32.
```
NEW:
```
  surviving negative-degree arrow is precisely the failure of the tilting
  condition `Hom(T, T[<0]) = 0` (Aihara–Iyama Definition 2.1(b); their
  Theorem 2.32(b) says that, for a mutation of a tilting object, this is
  equivalent to injectivity of `Hom(g, D)`).
```
Reason: audit row 9. (The "precisely" identification of arrows with Hom(T,T[<0]) is Oppermann's dg-quiver statement and was not re-checked here: UNVERIFIED, no Oppermann source in `sources/`.)

## Patch 4 -- `research/literature/2509.12983-chz-criterion-derived-equivalences.md`, line 70

The r055 correction sits at line 131, far from the statement; a reader of the quoted Corollary sees no caveat. Line 70 is exactly `> **Corollary 3.6** (` + `Λ = kQ/I` + `). ` (trailing space, then end of line).

OLD:
```
> **Corollary 3.6** (`Λ = kQ/I`). 
```
NEW:
```
> **Corollary 3.6** (`Λ = kQ/I`, as printed: no "monomial" hypothesis appears in the paper; its Example 3.4 needs one -- see the r055 correction below, and E-171). 
```
(The r055 correction is at line 131 of the same file.)

## Patch 5 (optional) -- `research/literature/rickard-morita-theory-derived-categories.md`, lines 8-10

I found no wrong statement at line 9 or 35: Prop 2.3 is already described there as "an algebraic triangulated category with a tilting object". Only provenance can be added.

OLD:
```
Aihara–Iyama, arXiv:1009.3370 Definition 2.1(b), Example 2.2(a) and Proposition
2.3, which *was* read; (b) an elementary computation done here from those
```
NEW:
```
Aihara–Iyama, arXiv:1009.3370 Definition 2.1(b), Example 2.2(a) and Proposition
2.3, which *was* read (numbers checked against the LaTeX in round 055, E-171); (b) an elementary computation done here from those
```

## Patch 6 (optional) -- `research/syntheses/001-rounds-001-052.md`

The synthesis is a dated checkpoint; use the project's annotation style, not a rewrite.

6a, lines 36-38. OLD:
```
  tilting exactly when some `J_i = Hom(S_v, e_i A) != 0` (E-124, E-128,
  E-130; conditional on Aihara-Iyama Thm 2.31, which could not be checked
  against the paper). `dim J_i <= d_i - 1` with `d_i = dim e_i A e_v` (E-128).
```
NEW:
```
  tilting exactly when some `J_i = Hom(S_v, e_i A) != 0` (E-124, E-128,
  E-130; conditional on Aihara-Iyama Thm 2.31 *[round 055, E-171: Thm 2.31 and
  Thm 2.32(b) checked against the LaTeX; the iff is 2.31 + 2.32(b)]*).
  `dim J_i <= d_i - 1` with `d_i = dim e_i A e_v` (E-128).
```
6b, lines 84-85. OLD:
```
derived equivalence rests on Aihara-Iyama 2.31 / Ladkani 2.3(c), cited but
unread, and the independent test is Cartan-level with generation assumed; the
```
NEW:
```
derived equivalence rests on Aihara-Iyama 2.31 / 2.32(b) (read, E-171) and Ladkani 2.3(c), and the independent test is Cartan-level with generation assumed (generation proved per loopless step, E-168); the
```
6c, lines 194-196. OLD:
```
3. Two citations carry weight and have not been read: Aihara-Iyama
   arXiv:1009.3370 (Thm 2.31/2.32) and CHZ arXiv:2509.12983 (Cor 3.6).
   `arxiv.org` is blocked by the environment's network policy.
```
NEW:
```
3. Two citations carry weight: Aihara-Iyama arXiv:1009.3370 (Thm 2.31/2.32)
   and CHZ arXiv:2509.12983 (Cor 3.6). *[Round 055, E-171: both read from
   the LaTeX; the numbers stand; Cor 3.6 states no "monomial" hypothesis and
   is an iff only for monomial `I`.]*
```
(Line 190 "conditional only on AI 2.31" may stay: 2.31 is a theorem and the wording remains true.)

## Not patched (audit rows already done or outside `literature/`)
Rows 16-19 (provenance of the two summaries): already edited in `1009.3370-silting-mutation.md` and `2509.12983-...md`, checked present at lines 8-16 / 131. EXPERIMENTS annotations (lines 245, 303, ...) were applied by the chair in round 055. `tests/test_gate_without_tilting.py:5` (row 2) is outside `research/`.
