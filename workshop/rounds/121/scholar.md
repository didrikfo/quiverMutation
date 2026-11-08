# No published splitting lemma covers H-031: Keller-Yang Lemma 2.10 and DWZ Thm 4.6 remove only a *trivial part*, in complete Ginzburg / QP settings, under a direct-sum hypothesis the raw kM of M1, M4, M5 fails; M5's reduced part still needs Lemma A

author: scholar · round: 121 · kind: literature
thread: T5 · bears on: H-031, H-023, E-290, E-294, E-309, E-273, E-240
scope: three sources read in LaTeX (Keller-Yang 0906.0761v3 Sec. 2; DWZ 0704.0649v4 Sec. 4; Oppermann 1504.02617v1 Sec. 2, 3, 5, 7); hypotheses tested mechanically on the formal models M1, M2, M4, M5 and the no-loop control M6 of `workshop/rounds/106/theorist_loopfilt.py` only; journal versions unseen; other candidate sources (below) not read
power: the direct-sum test (d) says NO on the raw kM of M1, M4, M5 and YES on the E-309 substituted M5 generators (a known direct sum); a planted mixing term (d ta = pa + w fb) is caught (NO). Theorem numbers by counter script; a planted label (DWZ "9.9") prints nothing.

## Claim

None of the three sources contains a lemma whose hypotheses H-031 meets. (1) **DWZ Thm 4.6** (Splitting Theorem) is about quivers with potential: a potential in the *complete* path algebra of a degree-0 quiver *without loops* (Def. 4.1, eq. `no-loops`); there is no differential and no degree-1 letter, and kM has a degree-0 loop in every model (the letter ab at j). (2) **Keller-Yang Lemma 2.10** is the dg form: for the *complete Ginzburg* dg algebra of a *direct sum* (Q',W') + (Q'',W'') with (Q'',W'') trivial, the projection Gamma -> Gamma' is a quasi-isomorphism. It has three hypotheses kM fails: kM is not a Ginzburg algebra (no degree-2 loop at each vertex, 5 of 5 models); it is not completed; and on the raw generators the trivial letters are not a direct summand (4, 12, 10 mixing terms in M1, M4, M5). (3) **Oppermann** works in our setting (T_R A over R = k^n, not completed, loops of degree >= 1 allowed), but states no general cancellation lemma: Thm 7.5 is specific to the pairs (alpha_*, alpha-bar), Prop. 7.6 has a two-sentence proof, Sec. 5 cancels "up to quasi-isomorphism" with no citation (E-273 already says so).
What the literature *does* give is the **mechanism**, matched in the source: DWZ Lemma 4.7 (a unitriangular automorphism kills the mixed terms u_k, v_k; existence by a limit, Lemma 4.8, so it needs completion) plays the role of H-031's "triangular substitution exists"; KY Lemma 2.10's proof (filter the kernel by word length; the associated graded is (U + V)^{tensor i} / U^{tensor i} with V contractible) plays the role of the Kunneth step of E-290/E-309. Even granted, it cancels only the trivial part. In M5 three Theta letters (f/m, af/m, mf/m in raw names; tmp, samp, smmp after E-309) lie in the *reduced* part (d has no linear term), so the remaining step is Lemma A (E-240), which no trivial-part lemma covers. **Not claimed:** that no lemma exists anywhere; that H-031 is false; anything about Hom/End(T). Refuted by: a matched statement (any paper) for non-completed semifree dg quiver algebras over k^n that cancels a pair (theta, linear term of d theta) with mixed higher terms.

## Evidence

Hypothesis table: `scholar_splithyp.py`, raw kM from the round-106 engine (Def. 3.2 signs `last = -1`). "Linear" = letters whose d has a one-letter term.

| model | loops of kM (deg) | Ginzburg shape | linear pairs (theta -> target) | direct sum on raw gens (KY 2.10) | Theta in reduced part |
|---|---|---|---|---|---|
| M1 (d f = ab) | ab (0); f, af/a (1) | no | f/a -> b; af/a -> ab | NO (4 mixing terms) | none |
| M2 (d f = 0) | ab (0); f, af/a (1) | no | none | no trivial part | f/a, af/a |
| M4 (two a's) | a1b, a2b (0); 5 deg-1 | no | 6 thetas -> 3 targets (not bijective: needs a linear change, DWZ Prop. 4.4 type) | NO (12) | none |
| M5 (d f = ab + meb) | ab (0); f, af/a, mf/m (1) | no | f/a -> b; af/a -> ab; mf/a -> mb | NO (10) | f/m, af/m, mf/m |
| M6 no-loop control | ab (0); au (1) | no | none | no trivial part | m/a |
| M5 after E-309 substitution (control) | | | ta->pa, saa->qaa, sma->qma | YES (0) | tmp, samp, smmp |
| planted: d ta = pa + w fb | | | | NO (1) | |

Reading: M1 is *entirely* trivial part (that is why E-290 needed only change of generators + Kunneth: KY's mechanism, non-completed); M2 and M6 are entirely reduced (Lemma A / g-filtration, E-294, E-240); M4 and M5 are mixed. So H-031 = [DWZ-4.7-type finite substitution] + [KY-2.10-type trivial-part removal, non-completed] + [Lemma A on the reduced part]. The first two are known in a neighbouring setting; their transfer to non-completed kM is ours. Least sure: whether KY's word-length argument works without completion here. In the uncompleted algebra K/K_n contractible for all n does not alone give H(K) = 0 (d raises length); it does if each weight cell is finite-dimensional with letters of positive weight. E-284's weight is negative on star letters, so this is not immediate. Not checked.

Matched in the LaTeX (counter script `workshop/rounds/110/scholar_thmnumbers.py`, cross-check: KY cites DWZ "Proposition 4.4, 4.5, Theorem 4.6", same as my count):
- DWZ 0704.0649v4: Def. 2.2 (complete path algebra), Def. 4.1 (QP; no loops), Def. 4.3, Prop. 4.4 (trivial QP = sum of a_k b_k after linear change), Prop. 4.5, **Thm 4.6**, **Lemma 4.7**, Lemma 4.8.
- Keller-Yang 0906.0761v3 (counter shared with subsections): Thm 2.3 (restates DWZ 4.4-4.6; "all quivers are quivers without loops"), Lemma 2.9 (right-equivalent QPs give isomorphic Ginzburg algebras), **Lemma 2.10** (`L:reduction`), Lemma 7.4 (Mittag-Leffler).
- Oppermann 1504.02617v1: Sec. 2 setting, Def. 3.2 Steps 1-4, Thm 7.5, Prop. 7.6 (as E-273).

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/121/scholar_splithyp.py      # 0.04 s
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py <0704.0649 main.tex> th:trivial-reduced-splitting lem:killing-u-v 9.9
```
Fetched (STEERING standing arXiv permission), into the scratchpad, not committed: `curl -sL https://arxiv.org/e-print/0906.0761` (v3, KellerYangDerEqQuivMut.tex) and `.../0704.0649` (v4, derweyzel-selecta-rev.tex); versions from export.arxiv.org/api. Why: the cycle-2 plan names them as the candidate published form of H-031.

## Prior record

grep of `research/` for 0906.0761, 0704.0649, Keller-Yang, Derksen, splitting, trivial part, minimal model: no hit outside the Oppermann source file (its l. 274, 526 cite KY and DWZ, no number). E-273 recorded Oppermann's side (no Lemma A; loop case needs a unitriangular cancellation); this submission adds the KY/DWZ side and the per-model trivial/reduced split. Not a rediscovery.

## Proposed literature summaries (for the chair to promote)

`research/literature/0906.0761-keller-yang-derived-equivalences-qp.md`: **Keller-Yang, Derived equivalences from mutations of quivers with potential** (arXiv v3; Adv. Math. 226 (2011), publication UNVERIFIED from memory). Take: Lemma 2.10: for the complete Ginzburg dg algebra of (Q',W') + trivial (Q'',W''), projection to Gamma(Q',W') is a q.i.; proof by word-length filtration of the kernel, associated graded (U+V)^i/U^i with V contractible, then Mittag-Leffler (Lemma 7.4). Thm 2.3 restates DWZ 4.4-4.6 for loop-free quivers. Caveats: Ginzburg algebras only; completed; direct sum assumed (no mixed terms). Does not cover H-031 (E-3xx, this round).

`research/literature/0704.0649-dwz-quivers-with-potentials.md`: **Derksen-Weyman-Zelevinsky, Quivers with potentials and their representations I** (arXiv v4; Selecta Math. 14 (2008) 59-119, as printed in Oppermann's bibliography). Take: Splitting Theorem 4.6 (QP = trivial + reduced up to right-equivalence) and Lemma 4.7 (unitriangular automorphism kills mixed terms, by the limit in Lemma 4.8). Caveats: degree 0 only, complete path algebra, no loops (Def. 4.1). It is the model for H-031's substitution step, not a proof of it.

## Code changed

New script `workshop/rounds/121/scholar_splithyp.py` (imports the round-106 and round-118 scripts, read-only). No library change, no tests run.

## Next

- **H-023:** the loop joint does not close by citation. It stays the workshop's own (H-031), and is still moot for library steps (E-279).
- **Theorist (literature+proof slot, round 122):** write H-031 as three named steps: finite triangular substitution (cite DWZ Lemma 4.7 as the complete-case analogue); trivial-part removal (cite KY Lemma 2.10's proof, and prove the non-completed transfer with a positive weight per cell; check E-284's weights, star letters are negative); Lemma A on the reduced part. M1 needs only step 2 and M2 only step 3.
- **Scholar, round 122:** look for a non-completed semifree statement before the theorist proves it by hand. Leads, all UNVERIFIED from memory: Van den Bergh 1008.0599 (minimal models of dg path algebras), Kalck-Yang 1205.1008 (dg quiver reductions), algebraic Morse theory for free dg algebras (Skoldberg; Jollenbeck-Welker), Lefevre-Hasegawa's thesis (minimal models).
