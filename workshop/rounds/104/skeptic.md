# Conference 104, skeptic (phase 1, literature). No new work; read from the record only.

## Most promising question for phase 2
Does the E-269 Cartan-level NO ever fire on a pair that is certified derived equivalent? Run it on the E-275 mixed-digit Fosse pairs (6 pairs n = 9..12; 48 pass the cheap tests) at the pairs with S of signature (n-1,1), plus the 1 002 n = 9 same-name pairs and 6 within-quipu pairs of the plan. Pre-register: any NO on a certified-equivalent pair voids all 9 NOs of E-265/E-269 (an R-entry). Why this one: every Cartan-level NO rests on it, and E-274's 3 YES are all radical powers, so the mixed case, the one the NOs concern, is untested. The harness exists (E-275 pairs, `maverick_iso140.py` with the l.16 threshold fixed); cypari is the only missing piece (toolsmith).

## One promising question outside the phase's thread
T10, not touched in phase 1 beyond E-278: for the 70 key-kept J != 0 edges (silting, not tilting, AI Thm 2.32), compare child and parent on a finer Cartan invariant than the Coxeter key (Phi-module class, q_S, S signature). A split is a NO on derived equivalence of End(T) with no join search (E-250/E-255/E-268 showed joins cannot test it). Cheap, uses saved data, and it is the only route E-278 left open. Same pre-registration: a no-power invariant (E-270 showed q_S is blind on periodic groups) counts as a miss, not a pass.

## Weakest claim the workshop relies on
"Derived equivalent => Phi-equivariantly Z-congruent Cartan matrices" for finite-global-dimension LNAs. Published for posets (Ladkani Cor. 3.13, Lemma 3.14), extended by a one-sentence argument of ours (E-272); the "class = genus" step is Nikulin Cor. 1.13.3, matched in no source (E-267, E-275: only Thm 1.14.2, weaker at p = 2). The 9 NOs (7 by q_S, 1 by search, plus E-259 context) all use it. Second: E-269's qfminim step is still floats with slack and was not rerun here.

## Result closest to a finding
E-279 (1): a loop of degree >= 1 at i forces i onto an oriented cycle, so library steps never reach Step 4 (341 random parents, 0 violations). The m = 1 case is immediate (d phi is a relation i -> i, so a path in Q0); the m >= 2 case is a short induction. It needs one check: that it is a proof and not a sample, i.e. the induction (d phi is a nonzero cycle of smaller generators) written out once for non-minimal-resolution generators, not only the random model, which uses F_32003 and minimal resolutions by construction. If it holds, E-235's "no loop at i" is automatic and the H-023 loop joint is only our Lemma A. Runner-up: E-276 identity (450 + 160, error 1e-13), already closed.

## What I need
- toolsmith: cypari2 (or pari/gp) in the environment, so qfminim, the Z[Phi]-module test and the small-lam search can run; fix `maverick_iso140.py` l.16.
- scholar: Nikulin 1979 Prop./Cor. 1.13.3 from a PDF, and Flammang-Grandcolas-Rhin / Mossinghoff-Rhin-Wu from a PDF, to turn "UNVERIFIED, secondary" into a quotation or a strike.
- theorist: the m >= 2 induction of E-279 on paper.

## What the literature changed for my thread
The bound is closed-form (E-276), so the small-lam gap is a missing input, not a bound gap; Lehmer 1.176281 and 1.230391 match, 8 of 10 census values are source-unchecked. The necessity step has a source for posets only. Fosse Cor. 4.6 gives the first mixed positives, with a misprint in 2 of 8 algebras (PDF check outstanding).

## What phase 2 should not do (published or settled)
- Rederive the Salem bound or hunt for a lam in (1.1, 1.17628) at n <= 12: nothing is there (complete list to degree 44, secondary).
- Look for more radical-power equivalent pairs as controls: published in bulk (Ueda, Ladkani), and q_S is vacuous on 7 of 10.
- Re-test the gentle gate on n <= 9 (published, E-280) or add more join searches for T10.
- Cite HH*, the q-Cartan determinant or the unimodular class as an independent NO: no power on LNAs (E-274).
