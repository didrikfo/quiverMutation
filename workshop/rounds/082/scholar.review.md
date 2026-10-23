# Review of workshop/rounds/082/scholar.md

referee: theorist · round: 082
verdict: minor revision

## Reproduction

`timeout 600 .venv/bin/python workshop/rounds/082/scholar_cyclo.py 11` took 19 s (the submission says about 1 min). The output matches the submission's table in every row: four groups of 535, 64, 24 and 15 (total 638), ord Phi 30, 72, 24 and 42, lcm phi(m) 8, 12, 4 and 6, and windows (15 odd, 30 even), (72 even), (24 even) and (21 odd, 42 even). It also reports radical powers L = [4, 5], all in the 535-group. I did not re-derive the 638-LNA census, which is upstream.

I fetched both LaTeX sources and checked the following.
- 1310.1557: Thm 2(a), 2(b) and Lemma 2.9 as stated. The 2(b) converse is only for "every component of [Gamma] is a tube". The FCY definition is `S^q=[p]`.
- 0911.5137: `c:rectline` and the unnumbered line corollary `(m-1)/(m+1)+(n-1)/(n+1)`. The sign argument checks out: [S] = -Phi and [1] = -1, so (-Phi)^a = (-1)^b, i.e. Phi^a = (-1)^(a+b) I. I also checked each window by hand. For example, at a = 21 in the 15-group, the primitive 2nd, 6th and 14th roots all go to -1.

## True?

Nothing false found. Three points are loose.
1. "The only FCY LNAs the papers hand us are the lines" overreaches. 0911.5137 Cor `c:fracCY` says Lambda being lambda-CY makes T_n(Lambda) (lambda + (n-1)/(n+1))-CY, and a tensor product of FCY algebras is FCY. The "only lines" claim needs a stated argument. The argument is the prime rank: 11 = rank(Lambda)·n forces a trivial factor. The submission never says this. As written, the claim is "the lines", not "the only FCY LNAs the papers give".
2. Cor 1.9 at n = 11 is not "nothing". Taking m = 1 or n = 1 gives A(11,2) ≃ kA_11, which is FCY but Dynkin and so outside the indefinite set. This is harmless, but "gives nothing" should read "gives only Dynkin".
3. "The 64 and 24 groups allow a = 72, 24 only" holds for the minimal a. Any multiple of the period also satisfies Phi^a = +I. The table is of minimal a, and the text should say so.

## New?

Mostly recorded already.
- The misprint (lcm phi(m) vs lcm m) is in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` lines 37-42, checked on 2026-09-19 against N_15(5). The submission admits this.
- "Cor 1.2 and 1.9 are silent for prime N" is in `research/literature/0911.5137-lines-rectangles-triangles.md` lines 137-139.
- "FCY implies periodic, not the converse" is in the 1310.1557 summary lines 54-57 and 118.
- The four chi factorisations are in E-190 / E-197 / E-209.

Not found by grepping research/ for `Phi\^(15|21|42|72|24|30)`, `a = 15|21|42` and "half period": the per-group Phi^a = ±I windows at n = 11 (a = 15 or 30, 72, 24, 21 or 42), and the observation that only the 535 and 15 groups admit odd a+b. E-190 has this only for n = 10. The theorem-number audit is also new but is only bookkeeping. That 535 contains kA_11/rad^4 and rad^5 is new at idea level. I did not look for a record of it.

## Evidenced?

The literature claims are specific and I verified them. The weak part is the headline "no theorem there reaches n = 11". It is negative evidence from two papers only. The title scopes it to those two, which is acceptable. The "Next" test is a proposal and is not evidence.

## Scope

The title is fine as a statement about two papers and does not generalise. It should say it covers the two papers only, not the literature. For example, Happel-Seidel and Lenzing-de la Peña are noted but not checked for FCY of rad^4 and rad^5 at n = 11. The submission marks this UNVERIFIED, which is honest.

## Required for acceptance

1. Add the prime-rank argument that excludes `c:fracCY` / tensor-product / T_n(Lambda) constructions at n = 11, or weaken "only" to "the lines (and, by primality, nothing else from the tensor and T_n constructions)". Do the same for the Dynkin case, A(11,2) ≃ kA_11.
2. State that the windows are minimal a (multiples also work) and that the 535 and 15 groups are the ones with an odd-a window.
3. Label the lcm-misprint finding as already recorded (literature summary lines 37-42), not new.
4. Fix the runtime ("about 1 min" is 19 s). Optional.
