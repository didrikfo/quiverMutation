# A second, independent nu (termwise nu + cone replacement; no bimodule resolution) reproduces E-346 at n = 14-15 object by object: kA_15/rad^3 nu^144 S_5 ~ S_5[160], P_2 ~ P_2[162]; kA_15/rad^8 160/162; kA_14/rad^5 70/72 (m = 60), and the uniserial certificates of N_14(4), N_14(6) and three n = 10-11 groups; GF(2^31-1) and GF(2) (n = 15 certificate objects also GF(3), GF(1000003))

author: maverick · round: 135 · kind: result + tool
thread: T3/T8 · bears on: H-025, E-346, E-339, E-245, E-246
scope: the 8 algebras of E-346's table and its 6 controls at the m of E-346. N_15(3), N_15(8) at m = 144 and N_14(5) at m = 60 on all simples and projectives. N_14(4), N_14(6) at m = 132, 050404400 at m = 24, 48, 040404000 at m = 72, 45050400 at m = 84, 34504030 and 50505000 at m = 9, 18, 334504030 at m = 21, all on every uniserial. Exact over GF(2^31-1). N_15(3) (all S, P, m = 144) and N_14(5) (all S, P, m = 60) also over GF(2), with identical output. N_15(3) S_5, S_11, P_2, P_14 at m = 144 also over GF(1000003), and S_5, P_2 over GF(3). Uniserials only over GF(2^31-1). Char 0 follows by E-346's base-change argument, which I did not re-check. m = 72·j, j > 2 and other primes not run. This checks the nu computation, not the deduction, which E-346's review already checked.
power: four controls. (1) Hypersurface N_14(3) (FCY is the truth). Every one of 27 simples and projectives and 39 uniserials returns with the single shift 26 at m = 24. N_11(4) gives 38 of 38 with shift 16. (2) Wrong nu = identity. It is a genuine autoequivalence and passes every internal assert, and the certificate says "one shift 0, consistent with FCY". The Serre-duality check rejects it at step 1, and the comparator reports 648 of 648 (object, k) mismatches. (3) Wrong injectives. "I_i one vertex short", "I_i one vertex long" and "injective length = projective length" each fail the module-map assert or d^2 = 0 at the first or a later step. (4) Non-functorial scalars: length-2 standard maps times 2. This fails d^2 = 0 on N_14(5). On rad^3 it is invisible, because it is gauge-equivalent there (see Claim).

## Claim

I wrote `maverick_nu2.py`, a second nu on D^b(mod A) for linear Nakayama A. It shares no code with `rounds/089/maverick_bimod.py` or `rounds/087/scholar_rankgrowth.py`, and it uses a different method. A minimal complex of projectives X is "interval types + scalar matrices", because Hom(P_i, P_j) has dim <= 1. nu is applied termwise: P_i -> I_i, and the standard map -> the standard map. This is nu up to rescaling the isomorphisms nu P_i ~ I_i, since every Hom is spanned by a composite of arrows. The complex of injectives is then turned back into projectives by the cone algorithm, at the level of representations: descending in degree, P^k is the projective cover of the cone cycles W^k = {(x, y) : dx = 0, fx = dy}. Gaussian elimination of iso entries makes it minimal. No bimodule resolution of D(A) and no bar resolution is used.

**Result.** For every (object, k <= m) on N_15(3) (m = 144), N_15(8) (m = 144) and N_14(5) (m = 60), the two implementations agree on the full table dim H^t(nu^k X)_v for every v and t. That is 4176 + 4176 + 1620 (object, k) pairs and 0 disagreements (`maverick_compare.py`). So the E-346 certificate objects are reproduced: S_5, S_11 return with shift 160 and P_2, ..., P_14 with 162 (N_15(3)); 160/162 (N_15(8)); 70/72 (N_14(5)). The same holds for the intermediate shifts 80/81 at m = 72. On every uniserial, the new code alone reproduces each row of E-346's uniserial table: the same returning objects with the same shifts. That is 154/156 for N_14(4) and N_14(6), 28/30 and 56/60 for 050404400, 80/90 for 040404000 and 98/108 for 45050400. The three no-certificate groups keep a single shift (10, 20; 10, 20; 24). The Serre-duality identity dim H^t Hom(X, P_i) = dim H^{-t}(nu X)_i needs no nu code. It holds at every step on N_15(3) (all 29 objects, k <= 144) and on N_14(5) (all 27, k <= 60). The "one implementation of nu" blocker in E-346 is therefore removed at n = 14-15, at GF(2^31-1) and GF(2) (GF(3), GF(1000003) on the n = 15 certificate objects).

**What a wrong nu looks like.** A termwise error on objects, such as wrong injective envelopes, cannot survive: the standard maps stop being module maps, or d^2 != 0 (controls 3). A scalar error on morphisms is either a gauge change (rad^3: every relevant complex is tree-like, so the output is unchanged) or it breaks d^2 = 0 (rad^5). The dangerous wrong nu is a different *coherent* autoequivalence F. F passes every internal assert, and the certificate then decides FCY-ness of F, not of the Serre functor (control 2: identity, a false "consistent"). The Serre-duality check is exactly the test that tells the Serre functor apart from such an F, and both implementations pass it. **Not claimed:** other primes with the new code, or uniserials at primes other than 2^31-1; m beyond those listed; a proof that the two codes cannot share a conceptual error. They share only the Kupisch series, which is asserted equal, and the convention "X ~ [s, e][c] iff H^* is the dimension vector of [s, e] in one degree". **Refuted by:** a third implementation, or a hand computation, that gives nu^144 S_5 of N_15(3) cohomology outside degree -160.

## Evidence

| algebra | m | objects | new nu: returns, shifts | vs. E-346 | object-by-object vs. nu_bimod |
|---|---|---|---|---|---|
| N_15(3) 3333333333330 | 72 / 144 | 29 S, P | 9: 80 (S5, S11), 81 (P2..P14) / 160, 162 | same | 4176 / 4176 (k <= 144), Serre ok |
| N_15(8) 8888888000000 | 144 | 29 S, P | 160 (P3, 6, 10, 13), 162 (S2, 5, 8, 11, 14) | same | 4176 / 4176 |
| N_14(5) 555555555000 | 60 | 27 S, P | 70 (S4, S11, P3, 4, 7, 8, 11, 12), 72 (S2, S13) | same | 1620 / 1620, Serre ok |
| N_14(4) 444444444400 | 132 | 50 U | 154 ([4,5], [10,11]), 156 (seven) | same | not compared (141 s) |
| N_14(6) 666666660000 | 132 | 69 U | 154 ([6,7], [6,9], [8,9]), 156 (seven) | same | not compared (275 s) |
| 050404400 | 24 / 48 | 43 U | 28/56 ([6,8]); 30/60 ([1,5], [2,3], [5,9], [9,11]) | same | -- |
| 040404000 | 72 | 42 U | 80 (four), 90 ([2,3], [9,10]) | same | -- |
| 45050400 | 84 | 39 U | 98 ([4,5], [4,7], [6,7]), 108 ([1,3], [3,8], [8,10]) | same | -- |
| 34504030 / 50505000 / 334504030 | 9,18 / 9,18 / 21 | 34 / 42 / 37 U | one shift: 10,20 / 10,20 / 24 | same | -- |
| control N_14(3), N_11(4) | 24, 15 | U and S, P | all return, 26 / 16 | same | 648 / 648 (N_14(3)) |
| wrong nu = id on N_14(3) | 12, 24 | 27 S, P | all return, shift 0 ("consistent") | -- | 0 / 648; Serre fails at k = 1 |
| I_i short / long / p-length | | N_14(3) | assert: not a module map / d^2 != 0 / not a module map | -- | -- |
| scalar plant (len-2 maps x2) | 24 / 60 | N_14(3) / N_14(5) | unchanged (gauge) / d^2 != 0 | -- | -- |

A bug found on the way: my first version of the d^2 check multiplied int64 matrices mod 2^31-1 and overflowed. It raised a false "d^2 != 0" on 050404400 at k about 30. It is fixed with a split multiply (`mulmod`), and all runs above are from the fixed code. The row-reduction itself never multiplies more than two residues, so it was unaffected. Anyone writing GF(2^31-1) code with numpy should watch for the same thing.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 3333333333330 72,144 S,P        # 43 s (54 s with --serre)
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 8888888000000 144 S,P           # 57 s
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 555555555000 60 S,P --serre     # 15 s
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 3333333333330 144 S5,P2,S11,P14 --prime 1000003   # 2 s; --prime 3 S5,P2: 1 s; --prime 2 S,P: 37 s (555555555000 60: 11 s)
timeout 10m .venv/bin/python workshop/rounds/135/maverick_compare.py 3333333333330 144           # 56 s; 8888888000000 144: 70 s; 555555555000 60: 16 s
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 444444444400 132 U              # 141 s; 666666660000 132 U: 275 s; 45050400 84 U: 556 s;
#   050404400 24,48 U: 238 s; 040404000 72 U: 299 s; 34504030 9,18 / 50505000 9,18 / 334504030 21 / 333333333330 24 / 444444400 15 (U): <= 20 s
#   all of these -> workshop/rounds/135/maverick_nu2_unis.out
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 333333333330 12,24 S,P --plant=id      # wrong nu; add --serre (m = 1) to see it rejected
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 333333333330 12,24 S,P --plant=trunc   # also long, kup, scale (< 3 s each)
timeout 10m .venv/bin/python workshop/rounds/135/maverick_compare.py 555555555000 60 S,P --plant=scale  # d^2 != 0
timeout 10m .venv/bin/python workshop/rounds/135/maverick_compare.py 333333333330 24 S,P --plant=id     # 0 / 648 agree, 1 s
```

## Prior record

E-346 (one nu implementation at n = 14-15; the review's bar-method attempt at n = 14 took > 5 min). E-231/E-245 (rounds 087/089): the bar method and the bimodule method were already cross-checked at n = 10-11, with a Serre-duality check (E-245 scope line). Not new: the shifts and verdicts themselves. They are E-346's. New: an implementation that is independent at n = 14-15, the full per-step agreement there, the Serre check at n = 15, and the wrong-nu controls (in particular that a coherent wrong autoequivalence is caught only by the Serre check, not by any internal assert).

## Code changed

New files: `workshop/rounds/135/maverick_nu2.py`, `maverick_compare.py` (imports rounds/089, 087, 134 unchanged), and the output `maverick_nu2_unis.out` (2.5 KB). No library edits, so no tests apply.

## Next

- Chair: E-346's blocker "one implementation of nu at n = 14-15" can be lifted, with H-025 -> CONFIRMED for those instances if the referee agrees.
- Toolsmith: `maverick_nu2` is about 10x slower than nu_bimod at n = 15, but it needs no precomputed resolution and runs at any n. Promote one of the two into the library, with the Serre check as a standing test.
- Theorist (idea): every planted error was either caught by an internal assert or was a gauge change, except a coherent autoequivalence. So "is this the Serre functor?" reduces to the Serre-duality identity on Hom-dims. Is the identity on all S_l and P_l at every step sufficient to pin nu up to isomorphism on these orbits?
