# E-323's F3 step at 0040033030·F3 R9 F3 is non-tilting at Hom level (Hom(T,T[-1]) = 2, at p = 2, 7, 101, 32003), and 0040033030 / 0400330030 are separated by Phi-conjugacy invariants that need no lattice theory; all of them sit at eigenvalue 1 mod 3 and mod 5 and go through the same necessity chain as SNF(S), one ((x - 1)^2 at p = 3) being finer than SNF(S)

author: skeptic · round: 129 · kind: result
thread: T10 · bears on: E-323, E-318, E-324, H-024, F-055, H-032, R-017, the F-candidate of conference 128
scope: n = 12, the one E-323 site (root 0040033030, parent P = root·F3 R9 F3, step F3, and its dual mirror R3); every gate-admitted F/R step at the root and at P (35 cases); the pair 0040033030 / 0400330030 and, as a cross-check, E-318's n = 11 pair. Invariants: rank over F_p (p <= 13) of f(Phi) for 23 polynomials f, plus charpoly. Not covered: the other 49 crossing steps of E-323 at Hom level; whether End(T) is isomorphic to the library child as an algebra (only its Cartan matrix, from Hom dimensions, is compared); any invariant not determined by the Cartan matrix.
power: Hom test: says NO on E-080's algebra abde = acde at d (known non-tilting; Hom(T,T[-1]) = 1) and YES on abd = acd and on all 33 J = 0 steps (3 walk steps, 2 last steps, 28 neighbours); it can and does say NO. Invariant: YES on 5 of 5 Hom-certified tilting steps of the walk (End(T) Cartan from Hom dims) and on 10 of 10 congruences r C r^T, r in GL_12(Z). NO-power: none independent of the necessity chain -- every certified NO on an LNA pair goes through the same chain; the planted det-2/3 base changes are caught only by the unimodularity check, not by the invariant (no power from them).

## Claim

(1) At P = 0040033030·F3 R9 F3 the Okuyama-Rickard complex T for F3 (and for the dual mirror R3) has Hom(T, T[-1]) of total dimension 2 and Hom(T, T[m]) = 0 for m = -2, 1, 2, computed by the E-327 Hom engine (mod 32003, independent of the gate, of perI and of `isTiltingAtVertex`). So the step is not a tilting mutation by direct Hom computation, not only by the injectivity test. Hom dimensions of T give exactly the library child's Cartan matrix, so the child has the Cartan data of End(T) for this non-tilting T. On all 35 gate-admitted F/R steps at the root and at P, J != 0 <=> Hom non-tilting <=> `isTiltingAtVertex` False (3 J != 0: F3 and R5 at P, mirror R3; R5 changes the key and is guard-refused). (2) The pair has equal charpoly(Phi), and 15 of the 126 (f, p) ranks of f(Phi) separate it. 13 of them are not f = x - 1. But each of the 13 reduces mod p to a multiple of (x - 1): x + 2, x - 4 at p = 3; x + 4 at p = 5; (x - 1)(x + 1), x^k - 1; x^2 + x + 1 = (x - 1)^2 mod 3. f = x + 1 (the skew form C - C^T) and every f prime to x - 1 mod p agree. So there is a second computation but no second argument, which is E-324's conclusion for E-318's pair, now at n = 12. (3) The NO direction does not use Nikulin 1.13.3. Derived equivalence gives C_B = r C_A r^T with r in GL_n(Z) (Ladkani 1001.4765 Lemma 3.5 + Aihara-Iyama Thm 2.27 + Rickard, E-308; math/0610685 Cor. 3.13, E-272). Then S_B = r S_A r^T, and Phi_B = r^-T Phi_A r^T, so rank_p f(Phi) and SNF(S) = coker(Phi - 1) are equal. Each step is sourced or elementary. Nikulin is the class = genus step and serves only YES answers. E-272 and E-308 already say no NO depends on it. The premise of the conference-128 blocker ("rests on Nikulin 1.13.3") is therefore a misreading. Refuted by: a Hom-tilting computation at P·F3, or a GL_12(Z)-conjugation between the two Phi.

Does not claim: that "non-tilting T" alone makes P and its child inequivalent. That still rests on the SNF / rank invariants and the necessity chain. Nor that End(T) is isomorphic to the child, nor anything about the other 49 steps.

## Evidence

Hom test (`skeptic_hom323.out`):

| case | J != 0 | isTiltingAtVertex | Hom-tilting | Hom(T,T[m]) != 0 | End(T) Cartan = child |
|---|---|---|---|---|---|
| E-080 abde = acde at d (NO control) | - | - | False | m = -1: 1 | - |
| abd = acd at d (YES control) | - | - | True | none | - |
| walk F3, R9, F3 from 0040033030 | F, F, F | T, T, T | T, T, T | none | T, T, T |
| **E-323 step F3 at P** | **True** | **False** | **False** | **m = -1: 2** | True |
| E-323 mirror R3 (dual side) | True | False | False | m = -1: 2 | True |
| P·R5 (key changes) | True | False | False | m = -1: 1 | True |
| other gate-admitted steps at root (14) and P (14) | all False | all True | all True | none | all True |
| last step X -> 0400330030 (F7, F9) | F, F | T, T | T, T | none | T, T |

35 cases + 2 last steps: 0 disagreements among J, injectivity, Hom. Mechanism at P: E-080's shape, read off the relations (one arrow 3 -> 8 out of 3; e_6 A e_3 = <6-7-3, 6-9-3>, input relation (6 7 3 8) + (6 9 3 8)). Read by eye, not classified by script.

Invariants (`skeptic_inv323.out`), ranks of f(Phi) over F_p, A = 0040033030, B = 0400330030:

| f | p | rank A | rank B | reduces to (x-1)-multiple mod p |
|---|---|---|---|---|
| x - 1 (= SNF(S)) | 3, 5 | 11, 10 | 10, 11 | itself |
| x + 2, x - 4 | 3 | 11 | 10 | yes |
| x + 4 | 5 | 10 | 11 | yes |
| x^2 - 1, x^3 - 1, x^4 - 1, x^6 - 1 | 3 / 5 | 9,9,9,7 / 8,8,8,6 | 8,8,8,6 / 9,9,9,7 | yes |
| x^5 - 1 | 3 | 11 | 10 | yes |
| x^2 + x + 1 | 3 | 10 | 8 | yes ((x-1)^2) |
| x + 1, x^k + 1, x^2 - x + 1, all f at p = 2, 7, 11, 13 | | equal | equal | |

Along the walk the 126 ranks are constant on the 3 tilting steps, change at 15 (f, p) at the F3 step, and are constant on the last tilting step: the J != 0 step is the only place the class changes. E-318 n = 11 pair: separated at p = 2 by x^2 +- 1, x^6 +- 1 (= (x - 1)^2 mod 2), consistent with tails (2,4) vs (2,2), which plain x - 1 ranks cannot see.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/skeptic_hom323.py   # 3 s, out skeptic_hom323.out
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/skeptic_inv323.py   # 5 s, out skeptic_inv323.out
```

## Prior record

E-323 and its review (round 125, by me) used `isTiltingAtVertex` and `tiltingPlus`. Both are the injectivity test (Ladkani 2.3(c)), and neither is a Hom computation. So (1) is new: the first Hom-level check of an E-323 or E-318 crossing step. E-327 showed injectivity = Hom on random parents, not on this one. (2) at n = 12 is E-324's "no logically independent Cartan-level second witness" (n = 11) repeated: a rediscovery, now explained (everything factors through (x - 1) locally). (3) is in E-272 and E-308 ("no NO result depends on Nikulin"). F-055's "What it does not say" puts Nikulin next to q_S and can be read the other way, and STATE inherits that reading.

## Code changed

None in the library. New scripts `workshop/rounds/129/skeptic_hom323.py`, `skeptic_inv323.py` and their `.out` files. No tests touched.

## Next

- Chair: the E-323 blocker as named is met. The step is Hom non-tilting, and inequivalence needs only the sourced chain, not Nikulin. A non-Cartan invariant does not exist here (HH trivial, E-324), so asking for one blocks permanently. Write the F-candidate with "inequivalent = Cartan matrices not GL_n(Z)-congruent (E-308 chain)" in its scope.
- Toolsmith: run `skeptic_hom323.check` on the other 49 E-323 steps (one line per step, seconds each). It gives Hom-level coverage of the whole finding and tests E-327's agreement on walk-reachable parents, the R-candidate's blocking check.
- Theorist (agenda item 1): at P, End(T) has the child's Cartan matrix while Hom(T,T[-1]) = 2. Is the library child End(T) of this silting-not-tilting T? Compare the algebras, not only their dimensions.
