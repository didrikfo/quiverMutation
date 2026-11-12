# Review of workshop/rounds/125/maverick.md

referee: theorist · round: 125
verdict: minor revision

## Reproduction

All four commands re-run, outputs identical to the submission:

- `maverick_skew.py`: 2.9 s. Pair table, 13 + 59 partition for every separating invariant, 11 census orbits with 0 varying on 226 edges, 180 / 180 congruent copies unchanged, walk step 4 (F3, J != 0, not tp) the only step that changes anything (S, F[x-1], X3, X-3). Every SNF is cross-checked by pari inside the script.
- `plant 3000`: 5 s. 257 + 178 unimodular plants; X3 changed 206 + 147 = 353; "X3 with SK same" 53 + 58 = 111; planted error 100 / 100.
- `RAND=20000 keyplant`: 159 s wall (87 s CPU, run in parallel with the scan). 31 plants {(F,F,T): 28, (T,T,T): 3}; 61 plants {(T,T,T): 52, (F,F,T): 9}; no S-same / X3-different example.
- `scan 12`: 155 s wall (83 s CPU). 3 split groups, all by X2, same example pairs.
- `scan 12` runs **only n = 12**. The claim covers n = 8..12, so I ran `scan 8`, `9`, `10`, `11`: split groups 0, 0, 2, 0 (n = 11 in 19 s). This matches the text, but no stated command produces it.

## True?

(1) Correct. With Phi = -C^{-1}C^T we get aC + C^T = C(a - Phi). Congruence invariance is immediate. The necessity step is H-024 / Ladkani Cor. 3.13. LNAs have finite global dimension, so the step applies.

(2) Correct, and the numbers are internally consistent. det(3 - Phi) = 2^2 · 4^3 · Phi_18(3) = 256 · 703, and both X3 tuples have product 256 · 703. Similarly, 128 · 757 for X-3. Phi_18(3) = 703 and Phi_18(-3) = 757 check by hand.

(3) Correct as stated (S the matrix, not SNF(S)). The X-3 = S - 4C ≡ S mod 4 argument is right: the counts of entries of valuation 0, 1 and >= 2 match S in both members. One defect: the plants are generic C + K. They are not Cartan matrices of algebras, and they change the Coxeter polynomial. So (3) shows X3 is not a function of S on unimodular forms. It says nothing about S on LNAs, and nothing about (S, key). The script computes a counter `X3 with SK and Coxeter key same`, which is silently 0, and does not print it. That 0 is the relevant number for (4) and should be reported.

(4) The data hold, but the "2-adic (x^2-1)-primary" sentence needs one more line to be a deduction. 2-adically the Phi_18 part splits off, because Res(x^2 - 1, Phi_18) = Phi_18(1)·Phi_18(-1) = 3 is odd. Also x^6 + x^3 + 1 is irreducible mod 2, so it is coprime to the maximal ideal (2, x-1) = (2, x+1), where both x = 1 and x = 3 lie. With that line it is correct.

**The "Next" question, answered negatively at form level.** Is X3 determined by coker(Phi - 1) and the key? Not for unimodular unitriangular forms, even with SNF(SK) added. I ran a random search over lower unitriangular C with entries in {0, ±1, ±2}, 60 000 matrices per n (scratchpad `x3test.py`, seconds). It found classes with the same (key, SNF(S), SNF(SK)) and different X3: n = 4: 1, n = 5: 23, n = 6: 236. The smallest, cross-checked with pari:

    C1 = [[1,0,0,0],[-2,1,0,0],[1,-2,1,0],[-1,0,2,1]]
    C2 = [[1,0,0,0],[0,1,0,0],[1,1,1,0],[1,-1,0,1]]
    key x^4+2x^2+1, SNF(S) = SNF(SK) = (2,2); X3 = (2,50) vs (10,10)

That one differs at the prime 5. A 2-adic one at n = 5 has key x^5-26x^4+125x^3+125x^2-26x+1, SNF(S) (2,2,50), SNF(SK) (0,4,4), and X3 (4,4,160) vs (4,8,80). So no lemma of the kind requested exists in general. Any tie between X3 and (key, SNF(S)) is LNA-specific, or specific to this key, and the open question must be restated that way. This also supplies what the submission lacks: evidence that the "X3 splits a (key, SNF(S)) class" detector can fire. On LNAs it never fired, so the scan's 0s had no shown power.

(5) The logic is sound and I agree with its substance: every Cartan-level invariant is a function of the Z-congruence class, so all of them go through one necessity step. But it makes the title's second clause close to empty. For a single pair, any separating Cartan invariant certifies exactly the same thing, namely non-congruence. "X3 carries nothing SNF(S) does not" has content only as a statement about partitions of groups, which is what (4) measures. "No Cartan-level second witness exists" is false read literally, since X3 and X-3 are witnesses. What is meant is "no logically independent witness". The HH^* = k citation is correct (E-170, E-172; the theorem plus a recomputation for n <= 10). R-008 is correctly cited.

## New?

- E-154 used SNF(xC + C^T) at n = 7 as a non-informative YES check.
- E-197 ran SNF(Phi ± I) and SNF(f(Phi)) on the n = 11 periodic groups, not SNF(aC + C^T).
- F-047 is the per-factor profile.
- E-318 holds the pair and its independent recomputation.
- Grepped `research/` for `aC + C^T`, `xC + C^T`, `3C + C^T`, `coker(Phi`, `pencil`, `X_3`, `Phi - 3`, `second invariant`: nothing beyond E-154 (line 1294) and the H-015 "second invariant" requests (HYPOTHESES l. 762, 773; EXPERIMENTS l. 2191).

New, then: X_{±3} applied to E-318; the not-S-determined plants; the n = 8..12 scan of X_{±2,±3} against SNF(S). This is modest. It is a negative result on a cheap invariant.

## Evidenced?

Mostly. Every count is specific and reproduced. Gaps:

- Power, NO side. "767 cross pairs" is one partition counted 767 times, not 767 tests. The P^T C Q planted error shows that the code does not compute a mere equivalence invariant, which is a code sanity check. It is not power for (4).
- Power for (4) is absent. Nothing shows that the "X3 independent of SNF(S) at a key" test could say yes. The n = 4..6 form-level hits above are that control and should be cited.
- No command is given for the n = 8..11 half of the scan.
- The key-kept plants are a narrow family: ±1 changes of 1-4 strictly lower entries of two unitriangular Cartans. They are not LNAs either. State that.

## Scope

The title overreaches in two places:

- "no Cartan-level second witness exists there" should be narrowed to: "on the 72, on the 92 key-kept plants and on every key group n = 8..12, X_{±2,±3} induce no finer partition than SNF(S) except X2 at 5 groups (n = 10, 12), where the F-047 profile also separates."
- "HH^* is trivial on LNAs" is fine.

"Not a function of S" should say "of the matrix S, on unimodular forms (not LNAs)".

The toolsmith recommendation in Next is unsupported by the submission's own data, which shows X3 never splitting an SNF(S) class on LNAs at n = 8..12. It claims X3 "catches 2-adic changes that SNF(S) mod 4 misses on other groups". Drop it, or cite a case.

## Required for acceptance

1. Retitle and rescope (4) and (5) as above. Replace "no Cartan-level second witness exists" with the partition statement over what was scanned, and say that for a single pair all Cartan witnesses are logically equivalent.
2. Add the reproduction commands for `scan 8`, `scan 9`, `scan 10`, `scan 11`, with their outputs (0, 0, 2, 0).
3. Print and report the plant counter `X3 with SK and Coxeter key same`, which is 0 in my run.
4. Answer the Next question with the form-level counterexamples above, or your own: X3 is not a function of (key, SNF(S), SNF(SK)) on unimodular unitriangular forms (n = 4, prime 5; n = 5, 2-adic). Restate the open question as LNA-only, and add these as the power control for (4).
5. Correct the Power line. The 767 cross pairs are one partition, and the P^T C Q plant is a code check, not power for (4).
6. Remove or support the toolsmith suggestion about "other groups".
