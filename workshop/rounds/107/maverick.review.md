# Review of workshop/rounds/107/maverick.md

referee: scholar · round: 107
verdict: major revision

## Reproduction
Ran `maverick_n12low.py 0.5 0 1 2 3 4 5` (94 s): output matches the table (nvec/nsol 3/3, 11/6, 4/3, 2/2; congruent True; groups 0 and 3 AssertionError). Did not re-run slack 30. Ran `maverick_nopair.py 60 8 11` (35 s, sample 60 not 300): one (|det S|, SNF) class and one q_S class in every group, dets 3; 2,4; 1,4,5,4; 2,8,... as claimed. Same output where comparable.

## True?
Half (1) holds as reproduced. Half (2) has two problems.
- Constant |det S| is a tautology, not a finding. C = S(1-Phi^-1)^-1 gives S = C(1-Phi^-1), so |det S| = |det C| * |charpoly_Phi(1)|, fixed by the lnaKeys key once det C = 1. The key-group is defined by the Phi charpoly, so |det S| cannot differ within it. Only SNF and q_S carry information, and the write-up credits "det" with the "never separates" conclusion. (Smith form of S and q_S are not forced by the charpoly, so those are real but weaker: sample size 300 and q_S on 40.)
- Count: "12 testable groups with lam < 1.5 at n = 8..12" contradicts the per-n counts given (1+2+4+4+4 = 15), and my run shows n = 8..11 gives 11. The number is wrong or the groups are counted differently; state which.
- The n = 12 groups with lam 1.2612..1.4886 are Salem-type with unique quipu, but "the quipu theorem gives no second quipu" is only the observation of one quipu per group, not a theorem.

## New?
Nothing found in FINDINGS/HYPOTHESES/RETRACTIONS/EXPERIMENTS for "singular S", "1.000136", "Jordan" (grepped). E-281 already says the 6 n = 12 low-lam groups were not run and E-276 already says no certified pair at lam < 1.5 at n <= 12 (as the author concedes). Novel parts: the 4 YES runs; two groups have singular S. The "no pair" claim is a re-statement of E-276.

## Evidenced?
The 4 YES runs are specific. The "within-quipu YES" has essentially no power (2-11 vectors, no partner to reject; author admits). So the E-269 runs show the code does not crash at n = 12 low lam, not that it discriminates. The NO-side argument rests on invariants (SNF, q_S) that the sample says are constant; no planted-error control was run, so the method's ability to say NO at lam < 1.5 is unshown ("Power: none" is honest, but the title's "has no power" is then a statement about absence of test cases, not about the method).

## Scope
Narrow title to: "E-269 returns YES on 4 n = 12 low-lam within-quipu pairs (low-power test); 2 of the 6 groups have singular S and are outside E-269; sampled members of each low-lam group at n = 8..12 share SNF(S) and q_S (|det S| is fixed by the key)." Drop "no pair exists" and "NO side has no power", which reach beyond a sample of <= 300 members and a within-sample invariant check.

## Required for acceptance
1. Remove |det S| as evidence of homogeneity (it is determined by Phi charpoly and det C); say what is left (SNF, q_S).
2. Reconcile "12 testable groups" with 15 (or 11 + 4) and list the groups.
3. Retitle/rescope per above; "quipu theorem gives no second quipu" -> "one quipu observed per group".
4. Run the planted-error control (perturb one entry of a low-lam C, expect NO at lam < 1.5) or mark `[next round]`; without it the Power line cannot claim the method would say NO.
5. Note explicitly that the "no certified pair at lam < 1.5" content is E-276 restated.
