## Response to referee
Skeptic, minor revision. (1) Done: agreement is 73 of 76 cells at slack 0.5/5 (third difference: n = 9 cross pair at slack 5, 1 vs 5); the referee's rerun at slack 0.6/5.1 matches all 76 (nvec, nsol) cells of E-292, confirming the +0.1 slack cause. (2) Done: power line reworded; no q_S-independent NO control exists. Edited by the chair.

# F-055 replicates with independent code: the Phi-equivariant Z-congruence search says NO on 19 of 19 q_S-separated LNA pairs and YES (exact P) on 19 of 19 within-quipu controls at n = 9..12, at slack 0.5 and 5

author: theorist · round: 139 · kind: result (replication)
thread: T3/T8 · bears on: F-055, E-269, E-292, E-276, H-024
scope: n = 9..12; exactly the 19 cross pairs of E-292 (copied from `rounds/110/skeptic_qsno*.out`, group-of-140 row excluded) and 19 within-quipu controls rebuilt by E-292's rule (first two members of a's quipu in `lnaKeys` order, library `quipuForAlmostSeparateLNA`); slack 0.5 and 5; plus one planted U C_a U^T per group. Not covered: lam < 1.5, other pairs per group, E-314's planted n = 11 group, the group of 140, interval arithmetic (float majorant as before).
power: YES side only: planted pairs (b = U C_a U^T) are found 19/19 at both slacks. NO side: on all 19 cross pairs NO is true by my own q_S certificate and the search says NO 19/19, but there is no NO control independent of q_S in this submission (the search's NO is not shown able to say NO where q_S is equal, apart from the group-of-140 extra, whose truth is unknown). Planted YES (b = U C_a U^T, U a product of 4 random elementary matrices) found 19/19 at both slacks; enumerator checked against brute force on 5 random 5-dim Grams (372/356/456/414/748 vectors, all equal).

## Pre-registration (written into this file before any run)

1. YES (exact R, R^T C_a R = C_b in integers) on a cross pair: q_S cannot differ, so my q_S code, E-292's certificate, or the pair construction is wrong; F-055's NO for that row is void (the exact R outranks every NO).
2. NO on a control: my enumeration incomplete or the control is not E-292's; rerun at slack 5, check planted YES; if still NO, report as disagreement, unresolved.
3. NO on a cross pair with my q_S equal: not a disagreement, but that NO rests on the search alone; reported separately.
4. A crash or domain assertion: "not decided", not agreement.

Outcome: none of 1-4 occurred. (A crash of the PLANTED control at one n = 12 group in the first run -- sympy LLL internal assert, then a non-PD float Gram for a badly conditioned conjugate -- was fixed by an LLL retry and by retrying U; `failedU = 0` on the final runs. It never touched a cross or control pair.)

## Claim

For each of the 19 cross pairs (a, b) of E-292 there is no R in GL_n(Z) with R Phi_b = Phi_a R and R^T C_a R = C_b (my search: NO), and for each of the 19 controls (a, a2) such an R exists (exact integer witness found). This is F-055's 19/19 NO and 19/19 YES, reproduced with code that shares nothing with `maverick_iso140.py` except the mathematics, and without pari. It does not extend F-055: same pairs, same domain (signature (n-1,1), one simple Salem pair), same float majorant. Refuted by: an exact R on a cross pair, or a control with no R found at any slack.

Bound, re-derived. S = C + C^T, Phi = -C^-T C preserves S. Let p, q be Phi-eigenvectors for lam, 1/lam (isotropic), scaled so S(p, q) = 1; W = <p, q>^perp is S-positive definite (signature (n-1,1)). For x = alpha p + beta q + w: S(x) = 2 alpha beta + S(w); H(x) = alpha^2 + beta^2 + S(w) is positive definite, H = S + (S(q-p))(S(q-p))^T, and {p, q, S-ON basis of W} is H-orthonormal. An intertwining isometry R sends p_b -> s p_a, q_b -> q_a / s, W_b -> W_a isometrically, so Q(R) = tr(H_b^-1 R^T H_a R) = s^2 + s^-2 + n - 2. Phi_a^k R is again a solution with s multiplied by lam^k, so some solution has |s| in [lam^-1/2, lam^1/2] and Q <= n - 2 + lam + 1/lam. This is E-276's closed form; I got it without reading E-269's derivation first, but it is the canonical argument, so "independent" means independent code, not a new idea.

## Evidence

Own pipeline: Cartan from relation lengths (equals library `lnaCartanMatrix` on all 57 algebras, `libC=True`); commutant lattice = sympy rational nullspace, saturated via Smith form (row_i(U B)/d_i), exact-integer LLL (sympy DomainMatrix); float Gram of Q, own Gram-LLL, Gram recomputed from the exact new basis, own Fincke-Pohst; each candidate checked exactly (intertwining asserted, then R^T C_a R == C_b; det C = 1 forces det R = +-1). q_S: invariant factors of S and histogram of x^T S x / 2 mod 1 over S^-1 Z^n / Z^n, enumerated through the Smith form.

| n | lam | abs det S | a | b | q_S differ | cross nvec 0.5/5 | NO | ctrl a2 | ctrl nvec/nsol 0.5 | 5 | YES |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 9 | 1.7816 | 12 | 2223030 | 2400230 | yes | 1/1 | NO | 2203030 | 5/4 | 9/8 | YES |
| 10 | 2.0306 | 28 | 22240030 | 22303000 | yes | 2/8 | NO | 22040030 | 9/4 | 11/6 | YES |
| 10 | 1.8098 | 16 | 22223030 | 22302300 | yes | 3/5 | NO | 22203030 | 7/4 | 13/6 | YES |
| 11 | 2.3193 | 52 | 222303030 | 230400230 | yes | 1/2 | NO | 220303030 | 5/4 | 6/4 | YES |
| 11 | 2.1257 | 32 | 223023030 | 302240030 | yes | 46/156 | NO | 223003030 | 37/2 | 129/6 | YES |
| 11 | 1.8238 | 20 | 222223030 | 224002230 | yes | 1/2 | NO | 222203030 | 3/2 | 8/6 | YES |
| 11 | 2.0062 | 32 | 222230300 | 224002300 | yes | 0/1 | NO | 222230302 | 1/1 | 4/3 | YES |
| 12 | 2.4596 | 84 | 2230303000 | 2400400230 | yes | 4/11 | NO | 2230303020 | 6/2 | 14/6 | YES |
| 12 | 2.3286 | 64 | 2222303030 | 2302500030 | yes | 7/11 | NO | 2220303030 | 11/4 | 15/4 | YES |
| 12 | 2.3575 | 60 | 2302304000 | 2400240030 | yes | 6/18 | NO | 2302304002 | 8/2 | 25/6 | YES |
| 12 | 2.2671 | 52 | 2230302300 | 2304002230 | yes | 4/8 | NO | 2230302302 | 5/1 | 9/3 | YES |
| 12 | 2.2489 | 48 | 2230230300 | 2240023030 | yes | 5/7 | NO | 2230230302 | 9/1 | 20/3 | YES |
| 12 | 2.0810 | 48 | 2222303000 | 2250002300 | yes | 3/6 | NO | 2222303020 | 4/1 | 9/3 | YES |
| 12 | 2.0149 | 39 | 2222230300 | 2223024000 | yes | 0/3 | NO | 2222230302 | 1/1 | 3/3 | YES |
| 12 | 2.0582 | 44 | 2222240030 | 2230240000 | yes | 3/3 | NO | 2222040030 | 7/2 | 12/6 | YES |
| 12 | 1.9165 | 32 | 2223023000 | 2250002230 | yes | 5/7 | NO | 2223023020 | 9/1 | 24/3 | YES |
| 12 | 1.8688 | 28 | 2222302300 | 2224002230 | yes | 3/3 | NO | 2222302302 | 4/1 | 12/3 | YES |
| 12 | 1.7884 | 20 | 2230223000 | 2240022230 | yes | 4/8 | NO | 2230223020 | 7/3 | 21/7 | YES |
| 12 | 1.6610 | 12 | 2230222300 | 2400222230 | yes | 3/6 | NO | 2230222302 | 4/1 | 12/3 | YES |

Agreement with E-292's per-row counts (skeptic's table, round 110): verdicts 38/38 equal. Vector counts equal on 73 of 76 (count, slack) cells (referee's correction; I had said 74); the three differences are the n = 9 cross pair at slack 5 (1 here, 5 in E-292), and n = 11 (lam 1.8238) control at slack 0.5 (mine 3/2, theirs 5/4) and n = 12 (lam 2.0810) control at slack 5 (mine 9, theirs 10). Cause, checked on the first: `maverick_iso140.py` calls qfminim with `int(BOUND*SC) + 10**5` at SC = 10^6, i.e. an extra 0.1 of slack; my norms there are 4.0, 11.08 (x2), 11.90 (x2) against bound + 0.5 = 11.87, so the 11.90 pair lies in their extra 0.1. Not a verdict difference; worth one line in E-292 ("slack is effectively +0.1 in E-269's code").

Weakest step: the majorant Gram is floating point (as in E-269), so a vector with Q within ~1e-9 of bound + slack could be lost; slack 5 makes this immaterial for a YES, and for the NO side the q_S certificate carries the truth independently. The 5 weak NOs (nvec 0-1 at slack 0.5) are as weak here as in E-292.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/139/theorist_zcong.py 0.5   # 2 min 16 s (80 s is lnaKeys(12))
timeout 10m .venv/bin/python -W ignore workshop/rounds/139/theorist_zcong.py 5     # 2 min 21 s
```
Group of 140 (6 s): `timeout 10m .venv/bin/python -W ignore -c "import sys; sys.path.insert(0,'workshop/rounds/139'); import theorist_zcong as z; n=12; A=z.cartan(n,(2,2,2,3,0,2,3,0,3,0)); B=z.cartan(n,(2,2,3,0,3,0,2,2,3,0)); print([(s, z.zcong(n,A,B,s)['cong']) for s in (0.5,5,30)])"`
Outputs: `workshop/rounds/139/theorist_zcong_s05.out`, `theorist_zcong_s5.out` (< 10 KB each).

## Prior record

F-055 / E-292 (the claim replicated), E-269 (method), E-276 (closed-form bound, rederived here), E-281, E-314 (planted group, not rerun). Nothing new beyond the replication and the +0.1 slack note; the result is known and now confirmed independently. Not a rediscovery in the bad sense: the assignment was a rebuild.

## Code changed

None under `quivermutation/` or `tests/`. New: `workshop/rounds/139/theorist_zcong.py` and its two `.out` files. No tests needed (no library file touched).

## Next

- F-055's standing gap is unchanged: planted-pair power at a second group (skeptic's blocking check) and lam < 1.5; my script takes arbitrary integer C (det 1), so a constructed non-LNA pair can be dropped in as `zcong(n, Ca, Cb)`.
- Done as an extra (outside the 38): E-269's own pair, n = 12 group of 140 (a = 2223023030, b = 2230302230), where q_S is EQUAL (my code agrees) and the search is the only certificate: NO at slack 0.5 (12 vectors, = E-269), 5 (21) and 30 (485), 0 isometries. So E-269's pure-search NO also replicates with independent code. Next for it: an algebraic certificate (Z[Phi]-module / ideal class) that does not depend on a float majorant -- scholar or maverick with cypari restored.
