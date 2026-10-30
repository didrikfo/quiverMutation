# The n = 12 key group of 140: the certified-inequivalent pair is NOT Z-congruent as Cartan matrices (9 of 9 cross pairs; 6 of 6 within-quipu pairs are), by an exhaustive finite search in the Phi-commutant lattice

author: maverick · round: 099 · kind: result
thread: T3/T8 · bears on: E-265, E-259, E-267 (and the "Cartan-level NO" use of q_S)
scope: n = 12, key group of 140 only (36 + 40 members, two quipus); the open pair reps (2,2,2,3,0,2,3,0,3,0) vs (2,2,3,0,3,0,2,2,3,0) plus up to 3 members per quipu (15 pairs). Floating-point majorant used to size the search (not a proof-grade bound, see Evidence); every accepted/rejected R is then exact integer arithmetic. Not covered: other groups, Z[Phi]-module iso alone (only the stronger isometry question), n = 13.

## Hit / miss criterion (stated first)
HIT = an integer P with det P = +-1 and P C_a P^T = C_b is found (explicit certificate), or an exhaustive search proves there is none. MISS = neither (search incomplete). Control required: the same code must find P for pairs known congruent.

## Claim

C_a and C_b of the open pair are not Z-congruent (nor is C_a congruent to C_b^T). So the Cartan matrices separate it, and E-259's list of Cartan-indistinguishable certified pairs is now empty at n = 11, 12 (q_S: 7 of 8, this: the 8th). Same result for all 9 cross pairs among 3+3 members of the 140 group. It does NOT say the S = C + C^T lattices are non-isometric: S_a, S_b are in one genus and the indefinite-class-number argument (not verified here) would make them Z-congruent; the obstruction is the Phi-equivariance. Refuted by: one integer P with P C_a P^T = C_b, which the code would also print.

## Evidence

Reduction. C = S (1 - Phi^-1)^-1 with S = C + C^T and det(1 - Phi) = -40 != 0 (E-265), and Phi = -C^-T C. Hence for R = P^T: R^T C_a R = C_b iff R Phi_b = Phi_a R and R^T S_a R = S_b. R is an integer matrix of the commutant lattice Lam = {R : R Phi_b = Phi_a R}, computed by `matkerint`: rank 14 (= 2^2 + 2 + 8, Phi semisimple at -1, char poly (x+1)^2 (x^2+1) g_8).
Finiteness. S has signature (11, 1) for both. g_8 has exactly one real pair (lam = 2.1451, 1/lam, a Salem number) and three unit-circle pairs; the hyperbolic plane of S is the Salem eigenplane, everything else is positive definite. So the S-isometries commuting with Phi are a finite set times Phi^Z. Majorant H = S + (S d)(S d)^T, d = q - p (p, q eigenvectors for lam, 1/lam, S(p,q) = 1). For an isometry R, tr(H_b^-1 R^T H_a R) = tr(H_b^-1 S_b) + z^T S_b H_b^-1 S_b z with z = t q_b - p_b/t, and R -> Phi_a^k R moves t into [lam^-1/2, lam^1/2]. This gives norm <= 13.11 (10 from the compact part, the rest from the hyperbolic plane). `qfminim` lists ALL vectors of Lam up to that bound; each is tested exactly (R^T S_a R == S_b, det, P C_a P^T == C_b).
Results (bound 13.11 and again 42.6 = slack 30; `maverick_iso140_open.out`): open pair 12 resp. 487 lattice vectors up to sign, 0 exact isometries. Table over the 140 group (`maverick_iso140_table.out`, 15 pairs):

| pair type | pairs | congruent found | lattice vectors <= bound |
|---|---|---|---|
| within quipu (control, derived-equivalent) | 6 | 6 (2 or 4 solutions each: +-P, and automorphisms) | 14 or 16 |
| across the two quipus (certified inequivalent) | 9 | 0 | 12 |

Controls: (i) a = a gives P = I (and the Phi-multiples inside the bound); (ii) C_b := Q C_a Q^T for a random product of 25 elementary matrices Q (`maverick_ctl.py`): found, 4 solutions at 13.1, 8 at 42.6; (iii) the 6 within-quipu pairs. An earlier version had the commutant equation transposed and gave false NOs on the controls; fixed (R Phi_b - Phi_a R, row-major vec = I (x) Phi_b^T - Phi_a (x) I), controls then pass. Opposite-algebra convention: C_a vs C_b^T also NO (`maverick_transp.py`; control a vs a^T found).
Least sure step: the bound. The derivation (window of t of width log lam, norm identity above) is checked by the controls (true solutions of norm exactly inside the bound) and by the NO persisting at 3x the bound, but the float majorant is not an interval-arithmetic proof; the qfminim Gram was rounded at 1e-6. A referee should re-derive the identity tr(H_b^-1 R^T H_a R) and perhaps rerun with slack 100 (the first slack-100 run overflowed my 150000-vector cap, not a result).
Why it fits: the pair has q_S equal (E-265) so genus-level obstructions vanish; what separates it is the Phi-compatibility, a conjugacy-type condition on Phi in O(S), consistent with E-265's note that C-congruence = S-isometry + Phi-equivariance.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/099/maverick_iso140.py 30        # open pair + control, ~5 s
timeout 10m .venv/bin/python workshop/rounds/099/maverick_iso140.py table 3   # 15 pairs of the 140 group, 85 s
timeout 10m .venv/bin/python workshop/rounds/099/maverick_ctl.py              # random-P control, ~5 s
timeout 10m .venv/bin/python workshop/rounds/099/maverick_transp.py           # transposed convention, ~5 s
timeout 10m .venv/bin/python workshop/rounds/099/maverick_setup.py            # signature (11,1), Phi spectrum
```
(cypari, python-flint were missing in the venv; installed with uv pip.)

## Prior record

E-265 / E-267: this pair open; E-267: Nikulin 1.13.3 unverified. E-215: definite-stage isometry test P C_a P^T = C_b (n = 10, 11 definite); E-232, E-242: YES at n = 11, 12 periodic groups. Grep of `research/` for "Salem" and "commutant" in this role found no recorded decision of this pair. New: first complete (not random) decision procedure, and a Cartan-level NO with certified-inequivalent pair that q_S could not see. Does not touch derived equivalence for infinite global dimension.

## Code changed

None in the package. New scripts in `workshop/rounds/099/`: maverick_setup.py, maverick_iso140.py, maverick_ctl.py, maverick_transp.py (+ two .out files). No tests run (no package files touched).

## Next

Skeptic: re-derive the bound and rerun with interval/larger slack; try to make a case where it wrongly rejects. Toolsmith: fold `test` into a routine for all key groups with one Salem pair (signature (n-1,1)) at n = 11..13, to see which q_S-equal pairs are separated. Chair: E-265 can be amended (eight of eight). Scholar: why Phi-equivariance, not genus, is the obstruction (spinor-genus / hermitian class of the K8-piece).
