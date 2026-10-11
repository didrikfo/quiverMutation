# On LNA walks Q(x) = x(adj S_ii - adj S_wi - adj S_iw); its x^2 coefficient is 1 exactly when two Coxeter-orbit numbers of the simple S_w vanish, which is proved for F e_w = e_i and only observed otherwise

author: theorist · round: 042 · kind: proof (partial) + result
thread: T5 · bears on: E-143, E-140, E-142, E-138, H-015

## Claim

Setting (E-138, as in 041). A has Cartan matrix C_A (acyclic, so unitriangular after a permutation), v is a gate-admitted vertex, J = e_i (every J != 0 step seen has dim 1, one vertex), C' = r C_A r^T, C_B = C' + E_{vi}. Order v last and put Z = C_A restricted to V \ v, Y = Z^-1, T = Y^T, F = Z Z^-T (minus the Coxeter matrix of Z), S(x) = xZ + Z^T, N = S^-1.

**(P1, proved, exact).** If (H1) C' e_v = e_v (column v of C' is the unit vector) and (H2) the row of v in C' is u = e_w - e_i off the diagonal, then C' = [[Z,0],[u^T,1]], C_B = [[Z,0],[e_w^T,1]], and by the Schur complement det(xC+C^T) = (1+x) det S - x u^T adj(S) u. Hence
 Q(x) = x B(x),  B = adj(S)_ii - adj(S)_wi - adj(S)_iw = det S * (N_ii - N_wi - N_iw)
(checked against the library's C_B: 716/716 n = 6, 120/120 n = 7 steps in scope). Writing N = sum_k (-x)^k T F^k and c_k := (T F^k)_ww = chi(S_w, F^k S_w): c_0 = 1 and c_{-k} = c_{k-1} (Serre, from F^T T F = T and T F^-1 = Y; checked), and **Q_1 = B_0, Q_2 = B_1** where, if F^s e_w = e_i, B_0 = c_0 - c_s - c_{-s}, B_1 = -(c_1 - c_{1+s} - c_{1-s}) + sigma_1 B_0 (sigma_1 = tr(T Z)). So "lowest term x^2 with coefficient 1" is the statement B_0 = 0, B_1 = 1: a statement about the first two moments of the F-orbit of the simple S_w.

**(P2, proved).** If in addition F e_w = e_i (s = 1; "Phi S_w = S_i"), then Y_wi = Y_ii = 1 forces w before i in the triangular order, so c_1 = T_wi = Y_iw = 0, B_0 = 0 (Q_1 = 0 for free), and B_1 = 1 + c_2 - c_1 = 1 + c_2 with c_2 = chi(S_w, F S_i) = sum_b Y_bw (F e_i)_b. So: **Q = x^2 + O(x^3) iff c_2 = 0**, and c_2 = 0 holds whenever F e_i = -e_m (then i before m since Y_im = -1, so Y_mw = 0), or more generally when F e_i is supported off the vertices that precede w.
Proof sketch of every step: F e_w = e_i <=> Y^T e_w = Y e_i, i.e. Y_wb = Y_bi for all b; b = i gives Y_wi = 1. Terms: W = T Z T = T F, W_iw = T_ii = 1, W_wi = (T F e_i)_w = c_2. Weakest step: none in P2 itself; the weak points are the hypotheses, below.

**Not proved, only computed (the exact gaps).**
(G1) H1 and H2 are observations about the library's walk steps, not consequences of the gate. H2 says C_B has row (v) = e_v + e_w: after the mutation the only new arrow at v is the old one. Counter-shape: in the n = 4 example of E-143 u = e_2 + e_4 (not e_w - e_i, u_i = 0) and Q = 0 identically, so the reduction does not apply there, consistently with E-143.
(G2) The orbit relation e_i = F^s e_w (s = 1 in 176/716 n = 6 steps and 120/120 n = 7 steps; s = -2 in the other 540 of 716 at n = 6; never absent when H1, H2 hold). For s = -2 the same moment formulas apply (B_0 = 1 - c_{-2} - c_2 = 1 - Y_iw - Y_wi, B_1 = c_{-1} + c_3 - c_1 = 1 + c_3 - c_1, verified numerically) but I have no argument for the numbers c_1 = 0, c_2 = 1, c_3 = 0 that give B_0 = c_0 - c_1 - c_2 = 0 and B_1 = 1; they are observed. Why e_i lies in the F-orbit of e_w at all is the open question; for s = 1 it reads "the relation i ~> v -> w = 0 is a relation between S_i and S_w in eAe, and Phi S_w = S_i".
(G3) c_2 = 0 in the s = 1 cases (176 n = 6, 120 n = 7) is observed; the termwise vanishing sum_b Y_bw (F e_i)_b = 0 holds in 130/176 (n = 6) and 120/120 (n = 7); in the other 46 F e_i has full support and the sum cancels.
(G4) Random unitriangular Z with Y_wi = 1 (so B_0 = 0) have B_1 anything (-37..13): the identity B_1 = 1 is not a matrix identity of Z; it needs F-orbit data that only realisable (Cartan) Z have. This is the "identity that fails to generalise".

## Evidence

| set | J != 0 steps | H1 & H2 hold | s = 1 | s = -2 | Q low term | exceptions to H1/H2 |
|---|---|---|---|---|---|---|
| n = 6 c0 guarded walk, 522 s | 766 (314 distinct (Z,w,i)) | 716 | 176 | 540 | x^2, coeff 1: 766/766 | 50: u = e_a + e_b - e_i (42) or H1 fails (8); Q still x^2 (1 + ...) |
| n = 6 c0 key guard OFF, depth 10, parent has class key | 82 | 82 | 38 | 44 | 82/82 | 0 |
| n = 7 c0 guarded walk, 521 s | 120 | 120 | 120 (s in {1, -4}) | 0 | 120/120 | 0 |

C_B = C' + H (E-138) holds on every gate-admitted record: 74 384/74 384 (n = 6), 44 761/44 761 (n = 7); this closes the round-041 review item for n = 7. Every J != 0 step has J = e_i (dim 1). Moment sequences are few: s = 1, n = 6: c_{-3..6} = 0 0 1 1 0 0 -1 -1 0 0; s = 1, n = 7: 0 0 1 1 0 0 0 1 1 0; s = -2, n = 6: 1 0 1 1 0 1 0 0 1 0. Newton form of the result: tr F and the x^1 coefficient agree for A and B; the e_2 coefficient of det(xC+C^T) differs by exactly 1, equivalently tr(Phi_B^2) = tr(Phi_A^2) - 2.
The n = 6 counts here (766) exceed E-142's 147 because the budget is 522 s and different walk order; the claim is about this sample, not a re-count of E-142.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/041/theorist_dump.py 6 0 520 /tmp/n6.pkl      # 522 s (pkl not kept)
timeout 10m .venv/bin/python workshop/rounds/041/theorist_dump.py 7 0 520 /tmp/n7.pkl      # 521 s
timeout 10m .venv/bin/python workshop/rounds/042/theorist_dump_off.py 6 0 520 /tmp/n6off.pkl 10   # 185 s, guard off
.venv/bin/python workshop/rounds/042/theorist_reduce.py /tmp/n6.pkl        # P1 formula, u shape, orbit s, moments; ~2 min
.venv/bin/python workshop/rounds/042/theorist_reduce.py /tmp/n6off.pkl onlyparentkey
.venv/bin/python workshop/rounds/042/theorist_lemma.py /tmp/n6.pkl         # P2 hypotheses, seconds
.venv/bin/python workshop/rounds/042/theorist_moments.py /tmp/n7.pkl       # c_k sequences and Serre symmetry
.venv/bin/python workshop/rounds/042/theorist_orbits.py /tmp/n6.pkl        # all F-orbit coincidences
.venv/bin/python workshop/rounds/042/theorist_exceptions.py /tmp/n6.pkl    # the 50 off-shape steps
.venv/bin/python workshop/rounds/042/theorist_cbcheck.py /tmp/n7.pkl       # C_B = C' + H
.venv/bin/python workshop/rounds/042/theorist_n4.py                        # E-143 n = 4 example: u_i = 0, Q = 0
.venv/bin/python workshop/rounds/042/theorist_zrandom.py 5 3000            # prints a LOT (head it); G4
```
(`theorist_explore.py`, `theorist_blocks.py` are exploratory.) The n = 6 sample is a time-limited walk prefix, so counts vary a little with machine speed.

## Prior record

E-143 states Q's lowest term x^2 (n = 6, 7) as observed, "no proof; x = 0, infinity, trace routes closed". New here: the exact block form Q = x B (P1), the reduction of Q_2 = 1 to Euler-form moments of the simple S_w under F, the proof for F e_w = e_i (P2), the orbit relation e_i = F^s e_w (s = 1 or -2) as the new regularity, the exceptions (50) which also have x^2 lowest term, and the n = 7 check C_B = C' + H. Grep of EXPERIMENTS/FINDINGS/RETRACTIONS/HYPOTHESES for "Serre", "moment", "one-point": nothing relevant; no retraction touches these. Bears on the key guard: R = 1 would need B = 0, which requires u not of the form e_w - e_i (as in the n = 4 example); so a proof of the walk law would need G1 (why the new vertex v of B has exactly one arrow) plus G2.

## Code changed

None in the library. New scripts in `workshop/rounds/042/`: `theorist_dump_off.py` (041 dump with guard off and depth cap), `theorist_reduce.py`, `theorist_lemma.py`, `theorist_moments.py`, `theorist_orbits.py`, `theorist_exceptions.py`, `theorist_cbcheck.py`, `theorist_n4.py`, `theorist_zrandom.py`, `theorist_explore.py`, `theorist_blocks.py`. No tests run (no library file touched).

## Next

- Theorist/skeptic: derive H2 and the orbit relation from the algebra. Conjecture to test: w, i are the endpoints of the single minimal relation r: i ~> v -> w = 0 (J_i = 1), and eAe (Cartan Z) has Ext^2(S_i, S_w) = k with S_i = Phi S_w up to sign; then c_1, c_2 follow from Serre duality in D^b(eAe). Needs a module-level (not Cartan-level) computation: build eAe and its simples' AR data at n = 6 for the 12 s = 1 shapes.
- Skeptic: independent derivation of P1 from C_B = C' (I + E_vi) (transvection), and a J != 0 step with H1 and H2 but c_2 != 0 (that would be a key-preserving candidate: Q_2 = 1 + c_2 - c_1 = 0 needs c_2 = -1).
- Experimentalist: the same tally for n = 8 (a J != 0 record needs H1, H2; report s and c_2); and any step where u has the form e_w - e_i but the orbit relation is absent (none at n = 6, 7).
- Search idea (cheap, matrix level): among all unitriangular Z with H1, H2, find F-orbit data making B = 0 (key-preserving); if such Z always violates a realisability constraint, that constraint is the theorem.
