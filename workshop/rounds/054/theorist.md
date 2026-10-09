# Generation of K^b(proj A) by the one-step complex T is automatic for every loopless step (so for all tiltingPlus steps on acyclic algebras), by a two-line cone argument; what is not automatic is "End(T) is the next algebra"

author: theorist · round: 054 · kind: proof
thread: T10 · bears on: E-155, E-158, E-161, E-164, E-165, E-159, E-128, H-015
scope: proof for any finite-dimensional A (quiver with relations) and any vertex v with no loop, one step at a time; hypotheses checked by script on all LNAs and duals n = 5, 6, 7 (2116 steps, every vertex with an out-arrow) and on the 13 E-161 edges + the 16 failing J != 0 steps (n = 7, class 1). Says nothing about End(T) = next algebra beyond E-165, nor about steps with loops or non-acyclic quivers in the code.

## Claim

Let A be the current algebra, v a vertex without a loop, T = (+_{i != v} P_i) + T_v with T_v = cone(f: P_v -> M), M = +_{b: v->h} P_h (all out-arrows, the map given by the arrows). Then thick(T) = K^b(proj A), whatever Hom(T,T[m]) is. Proof: no loop gives h != v for every arrow, so P_h is a summand of T for each h and M in add(T); the triangle P_v -> M -> T_v -> P_v[1] puts P_v in thick(M, T_v) in thick(T); the remaining P_i, i != v, are summands of T; K^b(proj A) = thick(A). The weakest step is none; it uses only that the arrows of v are all in f (they are: `arrowsOutOf`) and that v has no loop (acyclic quivers, which all walks and LNAs are).
Consequences: (1) every tiltingPlus step is a tilting complex as soon as Hom(T,T[m]) = 0 for m != 0, and the "generation assumed (Okuyama-Rickard)" caveat in E-159/E-164/E-165 can be dropped for the single step: it is proved, not assumed. (2) m = +1 never fails (data: Hom(T,T[1]) = 0 in all 2116 + 29 steps); J != 0 means Hom(T,T[-1]) != 0, i.e. T is silting but not tilting (E-128, AI 2.31 cited there; this proof is independent of AI). The 16 failing steps have dim Hom(T,T[-1]) = 1 (14 steps) or 2 (2 steps), all with generation. (3) Hence a failing step is still a genuine silting mutation but not a derived equivalence; the premise "J = 0" is exactly the vanishing of Hom(T,T[-1]), nothing about generation.
Not automatic, and not covered: (a) End(T) isomorphic to the algebra the code builds as `reducePathAlgebra(quiverMutationAtVertex)` (E-165: quiver level, label-preserving, only the 13 E-161 edges, class 1; the 3 E-158 paths and E-155 paths are Cartan level, E-164); (b) the iteration: a path A_0 -> ... -> A_N is a derived equivalence only if each End(T_k) is A_{k+1} as an algebra, not just in Cartan matrix; generation is then applied to A_k, a new algebra, step by step (fine, the proof is per step); (c) R steps are F steps on the dual algebra, so they need the same hypothesis on the dual (checked in path13: 7 R-steps, no loops). A counting argument (n summands, K0 class matrix of det +-1) is also satisfied but is vacuous evidence: det = +-1 holds by construction ([T_v] = sum [P_h] - [P_v]), so it does not test anything beyond the proof. Finite global dimension plays no role. Would be refuted by: a step whose quiver has a loop at v (the code would then not be tiltingPlus as described).

## Evidence

Script `theorist_gen.py` checks the proof's hypotheses and records the Hom signs (own Hom-complex code from rounds/050 skeptic_tilt.py, independent of tiltingPlus/perI):

| set | steps | acyclic, targets != v, T has n summands | det of class matrix | Hom(T,T[1]) = 0 | Hom(T,T[-1]) = 0 |
|---|---|---|---|---|---|
| LNAs + duals, n = 5 | 112 | 112 | 112 of +-1 | 112 | 70 |
| n = 6 | 420 | 420 | 420 | 420 | 252 |
| n = 7 | 1584 | 1584 | 1584 | 1584 | 924 |
| E-161 path, 13 edges | 13 | 13 | 13 | 13 | 13 |
| 16 failing steps (class 1, n = 7) | 16 | 16 | 16 | 16 | 0 |

(Dual algebras included in the LNA rows; every vertex with at least one out-arrow, not only gate-admitted vertices, so the J != 0 rows there are not all gate-admitted steps.) Failing-step Hom(T,T[-1]) dimensions: 1 x14, 2 x2. The failing set is rebuilt from the E-149 walk (20 000 expansions, 519 s) and agrees with E-165's 16.

## Reproduction

```
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 5      # 2 s
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 6      # 3 s
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 7      # 8 s
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # 519 s (mkdir -p /tmp/tsm first)
.venv/bin/python workshop/rounds/054/theorist_gen.py path13 > workshop/rounds/054/theorist_gen_path13_out.txt   # 2 s
```

## Prior record

E-128/E-126/E-136 use "AI 2.31: mutation is silting" (cited from memory; the literature summaries are marked UNVERIFIED in STATE T5) and E-159/E-164/E-165 say "generation assumed (Okuyama-Rickard)". So the content "T is silting, hence generates" is implicit in the record, but no entry gives a proof or a check of the hypotheses; this submission gives the elementary proof that avoids AI 2.31 and the check. It is therefore a closing of a caveat, not a new phenomenon. Not in RETRACTIONS.

## Code changed

None (new script `workshop/rounds/054/theorist_gen.py`, output `theorist_gen_path13_out.txt`, 2 KB). No tests touched.

## Next

- toolsmith/skeptic: the real residual gap of T10 is (a), End(T) = next algebra as algebras on the E-155/E-158 paths, class 2, and the 8 parallel-arrow failing steps (needs matrix-valued arrows); generation can now be struck from their caveat lists.
- skeptic: check that the code's mutation never meets a loop (acyclicity of every intermediate child on the walk; the script asserts it for the 29 edges only).
- Reword STATE T10(b)/(c): "Generation still assumed" -> "generation proved per step; premise = Hom(T,T[-1]) = 0 plus End(T) iso".
