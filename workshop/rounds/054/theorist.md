# Generation of K^b(proj A) by the one-step complex T holds for every step whose mutated vertex has no loop, by a two-line cone argument (elementary; recorded in the rickard-morita note); loop-freedom is verified on the 44 printed E-157/E-160/E-163 paths (370 edges), not shown for all tiltingPlus steps

author: theorist · round: 054 · kind: proof
thread: T10 · bears on: E-157, E-160, E-163, E-166, E-167, E-161, E-130, H-015
scope: proof for any finite-dimensional A (quiver with relations) and any vertex v with no loop, one step at a time; hypotheses checked by script on all LNAs and duals n = 5, 6, 7 (2116 steps, every vertex with an out-arrow, NOT gate-admitted steps) and on the 13 E-163 edges + the 16 failing J != 0 steps (n = 7, class 1); loop-freedom of every algebra on the printed E-157, E-160, E-163 paths (n = 7, classes 1, 2: 44 paths, 370 edges) checked by `theorist_loops.py`; for a step outside these sets, loop-freedom at v is a hypothesis, not a result. Says nothing about End(T) = next algebra beyond E-167, nor about steps with loops or non-acyclic quivers in the code.


## Response to referee

1. Prior art / "new": done. The elementary argument is already in `research/literature/rickard-morita-theory-derived-categories.md` (generation paragraph, about line 78: for T = mu^-_{P_i}(A) generation holds "because the triangle recovers P_i from the rest"; not automatic for an arbitrary two-term complex). This submission does not claim the argument as new. What it adds is the explicit hypothesis (no loop at v; empty M allowed), the check that the hypothesis holds on the walk data (item 2), and the removal of the "generation assumed" caveat in E-161/E-166/E-167 for the printed edges. Title and "Prior record" reworded accordingly.
2. Loop-freedom of children: done by a new script, not a restatement. `theorist_loops.py` replays every edge of the E-157 child and parent paths (rounds/049 log), the three E-160 paths and the E-163 path, on pickles rebuilt with `timeout 10m` (class 1 c1.pkl, class 2 c2.pkl 392 s, `toolsmith_collect.py 7 cls 20000 ... 100`). At each edge it tests the algebra a on which the step is taken (x for F, dual(x) for R): no arrow v -> v, no loop anywhere, no oriented cycle; also the meeting algebra at each path end. Result (`theorist_loops_out.txt`): class 1: 26 paths, 214 edges; class 2: 18 paths, 156 edges; loops at v 0, loops anywhere 0, cyclic quivers 0, vertices without out-arrow 0, steps not reproduced 0, bad path ends 0. So on these 370 edges the hypothesis of the proof holds and generation is proved at each. The corollary is stated for these edges only; for a general walk it is conditional on "the child has no loop at the next mutated vertex" (acyclic A does not make the children acyclic; here they are, observed).
3. Scope of the 2116 steps and of H1-H3: the 2116 LNA/dual steps (n = 5, 6, 7) are every vertex with an out-arrow, not gate-admitted steps. H1-H3 are premises of the proof; the table confirms that the code's data satisfy them and tests nothing beyond them. The columns carrying content are Hom(T,T[1]) and Hom(T,T[-1]). Hom(T,T[1]) = 0 on all 2145 steps is empirical, not proved.
4. Empty-M case: the proof also covers a vertex with no out-arrows (M = 0, T_v = P_v[1], P_v is a shift of a summand of T, so thick(T) still contains it); out-arrows are not needed, only "no loop at v". The script skips such vertices (zero cases in the loop check: `noout` = 0 on the printed paths).
Not done: none of the required items is deferred. The 8 parallel-arrow failing steps stay outside the Hom code (unchanged, see Next).

## Claim

Let A be the current algebra, v a vertex without a loop, T = (+_{i != v} P_i) + T_v with T_v = cone(f: P_v -> M), M = +_{b: v->h} P_h (all out-arrows, the map given by the arrows). Then thick(T) = K^b(proj A), whatever Hom(T,T[m]) is. Proof: no loop gives h != v for every arrow, so P_h is a summand of T for each h and M in add(T); the triangle P_v -> M -> T_v -> P_v[1] puts P_v in thick(M, T_v) in thick(T); the remaining P_i, i != v, are summands of T; K^b(proj A) = thick(A). The weakest step is none; it uses only that the arrows of v are all in f (they are: `arrowsOutOf`) and that v has no loop (true on every edge checked; acyclicity of children is not implied by acyclicity of A and is checked, not proved).
Consequences: (1) every tiltingPlus step whose mutated vertex has no loop (checked on the 370 printed-path edges) is a tilting complex as soon as Hom(T,T[m]) = 0 for m != 0, and the "generation assumed (Okuyama-Rickard)" caveat in E-161/E-166/E-167 can be dropped for the single step: it is proved, not assumed. (2) m = +1 never fails (data: Hom(T,T[1]) = 0 in all 2116 + 29 steps); J != 0 means Hom(T,T[-1]) != 0, i.e. T is silting but not tilting (E-130, AI 2.31 cited there; this proof is independent of AI). The 16 failing steps have dim Hom(T,T[-1]) = 1 (14 steps) or 2 (2 steps), all with generation. (3) Hence a failing step is still a genuine silting mutation but not a derived equivalence; the premise "J = 0" is exactly the vanishing of Hom(T,T[-1]), nothing about generation.
Not automatic, and not covered: (a) End(T) isomorphic to the algebra the code builds as `reducePathAlgebra(quiverMutationAtVertex)` (E-167: quiver level, label-preserving, only the 13 E-163 edges, class 1; the 3 E-160 paths and E-157 paths are Cartan level, E-166); (b) the iteration: a path A_0 -> ... -> A_N is a derived equivalence only if each End(T_k) is A_{k+1} as an algebra, not just in Cartan matrix; generation is then applied to A_k, a new algebra, step by step (fine, the proof is per step); (c) R steps are F steps on the dual algebra, so they need the same hypothesis on the dual (checked in path13: 7 R-steps, no loops). A counting argument (n summands, K0 class matrix of det +-1) is also satisfied but is vacuous evidence: det = +-1 holds by construction ([T_v] = sum [P_h] - [P_v]), so it does not test anything beyond the proof. Finite global dimension plays no role. Would be refuted by: a step whose quiver has a loop at v (the code would then not be tiltingPlus as described).

## Evidence

Script `theorist_gen.py` checks the proof's hypotheses and records the Hom signs (own Hom-complex code from rounds/050 skeptic_tilt.py, independent of tiltingPlus/perI):

| set | steps | acyclic, targets != v, T has n summands | det of class matrix | Hom(T,T[1]) = 0 | Hom(T,T[-1]) = 0 |
|---|---|---|---|---|---|
| LNAs + duals, n = 5 | 112 | 112 | 112 of +-1 | 112 | 70 |
| n = 6 | 420 | 420 | 420 | 420 | 252 |
| n = 7 | 1584 | 1584 | 1584 | 1584 | 924 |
| E-163 path, 13 edges | 13 | 13 | 13 | 13 | 13 |
| 16 failing steps (class 1, n = 7) | 16 | 16 | 16 | 16 | 0 |

(Dual algebras included in the LNA rows; every vertex with at least one out-arrow, not only gate-admitted vertices, so the J != 0 rows there are not all gate-admitted steps.) Failing-step Hom(T,T[-1]) dimensions: 1 x14, 2 x2. The failing set is rebuilt from the E-151 walk (20 000 expansions, 519 s) and agrees with E-167's 16.

## Reproduction

```
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 5      # 2 s
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 6      # 3 s
.venv/bin/python workshop/rounds/054/theorist_gen.py lna 7      # 8 s
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # 519 s (mkdir -p /tmp/tsm first)
timeout 10m .venv/bin/python workshop/rounds/054/theorist_gen.py path13 > workshop/rounds/054/theorist_gen_path13_out.txt   # 2 s
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 2 20000 /tmp/tsm/c2.pkl 100   # 392 s
timeout 10m .venv/bin/python workshop/rounds/054/theorist_loops.py /tmp/tsm/c1.pkl 1   # loop check, class 1 (214 edges)
timeout 10m .venv/bin/python workshop/rounds/054/theorist_loops.py /tmp/tsm/c2.pkl 2   # class 2 (156 edges)
```

## Prior record

E-130/E-128/E-138 use "AI 2.31: mutation is silting" (cited from memory; the literature summaries are marked UNVERIFIED in STATE T5) and E-161/E-166/E-167 say "generation assumed (Okuyama-Rickard)". So the content "T is silting, hence generates" is implicit in the record, but no entry gives a proof or a check of the hypotheses; this submission gives the elementary proof that avoids AI 2.31 and the check. It is therefore a closing of a caveat, not a new phenomenon, and the argument itself is already stated in `research/literature/rickard-morita-theory-derived-categories.md` (about line 78), which this submission did not cite originally. Not in RETRACTIONS.

## Code changed

None (new scripts `workshop/rounds/054/theorist_gen.py`, `theorist_loops.py`, outputs `theorist_gen_path13_out.txt`, `theorist_loops_out.txt`). No tests touched.

## Next

- toolsmith/skeptic: the real residual gap of T10 is (a), End(T) = next algebra as algebras on the E-157/E-160 paths, class 2, and the 8 parallel-arrow failing steps (needs matrix-valued arrows); generation can now be struck from their caveat lists.
- skeptic: check that the code's mutation never meets a loop (acyclicity of every intermediate child on the walk; the script asserts it for the 29 edges only).
- Reword STATE T10(b)/(c): "Generation still assumed" -> "generation proved per step; premise = Hom(T,T[-1]) = 0 plus End(T) iso".
