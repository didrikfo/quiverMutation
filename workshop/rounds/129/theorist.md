# Every one of E-323's 50 key-kept J != 0 steps is a non-tilting irreducible silting mutation with dim Hom(T,T[-1]) = J, and its child has the Cartan matrix of End(T); the equivalence itself is already a theorem on record (E-130, E-171), not new

author: theorist · round: 129 · kind: rediscovery (theorem) + result (check at new scope)
thread: T10 · bears on: H-032, H-023, R-017, E-323, E-327, E-328, E-130, E-138, E-161, E-168, E-171
scope: proof: any acyclic quiver with admissible ideal, one vertex k, one step. Check: the 50 recorded key-kept J != 0 steps of E-323 (n = 11, 12; 25 F, 25 R; from `toolsmith_jxclass_out.jsonl`, not re-searched), 2 replayed crossing walks (6 J = 0 steps), 4 hand controls, 3 planted errors. Not covered: cyclic parents (E-329), algebra-level (quiver + relations) identity of child and End(T), any step not in that file.
power: the Hom test says NO (not tilting, H[-1] = 1) on E-080 (abde = acde at d), where NO is known; YES on abd = acd and on the 6 in-family J = 0 steps; the silting test says NO on 2 planted non-approximations (H[1] = 1, 2); the End(T)-vs-child comparison says NO on a planted wrong child.

## Claim

(1) **Proved (relative to Aihara-Iyama 1009.3370 Def 2.30, Thm 2.31, Thm 2.32(b), as read in E-171; no new input).** Let a = kQ/I, Q acyclic, I admissible, k a vertex, f: P_k -> M' = (+)_{beta: k -> t} P_t given by the arrows out of k, and T = Y + (+)_{i != k} P_i with Y = cone-type two-term complex [P_k -f-> M']. Then f is a minimal left add(a/P_k)-approximation, so T is the irreducible one-summand silting mutation of a at P_k; Hom(T, T[m]) = 0 for m != 0, -1; and Hom(T, T[-1]) = (+)_{i != k} J_i, J_i = ker(p |-> (p beta)_beta) on e_i a e_k. Hence for every step of an acyclic walk: **J != 0 <=> T is silting and not tilting**, and dim Hom(T, T[-1]) = J. This is E-130 (M2) plus E-171's reading of AI 2.32(b); I re-derived it before finding them. Also already true by code: the library's J (`perI`) and its injectivity test (`isTiltingAtVertex`, `tiltingPlus`) compute the rank of the same matrix, so "J != 0 <=> injectivity test fails" is a tautology of the implementation, not something E-327 is needed for.
(2) **Checked, at new scope (n = 11, 12).** By an independent Hom-complex computation (chain maps modulo homotopy, exact rationals, no call to `perI`/`tiltingPlus`): on all 50 E-323 steps, f is an approximation, Hom(T,T[1]) = Hom(T,T[2]) = Hom(T,T[-2]) = 0, Hom(T,T[-1]) = J (1 on 48 steps, 2 on the F3 step at 0040033030 and its dual), Hom(Y,Y[-1]) = 0, and the matrix dim Hom(T_x, T_y) (Cartan of End^0(T)) equals the exact Cartan matrix of the walk's child (child^op for R steps) up to simultaneous permutation, 50 of 50 (never as the transpose).
(3) **Not claimed, and what this does NOT give.** (a) Child = End(T) as algebras: only Cartan level, unlabelled. (b) That a non-tilting step leaves the derived class: false in general (E-316's 19 J != 0 F-010 steps stay in class); non-tilting says only that *this* T gives no equivalence. The inequivalence of E-323 still rests on SNF(S) via the H-024 necessity step, so this does not discharge E-323's F-blocker by itself; it does supply its "Hom-level non-tilting check of F3 at 0040033030" (H[-1] = 2). (c) The agenda's corollary "option-on `refuseNonTilting` walks cannot leave a class (theorem relative to AI + E-327)" is **misstated**: it needs J = 0 => child ≅ End(T), which is H-023 (the J = 0 premise), not AI or E-327. Correct form: relative to AI 2.31/2.32 and H-023, on acyclic parents, an option-on walk is a chain of tilting mutations and stays in one derived class. Refuted by: a J != 0 step with Hom(T,T[-1]) = 0, or an acyclic admissible step whose arrows out of k fail to approximate.

## Evidence

Proof of (1), the step I am least sure of marked (*).
- Approximation: a map P_k -> P_j (j != k) is a class of paths k ~> j; Q acyclic, so each has positive length and begins with an arrow beta out of k, i.e. it is f followed by a map P_t -> P_j. Minimality (*): the components of f are distinct arrows, linearly independent modulo rad^2 because I is admissible; an endomorphism phi of M' with f phi = f is then the identity modulo the radical, hence invertible. (Standard; I did not find it written in AI for this case, and parallel arrows are covered only by this argument.) So T = mu_{P_k}(a) of AI Def 2.30 with #ind X = 1: irreducible by definition, silting by Thm 2.31.
- Negative degrees (degree bookkeeping, Y in degrees -1, 0): Hom(P_i, P_j[m]) = 0 for m != 0; Hom(Y, P_j[-1]) = 0 and Hom(T, T[<= -2]) = 0 for degree reasons; Hom(Y, Y[-1]) has only the component M' -> P_k, and Hom(P_t, P_k) = 0 for t a successor of k (acyclic); Hom(P_i, Y[-1]) = {g: P_i -> P_k : g f = 0} with no homotopy = J_i. Total: dim Hom(T,T[-1]) = sum_i J_i (E-130 (M2), acyclic case).
- Consequence for Cartan (E-138, recorded): with M the K_0 change of basis, C(End^0 T) = M C_a M^T + N, N supported in the (P_i, Y) entries with N_{iY} = J_i. So S' = M S M^T + N + N^T: the SNF change of E-323 is this rank <= 2 perturbation, and E-316 shows it can be invisible.

Checks (`theorist_jsilt.out`; columns: J = own copy of perI, lib = `isTiltingAtVertex`, appr, H[m] = dim Hom(T,T[m]), End0 ~ child):

| case | J | lib tilt | H[-2,-1,1,2] | silting | tilting | End0 ~ child |
|---|---|---|---|---|---|---|
| E-080 abde = acde @4 (J != 0 known) | 1 | F | 0 1 0 0 | T | F | T |
| abd = acd @4 (tilting known) | 0 | T | 0 0 0 0 | T | T | T |
| LNA A5/(abc) @2 / @3 | 0 / 1 | T / F | 0 0 0 0 / 0 1 0 0 | T | T / F | T |
| planted: 1->2->{3,4} @2, one arrow dropped from Y | - | - | 0 0 1 0 | **F** | F | - |
| planted: E-080 @1, one arrow dropped | - | - | 0 0 2 0 | **F** | F | - |
| planted: LNA @2 compared with child @3 | 0 | T | 0 0 0 0 | T | T | **F** |
| walk 0040033030 F3 R9 F3 (J = 0, 3 steps) | 0 | T | 0 0 0 0 | T | T | T |
| walk 0040033030 ... F3 (the E-323 crossing) | 2 | F | 0 2 0 0 | T | F | T |
| walk 003033030 F3 R8 F3 (3) then F3 | 0..0, 1 | T, F | 0 0 0 0, 0 1 0 0 | T | T, F | T |
| **E-323's 50 steps** (25 F, 25 R) | 1 (48), 2 (2) | F 50/50 | 0 J 0 0, 50/50 | 50/50 | 0/50 | 50/50 (transpose 0/50) |

The replayed walks rebuild the recorded paths from the LNA roots with the walk's semantics (gate, illegal-relation drop, reduce, op for R); every replayed step was gate-admitted. The "appr" column on planted rows refers to the full arrow set; the planted error is caught by H[1] > 0, as it should be (dropping an arrow kills the approximation property).

On the E-328 certificate: I did not run `scholar_irred.py` on these steps, and it cannot be run as written: it needs the complete list of indecomposables of D^b (available for B_n ~ D_n, type finite) to rule out a silting object strictly between; E-323's parents at n = 11, 12 have no such list. It is also not needed: E-328 had to certify exchanges found by search, whereas every walk step is by construction the mutation at one indecomposable summand with an explicit minimal approximation, which is AI's definition of irreducible (and a covering in the silting order by Thm 2.35).

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/129/theorist_jsilt.py all > workshop/rounds/129/theorist_jsilt.out   # 4 s
```
(`controls` or `steps` alone; reads `workshop/rounds/125/toolsmith_jxclass_out.jsonl`.)

## Prior record

Not new as mathematics: E-128/E-130 (Hom(T,T[-1]) = sum J_i, "silting-not-tilting iff some J_i != 0"), E-171 (AI 2.31/2.32(b) checked against the LaTeX; whole-T reading), E-138 and E-095 (Cartan defect N at (v, i) = J_i), E-161 and E-168 (independent Hom tests at n = 7: 25 and 16 failing J != 0 steps have Hom(T,T[-1]) = 1 or 2; generation). New: the same independent test at n = 11, 12 on E-323's 50 crossing steps; the End^0(T) = child Cartan match on J != 0 steps at this scope; the observation that the round-128 agenda's corollary rests on H-023, not on AI + E-327; that the E-328 certificate is inapplicable and unnecessary here.

## Code changed

None in the library. New: `workshop/rounds/129/theorist_jsilt.py` (Hom-complex engine for two-term complexes over an acyclic path category, controls, step replay) and its output `theorist_jsilt.out`. No tests touched.

## Next

- T10 is, as a mechanism question, answered by the record: J != 0 is exactly "the vertex silting mutation is not tilting". The live questions are (i) H-023: J = 0 => child ≅ End(T) at algebra level (toolsmith: an End(T)-presentation builder to compare quivers with relations, as E-167 did on 13 edges); (ii) whether a J != 0 child can be in-class (E-316 says yes), so whether `refuseNonTilting` costs connectivity (toolsmith: F-041 n = 10, 11 with the option on).
- Chair: amend the agenda's corollary to "relative to AI 2.31/2.32 and H-023, acyclic parents".
- Skeptic: E-323's remaining blocker is the second invariant on the pair; the Hom-level part is supplied here (H[-1] = 2 at F3, 0040033030).
- Cyclic parents: E-130 (M3) says Hom(N,N[-1]) may add; my engine asserts acyclicity, so the E-329 cycle-upstream cases need a different basis (non-simple paths).
