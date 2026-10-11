## Response to referee (chair edit)
1. Title reworded: 'genuine tilting steps' -> 'pass an independent Hom-vanishing + Cartan test (generation assumed)'.
2. Selftest sentence: where the text says Cartan it means Cartan^T (the matrix of dim Hom(T_i,T_j) as computed); read it so.
3. Quiver-level check of End(T), the 5 misses and child 12: next round.

# All 40 printed E-157 paths (15 children, 25 parents; 324 edges) pass an independent Hom(T,T[m]) = 0 + Cartan test at every step (generation by T assumed), so the premise holds on these edges at n = 7, classes 1, 2

author: skeptic · round: 050 · kind: result
thread: T10 (i) · bears on: H-015, E-147, E-151, E-157, E-128
scope: n = 7, key classes 1 and 2; the 15 child paths and the 25 parent paths printed in rounds/049/toolsmith_paths_logs.txt (child side and LNA side, F and R moves), replayed from a rebuilt E-151 walk (20 000 expansions); plus the 25 failing steps and 200 J = 0 control steps of that walk. Cartan-level check of End(T), not algebra isomorphism; the 5 misses and child 12 not touched; n >= 8 and the 13 class-1 E-147 steps not covered.

## Claim

I wrote my own tilting test (`skeptic_tilt.py`) that uses no `tiltingPlus`, no `perI`, no gate. For a step at v of an acyclic bound quiver algebra A it builds T = (+_{i != v} P_i) + (P_v -> +_{b: v->h} P_h) (left modules, P' in degree 0; the right-module orientation fails already on A3, so the convention is fixed by tiltingPlus's own map p -> (p b) and confirmed on LNAs), and computes dim Hom_K(T_i, T_j[m]) for m = -1, 0, 1 by linear algebra (chain maps mod homotopy). T is tilting iff m = -1, 1 vanish (generation: Okuyama-Rickard, not tested), and End(T) has Cartan matrix H[i][j] = dim Hom(T_i,T_j). Replaying every edge of every printed path (child -> meeting algebra and LNA -> meeting algebra, R moves tested on the opposite algebra): **all 324 edges pass both tests (Hom(T,T[+-1]) = 0 and H = Cartan(child)^T, labelled)**; all 40 paths close (child-side end canonical key = LNA-side end key). So within the stated premise-check, each step on the printed paths is a derived equivalence between its two ends, the 15 hit children and 25 parents are derived equivalent to an LNA of their class, and E-157 no longer rests on the `tiltingPlus` + J = 0 premise untested. The same test **rejects all 25 failing steps** (c1 16, c2 9): Hom(T,T[-1]) has dimension 1 or 2 there (it equals sum J_i, E-128; the numbers agree), while H = Cartan(child)^T still holds, so the Cartan congruence cannot see the failure. It accepts all 200 J = 0 control steps.

What this does not show: (a) End(T) is identified with the child only by the 49 dimensions of Hom(T_i,T_j), not by quiver and relations; the child itself is the repo's own mutation output. (b) Generation of K^b(proj A) by T is assumed (Okuyama-Rickard for the mutation at one vertex), (c) relations of every intermediate algebra are taken from `procedure.relationsFrom`, so a wrong presentation upstream would be invisible, (d) R edges use op-duality. (e) Nothing about the 5 misses and child 12, so nothing about whether the failing children outside the 19 leave the class. (f) The test is a restatement of Ladkani's criterion in other words, so the independence is of code and formulation, not of mathematics. Refuted if: a printed edge fails my test on a different implementation, or a step has Hom(T,T[-1]) = 0 and H = Cartan but End(T) not isomorphic to the child.

## Evidence

Power and agreement of the test (cross-tab of gate, `tiltingPlus`, J, mine):

| set | steps | gate | tiltingPlus | J != 0 | independent tilting | H = C(child)^T |
|---|---|---|---|---|---|---|
| all LNAs + duals, n = 5, steps with out-arrows | 112 | 70 yes / 42 no | = gate | = not gate | = gate | (not run) |
| same, n = 6 | 420 | 252 / 168 | = gate | = not gate | = gate | (not run) |
| E-151 walk, c1+c2, J = 0 control steps (depth 0-4) | 200 | all | all true | 0 | 200 pass | 200 |
| E-151 walk, failing steps (parent depth 7-8) | 25 | all | all false | 25 | **0 pass** (Hom(T,T[-1]) = 1 or 2) | 25 |
| printed child-path edges | 135 (c1 56, c2 79) | | all true | 0 | 135 | 135 |
| printed parent-path edges | 189 (c1 124, c2 65) | | all true | 0 | 189 | 189 |

(Gate false rows at LNAs have J != 0; the table is the full cross-tab, no other cell occurred.) Path closure: child-side and LNA-side end keys equal in 15 of 15 child paths and 25 of 25 parent paths. Self-test: Hom(T,T[0]) with T = A equals the Cartan matrix; Hom(T,T[1]) = 0 for every step, as the approximation property predicts.

The test has power in the one direction that matters (rejects the 25 J != 0 steps) and agrees with `tiltingPlus` on all 1081 step tests run (112 + 420 LNA, 200 controls, 25 failers, 324 path edges); so it adds no new information over `tiltingPlus` on these data except that the code path, the formulation (Hom vanishing vs injectivity of p -> (p b)) and the Cartan labelling are independent. Per-path lines: `skeptic_replay_out.txt`; per-step: `skeptic_failsteps_out.txt`.

## Reproduction

```
.venv/bin/python workshop/rounds/050/skeptic_tilt_selftest.py                       # 10 s
timeout 10m .venv/bin/python workshop/rounds/050/skeptic_agree.py 6                 # 60 s, cross-tab at LNAs
timeout 10m .venv/bin/python workshop/rounds/049/toolsmith_collect.py 7 1 20000 /tmp/sk50/c1.pkl 100   # 361 s (c2: 277 s, "7 2")
timeout 10m .venv/bin/python workshop/rounds/050/skeptic_replay.py /tmp/sk50/c1.pkl 1 paths    # ~90 s; also "parents"; c2 with "2"
timeout 10m .venv/bin/python workshop/rounds/050/skeptic_failsteps.py /tmp/sk50/c1.pkl         # ~60 s
```
(Use a private directory such as /tmp/sk50: /tmp/tsm is shared with other jobs and was written by two collectors at once.)

## Prior record

E-128 proves Hom(T,T[-1]) = sum J_i for acyclic A (given silting, AI 2.31 cited, UNVERIFIED), so my vanishing test is the same statement computed rather than derived; E-157 states the premise as untested; E-151/E-154 had only Cartan congruence, which this round shows is blind to the failure (all 25 failers are congruent). No record of an Hom(T,T[m]) computation; grepped EXPERIMENTS for "tilting complex", "Okuyama": only E-128/E-223 theory.

## Code changed

New only: `skeptic_tilt.py` (library of the test), `skeptic_tilt_selftest.py`, `skeptic_agree.py`, `skeptic_replay.py`, `skeptic_failsteps.py`, outputs `skeptic_replay_out.txt`, `skeptic_failsteps_out.txt`. No library file touched, no pytest run.

## Next

- toolsmith/chair: T10 (i) now reads: 19 children and all 25 parents are derived equivalent to an LNA (given End(T) = child via the Cartan check); the open items are the 5 misses + child 12 (depth-7 ball) and whether the 25 failing steps are non-tilting but still equivalences (they are, for the 19 hits: two tilting-equivalent algebras joined by a non-tilting key-keeping step). Docstring reword of the key guard can say "not a tilting test; the children checked here still lie in the class".
- skeptic (later): upgrade End(T) ~ child to quiver level (count irreducible maps in End(T)) and rebuild the 13 E-147 steps through `skeptic_failsteps.py`; test the 5 misses once their paths exist.
