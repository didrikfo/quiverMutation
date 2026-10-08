# Gate-admitted J != 0 steps with H1, H2 and x^2 coefficient != 1 exist, and at n = 7 classes 1 and 2 some are key-preserving (D = 0), so "J != 0 moves the key" is a class-0 statement

author: skeptic · round: 043 · kind: result (counterexample to E-140's reading; independent build)
thread: T5 · bears on: E-138, E-140, E-141, E-143, H-015

## Claim

(1) Independent of `theorist_*.py`: on guard-walks from the n = 6, 7 class-0 LNAs (316 and 300-317 distinct J != 0 steps) every step has H1/H2, an orbit relation e_i = F^s e_w (s = 1 or 6 = -2 mod 8 at n = 6; s = 1 or 3 at n = 7) and D(x) = det(xC_B+C_B^T) - det(xC_A+C_A^T) with lowest term exactly x^2, coefficient 1. E-143's observation reproduces; no c_2 != 0 inside class 0.
(2) Class 0 is not typical. At n = 7, key-guarded walks (500 s, depth <= 30) in class 1 (67 distinct J != 0 steps) and class 2 (64) contain **D = 0 steps: 13 and 9** (12 and 9 with H1, H2 in the strict sense). These are gate-admitted, J = e_i, key-preserving, and still occur when the walk is also restricted to tiltingPlus steps (12 and 9 of 66 and 64). So E-140's "none of the J != 0 steps keeps the class key" (c1: 2 steps, c2: vacuous) was a depth/time artefact; the law does not hold at n = 7 c1, c2.
(3) The orbit relation can fail: in c1, 44 of 66 H1/H2 steps have no s (F^s e_w = e_i has no solution for |s| <= 40), 52 of 52 in c2 with D low term x^3 (coefficient 1); in c2 the nine D = 0 steps do have it (s = 10). With random acyclic parents (m = 5, 6, 7) it fails in most strict H1/H2 steps (m = 6: 19 of 23 hits with D = 0 and 800+ others lack s) and then the lowest term is x^1 (coefficient mostly 1).
(4) Where the orbit relation holds on non-LNA-walk parents, x^2 coefficients 2, 3 and -1, -2 occur (random m = 6, 7: 33 strict H1/H2 steps with an orbit relation and lowest term x^2, coefficient != 1), i.e. c_2 != 0 is realisable at the matrix-plus-gate level; none of them has D = 0.
Not claimed: that any D = 0 child is derived equivalent to its parent. Every J != 0 step fails tiltingPlus by definition (all 700+ here: `(tiltingPlus, D == 0)` tables), the four D = 0 steps rebuilt by hand are Cartan-incongruent (R C_A R^T != C_B), as in E-137. The key guard is therefore useless as a test of derived equivalence off J = 0 steps; that is H-015's guard, not its conclusion.

## Evidence

Own code: `skeptic_c2.py` takes the library only for gate, mutation, reduction, Cartan matrices and J (the `perI` of rounds/033). It computes D(x) by exact interpolation of det(xC+C^T), tests H1/H2 on C_B (C_B = [[Z,0],[e_w^T,1]], Z = C_A off v), the strict version (w != i and det poly of C_B - E_vi equals that of C_A), and the orbit exponent from F = Z Z^-T (no use of theorist code).

| set | walk | distinct J != 0 | H1/H2 | D lowest term | orbit s |
|---|---|---|---|---|---|
| n=6 all 4 classes seeded, guard off, depth 10, 540 s | no key guard | 316 | 308 | (2,1) in 316 | 1: 88, 6: 220 |
| n=6 classes 1-3, key guard, depth 10, 520 s | key | 0 | - | - | vacuous |
| n=7 c0, key guard, 500 s | key | 300 (317 in another run) | 292 | (2,1) in 300 | 1: 276, 3: 16 |
| n=7 c1, key guard, 500 s | key | 67 | 65 | (1,-1) 2, (2,1) 48, (3,2) 4, D=0 13 | 1: 22, none 44 |
| n=7 c1, key + tiltingPlus guard | tilting | 66 | 65 | same, D=0 12 | same |
| n=7 c2, key guard, 500 s | key | 64 | 61 | (3,1) 52, D=0 9 | 10: 9, none 52 |
| n=7 c2, key + tiltingPlus guard | tilting | 64 | 61 | same | same |
| random acyclic m=5/6/7 (7 min each) | none | 582 / 6237 / 2638 | see script output | x^1 dominant | see script |

The c2 D = 0 example, rebuilt by hand from arrows and relations (`skeptic_zero.py`): n = 7, arrows 3->1, 3->4, 4->5, 5->6, 1->6, 6->2, 7->1; relations (3->1->6->2) = (3->4->5->6->2) and (7->1->6) = 0; v = 6 (out-degree 1); J_3 = 1; gate True, tiltingPlus False; child arrows 1->2, 2->6, 3->1, 3->4, 4->5, 5->2, 7->1, relations 1->2->6 = 0, 3->1->2 = 3->4->5->2, 5->2->6 = 0, 7->1->2 = 0; key of parent = key of child = (1,1,0,-2,-2,0,1,1); dim A 22, dim B 17; Cartan incongruent.
Traps checked: key comparison is on the child (both library key and my D agree on 21 steps); parents of the D = 0 steps are not LNAs (they are walk descendants), the seeds are all LNAs of the class and duals. The c1 pickle does not rebuild by `skeptic_zero.py` (parallel arrows 7->4 break the vertex-sequence relation format), so the c1 D = 0 steps rest on `skeptic_c2.py` alone.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/043/skeptic_c2.py 6 0 100 540 14           # guard off, 540 s
timeout 10m .venv/bin/python -u workshop/rounds/043/skeptic_c2.py 7 C C+1 500 30 guard     # C = 0, 1, 2; 500-520 s each
timeout 10m .venv/bin/python -u workshop/rounds/043/skeptic_c2.py 7 C C+1 520 30 tguard    # C = 1, 2 (adds tiltingPlus filter)
timeout 10m .venv/bin/python workshop/rounds/043/skeptic_zero.py workshop/rounds/043/skeptic_zero_n7c2.pkl   # 30 s, hand rebuild
timeout 10m .venv/bin/python -u workshop/rounds/043/skeptic_rand.py M 11 420                # M = 5, 6, 7
```
(Walks are time-capped, so counts move by a few percent with load; the D = 0 steps appeared in every run of c1 and c2: g, h, i, t runs.)

## Prior record

E-143 (H1/H2, orbit relation, c_2 = 0 observed, class 0 only); E-140 (278 steps, none key-preserving; n=7 c1 2 steps, c2 vacuous); E-138, E-141 ("J != 0 implies key moves", false for general parents). E-137 notes a key-preserving non-tilting step at n = 6 out of an E-134 hit (kernel dim 2). New: dim-1 key-preserving J != 0 steps with H1/H2 at n = 7 c1, c2 reachable from LNAs by key (and tilting) steps; failure of the orbit relation; the c_2 != 0 matrices. Grep of research/ for "key-preserving" and "orbit relation": nothing contradicting. RETRACTIONS: not touched.

## Code changed

None in the library. New: `skeptic_c2.py`, `skeptic_rand.py`, `skeptic_zero.py` and two small pickles in `workshop/rounds/043/`. No tests run.

## Next

- Theorist: P1 gives Q = xB; D = 0 needs B = 0 = 1 - c_s - c_{-s} with s = 10 in c2. Find which orbit data realise it; do the c2 steps need a long F-cycle (F of order > 8 here)?
- Experimentalist: n = 8 c1, c2 tally with `skeptic_c2.py 8` after `--plan`-style sizing (n = 7 c1, c2 took 500 s); is the key-preserving step reachable from the LNA by a path whose earlier steps are all J = 0 (tilting) -- the tguard run says yes at the key level; replay one path with Cartan congruence.
- Skeptic/referee: any "key guard excludes J != 0" claim (E-138 law, E-140) must say class 0; E-140's c1, c2 cells should be marked non-vacuous.
