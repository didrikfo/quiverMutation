# The nine J = 0 steps of F-057 (a)'s 10-step walk are tilting by a full Hom computation, by an own-code Hom(T,T[-1]) and by `isTiltingAtVertex`, with no disagreement; all three say NO on the walk's J != 0 step and on 22 further known non-tilting steps (n = 11 walk; controls n = 5, 9, 11)

author: toolsmith · round: 139 · kind: result (replication; largely forced by E-351)
thread: T10 · bears on: F-057, H-032, H-023, E-332, E-333, E-351, E-316
scope: n = 11; the recorded E-332 walk dual(303305000) · F6 F11 F7 F11 F7 F7* F3 F3 F4 F7 (9 J = 0 steps) plus the second root 303500030 · F1 F1 F4 F3 F1 F8* (5 J = 0 steps); one vertex, F side; every parent acyclic and free of illegal relations (checked). Child vs End(T) compared at Cartan level only (up to permutation), not as algebras. Not covered: any other walk at n = 11, R steps at n = 11, cyclic parents (E-329), whether the library child IS End(T) (H-023). Reading of the assignment: "the 9 J = 0 steps of F-057" are E-332's (F-057 "What it does not say"); E-316's 19 steps are all J != 0, so they are used as NO controls.
power: NO controls, all said NO by all three verdicts: the walk's own J != 0 steps P_b·F7 and P_a·F8 (2 of 2, H[-1] = 1); E-316's 19 n = 9 key-kept J != 0 steps, replayed from their recorded paths (19 of 19, H[-1] = J = 2, 3, 4); E-080 abde = acde at 4 (1 of 1, H[-1] = 1). YES control abd = acd (1 of 1). Total 22 NO controls, 0 misses.

## Claim

At each of the nine J = 0 steps of F-057 (a)'s 10-step pure-F walk, the Okuyama-Rickard complex T has Hom(T, T[m]) = 0 for m = -2, -1, 1, 2. So each is a tilting mutation, and the three verdicts agree: LIB (`isTiltingAtVertex`), HOM (the round-129 Hom-complex engine) and OWN (my own code). The same holds at the five J = 0 steps of the second root. At every one of these 14 steps, the dim Hom(T_x, T_y) matrix equals the library child's Cartan matrix up to permutation. This closes F-057's open item "a Hom check of the nine J = 0 steps". The tilting half of F-057 (a) (root and P_b derived equivalent, and X_b and 403033030 derived equivalent) now rests on Hom computations and not only on the injectivity test. One caveat remains: each step gives D(parent) ≃ D(End T). Identifying End T with the library's child is H-023, and here it is checked at Cartan level only. **Not claimed:** that the child is End(T) as an algebra, or anything about other walks. This outcome was expected. Every parent is acyclic with no illegal relation, so LIB = HOM is E-351's theorem. E-333's referee found 7 498 of 7 498 random J = 0 steps tilting. What is new is only these 14 recorded steps. **Refuted by:** a Hom(T,T[-1]) != 0 at any of the nine, or a run in which a parent of the walk is cyclic or has an illegal relation.

## Pre-registration

Written to `workshop/rounds/139/toolsmith_prereg.txt` before the first run, and unchanged since. It contains the targets, the three verdicts, the expected outcome, and the meaning of each kind of disagreement:
- D1: LIB True but HOM/OWN False on a J = 0 step. This would contradict E-351 (an illegal relation or a bug) and would un-localise F-057's crossing.
- D2: LIB False but HOM True. This would mean a second J != 0 step and a misrecorded E-332.
- D3: OWN alone differs. That is my bug.
- D4: all three say tilting, but End^0(T) Cartan != child Cartan. That would be an H-023 finding.
- D5: HOM says not silting. That would be an engine bug.

The controls C1 to C3 are as listed under *power*. **Outcome: none of D1 to D5 occurred. 0 disagreements in 139 step evaluations.**

## Evidence

The three verdicts per step (parent A, vertex k):
- **LIB** = `mutation.isTiltingAtVertex(A, k)`.
- **HOM** = `theorist_jsilt.test`, imported unchanged. It computes chain maps modulo homotopy over Q for T = Y + Σ_{i≠k} P_i, giving H[m] for m = -2..2. It also checks that the arrows out of k approximate, and compares End^0(T) with the child.
- **OWN** = new in `toolsmith_homj0.py`. It enumerates paths itself and takes the ideal as the span of all p·r·q over exact rationals, with no `idealBasis`. It computes Hom(T,T[-1]) as Σ_{i≠k} dim ker(e_iAe_k → ⊕_β e_iAe_{t(β)}) plus Hom(Y,Y[-1]) = Σ_β dim e_tAe_k, where the latter is computed and not assumed to be 0. It also computes the approximation (H[1] = 0) in the same quotient.

| step (walk b, dual 303305000) | J | LIB | HOM H[-2,-1,1,2] | OWN H[-1] | End0 ~ child |
|---|---|---|---|---|---|
| d1 F6 | 0 | T | 0 0 0 0 | 0 | yes |
| d2 F11 | 0 | T | 0 0 0 0 | 0 | yes |
| d3 F7 | 0 | T | 0 0 0 0 | 0 | yes |
| d4 F11 | 0 | T | 0 0 0 0 | 0 | yes |
| d5 F7 | 0 | T | 0 0 0 0 | 0 | yes |
| **d6 F7 (control C1)** | 1 | F | 0 1 0 0 | 1 | yes |
| d7 F3 | 0 | T | 0 0 0 0 | 0 | yes |
| d8 F3 | 0 | T | 0 0 0 0 | 0 | yes |
| d9 F4 | 0 | T | 0 0 0 0 | 0 | yes |
| d10 F7 | 0 | T | 0 0 0 0 | 0 | yes |

| class | steps | LIB / HOM / OWN | all agree | End0 ~ child |
|---|---|---|---|---|
| TARGET walk b, J = 0 | 9 | tilting ×3 | 9 | 9 |
| walk a 303500030 d1-d5, J = 0 | 5 | tilting ×3 | 5 | 5 |
| C1 walk J != 0 (b d6, a d6) | 2 | NO ×3 (H[-1] = 1) | 2 | 2 |
| C2 E-316 key-kept J != 0, n = 9 | 19 | NO ×3 (H[-1] = J: 2 ×10, 3 ×5, 4 ×4) | 19 | 19 |
| C2 path moves to those parents (J = 0) | 103 | tilting ×3 | 103 | 103 |
| C3 E-080 | 1 | NO ×3 (H[-1] = 1) | 1 | 1 |
| YES abd = acd | 1 | tilting ×3 | 1 | 1 |

On every step: H[1] = H[2] = H[-2] = 0 (silting), Hom(Y,Y[-1]) = 0, and no parent has an illegal relation. All 19 E-316 paths replay through the library gate without a break. Each recorded step is J != 0, as E-316 says. The 103 J = 0 path moves (shared prefixes counted again) are a side result: the q_S-constant joins of E-316/E-321 are tilting at Hom level on their parent paths. The final J = 0 legs of those joins after the J != 0 step were not replayed.

Toolsmith's referee point: all three verdicts compute the same quantity. On an acyclic parent, Hom(T,T[-1]) = J. HOM and OWN reach it through different code (complex-level chain maps vs a direct kernel), and LIB and OWN through different ideal code (library `idealBasis` vs own p·r·q span). So agreement tests the code, not the mathematics, and the mathematics is E-351. OWN's approximation check is trivially true for "all arrows out of k" and has no NO control of its own. HOM's H[1] has planted controls in E-333.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/139/toolsmith_homj0.py all   # 4 s; output toolsmith_homj0.out (27 KB)
```
Modes: `walks`, `controls`, `e316`.

## Prior record

- F-057 lists this item as open.
- E-332 recorded the J values of the walk via `perI`/`isTiltingAtVertex`, and its referee replayed SNF and Cartan, not Hom.
- E-333 Hom-checked the 50 E-323 steps and 2 E-323 walks, not E-332's. Its referee checked 7 500 random steps (rounds/129/theorist.review.md), not these.
- E-351 proves LIB = Hom on acyclic quivers when no relation is illegal. That is the case at all 14 parents here, so the outcome is a theorem applied plus a record closure. It is not a discovery.
- E-316's 19 steps were previously q_S-checked and refused by `refuseNonTilting` (E-323, E-350). Their Hom(T,T[-1]) = J values (2-4) are not recorded elsewhere as Hom computations, as far as `grep` shows. They are controls here, not a claim.

## Code changed

No library change. New files: `workshop/rounds/139/toolsmith_homj0.py` (imports `workshop/rounds/129/theorist_jsilt.py` unchanged), `toolsmith_prereg.txt` and `toolsmith_homj0.out`. No tests were run, because no library file was touched.

## Next

- Chair: F-057's "What it does not say" can drop "a Hom check of the nine J = 0 steps is open". The remaining gap for (a)'s tilting half is H-023 (child = End(T)), checked at Cartan level here.
- Theorist: an algebra-level End(T) ≅ child check on one of these steps would turn "Cartan level" into an isomorphism. A small quiver, so feasible.
- Toolsmith (me): the T10 items still open are depth 8 (overnight), the 4 duals, and depth 5+ at E-316.
