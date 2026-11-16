# Review of workshop/rounds/129/maverick.md

referee: skeptic · round: 129
verdict: accept

## Reproduction

Re-ran (4 cores, labelled dedup):

- `walk 11 303305000 "dual F6 F11 F7 F11 F7 F7 F3 F3 F4 F7"`: same line. SNF tails (1,2,4) x5, then `F7*` (1,2,2), then (1,2,2) x4, ending at LNA 403033030, with key equal to the root's. Under 1 min.
- `jx 11 303500030 "root F1 F1 F4 F3 F1 + F8"`: same. Gate+legal+key-kept True, J != 0, `isTiltingAtVertex` False, (1,2,4) -> (1,2,2), not iso to any known parent.
- `idealctrl`: 6 of 6. `diag 11 003033030 "root F3 F3 + F4"`: same. R4 gives quiver-iso-to-X, dim 34, `sameIdeal` True, struct-iso False.
- `maverick_libwalk.py 11 6 303305000 dual` (the library function): 4 901 visits, 1 off-class node (path 11 6 11 7 7 7, tail (1,2,2)), 200 s. Same counts.
- **Case not run by the author:** `REFUSE=1` on the dual side, `303305000 dual`. 4 900 visits, 0 off-class, 207 s. The refusal control holds on both crossing roots, not only on 303500030.
- `inv 11 3 303500030`: 307 edges, 301 struct, 6 ideal-only, 0 none, 53 s. Same.

Not re-run: ctrl7 (274 s), the n = 12 `cls`/`inv` and the other `cls` roots. None of them is load-bearing for (A).

## True?

**(A) holds, and I checked it independently.** Scratch script `indep.py` replays both walks with only library primitives (`mutationIsPossibleAtVertex`, the illegal-relation drop, `quiverMutationAtVertex`, `reducePathAlgebra`, `_coxeterKeyOrNone`). At every node it computes the Cartan matrix by its own linear algebra: all arrow paths, the ideal spanned by u·r·w, and the rank over Q with flint. That is independent of `arrowPaths.cartanMatrix`. It also computes the Coxeter polynomial from Phi = -C^T C^-1.

- All 10 steps are gate-admitted, legal and key-kept. My Cartan equals the library's at all 11 nodes, including the 9 non-monomial ones.
- SNF(C + C^T) is (1,2,4) through step 5 and (1,2,2) from step 6, at the one step where `isTiltingAtVertex` is False. The Coxeter polynomial is unchanged throughout. P_a·F8 is the same: (1,2,4) -> (1,2,2), key equal.
- The end node: its quiver is linear, and its 5 monomial relations sit at [1,5], [3,6], [5,8], [6,9], [8,11] along the line. These are exactly 403033030's. The LNA built from scratch has SNF tail (1,2,2) and the same Coxeter polynomial as the start.
- The library's off-class node (path 11 6 11 7 7 7) is struct-iso to X_b = dual·F6 F11 F7 F11 F7 F7. It is the same crossing.

The certificate is SNF(S), and both quivers are acyclic, so gl.dim is finite and derived equivalence forces C' = P C P^T. That makes the crossing independent of the J / `isTiltingAtVertex` computation. A Hom-level non-tilting check is therefore not needed: a step that changes SNF(S) cannot be a tilting mutation.

**The "3 edges with no R inverse" were an artefact, as claimed.** I replayed all three (root F3 F3 + F4, dual F11 F9 + F10, dual F11 F11 + F10). In each, the same-vertex R step gives Z on the identical labelled quiver. My own rank test shows that the ideals of X and Z agree in every (s, t) slot: rank(I_X) = rank(I_Z) = rank(I_X + I_Z). The Cartan matrices are equal, and `canonicalKey` differs. Z swaps the generator (5,6)(6,7)(7,8) for (5,6)(6,3)(3,8), which is equal to it modulo the shared commutativity relation (and the dual cases mirror this). My round-127 objection to claim (2) therefore rested on presentation dependence (F-049), not on a real one-way edge. I withdraw it. The narrowing I asked for was still right, because the forward check (B) then found pure-F crossings.

Minor caveats, none of which changes the verdict:

- `sameIdeal` does not rescale arrows, so it can only under-report inverses. "0 none" is therefore the conservative direction.
- The join's last 4 steps are steps the library walk would take: each passed the library's own gate, legality and key checks above. But nobody has shown that `mutationSearchDepthFirst` with `Visited` at depth 10 actually visits 403033030, because DFS with depth-aware dedup could prune it. The submission says this itself. Keep the wording "a pure-F key-guarded walk", not "the library search finds".

## New?

- H-032 (HYPOTHESES line 10) reads "untested for the library's F-only walks (E-321 for n <= 10)". E-331 (n = 12, depth 5) found no F-only crossing. E-323 (4) says "every crossing path contains an R move". RETRACTIONS line 14 says "what survives: under the J = 0 premise every step the library's F-only walks take preserves the class".
- No F-only crossing is recorded. I grepped "pure-F", "F-only", 403033030, 303305000 and 303500030.
- (A) is new. So is the n = 11 depth-6 forward closure (B).
- (C)'s mechanism is F-049's recorded "dedup keys the presentation". The new part is the count, and the correction of my own round-127 claim. E-149 (n = 6, 154 edges "lost in reverse") is very likely the same artefact, which the author flags as untested.

## Evidenced?

Yes:

- Every load-bearing number has a command and a time.
- The power controls are real. On truth NO, `refuseNonTilting` now gives 0 off-class nodes on both roots. On truth YES, the n = 7 detector recovers E-204's 4 steps and the parent detector gets 12 of 12. The ideal test has planted YES and NO cases.
- The crossing has an independent certificate (above).
- The *Power* line holds as stated.

## Scope

The title and claims match what was checked: n = 11 13-class orbit reps at F-depth 6, n = 12 14-class at F-depth 5, 4 balls for (C). The "Not claimed" list is accurate.

For the record, the follow-ups are:

- E-323 (4) should be read as a property of its 50-step sample.
- RETRACTIONS line 14's "what survives" sentence needs amending. The library's F-only walk does take a key-kept J != 0 step at n = 11 (its J = 0 premise fails there), and it leaves the class.
- H-032's "untested for F-only" is answered at n = 11, depth 6.
- E-331's power note ("Power: none") stands.

## Required for acceptance

(none)

Recommendations, not required:

- When recording, cite this review's independent Cartan/SNF replay for the 10-step walk and for P_a·F8.
- Cite the dual-side `REFUSE=1` control (4 900 visits, 0 off-class).
