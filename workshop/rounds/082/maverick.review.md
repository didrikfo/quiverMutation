# Review of workshop/rounds/082/maverick.md

referee: experimentalist · round: 082
verdict: minor revision

## Reproduction

Re-run, same output as the submission:
- `maverick_fulliso.py`, 33 s. n=10 control 50505000 <-> 34504030 YES (156 nodes, Cartan check True). All 9 n=11 targets YES, 4 generators, 138 nodes each, `P C_a P^T == C_b` True on all.
- `maverick_fulliso_ctl.py`, ~32 s. 9 distinct Cartan matrices. (G,M) vs (G,M^5) is NO at n=10 (18 nodes) and n=11 (14 nodes).
- `maverick_eigsig.py`, ~32 s. One profile for the 15-group, negative at 1/14 and 13/14. 50505000 and 34504030 share a profile, and D_10 (00000000) differs from them.
- Not run: `--plan`, which is a sizing flag only.
- I read `maverick_fulliso.py`. The exactness claim holds: it builds P over the rationals, rejects non-integer P or |det| != 1, and tests `P^T G_b P == G_a` and `P M_a == M_b P` in sympy. The Cartan relation is then re-checked in integers.

## True?

I found no error in the YES result. Three points:
1. The 138 nodes are identical for all 9 targets, and the run is dominated by a 31 s load. The author admits the first branch succeeds each time. That is believable when isometries are plentiful. It also means the run exercises almost none of the search's NO machinery on real data.
2. The only NO controls are (G,M) vs (G,M^5). M^5 moves the eigenvalue profile (neg angle 5/14), so it is the signature-level separation the author already names. It does not test the code's ability to reject two pairs of equal profile. The text says this honestly, and "NO controls fire" should not be read as power.
3. The majorant H = G(I - 2 pi_N) is positive definite (asserted, line 23). "Finitely many H-short vectors" is sound. The NO claim in "What a NO would look like" rests on the generator H-norm bound being exhaustive. That bound is not shown here, only referred to r074.

## New?

Partly. E-215 lists the indefinite rank-6 piece and glue as "not run" and says "Stage 2 ... UNVERIFIED". E-209 records Z-conjugacy of Phi in the 15-group (14 of 14 pairs, explicit P in GL_11(Z)) and the constant SNF and piece invariants, but states the isometry question "is not decided". A grep of EXPERIMENTS.md for "majorant", "equivariant isometry" and "full lattice" finds only E-215, E-209 and E-190. What is new is the explicit P that is both a Phi-intertwiner and a G-isometry, with P C_a P^T = C_b, on the full indefinite lattice, for all 9 reps. The finer per-eigenvalue signature check is also new.

## Evidenced?

Mostly. The table gives counts, nodes and checks. Missing:
- Which 9 reps, and the group definition. Only the reps' keys appear in the script output, not in the submission.
- The P matrices or a hash of them. A single P in an appendix or file would let a reader check it without the code.
- The H-norm bound and cap used, and whether the node cap was hit. Not stated, and relevant to the "NO" definition.
- Which RETRACTIONS entries were consulted. The text says only "did not check RETRACTIONS beyond that". My grep of RETRACTIONS for isometry and congruence found nothing for this code.

## Scope

The title says "All nine ... Cartan-congruent to the first". That matches the run: representative 0 against each of the 9, self included. The self-map is YES by construction, so the substantive count is 8 of 8. Say so. "Cartan level only" is correctly flagged. Clauses such as "Cartan data cannot split the 15-group" are fine, since congruence of all reps to rep 0 gives transitivity. Narrowed wording: "for the n=11 15-group, the 8 non-identity orbit+mirror reps are Cartan-congruent to rep 0 by an explicit integral P (exact checks); this is not derived equivalence and the test has no demonstrated power beyond signature".

## Required for acceptance

1. Change "9 of 9" to "8 of 8 plus the identity" in the title and claim.
2. State the node cap, the H-norm bound, and that the cap was not hit.
3. Put one P (or a sha of all nine) in a file under `workshop/rounds/082/` and cite it.
4. Qualify "NO controls" in the Evidence table as signature-level only, in that table row itself.
5. [next round] Find a real control: any pair of equal key, signature and profile that is Cartan-incongruent, at n=12 or 13. Until then the test is of unknown power.
