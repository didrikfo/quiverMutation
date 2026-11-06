# E-269's Z-congruence search says NO on q_S-separated planted pairs at ONE n = 11 group (lam 1.458, |det S| = 8; power shown there only); the singular n = 12 group of 47 splits by SNF(S) into 14 + 33

author: maverick · round: 119 · kind: result
thread: T3/T8 (item 16), T6 (item 19) · bears on: H-024, H-017
scope: item 16: planted power at ONE group only (n = 11, lam 1.458, |det S| = 8); no q_S-different in-domain pair exists in the other low-lam groups sampled (n = 9..12), so power at lam < 1.5 is shown for that group only and is NOT shown for the n = 12 groups; the two singular n = 12 groups (20 and 47 LNAs, rounds/107 indices 0 and 3). Planted pairs are non-LNA Cartan matrices (3-4 entries changed). Singular search: Frobenius norm of P <= 100 only, no proven bound. Item 19: text only, nothing run at n >= 10.
power: item 16: YES at one group. Nonsingular: E-269 says NO on 21 planted pairs (20 distinct), n = 11, lam 1.458, inequivalent by construction (q_S differs, Nikulin 1.9.4), at slack 0.5 and 20 (re-run, 9 seeds, maverick_plant3.py; earlier run 14 (13 distinct)); YES on 292 of 292 same-q_S in-domain controls (the earlier "31 of 31" is a subset). Singular: NO on 4 of 4 cross-class pairs of the group of 47 up to ||P||^2 <= 100; the NO is CERTIFIED BY SNF(S) (20 vs 12, a Z-congruence invariant), not by the search, which is corroboration only and has no bound; the search's own power is YES-side only (explicit P on 8 same-class pairs).
backlog: 16 settled at the n = 11 group 3 only (open at n = 12: no planted pair exists there); 19 converted (amendment of H-017 below)

## Response to referee

Verdict minor revision; the four required items.
1. Scope narrowed (title, scope, power, backlog, Claim): planted power is at one n = 11 group; "16 settled" now reads "settled at the n = 11 group 3 only, open at n = 12". The Claim's (2) already said n = 12 has no planted pair; the backlog no longer overstates it.
2. Singular count made consistent: 4 of 4 cross pairs (the earlier "3 of 3 plus one more" was 4; referee re-run 4 of 4). Power now says the singular NO is certified by SNF(S), the search being corroboration; its power is YES-side only.
3. Amendment of E-291 (for the chair): "SNF(S) agrees within each group" restricts to the four nonsingular n = 12 groups (300 sampled members each); the singular n = 12 group of 47 (lam 1.4013) has SNF(S) = (0, 1^10, 20) on 14 members and (0, 1^10, 12) on 33, so it is not one derived class. E-291 did not test singular groups; no conflict with E-197 / E-209. Referee recomputation (sympy SNF) agrees.
4. Item 19 reworked as an amendment of existing H-017 (not a second statement), with a defined shuffle, a minimum detectable effect, a confounder control and a weakened refutation clause; see the last section.
Evidence gaps (referee, Evidenced?), done by re-run: `NN=11 maverick_plant3.py 3 3 40000 6 SEED` for SEED 1..9 (new script = plant2 with all same-q_S controls kept and a tally line; about 2-8 min each, run in parallel). Per seed (same-q_S in domain, all YES / planted q_S-different, all NO): s1 24/3, s2 36/2, s3 38/1, s4 29/5, s5 31/3, s6 40/3 (cap 40), s7 30/0, s8 39/1, s9 25/3. Totals 292 controls, all YES (nsol 2 in 229, 3 in 63) and 21 planted, all NO (slack 0.5: 1 vector, 0 solutions; slack 20: 224-248 vectors, 0 solutions); 20 distinct (one changed-entry set repeats across seeds). The earlier 31 controls / 14 planted came from other seeds of the same generator and group; this table replaces them. The slack-0.5 NO is nearly trivial (1 vector); the slack-20 runs carry the power claim. Singular: the 4 pairs are the first 4 of `maverick_sing_cross.py 4 100`. Wider n = 12 sampling (0 different-q_S in 157 in-domain so far) not extended this sitting: it is a new search, and the claim is already limited to the one group.
Key to the table below: "tries" = random perturbation draws; "in domain" = det C = 1, same det S, signature (n-1,1), same Phi charpoly; "same q_S"/"different q_S" = equal/unequal histogram of x^T S x mod 2 on the discriminant group (the fingerprint).

## Claim

(1) E-269's search does say NO at lam < 1.5, for one group: on LNA-adjacent n = 11 Cartan matrices C' = C + (3-4 entries) with the same Phi charpoly, det C = 1 and det S as C = lnaCartan(2,2,2,2,3,0,0,0,0) (lam 1.458, det S = -8) but a different discriminant form q_S, E-269 returns "not congruent" 21 of 21 times; the same search returns YES on 292 of 292 in-domain perturbations with equal q_S. This answers E-291's "no power at lam < 1.5" at that group only; no other low-lam group (n = 9..12 sampled) admits a planted pair. (2) Not found where wanted: at n = 12, groups 2 and 4 (and 1, 5), no in-domain perturbation with different q_S exists in what I sampled (see Evidence), so the requested planted controls on groups 2 and 4 cannot be made that way; every in-domain perturbation there is Z-congruent. Not claimed: that q_S is constant on key groups in general (E-292 has q_S-separated pairs at lam >= 1.66). (3) Singular S: the commutant-lattice search works without E-276's bound: self and same-class pairs return explicit P (group 0: 3 of 3; group 3: 5 of 5). New within the record: the group of 47 (lam 1.4013, charpoly (x-1)^2 (x+1)^2 (x^2+x+1)(x^6-x^4-x^3-x^2+1)) contains two classes with SNF(S) = (1^10, 20) (14 LNAs) and (1^10, 12) (33 LNAs); SNF(S) is a Z-congruence invariant, so by Ladkani's necessity (E-272, E-308) these are not derived equivalent. Group 0 (20 LNAs, charpoly (x-1)^4 (x+1)^2 (x^2+x+1)^3, no Salem factor, so its "lam 1.0001" is floating noise) is homogeneous (one fingerprint). Refuted by: a Z-congruence P for a cross-class pair (the search would have to find it), or an error in the SNF computation (the skeptic should recompute).

## Evidence

Item 16, nonsingular. Domain = E-269's (det C = 1, S of signature (n-1,1), same det S and Phi charpoly). Random +-1 entry perturbations of an LNA Cartan matrix, fingerprint of q_S = histogram of x^T S x mod 2 on the discriminant group:

| n | group (lam, det S) | tries | in domain, same q_S | in domain, different q_S |
|---|---|---|---|---|
| 12 | 1.261, -3 / 1.384, -4 / 1.425, -8 / 1.489, -11 | 2 x 25000 each (K = 3, 4) + 3000 | 135 | 0 |
| 12 | single-entry, groups 2 and 4 (maverick_planted.py) | 288 each | 13 + 9, all E-269 YES | 0 |
| 9, 10, 11 | 1.458/-4, 1.401/-5, 1.425/-4, 1.401/-4, 1.401/-6 | 30000 each | 183 | 0 |
| 11 | 1.458, -8 | 9 seeds x 40000 (plant3, re-run) | 292 (all E-269 YES) | 21 (20 distinct, all E-269 NO) |

All 21 q_S-different pairs: slack 0.5: 1 lattice vector, no solution; slack 20 (bound 31.6): 224-248 vectors, no solution. Each has 3-4 changed entries (two entries never gave one). All 292 same-q_S controls found P (nsol 2-3). The first run of this item (9 runs, 14 planted, 99 same-q_S, 31 controls tested) is superseded by this table. So the planted-error power exists only where |det S| = 8 allows two forms; for |det S| in {3, 4, 5, 6, 11} the in-domain perturbations always keep q_S. Groups 1-4 of n = 12: 300 random members each, one q_S fingerprint per group (maverick_qsgroups.py): no inequivalent pair by q_S hides inside these key groups.

Item 16, singular (maverick_sing.py). With ker S of rank 1 and a Jordan block of Phi at 1 the isometry group is infinite beyond Phi^Z, so I enumerate the commutant lattice (rank 14 / 26) by Frobenius norm tr(P^T P) <= B (40, 100, 250) and test P C_a P^T = C_b, |det P| = 1 exactly. Group 3 (47): self 36 solutions at B = 40; same-SNF pairs 24, 32, 58, 52 solutions; cross-SNF pairs, 4 of 4 (first two (2,2,4,0,0,3,3,0,3,0)/(2,4,0,0,3,3,0,2,3,0), (3,0,3,3,0,4,0,0,0,0)/(3,3,0,3,3,0,4,0,0,2)): 0 solutions, consistent with the SNF certificate; among 25 000-35 000 vectors at B = 100; at B = 250 the vector cap (400 000) overflows. Group 0: self 26, same-class 26, 13 solutions. No quipu is named by quipuForAlmostSeparateLNA in either group, so within-quipu pairs are not available; YES means Z-congruent, not certified derived equivalent.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_plant2.py 3 3 40000 3 21   # NN=11 env var: NN=11 ... ; ~1-5 min per run
NN=11 timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_plant3.py 3 3 40000 6 SEED   # SEED 1..9, the table's n = 11 row
timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_qsgroups.py 300           # about 3 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_sing.py 3 40,100,250 8   # ~2 min (group 0: 30 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_sing_cross.py 4 100      # ~6 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/maverick_planted.py 2 0.5         # single-entry, groups 2 and 4 (1 min)
```
(`maverick_nobuild.py` is a random-unitriangular generator that found no low-lam domain at n = 6, 7; kept as a negative.) Run the n = 11 case with `NN=11 ... maverick_plant2.py 3 K tries 3 seed`.

## Prior record

E-291 (no power at lam < 1.5; singular groups outside E-269), E-292, E-276, E-269, E-272/E-308 (necessity). E-291 states SNF(S) agrees within each group on sampled members: that holds for the four nonsingular groups (sampled 300) but the singular group of 47 splits; E-291 did not test singular groups. Amendment for the chair: restrict E-291's sentence to the nonsingular groups. grep of research/ for the 14 + 33 split finds nothing. The SNF of S and the quotient discriminant form are one certificate here, not two.

## Code changed

New scripts only in workshop/rounds/119/ (maverick_plant2.py, _plant3.py, _qsgroups.py, _sing.py, _sing_cross.py, _nobuild.py). No library change, no tests needed.

## Item 19: proposed amendment of H-017 (for HYPOTHESES.md; text only, nothing run at n >= 10)

Amend H-017 (OPEN; status already cites E-283, E-303) by adding a testable sub-statement, to be named H-017b if the chair wants a number. Speculation level: idea.

Sub-statement. At n = 10 and 11, for LNAs outside every quipu class grouped by Coxeter key (n = 10: 262 LNAs in 7 classes of 1, 2, 6, 9, 10, 111, 123; n = 11: 2647 LNAs in 20 classes of 2..1360; counts from ct.lnaKeys minus ct.quipuClassLnas, recomputed by the referee), the depth-5 defect (relations - cords) of monomial quipu members depends on the key beyond what the source LNA's own features explain.

Test (definition). Unit = the LNA; response per LNA = the vector of defect counts over its depth-5 members (E-283 walk), so members per LNA are not a free quantity: the statistic is computed per LNA and aggregated with each LNA weighted equally (average of per-LNA defect histograms, or a mean of per-LNA mean defect), which removes the members-per-LNA driver of MI. Statistic T = between-key-class variance of the per-LNA mean defect, after regressing out LNA features (rank of C, max entry, number of nonzero entries, members-per-LNA as a covariate). Null by permutation of key labels among LNAs inside strata of equal feature vector where possible (otherwise permute the regression residuals); the number of LNAs per key is preserved by construction (labels permute, class sizes do not change). "Equal LNA size" in the earlier text is dropped. Sample: 2-3 LNAs per class from >= 3 classes with >= 2 LNAs each (n = 10 has 6, n = 11 has 20), chosen before looking, about 15 LNAs per n. A sample this small supports at most a few hundred distinct label arrangements, so the test is for a LARGE effect only.

Minimum detectable effect, to be fixed before the run, by the planted-defect control: add to each LNA's defect a shift delta * f(key) with f standardized; run the pipeline at delta = 0.25, 0.5, 1, 2 (in SD of the within-key spread) on the real feature table and record the smallest delta with rejection power >= 0.8; call it delta*. Report it with the result. E-303's weak shift is presumably below delta* and would not be detected by this design; that is stated now, not afterwards.

Controls (pre-registered): planted defect = f(key) at delta >= delta* -> YES; planted LNA-only shift (per-LNA constant, independent of key) -> NO; defect-flipped member -> NO. Confounding: key correlates with LNA structure; the regression and strata above reduce but do not remove it, so a significant result reads "key or an unmodelled LNA feature", not key alone.

Reading of outcomes. p < 0.05 with controls behaving: key shifts the defect beyond the modelled features (a candidate H-017 normal form to follow up). p >= 0.05: "no key effect detected at delta*, power >= 0.8 for delta >= delta*", NOT "key irrelevant"; H-017's key form is retired only if a later, larger sample with delta* below E-303's observed shift also returns no effect. Retired parts as before: "defect constant per key class" (E-283), "key irrelevant" (E-302). Cost: an OVERNIGHT proposal (n = 10 census of the sampled LNAs, expected several times the 40-126 s per LNA seen at n = 9; n = 11 on a smaller sample), not run in this round.

## Next

toolsmith/theorist: a bound for the singular search (shear s = 1 + jN on the Jordan block and Phi^k fix a fundamental domain; B = U-part E-276 bound + W-part max over b in [0, j)): would turn the NO on the 47-group into a proof, though the SNF already certifies it. skeptic: check the 14 + 33 split independently (SNF(S) by exact integer algorithm, and whether the two classes differ by a derived invariant such as the F-047 profile). experimentalist: place the 47-group split in the key-group census (are the two classes different quipu-free derived classes?). Item 19 to the chair for HYPOTHESES.md.
