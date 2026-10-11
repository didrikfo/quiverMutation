# Review of workshop/rounds/058/skeptic.md

referee: theorist · round: 058
verdict: minor revision

## Reproduction

Ran `CLS=1 timeout 10m .venv/bin/python workshop/rounds/058/skeptic_wrongalg.py`, 7 s. The output is identical to the claim: 280 (eqdim, truth=NOT, test=NO), 168 (eqdim, truth=iso, test=iso), total 448, and 47 of 47 drop-one cases "dim > End(T), test=iso". The seven steps and their killed points match. I did NOT rebuild `/tmp/tsm/c1.pkl` (about 10 min of collection); I used the pickle already on disk (timestamp 01:00-01:02, apparently from this round). The 20 000-expansion provenance of the pickle is therefore unverified by me.

## True?

I found no error. I checked the following.
- Ground truth: a label-preserving iso must act as scalars on the x-arrows and as GL2 on the parallel pair, so it is a Moebius map on the killed points. For at most 3 ordered points on P^1, the PGL2 orbits are exactly the coincidence patterns. The argument is sound. It holds only with 3 relations; the author states this.
- Reading `symcheck2`: it compares arrow counts per pair only, then solves for a surjection End(T)-from-child (matrices plus radical-square corrections). It never compares dim K Q/I' with dim End(T). The "dimension is an unchecked input" observation is correct.
- Gap in the 47 drops: the quotient is bigger than End(T) and `symcheck2` says "iso". That is consistent with surjection-only, and it shows "iso" is not a statement about the algebra K Q/I' unless `crels` is complete. The author says so. The count 47 covers all 8 parallel-arrow steps (step 12 included), while the main table covers 7. The text blurs this ("47 variants over the parallel-arrow steps"). State it.
- The "NO" on 280 might in principle come from a cause other than non-isomorphism, such as a solver or timeout artefact. The script maps any verdict not starting with "iso" to NO. Check that all 280 are genuine "NO (no solution)" and not an error or timeout string.

## New?

`grep E-174 research/`: E-174 (EXPERIMENTS.md:9, :16) lists "equal-dims wrong-algebra control with a parallel pair" as open. This is that control, so it is new as executed. E-174's 9/9 and 9/9 perturbations are not equal-dim-verified, as the author says. Nothing found in RETRACTIONS or HYPOTHESES for "equal-dim", "wrong algebra".

## Evidenced?

Mostly. The 2x2 table, the variant count (4^3 x 7 = 448), the independent dimension computation, and the reproduction are specific enough. Weak points:
- The family is a single type of error (killed lines) with no moduli. The author admits this. It shows the test separates coincidence patterns, not that it detects a continuous parameter. The claim "has power against relation errors of this kind" is fine; the title's "not vacuous on the relations" is slightly stronger than a one-family test supports.
- The 168 "iso" acceptances show only that the solver finds a solution. It does not show the solution is an isomorphism, since dim equality is not checked. In this family it is by construction (equal dim, checked externally), so this is fine here but not in general.

## Scope

Matches: n = 7, class 1, 7 of 8 parallel-arrow steps, label-preserving, equal-dim, killed-line family. The title says "280 of 280 ... parallel-arrow children". Narrow to "7 steps, one parallel pair, 3 line-type relations".

## Required for acceptance

1. Report the exact verdict strings for the 280 NO cases (count by string), to rule out error or timeout artefacts. One sitting.
2. Say in the Evidence that the 47 drops include step 12 while the 448 table excludes it.
3. Soften the title: "rejects 280 of 280 equal-dimension killed-line variants" rather than "non-isomorphic parallel-arrow children" in general.
4. Record that the pickle was not rebuilt in this review, or give its hash and expansions count. [optional]
5. Step 12 and class 2 stay `[next round]`, as the author lists.
