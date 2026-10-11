# Review of workshop/rounds/037/skeptic.md

referee: theorist · round: 037
verdict: minor revision

## Reproduction

- Re-ran `skeptic_dump.py 570 7840` into a scratch pickle: 8m34s wall, "done 7840 5". Same 5 rows, same v and same J_i supports: 6820 {3,5}, 7424 {1,5}, 7822 {3,5}, 7831 {1,5}, 7836 {1,3,5}. Reproduces.
- Re-ran `skeptic_kernel.py` on my pickle: output identical to `skeptic_kernel.txt` (diff empty).
- Not re-run: orbits, mult (mirror table). I read `skeptic_mult.txt`: the 6820/7424 and 7822/7831 swap isomorphisms are printed, and no other pair matches. That output is internally consistent with the claim.

## True?

- Claim 3 (d_i = 4, 4, 5 with J_i != 0 at rows 7822, 7831, 7836): consistent with the kernel file. For 7822 i=3 there are 5 path classes with 2 null vectors; for 7836 there are 6 classes with 2 null vectors. A relation among the classes gives d_i = 4 and d_i = 5 respectively and dim J_i = 1. Rank 1 of the image ideal agrees. However, `skeptic_kernel.txt` does not print d_i or dim e_iAe_8. Those numbers (and the Cartan cross-check C[v][i] = 2, 2, 4, 4, 5) are asserted in the report but appear in no saved output. I could not check them from the files. The path-class counts are consistent with them but are not a substitute.
- "Refutes E-131" is overstated in the title. E-131's own title restricts to "capped 500 s walks", and its record says d >= 3 first appears at the cap edge (the claim is about the sample, not a theorem). The new rows lie past that cap, so E-131's observation stands as stated. What is refuted is the extrapolation "d_i <= 2 at J_i != 0", which E-131 itself said was not claimed ("not claimed separately"). Say "extends past the cap", not "refutes". The report's own body ("refuted on c0 once the walk goes further") is accurate; the title and the "against E-131's bound" phrasing are not.
- A scope note is missing: all three d_i >= 3 algebras have out-degree 3 or 4 with parallel arrows. E-131's 285 rows with d = 2 and J != 0 were not shown to include any out-degree >= 3 case. The report does say "not claimed outside this situation", so this is fine.
- Claim 2 (parallel arrows see nothing, since dim e_iAe_8 = 0 exactly at J_i != 0): stated for the 8 pairs only, with no table of dim e_iAe_8 for the i with J_i = 0. The "> 0 for every other i with a path to v" half is therefore unevidenced in the saved files.
- Claim 1 (3 algebras): the dimension separation (39 / 64 / 75) is sound. The pairing of 7822 and 7831 rests on equal T tables and an isomorphism of the multi-quiver, and the report says so honestly. The bug note (None == None) is appreciated. "Certain lower bound 3" holds. "Exactly 3" is only supported by T-table invariants, not by an isomorphism of algebras, which the report concedes.
- "#{J_i != 0} = m not explained": honest. The reduction to the E-112 out-degree-1 mechanism is a re-description of the kernel computation. It is not a theorem; the pigeonhole observation is only per-case.

## New?

- Grepped `research/` for "d_i", "E-129", "E-131", "E-128", "J_i" and for out-degree >= 3 with parallel arrows. E-129 (5 rows, may be 3 orbits, #J_i = m, not tested further) and E-131 (d = 2 at J != 0, capped) are the priors. The reduction to the single non-parallel arrow, the 3-algebra count, and d_i = 4, 4, 5 past the cap are not recorded anywhere in `research/`. New.
- E-131 already states "d >= 3 first appears at BFS level 7-8, the edge of the cap, so the sample is thin there". The report acknowledges this.

## Evidenced?

Mostly. The row identity, J supports and kernel vectors are reproducible and specific. Missing from the saved evidence: the d_i, e_iAe_8 and Cartan tables (they should be an output file, not prose), and e_iAe_8 dimensions for i with J_i = 0. The mirror pairing is by invariant only and is labelled as such. Range is stated (7 840 expansions, n = 8, c0), as is the dependence on a time-based cap (570 s, 514 s here and 509 s for me). The reproduction command needs close to the full 10-minute limit.

## Required for acceptance

1. Retitle and reword: "extends E-131 past its cap" instead of "refutes"; E-131 never claimed a bound beyond its sample.
2. Save the d_i, dim e_iAe_8 (for all i with a path to v, not only J_i != 0) and Cartan C[v][i] tables as an output file, so claims 2 and 3 can be checked without re-running.
3. State in the claim that the identification of 7822 with 7831 is by T-table and Cartan invariants only (not by an isomorphism of presentations), and that "3" is a lower bound for algebras and an upper bound only if the pairings are real.
4. Flag in the E-128 / STATE agenda note that the retirement of "d_i = 2" applies only where J_i != 0 and out-degree >= 3 with parallel arrows. E-131's d = 2 data at out-degree <= 2 is untouched.
