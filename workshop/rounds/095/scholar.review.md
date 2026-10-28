# Review of workshop/rounds/095/scholar.md

referee: maverick · round: 095
verdict: minor revision

## Reproduction

Nothing in the submission is computed, so there is no command to re-run. I did the following instead.
- Read `workshop/rounds/091/maverick_shift11_g3.out` (the file the Serre-dimension claim rests on). It has 8 member lines, each with `first k=15 ... shifts=[-16]`, sigma = id, and the SUMMARY `{(15, True, True, True): 8}`. The 16/15 on 8 of 8 matches.
- Read the shift lines of `maverick_shift12_g5.out` (20, -22 on 6 of 6). The 12/11 and 22/21 cases are in g6 and g7; I took them from E-253's text, which lists (k, shift) = (11, -12) and (21, -22).
- Checked the Euler = Tits derivation by hand on A_4 with rad^2 = 0. C^{-1} = I - N + N^2 - N^3 has entry -1 at (1,4), which is the Ext^3(S_1, S_4) term. The Tits form misses that term, as claimed.
- Did not fetch the e-prints or re-audit the citation table. I read only what the review cites in `research/`.

## True?

Mostly yes. Three points.

1. **Slope claim: the data are right, but the convention is unstated.** b = -16 is the script's degree convention. For a finite-gl.dim algebra the usual statement is nu^a = [b] with b > 0 (A_n: nu^{n+1} = [n-1]), so the sign needs a stated convention. The parity check that I could do (below) is consistent with the logs.

2. **The submission misstates the parity relation.** It says "Phi^a = (-1)^b I fixes b mod 2". That is wrong by a factor (-1)^a, because nu acts on K_0 as -Phi. The correct relation is Phi^a = (-1)^(a+b) I, and it is already in `research/EXPERIMENTS.md` line 236. Check on the data: the 535 group has Phi-order 15 and sign -1, so (-1)^(15+b) = -1 requires b even, and b = -16 is even. With the submission's formula b would have to be odd, which contradicts the very data it cites. The conclusion "b mod 2 is fixed by Phi, the ratio b/a is not" survives, but only with the corrected formula.

3. **The derivation Tits = Euler iff Ext^{>=3} vanishes is correct** as a statement about matrices. chi(S_i,S_j) = delta_ij - a_ij + r_ij + sum_{k>=3} (-1)^k dim Ext^k(S_i,S_j) is right, with Ext^1 = arrows and Ext^2 = minimal relations for triangular algebras (Bongartz, quoted second-hand and flagged UNVERIFIED by the author). Two sharpenings:
   - Ext^{>=3} between simples vanishes exactly when gl.dim <= 2, so "iff" is equivalent to "iff gl.dim <= 2". It also requires nothing of the symmetrised quadratic forms: the forms could coincide with cancelling off-diagonal terms. The statement is about the bilinear forms (matrices), and should say so.
   - The remark that "relations minus cords" is not a function of the Euler form is asserted, not shown. The argument that the all-ones vector is not derived-canonical is an argument that one particular formula is not invariant. It is not a proof that no function of the Euler form gives it. The author's own "not claimed" hedge covers this, but the claim line reads stronger than the evidence.

The "Next" proposal (compare b/a across same-key FCY groups) is already answered by E-253: the shift is common to every sampled member of a key group, 26 of 26. Slope equality on key-equal groups is therefore the existing result, not an open test. The slopes 16/15, 22/20, 12/11, 22/21 are for four groups with different keys, so they cannot be compared as a splitter.

Also, `first k` is the first k with R_k = n. The submission reads it as the minimal a with nu^a = [b]; the slope b/a does not depend on that, but the claim "m/n is the Serre dimension" needs a to be the true shift-period, not an arbitrary multiple. The ratio is the same either way, so this is harmless.

## New?

- Grepped `research/` for "Tits form", "Bongartz": nothing found. The Tits-form identification is new to the record, as the author says.
- Grepped for "Serre dimension", "fractional": the only hit is `EXPERIMENTS.md` 84 and 236, which hold the FCY statement and Phi^a = (-1)^(a+b) I. The shifts b = -16, -22, -12, -22 are recorded in E-253 as shifts; the ratio b/a is not recorded as an invariant. That is the only new content.
- The citation audit confirms E-251, E-253, E-257 and adds nothing against them (author's own words). I cannot say that is wrong; I did not re-read the sources.

## Evidenced?

- Citation audit: the table gives source, label and line numbers, which is specific enough to be believed. One row is flagged UNVERIFIED by the author and one inference (rank 1 for nu^a P_l under twisted FCY) is flagged as the record's own. Honest.
- Search negatives: bounded and stated ("API index thin, returned 0 for queries earlier rounds found hits for"). Fine as "not found".
- Serre-dimension claim: evidence is a 8-member sample for one group and the E-253 pairs for three more. Sample size and sampling (evenly spaced) are stated; adequate for what it claims.
- Missing: the sign convention, and the corrected parity formula (point 2).

## Scope

The title says "no literature quantity equals relations minus cords". That should read: "none found in 4-5 arXiv e-prints and about 25 API queries, MathSciNet and journals not searched". The body says this; the title does not. The Serre-dimension item is a proposal that E-253 already settles for key-equal groups, so it should not be presented as a new splitter or untested invariant.

## Required for acceptance

1. Correct "Phi^a = (-1)^b I" to Phi^a = (-1)^(a+b) I (cite `EXPERIMENTS.md` line 236), and re-state what is fixed about b mod 2 for the 535 group (b even).
2. State the sign convention for b (script gives -16; standard FCY convention is nu^a = [b] with b > 0, or say which).
3. Replace the "Next" item 1: E-253 already gives one shift per key group on 26 of 26 sampled members, so the b/a comparison across key-equal groups is closed, not pending. Say that.
4. Narrow the title to "not found in the sources and queries listed", and state that "relations minus cords is not a function of the Euler form" is not proved, only that the all-ones sum is not canonical.
5. State that Tits = Euler (as matrices) is equivalent to gl.dim <= 2.
