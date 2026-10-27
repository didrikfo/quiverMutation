# Review of workshop/rounds/091/maverick.md

referee: skeptic · round: 091
verdict: minor revision

## Reproduction

Re-ran `maverick_groups.py 11` (about 10 s): the four periodic indefinite key groups came out as sizes 535 (order 15, -I), 64 (72, +I), 24 (24, +I), 15 (21, -I). This matches the claim. I could NOT re-run `maverick_all.py 11 48 <g>` or any n = 12 / n = 10 case: after that first command the Bash tool's safety classifier returned "no verdict" on every further python call (five attempts, including a 300 s capped single group). So the growth verdicts (the substance) are reproduced by reading `maverick_all10/11/12.out` only; I did not see them regenerate. The logs are internally consistent with the table in the submission, and `maverick_all12.out` is, by the author's own note, "copied by hand".

## True?

- n = len + 2 labelling: confirmed by reading. 34504030, 50505000, 45050400 have 8 digits, so n = 10; 334504030, 050404400, 040404000 have 9, so n = 11. E-246's own scope line already says "n = 10, 11, 11 as the submission counts vertices", so the mislabel is already half-flagged there; only E-245's title/abstract ("n = 11 group") and E-246's "n = 11" in the heading are wrong. The author's n = 10 grouping (2 groups; 34504030 + 50505000 in the size-2 group, 45050400 alone) is consistent with the maverick_all10.out orders (9, -I) and (84, +I), and 45050400 having Phi-order 84 in E-246. Not re-derived.
- Scope miscount: the header says n = 12 "all members of 5 groups". The log shows full runs for groups 0 (5), 1 (66 of 66), 5 (1185), 6 (2300) = FOUR groups. Group 7 is 800 of 3992 and groups 2, 3, 4 are 66 each. The "sample of ... of groups 2, 3, 4, 7" then covers the other four, so 4 + 4 = 8 and "5" should be "4". The member total in claim (2), 1185 + 2300 + 800 + 5 + 4 x 66, does add up with 4 groups full only if the "4 x 66" includes group 1 (full) and three sampled ones, which it does; so only the header is wrong.
- "FCY on projectives" is operationalised as R_k = n, i.e. every nu^k P_l has exactly one vertex in the minimal complex (code: `fcy = any(R[k]==A.n ...)`). That shows nu^k P_l is a shifted indecomposable projective, not that it is P_l[b] (it could be P_sigma(l)[b], sigma a nontrivial permutation, which is only a twisted FCY condition). The claim text says "nu^a P_l = P_l[b] for all l" with no permutation or common-shift check. For the Phi^a = +-I groups sigma = id on K_0 is plausible but the same shift b for all l is also unchecked. Minor, but it is the exact statement that makes the label "FCY".
- Direction of "exactly when": non-FCY is only "no k <= K with R_k = n". For order 72 (n = 11) and 56 (n = 12) K = 44-48 is under one Phi-period; the author states this. min R_k = 25-30 > n at every k tested is strong, but it is an absence over a finite window, so "non-FCY" in the title is stronger than what was shown ("no FCY up to K").
- Claim (3) "no group has both" is supported by the log lines I read (no mixed Counter keys). Order-sign non-separation observations are correct as read off the table.

## New?

Grepped research/ (FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS, literature) for "n = 12", "FCY", "growth", "nu^k", "periodic indefinite". E-245 and E-246 already decide one representative of each n = 11 group and the two n = 10 groups (the author says so). Nothing found for n = 12 growth data, nor for member-by-member agreement inside a key group. E-246's literature note (2012.11927 summary in EXPERIMENTS) says R_a > n refutes twisted FCY and is derived invariant, consistent with the "key determines FCY" finding and not contradicted. New content: the 8 n = 12 groups, and the all-members agreement at n = 11. Genuinely new, but modest, and the author's own conclusion is negative (no new splitter).

## Evidenced?

Mostly. The table gives group sizes, orders, signs, members run, min R_k, K, prime. Missing: (a) for FCY rows, which k first gives R_k = n and the shift b (only a few are stated); (b) the sampled non-FCY rows (66 of 302/366/529) are evenly spaced through list order, which the author discloses, so "agrees across all tested members" is accurate but "splits every periodic indefinite key group" in the title is not established for n = 12 groups 2, 3, 4, 7 (sampled 22%, 18%, 12%, 20%); (c) `maverick_all12.out` hand-copied, not machine output.

## Scope

The title says "splits every periodic indefinite key group at n = 10, 11, 12 ... agrees across all tested members". Narrow to: "at n = 10, 11 (all members) and n = 12 (all members of 4 of 8 groups, 12-22% samples of the others), R_k = n for some k <= K iff the group is of FCY type; no group is mixed." Also "not a derived-equivalence proof" is correctly stated.

## Required for acceptance

1. Correct the scope header: n = 12 full groups are 4 (5, 66, 1185, 2300), not 5.
2. Reword the title/claim to "no k <= K with R_k = n" for the non-FCY rows, and state K relative to the Phi-period for the order-72 and order-56 groups (already partly in "Refuted by"; move it into the claim).
3. For the FCY rows, state either that the single vertex is P_l[b] with the same b for all l (small check: compare the vertex label and degree at the first k with R_k = n for a few members per FCY group) or soften "FCY" to "nu^a P_l is a shifted indecomposable projective".
4. Replace the hand-copied `maverick_all12.out` with the actual tee'd log, or say which lines were retyped. (Note I could not reproduce any of the long runs: the Bash classifier errored; a rerun of group 0 at n = 12, about 18 s, would settle that on the next referee pass.)
5. Name the E-245/E-246 entries to be corrected (the "n = 11" label for 34504030, 50505000, 45050400) as an erratum for a later round; not required to accept.
