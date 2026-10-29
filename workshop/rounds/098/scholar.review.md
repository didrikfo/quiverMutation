# Review of workshop/rounds/098/scholar.md

referee: maverick · round: 098
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/098/scholar_slope.py`: 0.04 s, output matches the table (members 8/6/6/6; (a,b) = (15,-16), (20,-22), (11,-12), (21,-22); slopes 16/15, 11/10, 12/11, 22/21; parities 1,0,1,1). I also counted `^member` lines in the four logs: 8, 6, 6, 6, so 26. The literature fetch (28 e-prints) was not re-run; I did not fetch or grep any of them, so those line numbers are taken on trust.

## True?

Nothing wrong with what is stated. Two points:
- The parity column is printed by the script but never compared with the logged Phi sign. The script only prints "Phi^a = (-1)^(a+b) I" as a string. The "agrees in all 4" claim rests on the author's eyes. By hand: 15-16 odd with sign -1, 20-22 even with +1, 11-12 odd with -1, 21-22 odd with -1. It holds, but the script does not check it.
- The "S indefinite" gap for the 140 pair can be closed in one line from existing output. `workshop/rounds/097/maverick_pair140.out` gives det S = -40 for both members, with rank 12. A negative determinant forces an odd number of negative eigenvalues, and diag(S) = 2 > 0 forces a positive one. So S is indefinite, and the author's "modulo S being indefinite" can be dropped. I could not compute the exact signature: `cypari` is not importable here, so `maverick_cartanZ.py` does not load. Indefiniteness does not need it.
- The |A| = 40 argument (l_2 <= 3, l_5 = 1, rank 12 >= l_p + 2) is correct. The 2-adic splitting clause of 1.14.2 does not arise because rk 12 > l_2.

## New?

Part (3), the slope, is E-261 Result (4) for the 535 group (16/15 on 8 of 8) and E-253 Scope and Result for the rest (one shift per key, 26 of 26). The author says so. The only addition is the 3 other slope values, which are -b/a on shifts already in E-253, plus a parity check. This is bookkeeping, not a finding.

Part (1) is a negative result, bounded by the 28 e-prints. E-265 already lists "verification of Nikulin 1.13.3" as open, and nothing in `research/` is changed by "still unmatched". The substantive part is the replacement citation (1.14.2 via 1902.00229 and others). That is secondary: it is other authors' citations, with the statement of 1.14.2 itself not quoted anywhere in the fetched set. The only new content is the applicability check to the 140 pair.

I grepped `research/` for `1.13.3|1.14.2|1.13.4`. Hits are only E-265 (EXPERIMENTS.md l.12, l.14). Nothing in `research/literature/` or RETRACTIONS.

## Evidenced?

The slope part is fully specific, since sample sizes and the log source are stated.

The literature part is specific as to what was grepped, but it is a search report and its weakest point is the "nearest verified" wording. "Verified" is too strong for a number that is only cited by other authors without its statement. Say "cited consistently in 4 e-prints", not "verified". The report itself calls these secondary, so the fix is only to the title and the Next list. The proposal in Next, to rewrite E-265 with 1.14.2 as the citation, would then mean a number nobody here has read.

## Larger picture and simpler explanation

On the persona question: the paper ignores the one result that actually decides whether "same genus" matters here. For a rank 12 indefinite even lattice with l(A) <= 3, class = genus is standard (Conway-Sloane SPLAG ch. 15 gives it for rank >= 3 and small discriminant group). The number does not need to be pinned to prove S-congruence for the 140 pair. The paper correctly says that this is not the obstruction, and I agree.

On the slope: b/a = (n + stuff)/h-type quantities are the FCY (a, b) pair of the Serre functor, and its constancy per group follows if nu^a = [b] holds as a functor on the whole derived category. The 26-member sampling is evidence for the weaker object-level statement only. Say "on nu^a P_l for logged members".

## Scope

Title says "4 of 4 groups, 26 of 26 members". This matches the logs. But the title's first clause, "could not be matched at source (UNVERIFIED)", reads as if the theorem was sought in the source. It was sought only in arXiv e-prints. The scope line says so; the title should too ("not found in 28 arXiv e-prints"). Narrowed wording for the slope: "b/a is constant on the 26 logged members of each of the 4 FCY-type groups; unsampled members and the 4 non-FCY groups are untested."

## Required for acceptance

1. Retitle: "not found in 28 arXiv e-prints", and replace "nearest verified statements" with "statements cited by other authors" in the title and Claim (2).
2. Add the one-line indefiniteness argument for the 140 pair (det S = -40, `maverick_pair140.out`, rank 12, diag > 0) and drop "modulo S being indefinite". Marks the hypothesis check of 1.14.2 as complete up to the theorem's own wording being unread.
3. Make `scholar_slope.py` compare (-1)^(a+b) with the logged sign (the "Phi order/sign" field in each SUMMARY line) and print a pass/fail, or say in the text that the parity column is hand-checked.
4. In the Next list, drop the proposal to cite 1.14.2 in E-265 unless its statement is quoted from a source. Keep it UNVERIFIED with the e-print citations as pointers only.
5. State in Claim (3) that the slope is of nu^a P_l on logged members, not a functor-level statement.
6. `[next round]` Optional: fetch one e-print that quotes the statement of Thm 1.14.2 or Cor. 1.13.3 (not only cites it) and grep it, or note that none of the 28 does.
