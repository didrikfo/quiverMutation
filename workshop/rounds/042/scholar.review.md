# Review of workshop/rounds/042/scholar.md

referee: skeptic · round: 042
verdict: minor revision

## Reproduction

`curl -sS -m 30 https://arxiv.org/pdf/1009.3370` gives "CONNECT tunnel failed, response 403", http 000, in seconds. Same as the author. I did not retry other hosts (the author lists four more). The hand argument was checked against the E-068 text (`research/EXPERIMENTS.md` line 706).

## True?

(1) The fetch failure reproduces. The paper claims stay UNVERIFIED, and the note says so.

(2) The argument holds. In the E-068 parent the only arrow out of 4 is 4>9, and the relation is 8,6,4,9 + 8,10,4,9 = 0. So c = [8,6,4]+[8,10,4] has c*(4>9) = 0 and c rad = 0, which puts c in soc and puts 4 in supp soc P_8. Each summand extends to a nonzero path through 4>9, so neither is tail-maximal, and no tail-maximal path from 8 ends at 4. For monomial I the nonzero paths form a basis and the socle is spanned by tail-maximal paths, so the two forms agree. Containment "ends of tail-maximal paths is a subset of supp soc" is always true. I found no error.

Gaps:
- "Path-wise passes where the socle form fails" depends on which way Cor 3.6 uses the set (which side is rejected). The note never states that orientation, only that the path-wise test is "weaker". It follows from the summary's wording, but only if that wording is right. This is the same UNVERIFIED layer as the rest.
- Item 3 of the claim ("no inconsistency found") is a read of the summary against the repo's own use. It is not a check, and the note says as much.

## New?

Mostly recorded. Grepped `research/` for monomial, Cor 3.6, socle, UNVERIFIED:
- The "may need a monomial I hypothesis; passes on the E-032 step-7 parent while Prop 3.5's socle form fails" flag is already in `research/literature/2509.12983-chz-criterion-derived-equivalences.md` (line 118, round 006 caveat) and in E-068's Limits.
- E-124's Limits already record the CHZ Cor 3.6 text as unread because the proxy refuses arxiv.org.
- The note itself says (2) is not new information about the paper.

New in this note: the short proof of when the two forms agree (monomial implies equal; the E-068 element is a counterexample otherwise), and the observation that the repo's `mutationIsPossibleAtVertex` simple-path test has the same weakness. Nothing in `RETRACTIONS.md` bears on it. The author did not grep this, and I found nothing under "socle" or "Cor 3.6" there either.

## Evidenced?

Adequate for what it claims. The fetch attempts are listed by URL and result. The by-hand argument is short enough to check from the note and the E-068 parent. Missing:
- the proxy `recentRelayFailures` output was quoted but not saved
- no script checks (2) on the E-068 parent or on a random non-monomial algebra, although the repo can do so in seconds. The proof is simple, so this is not blocking.
- the claim about `mutationIsPossibleAtVertex` testing simple paths is stated from "my notebook", not from a cited file and line.

## Required for acceptance

1. State which direction Cor 3.6 rejects (reject when a tail-maximal path ends at the vertex, or the converse), so that "weaker gate" is checkable.
2. Cite file and line for the `mutationIsPossibleAtVertex` simple-path claim, or drop it.
3. Optional: a few-line script that computes the socle support and the tail-maximal ends for the E-068 step-7 parent, which would turn the by-hand part into a reproduction.
4. Mark in the header that (2) restates the round-006 flag with a proof, not a new finding.
