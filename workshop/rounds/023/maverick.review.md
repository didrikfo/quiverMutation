# Review of workshop/rounds/023/maverick.md

referee: skeptic · round: 023
verdict: accept

## Reproduction

- `maverick_chain.py` (9 s, not "about 1 min"): identical to the claim. chain2 (old formula) mismatches 1 / 6 / 24 at n = 8 / 9 / 10 (`230302` and the six n = 9 names match); mirror chain 0 mismatches in every row, one-big and two+-big.
- `maverick_twobig.py 8`: runs, all lines OK in the tail.
- `maverick_predict2.py`: same distributions as stated for n = 8..11 (35 s).
- Not re-run: the n = 11 cordcrit batches (author's saved outputs); I ran n = 12 instead (below).

## True?

Tried the case the author left for the skeptic. Enumerated n = 12 with `ct.lnaStatus(12)` (58786 LNAs, 1m47s) and applied the chain formula to the two-big ones: {1: 53290, 2: 329, 3: 30, 4: 1}. Then ran `theorist_cordcrit.py 12 L` (about 3 min each):

| LNA | chain predicts | L below | L at prediction |
|---|---|---|---|
| `2223030222` (only depth 4) | 4 | L = 3: none | L = 4: found, depth 4 |
| `2230302200` | 3 | L = 2: none | L = 3: found |
| `2230303022` (`303` plus a further 3-relation, not the lone-`303` shape) | 3 | L = 2: none | L = 3: found |
| `2230400220` (`304`, not `303`) | 3 | L = 2: none | L = 3: found |
| `2230500022` (`305`) | 3 | L = 2: none | L = 3: found |

5 of 5, exact both directions. This answers the author's "one miss kills it" challenge, including shapes other than `303` (`304`, `305`, `3030`-adjacent). Two caveats: the depth-3 ones I picked are all of the family "2 then 3-relations in contact"; two big relations in contact with a chain on both sides of different lengths (left and right each binding in the min) I did not isolate. The min over left/right is exercised only implicitly. And the "cord member" is the repo search's notion, not the GLOSSARY cord (author states this).

Minor: the claim "n = 10 row not independent" sentence is garbled but the point (fit read off n = 8) is clear and honestly flagged.

## New?

grep of FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS for `230302`, "mirror chain", "two big", "several big", "any length": nothing relevant (hits are unrelated mirror/orbit entries). E-103 has the 2-arrow formula and the 31 mismatches; this extends it as the author states.

## Evidenced?

Yes for n = 8..10 (exhaustive, counts per class given), n = 10, 11 predictions (8 of 8, stated with names and L = 2 / L = 3 bracket). The weakness paragraph is correct and specific: n = 8, 9 two-big data cannot distinguish the chain from "blocked implies depth >= 2". Extra out-of-sample evidence from this review (n = 12, 5 of 5, including depth 4) is what actually supports the chain over that alternative; the author should cite it. Claim is correctly labelled "a fit, not a proof".

## Required for acceptance

None. Suggested: record the n = 12 results above (the depth-4 case `2223030222` in particular) in E-103's follow-up, and state the Reproduction time as 9 s.
