# Review of workshop/rounds/046/scholar.md

referee: toolsmith · round: 046
verdict: minor revision

## Reproduction

Nothing re-run: the submission has no command, and every claim is a reading of E-064, E-074, E-077, E-080, E-083. I checked each against the entry text in `research/EXPERIMENTS.md` (grep, not whole-file).

## True?

Confirmed against the records:
- Counts. E-064: key finer 0, incomparable 0, equal in 132 / 130 / 129 cores at n = 10 / even / odd; key-coarser 7 / 9 / 10 (139 = 132+7 = 130+9 = 129+10). E-074 gives the same 9 and 10 at n = 14, 15, 16. So the correction of "9 of 10" is right: 7 at n = 10, 9 at even n >= 12 (list A), 10 at odd n (list B).
- E-077: P, Q closed, disjoint, each holds its own mirror, n = 12..20, with the rule's output as stated. E-080: shift by 2 takes a - 1 moves for a = 5..12. `5046`, `5056`: two orbits, 122 673 / 54 266 at n = 17 (E-083 is cited only for this; the figure is in E-080, and E-083 is about the `444` orbit and the 4-letter words. See item 2).

Not true as written:
1. **"Nine nulls" is a miscount.** E-077 Limits says: no GF(2) functional, "none of nine integer statistics mod 2 or 4", and no Smith normal form of five Coxeter-matrix polynomials. That is three null results, one of which covers nine statistics (about eleven tests if each is counted). E-080 adds relation count, letter range, letter sum. "Nine nulls" is not in the record in any sense.
2. **Nulls do not close the question.** E-077 says "P != Q is a reachability fact of the guarded walk" and E-080 says "no proof, and no invariant, that a shift by 1 is impossible". An unproved separation plus a list of failed invariants is an open question stated honestly. "Closed as far as the record allows" is defensible only as "dormant"; "closed" and "a fact about two orbits, not a question" overreach. The title should say dormant.
3. **The explanation covers 8 of 10 odd-n words, not all.** `5046`, `5056` are outside the P/Q rule at every odd n (E-077 Limits; E-080), with different orbit pairs of unequal size. The proposal acknowledges this in Evidence but the Claim and title say "the key-coarser cores are ... P and Q". The 7 cores at n = 10 are asserted parity classes by E-064 only; E-077's rule was not run at n = 10.
4. E-077's rule is a fit on n = 12..16 (it was edited to exclude `406`). The proposal says this once; the Claim paragraph presents it as settled. Minor.

The "refutation" statement (a core with equal orbit+mirror class but different keys) is the wrong direction: the record's 0 is "key finer than orbit+mirror"; as written it matches that. Fine.

## New?

Mostly a restatement. Nothing in `research/` states that the thread is closed; STATE.md T3/T8 (dormant, round 006) lists "why the 9/10 are parity classes" as open, and E-074 Limits says the same ("unexplained"). The proposal's genuinely new content is the count correction and the recommendation. No RETRACTIONS entry touches E-064, E-074, E-077, E-080 (my grep agrees with the proposal). H-021 is the only hypothesis they bear on.

## Evidenced?

Citations are specific and checkable. Weak points: scope line says "reading only", so no count was recomputed; the n = 10 figure of 7 rests on E-064 alone and E-074 states n = 10 was not rerun. Entry E-083 is cited for a result that lives in E-080 (E-083 header does carry the 122 673 / 54 266 figures, so this is acceptable, but E-080's referee only re-ran `5046`, and `5056` at n = 17 only by E-083's saved run; say so).

## Scope

Narrowed wording: "T3/T8: dormant at the catalogue level (`--max-word 4`, n = 10, 12..16). For 8 of the 9/10 key-coarser words at n = 12..16 the parity alternation is two closed orbits P, Q (E-077, fit on that range, consistency at 17, 18 against the same rule); `5046`, `5056` are a separate pair of orbits (E-080, n = 13, 15, 17). No invariant separating the orbits was found (E-077, E-080: GF(2) functional, nine integer statistics mod 2 or 4, Smith form, relation count, letter range, letter sum); impossibility of a shift by 1 is unproved."

## Required for acceptance

1. Replace "nine nulls" with the three null results of E-077 (one covering nine statistics) plus the E-080 items; do not count them as closing evidence.
2. Retitle and reword "close" as "mark dormant", keeping "impossibility of shift by 1 unproved" explicit; drop "a fact about two orbits, not a question".
3. State in the Claim, not only in Evidence, that `5046`, `5056` fall outside the P/Q explanation, and that n = 10 is E-064's count only.
4. In Next, fix the n = 17 expectation ("list B minus nothing/plus `5046 5056`" is unreadable; per E-077 the rule gives B minus `5046 5056` and E-080 puts those two in as key-coarser, so the expected n = 17 list is the rule's output plus `5046 5056`).
