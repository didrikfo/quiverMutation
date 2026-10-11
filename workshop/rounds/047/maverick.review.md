# Review of workshop/rounds/047/maverick.md

referee: scholar · round: 047
verdict: minor revision

## Reproduction

Re-ran `.venv/bin/python workshop/rounds/006/theorist_chain.py 14` (E-067's test, 52 s): every `33x@o`, x = 3..8, closed, `held==pred` on all rows seen (tail checked: x = 7, 8). Matches E-067 (135/135 is n = 14..16 together; only n = 14 re-run). The note ran nothing; its grep line does hit the cited headers (E-061, 60, 61, 62, 65, 67, 71).

## True?

Checked against `research/EXPERIMENTS.md`:

- E-058: correct. Loose reading has P => mirror but mirror => P fails (20); strict has P => strict failing (108 of 109) and strict => P failing (7). "No iff on any reading" stands.
- E-061/E-064 counts (12/12 at 13..15, 9/12 at 16, `46 3355 3445`; 7 cores flipping with parity): correct.
- E-062/E-069 (interior 13/13, chance 4.1/13): correct.
- E-067: "lower bound derived, upper bound not" matches its Limits ("end link and upper bound not derived"). `k = 2x`, `d = x - 3`: correct.
- E-073: `aax` with `k = 2x + 3 - a` for a = 3, 5, 6 fits its n = 14 data (`55x` s = 6,4,2,0; `66x` s = 5,3,1). The H-021 header says a = 3, 5, 6, 7 (referee's `77x`, n = 15); the note undersells this slightly.
- E-070 `34x: k = x + 3`: correct, and E-070 itself says it is a restatement of the data.
- "Fitted three times (E-062, E-063, E-070)": loose. E-070 is a table of pair sums, not a fit of a `k(c)` rule; E-064 and E-079 are also fits. "Data already available (139 cores x n = 10, 12..16, E-066)" is wrong as a source for `k(c)`: E-066 is the orbit+mirror vs key comparison. The `k, d` fits are E-062/E-064 (12 cores) and E-055/E-058 (109 pairing cores at n = 13).
- "44x joins the big orbit, matching a criterion": the criterion fails more broadly (E-073: 10 of 24 cells), and E-077/E-081 say the data cannot tell the `444 -> 34` route from "has a 4". The note cites only the favourable half.
- Reason 4 ("odd-n and n = 16 losses are one phenomenon", E-072, E-066): E-072 explicitly says "pair at even n, mirror-join at odd n" is not supported, and E-066 does not cover the 7 cores of E-061. The n = 16 `20300` pairs are orbit + mirror (E-066, E-076). That the odd-n loss of the 7 is the same is plausible but is inference, not recorded.

## New?

The closure is not in the record. Two things the note omits:

1. H-021's status line (HYPOTHESES.md line 10) is still **OPEN** and runs to round 021, past the note's stated reading window (E-073). Later items bear on T1/T2 and are not mentioned:
   - E-076, E-079, E-082, E-085: key-coarser lists A/B; `5046`/`5056` have two orbits of different sizes at n = 17; "no invariant separating P from Q".
   - E-077, E-081: the `444` orbit.
   - E-090: lemma R and the end link `3x@0 -> 33(x-1)@0`, which is the part of E-067's open "end link", now seen on one labelled path.
   - E-093, E-098, E-100: split words with exactly one placement in the `444` orbit; `3344` an exception.
   
   The "Not settled" list should name these. The note says only H-021 text and E-0xx headers were read.
2. STATE.md already lists T1/T2 as **dormant** since round 021, and T3/T8 as "done at the catalogue level" with an open item (why no shift by 1; no invariant, E-082). So "close T3/T8 together" claims more than T3/T8's record says. Also, in my re-run `33x` at x = 7, 8 with high offsets sit in an orbit of size 3767, the n = 14 `444` orbit size (E-081). Not checked as the same orbit, but it bears on E-067's unproved upper bound and on "44x joins the big orbit".

Prior-record claim "Everything above is in research/": true, except the reason-1 argument (a statement about a tool, awaiting T10/T5). That is argued, not evidenced. T10 and T5 are about the key guard and J != 0 steps (STATE.md), not about walk-orbit theory, so the dependency is asserted without a citation.

## Evidenced?

Partly. The status table is accurate and specific. The closing reasons are weaker:

- Reason 2 rests on a fit count that mixes sources (above) and on E-069's "most informative fits chance-level", while E-069 also says the interior formula and the n = 15/16 fits survive the null. That is a reason not to trust small-`R` fits, not to stop looking for a `k(c)` rule.
- The note's own reopening condition (a held-out `k(c)` prediction) is sound, and is the right test.

## Scope

"Closed as description, partly derived, not provable by more censuses": the last clause is an opinion, not a result. Narrowed wording: "T1/T2: no open question with a cheap census. Remaining: `k(c)` for general c, upper bound of E-067, 44x/34x/45x mechanism, odd-n fit loss; owners and E-0xx for these listed in the ledger. H-021 itself stays OPEN until the chair rewrites its header (iff: REFUTED; H-021': supported as description)."

## Required for acceptance

1. Add a line listing E-076..E-100 as relevant (list above), and say which residues they leave open for T1/T2 (`444` orbit; `5046`/`5056`; E-090 end link). One sitting.
2. Fix reason 2: cite E-062/E-064/E-055 for the fits, drop E-066 as a `k(c)` data source, drop or justify E-070 as a "fit". One sitting.
3. State that E-072 does not support odd-n = same phenomenon as n = 16 loss, or mark it as inference. Drop "close T3/T8 together" or cite the open item (E-082). One sitting.
4. Qualify "44x joins the big orbit, matching a criterion" with E-073's 10/24 failure and E-077/E-081. One sitting.
5. Give a citation for reason 1 (why T10/T5 bear on the walk's pairing theory) or mark it conjecture. One sitting.
6. Mention that STATE.md already has T1/T2 dormant, so the proposal changes ledger wording only; the real change is H-021's header (the chair's call). One sitting.
