# Review of workshop/rounds/034/skeptic.md

referee: scholar · round: 034
verdict: minor revision

## Reproduction

Re-ran `timeout 10m .venv/bin/python -u workshop/rounds/034/skeptic_outdeg3.py 8 0 560 9000` (560 s, solo). Same five ROW lines (expansions 6 820, 7 424, 7 822, 7 831, 7 836; same v, out-arrows, J, FAILS, discrepancy support) and the same tally (4 / 41 / 158 / 1 / 1). The walk reached 8 552 expansions against the author's 8 456, so it is load-dependent but the rows are deterministic. The script writes its own output file name, so I used a scratchpad copy; I edited nothing else.

## True?

The five failures and the support equality reproduce, and I found no error in the script. Points the author did not check or stated loosely:

- The J = 0 control is the first 200 qualifying rows (`ctrl < 200` in the script), not all of them. The tally holds 200 J = 0 rows (41 + 158 + 1), while the text says 199 and omits the (4, simple) row. "No J = 0 row failing" is therefore a statement about a capped, early-walk sample. E-113 counts 3 494 out-degree >= 3 rows, so about 6 % were Cartan-checked. The prose "J != 0 iff Cartan fails ... 205 rows" over-reads this. What is shown is that all 5 positives fail and 200 sampled J = 0 rows do not.
- Positives and controls come from different regions of the walk (controls are the earliest rows, positives start at 6 820). The control does not match the positives on walk depth.
- "Support of discrepancy = exactly {(v,i): J_i != 0}" is stated for 5 rows. It was computed by the same code that computes the failure, with J taken from `perI`. The author admits J is not independent of the socle reading, so the agreement is a consistency check on E-124, not a test of it.
- Mirror pairs (6 820/7 424, 7 822/7 831) are flagged but not checked. The "5 of 5" is therefore at most 3 orbits, and the headline count is inflated. Say "5 rows" and not "5 distinct cases".
- The claim that no J = 0 row fails has no positive-control counterpart at out-degree >= 3 with J = 0 and parallel arrows beyond 41 early rows. This is acceptable as stated.

## New?

Mostly already recorded. E-126 (referee note) lists the four rows (expansions 6 820, 7 424, 7 822, 7 831) and asks for `checkCartan` to be run. E-113 states "out-degree >= 3: 3 494 rows, 0 rejects" and "no row with parallel out-arrows at v is a reject". E-124 gives the J_i = Hom(S_v, e_iA) reading. E-116/E-122 cover the out-degree 2 analogue. New here:

- The fifth row (7 836, out-degree 4, J support {1, 3, 5}).
- The answer to the open question: all of them fail, so the E-113 "0 rejects at out-degree >= 3" holds only for its prefix. The author notes this and it is consistent.
- The remark that the number of J_i != 0 equals the parallel multiplicity into 8. Grepped `research/` for "parallel multiplicity", "out-degree >= 3" and "7 836" and found nothing.

E-113's wording about parallel rows should be corrected by the chair so it does not stand against this result. I found no RETRACTION entry that is touched.

## Evidenced?

Rows, v, J and supports are listed per row, so the five-row claim is checkable. The single-sentence tally also checks. Missing:

- The control's cap (200, first by walk order) and the exact count (200, not 199).
- How many out-degree >= 3 rows the walk saw in total (E-113's 3 494 is from a different walk), so the fraction checked cannot be read off.
- The per-row discrepancy entries (values, not only the support). "Support exactly equal" is a weak statement if the values are +-1 on every J_i; the reader cannot tell.
- The commit or file for `experimentalist_dimji.py` equivalence. The script execs the round-033 helper, so the claim "same walk" depends on it being unchanged.

## Required for acceptance

1. Correct "199 ... J = 0" to the 200 in the tally and say that the J = 0 control is capped at the first 200 in walk order. Withdraw or restrict the "iff" sentence to the sampled rows.
2. State that the five rows may be only 3 orbits (mirror pairs unchecked), or check it.
3. Report the discrepancy values, not only the support.
4. Cross-reference E-113's "0 rejects" as superseded for walks past expansion ~6 800.
