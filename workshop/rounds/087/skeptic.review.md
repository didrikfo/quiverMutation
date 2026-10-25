# Review of workshop/rounds/087/skeptic.md

referee: theorist · round: 087
verdict: minor revision

## Reproduction

Re-ran four of the author's commands, all inside the 10-minute limit.
- `skeptic_homotopy.py 10 5 run 0 -1` gave `{'OK': 1008}`, 14 s.
- `10 5 run 1 1` (mirror convention) gave `{'OK': 1008}`, 13 s.
- `10 5 ctl_cox 0 -1` gave `{'DIM': 1008}`, 23 s.
- `10 5 ctl_seq 0 -1` printed rejections (NOTLINE and DIM examples). I did not recount its total of 1 008 rejected, so the author's split of 826/98/84 is unchecked.

The counts match the author's. Not re-run: `ctl_after`, the calibration cells, and the w5 L11, w4 L10 and w4 L11 cells.

I also read `verdict()` (lines 352-398). It does what the claim says:
- Hom(T,T[k]) is tested for k in -4..4, k != 0.
- Every End(T) vertex must be local, with Hom dimension at most 1.
- Arrows are found by irreducibility against composites through the other summands, and must form a line.
- Hom dimensions must equal the predicted path-nonvanishing pattern.
- Consecutive-arrow composites must vanish exactly where the prediction says.

For a monomial LNA that determines the algebra, so "End(T) = predicted LNA" is justified.

## True?

I found no error and no counterexample. Gaps in the statement, not in the result:
1. The tilting test covers |k| <= 4 only. The text writes "0 < |k| <= 4" in the claim but the Evidence paragraph never says why that range suffices. An argument is available: the mutation sequences are short and the complexes have bounded length, so Hom(T,T[k]) = 0 for |k| beyond the span is automatic. That argument is not given. Hom(T,T[k]) for k > 0 and k < 0 are both tested, so the check is two-sided.
2. "T is tilting" also needs T to generate K^b(proj A). This is inherited from iterated silting mutation starting at A, which is the Aihara-Iyama citation the author marks UNVERIFIED. The claim should say so; the script does not test generation.
3. The convention (side and sign map) was fitted on other cells and then frozen. That is fair, and the mirror convention agrees, so the data fixes it up to the mirror symmetry. But the "rejected 100%" calibration lines are the only evidence that the wrong convention fails, and the calibration output is not in the saved `.out`. Only the cell output is saved, with the first control lines truncated by `tail`.
4. The claim uses `lnaMoves.VERIFIED_MOVES` and `anchorOf` as inputs. The title says "does not use the repo's mutation code", which is true, but the sequences themselves are the repo's. The result therefore audits "sequence implies predicted End", not "the rule table is the right table". The claim states this; the title does not.

## New?

Nothing found for an independent silting-mutation engine deciding End(T) = LNA, among these terms grepped in `research/`: independent, homotopy, silting, End(T).
- E-207 (HYPOTHESES/EXPERIMENTS line 270) re-implemented matching, application and the Coxeter test, but kept the repo engine.
- E-188 and E-159 use `verifyMove`.
- E-161 and E-184 use an independent Hom(T,T[±1]) test, but generation of the tilting step came from the repo.
- E-181 and E-174 compare dimensions and End(T) against the repo's rewrite.

So the independent engine plus the full algebra decision is new as a check. It is not new mathematics, as the author says. The author's prior-record paragraph is accurate. The ctl_cox control (a different LNA with the same Coxeter polynomial) is also a new kind of power check.

## Evidenced?

Mostly yes. The scope line gives widths, lengths, application counts and exclusions, and the control table has counts. Missing:
- The calibration output, and the full control outputs. Both are only quoted.
- The caveat on ctl_cox is honest, but "most probably inequivalent" is not a result. Its counts show only that the labelled End algebra differs. That is all the control needs for the stated purpose (power against Coxeter-equal alternatives).
- The extra cells (3 168, 528, 1 716) are listed with no timing or output file. Only the L10 W5 output is saved.

## Scope

The title says "H1 at width 5, length 10 (one cell)". This matches the cell: floating rules, 24 rules, 1 008 applications. The title also names "tilting with End(T) = predicted LNA", and the check does give that. The statement "H1 holds" overreaches slightly if H1 includes anchored rules or other widths, but the title says one cell and the scope line excludes the rest. Narrowed wording for the Prior record: "supports H-020 at this cell; adds no statement about the rule table's completeness".

## Required for acceptance

1. Add one sentence in Evidence stating why |k| <= 4 suffices for the tilting test (complex length bound), or extend the range to the actual maximum degree range of the complexes and say so.
2. State that generation of K^b(proj A) by T is inherited from silting mutation, not tested, and that this rests on the UNVERIFIED citation.
3. Save the calibration run output (and the untruncated control outputs) next to `skeptic_homotopy_L10_W5.out`, or cite the exact commands and counts so they can be reproduced. The commands are already listed; only the saved outputs are missing.
4. Add one line in the title or claim that the rule sequences themselves are `VERIFIED_MOVES`, so this audits sequence-to-End, not the table.
5. Give the output location for the w5 L11, w4 L10 and w4 L11 cells.
