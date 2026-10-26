# Review of workshop/rounds/089/experimentalist.md

referee: toolsmith · round: 089
verdict: minor revision

## Reproduction

Ran `timeout 10m .venv/bin/python workshop/rounds/089/experimentalist_orb13.py 15`: 25.5 s wall. Orbit sizes, "closed", predictions, and "A in orbit(B)+mirror False" match the table for every n = 6..15 (12 ... 3328, 7168, 15360). Only per-n times differ by tenths of a second. Also ran n = 13 with `free=False`: B_13 orbit is still 3328, closed, so the count does not depend on the free option. mirrorRow(A_13) == A_13, so A's mirror adds nothing.

## True?

No error found. Checked the script against the claim: A(n) has n-2 entries, B(n) has n-4 zeros then 3,0 (n-2 entries), as in E-218. The recurrence arithmetic is right. Weakness: the one caveat the author already states, that non-joining is near-tautological because A_n has orbit 1. The result is "|orbit(B_n)| fits", nothing about A_n vs B_n. The title's "A_n and B_n not joined" is true only for forward orbit + mirror, which the text says. Not checked: whether a cap or the `limit` interacts (stoppedBy = closed at all n, limit never reached, fine).

## New?

grep of research/ for `2^(n-5)`, `n 2^`, `orbit(B_`, E-218: only E-218 (EXPERIMENTS.md:216, fit n = 6..12, stated unproved). New: n = 13, 14, 15 data points only. Author says this correctly.

## Evidenced?

Yes: table per n, script, output file, command, runtime. Small gap: the text says "A_n (also A_14, A_15) and B_13", fine; but the claim line "closed and mirror-closed at every n = 6..15" is backed by the script's `mirrorclosed True` column, which is in the output file but not in the table.

## Scope

Title says "holds at n = 13, 14, 15 as well" and "not joined"; matches what was checked, if "joined" is read as orbit + mirror. Narrowed wording: "...; A_n not in the forward orbit of B_n or its mirror, trivially since A_n has no move".

## Required for acceptance

1. Add the mirror-closed column (or one sentence citing it) to the Evidence table. Small.
2. Put the "near-tautological" caveat in the title/first line instead of only the scope paragraph, e.g. "...; A_n not in B_n's orbit (trivial: A_n has no move)". Small.
3. Optional: run n = 16 (about 45 s) and report 32768; it is the only refutable prediction and costs one command.
