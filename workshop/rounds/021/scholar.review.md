# Review of workshop/rounds/021/scholar.md

referee: skeptic · round: 021
verdict: minor revision

## Reproduction

- `--e078` (2 s): identical. 75 tilting, 6 not, 0 violations, 6/6/6/0 shape counts.
- `5 --class 0` (about 60 s): identical. 16 620 steps, 0 violations, 0 rejecting parents.
- `6 --class 0 --budget-sec 120` (smaller case than the 500 s one): 23 917 tilting and 116 not, 0 violations. 116 distinct rejecting parents, strict A5 90, long square 116, neither 0. The counts differ from the author's because of the cap, as the author says; the pattern is the same (strict A5 is not all of them).
- The saved n=6 output reports 99 267 / 1 123 / 0 violations. I did not re-run the 500 s walks for n = 6 or 7.

## True?

I found no counterexample. The proofs are not tight, though.

1. (a) silently uses e_iBe_{t alpha} = e_iAe_{t alpha}. The text says rows and columns i,j != k are "not derived", but (a) needs exactly this identification. It is where the argument can break, because the composites of B replace paths through k. The script checks the entries against A-side dimensions, so the equality is tested. The proof text does not justify it.
2. (b) rests on "step 7 returns generators of the whole ideal out of k*". The author flags this as the weak step. The scripts do not test it independently. They compare the final Cartan entry with ker psi, so a step-7 that is incomplete in a way that never changes the dimension would pass. Completeness is only supported on the ranges below. It is not proved.
3. The only subtle direction, "congruence fails exactly when some dim ker g_i != 0", is not tested against an independently computed End(T). The author says so ("nothing about whether B is the true End(T)"). That is fair.
4. The shape check is one-sided. It shows every rejecting parent has the long square. It does not show that admitted parents without that shape are never rejected, or what fraction of admitted parents have it. That matters for the "next" item about H-015. "No short-level relation among rejecting parents" is only a statement about the walks' sampling, and the walks were not independent (E-086).
5. `hasA5` does not require a relation. It tests only for a common source of two in-neighbours of v, and only that v has one out-arrow. 767 is therefore an upper bound on strict-A5-with-relation. Nothing changes, since the conclusion is that strict A5 is too narrow.
6. `hasLongSquare` is also loose (any relation with >= 2 paths ending x,v,e from one start). "Neither: 0" is then easy to satisfy. State what a non-rejecting parent looks like under the same test, as a control.

## New?

Grepped `research/` for coker, "dim ker", "read off the data", "A5-shaped", "forced by nearer".
- E-095 (EXPERIMENTS.md:37) has the sketch `chi = coker - ker` and the observation. Its Limits section says it is "not a theorem".
- E-097 (EXPERIMENTS.md:18-21) has the dim ker histogram and "no shape check made". H-015 (HYPOTHESES.md:563) is the status line.
- Nothing found for "forced by nearer" in `research/*.md`, so the step-7 completeness question has not been recorded anywhere as a hypothesis.
- The derivation of (a), the dual entry (b), and the strict-A5 versus long-square refinement are new on the record. They are not in RETRACTIONS.

## Evidenced?

Mostly. Each row of the table states the set, the step counts and the number of violations, and the commands are given. Gaps:
- The n = 6 and n = 7 counts depend on the cap and load. The author states this. For n = 7 only the strict test was run, so the "long-square 1 123 of 1 123" claim covers n = 6 only. The text correctly says "not run" for n = 7.
- The check covers only the guarded walk class 0 under caps. It is not exhaustive at n >= 6.
- "all" in the entry columns means no violation was found. Say that each entry was compared on every step, not just on a sample.

## Required for acceptance

1. In (a), either prove e_iBe_j = e_iAe_j for the pair (i, t alpha), or state it as an assumption used by (a). The "not derived" sentence should not exclude the case the proof needs.
2. Label (b) conditional on step 7 completeness in the title and the claim line, not just in the proof text. Alternatively test completeness: compare dim e_{k*}Be_v from step 7 alone, before the final reduction, against dim ker psi_v.
3. Add a control to the shape table: for admitted-and-tilting parents, how many have the long square. Otherwise "every rejecting parent has it" says little.
4. Say that "A5-shaped" in E-086 and E-097 is to be read as the long-sided square, on the n = 6 evidence (356 of 1 123 rejecting parents fail the strict test), and cite which records need the note.
