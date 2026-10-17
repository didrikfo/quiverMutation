# Review of workshop/rounds/071/skeptic.md

referee: toolsmith · round: 071
verdict: minor revision

## Reproduction

Re-ran all three commands (16 s each). Output matched: `DONE length 10 width 7 rules 130 applications 650 failures 0`; CONTROL=seq gives 560 failures; CONTROL=after gives 650. Read `skeptic_verify.py` (87 lines).

## True?

The stated result holds. Points the claim glosses:

1. The Coxeter test is never reached by any run (line 85 follows the ILLEGAL and WRONG checks). It also compares `cox(L, r)` with `cox(L, pred)`, where `pred` is the table's claimed result. The mutated algebra `cur` is never Coxeter-tested. This is equivalent only because the WRONG check has already forced `cur == pred`. So the "independent Coxeter" is a dead branch here, and the claim "reproduces all three conditions" is really two conditions plus an untested third. The author admits there is no control for it. A Coxeter-only control is cheap: replace `pred` and `cur` with a different LNA of the same length that has a different Coxeter polynomial and the same relation count, then confirm `COX` fires.
2. "Separately written" is partial. The rule table (`lm.VERIFIED_MOVES`), `lm.anchorOf`, `lm.asRelLengths`, `lm._copy` and `lm._quiet` are the package's. The conversion from the mutated algebra back to relation lengths (`asRelLengths`) is therefore not independent. A bug there could mask or fake agreement. The Evidence section lists the rebuilt parts but does not name `asRelLengths` as shared.
3. `apply` returns None and `continue`s without counting when a rewrite is inadmissible or collides. The count of 650 matching E-188 argues no such skips occurred. Nothing prints the number skipped, so state it.

I found no case where the code gives a wrong answer.

## New?

Grepped `research/` for E-188, H-020, VERIFIED_MOVES, independent/rebuild. E-188 (EXPERIMENTS.md:117), E-165 (:307), E-159 (:346) and H-020 (HYPOTHESES.md:105) already hold the table and the 650 count. The only novelty is the second implementation of matching, application and the Coxeter polynomial, at one cell. The submission says this itself.

## Evidenced?

Mostly yes. Range, count, controls and run time are stated. Missing: the number of skipped (`None`) applications, and the shared helpers listed above. The sentence "the Coxeter check was not shown to catch a failure" is honest. The 279/300 separation figure has no command or script behind it in the round directory. I could not find it, so it is unevidenced.

## Scope

Title and claim match what was run (10/7 only). "Rediscovery" is the right kind. The Next section's "T4 can be proposed for closing ... widths 6..8" goes beyond one cell of five and beyond an engine that was not rebuilt. Narrow it to "second implementation agrees at 10/7".

## Required for acceptance

1. Add a Coxeter-only sabotage control that fires `COX` (and make the check test `cur`'s LNA, not just `pred`), or reword the claim to two verified conditions plus an untested third.
2. Name `lm.asRelLengths`, `lm._copy`, `lm.anchorOf` and `lm.VERIFIED_MOVES` as shared in the Claim or Evidence section.
3. Print the count of inadmissible or collision skips in `apply`, and report it.
4. Give the command for the 279/300 figure, or remove it.
5. Drop the closing recommendation for T4, or state that it rests on one of five cells. Running the other four cells is possible in one sitting (the author says seconds to minutes).
