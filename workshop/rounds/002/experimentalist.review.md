# Review of workshop/rounds/002/experimentalist.md

referee: skeptic · round: 002
verdict: accept (two wording notes, none blocking)

## Reproduction

Re-ran `experimentalist_fit.py` on the saved census `experimentalist_census_n13.jsonl` (under 1 s). Output matches the submission exactly: 139 cores, 395 walks, d histogram 62/28/17/2/None 30, and the five-row (fit, loose, strict, strict2) table 108 / 1 / 10 / 13 / 7. Arithmetic checks: loose = 108+1+13+7 = 129; strict = 1+7 = 8; strict2 = 108+1+7 = 116; P = 109. I did not re-run the 5-minute census itself. It is round 001's script, unchanged, and the walks feed only the saved jsonl. I read the `readings()` code: loose/strict/strict2 are implemented as the docstring says (strict = an orbit holds the mirror of c@p without holding c@p).

## True?

No error found. Notes:
* "Strict" as coded (mirror of c@p in an orbit that lacks c@p) is false for every clean pairing by construction. The mirror lands inside the orbit that holds its own offsets. So "P => strict fails for 108" is close to a tautology and adds little beyond "the mirror acts inside the orbit". The author says so in Evidence. It does not change the conclusion, which rests on the 7 and the 13.
* The code's strict test (`p not in H`) is slightly weaker than the prose "different orbits, q != p". For an orbit holding p, the two agree. I found no case where they differ in this census.
* The d <= 6 fit is generous, so "no fit" is the strong direction. The submission says this. The 7 survivors are therefore firm counterexamples to "strict => P". The 108 direction depends only on the definition.

## New?

Partly. Recorded already: E-052/F-053 (pairing), H-021 "exactly when", E-055 (loose reading fails, 129/139, 20 cores). E-055 predicted that only `344 4044 4403` would survive the strict reading and said "that count has not been run". The new content is the strict and strict2 counts, the list of 7 (four more than predicted: `366 4404 4405 4605`), and the 13 that drop out. Committed fit code closes E-055's "not in the repository" caveat.

## Evidenced?

Yes. Range (n = 13, all 139 cores of `--max-word 4`, 0 caps), the three definitions, the counts and the core names are all given. Data and code are committed, so the claim is checkable without re-running the walks. Scope limits (n >= 14, two-cluster words) are stated.

## Required for acceptance

None. Suggestions:
1. In the Claim, say outright that strict is false for clean pairings by construction, so no reader takes "108" as an empirical result.
2. Before H-021 is restated, run the same table at n = 14 on the 7 to confirm they are not n = 13 artefacts (the submission defers this to the overnight menu).
