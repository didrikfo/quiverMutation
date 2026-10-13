# Round 062 -- proceedings
kind: ordinary. Round 061 left no questions (none to settle). Merge of origin/main: already up to date. No suggested question taken.

## skeptic -- equal-dimension wrong algebras at step 12 / class 2 (T10)
Claim: `symcheck2` rejects 1392/1392 equal-dimension non-iso variants at class 1 step 12 and 18/18 distinct non-iso ideals at the 9 class-2 steps, accepts 594 iso ones, 0 disagreements. Referee toolsmith: minor revision (family qualifiers, assert no parallel arrows at class 2, dedupe class-2 variants by ideal and give hashes, widen the mu_4 truth search). Response did all four (class 2 is 36 distinct ideals, not 60; truth over the full torus mod 13). Decision: **accept**. Promoted as **E-186**. Closes the E-177 open cases; says nothing about J.

## theorist -- single-mutation LNA rule (T7)
Claim: proof of the E-182 rule for right mutation, both "only if" cases, and the Catalan factorisation, relative to the repo's `procedure`. Referee skeptic: minor revision (per-vertex check, "relative to procedure" in title, argue generation, note non-monomial n = 6 remark). Response did all four. Decision: **accept**. Promoted as **E-187**; H-010 stays SUPPORTED. The reading of the procedure against arXiv:2112.08129 is not done.

## experimentalist -- H-020 length gaps (T4, breadth)
Claim: widths 6, 7, 8 floating rules pass at the five skipped cells, 0 failures in 8 248 applications. Referee theorist: minor revision (claim wording of the length range, name the checks and commit). The chair made both wording edits in the submission; the optional third item was not done. Decision: **accept**. Promoted as **E-188**. Also found: `theorist_rulelen.py`'s `L >= w+5` filter silently skips cells.

## Consequences
None for any F-, H-, R- status. T7's "no proof" is replaced by a proof relative to the procedure; E-182's "proof only sketched" is superseded by E-187.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 061. No overnight run proposed.
