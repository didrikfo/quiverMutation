# Round 002 -- proceedings

kind: ordinary · chair: round session · referees: skeptic (x2), theorist
Answers from `STEERING.md` followed: round 001 q1 (strict reading first; done), q3 (no `tiltingPlus` promotion without a second non-monomial *gate-admitted* negative; not met). q2 was applied in round 001. No unanswered questions were left over from round 001.

## experimentalist -- T1 (revision)
* **Claim.** At n = 13 over 139 cores, H-021's mirror clause is no "exactly when" on the loose, strict or strict2 reading. Strict is false for 108 of 109 pairing cores (by construction); 7 cores hold a strict mirror with no reflection.
* **Referee (skeptic): accept**, fit reproduced exactly; census not re-run. Suggestions: say "108" is by construction; check the 7 at n = 14.
* **Decision: accept.** Wording point carried into E-058; the n = 14 check is a next question.
* **Promoted.** E-058; H-021 status line updated (still OPEN).

## skeptic -- T3 (revision)
* **Claim.** `3346` never pairs (key proof); `4056` pairs by key, but at n = 16 the walk gives `{0,3}` and two separate mirror-image orbits `{1}`, `{2}`, so "orbit = key class" fails.
* **Referee (theorist): minor revision**, rerun identical (237 s). Wanted: state `mirrorRow` is the F-026 dual and {1,2} equivalence is inferred; state n = 16 is one core at one length; flag E-056's superseded "7 / 123".
* **Decision: accept**, all three points are wording and I wrote them into the entry (overruling "revise" as the points are checkable and applied).
* **Promoted.** E-060; E-056 annotated in place with the corrected scan counts (309 / 155 / 121).

## scholar -- T5 (revision)
* **Claim.** On non-monomial parents at n = 5, 6, 7 the exact criterion rejects exactly the gate-refused vertices; no gate-admitted rejection beyond E-032 step 7; 0 skipped parents.
* **Referee (skeptic): accept**, all three runs reproduced.
* **Decision: accept.** Modest consistency check, recorded as such.
* **Promoted.** E-059. Not promoted: `isTilting` to the library (q3 not met in the strong sense).

Glossary: strict mirror, gate/guard. Merging `main` was a no-op.

## Decisions taken for the steering committee
None (no open question was unanswered).

## Questions for the steering committee
1. The scholar suggests promoting `isTilting` as a cross-check beside the gate, with pinned tests. **Recommend: not yet**; keep waiting for a gate-admitted rejection (the overnight audit). If unanswered, next chair takes this.
2. H-021 now fails on three readings. **Recommend**: next round the theorist restates it without the mirror clause (pairing up to an overhang, with `d(c)` tested against the H-020 head/tail difference) and the experimentalist runs the 7 survivors at n = 14. If unanswered, next chair takes this.
