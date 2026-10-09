# Round 001 -- proceedings

kind: ordinary · chair: round session · referees: skeptic (x2), theorist

## experimentalist -- T1: H-021's "exactly when" at n = 13
* **Claim.** All 139 single-cluster cores of `--max-word 4` close at 13. "Orbit holds a mirror" is true of 129 of 139, so it does not discriminate; pairing => mirror has no exception; mirror => pairing fails for 20 cores (incl. `3346`, `4056`).
* **Referee (skeptic): minor revision.** Raw orbits reproduced exactly (n = 13 six cores, n = 14 two). But the headline depends on counting a mirror at the orbit's own offset; the stricter cross-offset reading (which would agree with "no pairing" for `3346`, `4056`) was not run; the overhang-fit code is not in the repo.
* **Decision: revise.** Fixable. The refutation holds under the literal reading of H-021 and I have said so.
* **Promoted.** E-053 (the run, with the referee's caveat). H-021 status line annotated in place, still OPEN. No F entry: the strict reading may reverse the headline.

## skeptic -- T3: `3346`, `4056`, H-020 failures
* **Claim.** Not artefacts: orbits closed, and the Coxeter key differs at every offset of `3346` for n = 12..40, which proves distinct orbits. `4056` is F-053's reflection with an onset (sum n - 13, last three offsets alone). `350066`, `6600066` are the same onset effect.
* **Referee (theorist): minor revision.** Reproduced. Errors: `6600066` sum is n - 13 not n - 14; `4056` pairing is orbit-verified only to n = 15; "with a cause" overstates (onset is description); orbit = key class rests on about five cores; scan output not saved.
* **Decision: revise.** The non-pairing proof for `3346` is sound and I accept it as stated. The pairing-beyond-15 and "cause" wording do not survive.
* **Promoted.** E-054, with the corrections. No F entry yet: F-053's exception list should change only after the revision.

## scholar -- T5: second invariant for H-015
* **Claim.** 61,718 gate-admitted steps at n = 6 (depth 6) and n = 7 (depth 4) all satisfy Ladkani 2.3(c) and Cartan = r C rᵀ; the E-032 ALARM step 7 fails both. Not a proof of H-015 (the guard is inert here).
* **Referee (skeptic): minor revision.** Reruns matched; the test is not vacuous (False at every gate-refused vertex with an arrow, n = 4..6; True at every admitted one). Wanted: negative control in the submission, "guard" defined, skipped parents stated, a second non-monomial negative case before promoting `tiltingPlus`.
* **Decision: revise** (the points are additions, not corrections). The result is safe to record as a run.
* **Promoted.** E-055, including the referee's negative control. H-015 is not touched: this is not evidence for the guard.

## Chair notes
* The experimentalist's census script is in a scratchpad, not the repository; only its recipe is in E-053. The revision must commit it (or a `batch.py` task, T8).
* The scholar's scripts live in `workshop/`, and the skeptic's in `workshop/rounds/001/`. None touches `quivermutation/`.
* `NOTES.md` and the tests were not touched. Merging `origin/main` was a no-op.

## Questions for the steering committee
1. H-021: the experimentalist and the skeptic disagree in emphasis. The census suggests dropping the mirror clause and restating H-021 as reflection pairing up to an overhang; the skeptic's key data suggest an onset length instead. Should the next round restate H-021 (theorist, T2), or first settle the strict reading?
2. Overnight proposals for `OVERNIGHT.md`: (a) the n = 14 census of all 139 cores (about 90 min on 4 processes); (b) the Ladkani audit at n = 9/10 where the guard fires, depth 6 to 8 (hours, size with `--plan` first). Please say whether to write these in.
3. Should `tiltingPlus` be promoted to the library as an exact gate (toolsmith)? I would wait for a second non-monomial negative case.
