# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 1
next_round_kind: ordinary

## Open threads

Round 001 worked T1, T3, T5 (records E-053, E-054, E-055). Each line: id · question · suited to · status.

- **T1** · H-021 "exactly when": under the literal reading it fails (129 of 139 cores at n = 13 hold a mirror; 20 hold one with no reflection, E-053). Open: the *strict* reading (mirror of `c@p` in the orbit of `c@q`, `q != p`) over the 139 cores; the n = 14 census (not run); commit the census script. Restate H-021 once done. · experimentalist (revision), theorist · revising
- **T2** · Why offsets pair by a reflection (F-053): find the mechanism and a formula for the centre `s(c) = n - k(c)`. New data (E-054): `k(4056) = 13`, `k(45) = 8`, `k(3344) = 6`, `k(556) = 10`; onset length below which nothing pairs; `skeptic_scan.py` at n = 24 gives 123 of 585 words with no pair and 7 with a triple. Candidate: is `d(c)` = the H-020 head/tail difference? (untested) · theorist · open
- **T3** · `3346` never pairs (Coxeter key, n = 12..40): settled as real by key separation (E-054). `4056` is a reflection with an onset (orbit-verified to n = 15). Open: walk one `4056` orbit at n = 16; the other four H-020 failures (census ledgers are not in the repo); the 20-core no-fit list, esp. `3033 3034 3044 3303`. · skeptic (revision) · revising
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position, only anchored and edge moves see the ends (F-051). State and prove the lemma. · theorist · open
- **T5** · H-015: Ladkani 2.3(c) (arXiv:1001.4765) is the exact per-step criterion; it agrees with the gate on 61,718 steps at n = 6, 7 and catches the E-032 ALARM step (E-055). Not evidence for the guard (inert at n <= 7). Open: a second non-monomial negative case; the audit where the guard fires (n = 9, 10; overnight); promote `tiltingPlus` to the library as `isTilting`. · scholar (revision), toolsmith · revising
- **T6** · H-017, a quipu carries more relations than cords: the Euler-form count; test on recorded classes at small `n`. · maverick, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: a `batch.py` task or script printing E-052's orbit report for a core at a length, with a test pinning `45` at 13. Now also wanted: the Coxeter-key prefilter (partition offsets by key; key difference proves separation, equal key does not prove pairing) and the T1 census as a resumable task. · toolsmith · open
- **T9** · H-019, H-013: long runs, not for a round. Overnight proposals waiting on the human: the n = 14 census of all 139 cores (about 90 min, 4 procs); the Ladkani audit at n = 9/10 depth 6 to 8. · any · parked

Suggested next round (the chair may change it): the three revisions above fill `researchers_per_round`; the theorist on T2/T4 goes in the round after unless STEERING asks otherwise.

## Awaiting revision

- `rounds/001/experimentalist.md` · experimentalist · run the strict-mirror reading over the 139 cores at 13 and say which reading is refuted; put the overhang-fit code where it can be re-run; reword the `3346` sentence. (review: `rounds/001/experimentalist.review.md`)
- `rounds/001/skeptic.md` · skeptic · correct `6600066` to sum n - 13; separate orbit-verified (n <= 15) from key-only claims; walk one `4056` orbit at 16 or say it is unconfirmed; save the scan summary; drop "with a cause". (review: `rounds/001/skeptic.review.md`)
- `rounds/001/scholar.md` · scholar · add the referee's negative control; define "guard" and state the skipped parents; a second non-monomial negative case, or say there is one. (review: `rounds/001/scholar.review.md`)

## Requests between personas

- toolsmith: the Coxeter-key prefilter and a committed census task (skeptic, experimentalist).
- experimentalist: test the key partition against the orbit partition over the whole `--max-word 4` catalogue at 13 (skeptic; find a key class the walk splits).
- theorist: the centre `s(c)` formula from `skeptic_scan.py` data (skeptic).
- skeptic: check `tiltingPlus` on non-monomial parents (scholar).

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 001 | - |
| theorist | - | 001 |
| skeptic | 001 | 001 (x2) |
| scholar | 001 | - |
| toolsmith | - | - |
| maverick | - | - |
