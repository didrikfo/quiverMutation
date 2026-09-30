# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 2
next_round_kind: ordinary

## Open threads

Round 002 revised T1, T3, T5 (records E-056, E-057, E-058). Each line: id · question · suited to · status.

- **T1** · H-021 "exactly when": the mirror clause fails on all three readings at n = 13 (loose 20 counterexamples; strict false for 108 of 109 pairing cores by construction; 7 cores `344 366 4044 4403 4404 4405 4605` hold a strict mirror without a reflection; E-056). Open: restate H-021 without the mirror clause (pairing up to overhang); run the 7 at n = 14; is the `344`-type family a reflection with a defect? · theorist (restatement), experimentalist (n = 14 for the 7) · open
- **T2** · Why offsets pair by a reflection (F-053): mechanism and formula for the centre `s(c) = n - k(c)`; is `d(c)` the H-020 head/tail difference (untested)? Data: `skeptic_scan_n24.txt` (309 one centre / 155 key class >= 3 / 121 none). · theorist · open
- **T3** · `3346` never pairs (proved by key, E-054). `4056` at n = 16: `{1,2}` are two mirror-image orbits, key pairs them (E-058), so "orbit = key class" is false there; compare orbit-plus-mirror against key over the `--max-word 4` catalogue at 14..16. Open: other cores that split; the other four H-020 failures (ledgers not in repo). · experimentalist, skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position (F-051). State and prove. · theorist · open
- **T5** · H-015: Ladkani 2.3(c) agrees with the gate everywhere tested, incl. non-monomial parents at n = 5..7 (E-055, E-057); only gate-admitted rejection is E-032 step 7. `isTilting` not promoted (STEERING q3). Open: the audit at n = 9/10 (overnight, Menu 4). · scholar, toolsmith · waiting on overnight
- **T6** · H-017, a quipu carries more relations than cords: the Euler-form count; test on recorded classes at small `n`. · maverick, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: a `batch.py` task for E-052's orbit report with a test pinning `45` at 13; Coxeter-key prefilter (note: key class is not orbit, E-058; use orbit-plus-mirror); the T1 census as a resumable task. · toolsmith · open
- **T9** · H-019, H-013: long runs. Overnight waiting on the human: n = 14 census of the 139 cores (about 90 min, 4 procs; use `workshop/rounds/002/experimentalist_census.py` and `experimentalist_fit.py`); Ladkani audit n = 9/10. · any · parked

Suggested next round (the chair may change it): theorist on T1 restatement and T2; experimentalist on the 7 cores at n = 14 (if 10-min commands allow) or the orbit-plus-mirror vs key comparison; toolsmith on T8; maverick on T6.

## Awaiting revision

(none)

## Requests between personas

- experimentalist: orbit-plus-mirror vs key classes over the catalogue at 14-16 (skeptic).
- theorist: at n = 16 `4056` has a self-dual pair `{0,3}` and a mirror pair `{1,2}`: is the split a function of parity or of position relative to the middle? (skeptic)
- theorist: the centre `s(c)` formula (round 001 request, still open).
- toolsmith: committed census task; key prefilter with the E-058 caveat.

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 002 | - |
| theorist | - | 002 |
| skeptic | 002 | 002 (x2) |
| scholar | 002 | - |
| toolsmith | - | - |
| maverick | - | - |
