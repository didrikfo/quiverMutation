# Round 003 -- call
kind: ordinary
## Assignments
- theorist: Restate H-021 without the mirror clause (pairing of offsets up to an overhang, per STEERING round 002 q2), and test whether the shortfall d(c) equals the H-020 head/tail difference on the recorded data (`workshop/rounds/002/skeptic_scan_n24.txt`, `experimentalist_census_n13.jsonl`). Outputs: a statement (T1) and, if reachable, the formula for the centre s(c) (T2). (thread T1, T2)
- experimentalist: Run the 7 survivors (`344 366 4044 4403 4404 4405 4605`) at n = 14 with `workshop/rounds/002/experimentalist_census.py` (size with `--plan`; shard so no command exceeds the limit) and report pairing, reflection and strict-mirror status against n = 13. (thread T1)
- toolsmith: T8: add a resumable `batch.py` task (or committed script) for the orbit report of E-052 with a test pinning `45` at n = 13; note key class is not orbit (E-060), so any prefilter must use orbit-plus-mirror. Keep code small, run only the touched tests. (thread T8)
## Revisions due
(none)
## Referees
- theorist: refereed by skeptic
- experimentalist: refereed by theorist
- toolsmith: refereed by skeptic
