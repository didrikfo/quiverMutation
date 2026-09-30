# Experimentalist notebook (rewritten each round)

## What I now believe (after round 003)
- The 7 survivors (`344 366 4044 4403 4404 4405 4605`) are a parity effect. Swept n = 8..18 (one core = 1 to 3 s): at n = 12, 14, 16, 18 all 7 pair by reflection (d = 0 for 4, 1 for 2, 2 for `4405`), no strict mirror; at n = 13, 15, 17 none pairs and all 7 hold strict mirror. n = 9..11 mixed (1 to 5 fail), few offsets, weak.
- So T1's "pairing with a defect" question: at n = 14 no defect. The n = 13 defect is an odd-n one (singletons of equal size, e.g. `344` at 13: `{2}50 {4}50`; at 15 `{2}64 {6}64`).
- The n = 13 census (139 cores) is therefore one parity class only. Its "30 no-fit" and "13 mirror-without-fit" verdicts may not carry to even n; unknown.
- Round 002 still stands: H-021's mirror clause fails on every reading at n = 13; pairing itself 62 exact / 47 overhang / 30 none at 13.

## What I tried
- `experimentalist_census.py` has no `--plan`; I sized by running one core. `--cores a,b,c` works; all runs completed, 0 caps.
- Merged data: `workshop/rounds/003/experimentalist_census_7cores_n8_18.jsonl`; table script `experimentalist_table.py`.

## What I would do next
1. Cheap: run the full 139 cores at n = 12 (fewer offsets, faster than 14) to see whether the 30 no-fit and 13 mirror-without-fit change class; and a sample of the 62 exact cores at n = 15 to see if they lose pairing at odd n.
2. The overnight n = 14 census, then fit; compare class per core against 13.
3. Test d(c) against `|head - tail|` of H-020 (unchecked, four rounds running).
4. Toolsmith: `batch.py` task for the census with a `45` pin; add `--plan`.
- Watch: n = 9..11 fits use 2 to 5 offsets, easy to satisfy or fail; do not read them as the same phenomenon. Only 7 cores were swept, so "parity" is shown for these, not for all.
