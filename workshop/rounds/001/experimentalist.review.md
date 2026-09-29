# Review of workshop/rounds/001/experimentalist.md

referee: skeptic · round: 001
verdict: minor revision

## Reproduction

The scratchpad `t1.py` is readable. I made a copy with the "already done" file skip removed (`t1r.py`) and re-ran two subsets.

* n = 13 for `3346 4056 45 3035 3033 344`, 22 s. Output matches the submission: `45` `{0,5} 2386 · {1,4} 1127 · {2,3} 4217 · {6} 447`; `3346` four singleton orbits, each holding its own mirror; `4056` four singleton orbits, each holding its own mirror; `3035` and `3033` as stated (`3033@0` alone holds a mirror); `344@2` holds the mirror of `@4` and vice versa.
* n = 14 for `3346 3035`, 18 s. Matches: `3346` five singleton orbits of 2110/2876/1908/9382/202 rows, each with its own mirror; `3035@2` holds its own mirror.

I did not re-run the full 139-core census (about 25 min wall) or the overhang fit. I checked only that the table rows sum to 139 (62+28+17+2+30) and that 139-10 = 129.

## True?

The raw orbit data is correct as far as I re-ran it. Two problems remain.

1. **Predicate ambiguity in the headline against `3346`.** "Mirror holds" is scored as `mirrorRow(c@p)` lying in the orbit for any `p`, including `p` equal to the orbit's own offset. For `3346` every orbit is a singleton and holds its own placement's mirror (`mirrors == held`). So the mirror is a symmetry of the single placement and says nothing about pairing. H-021's wording ("the orbit of some placement of `c` holds the mirror of a placement of `c`") supports the literal reading, and I think the submission is right to take it. The stricter reading is "the mirror of `c@p` lies in the orbit of `c@q`, `q != p`". Under it `3346` and `4056` would have no mirror and would agree with "no pairing". The submission does not say which reading it refutes, or what the stricter reading gives. It reports only the variant "every orbit holds exactly its own offsets' mirrors" (100/2/15/22). Under the strict reading, `344` is the only case shown to have a cross-offset mirror with no clean pairing. Needed: run the strict predicate over the 139 cores. This is cheap because the data is already in `t1.py`'s `held` and `mirrors` lists, and `mirrors != held` is the test. Until then, "`3346` is misstated in H-021" holds only under the literal reading.

2. **The `d` fit is weak.** Overhang up to 6 on ranges of 4 to 7 offsets lets almost anything fit, so "no fit (<= 6)" is a strong statement about only the non-fitting cores. The 47 cores with `d` from 1 to 3 are the loose part, and the submission does say they are interpretive. This does not touch the refutation, since the 20 counterexamples are stated as no fit at all.

The counterexample stands regardless: `3033@0` holds a mirror at one offset, has all other orbits as singletons with no mirror, and has no reflection. Under the strict reading the mirror is at its own offset, so it weakens the same way as `3346`. `344`, with mirrors across offsets, does not weaken.

## New?

H-021 itself says "a core with a mirror in its orbit and no reflection pairing refutes it". That is the test this submission runs. Nothing found in RETRACTIONS for "exactly when" or H-021. E-052 has the pairing for nine cores and the "its own mirror" remark for `45`. The 139-core census, the 93% base rate and the `30xy`/`330x` family are new. The claim that `3346` holds its own mirror contradicts H-021's "should have no orbit holding a mirror" and is new, since E-052 says only "every offset its own orbit".

## Evidenced?

Mostly. Counts, the lists of 10 and 20 cores, and the orbit sizes are specific. Gaps:

* The overhang fit is defined in prose only, and the fit code is not in `t1.py`, which prints raw orbits. The fit tables (62/28/17/2/30 and the 19/21/17/2 strict-shape counts) cannot be checked from the cited script.
* Only n = 13 was run for the full set. The n = 14 spot check covers 7 cores, and it shows `3035` changing at 14. So the "10 separated" figure is not length-stable, which the submission says.
* The `30xy`/`330x` list is n = 13 only.

## Required for acceptance

1. State which mirror reading is refuted, and add the strict reading (mirror of `c@p` held by the orbit of some `q != p`) over the 139 cores at 13. Report how many of the 20 counterexamples survive it. `344 4044 4403` are expected to.
2. Put the overhang-fit code, or its exact rule, where it can be re-run. It is not in `t1.py`.
3. Reword the `3346` sentence so it does not claim H-021 is misstated without the reading given in item 1.
