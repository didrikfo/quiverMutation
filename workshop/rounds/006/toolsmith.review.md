# Review of workshop/rounds/006/toolsmith.md

referee: skeptic · round: 006
verdict: accept

## Reproduction

Ledgers for n = 10, 12..16 existed; I re-ran `toolsmith_orbitclass.py N` for each (resumes, all 484 units done, seconds) and `--same-orbit` (about 4.5 min, fresh walks). Every number matches: 139 closed cores at every n, 0 capped; key == orbit+mirror in 132/130/129/130/129/130 at n = 10/12/13/14/15/16; key-coarser 7/9/10/9/10/9; key finer 0, incomparable 0. The key-coarser core lists printed at n = 13 and 16 match the table. `--same-orbit`: eight walks of 20300, all closed; 28 pairs, 12 with 20300 shared rows, 16 with 0; in every disjoint pair the mirror of a's start is in b. Matches the claim exactly. The n = 10..16 ledger walks themselves were not re-walked (they are the author's ledger; the orbit rows are not re-derived here).

## True?

Checked the suspected weakness: `mirrors` in `batch.py` (l.387-413) is the loose reading, "offsets whose mirror row is in this orbit's walk". The script joins orbit A to the orbit of any offset o whose mirror lies in A. That is exactly the orbit-of-mirror relation, so the loose reading does not weaken the join: it is the correct definition of "orbit plus mirror", and the loose/strict distinction (E-058) only matters when a mirror is used as evidence of a reflection, which is not done here. No self-mirror orbit is joined to anything spuriously, since a self-mirror orbit only joins itself.

What the claim really says: orbit+mirror <= key in all 139 x 6 cases (key never splits an orbit+mirror class), and the 9/10 exceptions are cores where two self-mirror orbits (the parity classes `{0,2,..}{1,3,..}`) share one key. I confirmed the exceptions are all parity pairs at n = 13 and 16. This is not trivial: HYPOTHESES (l.41) records that under the reduced walk the key failed to be orbit-invariant in 6 of 1192 comparisons, so "key constant on orbit+mirror classes" was a live risk; it is now tested catalogue-wide. Caveat, the same as the author's: "refines" is an observation over 139 cores x 6 lengths, not a proof, and the author says so.

The 20300 claim holds; one wrinkle worth a sentence: the shared rows are across different (core, offset) pairs, e.g. `4056`@1 and `46`@3 have the same 20300 rows. The claim's "one orbit and its mirror" is right; the statement "the eight singleton orbits" is an orbit count across cores, not per core.

## New?

Grepped `research/` for "orbit+mirror", "orbit-plus-mirror", "parity class": nothing. E-060 already records `4056`@16 `{1,2}` as two mirror-image orbits and E-064 states that whether 20300 is one shared orbit was not checked. So new: the catalogue-wide comparison at n = 10..16 and the 20300 identification (12 shared, 16 disjoint). F-053's "each pair is one self-dual orbit" is already corrected by E-060.

## Evidenced?

Yes. Range (all 139 placed cores of `--max-word 4`, n = 10, 12..16), counts per n, the exceptional cores by name, and the closed/capped status are stated, and the scope limits (4-letter cores, n <= 16) are stated. Minor gaps, none affecting acceptance: n = 11 is skipped without comment; the key-coarser cores are not compared with E-061's seven (the author says so); the "sound prefilter" is correctly flagged unproved.

## Required for acceptance

None. Suggested, not required: say why n = 11 is omitted, and state in the claim that the refinement is observed, not proved.
