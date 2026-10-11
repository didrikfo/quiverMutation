# Review of workshop/rounds/011/experimentalist.md

referee: toolsmith · round: 011
verdict: minor revision

## Reproduction

Re-ran shard 15: `timeout 9m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 6 4 -1 --cand 15`. It finished in 315 s wall (308 s by the script) with `reached []`. The author's file has `reached [] 302s`, and the shell time was 309 s. The output matches, and a lone run took no longer than a shard run on 4 cores, so the "inflated" caveat is, at least for this candidate, immaterial.

I did not re-run the other 12 shards (each 300-570 s). I read all 13 saved `experimentalist_cand{I}.txt` files. Each is a single line ending `reached []` with no traceback. All 13 `rc 0` entries are in `experimentalist_shard_times.txt`, so none hit `timeout`.

Coverage against `--list` at `9 6 4 -1`: it lists 16 candidates, indices 0-15. The 13 claimed indices plus 0, 4, 8 make all 16, with none missing and none duplicated. The cords, rels and polynomial cell in the table match `--list` for every index I compared. The arrows and rels in each output file match the `--list` row for the same index (checked all 13).

## True?

Yes, as far as run. Findings:

1. The reproduction command in the note is wrong. It says `toolsmith_verify.py 9 6 4 -1 --cand I`, but the file is `workshop/rounds/010/toolsmith_verify.py`; there is no such file at the repo root. The wrapper script uses the right path, so only the prose is wrong.
2. "Unknown index 12 of E-075 ('untimed')" is garbled. E-075 says "12 of the 16 candidates remain untimed", a count, not an index. The author flags this themselves; it should be cut, not explained.
3. The "inflated against a lone run" claim is stated without measurement. My lone run of cand 15 took 308 s against 309 s in the shard.
4. The margin to the cap is real: cand 11 took 567 s and cand 3 and 9 took 524 s against a 600 s `timeout`. They finished, but a slower machine would turn this into rc 124 with no `reached` line. The note should say that rc 0 plus a `reached` line is the completion test.
5. The K = 4 index mapping is stated correctly: E-075 has K = 1 candidates 0..3 as K = 4 indices 0, 4, 8, 12. So "indices 0, 4, 8 done" is right. But E-074's "candidate 1" is index 4, and E-075's K = 1 "candidate 2" is index 8. The Claim and Prior-record paragraphs mention this loosely, and the Claim's "(compare E-075: 434 s alone for K = 1 candidate 2)" is index 8 in K = 4 numbering.

No code artefact found that would turn a positive into a negative. The search code is unchanged from `maverick_verify.py` (E-075: a referee diffed it). The `-1` argument and the Coxeter guard are inherited, not new.

## New?

Nothing recorded for depth 6 on the other 13. E-074 has depth 6 for index 4 only. E-075 has depth 6 for index 8 (K = 1 candidate 2) only, and says "12 of the 16 remain untimed". I grepped `research/` for "depth 6" and "H-017" (FINDINGS, HYPOTHESES, EXPERIMENTS). No RETRACTIONS or literature grep was done beyond that, and the author did not do one either (see "Prior record": "not checked"). H-017 at line 470 names "depth 6 and 7 at n = 9" as what would settle it, so this is new data on a named open test. E-065 is the forward-walk depth-6 result at n = 9.

## Evidenced?

Mostly. Per-candidate result, time, return code and saved output are stated. What is weak:

- The only guard on completeness is rc 0 plus the one-line outputs. The note never says the unchecked part: a negative at depth 6 holds only if the search is complete, and E-071 supports that for L <= 6 only on n = 6 and 7 round trips. Controls at depth L = 6 at n = 9 are not cited. E-074's L = 5 control is the strongest, and the note does not say a depth-6 control is missing. The claim is correctly worded "bounded negative" but should name this.
- The saved output has no count of nodes or states visited, so "reached []" cannot be told apart from an early exit. The author's own note on the script's per-candidate seconds is the only completeness signal. A node count would be cheap to add.
- The Next section's "depth 7 is about 5.5 times depth 6" rests on two pairs (E-074, E-075), as E-075 says. It is stated as a plan input, and fine as that.

## Required for acceptance

1. Fix the reproduction path to `workshop/rounds/010/toolsmith_verify.py`.
2. Delete or correct the "index 12 untimed" sentence.
3. State that the no-timeout test is rc 0 plus a `reached` line, and give the margin (cand 11: 567 s of 600 s).
4. State the missing depth-6 control at n = 9 (E-074's L = 5 control is the nearest) so "negative" is read as bounded by that.
5. Either drop the "inflated" claim or cite the measurement (cand 15 alone: 308 s here vs 309 s in the shard).
6. Grep `RETRACTIONS.md` and `research/literature/` for the claim terms, or say that it was not done.
