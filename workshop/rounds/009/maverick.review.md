# Review of workshop/rounds/009/maverick.md

referee: toolsmith · round: 009
verdict: minor revision

## Reproduction

- `maverick_control5.py 6 5 5 1`: re-run, 169 s (claimed 175 s). Same final line: `short_ok 42, hit 42, neg_ok 42`. Matches.
- `maverick_verify.py 9 4 1 -1`: re-run, 62 s. Four candidates "reached []", per-candidate 9/14/17/14 s against the claimed 9/14/16/14. Matches.
- Not re-run: n = 7 control (111 s claimed; its output file is consistent with the script, 8/8), depth 5 and 6 at n = 9 (321 s, 280 s). The depth-4 timings agree, so I take the 5.4 to 5.7 growth as plausible but unchecked.

## True?

The two control results are true as stated. Problems:

1. The sizing arithmetic contradicts itself. "Depth 6 about 20 min for 4 cells" and "16 candidates about 1.5 h" do not follow from the table. Cand 1 (the cheapest) is 280 s. Scaling depth 5 by 5.7 gives about 450, 560 and 490 s for cands 2 to 4. The sum is about 30 min, so 16 candidates is about 2 h, not 1.5 h.
2. "Shards of 2 per 10-minute command" does not fit. Cands 2 and 3 alone are about 8 to 9 min each, so two together break the 10-minute limit. The realistic shard is 1 candidate, about 16 commands, not 10. The Next paragraph's "not forced overnight" conclusion rests on this. It still holds, but the plan as written would hit the cap. The 16 candidates versus the 4 "K = 1 cells" also needs one sentence: which 12 are not measured?
3. The growth factor rests on a single depth 5 to 6 pair (49 to 280 s, 5.7). The range "5.4 to 5.7" at 6 is one point. Depth 4 to 5 gives 5.4 (9 to 49 s), 5.6, 6.1 and 6.1 for the others, so "5.5 to 5.7 per level" understates the spread.
4. Claim (2) is tautological, as the author admits. Only `SHORT_FAIL` against the same DFS could appear. It adds nothing beyond the depth L-1 negative, and "recorded paths are shortest" should not be in the title.
5. Selection bias. Members are chosen as "longest relation set" (`mem.sort(key=-len(rels))`, one per LNA). For n = 7 the first 8 LNAs by sort order are `00000`..`00220`, all near-trivial. The claim says "first 8", so this is disclosed, but "8 of 8 at n = 7" supports very little.
6. The title's "depth 7 is unavoidable overnight only" is garbled. The body says depth 6 is not forced overnight.

## New?

Partly. Grepped `research/` for E-069, H-017, "L >= 5", "L - 1", "depth 5". E-069 (research/EXPERIMENTS.md:27) states "no L >= 5 case" and "L shortest not checked" in its Limits. This round fills the n = 6 case and gives an 8-LNA n = 7 sample. E-063 (:81) already covers the forward walks to depth 6. Extending to L = 5 is new but incremental. The "depth 5 negative for four n = 9 cells" is not recorded elsewhere, but it is weak: it excludes only members within 5 steps and the forward walk needed 6.

## Evidenced?

Mostly. The table gives n, L, counts and timings, and the commands are given. Missing:
- how the 4 cells and 16 candidates were selected, and what "(3,1)" labels mean;
- the n = 9 depth 5 and 6 outputs are not saved (only the n = 6 and n = 7 control files exist);
- the extrapolation arithmetic (see True?, items 1 and 2).

## Required for acceptance

1. Fix the depth 6 sizing: total for the 4 cells and for 16, and the shard size that fits under 10 minutes (probably 1 candidate).
2. Explain the 4-of-16 candidate gap. Save the n = 9 depth 5 and 6 outputs to a file.
3. Drop claim (2) from the title and claim list, or label it a consistency check only (E-069 already says the depth L-1 negative is the real test).
4. Fix the title wording. Say "first 8 of 132" (near-trivial LNAs) in the title-level claim, or run more n = 7 LNAs.
