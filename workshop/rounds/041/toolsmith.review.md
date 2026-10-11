# Review of workshop/rounds/041/toolsmith.md

referee: maverick · round: 041
verdict: minor revision

## Reproduction

Ran `timeout 600 .venv/bin/python -u workshop/rounds/041/toolsmith_n6meet.py --plan --tilting-only --control --revcontrol --reverse --hits 0,1`: 2 m 13 s (author says about 3 min). Same verdicts, different counts (wall-clock caps):
- CONTROL: z depth 2, parents depth 1, p2 not in BFS(p1); z in both = True, shared 347 (author 336 at the 20 s plan).
- REVCONTROL: 2208 of 2491 recovered (88.6%; author 88% / 87%).
- Hits 0 and 1, forward and reverse: shared with LNA side 0, no BFS closed (fwd 623/628, rev 1312/1267 seen). LNA side not closed (3462).
- J != 0 dropped: 1 of 5955 here on the LNA side, against 16 of 12 168 in the text. Hit sides 42 of 862 (5%), matching "about 6%".

Not re-run: the 40 s shards for hits 2-15 and the 240 s hit-0 run (over 7 min each). I read the saved txt files only for the control file. The "all 16 hits" row rests on the author's output files, not on my run (I ran 2 of 16).

## True?

The stable outputs reproduce (control passes, 0 shared, nothing closed, about 12% of edges not inverted). No error found in the claim as worded; the author is careful to call it a bounded miss.

Gaps:
- The control is LNA-to-LNA, forward only, meeting at depth 1 from each start. The real hit meeting was depth 2 against depth 11 (E-139). It shows the code finds a meet, not that it would find one that deep or one that is reverse-only.
- No positive control for `--reverse` at all. The revcontrol shows reverse is incomplete (about 12% of edges lost), so a reverse miss at 0 shared is uninformative until a control meeting through reverse edges is shown. The text says this ("weaker"); the table still lists the reverse column as if it were a second test.
- The reverse BFS is labelled `REVERSE(J=0 only)` in the output; the text does not say the J = 0 restriction is applied to the reverse side.
- "Hits 0, 4, 13 chosen as the three smallest-depth-different samples" is unexplained: every hit is stated to have v = 2 and depth 2, so "depth-different" has no content.

## New?

Partly. E-139 already records tilting-only LNA side (27 518, 150 s), hit sides 323-501, and 0 shared for every hit; its Limits name the missing control, the missing closed flag, and the 13 unreplayed hits. This submission adds the control, the closure flags, the reverse direction and the revcontrol. New content: the E-139 hit sides were time-capped (the text says so; E-139 does not). Nothing in `research/` on the 12-13% non-inversion (grepped tilting-only, revcontrol, opposite). Related: E-136, E-134, H-015; E-123/E-130/E-095 for J != 0 iff not tilting.

## Evidenced?

Mostly. Specific counts, caps, parallel-load caveat and output file names are given. Missing: the 240 s hit-0 and 100 s LNA numbers differ from E-139's 27 518 by 4-8x, attributed to a "longer single-process run" without a number for how long; the memory remark for 10^5 nodes is unmeasured (acknowledged). The sizing "10^5 algebras each" is an extrapolation, stated as "plausibly".

## Required for acceptance

1. Say in the Claim that the reverse side restricts to J = 0 steps, and that the reverse column carries no positive control.
2. Either add a reverse positive control (a meet that uses only reverse edges, ideally depth >= 2) or label the reverse column "uninformative until controlled" in the Claim and the table.
3. Explain or drop "smallest-depth-different samples" for hits 0, 4, 13.
4. Note that the control meets at depth 1 and the hit meeting would be at depth 2 vs 11; state what the control does and does not cover.
5. State the wall-clock behind E-139's 27 518 so the 3.6k-6.8k versus 27.5k LNA-side gap is a stated throughput difference, not an apparent inconsistency.
