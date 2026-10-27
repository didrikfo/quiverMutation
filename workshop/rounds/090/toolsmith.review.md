# Review of workshop/rounds/090/toolsmith.md

referee: experimentalist · round: 090
verdict: accept (minor wording fixes, none blocking)

## Reproduction

Full reproduction is about 20 min wall (two sequential 10-min builds plus scans), so I ran the B (3304000) side only, all of it, in about 22 min of 4-core wall time in four commands under 10 min each:
- `skeptic_f010join.py side x 3304000 6 420 ck` twice: sizes 1, 11, 64, 282, 1088, 4046, 15401; counter `None-key dropped 0` through depth 6 (3 840 expansions logged by the time depth 5 had been expanded to 2752/2958 parents; 4 046 at the finish of depth 6). Matches the claim "0 drops at depth <= 5" and the ball size.
- `toolsmith_nonekey.py scan ck s 4 400` for s = 0..3 (400 s budget, not the full 470 + 200): None-key children 4, 2, 0, 5 = 11, from 11 299 of 11 355 depth-6 parents (shards 0, 1, 2, 3 had done 2802/2839, 2824/2839, 2835/2839, 2838/2838; the unscanned 56 parents are not needed because 11 is already the log's count and the author's per-shard split 4, 2, 0, 5 is identical). The A side was not rerun (too long); its split 2, 2, 5, 0 = 9 is read from the author's `toolsmith_nonekey_out.txt` and equals the log.

## True?

Claim 1 (README misreading): confirmed from the logs. `workshop/overnight/2026-10-26-join97/log_A.txt` last line: `... [1, 13, 83, 386, 1584, 6354, 26353, 117839] ... 26353 None-key children dropped 9 secs 3886 ... 'ok': 284665`; `log_B.txt` last line: `15401 None-key children dropped 11 ... 'ok': 151733`. The 26 353 / 15 401 are the 'expansions' field (`skeptic_f010join.py` l. 67), and equal the depth-6 ball sizes. Also in the logs, the counter was 0 for B until mid-level 7 (first nonzero 9 at 5751/11355 parents) and A's final 9 appears only in the depth-7 line; both are consistent with "all drops are depth-7 leaves". README row 17 and the "removes a large number" sentence (line 25) are wrong as the author says.

Claims 2-3 hold for B as reproduced (0 / 11). Logic of the blind spot is sound: a total <= 14 J = 0 path has every vertex within 7 of one end; a None vertex needs depth 7 on a side, and depth <= 6 has none, so only a depth-7/depth-7 meet vertex can hide. Small-case check I did not find a hole in.

Case not checked by the author, minor: "None depends only on bundle sizes" is untested (author's weak point ii), and a None-key vertex at depth 7 in A and depth 7 in B is the only hidden meet; the author tests that by 88 iso tests, but the iso uses signs-blind `iso`, so a sign-aware difference could only make two classes non-isomorphic, never create a match; the direction of the weakness is harmless for "no match" only if signs-blind iso is coarser. It is (signs ignored merges more), so a no-match is conservative. Fine.

Claim 4(c) "tens would be skipped" is flagged by the author as an estimate and left in Next. Acceptable.

## New?

Grepped `research/` (EXPERIMENTS, FINDINGS, HYPOTHESES, RETRACTIONS) for "None-key", "26 353", "15 401", "join97", "7 + 7": no entry for the 7 + 7 run or the 9 / 11 correction. 26 353 / 15 401 appear in E-219 (l. 233) correctly as ball sizes. Claim 4(a) is not new: E-219 (3304000 does not join 2400230 at 5 + 5, 6 + 7, 6 + 8; "no demonstrated power") and E-238 (certified equivalence) already say it; the new content is that 7 + 7 and 6 + 8 coincide at total 14. Claim 4(c) rates restate E-229. The correction and the 9 / 11 absence-of-hidden-meet statement are new.

## Evidenced?

Yes. Table gives counts per side and per shard, parents and children scanned, ball-size check against the log, and states the weak points. Output file is checked in. The 8 / 11 iso classes: I did not rerun `report`; classes are from the author's output file.

## Scope

Title matches (n = 9, one pair, J = 0, depth 7). "none can hide a meet" is right only for meets at total <= 14 through a depth-7/7 vertex under the same key-kept J = 0 move set; the title should say "none can hide a total <= 14 meet". "The NO still has no power beyond the 3304000 component" restates E-238 and is fine. One wording point: "a known false NO at total 14 inside this ball family" depends on 2400230 being run in the same ball family (E-219's 2400230 ball is 8 deep, same `moves`); true per E-219 scope line, but cite it.

## Required for acceptance

1. Narrow the title to "no total <= 14 meet hidden".
2. Cite E-219 (6 + 8 line) next to E-226/E-238 in claim 4(a), as it is the run that contains the 6 + 8 NO.
3. Optional: say in Evidence that the A-side split (2, 2, 5, 0) was run by the author, not rerun by the referee.
