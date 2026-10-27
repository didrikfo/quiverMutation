# Toolsmith notebook

## Round 003-065 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py` (`--jobs`, `--plan`, `--summary`); ledger unsafe with two processes on one file. Traps: foreground sleeps <= 115 s, `pkill -f` kills own shell, > 3 jobs on 4 cores inflates timings, check `ps` after a timeout.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED; meet up to sign-aware iso. End(T) iso at every step: only Hom(T,T[-1]) tests the premise.

## Round 067-087 (T10/T5)
- E-168/E-198/E-211/E-228: n = 7 failing steps joined by Hom-verified chains; Hom+Cartan has no NO power, the iso meet is the check.
- 083-087: key-kept perI J != 0 steps exist in both F-010 balls (3060000 depth 5: 10/4770; 3304000 depth 5: 1/2958, depth 6: 6/11355), all real, all children rejoin their own class at total 6-8 by J = 0 moves (E-234, E-239, E-244).
- Tools: `085/toolsmith_d5scan.py` (env KEYROOT, DEPTH), `087/toolsmith_d6scan.py` (resumable build/parents/scan). Size jobs by expansions, not nodes.

## Round 090
- Believe: the join97 README misreads its log. The log prints `expansions N None-key children dropped M`; 26 353 / 15 401 are expansions (= ball sizes at depth 6), drops are 9 (A) and 11 (B), re-found by re-expanding all depth-6 parents (19 999 / 11 355). All depth-7 leaves; only a None-key vertex at depth 7 on both sides could hide a total-14 meet; 8 + 11 iso classes, none A~B. The key gap does not limit the n = 9 NO.
- Power of the 7 + 7 NO: nil beyond the 3304000 move component. 2400230 ~ 3304000 is certified yet no meet at 6 + 8 = total 14, which for geodesic balls is the same statement as 7 + 7. YES controls reach total 6-8 only. The NO says nothing on the J = 0 premise.
- Tool: `090/toolsmith_nonekey.py` (scan/report, sharded, resumable). Trap: a script that execs the join97 skeleton must not name a variable ARGV or src (049's preamble overwrites ARGV: silent exit, no output, no error).
- Trap: background jobs with budget close to `timeout 10m` die silently with no output file; use budget 400-470 s (pickle load plus preamble about 60 s). Builds of both depth-6 balls in parallel about 10 min each; scans about 10 parents/s per process.
- Next: count key-kept perI J != 0 steps among 3060000's depth-6 parents (estimate tens, unmeasured); walk J != 0 steps as extra edges and retry 2400230 vs 3304000; calibrate J = 0 distance of equivalent key-equal pairs at n = 8 to learn what radius gives a NO power.
