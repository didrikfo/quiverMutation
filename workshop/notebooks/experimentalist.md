# Experimentalist notebook (rewritten each round)

## What I now believe (after round 001)
- At n = 13 all 139 single-cluster cores of `--max-word 4` (`batch._singleCores(4, 6, False)` with a placement) close under the reduced walk; biggest orbit 4217 rows. Whole census ~25 min wall on 4 procs; n = 14 ~ 3x.
- "Own mirror in the orbit" is true of 129/139 cores, so it cannot be half of H-021's "exactly when". Pairing => mirror holds trivially; mirror => pairing fails for 20 cores incl. `3346`, `4056` (both DO hold own mirrors, contradicting H-021's text).
- Only the `30xy` / `330x` words (10 of them) have no mirror; all singleton orbits.
- Reflection fits `o <-> s - o` with overhang d: d=0 62 cores, d=1 28, d=2 17, d=3 2, none 30. Unchecked: d vs |head - tail| of H-020.
- E-052's `45` numbers reproduce at 13 and 14; `3344`, `3346`, `4056` at 14 as well.

## What I tried
- Script `t1.py` in the scratchpad (not saved in the repo): orbits per offset, held offsets, held mirrors, JSONless (stdout lines). Sharded 4 ways; resumable by reading earlier outputs. Pitfalls hit: cores with no placement (skip), variable `S` clash, `pkill -f` kills own shell.
- Monitor "until" loops return immediately in this environment; use a bounded foreground for-loop with sleep 5.

## What I would do next
1. Test d(c) against head/tail from E-051 ledgers (logs/cores-n*.jsonl) for the 109 fitting cores.
2. n = 14 census of all cores (overnight proposal; ~90 min wall on 4 procs), then check that fit class (d) is length-stable.
3. Look at `3033`-type words (interior 0 blocks): singleton orbits, one mirror only.
4. Ask toolsmith (T8) to turn `t1.py` into a `batch.py` task with a test pinning `45` at 13 (`{0,5} 2386, {1,4} 1127, {2,3} 4217, {6} 447`).
- Watch: my fit with big overhang is weak evidence; do not over-read d >= 2.
