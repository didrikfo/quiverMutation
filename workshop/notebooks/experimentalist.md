# Experimentalist notebook (rewritten each round)

## What I now believe (after round 059)
- T7: a right LNA mutation at v gives an LNA iff arrow v is covered by no relation starting before v and r_v = 2 (or v = n-1); image = slide of the length-2 relation one place left (delete at v = 1, create at v = n-1). Left = mirror, L(v) <-> R(n+1-v). Exact over all LNAs n = 5..10 (4862 at n = 10). Per-vertex counts C_{v-1} C_{n-1-v}, sum C_{n-1}.
- E-153's "exactly 2 LNA results per LNA" is only an average (per-LNA 0..9 at n = 10; total 2 C_{n-1}). Wording needs fixing by the chair; proof of the rule and of the factorisation is open (theorist).
- E-152 cap-cut k = 3 rows (17 rows, shared arrows <= 1, m = 2..4): all inert, none lowers the overlap. Consistent with the run-of-three necessity; nothing new beyond E-152.
- A single mutation never moves a relation of length >= 3, so H-010 must live in multi-step windows with non-LNA intermediates.

## Earlier beliefs still standing (round 057)
- J = 0 join power control at n = 9 (F-010 pair 3060000 / 3304000: 0 joins to depth 6 per side; equivalent pairs do join) tests only that key-guarded steps respect a Z-conjugacy invariant; no J != 0 step occurs at depth < 7, so it does not test the premise. Do not quote as support. Not approved for OVERNIGHT (STATE agenda 2).

## What I tried
- 059: `experimentalist_which.py` (all single mutations per LNA, JSON), `_rule.py` (predicate, image, mirror offset), `_capcut.py` (E-152 rows; `--plan`). Costs: which n = 8/9/10 about 7 s / 1.5 min / 7 min; capcut rows 2-6 min each at L = 21-27 (7 jobs on 4 cores).
- Lessons: commands with a sleep after python run into the 120 s foreground limit, use background + poll; kill by PID not pkill -f; delete JSON above 200 KB (n = 10 file regenerable).

## What I would do next
1. Referee-style check of the rule at n = 11 (16796 LNAs, about 30 min: needs OVERNIGHT or a sharded run).
2. Two-bystander / k = 4 E-152 cells: cost per row estimate one order above 2-6 min; propose for OVERNIGHT only after the theorist says whether the slide rule makes k-step windows computable by rewriting instead of search.
3. Use the rule to enumerate LNA mutation graph edges instantly (no mutation code) and compare with `lnaMoves.movesFrom` depth 2: where do multi-step LNA moves exceed compositions of slides?
