# Experimentalist notebook (rewritten round 081)

## What I now believe
- Row convention: n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))`. Check digit count against n.
- 081: E-219's 0.6% counts gate candidates with perI true BEFORE the key test. Mutating them: 3304000 69/69 key-moved; 3060000 98 moved + 1 key-KEPT (parent depth 4, non-LNA parent, kind R v 8, child not in ball). So "no J != 0 step enters the balls" (E-175/E-204) is false for the E-212 move set at n = 9 from parent depth 4; the J = 0 NO is unchanged. Counts reproduce E-219 exactly.
- 081: 7 + 7 join of 3060000/3304000 sized: about 20 000 + 11 400 expansions beyond depth 6 at 5-9 /s (4 shared cores), ~2 h + unsized meet. Not completed. Driver `rounds/081/experimentalist_join97.py run --budget-hours H` (exit 2 when spent, resumable) is untested beyond the budget-0 path.
- 078: n = 10 90000000 / 00000030 tilting-only J = 0: MISS at 5 + 4 (bound). Ball growth x4.3/level, 23 ms/key. n = 15 scan: key 2.06 ms/row; layout n = 14 key table in 4 jobs then 24 shards.
- 075: A_n vs B_n family n = 6..12 not joined by orbit + mirror except n = 5. Coxeter key == class at n = 8, 9; no fixed i rule predicts deletion image (E-194, E-200).

## What I tried
- 081: `rounds/081/experimentalist_jkey.py` (J != 0 candidate mutation + key test), `_join97.py` (driver). Gotchas: running 4 processes drops per-process rate ~2x; the Bash tool blocks `sleep N; cat` chains, use Monitor/until loops; Monitor needs ToolSearch select first.
- 078: `_join10.py`, `_shardplan15.py` (exec'd tiltpath source overwrites ARGV; 10-min `run 5 5` is killed). 075: `_join9.py`, `_delpair.py`. 073/070/067/065 earlier.
- Lessons: grep EXPERIMENTS.md first; mirror needs len(row)+2; forward orbit absent != not equivalent; tee to a file.

## What I would do next
1. Hom-test and walk the one key-kept J != 0 step (parent in ball at depth 4): is the child in the parent class, and does adding it as a move join 3060000 to 3304000?
2. Overnight: join97 with --budget-hours 3, then 3304000 power beyond distance 6 (2400230 puzzle).
3. Pickle n = 10 A-ball for (6,4)/(5,5); is the radius-5 A-ball near the hereditary class?
4. S-1 n = 14 key table Stage A, if approved; items 2, 3 (H-020 head/tail at n = 12 -> 11).
5. Exact classifyLength(10) overnight.
