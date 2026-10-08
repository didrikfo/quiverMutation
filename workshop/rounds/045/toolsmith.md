# An added J = 0 check on a key-guarded walk costs about 8-9% of wall time per step (n = 7 and n = 8, class 1, capped BFS), and 16 of 80 978 (n = 7 c1) and 9 of 79 143 (c2) key-keeping steps fail tiltingPlus

author: toolsmith · round: 045 · kind: tool
thread: T10 (guard audit) · bears on: H-015, E-140, E-145, `search.mutationSearchDepthFirst` docstring
scope: n = 7 classes 1 and 2 (key classes sorted by (size, str(key)), 12 and 14 seed LNAs+duals), n = 8 classes 0 and 1 (2 and 8 seeds); BFS with canonicalKey dedup, key-guarded, capped at 20 000 expansions (n = 7, depth 9, frontier 18 879 / 6 999 left) or 1 200 (n = 8, depth 6-8). Not the library's DFS itself; walks are samples, not classes. Wall times are one machine, 3 jobs on 4 cores.

## Claim

(1) Cost. Adding "refuse the step when J != 0" (or when `tiltingPlus` fails) before the mutation costs 0.41 ms (J) / 0.35 ms (tiltingPlus) per gate-admitted step at n = 7 and 0.76 / 0.70 ms at n = 8, against 5.0 ms (n = 7) and 9.6 ms (n = 8) for the whole step with the key guard. Same walk, 2 500 (n = 7 c1) / 1 200 (n = 8 c1) expansions: 42.0 -> 45.8 s (+9%) and 42.8 -> 46.3 s (+8%) with J; tiltingPlus +9% / +11%. The check refused 2 of 9 112 steps at n = 7 c1 and 0 of 4 827 at n = 8 c1; the refused steps are ones the key guard also drops (key moved), so the walk is unchanged. The check is far cheaper than the key guard itself (step cost 3.9 / 7.6 ms includes mutation, reduction and `coxeterKey`).

(2) Tally. Key-keeping steps that fail `tiltingPlus`: n = 7 c1, 16 of 80 978 mutated steps (20 000 expansions); n = 7 c2, 9 of 79 143; zero in the n = 8 c0 (1 200 exp., 4 992 steps) and c1 (4 827) samples. All of them have J != 0; J = 0 coincides with tiltingPlus on all ~250 000 steps (no J = 0 step failed, no J != 0 step passed). So the docstring is wrong as stated: `coxeterGuard` does not keep the walk to tilting steps at n = 7 c1, c2. These are step counts, not distinct steps (E-145: 13 and 9 distinct J != 0 key-keepers in c1, c2; the 9 in c2 agree). Would refute (2): a rerun with a different walk order giving zero. Not claimed: that any such step leaves the derived class (E-145 says the Cartan matrices are incongruent, so those children are probably outside, but this script does not test it) nor that n = 8 has none (shallow samples).

## Evidence

| walk | exp. | mutated steps | J != 0 key moved | J != 0 key kept (fail tp) | J = 0, tp true | J = 0 iff tp |
|---|---|---|---|---|---|---|
| n 7 c1 | 20 000 | 80 978 | 56 | 16 | 80 906 | yes |
| n 7 c2 | 20 000 | 79 143 | 46 | 9 | 79 088 | yes |
| n 8 c0 | 1 200 | 4 992 | 10 | 0 | 4 982 | yes |
| n 8 c1 | 1 200 | 4 827 | 0 | 0 | 4 827 | yes |

Cost (mode A = key guard as is, with J and tiltingPlus timed on every step but not acting; J, T = the check acts), seconds total / ms per call:

| walk | A total | J mode | T mode | gate | step (mutate+reduce+key) | J per call | tp per call |
|---|---|---|---|---|---|---|---|
| n 7 c1, 2 500 exp. | 49.0 (42.0 without J, tp) | 45.8 | 45.8 | 0.54 | 3.90 | 0.41 | 0.35 |
| n 8 c1, 1 200 exp. | 50.0 (42.8 without) | 46.3 | 47.5 | 1.03 | 7.59 | 0.80 | 0.69 |

Long walks (20 000 exp.): J 0.54 / tp 0.46 ms vs step 5.1 ms (c1); J 0.41 / tp 0.40 vs 4.2 (c2) -- about 9% each. Contention noise: a first run with 6 jobs on 4 cores gave 13 ms/step and is discarded.

Where it goes. In `mutationSearchDepthFirst`, directly after line 435 (`if mutation.mutationIsPossibleAtVertex(pathAlg, vertex):`) and before `quiverMutationAtVertex`: the test needs only the parent and `vertex`, so a refused step is saved the mutation, reduction and key (about 10x the check). A new keyword `tiltingGuard = False` (default off, library behaviour unchanged) forwarded in the recursive call would carry it. Blocker: `perI` and `tiltingPlus` live in workshop scripts (`rounds/033/experimentalist_bothdie.py`, `rounds/001/scholar_h015.py`), not the library; `tiltingPlus` is the one to promote (it returns on the first failing i, `perI` computes all). STEERING q3 left `isTilting` unpromoted, so this is the human's call.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/045/toolsmith_guardaudit.py 7 1 20000 A   # 569 s
timeout 10m .venv/bin/python workshop/rounds/045/toolsmith_guardaudit.py 7 2 20000 A   # 456 s
timeout 10m .venv/bin/python workshop/rounds/045/toolsmith_guardaudit.py 8 1 1200 A     # 50 s (also J, T)
timeout 10m .venv/bin/python workshop/rounds/045/toolsmith_guardaudit.py 7 1 2500 J     # 46 s (also A, T)
```
Run the 20 000 jobs two at a time or fewer; timings in the tables are from 3 concurrent jobs.

## Prior record

E-145 (skeptic_c2.py, 500 s walks, depth <= 30) found the 13 / 9 distinct key-keeping J != 0 steps and that they fail `tiltingPlus` ("every J != 0 step fails tiltingPlus"); this submission reproduces the c2 count with a different walk (BFS by expansion cap), adds per-step counts and the per-step cost, and the n = 8 samples. E-140 (none keeps the key; 278 steps) predates E-145 and is superseded for n = 7 c1, c2. I did not grep for a recorded cost measurement; none was found under "coxeterGuard" in FINDINGS/EXPERIMENTS by the task's identifiers.

## Code changed

None in the library. New script `workshop/rounds/045/toolsmith_guardaudit.py`. No tests run (no library file touched).

## Next

- Chair/human: whether to promote `tiltingPlus` into `quivermutation` and add `tiltingGuard` (default False); then fix the docstring sentence ("what makes the walk a walk in one derived class") to say key-preserving, not tilting.
- Experimentalist (T10 b): rerun the recorded class merges with J-mode on; the cost is ~9%, so it is affordable on every n <= 8 job. Expect no change at n = 7 c0 and n = 8 c0/c1 samples (J != 0 steps there all move the key).
- Skeptic: do the 16 / 9 key-keeping children lie outside the derived class (e.g. via a different invariant than the Cartan congruence)?
- n = 8 c2 and depth beyond 8 were not sampled; the key-keeping J != 0 steps at n = 7 appear only after depth ~6-9.
