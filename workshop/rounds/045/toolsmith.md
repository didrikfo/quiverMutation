# An added J = 0 check on a key-guarded walk costs about 8-11% of a step (n = 7 and n = 8, class 1, capped BFS); key-keeping steps can fail tiltingPlus (16 of 80 978 at n = 7 c1, 9 of 79 143 at c2, all at parent depth 7-8), so key preservation does not imply tilting

author: toolsmith · round: 045 · kind: tool
thread: T10 (guard audit) · bears on: H-015, E-142, E-147, `search.mutationSearchDepthFirst` docstring
scope: key-keeping failures are n = 7 classes 1 and 2 only, from capped BFS of a key-guarded walk with no tilting filter (mode A; mode G on c1 agrees). Classes 1 and 2 (key classes sorted by (size, str(key)), 12 and 14 seed LNAs+duals), n = 8 classes 0 and 1 (2 and 8 seeds); BFS with canonicalKey dedup, key-guarded, capped at 20 000 expansions (n = 7, depth 9, frontier 18 879 / 6 999 left) or 1 200 (n = 8, depth 6-8). Not the library's DFS itself; walks are samples, not classes. Wall times are one machine, 3 jobs on 4 cores; the ratios are the claim.

## Response to referee

Verdict taken: minor revision. Required items 1-5 done. Script `toolsmith_guardaudit.py` gained a depth / taint / distinct tally and mode G (tguard: no expansion into children of a failing step). Reruns (3 jobs on 4 cores, under 10 min): c1 A 463.6 s, c1 G 461.4 s, c2 A 347.5 s.

1. E-086 reconciled (also in *Prior record*). E-086 (0 of ~1.3e6 admitted steps fail tiltingPlus) ran `scholar_walk.py --stop-on-reject`: each walk stopped at its first rejecting level (distance 5-7 at n = 7), before the E-091 rewrite fix (E-087); E-092 re-ran it to the same stopping depth only; its 1.3e6 sums walks of different length. My failures start at parent depth 7-8, beyond that horizon, so E-086's 0 does not extend there. No contradiction. F-038's guard entry ("costs 1.85x", whole guard) added; no recorded cost of an added J check, so the cost row stays new. E-086 class-index agreement is unchecked.
2. Depth reported. Parent BFS level of the key-keeping failures (child one deeper): c1 4 at 7, 12 at 8; c2 6 at 7, 3 at 8. None earlier (cap 20 000, 9 levels). Replaces "depth ~6-9". Both counts are also distinct (parent key, vertex) pairs in this walk; E-147's 13 / 9 distinct children is another measure, no clash.
3. Non-tilting parents. Mode A expands through failing key-keepers, but here it did not matter: 0 of the 16 and 0 of the 9 have a parent reached through an earlier failure ("parent tainted 0"). Mode G on c1 gives the same 16, same depths, 80 979 steps, J = 0 iff tp. At larger caps descendants could be counted; Claim (2) says so. Counts remain step counts.
4. Docstring restated. `search.py` l. 322 claims key preservation / one derived class, not tilting. Header no longer says "wrong". Claim (2): a key-keeping step can fail tiltingPlus, so key preservation does not imply tilting; whether the children lie in the derived class is untested. Next item on rewording is now conditional on the skeptic's test.
5. Timings are machine-dependent (c2 walk: 456 s author, 338 s referee, 347.5 s rerun); ratio is the claim. Rerun per-call: J 0.44 / tp 0.38 ms vs step 4.13 (c1: 10.7% / 9.2%); c2 0.33 / 0.30 vs 3.15 (10.5% / 9.5%). Title now says 8-11%, not 8-9%; wall-time deltas stay single runs.
Not done: n = 7 c1 2 500-expansion and n = 8 walks not rerun (referee reproduced the former); E-086 class-index mapping (needs its logs); out-of-class test (skeptic's).

## Claim

(1) Cost. Adding "refuse the step when J != 0" (or when `tiltingPlus` fails) before the mutation costs 0.41 ms (J) / 0.35 ms (tiltingPlus) per gate-admitted step at n = 7 and 0.76 / 0.70 ms at n = 8, against 5.0 ms (n = 7) and 9.6 ms (n = 8) for the whole step with the key guard. Same walk, 2 500 (n = 7 c1) / 1 200 (n = 8 c1) expansions: 42.0 -> 45.8 s (+9%) and 42.8 -> 46.3 s (+8%) with J; tiltingPlus +9% / +11%. The check refused 2 of 9 112 steps at n = 7 c1 and 0 of 4 827 at n = 8 c1; the refused steps are ones the key guard also drops (key moved), so the walk is unchanged. The check is far cheaper than the key guard itself (step cost 3.9 / 7.6 ms includes mutation, reduction and `coxeterKey`).

(2) Tally. Key-keeping steps that fail `tiltingPlus`: n = 7 c1, 16 of 80 978 mutated steps (20 000 expansions); n = 7 c2, 9 of 79 143; zero in the n = 8 c0 (1 200 exp., 4 992 steps) and c1 (4 827) samples. All of them have J != 0; J = 0 coincides with tiltingPlus on all ~250 000 steps (no J = 0 step failed, no J != 0 step passed). So a step that keeps the Coxeter key can fail tiltingPlus (n = 7 c1, c2), i.e. `coxeterGuard` does not keep the walk to tilting steps; key preservation does not imply tilting. The docstring claims key preservation / one derived class, not tilting, and whether the failing children lie in the derived class is untested here. First failures occur at parent depth 7-8 (4+12 in c1, 6+3 in c2); none of their parents descend from an earlier failure (mode G on c1 gives the same 16), but mode A expands through failing key-keepers, so longer walks could count their descendants. These are step counts, not distinct steps (E-147: 13 and 9 distinct J != 0 key-keepers in c1, c2; the 9 in c2 agree). Would refute (2): a rerun with a different walk order giving zero. Not claimed: that any such step leaves the derived class (E-147 notes incongruent Cartan matrices; this script does not test it) nor that n = 8 has none (shallow samples).

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
Run the 20 000 jobs two at a time or fewer. Seconds are machine- and load-dependent (c2: 456 s author, 338 s referee, 347.5 s rerun); tables are from 3 concurrent jobs; quote ratios. Mode `G` (tguard) is new; modes A and G print depth, distinct and tainted tallies.

## Prior record

E-147 (skeptic_c2.py, 500 s walks, depth <= 30) found the 13 / 9 distinct key-keeping J != 0 steps and that they fail `tiltingPlus` ("every J != 0 step fails tiltingPlus"); this submission reproduces the c2 count with a different walk (BFS by expansion cap), adds per-step counts and the per-step cost, and the n = 8 samples. E-142 (none keeps the key; 278 steps) predates E-147 and is superseded for n = 7 c1, c2. Added after review: E-086 (`scholar_walk.py --stop-on-reject`, 0 of ~1.3e6 guard-admitted steps fail tiltingPlus) stopped each walk at its first rejecting level (distance 5-7 at n = 7), predates the E-087/E-091 rewrite fix, and E-092 re-ran it to the same depth only; the failures here start at parent depth 7-8, so E-086's 0 does not extend there (superseded for n = 7 c1, c2 beyond the stopping depth, as E-147 said). H-015's SUPPORTED line quotes that 0; E-096 (n = 8 c2, depth 8, 0 key-moved) is consistent with my n = 8 zeros. F-038 guard entry: "costs 1.85x" for the whole guard; no recorded cost for an added J check found.

## Code changed

None in the library. New script `workshop/rounds/045/toolsmith_guardaudit.py`. No tests run (no library file touched).

## Next

- Chair/human: whether to promote `tiltingPlus` into `quivermutation` and add `tiltingGuard` (default False); and, only if the skeptic's test puts the children outside the class, reword the docstring sentence ("one derived class").
- Experimentalist (T10 b): rerun the recorded class merges with J-mode on; the cost is ~9%, so it is affordable on every n <= 8 job. Expect no change at n = 7 c0 and n = 8 c0/c1 samples (J != 0 steps there all move the key).
- Skeptic: do the 16 / 9 key-keeping children lie outside the derived class (e.g. via a different invariant than the Cartan congruence)?
- n = 8 c2 and depth beyond 8 were not sampled; at n = 7 the key-keeping J != 0 steps appeared at parent depth 7-8, the deepest levels reached.
