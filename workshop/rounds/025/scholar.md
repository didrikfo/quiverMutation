# At n = 8 the "reject iff long square" regularity breaks: the out-degree 2 rejects are D-type (commute into one arrow, killed into the other by zero relations) and all fail the Cartan test through the real rewrite

author: scholar · round: 025 · kind: negative · thread: T5 · bears on: H-015, E-066, E-097, E-100

## Response to referee

1. Walk table, n = 8: added (below). The walk statement "no D or G type algebra is reached" and the step-by-step iff "J != 0 <=> out-degree 1 and long square on alg.rels" are
   restricted to n = 5..7, class 0, 300 s caps. At n = 8 they are false: 55 of 59 J != 0 steps have out-degree 2 and no long square (referee saw 42 at 200 s; counts are cap-dependent).
2. The out-degree 2 rejects, classified (`scholar_n8rejects.py`, 300 s, 12 971 algebras): 61 distinct (parent, v) pairs (canonical key of parent + v), ALL out-degree 2.
   By my test (a relation of >= 2 paths all ending ..,v,e_beta, for each arrow beta out of v): every one commutes into exactly 1 of the 2 arrows ("D-part 1/2"); none into both, none into neither.
   So not G (not a presentation artefact: J != 0 is computed from the ideal, not from `alg.rels`), and not hand-case D. It is a third type, call it D': commutes into one arrow,
   and into the other the element c is killed by zero relations. First reject printed (v = 3, arrows 3->6, 3->8):
   `arrows (1,3)(2,4)(3,6)(3,8)(4,5)(5,1)(5,7)(6,8)(7,3)`; rels `1-3-8 = 0; 2-4-5-7 = 0; 4-5-1 = 0; 4-5-1-3 = 4-5-7-3; 5-1-3-6 = 5-7-3-6; 5-7-3-6-8 = 0; 7-3-6-8 = 7-3-8`.
   c = 5-1-3 - 5-7-3: c.(3,6) = 0 by the long square into 6; c.(3,8) = 5-1-3-8 - 5-7-3-8 = 0 - 5-7-3-6-8 = 0 (using 1-3-8 = 0, 7-3-8 = 7-3-6-8, 5-7-3-6-8 = 0). Checked by hand for this one only.
   Not a single path (gate admits); dim J = 1 (57 cases) or 2 (4 cases).
   Cartan test through the actual rewrite (`quiverMutationAtVertex` + `reducePathAlgebra`, compare Cartan(child) with R C R^T as rounds/018): all 61 fail (cong False), all tiltingPlus False. Consistent with r018/r021.
3. Next/chair wording rewritten: "long-sided square" is a shape of reached rejects at n <= 7 only; at n = 8 the intrinsic statement (J != 0 with no single path in J) is all that survives.
4. Claim line states class and caps.

## Claim

On guided class-0 walks (gate-admitted steps, 300 s caps per size) "J != 0 <=> out-degree 1 and a long square on alg.rels" holds step by step for n = 5, 6, 7 and FAILS at n = 8:
59 steps with J != 0, of which 55 have out-degree 2 and no long square; all 61 distinct out-degree 2 rejects found commute into one outgoing arrow only and fail the Cartan congruence through the real rewrite.
Only class 0, only these caps, counts are lower bounds on distinct objects and cap-dependent. Does NOT claim: classes 1-3, n >= 9, or that D' is the only extra type (my classifier tests only the "commute into arrow" shape; 61 are all D-part 1/2 but I read only one by hand).

## Evidence

Walk table (class 0, 300 s caps, steps; `scholar_longsquare_n{5,6,7}.txt` in rounds/023, `scholar_longsquare_n8.txt` here):
| n | steps | J != 0 | out-degree 1 and long square | J != 0, out-degree 2 | J = 0 with long square |
|---|---|---|---|---|---|
| 5 (6 240 algebras) | 16 620 | 0 | - | 0 | 0 |
| 6 | 101 045 | 1 139 | 1 139 | 0 | 0 |
| 7 | 66 049 | 162 | 162 | 0 | 0 |
| 8 (12 368) | 20 664 | 59 | 4 | 55 | 0 |
Referee's n = 8 run (200 s, 8 988 algebras) had 2 + 42, mine 4 + 55: the ratio is stable. The hypotheses of r023 (2)(a) "out-degree 1" is therefore not an invariant of reached rejects. Mechanism: at n = 8 the guarded
walk reaches algebras where zero relations make c die into the second arrow; at n <= 7 none was reached (the walk does not cover everything: 300 s caps).
Rejection table (61 distinct): (D-part 1/2, out 2, dim J 1, tiltingPlus False, Cartan cong False) 57; (same, dim J 2) 4.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/023/scholar_longsquare.py 8 --class 0 --budget-sec 300   # 300 s, writes the n = 8 row
timeout 10m .venv/bin/python workshop/rounds/025/scholar_n8rejects.py 8 --class 0 --budget-sec 300    # 300 s, classification + Cartan test; output scholar_n8rejects_n8.txt
```

## Prior record

E-066, E-097, E-100 as in r023; E-100's "every rejecting step has the long-sided shape" was counted on n <= 7 and is not contradicted there, but is contradicted at n = 8 (new). Nothing in RETRACTIONS.
The D/E/G split of r023 gains D' (one arrow by commutation, one by zero relations). grep of research/ for "out-degree" finds nothing at n = 8.

## Code changed

New `workshop/rounds/025/scholar_n8rejects.py`; no tests touched, nothing in `quivermutation/`.

## Next

- toolsmith: n = 9 class 0 walk with the reject classifier (needs > 10 min; overnight proposal, add a checkpoint). Does D' stay the only out-degree 2 type?
- theorist: the exact statement on walks is J != 0 with no path in J; prove or refute that D' arises from guarded steps only when zero relations enter.
- chair: for E-084/E-095/E-100 say "long-sided square" for n <= 7 class 0 only; do not state the iff for n >= 8.
