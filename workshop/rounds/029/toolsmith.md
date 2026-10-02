# The 16 unresolved n = 9 LNAs (and the 176 at n = 10) are all placed by the F-047 Smith profile: 16/16 into P^(1,4)_(1,0,1); 176/176 into one of the 4 quipu classes, in 15 s

author: toolsmith · round: 029 · kind: result (tool; partly rediscovery)
thread: S-1 validation · bears on: F-047, F-010, E-112, H-003

## Claim

Using the F-047 profile as the independent derived invariant (Smith normal form over Z of g(Phi) for each irreducible factor g of the Coxeter polynomial, plus SNF of C + C^T), every one of the 16 LNAs at n = 9 that `maverick_classes.classes` leaves '?' has the profile of P^(1,4)_(1,0,1) and not that of P^(1,2)_(1,1,2) (the other member of the cospectral key). At n = 10 all 176 get a unique class: 32 P^(1,5)_(1,0,1) and 104 P^(2,3)_(1,1,1) (key 1,1,-1,-3,-4,-4,-4,...), 16 P^(1,4)_(1,0,2) and 24 P^(3,3)_(1,0,1) (key 1,1,-1,-4,-7,-8,...). No LNA matched no class, none matched two.
Does not claim: that an LNA with profile of class A is *in* A (the profile is necessary, not sufficient: a third class in the same key with the same profile cannot be excluded by this alone; but the quipu theorem says an LNA class contains a quipu, and in these keys only the two quipus per key exist). It does not use the move orbits: it is an exclusion of the other quipu class, so the placement rests on "each LNA of this key lies in one of the key's quipu classes". It does not touch the n >= 11 cospectral LNAs.

## Evidence

Soundness: on every resolved LNA of the cospectral keys the profile is constant per class (all 2 classes at n = 9 and 4 at n = 10: 1 profile each, so "sound"); the two profiles in a key differ (the placement is unique for all 192 LNAs). The unresolved LNAs are all distinct as Cartan matrices (16, 176), so no caching helped; one profile costs about 0.03 s.

| n | unresolved | placed | P^(1,4)_(1,0,1) | P^(1,2)_(1,1,2) |
|---|---|---|---|---|
| 9 | 16 | 16 | 16 | 0 |

| n = 10 class | LNAs |
|---|---|
| P^(2,3)_(1,1,1) | 104 |
| P^(1,5)_(1,0,1) | 32 |
| P^(3,3)_(1,0,1) | 24 |
| P^(1,4)_(1,0,2) | 16 |

n = 9 members: 2223030, 2203030, 2023030, 2003030, 3030000, 3030200, 3030220, 3030222, 3030202, 3030020, 3030022, 3030002, 0223030, 0203030, 0023030, 0003030 (all P^(1,4)_(1,0,1)); the other 'unresolved' 3030000 x-family is thus one orbit-and-its-mirror family that the join of free/edge/double/mirror orbits did not connect to a seed.

Sizing (--plan): n = 9 0 s, n = 10 4 s estimate; actual 3.7 s and 11.7 s including the class tables. A mutation search was not needed. (The F-045/E-076 route, depth-6 searches of 25-30 min per candidate and `reached []` so far, would have been the expensive alternative.)

## Reproduction

```
.venv/bin/python workshop/rounds/029/toolsmith_snfresolve.py 9 --plan    # sizing, 1 s
.venv/bin/python workshop/rounds/029/toolsmith_snfresolve.py 9           # 4 s
.venv/bin/python workshop/rounds/029/toolsmith_snfresolve.py 10          # 12 s
```
It imports `workshop/rounds/027/maverick_classes.py` for the label set and unresolved list.

## Prior record

F-047 already says the profile splits the n = 9 pair, and names `3030000` (8 rows, with 3060000, 6000030, 2223030) as P^(1,4)_(1,0,1) (the orbit of 3030000 sits in F-047's table). The n = 9 result here is therefore a rediscovery, and consistent with F-047: the new content is that the profile is *computable on every LNA*, so the maverick gap is not a gap in the classification but in `classes()`: its orbit join (free + edges + doubles + mirror) does not link those 16 to a quipu seed, while the profile does. F-047 states the n = 10 split as "3 of 25 groups" at the orbit level; the 176 here are the LNAs of two such groups. I did not find the n = 10 per-class counts recorded (grep of 3030000 / the keys gave only F-047, E-076 context). Retraction check: nothing in RETRACTIONS.md touches F-047 that I found (not searched exhaustively).

## Code changed

None in `quivermutation/`. New: `workshop/rounds/029/toolsmith_snfresolve.py`. No tests run (no library change). Weakness: the profile is computed with sympy (slow for n >= 12; 0.02 s at n = 10).

## Next

- maverick/experimentalist: patch the labelling by adding a fallback to the profile for '?' LNAs (the script's `profile`) and redo the E-112 tables with all 1430 and 4862 LNAs; the K >= 3 count would then be over complete classes, especially n = 10 ('?' ends).
- Why the orbit join misses these 16 (theorist): a family 3030000 ... 0003030 of 16 with a '3 0 3 0' core; is it one orbit under a move not yet in `freeMoves` (H-003: orbits that ought to merge)?
- toolsmith: add the profile as an opt-in column of `coxeterTables`, with a test (pins n = 9 16 -> P^(1,4)_(1,0,1)); extend to n = 11 (cospectral groups there: 4) if wanted.
