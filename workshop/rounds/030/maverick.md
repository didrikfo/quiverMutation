# The free-end deletion rule (K >= 3) survives the corrected labels at n = 8, 9, 10 and fails at n = 11; K >= 4 holds at every n

author: maverick · round: 030 · kind: result
thread: S-1 · bears on: E-112, E-115, H-020, F-028

## Claim

Rule: delete the end vertex of a free end (all vertices left of every relation, or all right of every relation) of K free vertices; the image class depends only on the source class. With the E-115 class labels (the 16 n = 9 and 176 n = 10 cospectral LNAs placed by the F-047 profile, none dropped; resolved classes are now 16 / 36 at n = 9 / 10 instead of E-112's resolved subset) the K >= 3 statement **survives at n = 8, 9, 10** (3/3, 6/6, 12/12 source classes, 26 / 82 / 262 ends; E-112 had 3, 5, 10 classes and 4 coarse '?' images, now none). **At n = 11 it fails**: 1 of 21 source classes with a K >= 3 end has two image classes (34 / 60 ends; different n = 10 Coxeter keys, so the split is real). The failing source class lies in two free+edge+double orbits, and the images split already inside the larger orbit (943 members: 34 + 48 ends), so it is not a coarse-source artefact. **K >= 4 holds at every n = 8..11 (2/2, 3/3, 6/6, 10/10 classes; 252 ends at n = 11); K >= 5 too.** So the threshold is not 3 uniformly: K >= 3 at n <= 10, K >= 4 needed at n = 11 (the evidence suggests "K >= 3" was a small-n coincidence; it does not say the threshold is 4 for all n or that it grows with n).

Not claimed: a derivation; coverage of the 418 of 442 n = 11 LNAs unresolved by the profile (see Limits); that K >= 4 is a theorem at n = 11.

## Evidence

Image class is the corrected label of `piecewiseHereditary.removeVertex` at the end vertex. Classes with at least one end of free length >= K0, head and tail pooled:

| n | K0=1 | K0=2 | K0=3 | K0=4 | K0=5 |
|---|---|---|---|---|---|
| 8 | 8/10 | 6/6 | 3/3 | 2/2 | 1/1 |
| 9 | 8/16 | 10/11 | 6/6 | 3/3 | 2/2 |
| 10 | 21/36 | 21/23 | 12/12 | 6/6 | 3/3 |
| 11 | 32/69 | 32/41 | **20/21** | 10/10 | 5/5 |

(transported / classes having an end of free length >= K0.) K0 = 2 failures: 0 / 1 / 2 / 9 classes at n = 8 / 9 / 10 / 11 (E-112: 0 / 1 / 2 at n = 8..10).

n = 11 failure: class key (1,1,0,-1,-2,-3,-3,-2,-1,0,1,1), 1305 LNAs. Images: key (...,-2,-2,-2,...) 34 ends, key (...,-2,-3,-2,...) 60 ends. By orbit: orbit 15107 (943 LNAs) 34 + 48, orbit 15035 (362 LNAs) 12 ends, all to the (-3,-2) image. Free-run length K = 3 and 4 appear; K0 = 4 removes the failure, so the K = 3 ends carry it (read from the K0 runs, not a per-end table).

n = 11 labels: 418 of 442 cospectral-unresolved LNAs remain unresolved after the profile (E-115) and are dropped as sources; images are all at n = 10 (fully resolved), so no image is dropped. The n = 11 failure is therefore not a labelling artefact on the image side; on the source side, key-only classes could merge two derived classes, but the images split inside one orbit, whose members are equivalent by the free/edge/double moves, so the failure stands.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/030/maverick_endstrip2.py 8   # 25 s each at 8, 9, 10; 60 s at 11
timeout 10m .venv/bin/python workshop/rounds/030/maverick_endstrip2.py 11
timeout 10m .venv/bin/python workshop/rounds/030/maverick_fail11.py          # the n = 11 failure, ~2 min
```

## Prior record

E-112 (rule at n = 8, 9, 10 with the unresolved dropped), E-115 (profile placement; n = 11 gap). New: the rerun with corrected labels (statement stands), and the n = 11 counterexample at K = 3 with K >= 4 holding. Grep of `research/` found no other statement of the rule.

## Code changed

New `workshop/rounds/030/maverick_endstrip2.py` (imports `maverick_classes`, `toolsmith_snfresolve`) and `maverick_fail11.py`. No tests; no library code touched.
