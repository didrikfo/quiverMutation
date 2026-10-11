# The n = 11 K = 3 failure is a room-to-move effect: the same core word goes to image I1 when deleted from a free run of 3 and to I2 when deleted from a run of 4

author: maverick · round: 034 · kind: result
thread: S-1 · bears on: E-114, E-117, E-120, H-020, H-018

## Claim

Level: tested on small cases (one class, n = 11 only). In the single failing class of E-120 (key (1,1,0,-1,-2,-3,-3,-2,-1,0,1,1), 1305 LNAs, two orbits 15107 (943) and 15035 (362)), 82 ends have free run K >= 3 (66 with K = 3, 16 with K = 4; none K >= 5). Head and tail ends give identical counts (mirror). The image class is a function of (K, core word read away from the deleted end): 77 such keys, 0 with two images. The split is not between sources of different cores: 12 core words occur at both K = 3 and K = 4 and each goes to I1 at K = 3 and to I2 at K = 4. All K = 4 ends go to I2. So the failure is "K = 3 is one short": deleting one vertex leaves a free run of 2 at n = 10, the regime where E-120 already finds K = 2 failing (the core then has no room to move), whereas K = 4 leaves 3. This is S-1 question 3 realised once, as a table; it does not claim a threshold law for n = 12 or a derivation. It is not a labelling artefact on the image side (the two images have different Coxeter keys, E-120) nor on the source side (both images occur inside orbit 15107, whose members are related by free/edge/double moves).

## Evidence

Images (n = 10 classes): I1 key (1,1,0,-1,-2,-2,-2,-1,0,1,1); I2 key (1,1,0,-1,-2,-3,-2,-1,0,1,1).

| K | side | orbit | image | ends |
|---|------|-------|-------|------|
| 3 | head | 15107 | I1 | 17 |
| 3 | head | 15107 | I2 | 16 |
| 3 | head | 15035 | I2 | 5 |
| 3 | tail | 15107 | I1 | 17 |
| 3 | tail | 15107 | I2 | 16 |
| 3 | tail | 15035 | I2 | 5 |
| 4 | head | 15107 | I2 | 8 |
| 4 | head | 15035 | I2 | 1 |
| 4 | tail | 15107 | I2 | 8 |
| 4 | tail | 15035 | I2 | 1 |

(Orbit 15035 never reaches I1; orbit 15107 reaches I1 only at K = 3.)

By core length L (last relation end - first relation start + 1) at K = 3: L = 4: I1 2 / I2 0; 5: 2/2; 6: 4/4; 7: 8/8; 8: 18/28. At K = 4: I2 only, L = 4..7 (2, 2, 4, 10). Counts are dominated by the free run at the other end (K = 3 with 0 free at the other end: 18/28; 1: 8/8; 2: 4/4; 3: 2/2; 4: 2/0), so L is not what separates I1 from I2.

What separates them is the core word (the a_i string, relation i -> i + a_i). Oriented with the deleted end on the right and stripped of zeros, K = 3 words going only to I1 (25 of them) all end in a 3 (or are the single relation of length 7), e.g. 3, 23, 203, 223, 2003, 2023, 2223, 22223, 200203; those going only to I2 (37) include 2, 4, 5, 6 as last letters or a doubled 33: 32, 302, 333, 334, 2302, 22302, 62, 63, 66. The 12 words 3, 23, 203, 223, 2003, 2023, 2203, 2223, 20003, 20203, 22003, 22203 are in both K sets, K = 3 -> I1, K = 4 -> I2 (the orientation reversal for head ends was done on the digit string, a guess, but it gives 0 conflicts and head/tail counts agree exactly).

Reading: sources at K = 3 and K = 4 with the same core are in one class at n = 11 (they differ in the other end's run, a free move). After deleting, the n = 10 LNAs are (c, run 2) and (c, run 3): different classes, since with a run of only 2 the core cannot be slid to where I2 lives. That is exactly the "no room to move back" scenario of S-1 question 3, but at the image side only; I did not find a source whose class needs room that deletion removes in the other direction (simple LNA obtained only via room).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/034/maverick_endtable.py   # about 80 s
```

(imports `workshop/rounds/030/maverick_endstrip2.py` and `workshop/rounds/027/maverick_classes.py`.)

## Prior record

E-120 states the failure and that "K = 3 ends carry it ... read from K0 runs, not tabulated"; this tabulates it and adds the K = 3 vs K = 4 same-core comparison and the function-of-core-word statement. E-114 K = 2 failures are the same phenomenon one step down (not checked per end here). grep of research/ for "room to move" finds only H-018/H-020 context, no per-end table.

## Code changed

None (new script only). No tests run.

## Next

- Test the prediction at n = 12: the first failing K should be 3 again for the n = 11 failing core words (a length-12 class containing them), or 4 if the room needed grows; orbit-only (label-free) test, a source orbit whose K = 3 ends split.
- Compute for the 9 K = 2 failing classes at n = 11 whether they too split by (K, core word) with the same-core K = 2 vs 3 pattern; if yes, state "image is a function of (core word, K) and changes class exactly when the image run drops to 2".
- theorist: why does the last letter 3 (a relation of length 3 ending at the core's edge) put the K = 3 image in I1; derive from the H-020 rule table.
- skeptic: the head-end orientation of the digit string; redo with freeMoves.mirrorRow.
