# At n = 9, L = 5, the back search recovers 3033030 from 3 of 3 cord-bearing monomial tree members (all defect 2), and not at depth 4; no census cell has defect <= 0

author: maverick · round: 061 · kind: negative
thread: T6 · bears on: H-017, E-065, E-078, E-089, E-094, E-162
scope: n = 9 only; forward census depth 5 over the 9 outside LNAs; search-back control at L = 5 for 3 members of ONE LNA (3033030), all defect 2; plus one L = 6 member of the same LNA (also defect 2). No defect-1 source was used. Monomial tree quipus. Nothing at n >= 10, nothing at depth >= 7 run.

## Response to referee

1. Done: 32 cells; all nine defect-1 cell types listed in Evidence.
2. Not as asked: `maverick_cordcontrol_L6.txt` in the repository holds a finished L = 6 control, not only the header (line 2: `3033030 cords 3 rels 5 hit True neg True 376s`, then the table and `{'hit': 1, 'neg_ok': 1}`); the referee likely read it mid-run. I kept it, stated as ONE member of ONE LNA, defect 2, and dropped the "attempt" wording. This also does part of item 6 (L = 6 on 3033030, one of 14 candidates, 376 s); the defect-1 source was not run.
3. Done: title and Claim narrowed; all control members defect 2, no defect-1 source; round-trip-by-construction caveat added.
4. Done: E-103, E-094 Limits and E-065 cited as restatements under Prior record.
5. Done: Coxeter guard ON (default), stated in Reproduction.
6. Defect-1 source control not run (needs a new selection by defect across LNAs, over my budget this round); left to Next.

## Claim

(1) The search behind H-017 (`qr.reachedQuipuAlgebras` forward, `search.linesReachedFrom` back) already records the relation set of every tree-quipu member, so "a search that sees relations" exists; E-162 is a statement about invariants, not about the search. (2) The "no positive MONO control" worry attached to H-017 (E-089, E-094) is about *cycle* members (underlying quiver with an undirected cycle, always with a sum relation; GLOSSARY "Cycle member"). H-017's candidates are *trees* (quipu cords = branches) with monomial relations, so the MONO filter and those controls do not bear on them. (3) (round trip by construction: sources are forward-reached from 3033030, so this tests search completeness on reachable members, not detection of a defect <= 0 member) The control that does bear on them: take a forward-reached tree quipu with >= 1 cord and monomial relations at recorded path length L = 5, rebuild it as a path algebra, search back: the source LNA is returned at depth 5 and not at depth 4 in 3 of 3 (cords/rels 3/5, 3/5, 2/4, source 3033030, n = 9); and at L = 6 (the depth of E-078's negatives) 1 of 1: found at depth 6, not at 5 (cords 3, rels 5; 376 s including a 174 s depth-6 forward walk). The depth-6 forward walk of 3033030 still has min(rels - cords) = 2, as at depth 5. (4) Forward census at depth 5, all 9 outside LNAs: 0 members with relations <= cords; relations = cords + 1 is reached already at path length 2 (1 cord, 2 rels), 3 ((2,3)) and 4 ((3,4)), so the wall at defect >= 1 is reached cheaply and defect 0 is not, through depth 5. That is E-065's observation re-tallied with cords >= 1 only; it is not a new result and does not support a law beyond depth 5 (E-065 already saw minima fall with depth: 4444400 from 2 at depth 5 to 1 at 6).

What this does NOT claim: that H-017 is true; that depth 6 or 7 members are controlled (my control is L = 5 on one LNA); anything for n >= 10. A refutation would be a tree quipu with relations <= cords found in an outside class at any depth.

## Evidence

Forward census, n = 9, depth 5, 9 outside LNAs (532 s total, per LNA 38-88 s). Members (non-hereditary) per LNA 44-81, of which L = 5 with cords >= 1: 18, 39, 37, 41, 32, 30, 39, 41, 37. Table (path length, cords, rels): count is in `maverick_cordcontrol_plan9.txt` (32 cells; counts are certificate counts summed over LNAs, mirror pairs doubled). Defect-1 cells (path length, cords, rels): count: (2,1,2) 2, (3,1,2) 2, (3,2,3) 4, (4,1,2) 8, (4,2,3) 16, (4,3,4) 2, (5,1,2) 8, (5,2,3) 26, (5,3,4) 2 (summed over LNAs; per-LNA split not tabulated); no cell of defect <= 0.

Positive control, 3033030, L = 5, three members taken nearest the diagonal first (the sort puts defect 2 first, because no defect-1 member sits at L = 5 for this LNA):

| cords | rels | hit at depth 5 | hit at depth 4 | time |
|---|---|---|---|---|
| 3 | 5 | yes | no | 77 s |
| 3 | 5 | yes | no | 112 s (cumulative) |
| 2 | 4 | yes | no | 168 s (cumulative) |

So the back search is not blind to tree quipus with cords: it finds the source exactly at its path length, as E-071/E-074 showed without cords being tracked (they did not record cords for the members; this closes that).

Sizing for a design that would test more (from this round's timings and E-075/E-078: x5.5-5.7 per level): forward census depth 6 about 4 min per LNA, depth 7 about 22 min per LNA (9 LNAs, one per shard, about 3.5 h in 3 cores); back search from the 16 K = 4 candidates at depth 7 about 27 min each (E-074). Both are OVERNIGHT-size and neither can *confirm* H-017: a universal claim over an infinite class has no finite search proof, only refutation.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/061/maverick_cordcontrol.py 9 5 5 2            # census + table, 532 s (use run_in_background)
timeout 10m .venv/bin/python -u workshop/rounds/061/maverick_cordcontrol.py 9 5 5 3 1 --run  # control on 3033030, 168 s
```
```
timeout 10m .venv/bin/python -u workshop/rounds/061/maverick_cordcontrol.py 9 6 6 1 1 --run   # L = 6 control, 376 s
```
Coxeter guard: ON in `--run` (the script calls `search.linesReachedFrom` with defaults, `coxeterGuard = True` in `mutationSearchDepthFirst`).
Output: `maverick_cordcontrol_plan9.txt`, `maverick_cordcontrol_run.txt`, `maverick_cordcontrol_L6.txt`.

## Prior record

Restatements: claim (2) is E-103 (line 678: a "cord member" is not the GLOSSARY quipu cord) and E-094's Limits; claim (4) is E-065's depth 4-5 tally for all nine LNAs and depth 6 for 3033030/4444400 (minima fall with depth). Genuinely new: only the cords >= 1 control (3033030, L = 5, and the single L = 6 member).

H-017 header lists E-094 ("still no positive MONO control") and E-089 as open controls; E-094's own Limits say the candidates (monomial quipus with cords) and the controls (sum-relation cycles) "differ in kind". GLOSSARY now separates "cycle member" from "quipu cord". E-071/E-074/E-084 controlled the search on tree members but did not tabulate cords. E-065 has the forward (cords, relations) pairs; E-078 the 16 K = 4 backward negatives at depth 6. E-162 says the signature is Cartan-determined. E-154 (Cartan congruence) and the Smith-form filter already in `rounds/004/maverick_verify.py` mean no cheap invariant is left that sees relations: any Cartan-level invariant is blind to (cords, relations).

## Code changed

New: `workshop/rounds/061/maverick_cordcontrol.py` (adapted from rounds/009 `maverick_control5.py`, adds cords tally and cords >= 1 selection). No library code. No tests touched.

## Next

Proposal for the ledger (conference 064): do not close H-017 (it is not refuted), but retitle T6 "parked: search-only, refutation-only"; the header lines about "no positive MONO control" should be corrected to say the cycle-member controls are off topic and point to this control (chair/skeptic decide). Reopen only with (a) a proof route (theorist: a derived invariant that counts relations minus cords on a tree, e.g. from the Euler form plus tree structure; E-065 says gldim <= 2 gives only >= 1 or 2), or (b) an OVERNIGHT run: forward census depth 7, one LNA per shard, looking for any cell with defect <= 0. Experimentalist: run a defect-1 source (cells (5,1,2), (5,3,4), other LNAs) and widen the L = 6 control (one member of one LNA so far) to the other outside LNAs and to members nearer the diagonal (defect 1), which exist only at L <= 5 for 3033030.
