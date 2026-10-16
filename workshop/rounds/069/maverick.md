# T6 defect-1 source control passes (5 of 5 run: 3345000 x3, 3455000 x1 referee, 3505000 x1 with 3 cords; n = 9, L = 5); T6 may park as search-only

## Response to referee

Verdict was minor revision; all five required items done, plus the (5,3,4) cell and a six-LNA sizing.

1. Member identities: done. `maverick_defect1.py` now prints path, arrows and relation set (0-based vertices) for every selected member. Files: `maverick_defect1_members_3345000.txt` (the 3 run members, 57 s), `maverick_defect1_3455000_4444400.txt`, `maverick_defect1_scan_other5.txt` (3505000, 4550400, 5040330, 5050030, 5504030), `maverick_defect1_3505000_run.txt`. The walk is deterministic, so member order matches the old `maverick_defect1_3345000.txt` (cells (2,3), (2,3), (1,2)). The referee's 3455000 member is the first listed in the 3455000 file: path [2,5,1,1,3], rels [[0,1,3,4,5,6],[2,1,3,4,8],[3,4,5,6,7]].
2. Closure: reworded (Claim and Proposal). Not "controls complete"; "adequate for the one-, two- and three-cord cells measured". Closing T6 is now a proposal to PARK it as search-only, refutation-only, as E-185 did.
3. Filter and guard: stated in the Claim (cords >= 1; len(path) == L, shortest recorded path only; guard = `coxeterGuard=True`, the default of `search.mutationSearchDepthFirst`, which `quipuRelations.reachedQuipuAlgebras` calls without override; no argument in my script).
4. 4444400 relabelled "none within depth 5" in the table.
5. Referee's 3455000 member recorded: 83 s walk, 274 s total, member cell (2,3), hit True, neg True (referee's run, not mine).
6. (5,3,4) cell: located and run. It is cells (cords, rels) = (3,4), one member each in 3505000 and 5050030 at L = 5. Run on 3505000 (first member, path [7,5,1,2,3], rels [[3,1,0],[5,3,1,2],[5,3,4],[8,6,5,3,1]]): hit True, neg True, 187 s total (114 s walk). 5050030's (3,4) member was sized, not run. So a 3-cord member is now exercised.
7. Remaining LNAs: sized, not run on members: 4550400 (6 defect<=1), 5040330 (3), 5050030 (9), 5504030 (3), 3505000 (9). Member back-searches there are [next round]; the control is 5 member runs over 3 of 9 LNAs, not 9.

author: maverick · round: 069 · kind: result
thread: T6 · bears on: H-017, E-185, E-065, E-071
scope: n = 9 only, depth 5 (L = 5), monomial tree quipu members with >= 1 cord and shortest recorded path of length exactly 5, Coxeter guard on (search default). Run: 3 members of 3345000 (all its defect-1 members), the referee's 1 of 3455000, 1 of 3505000 (the 3-cord cell). Six other LNAs only sized, not run. Nothing at n >= 10 or depth >= 7.

## Claim

Filter, stated plainly: members are quipus reached in the forward walk to depth 5 that are not the hereditary line, with cords >= 1 (a cordless defect-1 member is excluded by design) and with recorded path length == L = 5 (a source with several members at different lengths contributes only its shortest-path ones to this test; fine for a completeness control). Coxeter guard: `coxeterGuard=True`, the default of `search.mutationSearchDepthFirst`, inherited through `quipuRelations.reachedQuipuAlgebras`; my script has no guard argument.

The one control E-185 left open is done on the cells measured. Three members of LNA 3345000 with defect (relations - cords) = 1, cells (cords, rels) = (2,3), (2,3), (1,2) at path length 5, were rebuilt as path algebras and searched back: the source is returned at depth 5 in 3 of 3 and at depth 4 in 0 of 3. Added after review: the referee's 3455000 member (2,3) passes (hit, neg; 274 s) and a 3505000 member of cell (3,4), a 3-cord member, passes (hit, neg; 187 s): 5 of 5. So the search finds defect-1 sources exactly as it found defect-2 ones (E-185). It does NOT show the search would detect a defect <= 0 member (none is known to exist; the sources are forward-reached, so this is completeness on reachable members). Refuted by: a defect-1 member whose back-search misses its source at depth L.

Proposal for the next ledger (narrowed): PARK T6 as "H-017 open, search-only, refutation-only; positive controls adequate for the one-, two- and three-cord cells measured (L = 5 and 6 at defect 2, L = 5 at defect 1, 5 members in 3 of 9 LNAs); the other 6 LNAs' members unrun; no deeper search without a proof route". The only thing that could reopen it is the already-listed overnight depth-7 census, which needs a human OVERNIGHT approval, not a thread slot. Reason to park (judgement, not a result): the cheap measurables are done, the Euler signature is a Cartan-matrix invariant and cannot see (cords, relations) (E-162), and further depth-6 negatives cost ~30 min per candidate (E-074).

## Evidence

Sizing (`--plan` equivalent, no `--run`; forward walk to depth 5, members with L = 5, cords >= 1, defect <= 1):

| LNA | members | L=5, cords>=1, defect<=1 | cells (cords,rels):count | walk time |
|---|---|---|---|---|
| 3345000 | 75 | 3 | (2,3):2, (1,2):1 | 51 s |
| 3455000 | 62 | 3 | (2,3):3 | 92 s (referee ran first member: 83 s walk, 274 s total, hit True, neg True) |
| 4444400 | 52 | 0 (none within depth 5) | - | 54-60 s |
| 3505000 | 81 | 9 | (2,3):6, (1,2):2, (3,4):1 | 114-127 s |
| 4550400 | 52 | 6 | (2,3):4, (1,2):2 | 109 s |
| 5040330 | 75 | 3 | (2,3):2, (1,2):1 | 55 s |
| 5050030 | 81 | 9 | (2,3):6, (1,2):2, (3,4):1 | 124 s |
| 5504030 | 62 | 3 | (2,3):3 | 79 s |

Run on 3345000 (517 s total): hit True, neg True for all three (cells (2,3), (2,3), (1,2)). Each member costs about 130-170 s (depth 5 plus depth 4 search). The (1,2) cell is the (5,1,2) cell named in STEERING. The (5,3,4) cell is in 3505000 and 5050030 (one member each); run on 3505000: hit True, neg True (187 s). 5050030's member and the other LNAs' members are [next round].

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/069/maverick_defect1.py 9 5 5 1 3 3345000 --run   # 517 s
timeout 10m .venv/bin/python workshop/rounds/069/maverick_defect1.py 9 5 5 1 2 3455000          # plan only, 93 s
timeout 10m .venv/bin/python workshop/rounds/069/maverick_defect1.py 9 5 5 1 1 3505000 --run   # 187 s, (3,4) member
timeout 10m .venv/bin/python workshop/rounds/069/maverick_defect1.py 9 5 5 1 1 3505000 4550400 5040330 5050030 5504030   # sizing + identities, 495 s
```
Output: `workshop/rounds/069/maverick_defect1_3345000.txt` (run), `..._members_3345000.txt`, `..._3505000_run.txt`, `..._scan_other5.txt`, `..._3455000_4444400.txt` (identities).

## Prior record

E-185 (round 061): defect-2 controls only, "no defect-1 source was used". E-071/E-074/E-083/E-089: the search returns the source iff a member lies within depth. H-017 stays OPEN; this adds a measured defect-1 control and no new evidence for or against the claim. Not a rediscovery: grep of research/ for "defect-1" / "(5,1,2)" finds only the census cell in E-185's plan and the STEERING request.

## Code changed

New `workshop/rounds/069/maverick_defect1.py` (now also prints member identities; variant of round 061's maverick_cordcontrol.py taking LNA names and a defect bound). No library change; no tests needed.

## Next

Chair: put "park T6 (controls adequate for measured cells)" on the ledger; keep the depth-7 census as an OVERNIGHT proposal only. Referee: re-run the 3505000 (3,4) member (~190 s). [next round]: member back-searches on 5050030, 4550400, 5040330, 5504030.
