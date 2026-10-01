# No monomial cord member exists in any LNA orbit searched (n = 4, 5 deep; n = 8 shallow), and every one of 2376 n = 8 cord members has its cord inside a commutativity relation; no positive MONO control at n = 8 could be built

author: maverick · round: 018 · kind: negative
thread: T6 · bears on: H-017, E-082, E-087, E-089

## Claim

Speculation level: tested on small cases. (1) A positive control for `MONO=1` at n = 8 was not found and, on this evidence, probably does not exist near the LNAs: monomial cord members (arrows >= n, no parallel arrows, all relations monomial) are absent from the orbit walks of all 5 LNAs at n = 4 (depth 8) and all 14 at n = 5 (depth 7), while sum-relation cord members appear at path length 1 (n = 4: from path length 1; n = 5: from path length 1). (2) The heuristic "a cord needs a sum relation" holds on all 2376 labelled cord members at n = 8, depth <= 5, from the six LNAs (indices 4, 9-13) that have any: every member has a sum relation, none is monomial-only, and in all 2376 every arrow lying on an undirected cycle lies on a path of some sum relation (cord = commutativity cycle). It does NOT claim a proof, nor that no monomial cord is derived equivalent to an LNA: depth is bounded, and LNAs 16..428 at n = 8 are not walked.
(3) Code-level control only: the `MONO=1` filter does return nonzero, including a member other than the seed, when the walk starts from a hand-built monomial cord algebra (not an LNA).

## Evidence

Producing relations (`maverick_producer.py 8 5 4 9 10 11 12 13`, 5 min total; same visitor as `toolsmith_cords.py`, both directions, members at path length <= 5):

| LNA (index) | members | 1 cycle | 2 cycles | monomial-only | cycle covered by sum relation |
|---|---|---|---|---|---|
| 000030 (4) | 448 | 448 | 0 | 0 | 448 |
| 000230 (9) | 114 | 114 | 0 | 0 | 114 |
| 000300 (10) | 444 | 432 | 12 | 0 | 444 |
| 000302 (11) | 226 | 225 | 1 | 0 | 226 |
| 000330 (12) | 464 | 464 | 0 | 0 | 464 |
| 000400 (13) | 680 | 615 | 65 | 0 | 680 |

Shapes of the sum relations (path lengths in arrows, members containing one, total): 2+2 1642, 3+2 358, 2+3 190, 4+2 114, 3+3 38, 5+2 25, 4+3 17, 2+4 10, 5+3 4, 6+2 1. So the producing relation is a commutativity of two paths, usually a square (2+2, 69%), never longer than 6+2; members with two independent cycles (78) are those from LNAs 10, 11, 13 (the "3"/"4" in position 4).

Monomial negatives (`MONO=1 toolsmith_cords.py N L 0 1 0 100 --plan`, MINREL 0): n = 4 L = 8: 0 members for all 5 LNAs (2 s each); n = 5 L = 7: 0 for all 14 LNAs (<= 6 s each). Without MONO the same quivers have cords from path length 1 (n = 4: LNA `30`, 12 members over L <= 6, all with 1 relation). Together with E-087 (n = 6, L = 5; n = 7, L = 5; all 0) and E-089 (n = 8, L = 6, six LNAs: 0), no monomial cord member has been seen at n = 4..8.

Weakness of the polynomial route (`maverick_monocord.py`): all unicyclic acyclic quivers with n arrows (n = 4: 9, n = 5: 54, up to isomorphism) with every monomial ideal (42 and 736 algebras): 15 and 190 of them share their Coxeter polynomial with some LNA. Only 2 distinct polynomials occur among the LNAs at n = 4, 5 (A_n-type-like classes), so the polynomial cannot support the heuristic; e.g. the triangle 0->1->2, 0->2 with zero relation 0->1->2 plus a leaf has the polynomial of `00`-class LNAs, and is nevertheless never reached. Not derived equivalence either way.

Code-level control (`maverick_filtercheck.py`, seed = triangle 1->2->3, 1->3, relation 1->2->3 = 0, leaf 3->4): `MONO=1` L = 5 returns 4 monomial cord members: 2 at path length 0 (seed and its dual) and 2 at length 1 (one with three monomial relations). So the filter and the dual both work. The walk prints "ILLEGAL RELATION" for several mutations of this seed (an oriented cycle created at a non-tilting step; harmless to the visitor), a sign the seed is not LNA-like.

Mechanism guess (idea, untested): a cord appears when a mutation composes an arrow with a path parallel to an existing path; the relation step 7 then keeps the two paths equal (a sum), and a monomial one would need the second path to vanish, which the Cartan determinant (unitriangular) does not forbid but which never occurred.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/018/maverick_producer.py 8 5 4 9 10 11 12 13   # 5 min; output maverick_producer_n8_L5.txt
MONO=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 4 8 0 1 0 100 --plan   # 8 s; maverick_mono_n4_L8.txt
MONO=1 timeout 10m .venv/bin/python workshop/rounds/015/toolsmith_cords.py 5 7 0 1 0 100 --plan   # 40 s; maverick_mono_n5_L7.txt
timeout 10m .venv/bin/python workshop/rounds/018/maverick_monocord.py 5     # 1 min; maverick_monocord_n{4,5}.txt
MONO=1 timeout 5m .venv/bin/python workshop/rounds/018/maverick_filtercheck.py 5   # seconds; maverick_filtercheck_L5.txt
```

## Prior record

E-087 states no monomial cord member at n = 6, 7 and that every walked cord member has a sum relation; E-089 states the n = 8 L = 6 `MONO` negative and calls "a cord needs a sum relation" a heuristic with a call for a positive control and the producing relation log. Not in RETRACTIONS. New here: the per-member log (2376 members, 100% sum, 100% cycle covered), the deeper n = 4, 5 negatives, the code-level positive control, and the Coxeter-polynomial non-test. Not new: the n = 8 sum-cord existence.

## Code changed

None in `quivermutation/`. New scripts only: `maverick_producer.py`, `maverick_monocord.py`, `maverick_filtercheck.py` (all in `workshop/rounds/018/`). No tests touched.

## Next

- toolsmith: `MONO=1` L = 5 `--plan` for LNAs 16..428 at n = 8 only matters if some LNA has cords at L = 5; the producer log says cords come from LNAs with a 3 or 4 early; first run the cheap non-MONO `--plan` over all 429 (one count per LNA) in shards and run MONO only where nonzero.
- theorist: prove or kill "cord = commutativity cycle": an invariant of the mutation step 7 (a new arrow parallel to a path is created together with its sum relation) would turn the heuristic into a lemma; try the n = 4 orbit (all cord members have one relation) by hand.
- skeptic: is a derived-equivalent monomial cord (e.g. discrete derived Lambda(1,3,m)) provably outside every LNA class? The polynomial cannot tell.
- consequence for H-017: the E-076 candidates are monomial quipus with cords (2-3 cords); the control members have sum cords, so the control class and the candidate class differ in kind, not only in size.
