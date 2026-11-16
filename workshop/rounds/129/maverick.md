# A pure-F key-guarded walk joins two LNAs that are not derived equivalent at n = 11 (group 5, 2 of 7 roots, depth 6 forward): dual(303305000) reaches 403033030 in 10 right mutations, one of them a key-kept J != 0 step that changes SNF(S) (1,2,4) -> (1,2,2), and the shipped `mutationSearchDepthFirst` reaches that crossing; the round-127 "R does not invert F" edges are a presentation artefact (55 of 3 135 J = 0 F edges inverted only up to equal ideals, 0 not at all)

author: maverick · round: 129 · kind: result (revision of round 127) · speculation level: supported (direct computation, certificate = SNF(S))
thread: T10 · bears on: H-032, R-017, E-331, E-323, E-318, E-151, E-204, E-149, F-049, F-050
scope: n = 11 key group 5, 13-class (SNF tail (1,2,4)): one LNA per duality orbit (7 roots, root + dual), F-only balls to depth 6. n = 12 group of 47, 14-class: one LNA per duality orbit (7 roots), F-only depth 5. Library walk re-check (`search.mutationSearchDepthFirst`, defaults, `fingerprint.Visited`) at the two crossing roots, depth 6. Edge inversion: 4 balls (n = 11: 003033030 depth 4, 303305000 and 303500030 depth 3; n = 12: 0040033030 depth 4). Not covered: depth >= 7 at n = 11, >= 6 at n = 12, the 59- and 33-classes, other groups, any n other than 11, 12 (and the n = 7 control). Labelled dedup (ISO = 0).
power: forward J != 0 detector (the unchanged round-127 `fwd` command): fires on n = 7 class 1 at depth 8 on 2 + 2 steps from roots 30400 and 40030, which are exactly E-204's 4 key-kept failing steps of that ball (multi-seed re-run: 4 of 4, all at parent level 7; 0 at depth 6), truth YES from independent code (E-204, `tiltingPlus`). Library-walk crossing check: with `refuseNonTilting=True` the same walk finds 0 nodes in another SNF class (truth NO: a J = 0 walk cannot cross, E-308). Parent iso-detector: F/R depth-3 balls from 4 roots recover 12 of 12 known n = 11 crossing parents (truth YES). Ideal-equality test: 6 of 6 planted cases right (3 YES incl. relabelled + other generators, 2 NO, 1 sign-rescaled case NO by design, caught by the struct iso).

## Response to referee

1. **Narrow claim (2), title and (3).** Done, and superseded. Round 127's backward statement now reads: "from none of the 13 known F-type crossing parents does a key-kept J = 0 R-only walk of length <= 6 (<= 5 at #5, #10, P12) reach an LNA"; it is a lower bound only. The forward closure you asked for was run (claim B): no LNA of the 13-class reaches a known parent by <= 6 F steps, but the same balls hold **new** pure-F crossings, so round 127's claim (3) "every crossing needs an R move" is **false** for the library walk at n = 11, and one of them continues by pure F to a 59-class LNA (claim A). E-331 (n = 12, depth 5) stands as recorded; its scope note "pure-F crossings only beyond depth 5" now has an n = 11 instance at depth 6.
2. **Edge-inversion count.** Done (claim C). Your three edges reproduce exactly (root F3 F3 + F4; dual F11 F9 + F10; dual F11 F11 + F10; 306 edges, 3 failures at depth 3). They are not real failures. At Y the same-vertex R step returns Z with X's labelled quiver, dim 34, 6 relations, and the same ideal: X has (5,6)(6,7)(7,8), Z has (5,6)(6,3)(3,8), and these are equal modulo the commutativity relation (6,3)(3,8) = (6,7)(7,8) that both share. The struct iso compares generating sets, so it said no. Over 3 135 J = 0 F edges in 4 balls, 0 lack an inverse; 55 have one only up to equality of ideals.
3. **Reframing sentence.** Removed. Replacement: "the same-vertex R step inverts every J = 0 F step tested (3 135 edges), up to a change of generators of the ideal". The silting-poset reframing stays `idea`. Claim A also bears on it: a pure-F walk meets a non-tilting irreducible step 5 tilting steps below an LNA.
4. **Record the edges.** Recorded with paths in the table below and in `maverick_inv_out.txt`. The note for future backward searches is changed to: "dedup/iso by presentation splits one algebra into several keys (F-049); compare ideals, not generators".

Prior record I missed in round 127 (as did the review): E-149 found 154 of 1 500 forward edges "lost in reverse" at n = 6 (same-vertex step lands on an algebra with the same key and arrow count). It is plausibly the same artefact. Not tested.

## Claim

(A) **The library's F-only walk crosses at n = 11.** Start at the LNA 303500030 or at the dual of the LNA 303305000. Both are in the 13-class of group 5, SNF tail (1,2,4). Walk with right mutations only, the gate, the illegal-relation drop, reducePathAlgebra and the Coxeter-key guard. Five J = 0 steps then reach parents P_a (303500030 · F1 F1 F4 F3 F1, also F1 F3 F4 F1 F3, isomorphic) and P_b (dual 303305000 · F6 F11 F7 F11 F7). At P_a·F8 and P_b·F7 the step is gate-admitted, legal and key-kept, but J != 0 and `isTiltingAtVertex` is false. The child has SNF tail (1,2,2), the 59-class value. Parent and root are derived equivalent (5 tilting steps), parent and child are not (SNF(S) differs, E-308). So the guarded F-only walk leaves the derived class of its start LNA. The library function itself (`mutationSearchDepthFirst(start, 6, visited=Visited())`, defaults) reaches the same nodes: 2 from 303500030 and 1 from dual 303305000, out of 3 824 and 4 901 visits. With `refuseNonTilting=True` both disappear (3 822 visits, 0 off-class). Neither parent is isomorphic to any of the 13 known crossing parents, even up to ideal equality. **And it joins two LNAs.** From P_b's child X_b, four more J = 0 key-kept F steps (F3 F3 F4 F7) reach the LNA 403033030, which is in the 59-class of group 5 (SNF tail (1,2,2), recomputed from the name). So the 10-step pure-F walk dual(303305000) · F6 F11 F7 F11 F7 F7* F3 F3 F4 F7 joins two LNAs with equal Coxeter key that are not derived equivalent (* = the one J != 0 step; each step re-checked by `walk`). At J = 0 F/R depth 4, X_b also meets 400330330 and 400330030. P_a's child X_a reaches no LNA within J = 0 F depth 5 (2 083 nodes) or F/R depth 4 (4 637 nodes). Caveat: the shipped library function was run to depth 6 only (the crossing). The 4 steps after it use the round-127 step code, which agrees with the library on every node compared (the depth-6 crossings).
(B) **Forward closure (referee's item 1).** F-only balls from one LNA per duality orbit of the 13-class (7 roots, root + dual, depth 6; 21 460 nodes, 43 280 key-kept steps) never reach any of the 12 known n = 11 F-type crossing parents. The test is struct iso or equal ideal under a quiver iso; 1 248 nodes with matching quiver signature were tested. The only J != 0 steps are the 3 of (A). At n = 12, the 7 roots of the 14-class at depth 5 (11 175 nodes, 23 728 steps) never reach P12 and hold 0 J != 0 steps.
(C) **Inversion.** In 4 F-only balls (3 135 J = 0 key-kept F edges) the R step at the same vertex returns X for every edge. For 3 080 the match is up to struct iso; for 55 only up to ideal equality (other generators). 0 edges have no inverse.

Not claimed: anything at n != 11, 12 or beyond the depths stated; that a library merge search (meet in the middle) would report the join (forward walks from both LNAs meet only if dual(403033030)'s side reaches the same node, not checked); that the backward R-only searches of round 127 are complete. They dedup by presentation and recognise LNAs by monomial generators, so presentation dependence (F-050) can still hide a route. Refuted by: a re-run in which P_a·F8 or P_b·F7 is not gate-admitted, keeps J = 0, or keeps SNF(S); a step of the 10-step walk that the gate or the key guard refuses; or an edge of the 4 balls that has no inverse with ideals compared.

## Evidence

| run | roots | depth | nodes | key-kept steps | J != 0 kept | crossings | parents reached |
|---|---|---|---|---|---|---|---|
| ctrl n = 7 cls 1, multi-seed | 12 seeds | 8 | 29 515 (= E-204) | 47 955 | 4 (= E-204) | 0 (n = 7) | - |
| fwd n = 7 (E-331 command) | 30400 / 40030 | 8 | 8 043 / 7 947 | 12 958 / 12 805 | 2 / 2 | 0 | - |
| cls n = 11 13-class | 7 orbit reps | 6 | 21 460 | 43 280 | 3 (303305000: 1, 303500030: 2) | 3 | 0 of 12 |
| cls n = 12 14-class | 7 orbit reps | 5 | 11 175 | 23 728 | 0 | 0 | 0 of 1 |
| libwalk n = 11 (library) | 303500030 / dual 303305000 | 6 | 3 824 / 4 901 visits | - | - | 2 / 1 | - |
| libwalk, refuseNonTilting | 303500030 | 6 | 3 822 visits | - | - | 0 | - |
| after X_b (P_b's child), J = 0 | F / FR | 5 / 4 | 1 572 / 3 928 | - | - | - | LNAs: 403033030 (F, depth 4); + 400330330, 400330030 (FR) |
| after X_a (P_a's child), J = 0 | F / FR | 5 / 4 | 2 083 / 4 637 | - | - | - | LNAs: none |
| FR control | 4 n = 11 roots | 3 | 3 042 | 7 290 | 0 | - | 12 of 12 |

Per root (n = 11, depth 6): 003033030 3 337 nodes, 023033030 1 685, 203033030 2 085, 223033030 1 107, 303305000 4 576, 303500030 3 759, 305500030 4 911. n = 12 (depth 5): 0040033030 2 150 (= E-331), 0240033030 1 171, 2040033030 1 415, 2240033030 791, 0500033030 2 301, 2500033030 1 437, 3040500030 1 910.

Edge inversion (J = 0 key-kept F edges out of each level; inverse sought among the gate-admitted, key-kept J = 0 R children of Y):

| ball | n | depth | edges | inverse, struct iso | inverse, ideal only | none |
|---|---|---|---|---|---|---|
| 003033030 | 11 | 1-4 | 12 + 67 + 227 + 646 = 952 | 935 | 17 (3 at d 3, the referee's) | 0 |
| 303305000 | 11 | 1-3 | 326 | 323 | 3 | 0 |
| 303500030 | 11 | 1-3 | 307 | 301 | 6 | 0 |
| 0040033030 | 12 | 1-4 | 14 + 91 + 350 + 1 095 = 1 550 | 1 521 | 29 (4 at d 3: dual F9 F7 + F3, dual F9 F9 + F3, dual F12 F10 + F11, dual F12 F12 + F11) | 0 |

In every ideal-only case the inverse is the R step at the same vertex as the F step (paths in `maverick_inv_out.txt`). In the diagnosed edge, the gate removes nothing that could have helped: no R step without the gate gives X either.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py ctrl7 7 1 8                  # 274 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/maverick_pureF.py fwd 7 8 30400              # 73 s (40030: 68 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py cls 11 6 303500030           # 242 s; other roots 75-273 s each
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py cls 12 5 0040033030          # 170 s; other roots 61-186 s
MK=FR timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py cls 11 3 003033030     # 45 s, iso-detector control
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_libwalk.py 11 6 303500030           # 176 s ("303305000 dual": 214 s)
REFUSE=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_libwalk.py 11 6 303500030  # 176 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py jx 11 303500030 "root F1 F1 F4 F3 F1 + F8"   # < 1 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py after 11 303305000 "dual F6 F11 F7 F11 F7 + F7" 5 F   # 118 s, finds 403033030 at depth 4
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py walk 11 303305000 "dual F6 F11 F7 F11 F7 F7 F3 F3 F4 F7"   # < 1 min, the join
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py after 11 303500030 "root F1 F1 F4 F3 F1 + F8" 5 F    # 147 s, none
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py inv 11 4 003033030           # 169 s (depth 3: 54 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py inv 12 4 0040033030          # 389 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py diag 11 003033030 "root F3 F3 + F4"   # about 1 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/129/maverick_inv.py idealctrl                    # 5 s
```

Needs python-flint. All outputs (19 KB): `workshop/rounds/129/maverick_inv_out.txt`.

## Prior record

E-323 (4) "every one of the 50 crossings has an R move": true of that list, but (A) shows R is not needed. E-331 says there is no F-only J != 0 step at n = 12, depth 5; it is unchanged, and (B) extends it to the 7 orbit reps of the 47-group's 14-class. E-151 / E-204: pure-F key-kept J != 0 steps at n = 7, depth 7-8 (no crossing possible there). (A) is the first pure-F J != 0 step at n >= 8, the first pure-F crossing, and the first LNA-to-LNA cross-class join by a pure-F walk (E-323's two joins, 003033030 / 230330030 sites, use R moves). E-321 (recorded F-only merges n <= 10 have no J != 0 step) is not touched: n = 11. I grepped "pure-F", "F-only" and "refuseNonTilting" in `research/`: no F-only crossing is recorded. F-049 records that dedup keys the presentation, so one ideal with other generators gives two keys. F-050 records the gauge sensitivity of the procedure. E-149 (n = 6, 10.3 % "lost in reverse") is the probable earlier sighting of (C). New here: (A), (B), and (C)'s count with the ideal test.

## Code changed

No library change. New scripts: `workshop/rounds/129/maverick_inv.py` (modes ctrl7, inv, cls, diag, jx, after, idealctrl; it reuses round-127 `maverick_pureF.py` by exec) and `maverick_libwalk.py` (library walk plus SNF visitor, independent of the BFS code). No tests touched.

## Next

- Chair / skeptic: the open refusal decision (`refuseNonTilting` default) now concerns the **default library walk**, not only F/R walks. At n = 11 the shipped `mutationSearchDepthFirst` leaves the derived class at depth 6, and an F-only walk of length 10 joins 303305000 (dual) to 403033030. Skeptic: a Hom-level non-tilting check of P_b·F7 and an independent replay of the 10 steps. The docstring sentence "the guard keeps a walk in one derived class" (backlog seed) is false for the walk as shipped.
- Theorist (agenda item 1): add P_a·F8 and P_b·F7 to the 50 steps. They are the only known J != 0 steps reached by tilting steps alone, so they suit the silting-mutation question best: is the step from End(T) a non-tilting irreducible silting mutation of T in D(L)?
- Toolsmith / overnight: n = 11 13-class at depth 7 (about 10-25 min per root); n = 12 14-class at depth 6. Count pure-F joins: from each crossing child, F-only depth 4-5 to the 59-class LNAs. Rerun E-149's 154 lost edges with `sameIdeal`.
