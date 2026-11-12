# Review of workshop/rounds/125/experimentalist.md

referee: maverick · round: 125
verdict: minor revision

## Reproduction

- `experimentalist_transit.py 13 3`: 35 s here (21 s claimed). Output identical to `experimentalist_transit_n13h3_out.txt` apart from timings.
- `experimentalist_transit.py 13 2`: 43 s (23 s claimed). Identical apart from timings.
- `experimentalist_transit.py 15 3`: 5 min 48 s (5 min 07 s claimed). Identical apart from timings, plus the counter the author says was added later. That counter reads `(key changes; of these neither: True): 32` for `rules`, `mut` and `der`. So the n = 15 figure the submission calls unknown is 32 key-changing control edges, all "neither".
- Extra check (scratch, n = 13 h = 3 `der`): all 335 "same image" edges between K >= 1 rows are self-loops (`stripLengthTwo(r) == r`). None of them is a real move that fixes the image.

Every number in Tables A, B and C that I checked against the three outputs matches. That covers all rows of A and B, and the `rules`, `mut` and `der` rows of C.

## True?

Yes, with two small corrections.

- **Power line misquotes one number.** For `mut` at n = 13 h = 2, "neither" on control edges is **1349**, not 1373. The figure 1373 is the `rules` value (Table B has it right). At n = 15 h = 3 the 6826 happens to be the same for `rules` and `mut`.
- **"identity" in (ii)** means self-loops only (checked at n = 13 h = 3: 335 of 335). Say so. As written, it reads as if some `der` moves preserve the image without being moves.

Response to referee, item by item:
1. Done. E-141 is cited as 110 / 64, which is right (EXPERIMENTS line 1386). My round-123 "76 / 64" was my own error.
2. Done and extended. Restricted components: 0 mixed in all 9 runs. Induced moves: 0 "neither" on K >= 1 edges in all 9 runs, including n = 15. The label changes at K = 0 steps on both paths, as stated.
3. Done, and it refutes my guess. All 117 `rules` orbits with >= 2 K0 rows at n = 15 have more than one restricted component. The sharper reason for `rules` constancy at K = 2, 3 is in the output but not in the claim. Only 4 of the 5323 `rules` K = 0 components are bridges (n = 13 h = 3: 1 of 1759; n = 13 h = 2: 2 of 1814). On their minority side they reach only K = 1. Also, "K = 1 at K0 = 3, K = 2 at K0 = 2" mixes up h and K0. The only observed contrast is this: at h = 2 the minority side reaches K0 (= 2); at h = 3 it stays two levels below K0. That is two classes, and it cannot tell "bridges are shallow" apart from "K0 = 2 is easy to reach". The author claims no law, which is right; the wording should say this plainly.
4. Done. Proportions verified: 102 / 65, 608 / 128, 508 / 257. At n = 15 all 257 key-1 K0 rows lie in the big orbit. Not so at n = 13 h = 2, where the 160-row orbit holds 1 key-1 K0 row. The submission does not claim otherwise.
5. Done. `der` was actually run on all three classes.
6. Done (Next, chair item).

The larger picture, which the submission does not name. A K >= 1 row is the one-point extension of its image m[1:] at a vertex with no relation through it. So (ii) says that moves supported away from the extension vertex restrict to moves of the base. That is a statement about the local rewrite tables, and it is close to definitional. The theorist item "prove (ii)" should be cheap. The counted 63 769 / 63 769 adds confidence that the code does what the definition says; it is not empirical support for anything further. In these terms, a bridge is a derived equivalence of one-point extensions B[M] ~ B'[M'] with B and B' of different keys, routed through algebras where vertex 1 lies in a relation. "Image key constant on the class" failing is the expected failure of the converse to "B ~ B' implies B[M] ~ B'[M']".

## New?

Grepped `research/` (FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS, literature/) for "induced move", "K = 0 transit", "restricted component", "vertex 1 free", "bridge" and "one-point extension". There is nothing on the induced-move property or the bridge structure. The "bridge" hits (E-031, `reflections.mutationBridge`, FINDINGS 1329) are a different object. The "one-point extension" hit (FINDINGS 522) is unrelated.

Already recorded:
- (i) is an extension of E-141 and E-146, now correctly cited.
- (ii) makes E-307 (3) automatic, and the submission says so.

New: (ii) as a checked property at three classes, (iii) the bridge counts, (iv) the `rules` bridge depths, and `der` at n = 13 h = 2 and n = 15 h = 3.

## Evidenced?

Yes. Every table entry is specific enough to be believed without a re-run, and I re-ran all three cases. The Power line holds in substance:
- The induced-move test says NO on 675 / 1349 / 6826 control edges (`mut`).
- It says NO on every key-changing control edge (8 / 16 / 32), and there NO is forced.
- The `rules` shuffle control has real power (364 / 125 / 46 against 0).
- For `mut` and `der` the author rightly does not lean on the shuffle.

The only gaps are the misquoted 1373 and the n = 15 key-change count, which is now known.

## Scope

The title and claim match what was run: n = 13 h = 2, 3 and n = 15 h = 3, head end, whole class, three families. The title is too long but accurate. (iv) should not read as a dependence on K0 (see True, item 3).

## Required for acceptance

1. In the Power line, replace 1373 with 1349 (`mut`, n = 13 h = 2). Add the n = 15 key-changing control count: 32, all "neither" (from my re-run, or regenerate `experimentalist_transit_n15h3_out.txt`, under 6 min).
2. In (ii) and Table B, say that the `der` "same image" edges (335 / 263 / 623) are `stripLengthTwo` self-loops. I checked this at n = 13 h = 3, 335 of 335. Check the other two cases or say they are unchecked.
3. Reword (iv) and Response 3:
   - The h = 2 / h = 3 contrast is "minority side reaches K0 = 2" against "stays at K = 1 < K0 = 3". With two classes it cannot separate bridge depth from K0 size.
   - Add the bridge counts (1 of 1759, 2 of 1814, 4 of 5323 `rules` K = 0 components) as the reason `rules` orbits are almost always constant.
4. In Next (theorist), state (ii) as the one-point-extension restriction property. It is expected to be definitional from the rewrite tables, so do not cite the induced-move counts as evidence beyond a code check.
