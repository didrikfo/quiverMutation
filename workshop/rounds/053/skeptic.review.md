# Review of workshop/rounds/053/skeptic.md

referee: theorist · round: 053
verdict: minor revision

## Reproduction

Used the author's pickles in /tmp/sk53 (not rebuilt; the 500 s / 381 s collection was not re-run). Re-ran all three replays and the key-13 check:
- c1 child 12: `RESULT edges 9 all-accept True meet True`, 1.7 s. Per-edge lines show hm1 0, hp1 0, cartan True, J False, tp True.
- c2 child 6: `edges 12 all-accept True meet True`.
- c1 child 5: `edges 12 all-accept True meet True`.
- `skeptic_key13.py`: keys equal, non-None, labelled quiver equal, parent (8,7),(8,7).
The output matches the claim. The move lists match E-160 and `rounds/050/toolsmith_depth7_logs.txt` lines 43, 44 and 50.

## True?

No error found in the arithmetic. Gaps, none of which the author hides:
- The claim is "Hom(T,T[±1]) = 0 and Cartan(End T) = Cartan(child)". That is the Cartan matrix of End(T), not End(T) ≅ child. The author says so (item 2). The title, "tilting-tested route", should not be read as "derived-equivalence route".
- Vanishing is tested only at m = ±1. For a tilting complex of two-term shape you need m ≠ 0 in general. For two-term T, Hom(T,T[m]) is nonzero only for m in {-1,0,1}, so ±1 is enough. The author writes "not argued"; it is argued by the length of the complex. Add that sentence.
- Generation is unchecked. A self-orthogonal two-term T with n summands does generate K^b(proj A) when it is a Bongartz/Okuyama-Rickard replacement, but this is assumed, not tested. The author lists it.
- "Every one of the 25 children is joined" holds under the J = 0 premise and rests on the paths being edge-tested. Child 13 and child 15 rely on key equality only (see item 4).

## New?

Not new mathematics; the author says so. E-160 printed the paths and noted they were not replayed. E-161 did the E-157 paths and E-163 did c1 14/15. This closes E-160's stated gap. Grepped EXPERIMENTS.md for E-160, replay3 and "33 edges": no prior replay of these three paths. No RETRACTIONS conflict found.

## Evidenced?

Yes for what was checked: the edge counts, the per-class tables and the reproduction commands are specific. The per-edge output files are not committed (about 3 KB each), so a reader cannot check them without re-running. This is a small defect.

One weakness: the author concedes the test is not independent in information, since J = False and tiltingPlus = True held at every edge and E-161 saw gate = test on 1081 steps. The gain is code independence only. The abstract's phrase "independent tilting test" is fair on code but should say "same verdict as the gate, independently coded".

## Scope

The title says "all 25 ... now have a tilting-tested route". Narrowed wording: "all 25 have a path of edges each passing the Hom(T,T[±1]) = 0 and Cartan-matrix test; 24 by their own path, c1 child 13 by key equality with child 5". That count is not what the title says: the title's "except" clause puts it at 24 plus 1, and the body supports exactly that. It also needs "c1 child 15 by key equality with child 14 (E-163)". Check: E-163 says child 15 joins by key equality. So two children (13 and 15) are covered by key equality, not one. The title mentions only 13 but its claim is about the three E-160 paths, so this is a title-versus-E-163 inconsistency, not an error in the replays.

## Required for acceptance

1. Narrow the title and Claim: say "Cartan matrix of End(T) equals the child's", not that the route is a derived equivalence. Add "under the J = 0 premise".
2. State that the m = ±1 test suffices because T is a two-term complex, or remove the "not argued" remark.
3. Commit the three per-edge output files, or paste their `RESULT` lines and per-edge tags.
4. State the count of children covered by key equality only (13 and 15, per E-163), so that "25" is read correctly.
5. [next round] End(T) ≅ child with relations, as the author's own Next list says. This is the open gap; it does not block this result.
