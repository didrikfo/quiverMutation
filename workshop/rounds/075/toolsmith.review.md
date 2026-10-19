# Review of workshop/rounds/075/toolsmith.md

referee: skeptic · round: 075
verdict: minor revision

## Reproduction
Re-ran `toolsmith_isojoin.py R 4 5 400000 500` for R = 0..3 (from repo root, four in parallel, ~90 s each, all uncapped). Output identical to the committed `toolsmith_isojoin_r{0..3}_out.txt` apart from timings: child totals 2, 2, 2, 3; parent totals 6, 6, 5, 6; level sizes and bucket counts (970 / 682) match. The search-1 table (`toolsmith_join.py`, 5-10 min each) was not re-run; its only load-bearing entry is the rep 0 child MISS, which is a negative and is not needed for the claim. Note the script must be run from the repo root (it opens `workshop/rounds/075/toolsmith_join.py` by relative path); the submission says "repo root" only in the docstring.

## True?
No error found. Read `iso`: a HIT requires a structure-graph isomorphism whose relabelled relations equal the target's gauge-canonical relations, so a false HIT would need a bug in `_relationsUnder`, not in the search. The `limit=3000` isomorphism cap can only cause misses, not false hits. Reversing the seed-to-node F/R path is legitimate (F and R are inverse tilting moves). Remaining unchecked points:
1. J = 0 is `tiltingPlus` only, as the author says; no Hom(T,T[-1]) check, so these are not Hom-verified chains the way E-212's are.
2. The path edge lists were not recorded, so the 073 verifier cannot be run on them yet.
3. "Seed" = `classes[base]`, i.e. the round-049 LNA list and duals for the key class; I did not independently confirm every seed is an LNA at n = 7, only that this is the round-049 construction used by E-157/E-163.
4. "The 4 reps are the E-213 unjoined classes": the 4 `NEW class rep` blocks in `rounds/074/skeptic_failkids_an_out.txt` (cls 1 v 7, cls 2 v 3/6/7) match E-213's "4 of 12 classes" in count and class split; I did not recompute that these 4 are the ones containing no E-151 window child.

## New?
Grepped `research/` for iso meet, isomorphism at meets, E-213, unjoined, the 4 classes. E-212 already uses isomorphism (not labelled key) at meets with Hom-verified chains for 8 of 9 class-2 window steps (class 1 not run); E-213 lists these 4 classes as having no recorded join; E-164 records that `canonicalKey` is labelled. The submission's new content is only: the 4 unjoined classes now have J = 0 (tiltingPlus-level) iso joins to seeds, and the labelled-key miss at rep 0 child is a relabelling artefact. Modest, but genuinely not recorded.

## Evidenced?
Mostly. Ranges, depths, caps and bucket counts are stated and reproduce. Gaps: (a) the search-1 table is not committed (stdout only), and the submission leans on it for the "labelled search missed" sentence in the title; (b) totals are first-meet levels, stated as upper bounds, fine; (c) no chain is recorded, so "at most 3 + 4 edges" cannot be audited per edge.

## Scope
The title says "reach an LNA by J = 0 F/R `tiltingPlus` paths of at most 3 + 4 edges". The data say: child of each rep is isomorphic to a node at distance 2, 2, 2, 3 from a seed; parents at 6, 6, 5, 6 totals (alg side 2, 2, 1, 2 + target side 4). Since sides are 0+2, 0+3 and <=2+4, "at most 3 + 4" is right. Narrowed wording: "...under the tiltingPlus (not Hom-verified) J = 0 premise, and joined to possibly different seeds for parent and child, so this says nothing about parent ~ child". The author already says the last part; the title omits it. "The labelled search missed the class-1 child at 6 + 5" should say it was a single run with the round-049 move set.

## Required for acceptance
1. Commit the search-1 stdout (or drop the 6 + 5 sentence from the title) so the labelled MISS is auditable.
2. Add "tiltingPlus-level, not Hom-verified" to the title/claim line, since E-212 sets the Hom-verified standard for the same kind of statement.
3. State that scripts must be run from the repo root in the Reproduction block.
4. [next round] Record the F/R move lists and run them through the 073 verifier, as the author proposes under Next.
