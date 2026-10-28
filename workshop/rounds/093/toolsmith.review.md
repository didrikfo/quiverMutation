# Review of workshop/rounds/093/toolsmith.md

referee: skeptic · round: 093
verdict: minor revision

## Reproduction

Re-ran the builds (A depth 8: first run exit 2 at budget, depth 7, 15 177 nodes; rerun 169 s to 20 210; B depth 6: 482 s, 15 401), then `lna`, `extra`, `meet --extra`, and the B-B control, all as given. Same output in every line I compared:
- ball sizes 20 210 / 15 401;
- A LNA census 1,1,3,4,2,19,28,27,75 with 9 classes, B 1,1,4,4,4,22,32 with 9 classes;
- 3 components touched on each side, 0 common;
- 2400230 component size 1;
- 8427 and 8063 children isomorphic but not labelled-equal;
- WL-bucket pairs 0 and iso tests 0 for A(8) against B(6+8427 child);
- B-B control: 95 449 pairs, 1 iso test, total 0.

I did not re-run the `meet` without `--extra`.

## True?

Nothing wrong found. Two points:
1. The union MISS is nearly forced by E-238 and adds little on top of the census. The jump can only be taken from an LNA node. A's LNA nodes lie in three components, none of which is 3304000's, and E-238 says the certificate connects them only through naming. A MISS for "F/R + one jump" is therefore implied by "common components = 0", and the 160-node scan does not test anything beyond that. The write-up says so ("as strong as the movable edge set"), but the title leads with the F/R+jump meet as if it were an independent test.
2. "2400230 isolated in the lnaMoves graph" depends on `movesByRule` being the complete edge set, with edges taken as undirected. This matches E-238's graph, so it is consistent. It is an observation about computed moves, not about equivalence.

## New?

Mostly new. Greps of `research/` found:
- E-226, E-238, E-244 and E-254 cover the F/R meets, the certificate, the 8063 and 11123 joins, and the 8427 proposal;
- no record of the union run, of 2400230's lnaMoves component (size 1), or of the 8427/8063 isomorphism (only the 087 `toolsmith_d6distinct_out.txt` grouping, which the author cites);
- E-254 line 13 says "8427 not run", so closing it as vacuous by isomorphism is new to the record.

## Evidenced?

Yes. Per-depth counts, component sizes, the control and the exact sizes are stated, and they match on re-run. The power limitation (no cross-component positive control) is stated honestly. The only soft spot is that the jump over-approximation is a single layer; this is declared in the scope.

## Scope

The title and claim match what was checked. One narrowing: "no meet on the union" should read "no meet with one jump layer, jumps from LNA nodes only". The title already says "one jump layer only". The 8427 point holds only at the level of "isomorphic child of 8063, whose join at total 6 is E-244".

## Required for acceptance

1. Add one sentence to Claim or Evidence saying the union MISS follows from the component census (0 common components) together with E-238. The F/R meet at 8 + 6 is the only new search content, and it repeats E-226/E-254 with the larger radius.
2. Put in the build instructions that the A8 build needs the rerun after exit 2. It is already in the comment, but a first run leaves a depth-7 pickle, and `lna` then gives depths 7 6 and different counts (59 at depth 8, not 75).
