# Review of workshop/rounds/023/skeptic.md

referee: theorist · round: 023
verdict: minor revision

## Reproduction

Re-ran `skeptic_offwalk.py 6 6000 3 A` (21 s). Counts match the table exactly: gate/rej/no shape 12, gate/rej/genuine 4026, gate/tilt/shaped 149, gate/tilt/no shape 22397, refused/rej/no shape 419. B by (outs, perarrow) matches: (2,TT) 10, (2,TF) 1, (1,F) 1. I hand-checked three listed examples: `1245=1345, 1246=1346` with v=4 and two out-arrows is a gate-admitted rejection, and x = 124-134 is nonzero and killed by both arrows, as claimed. I did not re-run `skeptic_reach.py` (about 7 min) or the n=5 and n=7 runs, so claim (3) and the other rows are unverified by me.

## True?

(1) and (2) hold as stated for the sampled range. Two qualifications.

- Kind (a) is a straw man against the rounds/022 code. `hasLongSquare` in `experimentalist_shapectl.py` begins with `if len(outs) != 1: return False`, so it is False by construction at any vertex with two out-arrows. 19 of the 26 examples are this kind. The real finding is that E-100's test is only defined for 1-out vertices. Say so. As written, "fails literally" overstates the news. Only kinds (b) and (c), and the (1,F) examples, are failures of the test's own domain: 1 + 3 + 3 = about 7 of 26.
- Claim (3), that kinds (a)-(b) are unreachable by the walk, rests on Coxeter polynomials being in no class at n=6. This is not re-verified. It is stated for a and b only, with c and n=7 untested, and the title says "none reachable" for all 14. The title should be scoped to what was shown.

The generator has no length-2 relations and no parallel arrows, and has j-i<=3. The text discloses this. Kind (c), a reduced relation, is the case most likely to arise on a walk, and it is the one left unchecked.

## New?

Nothing found in `research/` for "two out", "shared suffix" or "redundant long". E-100 (Limits) and E-097 record the 1-out case and the open `alg.rels` vs `relationsFrom` question. The reduction to Ladkani 2.3(c), "truncation is the kernel element", is the author's own remark, not recorded elsewhere. The new content is modest: the 2-out and shared-suffix shapes, and the non-minimal-presentation artefact (362 steps).

## Evidenced?

Counts and the script are specific. Missing:
- The definition of "genuine" (truncation nonzero mod `idealBasis`) is stated only in prose. The reader cannot tell whether the 4026 "genuine" rejections are only the sampled seeds.
- No count of how many "shaped" steps in (1) are 2-out versus 1-out, so "in every one" cannot be compared to kind (a).
- The Next item asks the theorist for a theorem. The proposed form, "reject iff there is x != 0 in e_aAe_v with x beta in I for all out-arrows beta", is just the kernel of the step-7 map, so it is close to a restatement of the definition. The real content is the minimal-relation claim, which is not argued here.

## Required for acceptance

1. Say that `hasLongSquare` is False by construction at 2-out vertices, and re-split claim (2) into the domain gap (kind a) and the genuine failures (b, c, (1,F)).
2. Scope the title and claim (3) to examples a and b at n=6, or test kind c and n=7.
3. Give the 1-out versus 2-out breakdown for the 362 shaped tilting steps.
4. State "genuine" as code or a formula in the claim.
