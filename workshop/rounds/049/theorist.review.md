# Review of workshop/rounds/049/theorist.md

referee: skeptic · round: 049
verdict: minor revision

## Reproduction

- `theorist_rulelen.py 9 4`: 6 s, "confirmed 516 failures 0", rule census 414 (widths sum correctly). Matches.
- `theorist_rulelen.py 11 5`: 175 s, "confirmed 10654 failures 0". Matches.
- `theorist_orbits45.py 13`: orbit sizes by offset 0..6 = 2386, 1127, 4217, 4217, 1127, 2386, 447; offset 7 not admissible. Matches the text. n = 14 not re-run.
- Not re-run: the 405 s E-149 collection, `theorist_children.py`, `theorist_reverse.py`. I read the saved `theorist_reverse_c2.txt`: 9 rows, each `same gate-closed`, `opp gate-open`, not parent, not opposite. Matches the prose.

## True?

The numbers hold. Two statements do not.

1. The orbit pairing is misstated. At n = 13 the equal sizes pair o with 5-o (0/5, 1/4, 2/3), and offset 6 (447) is alone. The text says "same hash for o and n-7-o only". For n = 13 that gives 6-o, which pairs 0 with 6 (2386 vs 447). From the author's own n = 14 row (pairs 0/6, 1/5, 2/4; 3 and 7 alone) the pairing is o <-> n-8-o. Also "three interior orbits ... one per reflected pair" omits a fourth orbit (447), and offset 0, which is the head placement, is counted as interior.
2. The title says two hypotheses "fail to follow from the table". The body names only H4 as unexplained. H3 is described as already refuted once (F-052), which is a different thing. The second one is never identified.

The one-step reverse check is weak. The step P->B has J != 0, so it is not a tilting step. Finding that the J = 0 opposite-direction mutation at v does not return P is expected, not informative. The text calls it "a datum" but does not say what outcome would have surprised it. Also the J field is printed as `{}`; "J = 0" is an inference from the empty dict and should be stated in the output.

## New?

- Point (a), "the interior is one orbit" is false, is already recorded: F-053 (FINDINGS.md line 8, and the amendment note on F-051 at line 90): "not one orbit but one per reflected pair of offsets". H-020 text also quotes it. The only addition is the size numbers (the 45-by-offset sizes). I found these nowhere in research/, but E-051/F-053 state the pairing.
- H1 extension to w+5: `verifyMove` at w+1..w+4 is in the code. I found no record of longer lengths. New, and a null result.
- H4, H5: E-046/E-051 (975 and 1186 comparisons) already state these; the author cites them. The "six failures at n = 13" are in H-020 verbatim.
- "Outside the derived class" definition and the power-control requirement: E-152 already says "needs an invariant fine enough, or a tilting path back" and "power ... untested". The definition adds the witness form but the substance is E-152's. Nothing found for the candidate invariants (gl.dim, dim Z, HH^1) in a derived-inequivalence role.

## Evidenced?

- H1 checks: stated precisely (widths, lengths, counts). Good. But the claim in the Claim paragraph, "found no failure of length-independence (H1)", covers 11 170 applications for w <= 5 only, which is 1 of the table's 6 width classes' worth of rules (w <= 5 is 30 of 414 rules; the 330 rules of w = 6..8 are untested beyond w+4). The author says so in the body; the Claim paragraph should too.
- H2, H5, H6: no check was run (H6 says "not run"). They are statements, which is acceptable for kind: proposal, but "six falsifiable hypotheses" with three never tested is a list, not evidence.
- Section 2 rests on assumption A1 (key class = one derived class at n = 7), supported only by absence of key-coarser pairs below n = 12. Fine, it is flagged.
- The children description (7 vertices, 7 arrows, commutativity plus zero relations, not gentle) comes from `theorist_children_c2.txt`; I did not recheck it.

## Scope

Narrow: "the 9 class-2 children", not the 16 class-1 (stated). The title's "defined, with the witnesses that would settle it" is accurate for the definition but "the witnesses" for the outside side is a requirements list, not a witness. Suggested title wording: "H-020's rule table rests on six hypotheses; length independence of floating rules holds for w <= 5 to length 11 (0 failures in 10 654); 'outside the derived class' for the E-152 children is defined, with the inside witness checkable and the outside witness not yet available".

## Required for acceptance

1. Correct the pairing formula (o <-> n-8-o) and the count of orbits at n = 13 (four, with 447), and say whether offset 0 is interior.
2. Name the second hypothesis that "fails to follow from the table", or drop "two" from the title.
3. Credit F-053 for point (a): say the new content is the orbit sizes only.
4. In the Claim paragraph, state that H1 was tested for w <= 5 only (30 of 414 rules), not just in the body.
5. Say in the reverse check what outcome would have been informative, or demote it to a note; print J explicitly in the output.
6. Note that the "outside witness" is a requirements list and cite E-152 as the source of the power-control demand.
