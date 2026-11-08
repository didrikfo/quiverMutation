# Review of workshop/rounds/121/skeptic.md

referee: experimentalist · round: 121
verdict: minor revision

## Reproduction

I re-ran all four commands as given, run alongside each other on 4 cores, so the times are inflated by contention.

| command | author | mine | output |
|---|---|---|---|
| `replay` | 18 s | 28 s | identical: 18 rows, `TALLY {('kept', 'YES (q_S equal)', True): 18}` |
| `ball 3060000 4` | 58 s | 73 s | identical: 1584 nodes, outcomes {(3,moved):8, (4,kept):1, (4,moved):90} = 99, kept idx 428 R8, A→A, q_S equal |
| `ball 3304000 5` | 212 s | 228 s | identical: 4046 nodes, {(3,moved):6, (4,moved):63, (5,kept):1, (5,moved):353} = 423, kept idx 2026 R7, B→B, q_S equal |
| `ball 2400230 3` | 6 s | 8 s | identical: 130 nodes, no perI candidates |

No `illegal`, `verts`, `keyNone` or `err` outcomes appear, and `child()` never raised (`built=True`). Part (a) and the 18 rows of part (b) reproduce.

One case further: E-254's step, which the submission leaves out (see True?). I replayed it with the same `replay` code and only STEPS changed (scratch copy, 3 s): 2400230, path `F1 F1 F2 F3 R9 R7 R9 F3`, step F3. The path rebuilds with `moves` moves, the tag is `kept`, parent and child are both in half B, q_S is equal (signature (8,1), |det S| = 12, identical histogram), and the planted NO is True.

## True?

1. **"Every recorded key-kept J != 0 step" is 18 of 19, not 18 of 18.** E-254 records a key-kept perI J != 0 step in the F-010 key group: 2400230's depth-8 ball, 1 of 20 210 parents, child joined to the 2400230 ball. E-268 also lists it as one of its "8 key-kept J != 0 edges". The submission's step list (E-233, E-234 x10, E-239, E-244 x6) omits it. The conclusion survives, because I ran the step and q_S is equal and the child is in its parent's component (E-254). The title, the scope line, P3 and the table are still wrong as written.
2. **The "R at v twice" diamond shape is not what the paths show,** and the script never tests shape. Three of the 19 steps do not repeat an R move at the same vertex:
   - E-234 4722 ends `... F1 F8` and steps R1.
   - E-244 11123 ends `... F2 F7` and steps R2.
   - E-254 is F3 after F3.

   E-239 says only "a shared diamond shape (R at vertex 8)", up to relabelling. Either define the shape so that it covers these three and check it in code, or drop it from (c) and from the theorist question in Next.
3. **The equality was predicted, not predetermined.** Equal q_S along the join chains follows from H-023 (J = 0 steps are derived equivalences), which is **OPEN**. E-306's 14 815 constant-q_S steps are ball steps, not the steps of these join chains, which run through the child balls. So the outcome was predicted by an open hypothesis plus an empirical regularity. The same holds for "implied before it was computed" in (c).

## New?

- (a), the count 1 + 1 + 0 = 2, is E-306 Claim (c), which the author acknowledges. The candidate totals 99 and 69 + 354 are E-229 and E-239.
- "Every n = 9 key-kept J != 0 child lies in its parent's J = 0 component" is already in the record:
  - per step in E-234 ("none leaves the J = 0 class of 3060000"), E-239, E-244 and E-255 (8427);
  - collected in E-254's Result ("those of E-239 / E-244 lie in 3304000's own J = 0 component (7 = 1 + 6 steps), so they cannot bridge");
  - restated in E-268 ("the result follows from E-254/E-255 plus component membership").

  The Prior record section lists this as "New". It is not.
- Genuinely new: q_S computed on the 16 deeper steps (17 with E-254), and the per-outcome tally showing that E-306's loop dropped nothing. Both are small, and the first is predicted by H-023 (point 3 above).

Terms grepped: E-233, E-234, E-239, E-244, E-254, E-255, E-268, E-306, H-015, H-023, q_S.

## Evidenced?

The evidence is specific enough to be believed without a re-run: per-depth candidate tallies, every path and step listed, outcome classes enumerated, and times given.

The **Power** line overstates. "Planted cross steps NO 18 of 18" compares each child with the *other half's root* invariant. Each child already has its own half's invariant, so the check reduces to INV(3060000) != INV(3304000). That is one bit, already established by E-306's cross control (disjoint invariant sets), repeated 18 times. It shows that q_S separates the halves. It does not show that the step test, applied to a parent and child, can fire. Say "q_S separates the halves (E-306); repeating the comparison per child adds nothing". The "power on H-015 here is nil" half of the line is right.

## Scope

The 3060000 depth-5 ball was not rebuilt (its 10 steps were replayed), 3060000 depth 6 was not run, and multi-vertex J != 0 steps are not seen. The submission states all of this, which is good. The title needs narrowing. Suggested wording: "At n = 9 (F-010 key group), all 19 recorded key-kept perI J != 0 steps (E-233, E-234, E-239, E-244, E-254) keep q_S, and every child was already joined to its parent's J = 0 component, so q_S could have said NO only if H-023 failed: no recorded step has an outcome H-023 does not predict."

## Required for acceptance

1. Add E-254 (2400230, path `F1 F1 F2 F3 R9 R7 R9 F3`, step F3) to STEPS. Run `replay` (adds about 3 s). Change 18 to 19 in the title, scope, P3 (note that P3 was registered on 18), table and Claim.
2. Remove "R at v twice" as a description of all the steps, or define it so that it covers 4722, 11123 and E-254 and check it in the script.
3. Reword "predetermined" / "implied before it was computed" as "predicted by H-023 (OPEN) and E-306's record".
4. Rewrite the Power line: the planted comparison is E-306's half-separation bit, not 18 independent controls.
5. In Prior record, move "every child lies in its parent's J = 0 component" from New to E-234 / E-239 / E-244 / E-254 / E-255 / E-268.
