# Round 123 -- proceedings
Ordinary round, cycle 2 phase 1 (round 3 of 4). Step 0.5: round 122 left no questions; STEERING.md has no unanswered ones.

## Submissions
**Skeptic -- T10 replay of recorded merges (n = 8-10).** Claim: the 19 n = 9 F-010 joins each contain exactly one J != 0 step (the audited one, by construction), q_S constant at every node; the library's LNA-to-LNA merges contain none, and that zero is forced (the depth-3 key-guarded F-041 balls hold no key-kept J != 0 edge); E-318's walk is flagged. Referee toolsmith: minor revision (zero is forced not informative; cite E-169; planted controls are one comparison; F-037 is 5 of 19; tied meetings none; wording of H-023 test). Response did all six items. **Accept (narrowed to the referee's reading).** Promoted: **E-321**.

**Theorist -- (LF) and (H1) for Construction-2.2 parents.** Claim: (LF) proved given (C) by a vertex shift; (H1w), (H3), H(I) = 0 machine-checked on 82 LNA instances; E-319's (H1) clause 1 not implied and never used. Referee scholar: minor revision (narrow to "given (C)", vertex-shift proof, E-258, tree H = A not checked, vacuous (LF) check). Response did all six. **Accept, narrowed.** Cancellation on no-loop steps is already E-235; the new part is (LF) and (H1w). Promoted: **E-322**; H-031 status line updated (still OPEN); GLOSSARY entry.

**Experimentalist -- H-026 class level at n = 13 (h = 2, 3), n = 15 (h = 3).** Claim: the image key is not constant on a mutation/derived orbit or the q_S class at K = K0. Referee maverick: **major revision** (reproduced and true, but it repeats E-141 and E-146, and it is automatic: moves among K >= 1 rows induce moves on the image, so the key can change only through K = 0 rows; E-307 (3)'s forward-component constancy is likewise automatic). I follow the referee. **Revise**; nothing promoted (the run waits for the revision to be recorded with the restricted-orbit result). The referee's scratch scripts are copied to `workshop/rounds/123/maverick_k1.py`, `maverick_ind.py`, `maverick_path.py` for the revision.

## Consequences
- E-321: no recorded LNA-to-LNA merge uses a J != 0 step, so the recorded F-only merges stand at n <= 10; the guard docstring claim stays false for R / J != 0 walks (E-318). Thread T10 now: a default J != 0 refusal in the library (prunes 0 edges in the F-041 balls, refuses E-318's walk) is the toolsmith's proposed next step; n = 11, 12 merges are not covered.
- E-322 corrects E-319's (H1) to (H1w) and makes H-023's loop joint rest on (H3) for real parents and Hom/End(T). H-023 stays OPEN.
- E-307 (3) (image key constant on forward-rewrite components) is a consequence of the induced-move property, not class-level evidence (referee, unverified at n = 15); noted on the board for the experimentalist's revision.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0). Chair: no overnight run proposed. Round 124 is the phase conference.
