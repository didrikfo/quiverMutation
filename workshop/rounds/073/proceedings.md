# Round 073 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. Round 072 left no questions; STEERING had no new answer.

## Toolsmith (T10/T5) -- E-212, accept, narrowed
Claim: at n = 7 class 2, 8 of the 9 E-151 key-kept J != 0 steps have parent and child joined by Hom-verified F/R chains with isomorphism at the meets. Referee (skeptic): minor revision; re-ran the collector and the chain mode for records 0-5, 7, 8 (all YES, same counts); no error found; NO power thin. Response did the required items (retitle, `parallel-at-meet` listed, E-157 in the prior record); the n = 9 F-010 5 + 5 NO control was not run (cap). Decision: accept, narrowed to the eight records; record 6 undecided and class 1 not run go on the board. Promoted as E-212 with the thin-NO-power caveat stated in the entry. Not promoted: any statement that the verifier "can say NO" at n = 7.

## Theorist (T5) -- E-211, accept
Claim: the guard does not force D(i) acyclic; the whole-quiver `noCycles` test does; one finite-gd cyclic parent (n = 4) on which the repo's step 7 drops rbar.b.c and the child is not derived equivalent. Referee (scholar): minor revision; re-ran (1.5 s), recomputed the Coxeter polynomials and the repo child's Cartan by hand; four required items (order label, gd finiteness, wording, narrowed title), all done. Decision: accept. Promoted as E-211.

## Experimentalist (S-1) -- E-210, accept
Claim: the E-190 pair is outside the D >= 3 rule's domain; unrestricted deletions land in P^(6)_(1,1) or the hereditary class; exact `classifyLength(10)` sized, not run. Referee (maverick): minor revision; reproduced the table, tested the exit-2 path at n = 8 four times; required retitle, convention, loop fix (exit 124), "unknown" cost: all done; item 5 (orbit + mirror join of 80000000 and 00000300) goes to the next round. Decision: accept. Promoted as E-210. (Referee wrote no file; the chair saved its returned text as `experimentalist.review.md`.)

## Consequences
- **Docstring / guard (T10):** `search.py` guard docstring: the guard is inert on cyclic input and a None key passes; the acyclic-parent property comes from `noCycles` (l. 377, 420). Docstring reword added to the housekeeping list. F-002's "the LNA search never reaches it" is consistent; no F-/H-/R- status changes.
- **E-199 (T5):** its identification "repo step-7 kernel = Oppermann kernel whenever D(i) acyclic" stands; the cyclic-D(i) case is now a counterexample, not open. Direct callers of `mutateAtVertex` (E-157-style joins) have no `noCycles` test; toolsmith to add a refusal or a check that the printed paths are unaffected (E-198 scanned them: acyclic).
- **E-157 "joined to an LNA of their class":** parent and child can reach different seeds (E-212); the T10 summary reads accordingly.
- No H-015 status change (stays OPEN).

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 072 (it left none).
- Overnight: none approved. Two proposals recorded, neither added to `OVERNIGHT.md`: exact `classifyLength(10)` (experimentalist, E-210; cost unknown, resolve step unsized) and the n = 9 F-010 5 + 5 NO control (`toolsmith_f010_overnight.py`, `--budget-hours`, exit 2 tested; the 4 + 4 control already ran). Reason: conservative default; each first needs a size estimate (resolve step) or a reason the bound would change a decision.
- Suggested question S-1 taken this round (experimentalist), not taken in 071 or 072.
