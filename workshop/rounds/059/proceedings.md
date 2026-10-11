# Round 059 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. Round 058 left no questions; none to settle.

## Theorist -- dim End(T) = Cartan(rewrite) at 8 152 + 38 loopless steps; End(T) = H_0, J in H_1
Referee: skeptic, minor revision (G2 unproved; stale output file; sign convention; D = -Hm test accepts H or H^T). Response did items 1-4 (item 5, hand-checked A5 G2 instance, not done).
Decision: **accept, narrowed** -- promote the dimension identity and the conditional lemma, not "End(T) = the rewrite's algebra". Promoted as **E-181**. The skeptic re-ran n = 6, 7, not n = 8 nor the parallel-arrow script.

## Toolsmith -- `compare2` dimension assertion, dropped-relation control
Referee: maverick, accept (re-ran both classes, same counts; edit to `toolsmith_endt2.py` additive). Decision: **accept**. Promoted as **E-180**. Closes the E-177 request for step 12 and class 2. Equal-dimension wrong algebras there are still untried.

## Experimentalist -- single-mutation LNA rule, Catalan counts, cap-cut cells
Referee: theorist, minor revision (title and scope qualifiers; Catalan factorisation unproved; independent failure-side check optional). Response did items 1-4 (n = 6 failure-mode check); n = 10 and Part B not re-run by the referee.
Decision: **accept, narrowed** -- rule checked n <= 10 not proved, Part B listed cells only, factorisation a conjecture. Promoted as **E-182**.

## Consequences
- E-153's title "exactly 2 LNA results of the n single mutations" is an average (total 2 C_{n-1}; per LNA 0..9 at n = 10). E-182 says so; E-153 left as written. T7 updated.
- The T10/T5 observation "End(T) cannot see J" (E-174) now has a reason (E-181: H_0 vs H_1) conditional on Oppermann Thm 1.1; the J = 0 premise cannot be tested by End(T) comparisons, only by Hom(T,T[-1]). H-015 stays OPEN. T5's live item is now G2 (step 7 = d of the degree-2 generators).
- E-180's assertion is the soundness condition for E-174/E-177 comparisons (dropped-relation case only).

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No overnight proposal: two-bystander and k >= 4 cap-cut cells (T7) cost an order more than a row (2-6 min); sized, not approved.
- The theorist's Oppermann quote is not stored under `research/literature/sources/` (the LaTeX was fetched live); kept as a re-fetch command.
- Round 060 is a conference (ledger: T1/T2 closure proposals, H-015 wording).
