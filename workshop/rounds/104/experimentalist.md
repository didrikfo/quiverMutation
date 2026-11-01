# Conference 104, experimentalist (phase 1, literature)
author: experimentalist · round: 104 · kind: conference · no new work, nothing run
## Most promising question for phase 2
Orbit versus class at n = 6: A_6 ~ B_6 is published (Fosse Cor. 4.6, table row D_6), yet A_6 is not in orbit(B_6) (E-247, n = 6..17). Is A_6 reached from B_6 by ANY tilting mutation (a J != 0 step, or a step the E-187 moves omit)? Why: it is the only question this phase produced that can fail either way at n = 6, where everything is cheap. A path means the move set is incomplete and "orbit = n 2^(n-5)" is not a class size. No path after an exhaustive small search means the gap is real. It also feeds T10 (guard audit) with a concrete pair that needs no join search over the J = 0 premise. Needs the theorist's or skeptic's silting check for the J != 0 steps, since my scripts only walk the moves.
## One question outside the main thread
T6/H-017 (defect constant per key class) on the saved n = 9 data: dormant since round 077, no census, and it is a plain count I can run in one round. It could fail (one key class with two defects). It is cheaper than any join search and independent of the literature.
## Weakest claim the workshop relies on
The T7 conclusion "two-bystander cells stay at overlap 2" (E-230, E-241, E-260, E-271). In E-271, 0 of 4 single-bystander controls could lower overlap, so the 9 cells could not have failed: no power. Close second: the K-threshold law K0 = 3, 4, 5 rests on three lengths at key level (E-146, E-263), and the E-280 literature search found nothing to support or contradict it. "Not addressed" at arXiv title/abstract level is not evidence for it.
## Result closest to a finding
E-280 finding 1: A_n and B_n are derived equivalent for every n >= 5 (Cor. 4.6, referee re-read the LaTeX), while orbit(B_n) is a proper subset of the class (n = 6..17 computed). It still needs one check: the n >= 18 closure argument depends on the moves being invertible and mirror-closed, which neither I nor the referee re-proved. Run A_n's move-out test at one n >= 18 (or prove invertibility) before it is stated for all n.
Runner-up: the gentle gate matches the Deng-Geng-Liu criterion on 3 168 cells. The orientation was fixed after seeing the match and the repo side is one line, so it is a sanity check, not a finding.

## What the literature changed for my thread
- E-218/E-247 "known only at n = 9, 10" is wrong in wording: equivalence holds for every n >= 5. Chair correction pending, by new entry.
- Gate: published only for gentle rows. For r_v >= 3 the single-mutation rule is ours alone. Catalan factorisation, K-threshold and n 2^(n-5): nothing found at arXiv abstract level; OEIS unreachable.
- Theorem numbers are LaTeX counts, PDFs unchecked (scholar).

## Phase 2 should not
- Search arXiv again for the factorisation, K-threshold or orbit count by title/abstract (done, 19 queries, no hit). Only full text or OEIS would add anything.
- Re-run the gentle gate check or extend it as a confirmation; only the r_v >= 3 rows (Ladkani 1001.4765 Prop 2.3(c)) are new.
- Re-derive Fosse Cor. 4.6 for A_n/B_n or the orbit at n <= 17; it is done.
- Add more two-bystander T7 cells before a control that can lower overlap exists.

## Needs from others
- scholar: PDF check of Thm 4.2/7.2, Cor. 4.6 (and its misprint), Prop 2.3; OEIS lookup of 12, 28, 64, 144, 320, 704 and 3, 7, 16, 36, 80, 176, 384.
- theorist/skeptic: an independent tilting/silting test for B_6 -> A_6, not tied to the move set.
- chair: correct the A_n/B_n wording (E-218, E-247, H-026) by new entry.
- toolsmith: nothing now.
