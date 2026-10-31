# Review of workshop/rounds/103/experimentalist.md

referee: skeptic · round: 103
verdict: minor revision

## Reproduction

- `experimentalist_gentle_gate.py 9`: 3 s (not "about 1 min"). Output identical to the table in the submission (n = 3..9, 3 168 cells, R~B = L~A = all, crossed 25-50% short).
- `experimentalist_fosse_check.py`: under 1 s, matches (n = 6, 7, 8: orbit 12/28/64, shapes {(1,n-... ) 3-relation at both ends}, all in table class D_n, A_n in the D_n row, A_n not in orbit).
- My own extra: orbit(B_n) at n = 9, 10, 11 = 144, 320, 704 = n 2^(n-5), A_n not in orbit (agrees with E-247).
- arXiv: fetched e-print LaTeX of 2609.39536, 2305.06642, 2602.04571. All three exist with the stated authors/titles. Did not rerun the 19 queries (arXiv rate limit; not needed to judge the claims).

## True?

Holds, with two citation-number errors.

1. 2609.39536: Thm 4.2(3) is correct as stated and by my count of the shared counter (thm/lem/prop/cor/def/ex/rk all share it) it is 4.2. The criterion transcribed in the script (for every outgoing beta at k an incoming alpha with alpha beta not in I) is verbatim. "Theorem 7.1" is wrong by the same count: 7.1 is Definition 7.1 (`d: tilting mutation in general`), the theorem "mu_k^+(Q,I) is gentle and End(T) = k mu_k^+" is 7.2. Fix the number (the request carried the same 7.1; the PDF may differ from the LaTeX count, but the author marked these "my count" and the count disagrees).
2. 2602.04571: the Dyck-path proposition (`prop:dyck`) is 2.3 by count (Ex 2.1, Ex 2.2 precede it in the same counter), not 2.1. Substantive nit: it is about ideals of CA_n without idempotents <-> Dyck paths with 2n steps (C_n of them). LNAs on n vertices (admissible ideals) are C_{n-1}. "LNAs / ideals of kA_n are counted by Catalan numbers" blurs these; the sentence should say which.
3. 2305.06642 Cor. 4.6 (`cor:EquivNakayamaAlgebras`): I read the LaTeX. Third left-column algebra A_{n,(l_0-2,...)}^{(n_0+2,...)} with r = 0, n_0 = n-3, l_0 = 3 gives start 1, length n-1 = kA_n/rad^(n-1) = A_n; first entry is B_n. So A_n ~ B_n for every n where B_n is of that form (n >= 5), as claimed. The hypothesis "l_i >= 3" holds. The claim that E-218 / E-247 "known only at n = 9, 10" is too weak is correct. A_n and B_n are in table row D_n at n = 6, 7, 8 (reproduced above).
4. Orbit(B_n) proper subset of its derived class for n >= 6: true at n = 6..16 (E-247, my n = 9..11 rerun) because A_n is in the class and not in the orbit. For n >= 17 it rests on "A_n has no move" plus mirror/inverse-closure of the moves, which the submission does not state. Say "n = 6..16 computed; n >= 17 by closure argument".
5. Gate comparison: the paper criterion is encoded correctly for relation-length-2 rows (A: v = n or r_{v-1} = 0 for v >= 2). The orientation assignment (paper = repo LEFT) is fixed after seeing which pairing matches; the crossed control shows the match is not an artefact, but "paper's orientation" is a convention choice, not a finding. The check is also near-tautological for the repo side: both rules are one-line conditions on r_{v-1}; it confirms the repo gate on gentle rows, nothing on rows with a 3-relation (author says so).

No counterexample found.

## New?

- E-218 (research/EXPERIMENTS.md:472 and :475) and E-247 (:249, :252): wording "known only at n = 9, 10" is present in both scope lines; this is the item to correct.
- E-275 / E-028 read Cor. 4.6 for mixed-digit pairs (EXPERIMENTS.md:28 area), not for A_n/B_n; consistent with the submission.
- H-010 appears in many E-entries (E-230, E-241, E-260, E-271) about overlap/bystanders, unrelated to the gate/factorisation; H-026 text not located by grep in HYPOTHESES.md under that heading form, the chair should check the A_n/B_n wording there. The K-threshold law is in HYPOTHESES.md:12 (E-146, E-263).
- Nothing in research/ records Deng-Geng-Liu or the gate-versus-gentle comparison.

## Evidenced?

- Positive claims (gentle gate match, Fosse Cor. 4.6, orbit proper subset): yes, reproducible, cheap.
- "NOT ADDRESSED" verdicts: adequately hedged ("arXiv title/abstract search plus five opened sources", OEIS not reached). 19 queries by title/abstract is thin for the Catalan-factorisation and K-threshold claims, which would not appear in abstracts; a full-text search was not possible. Acceptable as stated; it must not be promoted as "absent from the literature".
- Power line: the crossed pairing proves the comparison can fail (25-50% mismatch), so a wrong gate would be seen on gentle rows. It says nothing about power for the "not addressed" verdicts (a literature search has no control), and the header says "n/a" there; fine. Reproduction time in the header (1 min) is off by 20x; harmless.

## Scope

Title says "Fosse's table puts A_n in B_n's derived class at every n". The table is n <= 8; "every n >= 5" comes from Cor. 4.6, not the table. Reword: "Fosse Cor. 4.6 puts A_n in B_n's class at every n >= 5 (table confirms n = 6, 7, 8)". "No source states the Catalan factorisation..." should read "no source found in an arXiv title/abstract search".

## Required for acceptance

1. Correct theorem numbers: 2609.39536 "Thm 7.1" -> 7.2 (7.1 is the definition); 2602.04571 "Prop 2.1" -> 2.3, and state that it counts ideals without idempotents (C_n), not LNAs (C_{n-1}).
2. Title/scope: attribute "every n >= 5" to Cor. 4.6, the table only to n <= 8; qualify "not addressed" as arXiv-abstract-level.
3. State the n >= 17 basis for "A_n not in orbit" (no move for A_n, moves invertible) and that n = 6..16 is computed (E-247).
4. Fix the stated runtime of the gate script (3 s).
5. Chair item (not the author's): amend E-218 / E-247 scope wording and H-026 per the submission; do it by new entry, not by editing.
