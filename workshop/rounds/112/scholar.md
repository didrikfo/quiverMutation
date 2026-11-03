# Conference 112, scholar: phase 3 (research) view

author: scholar · round: 112 · kind: conference · no new work, nothing run. Read: my notebook, R110 submission and review (I have no R111 submission; my R111 work was refereeing toolsmith, E-295), E-289..E-296, STATE, plan. Formed before the board.

## Most promising question
Can a key-guarded walk leave a derived class through a J != 0 step (T10, H-015)? Is any of the 13 key-kept J != 0 steps of E-201 (n = 7 class 1) a derived equivalence? E-295 shows the Cartan route cannot decide: key-kept implies congruent, and the test has no NO case at n = 7. So the question is now a Hom question, not a lattice one. The literature already says what to expect: Aihara-Iyama 2.32/2.50(c) give silting-not-tilting steps, and for those End(T) need not be derived equivalent to A; nothing in the sources says a key-kept step must be tilting. Why this one: it is the only open item that decides the R-candidate "a guarded walk stays in one class" and the library docstring (search.py l. 322) hangs on it. A Hom-checked 13 (Hom(T,T[m]) = 0, m != 0, T generates) either gives a counterexample or removes the worry at n = 7.

## Outside the dominant thread (T10/T5/T3 got the phase)
H-025: Phi periodic but not FCY (converse of 1310.1557 Thm 2(a) is not a theorem). Dormant since the plan put it "literature only". It costs a reading slot, and a periodic-not-FCY LNA is either in the literature (Chan-Koenig-Liu, 2012.11927, 2512.19873 are on file) or a cheap search over the LNA data we already hold (Phi charpoly is cyclotomic, then test the dimension of Ext growth). A counterexample is a publishable fact, not bookkeeping.

## Weakest claim the workshop relies on
That a key-kept step is "same class": every "guard stays in class" statement, and the YES side of E-269, rest on class = genus (Nikulin Cor. 1.13.3), which I still cannot match in any primary (only 1404.5460 l. 336, secondary). The NO side is better off: Ladkani Cor. 3.13 (stated) plus Lemma 3.14's proof give equivariance, and E-292's 19 of 19 NOs sit on a q_S certificate. Also unchecked by me: LNAs have finite global dimension (standard for acyclic quivers with relations, but I wrote "UNCHECKED" and nobody has closed it).

## Closest to a finding
E-296: the K-threshold law, largest splitting K = (n-1)/2 - 2, now at n = 11, 13, 15, 17, 19 (lone-3 class, key level, head ends). One check it still needs: n = 19 (and the n = 17 K = 5 cell) re-run by someone other than the author, since both rest on one 518 s run; then the out-of-sample prediction K0 = 8 at n = 21. A planted-error control on the K0-1 split would give it a NO that can fail. Even then it is key level only, not class level, and not published.

## What I need
- toolsmith/theorist: for the 13 steps, a Hom-verified End(T) (E-290's Hom control shows the harness can say NO) so I can state the result against Aihara-Iyama.
- maverick: someone to re-run `run 19 7`; I can then check whether (n-1)/2 - 2 matches any published Nakayama vertex-deletion count.
- chair/human: PDFs of Nikulin 1979 (Cor. 1.13.3), Happel-Seidel, BGS, Anick 1982; arXiv does not have them, and four UNVERIFIED citations in STATE depend on this.
- I owe: Bondal-Kapranov / Huybrechts citation for Serre-functor intertwining; math/0702060 literature file.
