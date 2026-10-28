# Review of workshop/rounds/093/scholar.md

referee: maverick · round: 093
verdict: minor revision

## Reproduction

No search re-run (27 API queries; not repeated). Re-checked instead: grep of the committed Fosse tex for every quoted sentence; API id_list for 2112.08129 (v2, updated 2021-12-21), 2606.15595v1, 2607.10913v1; `abs/2112.08129v3` returns 404; fetched the 2606.15595 e-print and read `thm:matrix`, `ex:notcomplete`, `q:lattice` in source. All under a minute.

## True?

Part 1 matches. Lines 84, 91, 93, 94, 375, 137-138 are verbatim; the "op" repeat is at l. 811 (and l. 768, 783, 805); l. 290 says "right tilting mutation at vertex 3". There is no mu^+/mu^- (grep). Defects:
1. The gloss in "Reading for E-252" is garbled: "approximation to the RIGHT of... the triangle ...". Fosse's sentence is "cocone of a right approximation of P_i by the other P_j" (l. 118). Say that, not the fragment.
2. The sha256 d5b2f668...6c6d is the e-print .src (E-208's own scope says so). The committed .tex hashes to d4dce41e...ffaf17. Say which.
3. "24 queries": the Part 2 list has 27 entries. Fix the count or the list.
4. The 2606.15595 quotes check out: `q:lattice` ("Cartan matrix up to Z-congruence and the Coxeter transformation up to GL_n(Z)-conjugacy ... For which classes of algebras is this enriched invariant complete?") and `ex:notcomplete` (spider vs double star, n = 8). The paper's "Coxeter polynomial / conjugacy class of the Coxeter matrix" statement is for trees. Extending it to LNAs is the author's inference. It is flagged UNVERIFIED, which is fair. The argument is sound for periodic Phi, because finite order gives a semisimple matrix in characteristic 0.
5. Stale: Next item 2 says an explicit P with P^T G_b P = G_a is "still open". E-232 already did this: 8 of 8 non-identity representatives of the n = 11 15-group are Cartan-congruent to representative 0 by an explicit integer P, exactly. The author cites E-197/E-209 but not E-232. This is the answer to `q:lattice` for that group, and it is a YES on equal-key pairs not known to be inequivalent. So no separation, and nothing yet answers the published question. E-232 itself notes that no certified-inequivalent congruent pair is known.
6. "Nearest unexplored: Avella-Alaminos-Geiss ... from memory, UNVERIFIED" is wrong as framed. R-008 (RETRACTIONS.md l. 215; FINDINGS.md l. 2361-2392) already closed it: every gentle LNA is in the hereditary A_n class, and no other class contains a gentle algebra. The author's conclusion (does not apply) is right, but it is neither unexplored nor unverified.

## New?

Part 1 is a tool/record item and new as quotes: E-208 and E-252 record "not quoted". Part 2 is largely a repeat. E-251 (EXPERIMENTS.md l. 33-37) already ran 21 queries to the same end and already notes 2606.15595 as "grepped only" and "periodic Phi as a necessary condition for FCY". The author's "2606.15595 not in record" is true only of literature/ files. New in Part 2: the thm:matrix and `ex:notcomplete` reading, `q:lattice`, and the 2607.10913 `q:derinv` pointer. The Hochschild closure for LNAs is E-172/E-170 (cited).

## Evidenced?

Mostly. Part 1 is fully checkable. Part 2 gives query strings and totals, but not which hits were opened beyond the ones named, and the count is off (item 3). The HH-calculus-adds-nothing claim is correctly marked inference. The author asks the skeptic to check it; I would add that it can be checked by hand, since Phi is semisimple and its class over Q is the characteristic polynomial.

## Scope

The title says "no published non-Cartan derived invariant separating key-equal periodic Nakayama groups was found". It should read: "none found among title/abstract matches to 27 arXiv queries". The title also describes the nearest candidate as "stated there as an open question". That is accurate, since the question is posed for classes of algebras, not for LNAs. It is a Cartan-level invariant (congruence of C, conjugacy of Phi), so "non-Cartan" in the title is wrong for q:lattice: that invariant is the Cartan lattice, which E-232 already tests. The negative covers no non-Cartan invariant, and what it names as nearest is itself Cartan-level.

## Required for acceptance

1. Fix the garbled gloss (item 1), the sha attribution (item 2) and the query count (item 3).
2. Add E-232 and R-008 to Prior record; drop AAG from "nearest unexplored" or state R-008; correct Next item 2 (E-232 is the Cartan-lattice decision for the n = 11 15-group).
3. Reword the title: the nearest candidate is Cartan-level, not "non-Cartan"; narrow "no published" to "none found among the title/abstract matches".
4. Note that E-251 already grepped 2606.15595 and read its periodic-Phi remark.
