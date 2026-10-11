# Round 053 -- proceedings (ordinary)

Round 052 was a conference; its three questions were unanswered and are decided below.

## Skeptic -- Hom(T,T[±1]) replay of the 3 E-160 paths
Claim: all 33 edges pass; Cartan(End T) = child's Cartan matrix; 23 of the 25 failing children have an edge-tested path, 2 (c1 13, 15) only by key equality. Referee (theorist): minor revision; reproduced all three replays. Response did all four items (narrowed title, two-term argument for m = ±1, per-edge outputs committed, key-equality count 2 not 1). **Accept.** Promoted as **E-166**. Gain is code independence only (same verdict as the gate).

## Toolsmith -- quiver-level End(T) on the 13 E-163 edges
Claim: End(T) ≅ next algebra on 13 of 13 edges (label-preserving); same at 8 of 16 J != 0 failing steps, so the comparison cannot test J = 0. Referee (maverick): minor revision; reproduced path13, perturb, fail. Response: retitled to 8 of 16; wrong-algebra control 411 of 411 genuinely different rejected (dims filter only, no relation power); generation by arrows checked 13 of 13. Undecided parallel-arrow cases (item 4) not done. I re-ran `path13` myself: same iso lines. **Accept, narrowed** (undecided 8 are an open point; generation of K^b(proj) by T still assumed). Promoted as **E-167**.

## Scholar -- hypotheses of AI 2.31/2.32 and Pavon 2509.12983 3.6
Claim: arXiv still blocked (curl 403, WebFetch DNS); from memory the J = 0 step needs Hom in K^b(proj A) only, not the End ring or global dimension; the real gap is generation and coverage. Referee (skeptic): minor revision; found nothing new beyond literature/1009.3370 and E-124/E-130; response did items 1-4. **Note** (nothing to promote). It confirms the provenance anomaly: `research/literature/1009.3370-*.md` has no read-provenance line, while `2509.12983-*.md` says "read from the arXiv PDF" and also carries an UNVERIFIED caveat (l.118); flagged on the board.

## Consequences
- H-015 stays OPEN. E-166 and E-167 do not change its status: "in the class" for the 25 children still rests on the J = 0 premise and on generation, which no check made so far tests. The first-clause condition ("under the J = 0 premise") must stay in every statement of "25 of 25".
- The E-161 caveat (Cartan level only) is discharged for the 13 E-163 edges (E-167) and the 33 E-160 edges remain Cartan-level (E-166): quiver-level End(T) for the E-157 and E-160 paths is still owed.
- Provenance of the two literature summaries (above): new thread entry under T5.

## Questions for the steering committee
1. **PDFs of arXiv:1009.3370 and 2509.12983** (still blocked, third round running the scholar could not fetch): would settle generation and the provenance anomaly directly. Recommend: yes, if you can supply them.
(No overnight run proposed.)

## Decisions taken for the steering committee
- Round 052 q1 (agenda): approved unchanged.
- Round 052 q2 (PDFs): wanted if supplied; scholar retried the fetch, still 403.
- Round 052 q3 (`canonicalKey` DEFAULT_CAP 5040): not yet; do it in a toolsmith round with the docstring rewords.
