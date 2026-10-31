# Round 102 -- proceedings (ordinary, cycle 1 phase 1: literature)

Merge of origin/main: already up to date. Questions open from round 101: none (step 0.5: nothing to decide). Assignments, referees: `call.md`. All three reviews were minor revision; all three authors answered (step 3.5) and I checked each response against its review list.

## Submissions
- **scholar (T3/T8, H-024)** -- Mixed-digit derived-equivalent LNA pairs exist in print (among sources read) only in Fosse 2305.06642 Cor. 4.3/4.5/4.6; Cor. 4.6 has a misprint in 2 of 8 algebras; 6 pairs n = 9..12, 48 pass Cartan tests; "rank >= l_p+2" is Nikulin Thm 1.14.2. Referee theorist: minor revision (title overclaimed "only in print" and "not Cor. 1.13.3"; E-267 uncited). Response did all 4 items. Decision: **accept**. Promoted as **E-275**; misprint note added to the 2305.06642 literature summary.
- **skeptic (T3/T8, H-024)** -- E-269's bound is exactly n - 2 + lam + 1/lam + slack; small-lam gap is missing input; Salem lists secondary. Referee toolsmith: minor revision (iso140 threshold, "cheaper" unmeasured, title). Response did all 5 incl. 160 cross pairs. Decision: **accept**. Promoted as **E-276**.
- **toolsmith (T5/T10, H-023)** -- Ex. 3.1, 4.1, 4.2 reproduce; Ex. 3.2 does not (KeyError under `allowCyclic=True`); Sect. 6 only steps 1-2. Referee experimentalist: minor revision (16 vs 17, Sect. 6 power, stale F-002, hand-typed targets). Response did all 4 and gave `isoTo` a NO. Decision: **accept**. Promoted as **E-277**.

## Consequences
- **F-002** text and the `allowCyclic` docstring are stale (KeyError, not a loop): note added under F-002; docstring fix is a backlog item (toolsmith, phase 5).
- **H-024** stays SUPPORTED; status line now names E-275/E-276. Still not CONFIRMED: E-269 has not been run on the mixed-digit pairs, and NO at lam < 1.5 has no test case at n <= 12.
- Fosse 2305.06642 Cor. 4.6 misprint: if the PDF agrees, any use of the printed formula is wrong for 2 of 8 algebras (the library's `classPreservingOrbit` is unaffected).
- No H-023 status change: Oppermann prints no loop-of-degree->=1 example, so that joint stays a gap of our Lemma A (E-273).

## Questions for the steering committee
None. (Optional: PDFs of Happel-Seidel, Bobinski-Geiss-Skowronski and Nikulin 1979 would close the remaining unread primaries.)

## Decisions taken for the steering committee
Step 0.5: nothing to decide. No overnight run proposed.
