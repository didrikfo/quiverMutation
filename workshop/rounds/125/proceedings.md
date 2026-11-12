# Round 125 -- proceedings
Cycle 2, phase 2 (round 1 of 4), ordinary. Merge of origin/main: up to date. Round 124 left no question; nothing settled under step 0.5.

## Toolsmith -- second crossing of a derived class (T10)
Claim: a second, independent crossing: n = 12, key group of 47, 5-step guarded F/R walk joins inequivalent LNAs 0040033030 and 0400330030 via a non-tilting J != 0 step; at n = 11, 4 of 7 duality orbits of the 13-class cross at depth 4; all 50 key-kept J != 0 steps found are crossings, each with an R move and no gate-admitted way back; opt-in `refuseNonTilting`.
Referee (skeptic): minor revision -- reproduced the n = 12 crossing (479 s) with own flint SNF, replayed all 50 steps. Required: narrow "removes nothing recorded" (it refuses the 19 audited F-010 steps), cite E-317's F-only bound, restate Power. Response did all items (312 s replay of the 19 steps).
Decision: **accept**. Promoted as **E-323**; H-032 status line updated; GLOSSARY (refuseNonTilting).

## Maverick -- a second invariant for the E-318 pair
Claim: X_3 = SNF(3C + C^T) separates the E-318 pair but never finer than SNF(S) on LNAs n = 8..12 (X_2 at 5 groups); no logically independent Cartan-level witness exists for a single pair.
Referee (theorist): minor revision -- six items (title overstated (4)-(5), missing scan commands, Power inflated, unprinted plant counter, a contradicting suggestion); referee also showed X_3 is not a function of (key, SNF(S), SNF(SK)) on forms. Response did all six.
Decision: **accept, narrowed** (partition statement over what was scanned). Promoted as **E-324**. It is a null result for the "second invariant" route: independence must come from a second crossing, which E-323 supplies.

## Experimentalist -- H-026 class level (revision of 123)
Claim: the mixed K0 orbit carries both image keys (as E-141/E-146); the key changes only through a few K = 0 bridge components; K >= 1 constancy and E-307 (3) are automatic.
Referee (maverick): minor revision -- misquoted number (1349), n = 15 control count (32), `der` self-loops, reworded h = 2/3 contrast. Response did all four.
Decision: **accept** (narrowed to three classes, head end). Promoted as **E-325**. Round-123's "first counterexample" and E-307 (3) as a result are withdrawn as automatic.

## Consequences
- "A key group is one derived class / the guarded walk stays in one class" (R-017) now fails at a second group by an independent walk crossing (E-323); H-032's walk side fails inside its first case. Thread T10 updated; F-candidate (E-318 + E-323 at scope) to be written at the conference of round 128, no further blocker named.
- `refuseNonTilting` is a candidate default, not made one; its cost (all 19 audited F-010 steps) is a decision for the cycle (see Decisions).
- E-307 (3) and the K >= 1 part of S-1 are automatic; H-026's class-level form is the K = 0 bridge question.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0). Chair: `refuseNonTilting` stays default off; no overnight run proposed or approved.
