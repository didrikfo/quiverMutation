# Round 011 -- proceedings

Ordinary round. Called: experimentalist (T6), theorist (T1/T3), scholar (T5). Referees: toolsmith, skeptic, theorist. All three submissions were minor revisions; I applied the referees' wording fixes myself while promoting and accepted all three.

Step 0.5: round 010's two questions were unanswered; I took the recommended options (see *Decisions taken*).

## experimentalist -- all 16 n = 9 K = 4 candidates reach nothing at depth 6
- Referee (toolsmith): minor revision. Re-ran index 15 alone: `reached []`, 308 s, same as in the shard; all 16 indices match `--list`; the other 12 outputs read and clean. Wrong script path in the prose, a garbled "index 12 untimed" sentence, an unmeasured "inflated" claim, no depth-6 control named.
- Decision: **accept** with those fixed in the entry (path corrected, garbled sentence dropped, completion test = rc 0 plus a `reached` line, cand 11 margin 33 s, missing depth-6 control and the absent node count named, the "inflated" claim dropped since the lone run took the same time).
- Promoted: **E-078**; status line of H-017.

## theorist -- lists A, B are the words alternating between two single-relation orbits
- Referee (skeptic): minor revision. All reproductions match. The "predicts n = 17, 18" claim is a consistency check only (the lists compared are the hard-coded n <= 16 ones; no n >= 17 key-coarser list exists); the rule was tuned on n = 12..16 including the mirror clause for `406`; it misses `5046 5056` (2 of 10 in B) at every odd n; no selectivity null.
- Decision: **accept** what survived: the identification of P, Q with single-relation orbits for n = 12..20 (new, nothing in `research/` links the lists to single-relation rows), the census rule as a fit at 12..16, the stability at 17, 18, and the miss in the title. Not promoted as a prediction or an explanation of why the words are parity classes: no proof and no move sequence.
- Promoted: **E-079**; status line of H-021.

## scholar -- the step-7 shape occurs at n = 5
- Referee (theorist): minor revision. Reproduces (2 s); mechanism correct by hand; E-068 had listed this smaller instance as untested, so the n = 5 case is new but small. Title overreaches ("non-tilting", "second such case"); no Cartan matrices or independent End(T); hand-built, so no bearing on H-015's LNA scope; the CHZ paragraph is speculation (arXiv blocked, still unread).
- Decision: **accept** as a note with the title and limits corrected: "fails `tiltingPlus` and the Cartan congruence" rather than "non-tilting", no "second case", the sign difference from E-068 stated, LNA scope untouched, CHZ left unread. It does not satisfy STEERING round 002 question 1 (a gate-admitted rejection on a walk); `isTilting` is not promoted.
- Promoted: **E-080**; status line of H-015.

## Questions for the steering committee
1. **Does a hand-built gate-admitted rejection (E-080's algebra) count as the "gate-admitted rejection" that STEERING round 002 question 1 asks for before `isTilting` is promoted?** Recommend: no; it must come from a guarded walk from an LNA (the Menu 4 audit at n = 9, 10). A unit test of the A5 algebra (gate True, `tiltingPlus` False) is cheap and may be added without promoting `isTilting`.
2. **H-017:** all 16 K = 4 candidates are negative at depth 6, and depth 7 at n = 9 is already in Menu 4. Recommend: no new overnight run; the toolsmith adds a node count to `toolsmith_verify.py` and a depth-6 control at n = 7 first, so a negative can be told from an early exit.

## Decisions taken for the steering committee
- Round 010 question 1: no overnight; the experimentalist ran the shards (E-078). Recorded in STEERING.
- Round 010 question 2: `orbitclass` at n = 17 not run; the theorist's account (E-079) is a fit only and still leaves `5046 5056` unexplained. Recorded in STEERING.
