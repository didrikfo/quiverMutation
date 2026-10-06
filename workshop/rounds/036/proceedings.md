# Round 036 -- proceedings (conference)

All six personas wrote position statements (`rounds/036/<persona>.md`); no research, no referees, nothing promoted.

## Positions
- **experimentalist:** are the 28 rows with d_i >= 3 (all J_i = 0) one orbit or scattered, and does "J_i != 0 => d_i = 2" hold deeper than the cap? Weakest claim: that bound (capped data, no distinct-algebra count, no cyclic test). Needs: toolsmith to enumerate the 28 and walk past depth 8 with a checkpoint; theorist for a bound; skeptic to audit.
- **theorist:** why is d_i = 2 exactly when J_i != 0 on walks (E-129)? Weakest: E-123's "no LNA keys on circuit members" has no base rate. Needs: d_i histogram on long walks (experimentalist); a hand-built d >= 3 with J != 0, or proof d = 2 is forced (skeptic).
- **skeptic:** are the 5 out-degree 3/4 parallel rows (E-127) 2-3 orbits, and do Gamma_i and parallel multiplicity match J? Weakest: the d >= 3 => J = 0 pattern is empirical (E-126's gate bound is silent for d >= 3). Needs: theorist for a proof or an independent mechanism.
- **scholar:** same d_i question as the theorist, tied to long circuits, both-die squares at n = 8 and E-127. Weakest: E-123's key absence is circular (circuit members fail the gate, so carry no key). Needs: toolsmith for the cyclic L2 test; experimentalist for a circuit family with an LNA key.
- **toolsmith:** why does the n = 7 both-die BFS frontier stay at ratio 2.4-2.6 (E-130)? Weakest: "0/44 reached" assumes all 44 targets have non-None keys (unchecked). Needs: key diagnostic on the 44; theorist on the growth.
- **maverick:** why do the 12 core words ending in 3 go I1 at K = 3 and I2 at K = 4 (E-125)? Weakest: "room to move" is a hypothesis; head ends were reversed as digit strings, not with `mirrorRow`. Needs: theorist (H-020 row), toolsmith (n = 12).

Agreement: four of six converge on the d_i / J_i bound; two name the E-123 base rate as weakest. No open disagreement. One caution: the experimentalist's "needs" line reads "J_i = 0 force d_i >= 3" the wrong way round; the question itself is clear.

## Proposed agenda (proposed, round 036; replaces the round-032 agenda)
1. **T5: why d_i = 2 when J_i != 0** (theorist, scholar, skeptic). First question: from the gate bound dim J_i <= d_i - 1 and the walk structure, is d_i >= 3 with J_i != 0 possible on a walk; skeptic tries to hand-build one.
2. **T5: audit the 28 rows with d_i >= 3** (experimentalist, toolsmith). First: enumerate them with distinct algebras, orbits, classes, out-degrees; then a checkpointed walk past the cap.
3. **T5: n = 7 closure** (toolsmith, experimentalist). First: are all 44 targets keyed (non-None)? Then the reverse-direction search from them; frontier growth only if a structural reason is wanted.
4. **T5: the 5 parallel rows of E-127** (skeptic, theorist). First: orbit/mirror partition and Gamma_i against parallel multiplicity.
5. **S-1: last-letter-3 and K** (maverick, theorist, toolsmith). First: derive the I1/I2 split from the H-020 rule table; n = 12 only after `--plan`.
Carried, lower: positive control for E-123's key base rate (scholar, experimentalist); cyclic-quiver L2 case (skeptic, scholar).

## Questions for the steering committee
1. **Agenda:** approve the round-036 agenda above, or change it in `STEERING.md`. Recommend approve.
2. **Overnight:** none proposed; nothing needs more than 10 minutes per command. Recommend none.

## Decisions taken for the steering committee
- round 035, question 1 (agenda): keep the round-032 agenda until this conference -- decided by the chair of round 036; no answer from the human.
- round 035, question 2 (overnight): none -- decided by the chair of round 036; no answer from the human.
