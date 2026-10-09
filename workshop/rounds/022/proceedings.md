# Round 022 -- proceedings

Ordinary. Worked: experimentalist, theorist, toolsmith. Referees: theorist (on experimentalist), skeptic (on theorist), scholar (on toolsmith). All three: minor revision. I accepted all three with the scope fixes the referees asked for, written into the E-entries' Limits.

## experimentalist -- long-square control (T5)
Claim: 0 of 479 761 tilting steps at n = 5..7 have E-097's long-sided square; all 2 104 rejecting (parent, v) (n = 6 c0: 1 842, n = 7 c0: 262) do; strict A5 is not selective (0.5-1.3 % of tilting steps).
Referee (theorist): counts reproduce (n = 5 c0 identical; n = 6 c0 at `--maxexp 6000`: 0/22 140 and 102/102). Narrow the headline to the two runs with rejections; state the A5/long-square discrimination is n = 6 only; say whether `alg.rels` equals `relationsFrom(alg)`; save the n = 5 c0 output.
Decision: **accept**, narrowed. Promoted as **E-100**. The `alg.rels` question and the unsaved n = 5 output are left open in the entry; nothing was re-run by me.

## theorist -- cord criterion (T6)
Claim: depth-1 rule D1 (a relation of >= 3 arrows that is not blocked), peeling depth 1 + min(a, b) for single-big-relation LNAs; "within 3 steps" is false at n = 10 (`22230222`, depth 4); "cord" means a cycle member, not a GLOSSARY cord.
Referee (skeptic): n = 10 counterexample, D1 at n = 8, 9, peeling at n = 8 reproduce; an out-of-sample n = 11 sample (301 LNAs) gives 0 mismatches. D1 was selected among four variants by fit; the max-depth formula is exhaustive only for n <= 10; E-099's title is n = 8, so only the quoted depth clause is refuted.
Decision: **accept**, scoped as the referee asked. Promoted as **E-101**; glossary term "cycle member". No proof of either direction; the blocking step has no derivation. I did not edit E-099 itself (its title is n = 8); E-101 says it refines it.

## toolsmith -- replay and Cartan check (T5)
Claim: the 10 E-084 key-moved parents are recoverable from `rounds/014/scholar_walk_n8_c2.txt`; under the fixed library all 10 keep the key and pass the new opt-in `checkCartan` (env `QM_CHECK_CARTAN=1`); the pre-E-089 reduction fails on 10 of 10; cost about 15 ms per step, off by default; 2 tests.
Referee (scholar): replay, old-reduction failure, overhead (ratio 12 on the check alone) and 36 tests reproduce; the diff leaves the default path unchanged. Fix: the 1.5-2x versus 2.7x overhead wording, "closed" overreaches (Cartan matrices only), replay time is 3 s not 12 s.
Decision: **accept**. Library and tests committed (`procedure.py`, `tests/test_procedure.py`; `import os` is present at line 61, so the NameError the theorist hit with the uncommitted edit is not in this tree). Promoted as **E-102** with the corrected wording.

## Questions for the steering committee
1. **Agenda:** keep the round-020 agenda (recommend yes). Items 1 (T5) and 3 (T6) were worked this round; next round suits the skeptic on the "iff long square" claim off the walks, the scholar/theorist on a step-7 derivation of the long-square criterion, and the experimentalist on T1/T2 n = 15..17.
2. **Overnight:** none (recommend). Nothing this round needs more than 10 minutes per command except the n = 12 depth-5 cycle check (about 8 min).

## Decisions taken for the steering committee
- Round 021, question 1 (agenda): kept the round-020 agenda -- decided by the chair of round 022; no answer from the human.
- Round 021, question 2 (overnight): none; the non-MONO L = 5 walk is not worth 4.5 CPU-hours -- decided by the chair of round 022; no answer from the human.
