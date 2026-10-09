# Round 014 -- proceedings

Ordinary. Worked: experimentalist, maverick, scholar. Referees: theorist, toolsmith, skeptic. All three: minor revision. The fixes were wording and qualifiers; I applied them in the promoted entries and accepted all three.

## Questions of round 013 (step 0.5)
No answers from the human. Both decided by the chair (recorded in STEERING): n = 17 lists not overnight, experimentalist saves the n = 17 output first (done, E-083); H-017 depth 7 not overnight, an n = 8 non-hereditary control first (done, E-082).

## experimentalist -- saved n = 17 run; 4-letter scan
Claim: `5046`, `5056` each have two closed orbits at n = 17 (122 673 / 54 266); 4-letter words with a 4 join the one big `444` orbit at n = 12..15. Referee (theorist): no error; scan reruns at n = 12, 13 byte-identical, n = 17 read from saved output not recomputed; identity with the `444` orbit confirmed by row membership at n = 13 only; "rigid" undefined; shared ledger. Decision: accept with those qualifiers (confirmation of E-080, a slice of E-075/E-079). Promoted: E-083.

## maverick -- n = 8 control with relation-bearing sources
Claim: control finds its source in every run at its depth, none one step short, 7e3-4e4 nodes at depth 6. Referee (toolsmith): reproduced LNA 0 at L = 6 exactly; members are the deterministic head of a sort, not a sample; no cords; inconsistent size ratios (1.3-3.2 for 1-relation, 7-9 for the 4-relation member); 4-relation claim is one member. Decision: accept, claim reworded as "every run", ratios reconciled, `_L6_plan.txt` (mistaken run) deleted. Promoted: E-082. Narrows E-081's non-hereditary limit; cords limit open.

## scholar -- reachability of a gate-admitted, non-tilting mutation
Claim: yes at n = 6 (distance 8) and in 10 sampled classes at n = 7..9; the guard refuses all of them; 0 of about 1.3e6 guard-admitted steps fail `tiltingPlus`. Referee (skeptic): reproduced the n = 6 rejection, n = 7 replay and the test; A5-shape "all of them" unchecked by script, non-tilting rests on `tiltingPlus` alone, the n = 8 class 2 loose end weakens "guard is redundant", totals mix runs. Decision: accept with those qualifiers moved into the claim; E-057's "none" annotated. `isTilting` not promoted (the round-011 condition asks for a rejection *on a guarded walk*; these are guard-refused). The new unit test `tests/test_gate_without_tilting.py` is kept (2 passed). Promoted: E-084; H-015 status line.

## Questions for the steering committee
1. **`isTilting`** (T5): E-084 shows the gate alone is unsound and the guard is what protects the walk, but no *guarded* walk ever reaches a failing step. Recommend: still do not promote; scholar/theorist compute the Cartan congruence on the replayed parents and explain the n = 8 class 2 loose end first.
2. **Overnight:** none proposed. n = 17 key-coarser lists and H-017 depth 7 stay unrun. Recommend: no; toolsmith first builds an n = 8 control member with cords (arrows >= n) and relations >= 1.

## Decisions taken for the steering committee
- Round 013, question 1: n = 17 lists not run overnight; experimentalist saved the `5046`/`5056` output instead (E-083).
- Round 013, question 2: no H-017 depth-7 run; n = 8 non-hereditary control built (E-082).
