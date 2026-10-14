# Round 065 -- proceedings

kind: ordinary. Called toolsmith (T10), theorist (T5), experimentalist (S-1, suggested question, breadth slot). Referees: skeptic on theorist and experimentalist, theorist on toolsmith. Merge of origin/main: already up to date. No questions from round 064 to settle (the human has not changed the proposed agenda; ordinary rounds work from it).

## Toolsmith -- sizing a J != 0 join search (E-192)
Claim: no J != 0 key-kept edge in 128 945 n = 8 edges deep enough to bear on the n = 7 pattern (c0 depth 9, c2 depth 8); 0 in 52 190 shallow n = 9 edges; a search past that does not fit a command. Referee (theorist): minor revision: the "about 38 expected" used shallow classes (matched expectation about 18), the Poisson bound ignores clustering, the implicit n = 7 J != 0 path belongs in the claim. Response did all three. **Decision: accept, narrowed** (sizing on capped samples; says nothing deeper). Promoted as E-192.

## Theorist -- G2 past span (E-193)
Claim: for acyclic parents (parallel arrows allowed) step 7's kernel equals Oppermann's red(Z_1 e_i) as subspaces (Lemma L plus 7 479 random pairs, 0 differences); cycles open, repo's vertex-simple Cartan silently wrong there. Referee (skeptic): minor revision: title overreached (signs of steps 4-6, reading of Def. 3.2 untested), caveat already in E-128/E-171, rho_1 wording. Response did items 1-3 and the optional Cartan check (16 of 16 acyclic-child pairs; 60 cyclic unchecked). **Decision: accept, narrowed.** Promoted as E-193; H-015 status line gets pointers, status stays OPEN.

## Experimentalist -- S-1 item 1 at n = 9 -> 8 (E-194)
Claim: fixed-vertex deletion keeps the class for 37.5-48.9% of 147 681 pairs; no fixed i compatible. Referee (skeptic): reproduces, but **not new**: E-114 has the same rates. Response cited E-114 and stated the differences (exact-class labelling, per-class image counts). **Decision: accept as a refinement of E-114** (note, not a finding). Promoted as E-194 flagged "not new".

## Consequences
- E-193: `cartanMatrix` / `isInIdeal` / `_minimalGenerators` (vertex-simple paths) give wrong answers on cyclic quivers with no warning, while the docstrings say only "out of scope". Not a new caveat (E-128/E-171) but the instance (seed 12, i = 3, a gate-admitted child with a 2-cycle) shows it can happen inside the walk. New housekeeping item on the board (T5): refuse or warn on cyclic input; scan the printed paths for cyclic parents.
- E-192 withdraws nothing recorded; the "about 38 expected" figure from the submission is not used anywhere.
- No F-, no R-; no hypothesis status changed.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 064. No overnight run proposed: the n = 8 depth-9 / n = 9 searches are sized in E-192 as multi-hour or hopeless, and a J != 0 test has no independent endpoint certificate yet (E-192, obstacle 1), so nothing is worth a run.
- Suggested question S-1 taken once (item 1 only); not to be taken in round 066.
