# Round 060 -- conference proceedings

All six personas wrote position statements (`rounds/060/<id>.md`); nothing was run. Five of six (experimentalist excepted) name the same single most promising question: **G2**, whether repo step 7 of `procedure` equals d of Oppermann's degree-2 generators (1504.02617 Thm 1.1 / Ex. 4.2), checked by hand on A5 with two overlapping relations. All six name the same weakest claim in some form: "25 of 25 failing n = 7 children are in the derived class", resting on the J = 0 premise, since End(T) cannot see J (E-181), G2 is unproved, and no join test has returned NO on a path through a J != 0 step.

## Ledger
E-entries since the last conference (E-174 .. E-182):
- **H-015** (guard sufficient): stays OPEN. E-174/E-177/E-180 give the quiver-level End(T) check power against wrong algebras but E-181 shows it cannot test J; G2 and the premise are unproved; key-guard law is class 0 only.
- **H-010** (overlap reducible only at an end): stays SUPPORTED. E-182 (single-mutation LNA rule, n <= 10) refines the "2 per LNA" count (an average); no proof; no status change.
- **H-017**: stays OPEN; untouched since E-172 (no relation-sensitive handle).
- **H-020**: stays SUPPORTED; untouched. **H-021**: stays OPEN; untouched.
- **T1 (H-021')**: closed as stated (E-176: no linear letter rule for k(c); reopen only for a drift-based rule predicting a held-out core). **T2 (k(33x) = 2x)**: closed as stated (E-176; the 444 orbit belongs to H-020 / T3). No H- text rewritten this round (header rewrite of H-021 stays a housekeeping item).
- No new finding (F-), no retraction (R-).

## Proposed agenda (ranked)
1. **T10/T5 G2 by hand (scholar, theorist; skeptic replays).** Scholar: compare repo step 7 with Oppermann Section 3 / Ex. 4.2 on A5 with two overlapping relations, after re-fetching the LaTeX (Oppermann numbers UNVERIFIED until matched). Toolsmith supplies the A5 step-7 output. Theorist: is G2 derivable from Thm 1.1 alone; Cartan transpose convention (skeptic).
2. **Can the join test ever say NO across a J != 0 step?** Toolsmith: a module-level Hom(T,T[-1]) table on the 25 failing steps and the 13 E-163 edges (Cartan-level vs module-level); a positive control whose only short path uses a J != 0 step (skeptic's request).
3. **Outside T10: Candidate C (T3/T8).** Maverick, scholar reads 0911.5137 Cor 1.9 and 1310.1557 2.9 first: what the periodicity obstruction leaves open; Phi-orbit lengths of the n = 10 group. Level: idea. Skeptic notes it may be vacuous if Phi^a fixes a pair.
4. **T7 rule at n = 11 (experimentalist + toolsmith).** E-182 rule at n = 11 (16 796 LNAs, ~30 min sharded; `--plan` first); k >= 4 by composing slides; theorist proof/Catalan factorisation.
5. **Breadth/housekeeping:** T6/H-017 16 unplaced n = 11 LNAs (dormant since 050; no search that sees relations); S-1 n = 15 K0 = 5 scan stays parked (~33 min on 4 cores, not approved).

Disagreement: none on priorities; experimentalist alone puts T7 first, others G2. The agenda keeps G2 first and gives T7 the breadth slot.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- No unanswered questions from round 059. No overnight proposal (n = 11 rule run and S-1 n = 15 are each ~30 min, shardable under the 10-minute command limit; not approved).
- Proposed agenda set as open threads; ordinary rounds work from it until the human changes it in `STEERING.md`.
