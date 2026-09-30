# Round 005 -- proceedings (conference)

No new work; six position statements (`rounds/005/<persona>.md`). Nothing promoted to `research/`. Statements are from the haiku model and unchecked: treat their factual sub-claims (e.g. that the E-059 odd-n failure *is* an unmerged pair; "n = 13 breakdown" of the signature) as untested.

## Most promising question, by persona
- experimentalist: does an unmerged equal-size middle pair (E-062, `4056`) also explain the odd-n failures of E-059? (Untested guess.)
- theorist: why `k(33x) = 2x`, `d = x - 3`?
- skeptic: is orbit-plus-mirror equal to key class over the `--max-word 4` catalogue at 14..16?
- toolsmith: run `batch.py orbits` over the whole catalogue at 13..16 to give that comparison its data.
- scholar: why does the gate admit E-032 step 7 that Ladkani 2.3(c) rejects; does it repeat at n = 9/10?
- maverick: a Z-lattice/arithmetic reading of the Euler-form signature for H-017.

## Where they disagree
Each names a different "weakest claim": experimentalist and skeptic say parity, from 19 cores of 139; theorist and toolsmith say `k = 2x` is description without mechanism; scholar says H-015's support rests on n <= 7; maverick says "the Coxeter polynomial cannot see (cords, relations)" is an n = 9 observation. None contradicts another. Skeptic's and toolsmith's questions coincide (three personas ask for the orbit/key comparison), so it ranks first.

## Proposed agenda (ranked)
1. **Orbit-plus-mirror vs key, `--max-word 4`, n = 14..16** (T3, T8). Toolsmith first commits the whole-catalogue orbit task with `--jobs 4` and a fsync/append check; experimentalist runs it; skeptic referees. First question: for the 12 cores of E-060 and the 7 of E-059 at n = 14, 15, 16, do orbit-plus-mirror classes equal key classes, and is the size-20300 pair of `4056` the same orbit as the 20300 pair of `3355`/`46 3355 3445` at 16?
2. **Mechanism of `k(33x) = 2x`** (T2, T4). Theorist, with experimentalist supplying `k(c)` for `34x`, `44x`, `45x` at n = 14..17. First question: state T4 (the rule table acts the same at every interior position) precisely and derive `k` for `33x` from it, or say what fails.
3. **H-017 positive control** (T6). Toolsmith, maverick. First question: find a certified quipu-with-relations member from its LNA at its known path depth with `maverick_reached.py`'s search; until it does, "not reached" is not evidence. Skeptic to write the null test for `|R| <= 4` fits requested by the theorist.
4. **H-015 step 7 and the literature** (T5). Scholar, theorist. First question: does Aihara-Iyama Thm 2.32 (arXiv:1009.3370) or the CHZ criterion (arXiv:2509.12983) already decide E-032 step 7? Reading only; no run. The n = 9/10 audit stays parked in Menu 4 (not yet run by the human).
5. **Parity vs unmerged pair across the 127 other cores** (T1, T9). Parked until the n = 12 / 14 censuses of Menu 4 come back.

## Questions for the steering committee
None new. The maverick's Z-lattice idea is speculation, not on the agenda; an approver may promote it.

## Decisions taken for the steering committee
- Round 004, question 1 (H-017 overnight): approved the `n = 9` depth-7 run only; added to `OVERNIGHT.md` Menu 4. To satisfy the budget rule I gave `workshop/rounds/004/maverick_reached.py` a `--budget-hours` flag (stops between LNAs, exits 2; checked at n = 9 with a tiny budget: exit 2). `tests/test_overnight_doc.py` passes. The n = 10 run waits for a positive control.
- Round 004, question 2: round 005 kept as a conference.
