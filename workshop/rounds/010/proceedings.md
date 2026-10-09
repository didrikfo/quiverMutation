# Round 010 -- proceedings
kind: ordinary. Worked: toolsmith, skeptic, experimentalist. Referees: theorist, experimentalist, skeptic (one each). Questions of round 009 settled by the chair (below).

## toolsmith -- `toolsmith_verify.py` (tool)
Claim: `--list`, `--cand`, `--budget-hours` added to the round-004 search; E-074's depth-5 negative reproduces; a K = 1 candidate at depth 6 reaches nothing in 434 s. Referee (theorist): minor revision; reproduces (440 s), search code diffed unchanged; four wording fixes (K = 1 to K = 4 index map, machine-dependent budget test, `--cand` budget message, ratio rests on two pairs). Decision: **accept as note/tool**; I applied the four fixes in the entry. Promoted: E-075.

## skeptic -- neighbour-aware null for "a = 4 special" (negative, partial)
Claim: only `444` merges among `aaa`, a = 3..9, n = 12..15, but words with a 4 merge at 0.55 against 0.02 (no 4, no 2), so the `34` route is not singled out. Referee (experimentalist): minor revision; scan and stats reproduce exactly, n = 16 probe agrees. Required: orbit-sharing words collapsed, headline weakened, `34x` outside E-073's criterion. Decision: **accept with the weaker headline**; the collapse of merged words by orbit was not done (recorded as a limit and a request). Promoted: E-077; H-021's criterion stays out of its text (STEERING q2).

## experimentalist -- key-coarser lists at 14..16 (result)
Claim: lists A (9 words, even n) and B (10 words, odd n) repeat across n = 12..16; `348`/`349` size-20300 pairs at 16 are the `4056` orbit and its mirror (row sets identical). Referee (skeptic): accept; re-ran the script byte-identically and the ledger reads. Decision: **accept**. Promoted: E-076 (title carries the scope limit, as the referee suggested).

Status lines updated: H-021, H-017.

## Questions for the steering committee
1. **H-017 depth 6 for the other n = 9 candidates:** 12 untimed at K = 4 (indices 0, 4, 8 done or equivalent). Recommend: no overnight run; experimentalist runs the shards in chair-slot commands (`--cand I`, one per command) in round 011.
2. **`toolsmith_orbitclass.py` at n = 17 (about 2+ h) to test lists A/B further:** recommend not yet; the theorist first explains why those 9/10 words are parity classes, which may predict n = 17.

## Decisions taken for the steering committee
- Round 009 question 1: toolsmith adds the candidate-index argument and `--budget-hours`; no overnight (done, E-075).
- Round 009 question 2: the `aax` criterion stays out of H-021's text (supported by E-077).
