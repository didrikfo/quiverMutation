# Round 016 -- proceedings (conference)

No new work; six position statements (`rounds/016/<persona>.md`, haiku model). Nothing promoted to `research/`. Factual sub-claims are the personas' own and unchecked.

## Most promising question, by persona
- experimentalist: why do `3334` and `2455` sit in a different small orbit while the other 4-letter words with a 4 join the `444` orbit at n = 12..15?
- theorist: why does the rule table give P and Q a shift by 2 and never by 1 (E-082)?
- skeptic: does the row-set identity of E-088 hold at n = 16 (4-letter words) and n = 17?
- scholar: is the Cartan congruence the same criterion as `tiltingPlus`, independently of the library defect?
- toolsmith: does the walk ever build a monomial cord member from an LNA at n <= 8 (`MONO=1`, `--plan` first)?
- maverick: why do cord-bearing control members exist at n = 8 and not at n = 6, 7 or among the n = 9 candidates?

## Where they disagree
No contradictions. Weakest claims: skeptic and theorist both name `k(33x) = 2x` / lists A/B as fits without mechanism; scholar names E-086's n = 8 class 2 walk (run under the buggy rewrite); toolsmith names E-089's two-member search; maverick names E-065's range (n = 8..11). Toolsmith and maverick ask each other the same `MONO=1` question.

## Proposed agenda (ranked)
1. **Library fix and re-run** (T5). Toolsmith: patch `reduceAgainstPivots` with a unit test on the congruent pair of E-087. Experimentalist: re-run the E-086 n = 8 class 2 walk under the fix, say which counts change. Scholar/theorist: Cartan congruence against `tiltingPlus` as one map.
2. **Cords at n = 8** (T6). Toolsmith: `MONO=1 --plan` at n = 8, then run if sized. Maverick reads the result.
3. **`3334`/`2455` and `k(33x) = 2x`** (T1/T2/T4). Theorist; skeptic supplies the n = 16, 17 row-set test (`--plan` first) of the `444` orbit.
4. **Parity shift** (T3). Theorist: shift by 2 and not 1 in P/Q.

## Questions for the steering committee
1. **Overnight:** none proposed. Recommend: no, until `MONO=1` at n = 8 is sized.
2. **Agenda:** approve or change the ranking above in `STEERING.md`. Recommend: approve as is.

## Decisions taken for the steering committee
- Round 015, question 1 (`reduceAgainstPivots` fix): toolsmith patches it with a unit test in round 017; experimentalist re-runs the E-086 n = 8 class 2 walk -- decided by the chair of round 016; no answer from the human. (Library change for a defect with a referee-reproduced counterexample; `isTilting` stays unpromoted.)
- Round 015, question 2 (overnight): none; the n = 8 `MONO=1` plan is sized first -- decided by the chair of round 016; no answer from the human.
