# Round 012 -- proceedings (conference)

No new work; six position statements (`rounds/012/<persona>.md`, haiku model). Nothing promoted to `research/`. Factual sub-claims are unchecked. The maverick's "contradicted by Euler form working" is not shown: E-063 covers n = 8..11 and separates outside-every-quipu-class, not (cords, relations). The toolsmith's and skeptic's descriptions of E-077 repeat its own caveats.

## Most promising question, by persona
- experimentalist: why do `5046`, `5056` fail the alternating pattern at odd n, and can failures be predicted from orbit structure?
- theorist: does the move graph of the single-relation orbits P, Q force the parity alternation (move sequences `35@0 -> 35@2`, `4@0 -> 3@1`, obstruction for `5@0 -> 3@2`)?
- skeptic: is letter 4 special in the rule table, or is the 0.55-versus-0.02 merge rate a side effect of the `34` collapse?
- scholar: is there an LNA-reachable gate-admitted non-tilting mutation at n <= 9?
- toolsmith: make the H-017 bounded negatives legible (node count, depth-6 control at n = 7).
- maverick: is the Euler-form signature a lattice invariant that classifies quipu membership?

## Where they disagree
No contradictions. Weakest claims overlap: experimentalist, skeptic, toolsmith and theorist all say E-077's lists A/B are a fit with no mechanism and no real test at n >= 17; scholar says the "walk-reachable only" policy for `isTilting` is a policy, not a result; maverick says the Coxeter-polynomial negative is n = 9 only. Three of six (experimentalist, theorist, toolsmith) ask the theorist for the same move-sequence account, so it ranks first.

## Proposed agenda (ranked)
1. **Why lists A/B are parity classes** (T1/T3). Theorist with experimentalist. First question: move sequences `35@0 -> 35@2` in P and `4@0 -> 3@1` at n = 12, and why `5@0` cannot reach `3@2`; then `5046 5056` at odd n (`--plan` at n = 17 first). Skeptic referees.
2. **H-017 legibility** (T6). Toolsmith: node count and an n = 7 depth-6 control in `toolsmith_verify.py`. Maverick's Euler-form idea waits for it.
3. **Walk-reachable non-tilting mutation** (T5). Toolsmith with scholar; theorist on the one-map identity. A5 unit test allowed.
4. **Letter 4 versus the `34` collapse** (T2/T4). Theorist, skeptic. Merge rate stratified by `34` versus other 4-words; orbit-collapsed count of E-075.

## Questions for the steering committee
1. **Overnight: the n = 17 key-coarser lists** (`toolsmith_orbitclass.py`, 2+ h) would give the first real test of E-077 beyond n = 16. Recommend: not yet, until the theorist's account (agenda 1) says which words to expect.
2. **Agenda:** approve or change the ranking above in `STEERING.md`. Recommend: approve as is.

## Decisions taken for the steering committee
- Round 011, question 1: no; hand-built rejection does not count towards promoting `isTilting`; an A5 unit test is allowed. Recorded in STEERING.
- Round 011, question 2: no new overnight; toolsmith adds node count and n = 7 control first. Recorded in STEERING.
