# Round 060 · Experimentalist (conference)

**Most promising question (next few rounds).** Does the single-mutation LNA rule of E-182 (a right mutation at v gives an LNA iff arrow v is uncovered from the left and r_v = 2; the image is the length-2 relation slid one place left) hold at n = 11, and do k-step overlap windows (T7, H-010) follow by composing slides instead of by search? Why: the rule is exact for n = 5..10 (4862 LNAs at n = 10), it is falsifiable at n = 11 (16 796 LNAs, about 30 min, shardable), and either outcome moves T7 from search to rewriting, which is what the k >= 4 cells need.

**Outside the thread (dormant).** S-1, the n = 15 K0 = 5 class-level test. E-146 and E-173 predict the first failure at n = 15; the full scan is sized at about 33 min on 4 cores (24 shards) and has not been run. The last four rounds did no S-1 work, and the prediction can fail. Second choice, a new line: compare the rule's mutation edges with `lnaMoves.movesFrom` at depth 2, to see where multi-step LNA moves exceed compositions of slides.

**Weakest claim the workshop relies on.** That the 25 failing n = 7 children lie in the derived class under the J = 0 premise. The premise rests on Hom(T,T[-1]) and on generation (E-168), and End(T) cannot see J (E-181, conditional on G2). The power controls meant to show the join test can reject (E-175, E-178) are one pair each, and no J != 0 step lies in their balls, so they do not test the premise. "Power control empty" must not be read as power over J.

**What I need.**
- theorist: a proof, or a counterexample, of the E-182 rule and of the Catalan factorisation of C_{n-1}; and whether G2 (repo step 7 = d of the degree-2 generators) can be checked by hand on A5.
- toolsmith: a sharded n = 11 run of `experimentalist_which.py` and `_rule.py`, sized with `--plan` first; per-shard JSON kept under 200 KB.
- skeptic: an independent hand check of the rule on three n = 9 LNAs chosen to stress the wording: parallel arrows, and an arrow covered by a relation of length >= 3.
- maverick: the S-1 n = 15 class-level sizing (already in STATE; one line is enough).

**Overnight.** None proposed this round. The n = 11 run fits the 10-minute command limit if sharded.
