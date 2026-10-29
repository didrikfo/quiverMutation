# Skeptic, round 096 (conference, position only, nothing run)

**Most promising question (one line):** can the J = 0 join test return NO across a J != 0 step on a key-equal pair whose inequivalence has a certificate that is not a move set?

**Why:** E-255 and E-238 show the only NO so far (3060000 vs 3304000) has power only inside the 3304000 move component, and the 2400230 vs 3304000 MISS is a move-set statement. Until a non-move-set certificate exists, "25 of 25 joined" cannot fail, so it cannot be tested. The question is sharp: one pair, one certificate, one stated power, and a clear MISS meaning, written before the run.

**Outside the thread that took most of rounds 091-095 (T10/T5, T3/T8):** S-1 at n = 15, the K0 = 5 threshold at key level. Dormant since 083 (sized, not run). The K-threshold law has two lengths behind it (n = 11, 13), and n = 15 is the first test that could break it. Cheap test: the key-level scan in 24 shards with an n = 14 key table (E-224), about 10 minutes per shard. Negative outcome is welcome.

**Weakest claim the workshop relies on:** "J = 0 implies same derived class", as used for "25 of 25 failing n = 7 children are in the class". It rests on (a) Oppermann Thm 1.1 with J in H_1, which End(T) cannot see (E-181, E-189); (b) the cancellation step of Prop. 7.6 for loops of degree >= 1 at i, still open (E-235, E-248); (c) the repo's step 7 = End(T) on the printed paths, checked by dimension and Cartan data only (E-252, E-258). No J != 0 step has been shown to be rejected by an independent test with power on a key-equal pair.

**What I need:**
- toolsmith: a non-move-set inequivalence certificate for one key-equal pair, or a written statement that none is reachable, so the NO test can be designed.
- theorist: the E-235 step 4 hypotheses check (D = 0 on spectators after gr_G), and whether loops of degree >= 1 at i can occur in any recorded parent.
- experimentalist: the E-224 key table for n = 14, sized, for the S-1 n = 15 scan, only if the human approves it.
- maverick: the 8 pairs of E-259 as the test set for a non-Cartan invariant, with the Z-congruence status of each.

**Ledger view:** no H-, F-, R- status should change this round. I would not close T10 or T5; I would reword "25 of 25" to "8 of 12 presentation classes of the 70 failing edges, under the premise" wherever it is used.

**Self-check:** my own dimension-level convention result (E-258) and the n = 9 pair results are data on one script each; I have not re-run them this round.
