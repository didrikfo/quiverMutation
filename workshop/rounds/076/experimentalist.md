# Experimentalist, round 076 (conference, nothing run)

**Most promising question:** Can a tilting-only J = 0 join search reach `00000030` from `90000000` at n = 10 (equivalent, but not joined by orbit + mirror moves, E-218 style), and can the same search return NO on any pair whose inequivalence has an independent certificate?

Why: the J = 0 premise is the live question (T10/T5). Every "25 of 25" or "12 of 12" join so far is a YES; no join test has returned NO across a J != 0 step. The E-190 pair is a known equivalent control that the orbit + mirror moves miss, so it tests whether the join search reaches what is true. The NO half needs a certificate that does not share the Cartan-level blind spot. Size it with `--plan` first.

**Promising line outside the thread that took most of rounds 064-075:** The S-1 K-threshold prediction that K0 = 5 first fails at n = 15 (E-146; E-173 confirmed it at key level on a one-row witness). The full class scan was sized at about 33 min on 4 cores (24 shards) in round 055 and never run. It would turn a prediction into a table. Not approved; I propose it for `OVERNIGHT.md` Menu 4 as sharded runs, each under 10 minutes.

**Weakest claim:** "25 of 25 failing n = 7 children are in the derived class" (and its successors E-213/E-216: 12 of 12 representative classes joined), which rests on the J = 0 premise that a tiltingPlus step is a derived equivalence. End(T) cannot see J (E-181, E-189). The Hom(T,T[+-1]) test is Cartan-level. Second-weakest: the S-1 K-threshold law, which rests on two lengths.

**What I need:**
- toolsmith: `--plan` for a tilting-only join of `90000000` and `00000030` at n = 10, with the depth bound written down; then the search's answer on it.
- skeptic: a key-equal, inequivalent pair at n = 9 or 10 with an independent certificate that has power on a J != 0 step, or a statement that none is available.
- toolsmith or maverick: shard plan for the n = 15 K0 = 5 class scan, each shard under 10 minutes.
- theorist: whether the E-211 cyclic counterexample is excluded by anything beyond E-214's hypotheses (the premise needs it).

No file other than this one was edited. No run was made.
