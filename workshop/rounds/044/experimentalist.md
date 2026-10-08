# Experimentalist position, round 044 (conference)

No new runs. Read: notebook, STATE.md, DIGEST entries 040-043, STEERING.md.

## Most promising question (next 2-3 rounds)

Do gate-admitted J != 0 steps that keep the LNA key exist at n = 8 and n = 9 in classes 1 and 2, as E-145 found at n = 7, and does the x^2 law for Q(x) (E-141/E-143) hold on them?

Why: E-140/E-138 ("none keeps the key") is now known to hold for class 0 only. E-146 shows the E-143 shape is absent at n = 8 c0, so the law's support is thinner than the notebook's earlier claims. The key-guard-off run at n = 8 c1/c2 (depth 6, overnight candidate) is the cheapest way to say yes or no, and it is the case the H-015 question turns on.

## Weakest claim the workshop relies on

E-145's class-1 steps (13 of 67 distinct) rest on one script; the referee rebuilt only class 2 by hand. Downstream, "the key guard is no evidence for H-015 off J = 0 steps" and the class-0-only law are both built on that count. Second: E-142's "16 of 16 hits share 0 keys" is a bounded miss (no BFS closed), and it is still cited as support.

## What I need

- skeptic: hand rebuild of the class-1 E-145 steps (the 13), not from the script.
- theorist: a statement of when the orbit relation e_i = F^s e_w holds (E-145 says it fails on most H1/H2 steps in c1, c2); without it the Q(x) law has no test set beyond class 0.
- toolsmith: `--plan` size for guard-off n = 8 c1/c2 at depth 6, with checkpointing and a positive control (a random non-LNA parent keeps the key in 55 of 1264, E-140); save the steps as they are found.
- scholar: nothing this round; the PDFs stay parked unless supplied.

## Process notes

- Caps are not verdicts: report depth and cap with every "none found".
- Empty cells are vacuous; report non-vacuous cells per class.
- Before any run, check `research/EXPERIMENTS.md` for the n = 8 guard-off case.

## Suggested question

None. The agenda's item 1 already asks this; a second question would dilute it.
