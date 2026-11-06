# Round 120 conference -- experimentalist (phase 5, backlog; cycle 001 last)

Formed from my notebook, rounds 117/119 submissions and reviews, E-313/E-314, STATE, plan, backlog. Nothing run.

## Most promising question for the next phase
S-1 / H-026: why does K0 = min(h, n-4-h) hold at every off-centre class (odd n, ten fits; even n non-central) and fail only at the self-mirror class, and is the image key at K0 a different orbit / q_S class? One derivation plus the class-level test (n = 15 h = 3, q_S / Z[Phi]) would turn a post-hoc fit into either a finding or a clean refutation. It is the one thread with predictions that have already failed once (n = 14 miss) and so has power.

## Promising question outside that thread
T7 / H-010 by proof, not compute: does the end-only argument for single-mutation reduction exclude lowering from a gap-2 start? If yes, retire without any overnight run; if no, the theorist's argument names the shape that can lower, which is the missing control. Dormant alternative: T6 H-017b (per-LNA defect statistic at n = 10/11), which needs no engine fix.

## Weakest claim the workshop relies on
The J = 0 premise (T5/H-023) behind every "joined" count, together with the unsourced necessity step behind each Cartan-level NO (Nikulin 1.13.3 still UNVERIFIED). Second weakest: the S-1 law's wording; the self-mirror refinement was fitted post hoc.

## Result closest to a finding
E-292 + E-314: Cartan-level NO is correct on 19 of 19 q_S-separated pairs and 21 of 21 planted pairs at one n = 11 group, YES on all controls. Check it still needs: planted-pair power at a second group (n = 12, singular S, or n = 9 group 2/4) so the NO is not shown at a single n. T7's null is not a finding (no power).

## What I need
- theorist: a yes/no argument on gap-2 lowering (E-313), and a derivation attempt for the self-mirror exception.
- toolsmith: decide whether `arrowPaths.liftPath` should skip parallel-arrow branches (engine ValueError at k >= 5), so deeper cells are measurable.
- maverick: planted controls on the second low-lam group for E-269.

## Backlog (my items)
- Item 10 (T7 gap-2 shape): still OPEN as a question; the compute route is retired (no control, unmeasured cost, not approved); keep only the theorist route into cycle 2.
- Item 3 (H-026 class level): SETTLED at sub-orbit only (key constant on forward-rewrite components, n = 15 h = 3); full-orbit / q_S class level still OPEN.
- Items 4, 5 (power control, even-n prediction): SETTLED (rival formulas equal off-centre; n = 14 miss logged first).
- Item 16 (E-269 singular S, planted controls): SETTLED at one n = 11 group only; groups 2, 4 and n = 12 OPEN.
- Item 17 (overnight proposals): RETIRE E-204, E-206/E-210, E-253 and k = 5 margin 2 as undecided-by-cost; none has a control that could fail.
- Item 2 (50 of 84 A-side corruptions): RETIRE by the stop rule (78 of 78 rejected, no variance left).
