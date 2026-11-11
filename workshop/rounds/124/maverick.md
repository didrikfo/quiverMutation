# Maverick, round 124 (phase 1 conference, no new work)

Read: notebook, my r121 submission and review, E-314..E-322, STATE, plan. I have no r122 or r123 submission; my r123 files are referee scripts (maverick_ind/k1/path.py). Form below from the E-entries before the board.

## Most promising question
Does the E-318 walk shape (a gate-admitted, key-kept J != 0 F/R step that moves SNF(S) tail) exist inside the n = 12 group of 47 (14 + 33 by SNF(S))? E-317 showed my balls were all J = 0 tilting (crossing impossible by theorem). E-318 shows crossings exist at n = 11 in the 72-group (13 + 59) but toolsmith's search stopped at n = 11. Run `toolsmith_jxsearch.py` (F and R, depth 4) from the 3 + 2 roots of the 47-group, with E-318's planted crossing as the power control. A hit gives the T10 YES at n = 12 and a walk across the H-032 split. A NO at depth 4 is only a bound, but it is the first real test of this group.

## One question outside the T10 thread
Is the 14 + 33 split (and the 13 + 59 split) confirmed by a second, independent invariant? Plan item 6 asks for it and nobody has done it. Candidates: Coxeter charpoly of Phi, discriminant form q_S of the nonsingular parts. Every Cartan-level NO now flows through one invariant (SNF(S)). Same dormant line: H-017b (per-LNA defect statistic, text written r119, never run).

## Weakest claim
"SNF(S) differs, so not derived equivalent" (E-308 chain) carries E-314, E-317 and E-318 and so H-015's refutation. Its sources are Ladkani 1001.4765 Lemma 3.5, Rickard, Aihara-Iyama, matched in text only. Nikulin 1.13.3 is still UNVERIFIED (secondary), and the chain assumes finite global dimension for each LNA. The step is plausible, but nothing here exercises it on a pair known to be derived equivalent with a non-J = 0 route. Runner-up: "H-015 refuted" rests on one hit (one group, one root).

## Closest to a finding
E-318 (n = 11 crossing, 5-step walk, SNF tail (2,4) -> (2,2)). It was refereed with hand-built Cartan matrices and mod 4 and mod 8 solution counts, which all test SNF(S) again. One check it still needs: an invariant independent of S. Compute the Coxeter charpoly (or det and q_S) for 003033030 and 030330030. If they differ, the inequivalence has two legs. If they are equal, the refutation stands on SNF alone and the weakest-claim item above becomes the load-bearing risk.

## What I need
- toolsmith: the jx search tooling (`toolsmith_jxlib.py`, group JSON) generalised to n = 12 and the 47-group roots, plus a runtime estimate (my guarded balls cost 45 s at depth 4 per side, so F/R with 2 move kinds may exceed 10 min; segmented run).
- skeptic: a pre-registered prediction for the n = 12 search (hit or no hit at depth 4) before I run it, and a second referee computing the independent invariant on the E-318 pair.
- scholar: whether Avella-Alaminos-Geiss or Ladkani give any complete-invariant statement for LNAs (not gentle in general); my r121 literature slot was never used.
