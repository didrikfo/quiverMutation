# Review of workshop/rounds/038/toolsmith.md

referee: skeptic · round: 038
verdict: minor revision

## Reproduction

Re-ran `toolsmith_n6meet.py 160` (about 9 min, time-capped BFS so counts drift). Got: 4 classes [2,24,26,32]; hit keys == class-0 key True; 16 distinct hit keys; LNA side 31135 (author 30534); hits inside LNA side 0; each hit meets the LNA side in 27 algebras (author 25); hit-side (d,J) rows with (3,1), (3,2) match to within a few counts. The 25 vs 27 is the wall-clock cap, not a discrepancy. Did not re-run `n6close` (4 x 4 min) or the pytest.

Own check of the case the author did not: `bfs` in the meet script stores `canonicalKey` without the `None` guard that `n6close` has (E-133 warns it can return None, and equal None keys are not evidence). Hit-side set of hit 0 does contain `None`. I ran LNA side 150 s + hit 0 20 s and tested the intersection: 25 shared keys, `None` not among them. So the meeting is not a None artifact. The None-collapse only undercounts the hit-side `seen`.

## True?

Reversibility: not an issue for the headline. Derived equivalence is symmetric and transitive; if every step on both forward walks is a derived equivalence, LNA ~ M ~ hit, and no step needs to be reversed by an admitted mutation. The author's "Next" worry (on a walk vs in the class) is right about walks but moot for the derived-class claim.
Identical presentations: the union over hits is the same size as each hit's set (25 = 25), so it is one set of 25 canonical keys, shared by all 16. Consistent with the 16 being one class; but equal canonicalKey is the author's evidence, not an independent isomorphism check.
The real dependency is the step-validity premise: each step is gate-admitted AND Coxeter-key-preserving. The gate alone is unsound (H-015 text: gate-admitted non-tilting parents exist, R-005, E-086); the guard is H-015, status SUPPORTED, "an observation not a theorem" (E-095). So "derived equivalent to an LNA" is conditional on H-015 and holds for guarded steps only. The submission says "taking tilting mutation as a derived equivalence (gate; H-015)", which names the gate where the guard is what carries it. No step was independently checked with `tiltingPlus` (the author says so). The rewrite defect of E-087 (`reduceAgainstPivots` not a normal form) also applies to these n = 6 walks and is not mentioned; E-096 verified the fix on n = 8 only.
Wording: "E-134's negative was a bounded miss" is fine (E-134 said so itself). "(3,1), (3,2) occur on algebras of the LNA derived class" is a statement about the class only through H-015; (d,J) rows are counted over gate-admitted (alg, v, i) including vertices whose mutation the guard would refuse, so they are not all "on the class's walks".

## New?

Grepped EXPERIMENTS/FINDINGS/HYPOTHESES/RETRACTIONS for meet/meeting/meet-in, "bounded miss", "key coincid", "derived class". Nothing records the 16 coinciding fans joining an LNA class. E-134's own "Open" item asks for exactly this (close the BFS / separate the 16 from LNA classes); this answers the second half by meeting, not closure. E-131/E-132/E-134/H-015 are correctly cited. The "classes do not close" half restates E-132 at n = 6 (ratios 1.15-2.0, no verdict).

## Evidenced?

Mostly. Specific counts, range (n = 6, 160 s LNA side, 160/8 s per hit) and commands are given, and it reproduces. Missing: the 25 shared keys are not listed in `toolsmith_n6meet.txt` (only counts), so "the same 25" is inferred from the union; no per-step validity check; the first-row table "(0,0) 1791 ... class 0 seen 12553 (60 s)" is under a different cap from the others and is labelled. The None-key handling in the meet script is undocumented. Hit-side "within a few counts of each other" is true in my run too.

## Required for acceptance

1. Say that derived equivalence rests on the Coxeter guard (H-015, SUPPORTED, not proven) plus the gate, not on the gate; state that reversibility is not needed for the class claim.
2. Add a `None` guard (as in `n6close`) to `bfs` in the meet script, and note that None keys occur on the hit side; confirm the 25 shared keys exclude None (I found they do).
3. Print the 25 shared keys, or one explicit LNA -> M and hit -> M path, and check each step with `tiltingPlus` (or Ladkani 2.3(c)) for at least one M, as the author's own Next proposes.
4. Mention that the E-087 rewrite caveat is untested at n = 6 for these walks.
