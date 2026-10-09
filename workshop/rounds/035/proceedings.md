# Round 035 -- proceedings

## Scholar (T5): Hom(N,N[-1]) and the cyclic case
Claim: Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum b y_b = 0}, so "silting-not-tilting iff some J_i != 0" holds for every A and loopless v; on a cyclic quiver only the count changes. Referee (skeptic): accept; re-ran both scripts, all six rows match, re-derived M1 by hand. Decision: accept. Promoted: E-128.

## Experimentalist (T5): d_i at J_i != 0 on walks
Claim: 285 of 285 rows with J_i != 0 have d_i = 2; the 28 rows with d_i >= 3 all have J_i = 0. Referee (theorist): minor revision; no error, but a 150 s re-run did not reach d >= 3; "d = 2 exactly" is forced by L1; needs distinct-algebra counts and caveats. Decision: accept with qualifications (retitled around the real finding; the caveats are written into the entry; counts per algebra and a reproducible cap remain open). Promoted: E-129.

## Toolsmith (T5): size the n = 7 BFS closure
Claim: it does not fit one command; no hit in 31 000 expansions. Referee (skeptic): minor revision; numbers reproduce exactly at 120 s. Decision: note; promoted E-130 with "not falling" softened to a plateau of 2.4-2.5, the 8 h figure marked an extrapolation, memory and None-key drops marked unmeasured.

## Questions for the steering committee
1. **Agenda:** round 036 is a conference; the round-032 agenda stands until it. Recommend keep.
2. **Overnight:** the toolsmith proposes `toolsmith_closure.py run 1 --budget-hours 8`; it may not close and memory is unmeasured. Recommend no; the reverse-direction search is sized first.

## Decisions taken for the steering committee
- round 034, question 1 (agenda): keep the round-032 agenda -- decided by the chair of round 035; no answer from the human.
- round 034, question 2 (overnight): none; toolsmith sized the closure (E-130) -- decided by the chair of round 035; no answer from the human.
