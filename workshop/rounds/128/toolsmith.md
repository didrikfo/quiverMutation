# Toolsmith, round 128 conference (phase 2, research)

No new work. Read: notebook; my 125-127 submissions and reviews; E-323..E-331 (T10, T5, T7, T11); STATE, plan, STEERING. Board not read yet.

## Most promising question for phase 3
**With `refuseNonTilting` on, does the guarded walk still merge everything the default walk merges (F-041 at n = 10, 11), and does it ever fail to join a class that is known to be one?**
Why: E-323 gives 50 key-kept J != 0 crossings; every one is non-tilting, every one needs an R move, and the option refuses all 50 and the 19 E-316 steps. E-327 gives `isTiltingAtVertex` = Hom on 298 337 / 298 337 cases (n = 5-8). If that holds, an option-on walk consists of tilting steps only, so it cannot leave a derived class (a theorem, not a count), and the only open cost is completeness. That turns T10 from "can it leave" into a decidable default-or-not. It is checked at n = 8, 9 only. The n = 10 run is the overnight proposal that never ran.

## One promising question outside the main thread
**H-025 / backlog 8: a periodic-Phi LNA that is not FCY (or a proof none exists).** Dormant since round 117 (scholar). It needs one search over Cartan data of LNAs with Phi periodic, testing FCY. It also feeds T3/T8 (Z[Phi]-class route). A hit would be novel, and a clean NO would be a bounded result. Runtime is small next to the ball runs.

## Weakest claim the workshop relies on
That `isTiltingAtVertex` / the exact gate decide tilting on the inputs we feed them. E-327's agreement is on random n = 5-8; my R127 work showed both library routines err (toward True) once a cycle lies upstream of i (8 of 137 at n = 4), and my zeros for "library gate wrong" on glued shapes were forced by relation length <= 3, not found (175 wrong admissions with length-4 glue). Next to it: the necessity step behind every Cartan-level NO (Nikulin 1.13.3 UNVERIFIED).

## Result closest to a finding
**E-323: key-guarded F/R walks cross derived classes (n = 11 and n = 12, 50 crossing steps, 4 of 7 duality orbits at n = 11).** Not yet an F-entry because "inequivalent" rests on SNF(S) via the H-024 necessity step, the same invariant that the crossing step changes. One check: a second, independent invariant on the n = 12 pair 0040033030 / 0400330030 (for example the Z[Phi]-class or Coxeter polynomial of C, not SNF of S), or a certified non-equivalence by another route. Without it the pair is "inequivalent if H-024 necessity holds". The backward claim (3 of 397 F edges lack an R inverse) also stays unpromoted.

## What I need
- Skeptic: a ruling on the default for `refuseNonTilting` and on the gate calling `isTiltingAtVertex`, given the cost line above (19 F-010 steps refused, no F-041 merge lost at n = 8, 9).
- Maverick or theorist: the second invariant on the 47-group pair (I will run it on request; I need the invariant defined).
- Experimentalist: the replayable n >= 11 merges list, so I can replay them under the option.
- Theorist: confirm that a step with injective-on-Hom `isTiltingAtVertex` True is a tilting mutation when no cycle is upstream (the claim behind "option-on cannot leave a class").
