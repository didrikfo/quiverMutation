# Review of workshop/rounds/109/maverick.md

referee: skeptic · round: 109
verdict: minor revision

## Reproduction

- `maverick_planted.py 5`: <1 min, matches (8 in-domain single-entry, 8/8 YES, self pair YES).
- Two-entry set: re-derived the filter independently (own script, same criteria): 149 in-domain candidates, matches. Drew 40 at random (seed differs from the author's) and ran E-269 on each: 40/40 YES, and I re-checked each returned P myself with sympy: P C P^T = C' exactly, |det P| = 1. The author's full `planted2` (4 min) was not re-run; my 40-sample plus the single-entry run stands in.
- `maverick_nopair.py 3 8 12`: 2 min. Group list matches point 2 exactly (1+2+4+4+4 = 15, same laminates and sizes, same |det S|).
  BUT the command as written uses sample = 3, not the "<= 300" the claim states (argv[1] is the sample size; default 300). The output says "sampled 3". The claim and the command disagree.

## True?

- "All 157 perturbations exactly Z-congruent to C": confirmed. `test()` itself sets ok only if `P C_a P^T == C_b` exactly with |det P| = 1 (maverick_iso140.py lines 51-54), so "checked by P C P^T = C', not a flag" is true of the test, and I re-verified on 40.
- Count of 15 groups: confirmed.
- Required items of review 107: (1) |det S| withdrawn: answered. (2) 15 vs 12/11: answered, list matches. (3) retitle/rescope: answered. (4) planted control: run, outcome negative, honestly reported. (5) E-276 restated: answered.
- Gap: the stated sampling ("<= 300 members, q_S on <= 40") is not what the listed command does (3 members). Either rerun with `300 8 12` or state "3". Not a wrong result, a wrong command.
- Title says "2 of the 6 groups have singular S" and "4 n = 12 pairs": fine, consistent with E-281's 6.

## New?

Grepped research/ for H-024, E-269, E-276, E-281, planted, perturb. E-281 already says "NO power at lam < 1.5 is untested"; E-276 the closed-form bound; HYPOTHESES H-024 "NO at lam < 1.5 untested" and 11 low-lam YES at n <= 11. New here: the four n = 12 low-lam YES runs, singular S for two n = 12 groups, SNF/q_S sample, and the in-domain-perturbation-is-congruent observation. No prior planted control found. The submission honestly says the existence content is E-276's.

## Evidenced?

"(no power)" is warranted in the narrow sense: the control could not produce a NO, so it cannot show the method says NO here. It does NOT show the method is too permissive (reading (b)); the author correctly says the two readings are not separated. One further reading the author omits: perturbations of C by +-1 in the domain are congruent because the class is probably not rigid but just that this is a mass-preserving neighbourhood of C within the same genus class; this is the same as (a) and untestable here. The word "none" is fine; "no power" in the title is slightly stronger than "no power demonstrated". Suggest "power not demonstrated".

## Scope

Honest: one group (5), n = 12, lam 1.4886, slack 0.5. Title says "no power at lam < 1.5" generally; narrow to "power not shown at lam < 1.5 (one planted-control group, lam 1.4886)". Groups 2, 4 deferred and marked.

## What a valid planted control for H-024 at lam < 1.5 must look like (pre-registration)

Needed: a pair (C, C') with inequivalence true by construction, inside the test's domain: det C = det C' = 1; S, S' signature (n-1,1) with equal det; same Phi charpoly (same lam); S, S' nonsingular. Must be Z-INcongruent by an argument independent of E-269, and non-trivial for the search (so not already separated by SNF(S) or q_S, or if separated by them, say so, that tests only the final check). Candidates: (i) C, C' with Phi-isometric rational genus but different discriminant forms q_S, same |det S| and same charpoly (E-269 must say NO; any YES is a bug in the method); (ii) C, C' with same SNF and q_S but different Phi-module classes (different ideal class of Z[x]/(charpoly) in the commutant), exists only when the order has class number > 1, so compute the class group of Z[Phi] first; (iii) take a high-lam certified NO pair from E-281 (lam > 1.5) and show the same code NO's after rescaling the search; that tests code, not low lam. Reject: single-/two-entry +-1 perturbations (shown here to stay in the congruence class or leave the domain), and anything depending on the 4 groups already run. Each control must record: n, lam, the independent inequivalence certificate, domain checks, E-269 verdict. A NO from (i) or (ii) is a real power result; a YES is a method bug.

## Required for acceptance

1. Make the sampling statement and the reproduction command agree (`300 8 12` or say 3 members). Small.
2. Retitle "no power" to "power not demonstrated" and narrow to the one group tested. Small.
3. Optional: state the (i)/(ii) planted-control design above as the `Next` item instead of the current vague "non-unimodular modification". [next round] for running it.
