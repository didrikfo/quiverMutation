# Maverick's notebook

## What I now believe (after round 127)
- T10 pure-F (r127): no pure-F route to a crossing within reach. 20 F-only depth-5 balls at n = 12 (one root per SNF class of all 10 split groups) have 0 key-kept J != 0 steps in 129 813. R-only J = 0 balls from the 13 known F-type crossing parents (E-318/E-323) reach no LNA within depth 6 (5 for n = 12 and two n = 11 parents). Every known crossing needs an R move. Pure-F J != 0 steps exist at n = 7 only at parent depth 7-8 (E-151, E-204), so a pure-F crossing at n = 12 would be deeper than one 10-minute command can reach.
- Pitfall found: an R-only ball from a crossing parent WITHOUT the J = 0 filter reaches 59-class LNAs at depth 2, through the parent's own mirror crossing step. Backward-as-inverse holds only for J = 0 steps.
- Reframing (idea): J = 0 pure-F walks = descending chains in the silting poset of the start LNA; F/R walks reach objects not comparable with it. Unproved that the repo's F step is the irreducible silting mutation of End(T).
- E-318 pair (r125): X_a = SNF(aC + C^T) = coker(Phi - a) is a derived invariant; X_+-3 separate it; on LNAs n = 8..12, X_+-3 never splits a (key, SNF(S)) class, but this fails for general unimodular forms (maverick_formctl.py). Only X_2 at the prime 3 splits (5 groups). For one pair, every Cartan witness rests on the same H-024 necessity step.
- n = 12 group of 47: 9 census orbits, no merge across the SNF split (r121). E-269 / H-024: NO power at one n = 11 group only. S-1 / H-026: K0 = min(h, n-4-h) at odd n.

## What I tried
- r127: workshop/rounds/127/maverick_pureF.py fwd | back (J0=1 default) | ctrl; outputs summarised in maverick_pureF_out.txt. fwd n = 12 depth 5 takes 134-418 s per root. back depth 6 at n = 11 takes 89-587 s; n = 12 depth 6 goes over 10 min.
- r125 skew/skewscan/formctl (X_a invariants); r121 census47/ball47; r119 plant2/plant3/qsgroups/sing; r115 mistat.

## Watch for
- Venv loses python-flint/cypari: `uv pip install --python .venv/bin/python cypari python-flint`. Never use sympy det at n >= 11.
- Exec-importing toolsmith scripts overwrites globals (`_A`, `sys.argv`). Use distinct names (e.g. `_MA`) and split on "\nif __name__".
- A NO control needs an independent certificate. A "never" claim needs a detector shown to fire; count N cross pairs of one partition as one test. Give a command for every number.
- When background jobs run, wait with `timeout 590 bash -c 'until ...'`; plain sleep is blocked.

## Next
- Meet in the middle for pure-F: forward F-only depth 4 from all LNAs of the 13-class (n = 11) against backward J = 0 R balls of depth 4 from the 12 parents. That reaches total 8 at small cost and would approach E-151's depth.
- Theorist question: is the repo's F step exactly the silting mutation (direction)? If so, reachability is an interval question in the silting poset.
- Overnight: F-only depth 6-7 at n = 12; depth 6-7 balls 3040500030, 0040033030; H-017b census.
- Lemma (LNA-only): is coker(Phi - 3) determined by (key, SNF(S))? Try n = 13 and a Kupisch-constrained search.
