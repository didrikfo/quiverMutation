# Skeptic's notebook (after round 058)

## Believe now
- r058 (T10): equal-dims wrong-algebra control for E-174's `symcheck2`: 7 c1 steps with one parallel pair and 3 line-type relations, 64 variants each (killed line re-chosen among 4): 448 variants all equal Cartan and dim K Q/I' (own linear algebra); 280 non-iso all rejected, 168 iso all accepted. So the test has power against relation errors of this kind. But it takes dim as input: dropping a relation (47/47) gives dim K Q/I' > dim End(T) and "iso" anyway. All 280 NO verdicts are genuine "NO label-preserving iso" strings (no timeouts); the 47 drops span 8 steps incl. 12. 25/25 + 370/370 show End(T) = child as given (crels complete assumed), nothing about J.
- r053: replay of E-160 paths with Hom_K(T,T[m]) accepted; every failing child has a Hom-tested J = 0 path. Test is Ladkani restated. r050: Hom(T,T[-1]) != 0 at all 25 failing steps (= sum J_i). 200 J = 0 controls accepted.
- r046: invariants (Smith, signature, q mod m) have no power child vs parent. r043: class 0 n = 6 orbit law; c1, c2 have D = 0 J != 0 steps. r039: E-136 meeting not a derived-equivalence path. r037: J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r015 unit of evidence = orbit.

## Tried
- r058: skeptic_wrongalg (+ out). r053: replay3, key13. r050: tilt/selftest/agree/replay/failsteps. r051: replay13, k1415, sensitivity. r046: collect/inv/power/randpower/iso/rebuild. Earlier: r043, r039, r037, r034, r031, r029, r026, r025, r023; nulls r021..r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell (twice!); never `until pgrep -f name` (matches itself); private /tmp dir (/tmp/tsm; pickles vanish between rounds, rebuild c1 ~6 min with DEADLINE 520 then rerun to resume); right-module orientation is not tilting; python -u; sleep-then-cat blocked; 4 cores; 3 points on P^1 have no moduli (no continuous-parameter test).

## Not done
- Step 12 (two parallel pairs) and class 2 for the wrong-algebra table; a wrong algebra produced by the search itself rather than perturbed crels; dim assertion in compare2; 13 class-1 E-147 hand-rebuild; reverse question (non-tilting key-keeping steps between inequivalent algebras); Hom(T,T[-1]) on n = 10 E-169 edges; n = 8.

## Next
1. Referee any claim that key-keepers leave the class. 2. Hom(T,T[-1]) on the n = 10 edges. 3. Wrong-algebra table with a modulus (4 points on P^1).

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test; an invariant with no demonstrated power is not evidence; a test that agrees everywhere with the thing it audits adds only code independence; check what a test takes as INPUT (dims, relation sets) before trusting its output.
