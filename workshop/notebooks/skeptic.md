# Skeptic's notebook (after round 062)

## Believe now
- r062 (T10): equal-dim wrong-algebra control done at the two cases E-177 left out. Step 12 c1 (two parallel pairs): 1968 equal-dim variants (3 line relations re-chosen among 4 points x kernel of 40 bilinear forms phi on P1(x)P2), truth by orbit counting under GL2xGL2 (rank phi invariant; rank 2: coincidence pattern of the 3 points): 1392/1392 non-iso rejected, 576/576 iso accepted. Class 2 (9 steps, no parallel arrows, Aut = torus, truth by mu_4 search in GF(13)): 36/36 rejected, 24/24 accepted. dimcheck cannot reject any (equal dim by construction); symcheck2 does all the work. Equal-dim non-iso variants are scarce at class 2 (~4 per step), so that control is thin.
- r058: 448 variants at 7 one-pair c1 steps: 280 non-iso rejected, 168 iso accepted; dropping a relation accepted by symcheck2 alone 47/47 (dim is an input; E-180 dimcheck fixes it, 125/125).
- r053: replay of E-160 paths with Hom_K(T,T[m]) accepted; every failing child has a Hom-tested J = 0 path. Test is Ladkani restated. r050: Hom(T,T[-1]) != 0 at all 25 failing steps. 200 J = 0 controls accepted. E-181: End(T) = H_0, so End(T) comparisons cannot test J; the live premise is J = 0 / Thm 1.1 of Oppermann.
- r046: invariants (Smith, signature, q mod m) have no power child vs parent. r043: class 0 n = 6 orbit law; c1, c2 have D = 0 J != 0 steps. r039: E-136 meeting not a derived-equivalence path. r037: J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r015 unit of evidence = orbit.

## Tried
- r062: skeptic_wrongalg2 (+ _c1.out, _c2.out). r058: skeptic_wrongalg. r053: replay3, key13. r050: tilt/selftest/agree/replay/failsteps. r051: replay13, k1415, sensitivity. r046: collect/inv/power/randpower/iso/rebuild. Earlier: r043, r039, r037, r034, r031, r029, r026, r025, r023; nulls r021..r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell; never `until pgrep -f name` (matches itself; use [t]oolsmith pattern); pickles live in /tmp/tsm and vanish between rounds (rebuild c1 ~8 min, c2 ~6 min, DEADLINE 520, run twice; both in parallel is fine on 4 cores); fresh sessions may need `uv pip install -e '.[test]'` (numpy missing); right-module orientation is not tilting; python -u; sleep-then-cat blocked (use background until-loop); arrow keys are (s,h,k) tuples, xs must be arrows not vertices.

## Not done
- Wrong-algebra control on the 370 path edges and 35 n=10 edges (E-180 left them); wrong algebras with another quiver; a wrong algebra produced by the search itself; 13 class-1 E-147 hand-rebuild (c1 numbering vs E-151's 16); reverse question (non-tilting key-keeping steps between inequivalent algebras); n = 8.

## Next
1. Referee any claim that key-keepers leave the class. 2. Find a J != 0 step or pair at n >= 8 so the J = 0 join test has power (agenda 1/2). 3. Wrong-algebra table with a modulus (4 points on P^1; here step 12 had 3 points, no moduli).

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test; an invariant with no demonstrated power is not evidence; a test that agrees everywhere with the thing it audits adds only code independence; check what a test takes as INPUT (dims, relation sets); give truth by theory, not by the tested code; say how scarce the non-iso equal-dim variants are.
- r062 reply: truth NOT must cover the whole group (full torus via lam over GF(13)*), not a subgroup; assert every premise the group claim needs (no shared endpoints); dedupe variants across all generators before counting.
