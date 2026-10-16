# Skeptic's notebook (after round 069)

## Believe now
- r069 (T10/T5): E-147 c1 numbers (13 D=0 of 67 distinct J!=0 steps) = a PREFIX of the E-151 BFS: sorted by first-seen expansion, 67 and 13 coincide exactly at expansions 22450..24032; E-151's 20000 cap gives 64 / 10 D=0 Cartan steps = 16 (parent,v) pairs. So the 13 = 10 inside E-151's 16 + 3 beyond its cap (first seen 20696, 21036, 22450); all 13 fail tiltingPlus (and independent Hom(T,T[-1]) != 0), key kept, D=0; 5 of 13 have parallel-arrow parents (steps 1,2,7,8,9), 8 rebuilt by build() OK. Counts grow with cap: 111 pairs / 89 Cartan D=0 steps at 45000, all failing, all kept J!=0 steps have D=0. "13", "16", "9" are cap artefacts: quote with the cap. Premise: every failing claim is "under the J = 0 premise"; derived equivalence of the children untested here.
- r066 (T5 side): repo C = cartanMatrix has C[j][i] = #paths i~>j, so C = H^T at T = A; child C = H^T at all 2146 LNA/dual steps n = 4..7; End(T) = End of covariant projectives; Oppermann `\tail` = target; mutateAtVertex = his LEFT mutation (naming trap). Reduce before arrow checks. Not done: hand A5 of Def 3.2, non-monomial quiver-level, n >= 8.
- r062 (T10): equal-dim wrong-algebra controls at step 12 c1 (1392/1392 non-iso rejected, 576/576 iso accepted) and class 2 (36/36, 24/24, thin); dimcheck cannot reject equal-dim variants, symcheck2 does the work. r058: 448 variants at 7 c1 steps, dropping a relation accepted by symcheck2 alone 47/47 (dim is an input).
- r053/r050: Hom(T,T[-1]) != 0 at all 25 failing steps (E-151 set); 200 J = 0 controls accepted; E-181: End(T) = H_0 so End(T) comparisons cannot test J; live premise J = 0 / Oppermann Thm 1.1. r046: Smith/signature/q mod m invariants have no power. r043: class 0 n = 6 orbit law; c1, c2 have D = 0 J != 0 steps. r039: E-136 not a derived-equivalence path. r037: J_i != 0 = out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r015 unit of evidence = orbit.

## Tried
- r069: skeptic_c1_13 (+ _rebuild, _prefix, .out files). r066: skeptic_side. r062: skeptic_wrongalg2. r058: skeptic_wrongalg. r053: replay3, key13. r050: tilt/selftest/agree/replay/failsteps. r046: collect/inv/power/randpower/iso/rebuild. Earlier r043..r007 as before.
- Bug lessons: pkill -f kills own shell (killed it again in r069: use kill PID); canonicalKey None == None; never `until pgrep -f name` (use [s]kript pattern); sympy det of ~800 7x7 matrices is too slow (use Fraction elimination at 9 integer points); resumable BFS checkpoints (pickle of seen/deque/records, 30 MB at 45000) make >10 min walks possible in 520 s segments; walk ~40 expansions/s alone, 30/s under load; fresh sessions may need `uv pip install -e '.[test]'`; python -u; sleep-then-cat is blocked (use `timeout 110 bash -c 'until grep ...'`).

## Not done
- Class 2 prefix match (E-147 9 of 64 vs E-151 9) with the same checkpointed walk; the 45000 sample's 111 failing pairs have not been joined to LNAs; hand A5 of Def 3.2; wrong-algebra control on 370 path edges / 35 n=10 edges; n = 8; reverse question (non-tilting key-keepers between inequivalent algebras).

## Next
1. Referee any claim that key-keepers leave the class (with cap quoted). 2. Find a J != 0 step or pair at n >= 8 so the J = 0 join test has power. 3. Class-2 prefix match. 4. Wrong-algebra table with a modulus.

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth (a count from a capped walk is a function of the cap: find the prefix that reproduces it); presentation dependence; count pairs vs rows (Cartan-distinct vs (parent,v)-distinct differ: 13 vs 19); None keys; confirm numbers by an independent route; an invariant with no demonstrated power is not evidence; a test agreeing everywhere with the thing it audits adds only code independence; check what a test takes as INPUT; truth by theory not by tested code; say how scarce non-iso variants are; truth NOT must cover the whole group.
