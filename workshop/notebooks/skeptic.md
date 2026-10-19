# Skeptic's notebook (after round 074)

## Believe now
- r074 (T10): the 70 key-kept failing edges of the complete n = 7 balls (parent depth <= 8, c1 34, c2 36, nothing budget-stopped) have 60 distinct canonicalKeys, 12 distinct children up to vertex/arrow relabelling (7 c1, 5 c2; sign-aware check 0 failures). The E-151 20000-expansion window is reproduced by position: 16 c1 + 9 c2 = E-157's 16 / 9, 10 + 9 child keys; it covers 8 of the 12 classes; 4 classes (14 edges, parents at depth 8) are in no recorded join. canonicalKey over-counts by x5. Hand check (c2, v = 6, new class): Hom(T,T[-1]) = 1 at (3,6) = kernel of right mult by 6->2 on the two paths 3~>6, End(T) Cartan = child, Coxeter polys equal. Join of the 4 classes not searched.
- r071 (T4): own rebuild of verifyMove (skeptic_verify.py) gives w7 L10 130 rules / 650 applications / 0 failures = E-188; controls flag bad sequence 560/650; Coxeter test has no isolated power. Engine not independently rebuilt. T4 can be closed.
- r069: E-147 c1 numbers are a prefix of the E-151 BFS; 13, 16, 9 are cap artefacts: quote the cap. Every failing claim is "under the J = 0 premise".
- r066: repo C = cartanMatrix has C[j][i] = #paths i~>j; mutateAtVertex = Oppermann LEFT mutation. r062/r058: dimcheck cannot reject equal-dim variants, symcheck2 does. r053/r050: Hom(T,T[-1]) != 0 at the 25 failing steps (skeptic_tilt.py needs S = -1; S = 0 gives wrong End(T)); E-181: End(T) = H_0 cannot test J. r046: Smith/signature/q mod m invariants have no power. r043, r037, r034, r015 unit of evidence = orbit.

## Tried
- r074: skeptic_failkids (+checkpointed shards), _an (iso classes, window, sign-aware), _hand. r071: skeptic_verify. r069: skeptic_c1_13*. r066: skeptic_side. r062/r058: skeptic_wrongalg*. r053: replay3, key13. r050: tilt/selftest/agree/replay/failsteps. r046: collect/inv/power/randpower/iso/rebuild.
- Bug lessons: pkill/pgrep -f kills own shell (even the [x] trick when the pattern is in the same command line; use ps|grep|awk|xargs kill in its own call); checkpoint paths must not depend on the shard's output name; toolsmith_collect.py 7 1 20000 exceeds 10 min (class 1 slow, killed twice); canonicalKey None == None; sympy det of many 7x7 too slow; ~40 expansions/s; python -u; LNA row has n-2 entries.

## Not done
- J = 0 join of the 4 new classes and their depth-8 parents; are the 4 classes opposites of window classes (WL on op algebras); class 1 window via toolsmith_collect (only positional). T4 other cells, Coxeter-only sabotage, engine-independent mutation. 45000-sample 111 pairs not joined to LNAs; hand A5; wrong-algebra control on 370 path edges / 35 n=10 edges; n >= 8; non-tilting key-keepers between inequivalent algebras.

## Next
1. If toolsmith joins the 4 classes: referee with Hom test + isomorphism at meets. If a class does not join, try a NO control of demonstrated power (n = 9 F-010 5+5 is capped). 2. J != 0 step at n >= 8. 3. Wrong-algebra table with a modulus. 4. Close T4.

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth; pairs vs rows; labelled key vs isomorphism class (count both); None keys; confirm by an independent route (here: a positional reconstruction agreeing with E-157's counts); an invariant with no demonstrated power is not evidence; say which parts of a test are shared with the thing audited; check what a test takes as INPUT; truth by theory not tested code; say how scarce non-iso variants are.
