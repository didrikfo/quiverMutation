# H-020's rule table rests on six hypotheses; length independence of floating rules holds for w <= 5 to length 11 (0 failures in 10 654); "outside the derived class" for the E-154 children is defined, with the inside witness checkable and the outside witness not yet available

author: theorist · round: 049 · kind: proposal (statements) + small checks
thread: T4, T10 (i) · bears on: H-020, F-051, F-053, H-010, E-151, E-154, E-149
scope: statements are for the floating/anchored rule table (`lnaMoves.VERIFIED_MOVES`, 414 floating rules) and H-020's single-cluster slides; checks: floating rules of window width <= 5 re-verified at lengths w+5..11 (width <= 4: w+5..9 too); orbit sizes of core `45` at n = 13, 14 (reduced walk); the 9 failing children of the n = 7 class-2 E-151 walk rebuilt (the 16 class-1 children were NOT rebuilt). Nothing here is proved.

## Response to referee

1. Pairing and orbit count: corrected. From the saved output (`theorist_orbits45.txt`), equal orbit sizes pair o <-> n-8-o: n = 13, offsets 0..5 pair 0/5 (2386), 1/4 (1127), 2/3 (4217), and offset 6 = n-7 (447) is alone; n = 14, 0/6, 1/5, 2/4, offset 3 pairs with itself (2179), offset 7 = n-7 (272) is alone. "n-7-o" was wrong. So n = 13 has four distinct sizes, not three. Offset 0 is the head placement, not interior; offset n-7 is the tail placement. Among the interior offsets 1..n-8 (1..5 at n = 13) there are three sizes (1127, 4217, 2386), and the size 2386 recurs at the head offset 0, so the size symmetry runs across the head and is not an interior reflection. Text fixed in the Claim and H4.
2. Title: "two ... fail to follow" dropped. Only H4 is unexplained by the table; H3 is incomplete (F-052), a different statement. New title uses the referee's wording.
3. F-053 credit: point (a) was already recorded in F-053 and the F-051 amendment. The new content is the orbit sizes by offset only (n = 13, 14). Claim reworded.
4. Claim paragraph now states H1 was tested for w <= 5 only (30 of 414 rules; 330 rules of w = 6..8 untested beyond w+4); also H2, H5, H6 were not run (H6 is an unrun ablation, H2 is part of H1).
5. Reverse check: re-run (`theorist_reverse.py` on the rebuilt 405 s pickle) with `J_is_zero` printed explicitly; 9 of 9 `J {}`, `J_is_zero True`. Demoted to a note: the step P -> B has J != 0, so a J = 0 reverse step not returning P is expected. An informative outcome would have been a J = 0 reverse step landing on P or its opposite (B then inside the parent's class); that did not happen, and it says nothing about outsideness.
6. Outside witness: now called a requirements list, not a witness. The power-control demand is E-154's ("needs an invariant fine enough, or a tilting path back", "power ... untested"); the new part is only the witness form and the candidate list.

Not done: H1 for w = 6..8 at lengths >= w+5 (needs lengths 12-13, over 10 minutes), H6 ablation, class-1 children rebuild (previously flagged).

## Claim

H-020 ("placeability depends on the distance to the two ends") is a statement about the *move set*, so it inherits every hypothesis about the table. I list them as H1-H6 below, each with the check that could fail and the smallest case. Two points: (a) F-053 already records that F-051's "the interior is one orbit" is false; the only addition here is the orbit sizes: at n = 13 core `45` has orbit sizes by offset 0..6 of 2386, 1127, 4217, 4217, 1127, 2386, 447 (equal sizes pair o <-> n-8-o; offset 0 is the head, offset 6 the tail), so equal verdicts across the interior pairs (H4) are not explained by the table; (b) "outside the derived class" for a child B of the E-151/E-147/E-154 steps is defined in section 2, with a positive witness (a tilting path back, checkable) and a requirements list for a negative one (the power-control demand is E-154's; no outside witness exists). H1 was re-verified for floating rules of width w <= 5 only (30 of 414 rules; the 330 rules of w = 6..8 are untested beyond w+4): 0 failures in 11 170 new applications. H2, H5, H6 were not run.

## 1. Hypotheses behind the rule table (each: statement / check that can fail / smallest failing case)

**H1 length independence.** A floating rule `(w, before, after, seq)` verified at lengths w+1..w+4 (tests) holds at *every* length >= w+1 and wherever `matchesAt` fires (window exactly the pattern, no relation straddling the window, anything else outside). Check: `verifyMove` at larger lengths. Fails if some rule gives a wrong result, an illegal step or a moved Coxeter polynomial at length >= w+5; smallest: w = 3 at length 8 (two rules), w = 4 at 9 (four). **Done: 0 failures in 516 (w <= 4, lengths w+5..9) and 10 654 (w <= 5, to 11) applications.** Not done: w = 6..8 (the 330 rules that are most of the table) beyond w+4; w = 9..11 beyond w+2. One would need length 12-13 (about 2e5 LNAs): `OVERNIGHT.md` proposal below.

**H2 spectator independence.** The context outside the window never matters beyond "no straddling relation". It is part of H1 (verification runs over all LNAs), so it is tested exactly as far as length permits: a context needs w+5 or more positions before its first new pattern. Falsified for *non-window* effects in the sense of H-010 only if a rule's *sequence* is legal for one outside context and illegal for another; `verifyMove` records that as `illegal mutation`. None seen.

**H3 completeness vs soundness.** The table is *sound* (each rule verified, by gate, Coxeter) but H-020's "outside" needs *completeness*: outside means "no almost-separate row in the closed orbit of the moves (table + free/reduced + edges + doubles)", not "not derived equivalent to one". Check: enlarge the move set (a mutation search to depth d from an outside row) and see whether any outside placement acquires an almost-separate row. Fails at the first such row. This has already failed once: F-052, 60 placements at n = 11, 12 moved outside -> inside when the plain walk became the reduced one. Smallest case: those (n = 11, 12); the open question is whether the reduced walk is complete at n >= 13.

**H4 one verdict across interior pairs.** For a single-cluster core c and n >= 13 every interior offset has the verdict `o` (H-020 amended). The table is translation invariant (`windowStartsFor`), but *orbits of c@o and c@(o+1) are not shifts of each other*: n = 13, core `45`, reduced walk, orbit sizes by offset 0..6 are 2386, 1127, 4217, 4217, 1127, 2386, 447 (n = 14: 3767, 1636, 11820, 2179, 11820, 1636, 3767, 272; equal sizes pair o <-> n-8-o, offset 3 self-paired, offset 7 alone). So equal verdicts are data, not a consequence of H1. Check: a single-cluster core with an inside offset between two outside offsets at n >= 13 refutes it; E-051 found none in 1186 of 1192 comparisons, the six failures at n = 13 with slides of 1-2 offsets. Smallest failing case is at n = 13 (the six); a counterexample at n >= 18 would be the first beyond checked range.

**H5 stabilisation shape.** The slide at length n >= N(c) is `P(c) v^m S(c)` with fixed words P, S and v = o (amended statement). Falsifier: a core whose P or S changes between n and n+1 for n >= 13 (E-051: 0 of 1186 among complete pairs). N(c) is not known as a function of c; "three offsets or more" is an empirical reading. A falsifiable form: N(c) <= width(c) + 4 where width = arrows spanned; test on the six n = 13 failures (E-051: words of six or seven letters; spans not computed). Not tested here.

**H6 anchored rules see only the ends.** `windowStartsFor` puts an anchored rule at window 1 or n-w only, so only head and tail depend on anchored rules. Falsifier: a placement whose verdict changes when the anchored rules are removed (`overlaps.py --no-rules` style) while its core stays >= w from both ends. This is a code-level statement and should be checked as an ablation: smallest case n = 13, `45` at offsets 2, 3 with ANCHORED removed (not run).

Of these, H3 is the weakest (refuted once for the plain walk, F-052), H4 the least explained by the table.

## 2. What "outside the derived class" must mean for the 25 children

Setting. Parent P is derived equivalent to the LNAs of its key class K (n = 7, classes 1, 2) by a walk of forward steps all believed to be tilting; the step P -> B at vertex v is gate-admitted, keeps the Coxeter key, and has J != 0 (`tiltingPlus` false). **Assumption A1:** at n = 7 a key class is one derived class of LNAs (support: E-079 key-coarser pairs start at n = 12; it is a mutation-class fact for LNAs, not a proof of derived equivalence).

**Definition.** B is *inside* iff there is a tilting complex T over B with End(T) isomorphic to an LNA of K (equivalently D^b(B) ~ D^b(P)); *outside* iff there is none. The key being kept shows nothing: it is a necessary condition, and E-154 shows Cartan congruence (Smith form, signature, q mod m, an integral P) cannot separate even two LNAs of K.

**Inside witness (finite, checkable).** A sequence B = B0 -> B1 -> ... -> Bk with each step a gate-admitted mutation with J = 0 and `tiltingPlus` true (either direction, by the dual), and Bk (or its dual) an LNA whose key class is K. Verify by (1) independent rebuild of each Bi from arrows and relations, (2) Bk's canonical key equal to a listed LNA's. This is a proof of insideness (Ladkani 2.3(c) as used in `skeptic_back.py`). A *failed* search is not evidence of outside: E-149 loses 10.3% of reverse edges and the E-154 control missed too.

**Outside witness (a requirements list, not a witness; the power-control demand is E-154's).** A derived invariant I, computable from the presentation, such that I(B) differs from the common value I(A) on K. Requirements: (i) I is a derived invariant by a cited theorem; (ii) I takes one value on all 66 / 91 LNA pairs of K; (iii) **a power control**: some pair of algebras at n = 7 known to be derived inequivalent with equal Cartan data on which I differs (E-154: no known out-of-class equal-key pair at n = 7, so (iii) forces a smaller or non-LNA control, e.g. a hereditary algebra vs a tilted algebra of a different type with the same Coxeter polynomial, which Cartan data separates or not). Candidates: the finiteness of global dimension (derived invariant; all children have Cartan det +-1 so it is not excluded, but unmeasured); dim Z(B) = HH^0, dim HH^1, HH_1; whether B is piecewise hereditary of the type K's LNAs are (Happel: LNAs of K derive to a fixed hereditary or canonical type when Euler form is definite; F-048 is the indefinite criterion); for gentle B the Avella-Alaminos-Geiss invariant (the children carry commutativity relations with paths of length 2 and 3 and 5-6 arrows, so they are not gentle: it does not apply).

**What the children look like (class 2, 9 children, `theorist_children_c2.txt`).** All have 7 vertices and 7 arrows (one unoriented cycle, no parallel arrows), a commutativity relation (a path of length 2 or 3 equals another), plus 1-3 zero relations of length 2, dim 15 or 17 against 19 or 22 for the parent. So a child is a tree-with-one-square-like algebra, not an LNA, and not obviously monomial.

**Note: one-step reverse check (`theorist_reverse_c2.txt`, J printed as `J_is_zero`).** For each of the 9 children, mutating B at the same vertex v in the original direction is gate-closed in 9 of 9; in the opposite direction it is gate-open with J = 0 in 9 of 9, but the result is neither the parent nor its opposite (canonical key) in 9 of 9. So each child has a genuine tilting neighbour that is not P. This is expected, not informative: the step P -> B has J != 0, so it need not be inverted by a J = 0 step. An informative outcome would have been a J = 0 reverse step landing on P or its opposite (B then in the parent's class); none occurred in 9 of 9, which proves nothing about outsideness.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/049/theorist_rulelen.py 9 4      # 4 s
timeout 10m .venv/bin/python workshop/rounds/049/theorist_rulelen.py 11 5     # 183 s
timeout 10m .venv/bin/python workshop/rounds/049/theorist_orbits45.py 13      # seconds (n = 14 likewise)
timeout 10m .venv/bin/python workshop/rounds/046/skeptic_collect.py 7 2 20000 /tmp/c2.pkl 20   # 405 s
.venv/bin/python workshop/rounds/049/theorist_children.py /tmp/c2.pkl         # seconds
timeout 10m .venv/bin/python workshop/rounds/049/theorist_reverse.py /tmp/c2.pkl   # seconds
```

## Prior record

H-020 / F-051 / F-053 (orbits per reflected pair), F-052 (plain vs reduced walk moved 60), E-046/E-051 (975 and 1186 comparisons), E-154 (Cartan invariants non-informative; asks for a finer invariant or a tilting path), E-149 (reverse loss), E-151 (16 + 9 failing steps). `verifyMove` already ran w+1..w+4; the longer lengths and the orbit sizes by offset are new data, the list of hypotheses and the definition are new statements. I did not find H1's extension or the 45-by-offset orbit sizes in `research/`.

## Code changed

None under `quivermutation/`. New scripts only: `workshop/rounds/049/theorist_{rulelen,orbits45,children,reverse}.py` and outputs.

## Next

- toolsmith: H1 for w = 6..8 at lengths w+5..12 is an `OVERNIGHT.md` item (about 2e5 LNAs at 13; one rule-batch per core); H6 ablation (remove ANCHORED_MOVES) at n = 13, 14 on `45`.
- skeptic: for T10 (i) fix the power control first (an out-of-class pair at n = 7 or a non-LNA control with equal Cartan data) and only then compute global dimension, dim Z, HH^1 on the 25 children; rebuild the 16 class-1 children for my reverse check and follow each tilting neighbour (the reverse step in section 2) with `skeptic_back.py`, which would test insideness with a start one step nearer.
- theorist (me): a shift-correspondence between orbits of c@o and c@(o+1) (H4); why the pairs differ in size but not verdict.
