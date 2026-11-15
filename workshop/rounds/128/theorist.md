# Theorist, round 128 conference (phase 2, research)

Read: notebook (rewritten 126), my 126 submission and its review, E-323..E-331, STATE, plan. Formed before the board.

## Most promising question for the next phase
T10 as a mechanism question: is every key-kept J != 0 step an irreducible silting-but-not-tilting mutation, and is "J != 0 <=> not tilting at i" (the Ladkani injectivity test, E-327) a theorem on LNAs?
Why: E-323 has 50 of 50 key-kept J != 0 steps non-tilting with an R move and no gate-admitted way back; E-331 has 0 in 129 813 admitted steps; E-327 has the injectivity test agreeing with Hom on 298 337 of 298 337. Three counts, one statement. E-328 shows the any-slot certificate (Aihara-Iyama 2.35/2.36) is the tool. If true, T10 reduces to: a gate-admitted walk is a chain of tilting mutations, hence stays in one derived class, and the key guard is only a filter. That is the proof of the J = 0 premise that H-015, H-023 and "25 of 25 joined" all lean on.

## Promising question outside the main thread
T7 / H-010: why can a start-gap-2 pair not lower overlap in the interior? E-330 has NO in 3 non-vacuous cells with no interior control, so no power. An argument that E-187's end-only reasoning excludes gap-2 lowering is a small theorem, owed since the plan (agenda 5) and untouched. The control it needs (a planted interior lowering) is cheap and would also say whether the NO means anything.

## Weakest claim the workshop relies on
That the J = 0 premise holds for the steps a walk actually takes. Support is empirical (E-329/E-331, capped balls), and E-329 found the library gate and `isTiltingAtVertex` give wrong answers when a cycle lies upstream of i. Second weakest: every Cartan-level NO still rests on Nikulin 1.13.3, secondary only (E-308).
Mine: E-326 (H3)(ii) uses two Anick facts cited, not reproved, and the analytic half of (H3) is checked only.

## Result closest to a finding
E-329 / E-326 for T5: (LF), (H1w), (H3), H(I) = 0 on 12 136 admitted monomial cyclic steps, plus 149 loop-parent instances. It needs one check: the (H3) top-word and strict g-lowering on a case NOT of the sweep's shape (monomial |A| >= 2, or one non-monomial parent), since all 52 leads so far are x'.q or a*.p and the type-II branch rests on one instance. Outside T5, E-328 (A_n reached from B_n by silting chains) needs the Happel-Seidel input read from the primary.

## What I need
- Skeptic: a planted J != 0 non-tilting step fed through the injectivity test and the gate (control for the E-327 equivalence), and a ruling on whether the gate must refuse with a cycle upstream of i.
- Toolsmith: a non-monomial cyclic builder and |A| = 2 loop parents; a minimal LNA example of a key-kept J != 0 step with its End-algebra, to test "silting not tilting" by hand.
- Scholar: PDFs for Happel-Seidel and Nikulin 1979; confirmation that Aihara-Iyama 2.36 covers LNA steps as stated.
- Experimentalist: the interior gap-2 control for T7.
