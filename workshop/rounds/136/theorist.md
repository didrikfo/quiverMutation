# Theorist, conference 136 (phase 4 close, looking at phase 5)
Read: notebook, my 134 submission and review (none from me in 133, 135), E-341..E-349. Formed before the board.

## Most promising question (phase 5)
H-023 / T5: is the library child isomorphic to End(T) when the CHILD is cyclic, and what does the library refuse at
the 22-27 gate-admitted steps where `child_of` returns None? E-345 (587 / 587 ISO, n <= 6 acyclic) and E-348 (341
planted same-Cartan wrong children all NO) make the acyclic case a computation. H-023 as stated is not restricted to
acyclic children, and E-193 says the library's Cartan matrix is wrong on cyclic children (60 unchecked): the open part
is where the claim is likeliest to fail. My End(T) engine takes any acyclic parent; only a cyclic-safe child engine is
missing. ISO extends the statement; a counterexample gives it a precise hypothesis (E-340 already hints at one).

## Promising question outside the phase's main thread
T7 / H-010: why do clean orbits avoid every overlap-lowering row? E-343: lowering rows contain a gap-1 (3,3) pair in a
run of exactly three overlap->=2-linked relations (F-022); (8,4)(10,3) cells never reach one. Is "clean" equivalent
to an invariant of `VERIFIED_MOVES` (no (3,3) pair reachable)? The end-only gap-2 argument is still owed by me. A proof
would turn E-330 / E-343 from harness nulls with no clean YES control into a theorem.

## Weakest claim
E-344 / H-032: "refusing non-tilting steps costs 0 in-class nodes". It rests on 3 key-kept refusals (2 algebras up to
iso), all at step 6, depth 6 + 3-4, on 2 of 7 E-332 roots, with no in-class loss control; the zero is an upper bound
only by the author's accounting. Second: every NO in E-345 / E-348 is a sympy Groebner verdict (one hand-checked,
one rank-argument NO).

## Closest to a finding
E-347 (H-025 / T3): kA_15/rad^3 etc. are not fractionally CY; two independent nu codes agree on 15 700 pairs with 0
disagreements (GF(2), GF(2^31-1)), and 160a = 144b = 162a needs no Fuchsian input. One check: reproduce the shifts at
2m = 288 on both codes (expect 320 / 324 if genuine orbit data, not a window artefact) and re-run the char-0
base-change step the author did not re-check. Then CONFIRMED. (Second: E-345, which needs the cyclic-child run first.)

## What I need
- Toolsmith: a cyclic-safe child dimension / Cartan engine, or the `child_of` None reasons, to feed `theorist_endt.py`.
- Skeptic: a planted wrong child with same quiver and Cartan that is neither monomial nor KILL/SCALE (E-186 modulus
  family at n <= 6), to test the iso step beyond its shown power.
- Maverick / scholar: the m = 288 run for E-347.
