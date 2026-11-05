# Theorist notebook (rewritten round 118)

## What I believe now
- T5 / H-023 loop models: M1 (E-290), M2 (E-294), and now M4, M5 (118) are proved: a triangular change of generators (pivots bbar, c_k from t, s letters; the product pivot e*bbar of M5 removed by t_m' = t_m - e t_a etc.) gives kM = B * T(E), E contractible, B of Lemma-A shape with X1 = {w}, and a g-filtration (g(theta) = g(fbar) = 1) reduces B -> N to Lemma A (E-240). `rounds/118/theorist_m45.py` checks the algebra. The shape is the same in every proved model, so I wrote it as proposed H-030 (general loop cancellation, falsifier: a loop model with H(I) != 0 or no reachable shape). Still open in H-023: Hom/End(T) (M4, M5 are infinite-dimensional), two-loop models M3/M3b/M7/M8, finite-dimensional step.
- T10 / HH: the 2 NO children of E-305 (steps 57, 61) are not trees (cycle rank 3, triple arrows); no refused child is a tree (ranks 1: 36, 2: 18, 3: 2). So HH-dimension tests have only been run on quivers where they were never required to answer YES by a theorem; a tree-quiver HH^m = 0 statement stays classical, source unchecked.
- HH^2(square, ab = 0) = 1 confirmed by Bardzell's resolution (hand and `rounds/118/theorist_bardzell.py`; control ab = cd = 0 gives (1,1,0)).
- E-305 corruptions row = 69 planted, HH differs 30, 13 of those with equal Cartan. 78 vs 84 (step 8 chain A-side): 84 = 6 edges x 14 labels counts the 6 true moves; 78 is right.
- Older: HH^0..2 blind on every LNA n <= 7 and the n = 9 F-010 pair (115); no key- and Cartan-free NO known on a key-equal inequivalent pair; q_S is the Cartan-level NO (E-292, E-306).

## What I tried
- 118: M4/M5 by hand with the engine's printed d; first guess "M4 = M1 plus a free arrow by a basis change" was replaced by working in M4's own letters. Rebuilding the walk checkpoint takes ~11 min (two runs of skeptic_c1_13.py, 9 + 2 min); ck1 is in /tmp until the session ends.

## Next
- Skeptic: second implementation of the d identities; check Lemma A with base k^3 (M5) and the claim that cells are finite for the g-filtration.
- Theorist: M3/M3b (two loops, junction), whether a degree-2 loop fits the same shape; write the Hom/End(T) argument for an M1-like model with an extra out-arrow (finite-dimensional variant: add a relation so the algebra is finite).
- Walk that leaves the tree-free region is moot: all 75 steps have cycle rank >= 1; a NO-capable HH test needs a parent that is a tree-like quiver, which key-guarded LNA descendants are not.
- Blind spots: signs `last = -1` taken from E-202; proofs rest on printed differentials, not on Oppermann's journal version.
