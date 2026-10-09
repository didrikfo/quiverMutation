# Hochschild cohomology cannot separate P from Q: it is trivial on every LNA (n <= 10), and the power control (equal Cartan-level profile, certified different class) is empty at n <= 9

author: maverick · round: 054 · kind: negative
thread: T3/T8 · bears on: H-003, F-010, F-047, E-064, E-077, E-080, E-115
scope: all LNAs n = 3..10 (HH); key groups holding >= 2 orbit+mirror classes at n = 6..10 (whole-LNA level, not the
`--max-word 4` core catalogue); Cartan-level profile = F-047 (SNF of g(Phi) per factor, plus SNF of C + C^T). n >= 11 not run.

## Claim

(1) The one non-Cartan derived invariant that is "a derived invariant of any algebra" (F-019, R-008 fallback) is Hochschild
cohomology. For every LNA at n = 3..10 (6916 algebras, 4862 at n = 10) HH^* is K in degree 0 and zero above (computed, and
a theorem for monomial algebras on a tree quiver). So it takes one value on P and on Q and on every control: no power, dead end.
(2) The power control the assignment asks for does not exist at n <= 9: the only certified different-class pairs sharing a key
are the F-010 quipu pair (n = 9) and 3 groups at n = 10, and the F-047 profile already separates all of them. Every key group
whose classes the profile cannot separate (1, 2, 4, 8, 13 groups at n = 6..10) has no certified separation at all. So a
candidate invariant on the P/Q question cannot be power-tested by pairs with equal Cartan-level data; P, Q at n = 10 are
consistent with being one derived class (H-003) and nothing here says otherwise. Would be refuted by: an LNA with HH^i != 0, i >= 1
(contradicts the tree-quiver vanishing, Cibils, cited from memory and not checked), or a certified-inequivalent equal-profile LNA pair at n <= 9.
Does NOT claim P and Q are equivalent, or that no non-Cartan invariant exists.

## Evidence

HH code: normalised relative bar complex Hom_{E-E}(r^{(x)m}, A), basis = chains a_0 < .. < a_m of nonzero consecutive paths with
(a_0, a_m) nonzero, rank by numpy, `d^2 = 0` asserted on every algebra. Positive control for the code (incidence algebras of
posets, HH^m = H^m(order complex)): crown (circle) gives (1,1), commutative square (1), K_{2,3} gives (1,2): all as expected,
so the code does see nontrivial HH when the quiver is not a tree.

| n | LNAs | HH profile (dims of HH^0, HH^1, ..) |
|---|---|---|
| 3..10 | 2, 5, 14, 42, 132, 429, 1430, 4862 | (1) for every one |

Why: Hom_{A^e}(P_m, A) = sum over Bardzell associated paths p of e_{s(p)} A e_{t(p)}; on a linear quiver the parallel path to
p (m >= 2) is p itself, which contains a relation, so the term is 0.

Key groups with >= 2 orbit+mirror classes (orbits = tables + free + edges + doubles, joined with the mirror):

| n | key groups | with >= 2 classes | profile separates | profile equal (unresolved) |
|---|---|---|---|---|
| 6 | 4 | 1 | 0 | 1 |
| 7 | 6 | 2 | 0 | 2 |
| 8 | 11 | 4 | 0 | 4 |
| 9 | 19 | 9 | 1 (F-010) | 8 |
| 10 | 40 | 16 | 3 | 13 |

All class sizes and the per-class profile counts are in `maverick_pq_n10.txt` (HH profile (1) in every row; each class has one
SNF profile). In the 13 equal-profile groups at n = 10 the classes are exactly the unresolved ones of E-115's setting; none is a
certified inequivalence, so none can serve as a control.

## Reproduction

```
.venv/bin/python workshop/rounds/054/maverick_hhsweep.py 10      # about 20 s; output maverick_hhsweep_out.txt
timeout 10m .venv/bin/python workshop/rounds/054/maverick_pq.py 6 7 8 9    # about 15 s
timeout 10m .venv/bin/python workshop/rounds/054/maverick_pq.py 10 > workshop/rounds/054/maverick_pq_n10.txt   # 43 s
```

## Prior record

F-019 / R-008 name Hochschild cohomology as the remaining fallback after the AAG route closed; HYPOTHESES (H-013, around line 603)
lists it as a "would settle" tool; nobody computed it. This closes it for LNAs: it is vacuous (tree quiver, Cibils-type result, cited from memory;
here computed, and the argument above is the Bardzell one). The profile counts at n = 9, 10 (1 and 3 splits) agree with F-047. E-080's "no invariant found" is
extended by one invariant family, with the reason.

## Code changed

None in `quivermutation/`. New: `workshop/rounds/054/maverick_pq.py`, `maverick_hhsweep.py`, outputs `maverick_hhsweep_out.txt`,
`maverick_pq_n10.txt`. No tests touched.

## Next

- Any candidate must be non-vacuous on a tree quiver: derived invariants beyond K_0 + Euler form there are Z-conjugacy of Phi
  (Cartan-level, F-047's SNF is its shadow, already equal on P/Q) and structure of D^b itself (tau-orbits / fractional CY data);
  theorist: is the Euler form's integral congruence class (not only Phi's conjugacy) decided on P versus Q at n = 10?
- The honest route to a control is outside LNAs (gentle algebras with equal Cartan, different AG invariant, n <= 9); but HH-type
  invariants then do not transfer to P/Q since HH is blind on LNAs. Not run.
- A search for a merge (E-064 style) remains the only way P, Q at n = 10 can move; separation has no instrument.
