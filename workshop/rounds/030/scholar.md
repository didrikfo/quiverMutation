# Thin Hom is forced at depth 1 from an LNA and nowhere beyond: dim e_iAe_v <= (out-degree)^depth, and the E-113 data (dim 2, 3) show thinness is false, so it cannot be what bounds the circuit components

author: scholar · round: 030 · kind: result (derivation, short)
thread: T5 · bears on: E-110, E-111, E-113, H-015

## Claim

Let A' = End(T) be any algebra in a walk, P_j = A'-projectives viewed as the summands T_j of the tilting complex T = A'
(over the parent), and mutate at v: triangle T_v -f-> B -> C -> T_v[1], B = (+)_k T_k^{m_k} the minimal left
add(A'/P_v)-approximation (Aihara-Iyama Thm 2.32 setting, as in `research/literature/1009.3370-silting-mutation.md`).
(1) Off-diagonal entries of the child satisfy, with M(v) = sum_k m_k (= number of arrows at v in the relevant
direction, parallel arrows counted; convention UNVERIFIED against the code),
  dim Hom(T_j, C) <= M(v) * max_k dim Hom(T_j, T_k),  dim Hom(C, T_j) <= M(v) * max_k dim Hom(T_k, T_j).
So d(child) <= max(d(parent), M(v) d(parent)) off the diagonal at the new vertex; unchanged vertices keep their Hom.
(2) Hence from an LNA (d = 1, every v has M = 1 on a linear quiver) ONE mutation gives dim e_iAe_v <= 1: thinness
is forced at depth 1. At depth t, d <= prod of M(v_s) over the steps, which for out-degree <= 2 mutations is 2^t.
(3) It is not forced beyond depth 1, and it is false in the data: E-113 records dim e_iAe_v in {2, 3} at admitted
out-degree 2 vertices (c0: 261 rows with 2, 3 with 3; c1: 62 with 2). The tilting-complex structure alone gives only the
geometric bound in (1); nothing in AI/Ladkani/CHZ forces thin Hom for a derived-equivalent algebra of non-PH class.
(4) Therefore the bound does not explain "components of at most 2 edges". It does give: a circuit with k edges needs
dim e_iAe_v >= k (E-113, easy), so k >= 3 needs d >= 3, hence at least two mutations at out-degree >= 2 vertices (or one at
out-degree >= 3) from an LNA. At n = 8 c0/c1 the walks reach depth well beyond that, so the bound excludes only the
shallowest part of the walks. Not claimed: any proof of "no circuit >= 3" or "no nn 2-cycle".

## Evidence

Derivation (the step I am least sure of is the exact form of m_k, marked): apply Hom(T_j, -) to the triangle.
Exactness of Hom(T_j,B) -> Hom(T_j,C) -> Hom(T_j,T_v[1]) and Hom(T_j,T_v[1]) = 0 (T tilting: no Hom(T,T[1])) give
dim Hom(T_j,C) <= dim Hom(T_j,B) = sum m_k dim Hom(T_j,T_k). Dually Hom(C,T_j): Hom(T_v[1],T_j) = Hom(T_v,T_j[-1]) = 0
(negative degrees vanish for a tilting complex) and Hom(B,T_j) -> ... gives dim Hom(C,T_j) <= sum m_k dim Hom(T_k,T_j).
Diagonal Hom(C,C) not bounded here (only off-diagonal entries matter for the Gamma_i circuits).
Depth 1 from an LNA: a Nakayama parent has linear quiver, so the minimal left approximation of P_v is the single
irreducible map P_v -> P_{v+-1} (m = 1), and all dim Hom(P_a,P_b) <= 1 (uniserial, multiplicity-free), so d <= 1.
Not tested numerically this round (no script); the observed maxima (3 at c0) are consistent with d <= 4 only if some
mutation had out-degree >= 2, which E-113 shows exists. A test would be: record (depth, M product, max dim) per row.

Literature counterexamples: none needed to refute forcing; for the PH classes thinness holds because iterated tilted
algebras of type A_n are gentle with unique paths per pair outside oriented 3-cycles (Assem-Happel; recalled, UNVERIFIED,
paper not in `research/literature/`), but c0/c1 at n = 8 are not of that kind, so the PH fact does not transfer.

## Reproduction

None (derivation). Data quoted: E-113 line "observed dim e_iAe_v at admitted out-degree 2 vertices".

## Prior record

`grep -in "thin\|dim e_i" research/FINDINGS.md research/RETRACTIONS.md research/HYPOTHESES.md`: no entry; E-113 states the
dimension counts and the k-edge circuit needs dim >= k, but not the depth bound (1). Out-degree >= 3 has 0 rejects (E-111)
and 3 appears as dim only 3 times at c0. Not in RETRACTIONS. The inequality (1) is the standard cone estimate; I did not find
it in 1009.3370 / 2112.08129 / 2509.12983 summaries (those state tilting criteria, not Hom bounds).

## Code changed

None.

## Next

- experimentalist: per row log (depth, product of out-degrees of the mutated vertices along the walk, max dim e_iAe_v) at
  n = 8, 9: tests (1) and shows how loose it is; requested already as "max dim by depth", this sharpens it.
- theorist: the useful target is not a dimension bound but a sign/shape one: show the product p1 b1 = p2 b1 with
  b2 killing p2 (half-W) is the only way two classes share a v, using that C's Homs are cokernels/kernels of a single map
  (g_i of E-110) rather than sums; a circuit of k edges then needs k independent classes in ONE kernel of a map between spaces of
  dimension <= M d.
- scholar: read Assem-Happel / Buan-Vatne for thin Hom on iterated tilted A_n (outside the repo, network permitting).
