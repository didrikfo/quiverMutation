# Of the two papers 0911.5137 and 1310.1557 only (not the wider literature), neither decides Candidate C at n = 11: both give only necessary conditions (Phi periodic, Phi^a = (-1)^(a+b) I)

author: scholar · round: 082 · kind: rediscovery
thread: T3/T8 · bears on: E-190, E-197, E-209, E-215, F-045, F-047
scope: two papers only; Happel-Seidel, Lenzing-de la Peña and others not checked (UNVERIFIED). Literature check against the arXiv LaTeX source (fetched this round: 0911.5137 `tensor.tex` + `.bbl`, 1310.1557 `Algebras_cyclotomic_type_05.2013.tex`); numerics at n = 11 only (all 638 periodic-Phi indefinite LNAs, exact integers). No object-level S^a was computed.

## Response to referee

1. Prime-rank argument (accepted). Claim 2 now reads: the only FCY LNAs these papers give *that can occur at n = 11* are the lines. Argument: every FCY construction in 0911.5137 that yields an LNA (Cor 1.2 tensor kA_m (x) kA_n, Cor `c:fracCY` T_n(Lambda), tensor products of FCY algebras) has rank multiplicative in the factors (rank = rank(Lambda)*n, or m*n). 11 is prime, so one factor has rank 1 (trivial). The only nontrivial output is then kA_11 itself, i.e. A(11,2) = kA_11 (m = 1 or n = 1 in Cor 1.9), which is FCY but Dynkin, hence outside the 638 indefinite LNAs. So "Cor 1.9 gives nothing at n = 11" is corrected to "gives only Dynkin". This is a statement about those constructions, not about FCY LNAs in general.
2. Windows (accepted). The table lists the minimal a with Phi^a = +-I; every multiple of the period also gives Phi^a = +I (a+b even). Only the 535 and 15 groups have an odd-a window (a = 15 and 21, the half periods, Phi^a = -I); 64 and 24 admit only even a+b.
3. Misprint (accepted). The lcm phi(m) vs lcm m misprint is already recorded in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` lines 37-42 (checked 2026-09-19 against N_15(5)). Not new here; only the LaTeX confirmation is added.
4. Runtime (fixed): 19 s, not 19 s. Referee's re-run matches the table. Also noted: the claim of the title is now scoped to the two papers.

## Claim

Candidate C is "S^a(A) ~ A[b]" (Serre functor S = nu), i.e. fractionally Calabi-Yau (FCY) with the derived-invariant exponent (a, b). What the two papers give for it:

1. **Only a necessary condition.** 1310.1557 Thm 2(a) (intro l. 221): FCY => Phi periodic. Its Lemma 2.9 gives phi^(2q) = 1; the sharper form used by E-190, Phi^a = (-1)^(a+b) I, follows from [S] = -Phi and [1] = -id on K_0 (both papers define S = tau[1], = nu). **No statement in either paper goes from periodic Phi to FCY.** Thm 2(b) is the converse only for "every component of the class quiver [Gamma] is a tube", an AR-quiver property we cannot compute. So periodic Phi (the 535/64/24/15 groups) neither proves nor refutes Candidate C.
2. **The only FCY LNAs the papers hand us (and, by primality of 11, the only ones their constructions can give at n = 11, besides Dynkin) are the lines A(mn, m+1) = kA_N/rad^(m+1)**, 0911.5137 Cor 1.9: ((m-1)/(m+1) + (n-1)/(n+1))-CY, from Cor 1.2 (A(mn,m+1) ~ kA_m (x) kA_n). It needs N = mn composite. **n = 11 is prime, so Cor 1.9 gives only the Dynkin kA_11 at n = 11 (m = 1 or n = 1)**; at n = 10 it covers m = 2, n = 5 (L = 3, 6) but not 90000000 = A(10,9), whose S^9 = [8] comes from D_10 (E-190, Happel-Seidel), not from this paper.
3. **Nothing non-Cartan.** The papers' invariants are chi_A, the form h_A (1310.1557 Sect. 2.6, 6.1) and, in 0911.5137, the CY fraction. The fraction is the object-level content of Candidate C; Remark 1.14 shows how a paper separates fractions (orders of products of pairs of roots of the two Coxeter polynomials), which is again Cartan-level here, all four groups having one chi.
4. **Misprint confirmed in the LaTeX and on our data.** 1310.1557 states "the roots of Phi_n are primitive phi(n)-roots of unity" (Sect. 2.2) and Thm 2(c),(d) / Thm 2.10(3) use lcm{phi(m)}; the proof of 2.10 "(1)=>(3)" uses lambda^phi(m) = 1, false. The period is lcm{m}. At n = 11 all four groups agree with lcm{m}, not lcm{phi(m)} (table). Already recorded in `research/literature/1310.1557-algebras-of-cyclotomic-type.md` lines 37-42; not new, only checked against the source.

New small facts at n = 11 (idea level; the factorisations were already in round 066's table): the 535-group contains the radical powers kA_11/rad^4 and kA_11/rad^5 (found by digit string str(L)*(11-L)+"0"*(L-2)), so it is the one group with a Happel-Seidel / Lenzing-de la Peña anchor; whether they are FCY is **not stated in the repo's literature files** and not in Cor 1.9 (UNVERIFIED).

## Evidence

Theorem numbers matched to the LaTeX counters (`\newtheorem{lemma}{Lemma}[section]`, `prop` and `cor` share it; the introduction of 0911.5137 is Sect. 1 "The results"):

| summary says | source | status |
|---|---|---|
| 0911.5137 Thm A / B / C | `t:tensor` l. 272 (A), `t:ade` l. 363 (B, tridiagonal), `t:endo` l. 526 (C, full triangular); `\thethm = \Alph` | matches |
| Cor 1.2 (rectline), 1.3 (replica, Gorenstein), 1.4 | l. 439, 488, 510 (count: 1.1 itertilt, 1.2, 1.3, 1.4) | matches |
| Cor 1.7 (F-form), 1.9 (line CY), 1.10 (preproj), 1.11 (Aus), 1.12 (stable Aus), 1.13 (rectriang), Remark 1.14 | l. 588, 671, 695, 770, 851, 967, 972 (remark shares the counter: 1.5 is the remark at l. 544) | matches |
| 0911.5137 Prop 3.1 (RA_n tilting) | `p:RAntilt` l. 1395, first numbered item of Sect. 3 | matches |
| 1310.1557 Thm 1 (six readings), Thm 2 (a)-(e) | l. 182, 217; (c),(d) misprinted as above | matches / misprint |
| 1310.1557 Lemma 2.9, Thm 2.10, Thm 6.1, Sect. 6.2 | l. 407, 425, 1103, 1138 | matches; 2.9 gives only 2q |
| "h_A pos. def. iff derived Dynkin, citing [33]" | intro/Sect. 6 text; the Sect. 6.2 pair (A_n, B_n) has h equal, one finite, one wild | the summary's warning stands |

Fractional-CY definition in the two papers: 0911.5137 l. 628 `nu^e = [d]`, "common factors cannot always be cancelled"; 1310.1557 l. 399 `S^q = [p]`, minimal q. Same notion.

n = 11 table (`scholar_cyclo.py`, factor via sympy, order by exact powers; 19 s):

| group | chi (cyclotomic index, mult.) | ord Phi = lcm m | lcm phi(m) (misprint) | Phi^a = +-I (E-190 convention) |
|---|---|---|---|---|
| 535 | Phi_2 Phi_6 Phi_30 | 30 | 8 | a = 15: a+b odd; a = 30: a+b even |
| 64 | Phi_2 Phi_8 Phi_18 | 72 | 12 | a = 72: a+b even |
| 24 | Phi_2 Phi_4 Phi_8 Phi_12 | 24 | 4 | a = 24: a+b even |
| 15 | Phi_2^3 Phi_6 Phi_14 | 42 | 6 | a = 21: a+b odd; a = 42: a+b even |

Table entries are the minimal a; multiples of the period also work (a+b even). Only the 535 and 15 groups admit a *half* period (Phi^a = -I) hence an odd a+b window; 64 and 24 allow even a+b only. These are necessary windows, not evidence for FCY.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/082/scholar_cyclo.py 11     # 19 s, prints the table above
```

LaTeX source: `curl -L https://arxiv.org/e-print/0911.5137` (tar.gz) and `.../1310.1557` (single gzipped tex); the second answered "Rate exceeded" once, a retry after 20 s worked. Not committed.

## Prior record

E-190 (Cartan table at n = 10, four n = 11 groups), E-197, E-209, E-215 and round 066's group table already hold the factorisations and the "Cartan probes exhausted" verdict. Literature summaries of both papers already exist and are correct; I found no error in them beyond what they flag themselves. New: the theorem-number audit, the confirmation that **neither paper contains a periodic => FCY statement** (so no shortcut for the 15-group), and the observation that Cor 1.9 is silent at prime n.

## Code changed

New `workshop/rounds/082/scholar_cyclo.py` only. No tests apply (not a library file).

## Next

**One concrete test for the maverick (object level, can say NO).** For LNA A with Kupisch data, compute the *minimal* complexes of projectives X_k = nu^k(P_i), k = 1..a, i over all vertices, by iterating nu = - (x)^L DA (nu P_i = I_i, a module with a short uniserial projective resolution; apply termwise and resolve the cones). Candidate C holds with exponent (a, b) iff X_a = P_i[b] for every i; the *growth of the total rank of X_k* is a built-in NO signal, since FCY forces periodic ranks. Controls: kA_10/rad^9 (a = 9, b = 8, a YES), kA_10/rad^3 (Cor 1.9, m = 2, n = 5, YES with fraction 1/3 + 4/6 = 1), and a non-FCY LNA. Stage on n = 11: kA_11/rad^4 and rad^5 in the 535-group (a = 15 or 30), then one member of the 15-group (a = 21 or 42). Size with `--plan` first; the rank growth tells quickly whether a = 42 is affordable. A NO on any member plus a YES on another member of one group would be the first object-level separator.
Theorist: is ranks bounded for a *periodic but not FCY* Phi realisable by an LNA? Skeptic: re-derive Phi^a = (-1)^(a+b) I with the repo's Phi = -C^-T C.
