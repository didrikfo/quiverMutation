# Hom(N,N[-1]) vanishes whenever all J_i vanish, cyclic or not: so "T silting, not tilting iff some J_i != 0" holds unconditionally; only the count changes on a cyclic quiver

author: scholar · round: 035 · kind: proof (with computed checks)
thread: T5 · bears on: E-122, E-126, E-095, E-121

## Claim

Setting as E-126: right modules, v loopless with out-arrows b, D' = (+)P_t(b), g = left multiplication by the b, N = [D' -g-> P_v] (degrees 0, 1), D = (+)_{j != v} P_j, T = D + N (AI 2.31, cited: silting), J_i = {x in e_iAe_v : x b = 0 for all b} (E-122).

(M1) Hom_K(N,N[-1]) = { (y_b)_b : y_b in J_t(b), sum_b b y_b = 0 }, a subspace of (+)_b J_t(b). No acyclicity used.
(M2) Hence Hom(T,T[-1]) = (+)_i J_i  +  Hom(N,N[-1]), and Hom(N,N[-1]) != 0 forces some J_t(b) != 0. So Hom(T,T[-1]) != 0 iff some J_i != 0, for every finite-dimensional A and loopless v. With T silting this gives: the AI 2.32(b) mutation is silting-not-tilting iff some J_i != 0 (E-126's L2 hypothesis "no path from an out-neighbour back to v" is not needed for the iff).
(M3) What changes on a cyclic quiver is the size, not the iff: Hom(T,T[-1]) = sum_i dim J_i + dim Hom(N,N[-1]) with the second term possibly positive (it is 0 on acyclic Q since e_tAe_v = 0 there, E-126 L2). So "Hom(T,T[-1]) = sum J_i" is acyclic-only. Also E-126's L1 (dim J_i <= d_i - 1) needs the gate to test ALL paths into v; the repo gate enumerates simple paths only (`arrowPaths.allPathsBetween`), so on cyclic quivers L1 is not guaranteed (not investigated).
Not claimed: that T is silting (AI 2.31 cited, not re-derived; computed Hom(T,T[1]) = 0 below is only a check); anything about walks or LNAs (their algebras are acyclic).

## Evidence

Proof of M1. A chain map f: N -> N[-1] has only the component f^1: P_v -> D' (f^0 lands in N^{-1} = 0, f^2 starts at N^2 = 0); there are no homotopies (they would be N^2 -> N^0). Write f^1 = (y_b), y_b in e_t(b)Ae_v (Hom(P_v,P_t) = e_tAe_v). The two chain conditions are f^1 g = 0 and g f^1 = 0. First: for all b, b', y_b b' = 0 in e_t(b)Ae_t(b'), i.e. y_b b' = 0 for every out-arrow b', i.e. y_b in J_t(b) (this is E-126's Hom(N,D[-1]) computation, applied to each summand of D'). Second: g f^1 = sum_b b y_b in e_vAe_v. QED. Least certain step: signs in the Hom-complex differential (irrelevant, both conditions are "= 0" on a single component); and that the Hom complex in K^b(proj) is the derived one (true for finite-dimensional A of finite global dimension; in general K^b(proj) is the category of AI's silting theory).
Proof of M2. Hom(D,N[m]) for m = -1 is 0 (N[-1] has no degree-0 term); Hom(D,D[-1]) = 0 (modules); Hom(N,D[-1]) = (+)J_i (E-126). Add M1. If every J_i = 0 (i != v) then every y_b in J_t(b) is 0 (t(b) != v as v is loopless), so Hom(N,N[-1]) = 0.

Computed check (`scholar_hom_nn.py`, own code, truncated path algebra mod p = 32003, Hom-complex homology; independent of `tiltingPlus` and the repo mutation). Columns: dim Hom(N,D[-1]) = sum J, Hom(N,N[-1]), Hom(T,T[-1]); the script asserts Hom(T,T[-1]) = Hom(N,D[-1]) + Hom(N,N[-1]), Hom(N,N[-1]) = M1's formula, Hom(N,D[-1]) = sum J_i, and Hom(T,T[1]) = Hom(T,T[-2]) = 0.

| case | cyclic | sum J | Hom(N,N[-1]) | Hom(T,T[-1]) | Hom(T,T[1]) |
|---|---|---|---|---|---|
| E-078 square, v = d (control) | no | 1 | 0 | 1 | 0 |
| C2 (v->t->v), rad^2 = 0 | yes | 1 (J_t = <c>) | 1 | 2 | 0 |
| C2, rad^3 = 0 | yes | 0 | 0 | 0 | 0 |
| C3 (v->t->x->v), rad^2 = 0 | yes | 1 (J_x) | 0 | 1 | 0 |
| cyclic case 4 (below), L = 6 | yes | 2 | 2 | 4 | 0 |
| case 4 without the relation b y = 0 | yes | 5 | 2 | 7 | 0 |

Explicit cyclic test case (smallest): Q = v -b-> t -c-> v, relations bc = 0 = cb (self-injective Nakayama, dim 4). The vertex v is NOT gate-admitted (c b = 0 is a nonzero path into v killing the only out-arrow), but M1 to M2 do not use the gate. Here J_t = <c>, N = [P_t -b-> P_v], the chain map f^1 = c: P_v -> P_t satisfies c b = 0 and b c = 0, so Hom(N,N[-1]) = k and Hom(T,T[-1]) = 2 while sum J = 1. Cyclic and gate-admitted (case 4, 4 vertices): arrows b: v>t, c: t>x, g: x>v, c2: t>y, d: y>v; y = cg - c2 d; relations y b = 0 and b y = 0, and all paths of length >= 6 zero. `scholar_gate_cyc.py`: the repo gate admits v (True; it sees simple paths only and not the truncation). Computed: J_t = 2 (y and one truncation-induced element), Hom(N,N[-1]) = 2, Hom(T,T[-1]) = 4 > sum J = 2. This is a case where "Hom(T,T[-1]) = sum J_i" is false although the gate admits v; the iff of M2 still holds (both sides nonzero). The truncation makes J_t = 2 (case is not minimal); I did not find a gate-admitted case with J_t of dimension 1.

## Reproduction

```
.venv/bin/python workshop/rounds/035/scholar_hom_nn.py     # about 5 s, prints the table
.venv/bin/python workshop/rounds/035/scholar_gate_cyc.py   # 2 s, repo gate on case 4
```

## Prior record

E-122 left "Hom(N,N[-1]) for silting-not-tilting iff some J_i != 0" open; E-126 (L2) closed it for no-return-path (acyclic) only and the referee asked for a cyclic case or a restriction. grep of `research/` for "Hom(N,N" and "Hom(T,T" finds no entry stating M1 or M3; RETRACTIONS has nothing on it. New, small: M1 (the cyclic description) and M2 (the unconditional iff). Elementary; the same kind of observation as E-126.

## Code changed

None in the library. New: `workshop/rounds/035/scholar_hom_nn.py`, `workshop/rounds/035/scholar_gate_cyc.py`. No tests apply.

## Next

- chair/theorist: when promoting E-126, replace "(hypothesis: no path from t back to v)" by M2 and keep "Hom(T,T[-1]) = sum J_i" acyclic-only; note the gate's simple-path limit for L1 on cyclic quivers.
- skeptic: break M1 (it is a short linear-algebra identity; the AI 2.31 silting step is the cited part, so a non-silting output of the computation would be the news: none seen, Hom(T,T[1]) = 0 in all six cases).
- toolsmith (optional): `isMutable` on a cyclic quiver ignores non-simple paths; either document it or switch to paths in the truncated algebra.
