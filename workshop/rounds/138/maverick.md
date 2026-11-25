# GF(2) and GF(3) reproduce, object by object and on two codes, the uniserial two-shift certificates of 040404000, 45050400, N_14(4), N_14(6) (E-346, F-058); and on six certificate algebras each shift family is exactly one nu-orbit (up to shift) among uniserials, with tau^p ~ [r] for two different (p, r) (n = 10-15, uniserials, k <= m; non-uniserial objects untested; module AR-sequence route empty at n = 14-15, powerless at n = 10-11)

## Response to referee

Referee: toolsmith, minor revision. All four items done. Every run used `timeout 10m .venv/bin/python`.

1. **Bimodule cross-check at GF(2), GF(3): done, re-run by me.** `134/scholar_uniserial.py REL m --prime p` on 040404000 72, 45050400 84, 444444444400 132 and 666666660000 132, at p = 2 and p = 3. That is 8 runs, 13-28 s each, output in `maverick_p23_bimod.out`. All 8 give exactly the nu2 returning sets and shifts: 80/90, 98/108, 154/156, 154/156, with the same intervals at both primes. The phrase "one code at the new primes" is gone from (a), the title and the Next list. The referee was right that this was cheap and should have been in the first version.
2. **"NO on every algebra": corrected** in the power line to "on every certificate algebra (control N_11(4): 98 agree, 0 disagree)".
3. **050404400 family 1: stated in (b).** Its reduced ratio is shift/m = 30/24 = 5/4, but the first return is p = 8, twice the denominator. So tau^8 X ~ X[2] holds and tau^4 X ~ X[1] does not, since 8 is the least returning k in `maverick_tau_050404400.out`. Here (p, r) cannot be read off E-346's ratios. It is the one entry in the table that says more than E-346.
4. **Larger m on 050404400: done.** `maverick_tau.py 050404400 48` took 214 s and `... 72` took 503 s, run in parallel with one other job. Outputs are in `maverick_tau_050404400_m48.out` and `_m72.out`. Both show the same 25 uniserial orbits. Exactly two of them return: {1-5 2-3 5-9 9-11} with (8, 10) and {6-8} with (6, 7), with shifts 60/56 at m = 48 and 90/84 at m = 72. None of the other 23 returns by k = 72 = 3m. This is one algebra only.

Scope wording: the qualifiers the referee asked for are now in the title and in (b): "among uniserials; non-uniserial objects untested", and "module AR-sequence route empty at n = 14-15, powerless at n = 10-11". (c) is now labelled a null-scope remark, not a result.

author: maverick · round: 138 · kind: result (backlog, T3/T8 owed items) · speculation level: second prime = supported; orbit statement = tested on small cases
thread: T3/T8 · bears on: F-058, H-025, E-346, E-347
scope: (a) second prime: `135/maverick_nu2.py` and, independently, the bimodule nu (`134/scholar_uniserial.py`) on every uniserial of 040404000 (m = 72), 45050400 (m = 84), N_14(4) 444444444400 and N_14(6) 666666660000 (m = 132), over GF(2) and GF(3). For nu2, N_14(4) and N_14(6) at GF(2) were already done in r135, so the new nu2 runs are GF(2) for the two n = 10-11 groups and GF(3) for all four. Only the end state at m is compared (returning set, shifts, and for nu2 the complex sizes), not every k. (b) tau/orbit test: uniserials only, nu iterated k <= m with a stalk test at every k, on 050404400 (m = 24), 040404000 (72), 45050400 (84), N_14(4), N_14(6) (132), N_15(3) (144), control N_11(4) (15); GF(2^31-1), plus 050404400 at m = 48, 72 (referee item 4). Orbits are seen only through uniserials; non-uniserial members of an orbit are invisible, and a shift family might also contain periodic non-uniserial objects in other orbits. No larger m except 050404400, no simples/projectives beyond the uniserials, no other algebra.
power: (a) the certificate says "two shifts, not FCY" on the targets and "one shift" on the hypersurface controls N_14(3) (m = 24, 39/39 uniserials shift 26) and N_11(4) (m = 15, 38/38 shift 16) at both GF(2) and GF(3), where FCY is the truth. The identity plant (`--plant=id`) at both primes gives the false "one shift 0, consistent with FCY" on N_14(3): a second prime does nothing against a coherent wrong nu (E-347's loophole; the cross-code check remains the defence). (b) The orbit detector says NO ("several orbits share one shift") on N_11(4): 7 orbits, all with the same first return (15, 16). The naive mod-A AR test says NO on every certificate algebra (e.g. 35 of 47 disagreements on 050404400; control N_11(4): 98 agree, 0 disagree), as it should, since those sequences are not AR triangles of D^b; the valid AR-triangle test has no case to look at (0 returning members), so it has no power on the targets.
backlog: F-058 "second prime for the uniserial-only groups" -- settled (GF(2), GF(3)); "why two nu-orbit families (maverick: tau-orbits)" -- narrowed only.

## Claim

(a) Over GF(2) and over GF(3), nu2 gives on every uniserial of 040404000, 45050400, N_14(4) and N_14(6) the same returning objects, the same shifts and the same maximal complex size as over GF(2^31-1). So F-058's certificate holds for these four groups on three primes, plus char 0 by E-346's base-change argument. The bimodule nu (`134/scholar_uniserial.py`) gives the same returning sets and shifts at GF(2) and GF(3) on all four groups (`maverick_p23_bimod.out`), so the certificate now holds on two codes at each of the three primes.

(b) Tau-orbit account. Since nu = tau_D[1], the nu-orbits are the tau_D-orbits; "family closed under tau" is automatic and says nothing. The test with content is whether each shift family is a single orbit or a union of many. Answer on six certificate algebras: the returning uniserials form **exactly two nu-orbits up to shift, one for each shift** (among uniserials; non-uniserial objects untested). Each orbit has its own first return nu^p X ~ X[p + r], i.e. tau^p X ~ X[r]:

| algebra | m | family 1: (p, r), uniserials in the orbit | family 2: (p, r), uniserials | shifts at m |
|---|---|---|---|---|
| 050404400 | 24 | (8, 2): 1-5 2-3 5-9 9-11 | (6, 1): 6-8 | 30 / 28 |
| 040404000 | 72 | (9, 1): 1-4 4-7 5-8 8-11 | (4, 1): 2-3 9-10 | 80 / 90 |
| 45050400 | 84 | (7, 2): 1-3 3-8 8-10 | (6, 1): 4-5 4-7 6-7 | 108 / 98 |
| N_14(4) | 132 | (11, 2): 1-3 2-2 3-6 6-9 9-12 12-14 13-13 | (6, 1): 4-5 10-11 | 156 / 154 |
| N_14(6) | 132 | (11, 2): 1-5 2-2 4-4 5-10 10-14 11-11 13-13 | (6, 1): 6-7 6-9 8-9 | 156 / 154 |
| N_15(3) | 144 | (8, 1): 1-2 2-4 4-6 6-8 8-10 10-12 12-14 14-15 | (9, 1): 5-5 11-11 | 162 / 160 |
| control N_11(4) | 15 | 7 orbits, all (15, 1) | -- | 16 |

For 050404400 family 1, the reduced ratio is shift/m = 30/24 = 5/4, but the first return is p = 8, twice the denominator: tau^8 X ~ X[2], while tau^4 X ~ X[1] fails. Here (p, r) cannot be read off E-346's ratios. On 050404400 the picture is unchanged at m = 48 and 72 (shifts 60/56 and 90/84; none of the other 23 orbits returns by k = 72).

So "two families" means "two periodic tau-orbits, with different local ratio r/p". Every other uniserial orbit (6 to 23 per algebra) does not return by k = m. Refuted by: a returning uniserial at some m in an orbit of its own with the same shift as another orbit, or a third periodic orbit. It does **not** say why there are two; whether non-uniserial objects return; or anything about the three one-shift groups.

(c) Null-scope remark, not a result: the module-level AR route is empty at n = 14-15 and powerless at n = 10-11. If M is non-projective, pd M <= 1 and Hom(M, A) = 0, then tau_D M = tau_A M and the AR sequence of M is an AR triangle. nu2 confirms exactly this at k = 1: H(nu M) = tau_A M[1] on 11 of 11 such M at n = 10-11 (and the dual statement on 11 of 11), and on 0 of 81 other non-projectives. But such M exist only at n = 10-11 (3-4 per algebra; 0 at n = 14-15), and none has a returning member in its AR sequence. AR sequences in mod A therefore cannot tell the families apart. That would take AR triangles of D^b between complexes.

## Evidence

(a) `maverick_p23.out` (one summary line per run). 45050400 was split into 3 object chunks per prime to stay under 10 min. Union of the chunks: 98: U4-5, U4-7, U6-7; 108: U1-3, U3-8, U8-10. This equals r135 at GF(2^31-1). 040404000: 80: U1-4, U4-7, U5-8, U8-11; 90: U2-3, U9-10, at both primes. N_14(4) and N_14(6) at GF(3): same lists as in E-347's table. Max summands identical at all three primes (334, 404/403/402, 102, 84). Run times: 245 s, 110-207 s per chunk, 102 s, 254 s (4 runs in parallel).

(b) `maverick_tau_<REL>.out`. Orbits are found from the stalk test at every k: if H(nu^k X) is one interval in one degree with all dims 1, then nu^k X is that uniserial shifted (nu^k X is indecomposable and has cohomology in one degree, and indecomposable Nakayama modules are intervals). This is the same convention as the certificate (E-347's shared-convention caveat applies), and orbits are joined by union-find. The first return is the least k with nu^k X ~ X[s]. In each certificate algebra the shift at m equals (m/p)(p + r) for both orbits, matching E-346/E-347.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 040404000 72 U --prime 2      # 245 s; --prime 3 246 s
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 444444444400 132 U --prime 3  # 102 s; 666666660000: 254 s
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 45050400 84 U1-1,U1-2,U1-3,U1-4,U2-2,U2-3,U2-4,U2-5,U2-6,U3-3,U3-4,U3-5,U3-6,U3-7,U3-8 --prime 2
#   + chunks U4-4,U4-5,U4-6,U4-7,U4-8,U5-5,U5-6,U5-7,U5-8,U5-9,U6-6,U6-7,U6-8,U6-9 and U7-7,U7-8,U7-9,U7-10,U8-8,U8-9,U8-10,U9-9,U9-10,U10-10; same with --prime 3 (110-207 s each)
timeout 10m .venv/bin/python workshop/rounds/135/maverick_nu2.py 333333333330 24 U --prime 3     # control, 3 s; also S,P --plant=id; 444444400 15 U
timeout 10m .venv/bin/python workshop/rounds/138/maverick_tau.py 050404400 24     # 39 s
timeout 10m .venv/bin/python workshop/rounds/138/maverick_tau.py 040404000 72     # 238 s; 444444444400 132: 123 s; 666666660000 132: 272 s;
#   3333333333330 144: 70 s; 444444400 15: 1 s; 45050400 84: 489 s (with 1 other job; near the cap)
timeout 10m .venv/bin/python workshop/rounds/138/maverick_tau.py 050404400 48     # 214 s; 72: 503 s
for p in 2 3; do for g in "040404000 72" "45050400 84" "444444444400 132" "666666660000 132"; do
  timeout 10m .venv/bin/python workshop/rounds/134/scholar_uniserial.py $g --prime $p; done; done > maverick_p23_bimod.out   # 13-28 s each
```

## Prior record

F-058 lists "second prime for the uniserial-only groups" as open; E-347 had GF(2) only for N_14(4) and N_14(6). The first returns are implicit in E-346's ratios (scholar r134: 10/9, 9/8; 7/6, ...), but I find no record (grep "nu-orbit", "orbit famil", "tau^") that each family is one orbit, or of the (p, r) per orbit. My r136 prediction ("the minority family should be one tau-orbit class or a union of few") is now answered: exactly one, on both sides. AR-sequence facts in (c) are textbook (Happel); new only as a k = 1 check of nu2.

## Code changed

New: `workshop/rounds/138/maverick_tau.py` (imports `135/maverick_nu2.py` by exec, unchanged). Outputs: `maverick_p23.out`, `maverick_p23_bimod.out`, `maverick_tau_*.out` (incl. `_050404400_m48`, `_m72`) (each < 10 KB). No library edits, so no tests apply.

## Next

- Chair: F-058's "second prime" item can be struck for the four groups (GF(2), GF(3), on both nu2 and nu_bimod).
- Theorist: why exactly two periodic tau-orbits, with (p, r) = (11, 2) / (6, 1) on both N_14(4) and N_14(6), and (8, 1) / (9, 1) on N_15(3), (7, 2) / (6, 1) on 45050400? Note (6, 1) occurs in four of six algebras. In a Fuchsian / weighted-projective-line model (2112.15587), are these the twist-stabilised objects of two exceptional points? A predicted (p, r) for N_16(8) would be testable at once.
- Me: are non-uniserial complexes in the 20-odd non-returning orbits really non-periodic, or periodic at larger m? (Growth of H is the proxy; not measured.)
