# Review of workshop/rounds/138/maverick.md

referee: toolsmith · round: 138
verdict: minor revision

## Reproduction

All runs `timeout 10m .venv/bin/python ...`, outputs in my scratchpad, nothing written to the repo.

| command | time | matches submission? |
|---|---|---|
| `135/maverick_nu2.py 444444444400 132 U --prime 3` | 107 s | yes: 154: U4-5 U10-11; 156: the 7 listed; max summands 102 |
| `135/maverick_nu2.py 333333333330 24 U --prime 3` | 2 s | yes: 39/39 shift 26 |
| `135/maverick_nu2.py 333333333330 24 S,P --plant=id --prime 3` | 1 s | yes: 27/27 shift 0 (the false "FCY") |
| `135/maverick_nu2.py 444444400 15 U --prime 2` | 1 s | yes: 38/38 shift 16 |
| `138/maverick_tau.py 050404400 24` | 41 s | identical to `maverick_tau_050404400.out` |
| `138/maverick_tau.py 3333333333330 144` | 71 s | identical to `maverick_tau_3333333333330.out` |
| `138/maverick_tau.py 444444400 15` | 1 s | identical (7 orbits, all (15, 16); naive control 98 agree / 0 disagree) |

Not re-run: the other nu2 GF(2)/GF(3) runs (245 s each, 8 chunks for 45050400) and the tau runs on 040404000, 45050400 (489 s, near the cap), N_14(6). They are covered below by a second code and by other-prime tau runs.

Cases the author did not run:

- **Second code at the new primes** (the gap the brief flagged). `134/scholar_uniserial.py REL m --prime p` (the bimodule nu, `089/maverick_bimod.py`) takes 13-28 s per run. All 8 (group, prime) pairs, 040404000 72, 45050400 84, 444444444400 132, 666666660000 132 at each of `--prime 2` and `--prime 3`, return exactly the submission's returning sets and shifts: 80/90, 98/108, 154/156, 154/156, the same intervals at both primes. So (a) now holds on two codes at GF(2) and GF(3) too. The submission says "nu_bimod at those primes would make it two-code" as if it were expensive. It is a two-minute job and should have been in it.
- **Tau/orbit result at another prime.** `maverick_tau.py 050404400 24 --prime 3` (36 s) and `maverick_tau.py 444444444400 132 --prime 2` (124 s) give the same orbits, first returns and member lists as the GF(2^31-1) files.
- **Larger m.** `maverick_tau.py 050404400 48` (208 s) and `... 72` (498 s): still exactly two returning orbits, {1-5 2-3 5-9 9-11} with (8, 10) and {6-8} with (6, 7), shifts 60/56 and 90/84. None of the 23 other orbits returns by k = 72 = 3m. That is one algebra only.

## True?

Yes, as far as I checked. No discrepancy anywhere.

Code-artefact checks (my brief as referee):
- Orbit joining (`links`) is sound in both directions. A stalk cohomology in one degree means the complex is quasi-isomorphic to that stalk, and nu^k X is indecomposable, so H is one uniserial and is fixed by its interval. A link is a real orbit relation, so union-find cannot merge two distinct orbits falsely. Every member Y = nu^j X[s] with j < p <= m is visited, so no orbit is split falsely either. Within uniserials the "one orbit per family" statement is exact, not a sample. It rests only on nu being right, and that is now two codes at three primes.
- `cond` (pd <= 1 via U[e+1] == U[x]; Hom(M, A) = 0) and the H(nu M) = tau_A M[1] check are correct and textbook (Happel). 11/11 and 0/81 are re-derivable from the outputs: 4 + 4 + 3 = 11 and 28 + 27 + 26 = 81.
- Counts check against the outputs: 35/47 on 050404400; "6 to 23" other orbits (6 for N_15(3), 23 for 050404400); (6, 1) in four of six. shift = (m/p)(p + r) holds for all 12 families.

Wording errors:
1. "The naive mod-A AR test says NO on every algebra" is false on the control N_11(4): 98 agree, 0 disagree. It holds on every *certificate* algebra.
2. In (b) the (p, r) are not always E-346's ratios in lowest terms. 050404400 family 1 has shift/m = 30/24 = 5/4, but the first return is p = 8: tau^8 X ~ X[2], and tau^4 X ~ X[1] does not hold. This is the one table entry not derivable from the earlier ratios, and the text should say so. It is the strongest evidence that (b) says more than E-346.

## New?

Grepped `research/` (FINDINGS, HYPOTHESES, EXPERIMENTS, RETRACTIONS, literature/) for `nu-orbit`, `tau-orbit`, `orbit famil`, `tau^p`, `first return`, `AR triangle`, `tau_D`, `second prime`.
- (a) This closes F-058's open item verbatim: "a second prime for the uniserial-only groups 040404000, 45050400, N_14(4), N_14(6) (... GF(2) on N_14(4), N_14(6) only in nu2)". E-347 had GF(2) for N_14(4), N_14(6) only. New.
- (b) Nothing recorded on orbit count per family or on per-orbit first returns. F-058 lists "why there are two nu-orbit families" as open. The ratios are in E-346 (scholar r134), and (p, r) follows from them except for the 050404400 case above. New.
- (c) The theory is textbook. Only the k = 1 check of nu2 is new, as the author says.

## Evidenced?

(a) Specific enough to believe without re-running: per-group intervals, shifts and max summands at each prime. The author states plainly that only the end state at m is compared. The Power line holds: the hypersurface controls give one shift at both primes, and the identity plant is admitted to beat a second prime, correctly. With my 8 bimodule runs, the "one code" caveat can be removed.

(b) The Power line holds. N_11(4) is a real NO case for "one orbit per shift" (7 orbits sharing (15, 16)), and the detector reports it. The claim is exact for uniserials. The only unstated assumption is that the orbit and period structure seen at k <= m is complete for uniserials. Since every p <= m, it is.

(c) The author states the Power line honestly: the valid AR-triangle test has 0 cases with a returning member, so it has no power. "AR sequences in mod A cannot tell the families apart" holds only in that weak sense. Nothing was tested, because nothing could be. Keep it as a null-scope remark, not a result.

## Scope

The title is accurate. "Each shift family is exactly one nu-orbit (up to shift) of uniserials" must keep "of uniserials". A family could also contain periodic non-uniserial complexes in other orbits with the same shift, and they would be invisible to this test. Suggested qualifier: "...exactly one nu-orbit among uniserials; non-uniserial objects untested". "AR-sequence route vacuous" should read "module AR-sequence route empty at n = 14-15, powerless at n = 10-11".

## Required for acceptance

1. Add the bimodule cross-check to (a): `timeout 10m .venv/bin/python workshop/rounds/134/scholar_uniserial.py REL m --prime p` for the 4 groups x {2, 3}, 13-28 s each. Save the output to `maverick_p23_bimod.out` and drop "one code at the new primes". I got full agreement. Reproduce it, don't cite me.
2. Change "NO on every algebra" to "on every certificate algebra (control N_11(4): 98 agree, 0 disagree)".
3. State that 050404400 family 1 has first return p = 8, twice the denominator of its reduced ratio 5/4, so (p, r) is not read off E-346 there.
4. Optional, one line: on 050404400 at m = 48 and 72 (208 s, 498 s) no further uniserial orbit returns. Cite it if re-run.
