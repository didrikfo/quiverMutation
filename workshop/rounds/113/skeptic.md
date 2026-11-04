# Pre-registered NO rule for E-201's 13 steps, and planted controls: tested on small cases (controls n <= 6 LNAs), a Hom_K(T,T[m]) test says NO to a non-tilting T, and a K0 check is needed beside it; no control here speaks to "leaves the class"

## Response to referee

Verdict minor revision; all four required items done, numbers re-run after the edit.

1. P3 positive control printed (script now reports `Rcart`/`Rcartmatch`, the true child k2 == k compared with End(T_k)): n=4 20/20, n=5 70/70, n=6 64/64 match (referee's n=5 70/70 confirmed). So the C.T convention is right and the P3 mismatches (40/40, 210/210, 256/256) are not a convention artefact.
2. Accepted: P3's child Cartan is the repo's own `mutation` + `reduction` output, not independent. P3 therefore tests only the convention and labelled matching of the Hom side against the repo, and stays a weak control.
3. Accepted: K0-unimodularity is already part of the tilting test in E-287 ("unimodular K0"; generation unverified). P2 is not a newly found gap: it is a measured false-YES rate (262/262) of a Hom-only test, i.e. a control for the K0 check E-287 already has. Credited below under Prior record.
4. Dead `Tz` / first `Tz2['k']` lines deleted; re-run gives the same table (P1 30/112/120, P2 30/112/120, Rtilt 20/70/64, Rk0 30/112/120).
Also taken: P1-P3 test the detector, not the claim "13 of 13"; the title no longer claims power statements about "leaves the class" beyond the absence of a key-equal inequivalent pair at n <= 7 (E-212), not shown by this script. Newness: referee is right that the rule and the "Hom failure is not NO" point restate E-201/E-212; the new part is the planted controls only.

author: skeptic · round: 113 · kind: proposal
thread: T10 · bears on: H-015
scope: tested on small cases only, controls n <= 6 LNAs (not the 13 steps; P1-P3 test the detector, not the claim). Rule written before seeing any toolsmith output (I read none). Controls: all LNAs of length 4 and 5 with relation duals (10 and 28 algebras), 12 LNAs + duals at length 6 (24 algebras), every vertex with an out-arrow (30, 112, 120 steps); the two-term T of `rounds/050/skeptic_tilt.py` only. Not run: n = 7, the 13 steps, generation, End(T) as an algebra.
power: Controls where the Hom test says NO and NO is true by construction: P1 zero differential 262 of 262 NO; P3 End(T_k) against the Cartan of the child at another vertex 506 of 506 NO. Where it must say NO and does NOT: P2 duplicated summand (not a generator) 262 of 262 pass Hom; only the K0 check catches it (262 of 262). For "a guarded walk leaves its class" no control exists at n <= 7 (no key-equal inequivalent pair, E-212); the only candidate is the n = 9 F-010 pair, 5 + 5 hit the cap (E-212).

## Claim

Pre-registered, before any output:

**What a single-step Hom test can and cannot say.** The 13 steps all have dim Hom(T_v,T_v[-1]) != 0 (E-201, independent Hom test): that is already known and says "this two-term T is not a tilting complex", i.e. the step is not an exchange by T. It is NOT a NO for derived equivalence of parent and child (E-212: 8 of 9 class-2 steps with Hom(T,T[-1]) = 1 are joined by Hom-verified chains). So any rule that reads "Hom(T,T[m]) != 0 at the step" as NO is refuted in advance, and I will not accept it.

**YES (child in the parent's class), all of:** an explicit chain parent = X_0, ..., X_r (r >= 1 each side allowed, meet in the middle) in which every edge passes (a) Hom(T,T[-1]) = Hom(T,T[1]) = 0 for the edge's T, (b) K0 classes of the n summands have determinant +-1 (generation proxy; Okuyama-Rickard completes it, stated as assumed), (c) dim Hom(T_i,T_j) = Cartan of the next algebra, labelled, and the two sides meet in an explicit isomorphism of presented algebras (not only equal `canonicalKey`). Per step; the claim "13 of 13" needs all 13, with a record number for each.

**NO (a key-guarded walk leaves its class), only one of:** (i) a derived invariant that differs between a key-kept parent and its child and has a recorded power control: Cartan Z-congruence fails (already 0 of 89, E-295); Coxeter polynomial, det, SNF, signature differ (equal by the key, so cannot fire); a finer invariant (e.g. Hochschild/gl.dim-free) with its own YES/NO controls; or (ii) a complete finite proof, e.g. an exhaustive classification of tilting chains up to a stated depth together with a proof that the depth bounds the class. "Chain search found no meet to depth d" is NOT NO: it is a bound (E-156, F-010 4 + 4).

**Inconclusive:** a chain search ends without a meet, within budget, while all invariants agree. At n = 7 I expect this for any step the toolsmith cannot join; it leaves T10 open, and must be written as "no join to depth d + d", never as "outside".

**What the toolsmith's test must pass to count as having power** (all three, the first two on a set that includes n = 5 and 6 and this script's controls, which I put in `skeptic_homctl.py`): (1) reject P1 (zero differential, 262 of 262 here); (2) reject P2 (duplicated summand): a Hom-only test accepts it (262 of 262 false YES in my run), so a K0-determinant (or generation) check must be in the pipeline and shown to reject it; (3) reject a wrong labelled child (P3, 506 of 506 here); and (4) for the chain-search NO, a pair at the same budget that is inequivalent with all Cartan-level invariants equal, on which it returns no meet and a YES pair on which it meets: only the n = 9 F-010 pair is such, and it is capped. Without (4), a chain-search failure carries no "outside" weight.

## Evidence

P3 positive control (true child, k2 == k, labelled Cartan equals End(T_k)): n=4 20/20, n=5 70/70, n=6 64/64.

| set | steps | real T Hom-tilting | real T K0 det = 1 | P1 NO | P2 Hom says tilting | P2 K0 det = 0 | P3 mismatch |
|---|---|---|---|---|---|---|---|
| n=4 (10 alg) | 30 | 20 | 30 | 30/30 | 30/30 | 30/30 | 40/40 |
| n=5 (28 alg) | 112 | 70 | 112 | 112/112 | 112/112 | 112/112 | 210/210 |
| n=6 (24 alg, capped) | 120 | 64 | 120 | 120/120 | 120/120 | 120/120 | 256/256 |

The real-step column is only context: 42 of 112 n = 5 steps are not Hom-tilting, and those are the J != 0 or gate-refused ones (E-161: agrees with `tiltingPlus` on 1 081 tests); it is not a NO control for equivalence. P3 is weak as a control (different vertices usually give different labelled Cartans), and the child Cartan is the repo's mutation output, not independent. P1's truth is by construction (P_h and P_h[-1] are both summands, the identity gives a nonzero Hom in degree -1/+1).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/113/skeptic_homctl.py 5      # 3 s
timeout 10m .venv/bin/python workshop/rounds/113/skeptic_homctl.py 6 12    # 3 s
```

## Prior record

E-161/E-163/E-212 give the Hom test and the chain verifier; E-212 already says its NO power is thin (1 Hom rejection in 195 corruptions, no key-equal inequivalent control at n = 7). E-295: Cartan route cannot decide. New here: the pre-registered rule, the point that a Hom failure at the step is expected and not NO, and the planted controls incl. P2 as a measured false-YES rate for Hom-only (the K0 check of E-287 rejects it). K0-unimodularity is already in the E-287 tilting test (generation unverified); P2 only measures what a Hom-only test would miss. Generation assumed: E-168. Not a new retraction.

## Code changed

None outside the new file `workshop/rounds/113/skeptic_homctl.py` (revised: positive control added, dead lines removed; imports `rounds/050/skeptic_tilt.py` unchanged).

## Next

Toolsmith: add the K0 check to the chain verifier and run P1-P3 through it before reading the 13. Experimentalist: a key-equal inequivalent pair at n <= 9 (any) so (4) can be met; until then the answer for T10 at n = 7 is "13 of 13 joined" (YES) or "open", never NO.
