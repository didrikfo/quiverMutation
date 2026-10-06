# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 36
next_round_kind: ordinary

## Open threads

Round 006 (ordinary: toolsmith, theorist, scholar) recorded E-064..E-066; threads T3/T8, T2/T4, T5 are updated below. Round 005 was a conference; its ranked agenda (below, `rounds/005/proceedings.md`) is **proposed, round 005**, pending the human's approval in STEERING.md; ordinary rounds work from it. Priority order: T3 (with T8) first, then T2/T4, then T6's control, then T5 reading; T1 and T9 wait on overnight. Round 004 recorded E-061..E-063. Each line: id · question · suited to · status.

- **T1** · H-021: restated without the mirror clause (H-021': a pairing core has `s = n - k(c)`, `k`, `d` independent of n; covers only cores that pair). The 7 cores of E-059 pair at even n, fail at odd; the 12 cores of E-060 keep `k`, `d` at 15 and 16 and pair at 13, 14, 15; at 16 `46 3355 3445` lose the fit through an unmerged equal-size middle pair (size 20300, as `4056` in E-058), so parity is not the whole story (E-062). Open: is 20300 one shared orbit; the 3 at n = 17, 18; the other 127 cores (census at 12, 14 in Menu 4, overnight); why the outside block folds. · experimentalist, theorist · open
- **T2** · (proposed rank 2, with T4: theorist derives `k(33x) = 2x` from the rule table; experimentalist gives `34x`, `44x`, `45x`) Centre `s(c) = n - k(c)`. Interior blocks: `s` = first + last outside offset (13/13, E-060, E-061); end-touching failures `4045 3556 4506 4556` pick the smaller consistent centre for three and the larger for `4506` (no min/max rule); all-outside `4046 5046 5056` maybe a period-2 translation, not a reflection (unverified). `k(33x) = 2x`, `d = x - 3` for x = 3..7 (E-061). Open: why `k = 2x`; `k(c)` for other prefixes (`34x`, `44x`, `45x`); a null test for `|R| <= 4` fits; T4 may explain the top `x - 3` singleton offsets. · theorist, skeptic · open
- **T3** · (proposed rank 1, with T8: toolsmith commits the whole-catalogue orbit task, experimentalist runs 14..16, skeptic referees; first check: is the 20300 pair of `4056` the same orbit as that of the E-060 cores at 16) `3346` never pairs (proved by key, E-054). `4056` at n = 16: `{1,2}` are two mirror-image orbits, key pairs them (E-058), so "orbit = key class" is false there; compare orbit-plus-mirror against key over the `--max-word 4` catalogue at 14..16. Open: other cores that split; the other four H-020 failures (ledgers not in repo). · experimentalist, skeptic · open
- **T4** · H-020 as a theorem: the rule table acts the same at every interior position (F-051). State and prove. · theorist · open
- **T5** · (proposed rank 4: scholar reads Aihara-Iyama Thm 2.32, arXiv:1009.3370, and CHZ, arXiv:2509.12983, against E-032 step 7) H-015: Ladkani 2.3(c) agrees with the gate everywhere tested, incl. non-monomial parents at n = 5..7 (E-055, E-057); only gate-admitted rejection is E-032 step 7. `isTilting` not promoted (STEERING q3). Open: the audit at n = 9/10 (overnight, Menu 4). · scholar, toolsmith · waiting on overnight
- **T6** · (proposed rank 3: toolsmith builds the positive control first; skeptic: null test for `|R| <= 4` fits) H-017: not refuted at n = 9 to depth 6; the Coxeter polynomial and Euler form do not predict (cords, relations); the Euler signature `pos <= n-2` separates outside-every-quipu-class for n = 8..11 (E-063). Open: n = 12 signature (about 15 min in shards); an LNA in the n = 13 class of `P^(1,0,3,0,1)_(1,1,1,1)`; a positive control for the mutation search; a tree eigenvalue argument for `relations >= cords + 1`; depth 7 at n = 9 (overnight, pending question). · maverick, toolsmith, theorist · open
- **T7** · H-010, overlap reducible only at an end: a proof from step 7 of arXiv:2112.08129. Do not run `probe.py --steps 7`. · theorist · open
- **T8** · Tooling: `batch.py orbits N` exists (round 003; `45`, `344` pinned). Open: run it over the whole catalogue (`--jobs 4`), commit the fit/slide step as a task, orbit-plus-mirror prefilter (not key). · toolsmith · open
- **T9** · H-019, H-013: long runs. Overnight waiting on the human: n = 12 and n = 14 censuses of the 139 cores (Menu 4; `workshop/rounds/002/experimentalist_census.py` and `experimentalist_fit.py`); Ladkani audit n = 9/10. · any · parked

Next round (007) is ordinary; 008 is the next conference. Overnight run added in round 005: H-017 depth 7 at n = 9 (Menu 4). Still waiting on the human: n = 12/14 censuses, n = 14 census, Ladkani audit n = 9/10.

## Round 036 updates (supersede everything below where they differ)

Round 036 was a conference (all six personas; `rounds/036/`). **Round 037 is ordinary; 040 is the next conference.** Ranked agenda **(proposed, round 036)**, replacing the round-032 agenda, pending the human's approval in STEERING.md; ordinary rounds work from it:
1. T5: why d_i = 2 when J_i != 0 on walks (E-129); is d_i >= 3 with J_i != 0 possible; hand-build one. theorist, scholar, skeptic.
2. T5: audit the 28 rows with d_i >= 3 (distinct algebras, orbits, classes, out-degrees); checkpointed walk past the cap. experimentalist, toolsmith.
3. T5: n = 7 closure: are all 44 targets keyed (non-None)? then reverse-direction search. toolsmith, experimentalist.
4. T5: the 5 parallel rows of E-127: orbits, Gamma_i vs parallel multiplicity. skeptic, theorist.
5. S-1: derive the last-letter-3 I1/I2 split (E-125) from H-020; n = 12 after `--plan`. maverick, theorist, toolsmith.
Lower: positive control for E-123's key base rate; cyclic-quiver L2 case.
Weakest claims named: the d_i <= 2 bound (empirical, capped); E-123 key absence (no base rate); "0/44 reached" (keys unchecked); "room to move" (hypothesis).

## Round 035 updates (supersede everything below where they differ)

Round 035 (ordinary: toolsmith, experimentalist, scholar; referees skeptic, theorist, skeptic; scholar accepted, the others minor revision, promoted with qualifications) recorded E-128..E-130. **Round 036 is a conference.** The round-032 agenda stands (proposed, pending the human). Rounds 029-031 updates dropped as superseded (see their proceedings).

- **T5 (E-128):** Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum b y_b = 0} for any A and loopless v, so "silting-not-tilting iff some J_i != 0" holds without acyclicity (given AI 2.31); cyclic quiver changes only the count (v -> t -> v, rad^2 = 0: sum J = 1, Hom(T,T[-1]) = 2). The repo gate sees simple paths only, so E-126's L1 is unguaranteed on cyclic quivers (not investigated).
- **T5 (E-129):** capped walks (n = 8 c0, c1, c2; n = 9 c0): 285 rows with J_i != 0 all have d_i = 2, dim J_i = 1; the 28 rows with d_i >= 3 all have J_i = 0. Thin where d >= 3 lives (BFS level 7-8, cap edge); no distinct-algebra counts. Open: deeper walks with d >= 3; why J_i = 0 at d = 3.
- **T5 (E-130):** n = 7 BFS closure of the both-die classes does not fit one command (class 1 ratio about 2.4, class 3 falling from 2.6 to 1.9); 0 of 44 reached in 240 s per class. No overnight approved. Open: E-124's reverse direction (mutate each target, look for an LNA child); measure memory and None-key drops.
- Requests added: toolsmith: reverse-direction search from the 44 targets. experimentalist: d >= 3 rows deeper (checkpointed, distinct algebras); n = 7 class 3 plan. skeptic: a cyclic gate-admitted case where a non-simple path kills dim J_i bound (L1). theorist: why J_i = 0 at d_i = 3; Gamma_i / parallel multiplicity for the 5 rows (still open).

## Round 034 updates (supersede everything below where they differ)

Round 034 (ordinary: theorist, skeptic, maverick; referees experimentalist, scholar, skeptic; all minor revision, all accepted with qualifications) recorded E-125..E-127. **Round 035 is ordinary; 036 is the next conference.** The round-032 agenda stands (proposed, pending the human).

- **T5 (E-126):** on a gate-admitted v, dim J_i <= d_i - 1 (the gate tests single paths), so J_i != 0 needs d_i >= 2; with d_i = 2, dim J_i = 1. For acyclic Q, Hom(N,N[-1]) = 0 and Hom(T,T[-1]) = sum J_i (conditional on T silting, AI 2.31). A layered algebra has d = 3, dim J = 2, so E-124's "dim J_i = 1 on walks" still needs "d_i <= 2 at J_i != 0 on walks" (empirical; E-116 allows up to 8). Open: that bound on walks; a cyclic-quiver test of L2; compare with the AI text.
- **T5 (E-127):** the 5 out-degree 3/4 parallel rows with J != 0 at n = 8 c0 are real Cartan failures; defect support = J support; the gate is blind at out-degree >= 3 as at 2; E-111's "0 rejects at out-degree >= 3" holds for its prefix only. Open: are the 5 rows 3 orbits (mirror pairs); discrepancy values; other classes and n = 9; why the number of J_i != 0 equals the parallel multiplicity; Gamma_i with three-term relations.
- **S-1 (E-125):** at n = 11 in the failing class, the image is a function of (K, core word); same-core words go I1 at K = 3, I2 at K = 4. "Room to move" (image free run 2 against 3) is a hypothesis; the tail table carries the claim. Open: head orientation via `mirrorRow`; K = 2 failing classes at n = 11; n = 12; derive from the H-020 rule table (why a last letter 3 gives I1).
- Requests added: theorist: Gamma_i / parallel multiplicity for the 5 rows; the last-letter-3 rule. toolsmith: out-degree >= 3 rejects on c1 and n = 9 c0 with a checkpointed walk; size the n = 7 BFS. experimentalist: d_i histogram at gate-admitted (v, i) on long walks, any d >= 3 with J != 0 in a derived class. skeptic: a cyclic-quiver case for L2; a hand-built dim J_i >= 2 row. maverick: redo head ends with `mirrorRow`.

## Round 033 updates (supersede everything below where they differ)

Round 033 (ordinary: scholar, toolsmith, experimentalist; referees skeptic, skeptic, theorist; all minor revision, all accepted with qualifications) recorded E-122..E-124. **Round 034 is ordinary; 036 is the next conference.** The round-032 agenda stands (proposed, pending the human).

- **T5 (E-122):** AI 2.32(b) at v is exactly J_i = Hom(S_v, e_iA) (no monomial hypothesis); E-121's "H^{-1}(cone)" means H^{-1} RHom(cone, A). Open: Hom(N,N[-1]) for "silting-not-tilting iff some J_i != 0"; a case with coker != 0 or dim J >= 2; compare with printed AI text.
- **T5 (E-123):** key test has no base rate: walk parents carry the c0 key by construction; in E-121's layered family 0 of 2 704 have an LNA key, even circuit-free and pure-W. E-121's key absence is not evidence. Open: a positive control (sinks attached elsewhere); reconcile 1 593 vs 900.
- **T5 (E-124):** both-die squares exist at n = 6 (gate-admitted, J != 0) but none of 42 has an LNA key; 48 at n = 7, 1 408 at n = 8; capped n = 7 walks found none of the 44. dim J_i = 1 always, inside a 2-dim e_iAe_v. Open: close the n = 7 BFS or reverse-mutate the 44 toward an LNA (size with `--plan` first); why arrow-side cores never get an LNA key; out-degree 3 parallel rows with J != 0 under `checkCartan=True`; save the n = 8 enumeration output.
- Requests added: toolsmith: size the n = 7 BFS closure. theorist: why dim J_i = 1 inside a 2-dim e_iAe_v; Hom(N,N[-1]). skeptic: mutate the out-degree 3 parallel rows with checkCartan.

## Round 032 updates (supersede everything below where they differ)

Round 032 was a conference (all six personas, position statements; `rounds/032/`). **Round 033 is ordinary; 036 is the next conference.** Ranked agenda **(proposed, round 032)**, replacing the round-028 agenda, pending the human's approval in STEERING.md; ordinary rounds work from it:
1. T5: why no long circuits on walks; invariant separating LNA-derived from circuit-carrying algebras; check Aihara-Iyama 2.31/2.32 socle reading against step 7 on one mutation. theorist, scholar, skeptic.
2. T5: base rate of LNA Coxeter keys on random gate-admitted out-degree 2 walk rows at n = 8 c0. toolsmith, skeptic.
3. T5: why both-die rows first appear at n = 8; hand-built n = 6 both-die square; dim J_i on walk J != 0 rows. experimentalist, theorist.
4. T5: prove the Gamma_i two-edge bound. theorist, toolsmith.
5. S-1: core property fixing the deletion threshold K; K >= 4 at n = 12. maverick, theorist.
Requests: theorist: restate AI 2.31/2.32, invariant, why capped walks lack both-die rows; toolsmith: key base rate, per-(orbit, vertex) table for n = 11; experimentalist: dim J_i, parallel-arrow walk controls; skeptic: refute the socle reading on one mutation.

## Round 028 (conference) -- proposed agenda (proposed, round 028; supersedes the round-024 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/028/proceedings.md`. Round 029 is ordinary; next conference 032. T1-T9 and the round updates below stay current.
1. **No long circuits / converse of W** (T5): theorist + scholar prove or refute "no nn 2-cycle, no circuit >= 3 on LNA-derived algebras"; can D, G, H occur inside a walk? First: the skeptic's 22 n = 8 c0 rows (both J nonzero, empty intersection).
2. **Controls for W** (T5): experimentalist + toolsmith hand-build a parallel-arrow positive control and a "cancels" example; skeptic checks gate admission of the 61 D' rejects.
3. **Resolve the 16 n = 9 and 176 n = 10 cospectral LNAs** (S-1 validation): toolsmith + experimentalist; size with `--plan` first.
4. **S-1: derive K >= 3** (maverick + theorist): why K = 2 fails; relation to H-020; then n = 11, S-1 question 3.
5. Carried, lower rank: mirror chain (T6), `k(33x) = 2x` (T2/T4), n = 15..17 (T1/T2).

## Round 027 updates (supersede everything below where they differ)

Round 027 (ordinary: scholar, experimentalist, maverick) recorded E-110..E-112, all accepted with caveats. **Round 028 is a conference** (multiple of 4). The round-024 agenda stands.

- **T5 (E-110):** "reject => W" is not derivable from step 7 (J is defined by the parent alone). For monomial + two-term relations (scalar 1) J != 0 iff the circuit graph Gamma_i has a circuit; W is its length-2 ground-path case. Hand-built D (cancels / nn), G (out-degree 1, no long square), H (3-term) are gate-admitted rejects W misses. On walks J != 0 is always a length-2 ground path. Open: why no nn or longer circuit on LNA-derived algebras; balanced circuits for scalars != 1; >= 3-term relations.
- **T5 (E-111):** W has 0 mismatches on 32 132 fresh out-degree 2 rows (n = 8 c0, 2, 3; n = 9 c0 prefix; 1 135 parallel) but only 61 positives (n = 8 c0). Out-degree >= 3: 0 rejects. Open: parallel-arrow positive control; classes 4..10 at n = 8; n = 9 c0 beyond the prefix (no overnight adopted).
- **S-1 (E-112):** no positional rule for deleting a vertex (same-class pairs kept 0.43 vs 0.19 chance at n = 9); deleting a free end vertex with K >= 3 free vertices there gives an image class depending only on the source class at n = 8, 9, 10 (3, 5, 10 resolved classes). Open: derive it (F-028/H-020?); '?' ends at n = 10; independent validation of the class labels; n = 11; the 16 unresolved n = 9 LNAs; question 3 of S-1 (core with room to move).
- Requests added: theorist: why no nn/longer circuit on LNA-derived algebras; the balanced-circuit lemma; why K >= 3. toolsmith/experimentalist: resolve the 16 n = 9 and 176 n = 10 cospectral LNAs; hand-built parallel and "cancels" controls for W. experimentalist: n = 11 for the K >= 3 rule.

## Round 024 (conference) -- proposed agenda (proposed, round 024; supersedes the round-020 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/024/proceedings.md`. Round 025 is ordinary; next conference 028. T1-T9 and the round updates below stay current.
1. **The 42 n = 8 class-0 rejects** (T5): scholar classifies them (revision due); skeptic: reachable by key-preserving walks?; toolsmith: Cartan test through the rewrite; experimentalist: capped walks at n = 8 classes 1, 3 and n = 9 class 0.
2. **Mirror chain** (T6): maverick sorts the 24 failing LNAs by big blocker; experimentalist n = 11/12 deeper shapes; theorist a derivation.
3. **Split words / shuttle / `k(33x) = 2x`** (T1/T2/T4): experimentalist n = 15..17; theorist shuttle row.
4. **Cord = commutativity cycle** (T6): theorist, toolsmith.

## Round 020 (conference) -- proposed agenda (proposed, round 020; supersedes the round-016 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/020/proceedings.md`. Round 021 is ordinary; the next conference is 024. T1-T9 and the round updates below stay current.
1. **A5-shape check and kernel/cokernel identity** (T5): experimentalist checks all rejecting parents are A5-shaped; toolsmith replays the 10 E-084 parents and adds the Cartan assertion; scholar/theorist derive "(k,i) = dim coker" from step 7.
2. **Split words, gap, shuttle** (T1/T2): theorist, shuttle as a rule-table row and why the gap is word-only; experimentalist, R-terminal and right gap for `3334 2455 3335` and 4-letter words at n = 12..17; skeptic, selectivity null.
3. **Cords** (T6): toolsmith, non-MONO `--plan` over the 429 n = 8 LNAs; theorist, cord = commutativity cycle.
4. **`k(33x) = 2x`** (T2/T4): theorist, a mechanism.



(Round 023 updates dropped as superseded; see `rounds/023/proceedings.md`. Round 015-017 updates and the round-016 agenda are dropped from this file as superseded; see `rounds/015..017/proceedings.md` and E-085..E-090.)

## Round 006 updates (supersede the thread text above where they differ)

- **T3/T8 done at the catalogue level (E-064):** orbit-plus-mirror refines the key in all 139 cores at n = 10, 12..16; equals it in 129-132; the rest are key-coarser cores where the key merges the two parity classes (9 at even n, 10 at odd n). The 20300 pairs of `4056` and the 3 cores of E-062 are one orbit and its mirror, so the "unmerged middle pair" is resolved by the mirror join. Open: why those 9/10 cores are parity classes, and whether they are the 7 of E-059; `--max-word 5` at n = 14 (needs `--plan`); a prefilter by key difference (true in 139 x 6 cases, not proved); `toolsmith_orbitclass.py` has no `--budget-hours`.
- **T2/T4 (E-065):** `k(33x) = 2x`, `d = x - 3` from the drift of the double mutation plus seed `333` (lower bound derived; end link and upper bound computed only). Fails for `44x` (one orbit with `333@0`), `34x` only `345` explained, `45x` no drift. Open: the full-orbit size comparison for 33x; the end link as a proof; n = 17 for 34x/45x; n = 18 rule-table claim has no output file; x = 9 excluded by a script cap.
- **T5 (E-066):** step 7 is rejected at a commutativity element; AI 2.32(b), Ladkani 2.3(c), `tiltingPlus` appear to be one map (no derivation); CHZ Cor 3.6 path-wise wording may need "monomial" (UNVERIFIED, PDF unread, arXiv blocked by the proxy this round). Open: derive the one-map identity; build a smaller instance (commutative square into a vertex with one outgoing arrow) at n = 5..7 and run `tiltingPlus`; compute End of the two-term complex directly for an independent check.
- **T7 (H-010):** step 7 is a concrete non-monomial case to test the "overlap reducible only at an end" proof on.

Requests added: toolsmith to theorist: why the key-coarser cores are the parity classes. scholar to theorist: the n = 5..7 instance for T5. scholar/toolsmith: read the CHZ PDF where the network allows.

## Awaiting revision

(none)

## Requests between personas

- theorist (to experimentalist): `k(c)` for `34x`, `44x`, `45x` at n = 14..17; the 4046-type cores at 14 (period-2 translation?).
- maverick (to toolsmith): a positive control for the H-017 mutation search.
- maverick (to experimentalist): n = 12 signature census in shards.
- experimentalist (to theorist): what separates unmerged `{4}{5}` of `3355` at 16 from merged `{3,4}` at 14.
- experimentalist: orbit-plus-mirror vs key classes over the catalogue at 14-16 (skeptic).
- theorist: at n = 16 `4056` has a self-dual pair `{0,3}` and a mirror pair `{1,2}`: is the split a function of parity or of position relative to the middle? (skeptic)
- theorist: the centre `s(c)` formula (round 001 request, still open).
- toolsmith: committed census task; key prefilter with the E-058 caveat.

## Rota

<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 035 | 035 (theorist) |
<!-- 036 was a conference: all six wrote position statements -->
<!-- 012 was a conference: all six wrote position statements -->
<!-- 008 was a conference: all six wrote position statements -->
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 034 | 035 (experimentalist) |
| skeptic | 034 | 035 (scholar, toolsmith) |
| scholar | 035 | 034 (skeptic) |
| toolsmith | 035 | 027 (maverick) |
| maverick | 034 | 015 (toolsmith) |
