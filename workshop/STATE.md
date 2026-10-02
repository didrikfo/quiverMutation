# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 22
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

## Round 022 updates (supersede everything below where they differ)

Round 022 (ordinary: experimentalist, theorist, toolsmith) recorded E-100..E-102. **Next round (023) is ordinary; 024 is a conference.** The round-020 agenda stands (items 1 and 3 worked once more).

- **T5 (E-100):** among tilting steps on guarded walks (n = 5..7, 479 761 steps) none has the long-sided square; all 2 104 rejecting (parent, v) (n = 6, 7 class 0) do; strict A5 is in 0.5-1.3 % of tilting steps. Open: is "reject iff long square" a theorem (step 7 / E-066 commutativity element); a long-square tilting step off the walks (E-078 family); is `alg.rels` = `relationsFrom(alg)`; classes 2-3 at n = 6, 7; the two n = 8 c2 rejections.
- **T5 (E-102):** the 10 E-084 parents all keep the key under the fixed library and pass the opt-in Cartan check (`checkCartan` / `QM_CHECK_CARTAN=1`, ~15 ms per step); the pre-E-089 reduction fails 10 of 10. Open: whole-walk overhead; E-094's +7 distinct algebras.
- **T6 (E-101):** "cord" = cycle member (arrows >= n, no parallel arrows). Depth-1 rule D1 (some relation of >= 3 arrows not blocked) fits n = 6..10 and an n = 11 sample (fitted among 4 variants); peeling depth 1 + min(a, b) for single-big-relation LNAs; E-099's "within 3 steps" is true for n <= 9 only (`22230222` at n = 10 needs 4). Open: derive the blocking step and the peeling step; LNAs with two big relations (no formula); max depth over all LNAs for n >= 11; the quoted "within 3" in H-017 context.
- Requests added: theorist (to skeptic/scholar): derive "reject iff long square" from step 7. skeptic: a long-square tilting step off the walks. theorist: blocking-rule derivation, two-big-relation LNAs. experimentalist: n = 6 c2-3 and n = 7 c1+ with `--maxexp`.

## Round 021 updates (supersede everything below where they differ)

Round 021 (ordinary: scholar, skeptic, maverick) recorded E-097..E-099. **Next round (022) is ordinary; 024 is the next conference.** The round-020 agenda stands.

- **T5 (E-097):** rewrite entries derived: (k,i) = dim coker g_i (steps 4, 6), (i,k) = dim ker psi_i (step 7, conditional on step-7 completeness); congruence fails iff some dim ker g_i != 0. Rejecting parents are long-sided squares (strict A5 only 767 of 1 123 at n = 6); "A5-shaped" in E-084/E-095 reads so. Open: control (how many tilting parents have the long square); justify e_iBe_{t alpha} = e_iAe_{t alpha}; test step-7 completeness before final reduction; replay of the 10 E-084 parents and the Cartan assertion (toolsmith, carried).
- **T1/T2 (E-098):** "exactly one placement in the `444` orbit" is not enriched (below binomial) and holds for no-4 four-letter words as often or more (6/13/22 vs 6/9/12 at n = 12/13/14); g <= 1 beats a position null but is a property of the orbit's right-end shapes. Open: n = 15..17 (seconds per n); orbit closure of S at each n; do the no-4 exactly-one words reduce by lemma R into S.
- **T6 (E-099):** at n = 8, cord within depth 3 iff a relation has >= 3 arrows (365 / 429); the 64 rad^2-zero LNAs have none to depth 5 (E-087: depth-6 controls exist, so not "cordless"); no monomial cord at depth 3. Open: L = 4 for LNA 205 (first cord at depth 3); n = 9; MONO positive control; proof of the criterion; no overnight run adopted.
- Requests added: experimentalist: long-square control over tilting parents; skeptic: n = 15..17 for the null with |S| closure; theorist: prove or refute the cord criterion, shuttle row (carried); toolsmith: replay of the 10 E-084 parents and the Cartan assertion (carried).

## Round 020 (conference) -- proposed agenda (proposed, round 020; supersedes the round-016 agenda; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/020/proceedings.md`. Round 021 is ordinary; the next conference is 024. T1-T9 and the round updates below stay current.
1. **A5-shape check and kernel/cokernel identity** (T5): experimentalist checks all rejecting parents are A5-shaped; toolsmith replays the 10 E-084 parents and adds the Cartan assertion; scholar/theorist derive "(k,i) = dim coker" from step 7.
2. **Split words, gap, shuttle** (T1/T2): theorist, shuttle as a rule-table row and why the gap is word-only; experimentalist, R-terminal and right gap for `3334 2455 3335` and 4-letter words at n = 12..17; skeptic, selectivity null.
3. **Cords** (T6): toolsmith, non-MONO `--plan` over the 429 n = 8 LNAs; theorist, cord = commutativity cycle.
4. **`k(33x) = 2x`** (T2/T4): theorist, a mechanism.

## Round 019 updates (supersede everything below where they differ)

Round 019 (ordinary: toolsmith, theorist, experimentalist) recorded E-094..E-096. **Next round (020) is a conference.** The round-016 agenda stands.

- **T5 (E-094):** `rounds/019/toolsmith_walk.py` (checkpoint/resume walker). n = 8 class 2 depth 8 completes under the fix: 24 316 expansions, 2 rejections, 0 key-moved steps; E-084's 10 key-moved steps equal in count to the steps that now keep the key (not replayed). Open: replay the 10 E-084 parents; inspect the second rejection (path (17,8,5,6,8,8,2,5)); depth 9 (checkpoint 38 MB, outside repo); n = 8 slice-vs-uninterrupted check; Cartan-congruence assertion in `mutateAtVertex` (toolsmith, carried).
- **T5 (E-095):** histogram of dim ker on walks: exactly one bad vertex per rejecting parent, dim ker 1 (907 + 143 distinct parents at n = 6, 7); no non-tilting step has dim ker 0. Conditional on A5-shaped parents (unchecked). Open: A5-shape check of the 1 050 parents; walks of other classes; derive "(k,i) entry = dim coker" from step 7.
- **T1/T2 (E-096):** in-orbit placement of split words has right gap g = 0 (`2224 2334 4556..4889`) or 1 (`224x`, `344x`), n-independent, `3344` the exception; lemma R alone never reaches `333@0` (0 of 45), it ends at boundary shapes in the orbit and the rest runs through the `34@k <-> 403@(k-1)` shuttle. Open: prove why short shapes are in the orbit only at the boundary; `3344` (R valid only when interval does not swallow a left relation); n = 15, 17 R-chains; `3334`, `2455`.
- Requests added: theorist: derive the shuttle as a rule-table row; `3344`. experimentalist: A5-shape check of rejecting parents. toolsmith: replay of E-084's 10 `M` parents; Cartan assertion; non-MONO cord `--plan` over n = 8 LNAs (carried).

## Round 018 updates (supersede everything below where they differ)

Round 018 (ordinary: skeptic, scholar, maverick) recorded E-091..E-093. Next round (019) is ordinary; 020 is the next conference. The round-016 agenda stands.

- **T1/T2 (E-091):** at n = 12..17, 6/9/12/15/18/18 four-letter words (letters <= 9, a 4, >= 4 placements) are split across the `444` orbit, each with exactly one placement in it (last or second-to-last for 17 of 18 at n = 16); `3334`, `2455` outside it. E-086's "0 partial" is vacuous for merged words. Open: why exactly one placement (does it reduce by lemma R to `333@0`?); where the other placements go; words with letters >= 10 or < 4 placements.
- **T5 (E-093):** Cartan congruence = `tiltingPlus` read through the rewrite on all gate-admitted steps at n = 5..7 (807 non-tilting steps, 0 disagreements); the discrepancy is row k, off-diagonal, minus dim ker. Observation, not a theorem. Open: histogram of dim ker and number of distinct parents; derive "rewrite's (k,i) entry = dim coker" from step 7; diagonal/column-k; Cartan congruence as an assertion in `mutateAtVertex` (toolsmith); the n = 8 class 2 depth-8 walk (checkpoint first).
- **T6 (E-092):** 2376 of 2376 n = 8 cord members (depth <= 5, six LNAs) have a sum relation; no monomial cord at n = 4, 5, 8; no positive MONO control at n = 8 (n <= 5 negatives informative for 6/14 and 1/5 LNAs; depth shallower than E-087). The E-076 candidates are monomial with cords, the controls have sum cords: they differ in kind. Open: an n = 8 MONO run at L = 6 on an LNA outside the six; non-MONO `--plan` over all 429 n = 8 LNAs to find which have cords; theorist: cord = commutativity cycle as a lemma; a derived-equivalent monomial cord outside every LNA class?
- Requests added: theorist: one placement of a split word in the `444` orbit; the (k,i) = coker step in step 7; cord = commutativity cycle. toolsmith: Cartan-congruence assertion; checkpoint for `scholar_walk.py`; non-MONO cord `--plan` over n = 8 LNAs. scholar/experimentalist: dim ker histogram and distinct parents.

## Round 017 updates (supersede everything below where they differ)

Round 017 (ordinary: toolsmith, experimentalist, theorist) recorded E-088..E-090. Next round (018) is ordinary; 020 is the next conference. The round-016 agenda stands (items 1-3 worked once).

- **T5 (E-089, E-090):** `reduceAgainstPivots` is fixed in the library with a test. The E-084 walks re-run under the fix give the same counts wherever comparable (n = 7, 8 classes 0-1, 9); the n = 8 class 2 walk stopped at 17 058 of 20 899 expansions, before the 10 key-moved steps, so the loose end rests on E-085's replay. Open: complete depth 8 of n = 8 class 2 (about 15 min; needs a checkpoint/resume in `scholar_walk.py`); Cartan congruence vs `tiltingPlus` as one criterion; unit test on a sum-relation ideal with a non-pivot head.
- **T6 (E-089):** `MONO=1` at n = 8, L = 6 gives 0 monomial cord members for six walked LNAs (4, 9-13); no positive `MONO` control, LNAs 16..428 not walked. Heuristic: a cord needs a sum relation (untested). Open: positive `MONO` control; L = 5 `MONO` `--plan` over LNAs 16-428 in shards; log the producing relation of each sum cord.
- **T1/T2/T4 (E-088):** lemma R: `(a,b,b,d) -> (a-1,b,d+1)`; `3334`, `2455` are one double mutation from `35`, in class `J = {1, n-7}`, not the `444` orbit (`J = {0, n-6}`); class of `3x` is `x - 4`; n = 12..17. Lemma checked, not proved. Open: lemma R as a rule-table row and a proof; the anchored `3x@0 -> 33(x-1)@0` rule; 4-letter words with a 4 against "apply R, read `x - 4`" at n = 13, 15; x >= 10, letters >= 6; the shift by 2 (T3).
- Requests added: toolsmith: checkpoint/resume for `scholar_walk.py`; `MONO` positive control. experimentalist: 4-letter-word check of `J` (n = 13, 15). theorist: rule-table row for lemma R and the `3x` link; shift by 2. skeptic: n = 16 4-letter row sets (carried).

## Round 016 (conference) -- proposed agenda (proposed, round 016; supersedes the round-012 and round-008 agendas; ordinary rounds work from it until the human changes it in STEERING.md)

Details: `rounds/016/proceedings.md`. Next conference: 020. Threads T1-T9 and the round updates below stay current. Round 017 is ordinary.
1. **Library fix and re-run** (T5; STEERING decision for q 015.1): toolsmith patches `reduceAgainstPivots` with a unit test on the congruent pair (`rounds/015/theorist_fix.py`); experimentalist re-runs the E-084 n = 8 class 2 walk under the fix and says which E-084 counts change. Scholar/theorist: Cartan congruence vs `tiltingPlus` as one criterion.
2. **Cords at n = 8** (T6): toolsmith sizes `MONO=1` at n = 8 (`--plan` first); is there any monomial cord member, and why do none appear at n = 6, 7? Maverick consumes the answer.
3. **`3334`, `2455` in the small orbit / `k(33x) = 2x`** (T1/T2/T4): theorist, with the skeptic's n = 16, 17 row-set test of the `444` orbit as the data.
4. **Parity shift** (T3): theorist, why the rule table gives shift by 2 and never by 1 (E-080).

## Round 015 updates (supersede everything below where they differ)

Round 015 (ordinary: toolsmith, theorist, skeptic) recorded E-085..E-087. Next round (016) is a **conference**. The round-012 agenda stands until then.

- **T5 (E-085):** the n = 8 class 2 loose end of E-084 is a defect of the mutation rewrite: `arrowPaths.reduceAgainstPivots` is not a normal form, so step 7 of `procedure.mutateAtVertex` can drop a relation (child one dimension too big). `tiltingPlus` is right there; the Cartan congruence fails on all 11 replayed rejecting parents. Library untouched; fixed only by monkeypatch (`rounds/015/theorist_fix.py`). Open: apply the fix with a unit test (toolsmith); re-run the E-084 n = 8 class 2 walk and see whether E-084's counts change (experimentalist); is the rank step of `tiltingPlus` affected; how often the defect fires.
- **T6 (E-087):** n = 8 controls with cords (8-9 arrows, sum relations) are found at depth 6 and not 5, at 5.7e4-6.2e4 nodes. `reachedQuipuAlgebras` keeps only monomial quipu trees, hence E-082's "no cords". No monomial cord member at n = 6, 7. Open: `MONO=1` at n = 8; a monomial control; more members (only 2 run).
- **T1/T2 (E-086):** at n = 12..15 all merged words of the earlier scans have all placements in the `444` orbit's row set or none; E-075's "20 of 25" is 11 of 25. Open: row sets at n = 17 (sizes only); n = 16 4-letter scan; why `3334`, `2455` sit in the small `235/255/455` orbit.
- Requests added: toolsmith: apply the `reduceAgainstPivots` fix with a test; `MONO=1` at n = 8. experimentalist: re-run E-084 n = 8 class 2 walk under the fix. theorist: `3334`, `2455` in the small orbit.



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
| experimentalist | 022 | 021 (skeptic) |
<!-- 012 was a conference: all six wrote position statements -->
<!-- 008 was a conference: all six wrote position statements -->
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 022 | 022 (experimentalist) |
| skeptic | 021 | 022 (theorist) |
| scholar | 021 | 022 (toolsmith) |
| toolsmith | 022 | 014 (maverick) |
| maverick | 021 | 015 (toolsmith) |
