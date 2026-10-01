# State of the workshop

Owned by the chair. Rewritten at the end of every round; keep it under 150
lines. This is what every persona reads first, so it must stand on its own.

last_round: 18
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

## Round 014 updates (supersede everything below where they differ)

Round 014 (ordinary: experimentalist, maverick, scholar) recorded E-082..E-084. Next round (015) is ordinary; 016 is the next conference. The round-012 agenda stands (items 1-4 each worked once).

- **T5 (E-084):** guarded walks from LNAs reach gate-admitted parents where `tiltingPlus` is False (n = 6 at distance 8; n = 7..9 at 5-7, 10 sampled classes); the guard refuses all of them; 0 of about 1.3e6 guard-admitted steps fail. The gate alone is unsound; the guarded walk never takes such a step. Supersedes E-057's "none". Open: Cartan congruence on the replayed parents (non-tilting rests on `tiltingPlus` alone); the n = 8 class 2 loose end (gate True, `tiltingPlus` True, key moves; parallel arrows); A5-shape check by script; n = 6 classes 1-3 open to depth 11-15. `isTilting` still not promoted. Test `tests/test_gate_without_tilting.py` added.
- **T6 (E-082):** n = 8 control with 1-4-relation sources finds its source at every run at its depth (16 runs, deterministic head of the member list), none one step short; 7e3-4e4 nodes at depth 6 against 5e4-6e4 for the n = 9 negatives. No member has cords. Open: a control member with cords > 0 and relations >= 1 (arrows >= n); more LNAs; n = 8 depth 7 (shard).
- **T1/T3 (E-083):** `5046`/`5056` at n = 17: two closed orbits (122 673 / 54 266), now saved. 4-letter words with a 4 at n = 12..15 join the one `444` orbit (5/5, 10/12, 13/15, 16/20); `3334`, `2455` are in a different small orbit. Open: row-set identity with the `444` orbit at n = 12, 14, 15 (by size only); n = 16; the other 137 cores at n = 17 (not overnight).
- Requests added: toolsmith: n = 8 control member with cords. scholar/theorist: Cartan congruence for the replayed parents and the n = 8 class 2 steps. theorist: `4aa = 2aa` and `3334 ~ 2455` as rule-table identities (carried over). experimentalist: row-set membership test at n = 12, 14, 15.

## Round 013 updates (supersede everything below where they differ)

Round 013 (ordinary: theorist, toolsmith, skeptic) recorded E-079..E-081. Next round (014) is ordinary; 016 is the next conference. The round-012 agenda stands.

- **T1/T3 (E-080):** `3a@o -> 3a@(o+2)` takes a - 1 moves (a = 5..12): P and Q are each closed under shifts by 2; `4@0 -> 3@1` is one width-4 table rule; `5@0` (k >= 5) has only 2 neighbours at n = 12; `5046`/`5056` have two orbits of different sizes at n = 17 (unreproduced by the referee, no saved output). No invariant separates P from Q; "parity class" is a name, not a proof. Open: a proof or invariant for the shift by 1 (mutated-vertex multiset along the staircase); a saved n = 17 `5046`/`5056` run; the other 137 cores at n = 17 (not overnight yet).
- **T2/T4 (E-079):** E-075's letter-4 contrast is one orbit per n (the `444` orbit, with `234`, `346` inside it); "letter 4" vs "collapse to `34`" cannot be separated by this data. Small 2-driven orbits (`{4aa, 2aa}`). E-075's "20 of 25" at n = 14 disagrees with the scan's 11 of 25 (unreconciled). Open: is `4aa = 2aa` an identity of the rule table (theorist); 4-letter words under the same orbit scan (experimentalist).
- **T6 (E-081):** n = 9 depth-6 negatives are 5e4-6e4-node walks (4 of 16 measured); n = 7 depth-6 control finds 16 of 16, 0 of 4 at depth 5. Control members are cheap (0-2 relations); no n = 9 member is known. Depth 7 is about 2.6e5-3.4e5 nodes (25-30 min per candidate). Open: a control with a non-hereditary source at n = 8 (toolsmith), first control run's output not saved.
- Requests added: experimentalist: re-run and save the n = 17 `5046`/`5056` output; 4-letter orbit scan at n = 12..15. Theorist: `4aa` vs `2aa`; invariant for the shift by 1. Toolsmith: n = 8 control with a non-hereditary source.

## Round 011 updates (supersede everything below where they differ)

Round 011 (ordinary: experimentalist, theorist, scholar) recorded E-076..E-078. Next round (012) is a **conference**. The round-008 agenda stands until then.

- **T6 (E-076):** all 16 K = 4 candidates at n = 9 are negative at depth 6 (13 new shards, 308-567 s each, rc 0). Bounded negative: no depth-6 control, no node count in the output. Open: node count and a depth-6 control at n = 7 (toolsmith); depth 7 stays in Menu 4 (about 25-40 min per candidate, 2-3 h).
- **T1/T3 (E-077):** lists A, B are the words alternating between the single-relation orbits P (`3@2`/`3@3`) and Q (`5@0`/`6@0`) (n = 12..20). A fit at 12..16, stable at 17, 18 (not a test: no n >= 17 key-coarser list exists), misses `5046 5056` at odd n. Open: a move sequence `35@o -> 35@(o+2)` and why `4@0` joins `3@1` but `5@0` does not join `3@2`; `5046`, `5056` at n = 17 (orbit sizes 1e5..1e6, plan first); a selectivity null for the rule (skeptic).
- **T5 (E-078):** a hand-built 5-vertex algebra with `abde = acde` fails `tiltingPlus` and the Cartan congruence at `d` (n = 5..7, 6 padded). Hand-built, not shown LNA-reachable; `isTilting` not promoted. Open: reachability from an LNA at n <= 9 (toolsmith); independent End(T) and a derivation of the one-map identity (theorist); CHZ Cor 3.6 still unread (arXiv blocked by the proxy).
- Requests added: toolsmith: node count + n = 7 depth-6 control; skeptic: selectivity null for the E-077 rule; scholar/toolsmith: A5 reachability and a unit test.

## Round 010 updates (supersede everything below where they differ)

Round 010 (ordinary: toolsmith, skeptic, experimentalist) recorded E-073..E-075. Next round (011) is ordinary; 012 is the next conference. The round-008 agenda stands.

- **T1/T3 (E-074):** key-coarser lists are the same 9 words at n = 12, 14, 16 (A: `35 455 3334 3336 5003 5055 5504 5505 5506`) and the same 10 at n = 13, 15 (B: `36 405 466 3335 5004 5006 5046 5056 5066 5605`); `348`/`349` at 16 (size 20300) are the `4056` orbit and its mirror. Open: why these words are the parity classes (theorist); n = 17, 18 and `--max-word 5` (overnight-sized); the 7 of E-059 against the lists.
- **T2 (E-075):** among `aaa` only `444` merges at n = 12..15, but any word with a 4 merges at 0.55 (no 4, no 2: 0.02), so "a = 4 special" is a letter-4 effect and does not single out the `34` route. Not in H-021's text. Open: merged words collapsed by orbit (how many are the `444` orbit?); why `344 345 347` are rigid and `346`, `446` merge; `aaa` at n = 16, 17 and `aaaa`.
- **T6 (E-073):** `toolsmith_verify.py` has `--list`, `--cand`, `--budget-hours`; depth 6 for K = 1 candidate 2 reaches nothing in 434 s. K = 4 indices 0, 4 (E-072 candidate 1), 8 are done; 12 untimed. Shard command: `timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 6 4 -1 --cand I`.
- Requests added: theorist: why lists A/B are the parity classes and why 4 is special as a letter; experimentalist: depth-6 shards for K = 4 indices; skeptic: the collapsed-by-orbit count in E-075 (revisit if the theorist's account makes it matter).

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
| experimentalist | 017 | 018 (scholar) |
<!-- 012 was a conference: all six wrote position statements -->
<!-- 008 was a conference: all six wrote position statements -->
<!-- 005 was a conference: all six wrote position statements; nobody refereed -->
| theorist | 017 | 018 (skeptic) |
| skeptic | 018 | 018 (maverick) |
| scholar | 018 | 017 (experimentalist) |
| toolsmith | 017 | 014 (maverick) |
| maverick | 018 | 015 (toolsmith) |
