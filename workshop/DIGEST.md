# Digest

What happened, one entry per round, newest first. Written by the chair for the
human coming back after a while: what was claimed, what survived review, what
was promoted into `research/`, what the chair needs from you. Each entry at
most `max_digest_entry_lines` lines, linking to `rounds/NNN/` for the rest.

---

## Round 015 -- 2026-10-01 -- ordinary

Worked: toolsmith (T6), theorist (T5), skeptic (T1/T2). Referees: maverick, skeptic, experimentalist. Three minor revisions; I applied the scope fixes in the entries and accepted all three. Details in `rounds/015/`.

- **theorist:** the n = 8 class 2 loose end of E-084 is a bug in the mutation rewrite, not in `tiltingPlus`: `reduceAgainstPivots` is not a normal form, so step 7 can drop a relation. The Cartan congruence fails on all 11 replayed rejecting parents, agreeing with `tiltingPlus`. Referee reproduced a congruent pair with different residues and the patched runs; the library is untouched, and E-084's counts were made with the buggy rewrite and are not re-run.
- **toolsmith:** n = 8 controls with cords (8-9 arrows) are found at depth 6 and not 5, at the n = 9 negatives' size (5.7e4-6.2e4 nodes). Caveat: every cord member has a sum relation, the n = 9 candidates are monomial, and none was found at n = 6, 7; 2 members only.
- **skeptic:** at n = 12..15 every merged word of the earlier scans lies in the `444` orbit's row set or outside it (0 partial); E-075's "20 of 25" is really 11 of 25.

Promoted: E-085, E-086, E-087; H-015, H-017 status lines; E-075, E-084 annotated.

**Questions for you (the chair takes the recommended option if unanswered):** (1) the `reduceAgainstPivots` fix: recommend toolsmith patches it with a unit test after the conference, and the E-084 n = 8 class 2 walk is re-run. (2) Overnight: recommend none. Decided for you (round 014): `isTilting` still not promoted; no overnight. Round 016 is a conference.

---

## Round 014 -- 2026-10-01 -- ordinary

Worked: experimentalist (T1/T3), maverick (T6), scholar (T5). Referees: theorist, toolsmith, skeptic. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/014/`.

- **scholar:** yes, a guarded walk from an LNA reaches a gate-admitted parent where `tiltingPlus` fails: n = 6 at distance 8, n = 7..9 at 5-7 (10 sampled classes). The Coxeter guard refuses every one, and 0 of about 1.3e6 guard-admitted steps fail. So the gate alone is unsound; the guarded walk is not. Referee reproduced the n = 6 and n = 7 cases and the new test; non-tilting rests on `tiltingPlus` alone and an n = 8 loose end (parallel arrows) is open. E-057's "none" was depth-limited and is annotated.
- **maverick:** an n = 8 control with relation-bearing sources finds its source in every run at its depth and none one step short, at 7e3-4e4 nodes at depth 6 (n = 9 negatives: 5e4-6e4). Members are the head of a sorted list and have no cords, so it controls the walk, not coverage of classes with cords.
- **experimentalist:** `5046`/`5056` at n = 17 have two closed orbits (122 673 / 54 266), now saved (confirms E-080). 4-letter words with a 4 at n = 12..15 join the single `444` orbit (identity by row membership at n = 13 only, by size elsewhere).

Promoted: E-082, E-083, E-084; status lines of H-015, H-017, H-021; E-057 annotated.

**Questions for you (the chair takes the recommended option if unanswered):** (1) `isTilting`: recommend still not promoting; Cartan check on the replayed parents and the n = 8 loose end first. (2) Overnight: recommend none; toolsmith builds an n = 8 control with cords first. Decided for you (round 013): no n = 17 overnight, no H-017 depth 7.

---

## Round 013 -- 2026-10-01 -- ordinary

Worked: theorist (T1/T3), toolsmith (T6), skeptic (T2/T4). Referees: skeptic, theorist, experimentalist. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/013/`.

- **theorist:** the two orbits behind lists A/B are each closed under offset shifts by 2, by an explicit staircase of `a - 1` double mutations (`3a@o -> 3a@(o+2)`, a = 5..12; referee extended it to 10..12); `4@0 -> 3@1` is one table rule, `5@0` has only 2 neighbours. No invariant explains why a shift by 1 is impossible: "parity class" names two orbits. `5046`/`5056` have two orbits of different sizes at n = 17 (the referee later reproduced the `5046` run; `5056` not re-run; no output saved).
- **toolsmith:** the n = 9 depth-6 negatives are walks of 5e4-6e4 nodes (4 of 16 measured), and an n = 7 depth-6 control finds its target 16 of 16 (0 of 4 at depth 5). The control members are cheap, and no n = 9 member is known, so it shows the search ran, not that it was enough. The depth-7 cost is 25-30 min per candidate.
- **skeptic:** counted by orbit, E-075's letter-4 contrast is one orbit per n (the `444` orbit), which holds both `34`-words and 4-no-`34` words, so letter 4 and collapse to `34` are not separable. Referee re-ran identically; E-075's "20 of 25" at n = 14 does not match the scan's 11 of 25 (unreconciled).

Promoted: E-079, E-080, E-081; status lines of H-021 and H-017; glossary "Staircase".

**Questions for you (the chair takes the recommended option if unanswered):** (1) n = 17 key-coarser lists overnight: recommend not yet; save and reproduce the `5046` n = 17 run first. (2) H-017 depth 7 overnight: recommend no; the toolsmith sizes an n = 8 control with a non-hereditary source first. Decided for you (round 012): no n = 17 overnight; round-012 agenda approved unchanged.

---

## Round 012 -- 2026-10-01 -- conference

Six position statements (`rounds/012/`); no new work, nothing promoted. Next round (013) is ordinary; 016 is the next conference.

- Weakest claim named by four personas: E-077's lists A/B are a fit with no mechanism and no real test at n >= 17 (and miss `5046 5056` at odd n). Others: the "walk-reachable only" rule for `isTilting` is policy, not a result (scholar); the H-017 negatives cannot be read without a node count and control (toolsmith, maverick). Maverick's Euler-form claim is unchecked (E-063 covers n = 8..11 only).
- **Proposed agenda (ranked):** 1. why A/B are parity classes: move sequences in P/Q, then `5046 5056` (theorist, experimentalist); 2. H-017 node count + n = 7 depth-6 control (toolsmith); 3. walk-reachable gate-admitted non-tilting mutation (toolsmith, scholar); 4. letter 4 versus the `34` collapse (theorist, skeptic).

**For you: approve or change this agenda in `STEERING.md`; until then rounds follow it.** Question: n = 17 key-coarser lists overnight (2+ h) -- recommend not yet. Decided for you (round 011): hand-built rejection does not count towards `isTilting`; no new H-017 overnight.

---

## Round 011 -- 2026-10-01 -- ordinary

Worked: experimentalist (T6), theorist (T1/T3), scholar (T5). Referees: toolsmith, skeptic, theorist. Three minor revisions; I applied the fixes and accepted all three. Details in `rounds/011/`.

- **experimentalist:** all 16 K = 4 candidates at n = 9 reach nothing at depth 6 (13 new shards, 308-567 s). Referee re-ran one shard identically. A bounded negative: no depth-6 control, no node count.
- **theorist:** the key-coarser lists A, B are the words alternating between two single-relation orbits (`3@2`/`5@0` even n, `3@3`/`6@0` odd n), n = 12..20. Referee: a fit at 12..16 and a consistency check at 17, 18 (no n >= 17 list was computed), misses `5046 5056` at odd n; the single-relation identification itself is new.
- **scholar:** a hand-built 5-vertex algebra with `abde = acde` is gate-admitted at `d` yet fails `tiltingPlus` and the Cartan congruence (n = 5..7), so E-066's step-7 shape is not special to n = 10. Hand-built, not LNA-reachable; CHZ still unread (arXiv blocked).

Promoted: E-076, E-077, E-078; status lines of H-017, H-021, H-015. Round 012 is a conference.

**Questions for you (chair takes the recommended option if unanswered):** (1) does a hand-built gate-admitted rejection count towards promoting `isTilting`: recommend no, it must come from a walk; (2) H-017: recommend no new overnight run, toolsmith adds a node count and depth-6 control first. Decided for you (round 010): no overnight; `orbitclass` at n = 17 not run.

---

## Round 010 -- 2026-10-01 -- ordinary

Worked: toolsmith (T6), skeptic (T2), experimentalist (T1/T3). Referees: theorist, experimentalist, skeptic. Two minor revisions, one accept; I applied the fixes and accepted all three. Details in `rounds/010/`.

- **toolsmith:** `toolsmith_verify.py` now takes `--list`, `--cand` and `--budget-hours`; E-072's depth-5 negative reproduces and a second n = 9 candidate reaches nothing at depth 6 (434 s, fits one shard). Referee: a reproduction plus tooling; wording fixes only.
- **skeptic:** the neighbour-aware null for "a = 4 special": only `444` merges among `aaa` at n = 12..15, but any word with a 4 merges far more often (55/100 against 3/121 with neither 4 nor 2), so the claim is a letter-4 effect and does not single out the collapse to `34`. Referee reproduced it; the orbit-collapsed count was not done.
- **experimentalist:** the key-coarser cores are the same 9 words at n = 12, 14, 16 and the same 10 at n = 13, 15; `348`/`349` size-20300 pairs at 16 are the `4056` orbit and its mirror. Referee re-ran it byte-identically.

Promoted: E-073, E-074, E-075; status lines of H-021, H-017.

**Questions for you (chair takes the recommended option if unanswered):** (1) remaining H-017 depth-6 candidates: recommend chair-slot shards in round 011, no overnight; (2) `orbitclass` at n = 17 (2+ h): recommend not until the theorist explains why those words are parity classes. Decided for you (round 009): shard tooling added, no overnight; `aax` criterion kept out of H-021.

---

## Round 009 -- 2026-09-30 -- ordinary

Worked: experimentalist (T1/T3), theorist (T2/T4), maverick (T6). Referees: skeptic, experimentalist, toolsmith. All minor revision; I applied the fixes and accepted all three. Details in `rounds/009/`.

- **experimentalist:** the equal-size singleton pairs of `344 348 349` (n = 15..17) and `4046` (n = 14..16) are one orbit plus its mirror, orbit+mirror = key in 12/12 cells; the key-coarser cores at n = 12, 13 are not the 7 of E-059. Referee: mostly already in E-064 (only n = 17 new); the "pair at even n, mirror-join at odd n" summary contradicted the table and was withdrawn.
- **theorist:** `aax` drift families close with `k = 2x + 3 - a` for a = 3, 5, 6 (7 by the referee) and `44x` does not, because the seed `444` collapses to `34`, which reaches the slider `44`. Predicted before running for `55x`, `66x`; referee's extra runs agree. Computed criterion, not proved (one positive datum).
- **maverick:** L = 5 control passes (42/42 at n = 6; 8/8 near-trivial LNAs at n = 7); four n = 9 candidates reach nothing at depth 5. Referee corrected the sizing: depth 6 for 16 candidates about 2 h, one candidate per 10-minute shard.

Promoted: E-070, E-071, E-072; status lines of H-021, H-017.

**Questions for you (chair takes the recommended option if unanswered):** (1) H-017 depth 6 for all 16 n = 9 candidates: recommend toolsmith adds a candidate-index argument first, no overnight run (only depth 7 is overnight, already in Menu 4); (2) write the `aax` criterion into H-021's text: recommend not until a skeptic's null. No questions were outstanding from round 008.

---

## Round 008 -- 2026-09-30 -- conference

Six position statements (haiku), no new work, nothing promoted. Details in `rounds/008/`. Their factual sub-claims are unchecked.

- **Converging:** experimentalist, theorist and toolsmith all want the mirror-join / parity-class check of the equal-size singleton pairs (`344`, `348`, `349`, `4046`; the 9-10 key-coarser cores of E-064 vs the 7 of E-059). Skeptic's weakest claim: the centre formula was fitted at n = 13 and "confirmed" on pre-selected cores, so it needs a fresh random sample.
- **Proposed agenda (ranked):** (1) mirror join + parity classes; (2) fresh-sample centre test; (3) why the `33x` drift closes and `44x` does not; (4) H-017 L = 5 control, sizing before any overnight; (5) H-015 one-map identity with a small commutative instance. Maverick's Coxeter-spectrum question is parked.

**Please approve or change the proposed agenda in `STEERING.md`; until you do, ordinary rounds work from it.** No new questions. Decided for you (round 007 questions): no overnight depth 5-6 rerun of H-017 yet; no `34x` at n = 18, 19 until the mirror join is done.

---

## Round 007 -- 2026-09-30 -- ordinary

Worked: experimentalist (T2), skeptic (T2/T6), maverick (T6). Referees: skeptic, theorist, experimentalist. All minor revision; I applied the referees' wording fixes and accepted all three. Details in `rounds/007/`.

- **experimentalist:** `34x` offsets pair `o <-> hi - o` (`k = x + 3`, not `2x`) at n = 14..17 for x = 4, 5, 7, 8, 9 (20/20 cells); `346` one orbit; `45x` no reflection; `4046` is a (size-paired) reflection. E-060's `4046@13` line did not reproduce. Caveat: `k = x + 3` is just `s = hi`, and singleton pairing is by size only.
- **skeptic:** null test. 39 of 109 n = 13 fits are vacuous and informative ones are mostly chance-level (73%), so E-061's "45/62, 10/13" counts are padded; the interior-core centre formula (13/13 vs 4.1 expected) and the E-060 cores at n = 15/16 survive (joint p about 1e-3..1e-4 after the referee's correction).
- **maverick:** the H-017 search passes a positive control (273/273 round trips at n = 7 on 91 of 132 LNAs; 0 at one level too shallow), so the round-004 depth-4 negative only excludes members within 4 steps.

Promoted: E-067, E-068, E-069; status lines of H-021 and H-017. Round 008 is a conference.

**Questions for you (chair takes the recommended option if unanswered):** (1) overnight depth 5-6 rerun of the 16 H-017 candidates at n = 9: recommend not yet, size it first; (2) `34x` at n = 18, 19 overnight: recommend no, mirror-join check first. Decided for you (round 006 questions): agenda approved unchanged; no overnight run.

---

## Round 006 -- 2026-09-30 -- ordinary

Worked: toolsmith (T3/T8), theorist (T2/T4), scholar (T5). Referees: skeptic (x2), experimentalist. Toolsmith accepted; theorist and scholar minor revision, which I applied myself. Details in `rounds/006/`.

- **toolsmith:** over all 139 placed cores at n = 10, 12..16, orbit-plus-mirror classes refine the key classes (never finer, never incomparable; equal in 129-132); the exceptions are the parity-merged cores. The 20300 pairs at n = 16 (`4056`, `46`, `3355`, `3445`) are one orbit and its mirror, so E-058/E-062's "unmerged middle pair" is resolved. Referee reproduced it.
- **theorist:** `k(33x) = 2x`, `d = x - 3` come from a drift `33x@o -> 33(x-1)@(o+1)` of the double mutation (label `x + o` fixed) plus the self-dual seed `333`: 135/135 at n = 14..16, plus the referee's n = 17. Lower bound derived, upper bound computed only; the same argument is false for `44x`.
- **scholar:** E-032 step 7 is rejected at an explicit commutativity element; mostly known, and the literature is not an independent test of the code. CHZ "monomial only" caveat is UNVERIFIED (arXiv blocked).

Promoted: E-064, E-065, E-066; H-021 status line; GLOSSARY (two terms); UNVERIFIED flag on `literature/2509.12983`. No open questions were left by round 005.

**Questions for you (chair takes the recommended option if unanswered):** (1) approve the round-005 agenda (recommend yes, unchanged; T3 is now mostly answered); (2) no overnight run proposed yet, toolsmith to size `--max-word 5` at n = 14 first.

---

## Round 005 -- 2026-09-30 -- conference

Six position statements (haiku), no new work, nothing promoted. Details in `rounds/005/`.

- **Convergence:** experimentalist, skeptic and toolsmith all want the same test: orbit-plus-mirror vs key classes over the `--max-word 4` catalogue at n = 14..16 (T3/T8). Theorist wants a mechanism for `k(33x) = 2x`; maverick a positive control for the H-017 search; scholar to decide E-032 step 7 from the literature (Aihara-Iyama, CHZ).
- **Weakest claims named:** parity (19 of 139 cores tested), `k = 2x` (description, not mechanism), H-015 support (n <= 7), "Coxeter polynomial cannot see (cords, relations)" (n = 9 only). No disagreement between personas.
- **Agenda proposed** (ranked, in `STATE.md` as proposed, round 005): 1 orbit-vs-key at 14..16; 2 `k = 2x` from the rule table (T4); 3 H-017 positive control; 4 H-015 step 7 by reading; 5 parity across the other 127 cores waits on the overnight censuses.
- **Decided for you (round 004 questions):** approved the n = 9 depth-7 H-017 run only (in `OVERNIGHT.md` Menu 4; I added `--budget-hours` to `maverick_reached.py`, tested, `test_overnight_doc` passes); round 005 kept a conference. Overturn in `STEERING.md`.

**Please approve or change the proposed agenda in `STEERING.md`; until then rounds 006+ work from it.**

---

## Round 004 -- 2026-09-30 -- ordinary

Worked: theorist (T2), experimentalist (T1), maverick (T6, first time). Referees: skeptic (x2), scholar. All minor revision; I accepted all three after applying the referees' wording and one check myself. Details in `rounds/004/`.

- **theorist:** the interior/end-touch split is now a committed column and explains none of the 7 failures (the slide never picks a centre, and the orbits take the smaller or larger consistent one at random). For `33x`, `k = 2x` and `d = x - 3` at n = 13..17 for x = 3..6; I added `337` (holds). A description, not a mechanism.
- **experimentalist:** the 12 cores of E-060 keep `k` and `d` at n = 15, 16, and pair at 13, 14, 15; at 16 three lose the fit through an unmerged equal-size middle pair (as `4056` in E-058). So parity is not the whole story.
- **maverick:** H-017 survives to depth 6 at n = 9, but neither the Coxeter polynomial nor the Euler form predicts (cords, relations). The Euler-form signature does separate "outside every quipu class" for n = 8..11 (small, new to the record).

Promoted: E-061, E-062, E-063; status lines of H-021 and H-017. Step 0.5: decided round 003's question (yes) and added the n = 12 and n = 14 censuses to `OVERNIGHT.md`. Round 004 should have been a conference by the config; round 005 will be.

**Questions for you (chair takes the recommended option if unanswered):** (1) H-017 overnight: approve the n = 9 depth-7 run, not the n = 10 run until the search has a positive control? Recommend yes to n = 9 only. (2) Keep round 005 a conference? Recommend yes.

---

## Round 003 -- 2026-09-30 -- ordinary

Worked: theorist (T1/T2), experimentalist (T1), toolsmith (T8). Referees: skeptic (x2), theorist. All minor revision; I accepted all three after applying the referees' wording and test points myself. Details in `rounds/003/`.

- **experimentalist:** the 7 "mirror without reflection" cores of E-056 pair by a reflection at n = 14, and at every even n 12..18; they fail at n = 13, 15, 17. So that defect is an odd-n effect **for these 7** (chosen because they show it at 13); nothing known yet for the other 132.
- **theorist:** H-021 restated without the mirror clause, covering only cores that pair. For cores whose outside block is interior the shortfall is exactly the H-020 tail minus head (13/13 at n = 13, 3/3 at 14); with the block at an end 17/21. Near-tautological, F-053 has it for `45`; the interior/end split is not yet a committed column.
- **toolsmith:** `batch.py orbits N` (resumable E-052 orbit report), tests pin `45` and `344` at 13.

Promoted: E-059, E-060; H-021 status line (still OPEN); two glossary terms. `isTilting` still not promoted (round 002 q1). Merging `main` was a no-op. No unanswered questions to decide.

**Question for you (chair takes the recommended option if unanswered):** add the full n = 12 and n = 14 censuses (139 cores, about 90 min each, 4 procs) to `OVERNIGHT.md`, to see whether parity governs the other cores? Recommend yes.

---

## Round 002 -- 2026-09-30 -- ordinary

Revisions by experimentalist (T1), skeptic (T3), scholar (T5). Referees: skeptic (x2), theorist. All accepted (skeptic's with three wording points, which I applied). Details in `rounds/002/`.

- **experimentalist:** H-021's mirror clause fails on every reading at n = 13. Strict reading: false for 108 of 109 pairing cores (by construction); 7 cores (`344 366 4044 4403 4404 4405 4605`) hold a strict mirror with no reflection. `3346` and `4056` hold no strict mirror, so H-021's text on `3346` is right under that reading.
- **skeptic:** first counterexample to "orbit = key class": `4056` at n = 16, offsets `{1,2}` are two mirror-image orbits though the key pairs them (derived equivalent via F-026). Also corrected the n = 24 scan: 309 / 155 / 121, not "7 / 123".
- **scholar:** on non-monomial parents (n = 5, 6, 7) Ladkani 2.3(c) rejects exactly the gate-refused vertices; still no gate-admitted rejection except E-032 step 7. Negative control reproduced.

Promoted: E-056, E-057, E-058; E-054 corrected in place; H-021 status line updated (still OPEN); two glossary terms. Merging `main` was a no-op.

**Questions for you (chair takes the recommended option if unanswered):** (1) Promote `isTilting` to the library as a cross-check? Recommend not yet. (2) Restate H-021 without the mirror clause next round (theorist), with the 7 survivors run at n = 14 (experimentalist)? Recommend yes.

---

## Round 001 -- 2026-09-29 -- ordinary

Worked: experimentalist (T1), skeptic (T3), scholar (T5). Referees: skeptic (x2), theorist. All three came back **minor revision**; all three are due for revision next round. Details in `rounds/001/`.

- **experimentalist:** all 139 single-cluster cores of `--max-word 4` close at n = 13. H-021's "exactly when" fails under the literal reading: 129 of 139 cores hold a mirror, so it does not discriminate, and 20 hold one with no reflection (including `3346`, which holds its own mirror). A stricter reading of "mirror" was not run.
- **skeptic:** `3346` and `4056` are not artefacts. The Coxeter key differs at every offset of `3346` for n = 12..40, which proves distinct orbits. `4056` is F-053's reflection with an onset (orbit-verified to n = 15). Two H-020 failures fit the same onset pattern. Referee: pairing beyond 15 is key-only and "cause" is overstated.
- **scholar:** all 61,718 gate-admitted steps at n = 6, 7 pass Ladkani's exact tilting criterion (arXiv:1001.4765, Prop. 2.3(c), unused in `research/` until now), and the E-032 ALARM step fails it. This is not evidence for H-015: the guard is inert at those sizes.

Promoted: E-053, E-054, E-055 (run records with the referees' caveats); H-021 status line annotated, still OPEN; three glossary terms. No F entry yet, pending the revisions. Merging `main` was a no-op.

**Questions for you:** (1) Restate H-021 now (overhang or onset), or settle the strict mirror reading first? (2) Write into `OVERNIGHT.md`: the n = 14 census (about 90 min on 4 processes) and the Ladkani audit at n = 9/10? (3) Promote `tiltingPlus` to the library as an exact gate? I would wait for a second non-monomial negative case.

---

