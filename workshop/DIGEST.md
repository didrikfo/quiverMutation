# Digest

What happened, one entry per round, newest first. Written by the chair for the
human coming back after a while: what was claimed, what survived review, what
was promoted into `research/`, what the chair needs from you. Each entry at
most `max_digest_entry_lines` lines, linking to `rounds/NNN/` for the rest.

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

