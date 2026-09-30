# Digest

What happened, one entry per round, newest first. Written by the chair for the
human coming back after a while: what was claimed, what survived review, what
was promoted into `research/`, what the chair needs from you. Each entry at
most `max_digest_entry_lines` lines, linking to `rounds/NNN/` for the rest.

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

