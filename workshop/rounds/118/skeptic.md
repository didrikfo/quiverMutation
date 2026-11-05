# Backlog 13-15: the H-020 class-B relation delta is a per-rule constant with no formula found; F/R convention holds for parallel-arrow and cyclic parents but the repo's gate and F/R fail on non-minimal relation presentations; E-305's walk output reproduces

author: skeptic · round: 118 · kind: backlog
thread: T4, T5, T10 · bears on: H-020 (H3/H4), E-258, E-286, E-266, E-305, E-080, E-211
scope: (14) the 14 / 54 / 94 notLNA floating class-B rules of widths 5 / 6 / 7 of E-286, lengths <= 8 only (width 7: length 8 only, one match per rule); (15) random parents n = 4..6: parallel-arrow monomial acyclic (P, 1 115 parents), non-monomial acyclic with commutativity relations (N, n = 5, 6 only for the table: 687 parents), cyclic monomial (C, 923 parents); seeds 0..N of a fixed generator, a capped sample; relations/multiplication of the child not compared, arrow length >= 5 not covered; (13) the two scripts as written.
power: (15) the other side matches in 0 of every "power" step (about 7 600 F and 7 600 R tilting steps with T_R != T_L data); NO is reachable: the same algebra given with one redundant relation makes the repo fail 8 of 8 vertex/side cases and the minimal presentation 0 of 8 (seed 1614, n = 5). (14) a coarse-key control says NO: keying the delta by len(old window) is functional on 0 of 2, 0 of 3, 1 of 4 groups; the feature-search null reaches the same score as the best feature (see Claim), so the formula search has NO power to say yes. (13) none needed, a reproduction.
backlog: 13 referee part settled (output equals the recorded one; 3 open points of E-305 still open, see below) · 14 converted (proposed H-029; H-020 H3/H4 retired as stated for B rules) · 15 settled for parallel and cyclic parents, converted for non-monomial (proposed H-030)

## Claim

**Item 14.** For the 162 notLNA class-B rules of widths 5-7 the (arrows, relations) delta of the first step is a function of the rule alone: one value over all matching rows and window positions for 14 of 14, 54 of 54 and 94 of 94 rules, and a delta read off lengths <= L-1 predicts the delta at length L in 82 of 82 rule-length pairs (widths 5 and 6; width 7 has one match per rule, so its 94 of 94 is vacuous). The delta takes only two values, (1,1) and (1,2). I found no formula in the rule data: of 120 one- and two-feature predictors fitted on widths 5+6 (68 rules) the best scores 0.872 on width 7 (baseline 0.66), but label-shuffled data reach 0.872 as well (max of 20 shuffles), so the pre-registered 0.90 cut and the null both reject it. So E-286's "narrowed only" becomes: H3 (completeness of the rule table, beyond which the delta needs the host row) is false for the first-step delta; H4 (a mechanism computing the delta from the rule) is not found. I retire H3/H4 as stated for class B and convert the rest to H-029 below. Refuted by a width 8+ rule with two deltas, or by a rule-formula that predicts held-out width 8.

**Proposed H-029 (next cycle).** *For a class-B floating rule d (notLNA at step 1), the first mutation changes (arrows, relations) by (1, 1) or (1, 2), a number fixed by d alone (82 of 82 length transfers at widths 5, 6). It is not a function of len(added) - 1 (E-286: 14/14 at width 5, 84/162 at width 6) and no 1-2-feature rule of the window lists transfers to width 7 above chance.* Test: widths 8-9, lengths 9-10, rule-keyed table from widths <= 7 predicting delta (needs the rule table for widths 8, 9, about 10 min per width); falsified by one rule with two deltas.

**Item 15.** The side convention (repo F child = End(T_R)^op, R child = End(T_L)^op, labelled arrows and, for acyclic children, labelled Cartan; the other side never) and gate = tilting hold for P (0 failures in 2 882 F and 2 881 R tilting steps, n = 4..6, gate = tilting in every pair) and for C (0 failures in 692 F and 660 R tilting steps with a repo result, gate = tilting in every pair; the repo refuses 1 212 F and 1 233 R further tilting steps with `CyclicDescendantsError`, which is the documented E-211 refusal, so nothing is tested there). For N (non-monomial) it fails, and the failure is a precondition: with a minimal relation set and commutativity relations of arrow length 2 there are 0 anomalies in 360 parents (2 869 gate/tilting lines, all matches); with a minimal set that has an arrow-length-3 commutativity relation 7 anomalies in 33 parents, all of the E-080 kind (gate True, tilting False; never the reverse; never a mismatch at a tilting step); with a NON-minimal relation set that contains a redundant arrow-length-3 relation, 74 anomaly lines in 206 parents, both directions (gate False with tilting True, gate True with tilting False) and child mismatches at tilting steps. The same algebra with the redundant relation deleted gives 0 anomalies (seed 1614: 8 of 8 vs 0 of 8). Not claimed: that the repo is wrong in a supported setting: the procedure's hypothesis (minimal relations, as the paper) is not stated at the input of `isMutable` / `mutateAtVertex`, only `reduce` restores minimality on the output; I did not check the docstrings beyond `mutateAtVertex` (no mention). Refuted by a minimal length-2 commutativity parent with an anomaly. **Proposed H-030:** *on acyclic parents the repo's gate equals tilting and F/R equal End(T_R)^op / End(T_L)^op exactly when the relation set is minimal and has no commutativity relation of arrow length >= 3 (E-080 shapes); redundant relations alone, at arrow length 2, do no harm.* Test: n = 7, arrow length 4, and a minimalising wrapper as control.

**Item 13 (referee).** After rebuilding the checkpoint with `rounds/069/skeptic_c1_13.py 1 540 23000` twice, `theorist_steps.py` and `theorist_allJ.py` give the recorded output: 13 of 13 steps YES (HH parent = child = (1,0,0)); corruptions 69, HH differs in 30, of those Cartan equal to the child in 13 (the numerators E-305 lacked); allJ: 75 J != 0 records at nexp 23 000, kept 19 HH equal, refused 56 (54 equal, 2 differ, child (1,1,0) at v = 4). Equal to E-305/`theorist.md` in every number. Caveat: another process ran the same 069 rebuild on the shared `/tmp/skeptic_c1_13_ck1.pkl` at the same time (two such processes seen in `ps`); the walk is deterministic so the checkpoint equals, but the first run ended at a wall-clock deadline, so "equal" is shown by the printed counts (`distinct J != 0 steps 67`, 75 in all, `D == 0: 13`), not by a hash.

## Evidence

Item 14 (`skeptic_ruledelta.out`): matches per rule 8 (w5), 3 (w6), 1 (w7); delta distribution w5 {(1,2): 6, (1,1): 8}, w6 {(1,2): 26, (1,1): 28}, w7 {(1,1): 32, (1,2): 62}. Coarse keys as functions: `len(added)` 2/2, 2/3, 2/4 groups (it fits only the groups it was built from), `len(old window)` 0/2, 0/3, 1/4, width alone 0/1. Feature search (`skeptic_rulefeat.out`): top four feature pairs all 0.872 on width 7, shuffled-label max 0.872, mean 0.738.

Item 15, counts (`skeptic_conv3_*.out`; generator parameters in the script; seeds P4 1500, P5 300, P6 600, C4 1500, C5 300, C6 400, N5 4000, N6 3000, N4 4000 with 61 parents and 0 exceptions before tagging):

| kind | n | parents | F tilting steps matched / R | other-side match (power) | gate != tilting | mismatches |
|---|---|---|---|---|---|---|
| P parallel monomial | 4, 5, 6 | 607, 142, 366 | 1 293+377+1 212 / 1 299+381+1 201 | 0 | 0 | 0 |
| C cyclic monomial | 4, 5, 6 | 587, 153, 183 | 417+147+128 / 403+131+126 (arrows only when child cyclic: 144+49+54 / 133+42+46) | 0 | 0 | 0 |
| N, minimal, comm. rel. length 2 | 5, 6 | 143, 217 | all | 0 | 0 | 0 |
| N, minimal, length 3 | 5, 6 | 7, 26 | all but 7 gate-over-admissions | 0 | 7 (gate True, tilting False) | 0 |
| N, non-minimal, length 2 | 5, 6 | 24, 64 | all | 0 | 0 | 0 |
| N, non-minimal, length 3 | 5, 6 | 48, 158 | 74 anomaly lines in 1 849 steps | 0 | both directions | yes |

The seven minimal anomalies: R 871 v1 (n = 5), F 514 v4, R 1629 v2, F 1806 v4, F 1944 v3, R 2257 v1, R 2794 v1 (n = 6). Seed 871: relation (0,1,3,4) = (0,1,2,4) on 0 -> 1 -> {2,3} -> 4: a prefix square, at vertex 1 on the R side; that is the mirror shape of E-080's square-with-suffix (not checked as the mirror). My engine is presentation-independent (quotient by the ideal, row reduction mod 32003), so the redundant-presentation failures are the repo's, with the caveat that I did not locate the failing step inside `procedure` or `reduce`. Engine check: the monomial P and C runs reproduce the E-258 verdict; hand self-test (`selftest`) dimensions; no test against an external system.

Not covered from E-258: multiplication of the child (relations as such); arrow length >= 5; n >= 7; relations with coefficients other than 1, -1; non-acyclic non-monomial.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/118/skeptic_ruledelta.py 5 8   # also 6 8, 7 8 (about 5 s each; writes the json for the next)
.venv/bin/python workshop/rounds/118/skeptic_rulefeat.py                               # 2 s
.venv/bin/python -W ignore workshop/rounds/118/skeptic_conv3.py selftest
for a in "P 4 1500" "P 5 300" "P 6 600" "C 4 1500" "C 5 300" "C 6 400" "N 5 4000" "N 6 3000"; do set -- $a; ALLEXC=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/118/skeptic_conv3.py sweep $1 $2 $3; done   # 1 to 8 min each
timeout 10m .venv/bin/python -W ignore workshop/rounds/069/skeptic_c1_13.py 1 540 23000   # twice (9 min, then 2 min)
timeout 10m .venv/bin/python -W ignore workshop/rounds/115/theorist_steps.py; timeout 10m .venv/bin/python -W ignore workshop/rounds/115/theorist_allJ.py   # about 2 min each
```

## Prior record

E-286 has "rule-keyed table 70/0 on length 8" with seen rows; the new points are the cross-length transfer (82/82), the two-valued delta, and the feature null. E-258 lists non-monomial, parallel, cyclic parents as open: filled here; E-080 is the known non-monomial gate/tilting difference (length 4, suffix square) and the 7 minimal anomalies are the same kind; E-211 is the cyclic refusal. grep of EXPERIMENTS/FINDINGS/RETRACTIONS for "non-minimal / redundant relation": nothing about redundant input relations (only `reduce`'s output minimality, procedure.py l. 604-625).

## Code changed

None in the library. New scripts in `workshop/rounds/118/`: `skeptic_ruledelta.py`, `skeptic_rulefeat.py`, `skeptic_conv3.py` (imports `rounds/091/skeptic_convention.py` and overrides its Dmat and comp_classes for vector-valued composition). No tests (scripts only).

## Next

toolsmith: decide whether `isMutable` / `mutateAtVertex` should refuse or minimalise a non-minimal relation input (the 1614 case is a ready regression test: 5 vertices, arrows 0>1 0>2 0>4 1>3 1>4 2>3 3>4, relations 0>1>3 - 0>2>3 and the length-3 consequence), and add it to the docstring; theorist: E-305 open points still open: HH^2 = 1 for the ab = 0 square (unchecked), tree status of the two refused children (HH^1 = 1 forces a non-tree underlying graph if HH^1 = 0 for tree algebras, as the submission cites; not checked); chair: H-029 and H-030 (non-minimal presentations make the F/R gate presentation dependent) as next-cycle entries.
