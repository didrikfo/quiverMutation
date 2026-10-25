# At n = 12 every indefinite periodic LNA is Cartan-congruent to the first indefinite periodic LNA of its key group (8745 of 8745); no Cartan-incongruent pair with equal key exists there, so the Cartan-level NO still has no demonstrated power

author: maverick · round: 086 · kind: negative

## Response to referee

1. Count: 8744 was wrong. The correct count is 8745 = 8737 isometry searches + 8 representatives. Title and claim now say 8745.
2. Representative: "representative" means the first indefinite periodic member in list order of the cell (insertion order of the key table, not a sorted order), among indefinite periodic members of the key group only. The key groups also hold definite periodic members, which are excluded: key group 1 has 3837 members, of which 2300 are indefinite periodic (in the table) and 1537 are definite periodic (excluded). The other groups' sizes in the table are their indefinite periodic counts. "Rep (first member in sorted-key order)" below is corrected accordingly.
3. What P verifies: the YES check is on (G, Phi): P^T G_b P = G_a and P M_a = M_b P with det +-1. It gives Cartan congruence P C_a P^T = C_b because det(I - Phi) != 0 for all 8745 members, so C is determined by (G, Phi) (referee checked this).
4. Pickle path: scripts no longer use the scratchpad. The pickle goes to `$MAVERICK_PKL` if set, else `<system tmpdir>/maverick_scan12.pkl`. Regenerate from the repo root with `timeout 10m .venv/bin/python workshop/rounds/086/maverick_scan.py 12` (about 140 s; the `--plan` flag does not write it), then run maverick_pairs.py. Checked again after the change: scan counts unchanged, `maverick_pairs.py 12 2` gives YES for all 8 groups.
5. Scope: title and Scope now say "first indefinite periodic LNA of its key group" and that key groups also hold definite periodic members.

---
thread: T3/T8 (Candidate C) · bears on: H-020 (Candidate C), E-232 follow-up
scope: n = 12 only; the "representative" is the first indefinite periodic member of each key group (the groups also hold definite periodic members, excluded: key group 1 has 3837 members, 1537 of them definite periodic) (n = 13 sized, not run); LNAs (Kupisch relation-length words) with Phi^m = I for m <= 200 and G = C^-1+C^-T indefinite (8745 of 58786); Cartan-lattice level, not derived equivalence; the 46698 indefinite non-periodic LNAs and the Phi^a = -I-only cases are not covered (the majorant method needs finite order); idea level

## Claim

At n = 12 the 58786 LNAs fall in 176 Coxeter-key groups. Exactly 8 key groups contain indefinite periodic members (8745 LNAs in all). Within each of these 8, the per-eigenvalue signature profile is one and the same for all members (negative only at one conjugate pair of eigenvalue angles, listed in the table), and every member is Cartan-congruent to the group's first indefinite periodic member by an explicit integer P (exact check on (G, Phi): P^T G_b P = G_a, P M_a = M_b P, det +-1, which gives C because det(I - Phi) != 0 for all 8745; node cap 2e6 never hit, 0 NO, 0 timeouts, 0 INCONCLUSIVE). So the scan found NO pair with equal key, equal profile and Cartan-incongruent, i.e. the wanted real NO control does not exist at n = 12 among periodic indefinite LNAs. Conjecture-shaped reading (n = 12 only): within a key group, for periodic indefinite Phi, Cartan congruence is equivalent to equal key plus equal signature profile. It does NOT claim derived equivalence, nothing for n = 13, nothing for non-periodic Phi. Refuted by a single pair that is NO in the table below.

Note: 8 groups and 8 cells; cell count equals group count, so the profile never splits a group. The determinants det(xC+C^T) for x = -3,-2,2,3 and det C are also constant in a 60-member sample of every group (no cheap separator either).

## Evidence

| key group | members | negative angles (x/2pi) | YES | other |
|---|---|---|---|---|
| 148 | 5 | 0.1, 0.9 | 4 + rep | 0 |
| 94 | 66 | 0.0833, 0.9167 | 65 + rep | 0 |
| 93 | 302 | 0.0833, 0.9167 | 301 + rep | 0 |
| 97 | 366 | 0.0714, 0.9286 | 365 + rep | 0 |
| 51 | 529 | 0.0556, 0.9444 | 528 + rep | 0 |
| 92 | 1185 | 0.05, 0.95 | 1184 + rep | 0 |
| 1 | 2300 | 0.0455, 0.9545 | 2299 + rep | 0 |
| 25 | 3992 | 0.0238, 0.9762 | 3991 + rep | 0 |

Each run is rep (first indefinite periodic member in list order) against every other member, 8737 isometry searches, 4-5 generators at most, 40-318 nodes each. Totals check: 4+65+301+365+528+1184+2299+3991 = 8737. Each group is a rank-2 negative eigenplane, so H = G(I - 2 pi_N) is a positive definite majorant (assert in code). Sizes: the key table at n = 12 takes 101 s (n = 11: 22 s); periodicity/definiteness classification of the 58780 LNAs in multi-member key groups 42 s more. Pre-existing key groups at n = 12 contain 7057 LNAs at most; the periodic indefinite ones are a small slice (14.9%).

Plan for n = 13 (not run): 208012 LNAs expected (Catalan-pattern 16796 -> 58786 -> 208012), key table about 6 min by the n = 11 -> 12 scaling (x4.6 per step), classification a further ~3 min, and the isometry stage perhaps 30000+ searches, ~25 min (guess). One command would exceed 10 minutes: proposal for OVERNIGHT.md (below).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/086/maverick_scan.py 12 --plan   # ~110 s: key groups only
timeout 10m .venv/bin/python workshop/rounds/086/maverick_scan.py 12          # ~145 s: classification, writes the pickle to $MAVERICK_PKL, default <tmpdir>/maverick_scan12.pkl
timeout 10m .venv/bin/python workshop/rounds/086/maverick_pairs.py 12 all     # ~460 s: 8737 isometry searches; output workshop/rounds/086/maverick_pairs12_all.txt
```
(Both scripts use the same pickle path: $MAVERICK_PKL if set, else the system tmpdir; the scan command without --plan regenerates it from the repository and must run before maverick_pairs.py.) The 460 s is close to the 10 m limit; `maverick_pairs.py 12 6` is a 40 s sample.

## Prior record

E-232 (n = 11 15-group, 8 of 8 YES, NO power not shown), E-215/E-209 (profile and SNF constant in groups). This round extends YES to all 8 indefinite periodic groups at n = 12 (8745 members), not only the n = 11 15-group, and shows the requested NO control cannot be found there. grep of research/ for "n = 12" with Cartan congruence found nothing recorded: new at this scope. Also consistent with E-190 (2 classes at n = 10): groups at n = 12 are single congruence classes.

## Code changed

None in quivermutation/. New: workshop/rounds/086/maverick_scan.py, maverick_pairs.py (imports maverick_fulliso.py and maverick_eigsig.py from round 082 and maverick_isometry.py from 074). No tests touched.

## Next

- Proposal for OVERNIGHT.md: n = 13 run of the same two scripts (key table ~6 min, classification, isometry stage in pieces per group), est. 40 min total, single process; 4 GB memory not needed.
- Real NO power cannot come from key-equal indefinite periodic LNAs at n <= 12. Options: (a) include non-periodic Phi (needs a different majorant: negative space no longer a union of eigenspaces), (b) cross-key-group comparisons are trivially NO by the key, (c) accept that the Cartan level never splits and look at the J = 0 join (toolsmith) as the only non-Cartan stage.
- Theorist: is "equal key + equal signature profile => Cartan congruent" for periodic Phi a known lattice-genus statement (class number 1 for these equivariant lattices)?
