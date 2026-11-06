# Round 119 -- proceedings
kind: ordinary · cycle 1, phase 5 (backlog), third round. Merge of origin/main: up to date. No open questions from 118 (nothing to decide).

## Submissions
- **Toolsmith (E-312, R-016), backlog 15 remainder -- accepted (referee skeptic: minor revision; both items answered).** Seed 1614's failure was the round-118 harness encoding -1 as 32002 for a rational-arithmetic library; with -1 the repo agrees with Hom-verified tilting (8 of 8), N5/N6 failures fall 71 -> 8, all 8 computed minimal (E-080 kind). Library untouched; 5 regression tests. The "reproduces E-310" claim was dropped in the response. Refuse-or-minimalise retired. The skeptic's own harness error is a retraction (R-016) of E-310's non-minimal part and H-030's wording.
- **Experimentalist (E-313), backlog 10 -- note (referee theorist: minor revision, answered; no power).** No lowering in 3 non-vacuous cells (k = 4: 2, k = 5: 1); k >= 5 otherwise capped, ValueError diagnosed as an engine limit with parallel arrows. Referee: E-288 already has this null and declined the overnight run. Item 10 converted (theorist question, no overnight).
- **Maverick (E-314), backlog 16, 19 -- accepted narrowed (referee scholar: minor revision, answered).** NO power of E-269 at lam < 1.5 shown at one n = 11 group only (21 of 21 planted NO, 292 of 292 controls YES); none buildable at n = 12. Singular n = 12 group of 47 splits 14 + 33 by SNF(S), independently recomputed by the referee; E-291 amended. Item 16 settled at that one group; item 19 converted (H-017 amendment, H-017b text, overnight, not run).

## Consequences
- R-016 retracts the "presentation dependence" part of E-310/H-030; H-030 narrowed to the E-080 over-admission on minimal length-3 presentations (8 of 936). E-310 and E-291 carry amendment notes. T11 stays open at that narrower form.
- The 47-group split means that group is not one derived class; E-197/E-209 unaffected. No F or R-entry for T10: the R-candidate is still blocked.
- Lesson recorded in R-016: check a harness's encoding against the library's number type before claiming a library failure.

## Backlog marks after 119
Settled/converted: 1, 3-7, 9, 11-14 (as before), 15 (non-minimal part retired by R-016; minimal length-3 anomalies remain in H-030), 16 (settled at one group, open at n = 12), 19 (converted), 10 (converted, no power). Open: 2 (author-only), 8 (H-025 search), H-029 widths 8-9, E-080 shape of the 8 anomalies.
Item 17 (overnight proposals E-204 depth 8+8, exact classifyLength(10), n = 12 growth sweep, k = 5 margin 2): **decided by the chair: none approved; retired for cycle 1.** No result they would produce is blocking an F/R entry; the T10 control that blocks is a search for a crossing step, not these runs. Nothing is added to OVERNIGHT.md.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
(none from step 0; item 17 above is a chair decision)
