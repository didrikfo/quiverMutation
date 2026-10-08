# The arXiv is unreachable, so 2.31/2.32 and CHZ 3.6 stay UNVERIFIED; the path-wise CHZ Cor 3.6 needs a "socle spanned by paths" hypothesis, which "monomial" supplies

author: scholar · round: 042 · kind: negative
thread: T5 · bears on: H-015, E-066, E-122, E-066's CHZ flag

## Claim

(1) I could not read either paper. The proxy answers 403 to CONNECT for arxiv.org, export.arxiv.org and ar5iv, and WebFetch fails DNS (`ENOTFOUND arxiv.org`). The standing permission does not help while the proxy denies the host. Every statement below about what AI or CHZ say is therefore UNVERIFIED and rests on the repo's summaries (`research/literature/1009.3370-silting-mutation.md`, `2509.12983-chz-...md`) and on my memory.
(2) What I can settle without the papers: CHZ Cor 3.6 in its path-wise wording is not equivalent to the module-theoretic Prop 3.5 unless supp soc P_i is the set of ends of tail-maximal *paths*. This holds for monomial I and fails for E-066's parent. So the word "monomial" (or an equivalent hypothesis) is needed in the path-wise corollary, or the corollary is stated for I admissible in a weaker sense. Whether the paper says so is UNVERIFIED.
(3) The AI citations in E-066, E-122 and the `tiltingPlus` gate are coherent with the summary's statement of 2.32(b). I found no inconsistency, but "coherent with a summary" is not a check against the printed theorem.

## Evidence

Fetch attempts (Bash, WebFetch): `https://arxiv.org/pdf/1009.3370`, `.../pdf/2509.12983`, `export.arxiv.org/abs/1009.3370`, `arxiv.org/abs/2509.12983`, `ar5iv.labs.arxiv.org/html/1009.3370`. All: "CONNECT tunnel failed, response 403". The proxy status endpoint lists `arxiv.org:443` under `recentRelayFailures: connect_rejected`. Nothing downloaded. No file was added under `research/literature/`.

**Argument for (2).** Prop 3.5(2), as the summary gives it, is stated through Phi+({i}) = supp soc P_i. Write soc P_i = {x in e_i Lambda : x rad Lambda = 0}. This is a subspace of e_i Lambda, and a nonzero element need not be a path class. Cor 3.6 is stated through "nonzero path p with no nonzero extension pq". Call a path *tail-maximal* if p·a = 0 in Lambda for every arrow a. Then:
- supp soc P_i contains {ends of tail-maximal paths}, always (such a path class lies in the socle);
- equality holds when I is monomial, because then the nonzero paths form a basis of Lambda and the socle is spanned by tail-maximal paths;
- equality fails in general: E-066's parent has c = [8,6,4] + [8,10,4] in e_8 Lambda e_4 with c·(4>9) = 0, while neither summand is tail-maximal, since each extends to a nonzero path. So 4 is in supp soc P_8, but no tail-maximal path from 8 ends at 4.
Hence the path-wise Cor 3.6 can pass where Prop 3.5's socle form fails. That is exactly the E-066 flag. Only the direction "path-wise passes, socle fails" is shown. The path-wise test never rejects something the socle test accepts, so it is a weaker gate on non-monomial I. The paper may restrict to monomial I, or may define "path" or "maximal" in a way that avoids this; UNVERIFIED.

**The repo's gate.** My notebook (r033, r038) has `mutationIsPossibleAtVertex` testing simple paths, with the cyclic case unchecked. That is the same weakness (path-wise on a possibly non-monomial parent), independent of the CHZ wording. The tests (`tests/test_gate_without_tilting.py`) pin a gate-admitted, non-tilting step, which is the expected behaviour of a necessary-only gate.

**Citation check (all against the summary's text, UNVERIFIED against the paper):**
| citation | used for | status |
|---|---|---|
| AI Thm 2.32(b), M tilting, D contravariantly finite, mu^-(M;D) tilting iff each M has a right D-approximation g with Hom(g,D) injective | `tiltingPlus` (E-122), E-066 | statement matches the project's use at M = add A, D = sum of P_j (j != v); the condition sits at M = P_v. E-122's reduction to J_i = Hom(S_v, e_iA) is a separate derivation and does not depend on the paper beyond this statement. |
| AI Thm 2.31 "mutation of a silting subcategory is silting" | E-121, notebook r035 | The summary says "any mutation"; from memory the theorem requires D to be functorially finite. In K^b(proj A) with A finite-dimensional this is automatic, so the use is safe either way. The number "2.31" vs Prop/Def numbering is UNVERIFIED (assignment text and notebook disagreed). |
| "right mutation = repo left mutation" | convention note | convention only; mu^- in the summary, mu^+ is the dual statement. Not checkable without the paper. |
| CHZ "Cor 3.6 is iff, so for LNAs the gate is exact" (summary line ~118) | T5 | conditional on the repo's right tilting mutation being the HRS tilt at (filt S_v, S_v^perp), which the summary itself says it does not prove. For LNAs (monomial) the path-wise test is the socle test, so the conditional is the only gap. |

## Reproduction

```
curl -sS -m 60 -o /dev/null -w "%{http_code}\n" https://arxiv.org/pdf/1009.3370   # 000, CONNECT 403, seconds
curl -sS "$HTTPS_PROXY/__agentproxy/status"                                        # lists arxiv.org:443 connect_rejected
```
No script; the argument for (2) is by hand.

## Prior record

E-066 raised the "monomial" flag and the 2.32(b) identification; E-122 derived J_i from 2.32(b). My (2) turns the flag from "may need" into a short proof of when the two forms agree, but it is not new information about the paper. Nothing in `research/RETRACTIONS.md` concerns this (not grepped beyond E-066/E-122; see Next). The fetch failure is the same as round 006 and round 038.

## Code changed

None.

## Next

- Chair/human: arXiv access is the blocker for item 3 for the third time. Either whitelist arxiv.org in the proxy or paste the statements of AI Def 2.30, Thm 2.31, Thm 2.32 and CHZ Prop 3.5, Cor 3.6 (about 30 lines). Until then, T5 literature claims should carry UNVERIFIED.
- theorist: the only live gap is the HRS-tilt identification for the repo's mutation, not the citations.
- toolsmith: make the gate test all path combinations (socle of e_i A), not simple paths; this closes the E-066-type and cyclic gap on non-monomial parents. Not done here.
