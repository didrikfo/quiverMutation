# Literature check for T7 and S-1: no source found (arXiv title/abstract search) states the Catalan factorisation, the K-threshold law or the orbit count; the gentle fragment of the single-mutation gate is published and reproduced, and Fosse Cor. 4.6 puts A_n in B_n's derived class at every n >= 5 (table confirms n = 6, 7, 8)

author: experimentalist · round: 103 · kind: literature
thread: T7, S-1 · bears on: H-010, H-026 (and E-218 / E-247 wording)
scope: arXiv only (export.arxiv.org API and e-print LaTeX; PDFs not read, so theorem numbers below are counted from `\newtheorem` counters in the LaTeX, re-counted in round 103 with commented-out environments excluded, and marked as such). "Not addressed" means: not found by arXiv title/abstract search plus five opened sources; it is not a claim of absence from the literature (no full-text, zbMATH or OEIS search). OEIS is NOT reachable (HTTPS tunnel 403), so no OEIS lookup was made. 19 queries, about 75 titles seen, 5 sources opened. Computation (plus the n = 17 orbit run) only to check published statements against the repo (n = 3..9 gentle gate; n = 6..8 orbit vs table).
power: n/a (literature). The gate comparison has its own control: the crossed pairing (paper criterion vs the wrong-side repo gate) mismatches on 25-50% of cells.
backlog: none taken

## Response to referee

Verdict received: minor revision. Items 1-4 done; item 5 is the chair's.

1. Theorem numbers: done. I recounted the LaTeX counters with commented-out environments excluded. 2609.39536: Thm 4.2(3) confirmed (it is the criterion), "Thm 7.1" corrected to 7.2 (7.1 is Definition 7.1). 2602.04571: Prop 2.1 corrected to 2.3 (Ex 2.1, 2.2 precede it); it counts ideals of CA_n without idempotents (Dyck paths with 2n steps, C_n), not LNAs (C_{n-1}); text now says so. 2305.06642 Cor. 4.6 recounted: Lemma 4.1, Thm 4.2, Cor 4.3, Lemma 4.4, Thm 4.5, Cor 4.6, matches.
2. Title and scope: done. "Every n >= 5" is attributed to Fosse Cor. 4.6, the table only to n = 6, 7, 8; "no source states" became "no source found in an arXiv title/abstract search", and the scope line says this is not a claim of absence from the literature.
3. n >= 17 basis: done in the text (n = 6..16 computed, E-247; n >= 17 by closure: A_n has no move, moves invertible/mirror-closed). Extra: I ran E-247's script to n = 17 (`timeout 10m .venv/bin/python workshop/rounds/089/experimentalist_orb13.py 17`, 74 s for n = 17, output `workshop/rounds/103/experimentalist_orb17_out.txt`): orbit size 69 632 = 17 * 2^12 as E-247 predicted, closed, mirror-closed, A_n not in it. So computed range is now n = 6..17; n >= 18 rests on the closure argument.
4. Runtime: done, gate script about 3 s.
5. Chair item (E-218 / E-247 / H-026 by new entry): not mine.
Referee's point 5 on the gate comparison (orientation paper = repo LEFT is a convention fixed after seeing the match; check is near-tautological on the repo side and says nothing for rows with a 3-relation) is accepted: it is what the Claim table already says (gentle rows only); I do not claim the orientation as a finding.
Referee's remark on power of the 19 queries accepted: "not addressed" stays arXiv-abstract-level and must not be promoted as "absent from the literature".

## Claim

Verdict per claim, against what I could find and read.

| claim | verdict | source and match |
|---|---|---|
| Single-mutation rule, gate part (E-187 gate: right mutation of an LNA passes iff no relation ends at v+1) | SUPPORTED for gentle LNAs only (all relations of length 2, row entries 0 or 2); NOT ADDRESSED for r_v >= 3 | Deng-Geng-Liu arXiv:2609.39536v2, Theorem 4.2(3) (existence iff for every outgoing beta at k there is an incoming alpha with alpha beta not in I) and Theorem 7.2 (result is again gentle; 7.1 is the definition of mu_k^+). Numbers are my count, PDF unchecked. Repo check: all 2^(n-2) gentle rows x n vertices, n = 3..9, 3 168 cells: paper criterion = repo LEFT gate (paper's orientation), mirror = repo RIGHT gate, 0 disagreements |
| Single-mutation rule, image part (LNA iff r_v = 2 or v = n-1, result (r_{v-1}, r_v) := (2, 0)) | NOT ADDRESSED | no source covers LNAs with relations of length >= 3; the gentle sources give the result quiver and relations of the mutation but nothing I read specialises it to a line |
| Catalan factorisation C_{v-1} C_{n-1-v} per mutation vertex | NOT ADDRESSED as a mutation count. Adjacent published fact: ideals of kA_n without idempotents are counted by C_n, LNAs on n vertices by C_{n-1} | Calvo Cortes-Frost arXiv:2602.04571, Proposition on ideals <-> Dyck paths (label `prop:dyck`, my count 2.3; ideals of CA_n without idempotents <-> Dyck paths D_n with 2n steps, C_n of them; LNAs on n vertices = admissible ideals are C_{n-1}, so the proposition counts ideals, not LNAs); Bekkert-Giraldo-Velez-Marulanda 1910.01494 enumeration of LNAs 3 <= n <= 13 (already in the library summary). The repo's sum over v of C_{v-1} C_{n-1-v} = C_{n-1} is the Catalan convolution, so the total is consistent with these, but they say nothing about which vertices admit a mutation |
| K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15 (E-146, E-263), D >= 3 rule (E-206) | NOT ADDRESSED | no query found a vertex-deletion / stabilisation threshold for Nakayama algebras. The only deletion result in the library is Cor. `removevertex` of 2310.08346 (piecewise heredity only, already recorded) |
| |orbit B_n| = n 2^(n-5) (E-247) | NOT ADDRESSED; OEIS not checked | no count of tilting-mutation or move orbits of Nakayama algebras found. Fosse's table (below) gives only the class members with almost separate relations, a different object |

Two findings that bear on how these records are worded.

1. Fosse arXiv:2305.06642v1, Corollary 4.6 (my count; label `cor:EquivNakayamaAlgebras`) with r = 0, n_0 = n-3, l = 3 gives kA_n/rad^(n-1) = A_{n,(1)}^{(n-1)} derived equivalent to B_n = A_{n,(n-3)}^{(3)} for every n >= 5 (this is the source of "every n"; the table only confirms n = 6, 7, 8) (third displayed algebra of the left column: starts (l-2), length n_0+2; this is the column not affected by the E-275 misprint, which is in the right-hand column). The printed table (Section 5, rows D_6, D_7, D_8) lists both. E-218's "equivalent only at n = 9 and 10" and E-247's "known only at n = 9, 10" are therefore too weak: A_n and B_n are derived equivalent at every n >= 5 by the published corollary. The orbit of B_n is thus a PROPER subset of its derived class: n = 6..17 computed (E-247 to 16, n = 17 here; n = 9..11 re-run by the referee), n >= 17 by closure argument (a path from B_n to A_n would end in a move into A_n, whose inverse is a move out of A_n, and A_n has none; this needs the stated invertibility of the moves, which I did not re-prove here). So n 2^(n-5) is a move-orbit size and cannot be a derived-class size, and no published count is expected to equal it.
2. Reproduction of a published object (the request's "one published count or example"): the n <= 8 table. Orbit of B_n under orbit + free + mirror moves, keeping rows whose relations of length >= 3 are almost separate (consecutive ones share <= 1 vertex): at n = 6, 7, 8 the shapes are exactly {(1,3), (n-3,3)}, all in table class D_n, none outside any class, A_n absent. This is 2 of the 3 members of the table's D_n row, the third being A_n. Also, the table row counts (rows per n = 1..8) are 1, 1, 1, 2, 2, 4, 6, 11, equal to the library's `allQuipusOfOrder` counts (already matched in E-028 against the n = 8 table, `tests/test_checkpointing.py`; I counted the rows in the LaTeX again: 11 for n = 8).

What this does NOT show: that the gate is correct beyond gentle rows; that any of the three unaddressed claims is false. A sharper search (Math Reviews, zbMATH, OEIS, MathSciNet, PDFs of journals) could still find them. "Not addressed" is a statement about arXiv full-text search by title and abstract plus five opened sources.

## Evidence

Queries (all `export.arxiv.org/api/query`, title / abstract / all fields; list in `experimentalist_arxivq.py`, `experimentalist_arxivabs.py`): Nakayama x derived equivalence (x mutation, x classes, x Catalan, x vertex, x orbit, x Kupisch, x classification), "tilting mutation" x Nakayama, "silting mutation" x Nakayama, "tilting graph" x Nakayama, "Coxeter polynomial" x Nakayama, gentle x derived x counting. Rate limit 429 forced 6 s spacing and retries; one query returned only after retry.

Opened (LaTeX e-print, versions as of 2026-10-31): 2609.39536v2 (gentle tilting mutation; no mention of Nakayama), 2602.04571 (Dyck paths, LNAs), 2406.00254 (Sen, syzygy filtration, defect; nothing on mutation or counts), 2305.06642v1 (already in library; re-read Sect. 5, Cor. 4.6, Obs. 5.1), 2607.28767 / 2609.38456 abstracts only (mutation of gentle ideals; polytope interpolation: neither addresses derived classes of Nakayama algebras). Seen, not opened: 2610.25312 (derived classification of multifold extensions, affine type), 2002.02990 and 1309.2216 (counts of tau-tilting modules over Nakayama algebras: a different count, of modules, not of algebras reached by mutation), math/0607348 and 1904.02555 (gentle derived invariants).

Gate comparison output, columns rows x vertices, matches of (right gate, mirror B) and (left gate, paper A), then the crossed controls:
```
n  cells  R~B  L~A  R~A  L~B     n  cells  R~B  L~A  R~A  L~B
3     6     6    6    2    2     7   224   224  224  160  160
4    16    16   16    8    8     8   512   512  512  384  384
5    40    40   40   24   24     9  1152  1152 1152  896  896
6    96    96   96   64   64
```

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/103/experimentalist_gentle_gate.py 9    # about 3 s, output experimentalist_gentle_gate_out.txt
timeout 10m .venv/bin/python workshop/rounds/103/experimentalist_fosse_check.py      # under 1 s
timeout 10m .venv/bin/python workshop/rounds/103/experimentalist_arxivq.py           # queries, about 3 min with spacing
```
`tests/paper_classification.py` supplies the table; the e-prints are not committed (small, re-fetchable by id).

## Prior record

E-218 and E-247 (scope lines) say A_n / B_n equivalence is known at n = 9, 10 only; this is the item to correct (finding 1). E-275 already read Cor. 4.6 for the mixed-digit pairs, and found its right-column misprint; it did not apply it to A_n / B_n. E-187 / E-205 / E-208: the gate is stated as the repo's reading of steps 1-7; the new fact is an independent published iff on the gentle fragment. No `research/RETRACTIONS.md` item involved. Summary for promotion: a `research/literature/2609.39536-gentle-tilting-mutation.md` entry would hold Thm 4.2 (iff gate for gentle algebras, loops included), Thm 7.2 (result gentle), and the caveat that it covers no non-gentle LNA; 2602.04571 earns two lines (ideals without idempotents = Dyck paths, Prop 2.3; LNAs are C_{n-1}).

## Code changed

None in the library. New files in `workshop/rounds/103/`: the four scripts above and `experimentalist_gentle_gate_out.txt`. No tests run (nothing under `quivermutation/` touched).

## Next

- scholar: PDF check of the theorem numbers marked "my count" (2609.39536 Thm 4.2 / 7.2; 2305.06642 Cor. 4.6; 2602.04571 Prop. 2.3; the LaTeX count was redone in round 103 and agrees with the referee's); OEIS lookup of 12, 28, 64, 144, 320, 704, 1536 and of 3, 7, 16, 36, 80, 176, 384 (n 2^(n-5)) from a machine with access (UNVERIFIED from memory: I do not recall a match).
- chair: correct the A_n / B_n wording in E-218, E-247 and H-026 text ("derived equivalent for every n >= 5 by Fosse Cor. 4.6; not joined by moves"); the orbit is a proper subset of the derived class.
- theorist: the orbit-versus-class gap is itself a question: is A_n reachable by any tilting mutation outside the E-187 moves (this is T10's guard question at n = 6)?
- experimentalist (later): extend the gate check to rows with one relation of length 3 and rows with 2s next to it, using the published Ladkani 1001.4765 Prop 2.3(c) as the source (already in the library) instead of an arXiv search.
