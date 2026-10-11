# The LNA-key test has no usable base rate: 100% on walk rows by construction, 0% on every category of the layered family (circuit-free and pure-W included)

author: toolsmith · round: 033 · kind: negative
thread: T5 · bears on: E-123, E-115, E-116

## Claim

Agenda item 2 asked for the fraction of gate-admitted out-degree 2 walk rows at n = 8 c0 that carry an LNA Coxeter key. Answer: 100% (3 585 of 3 585 out-degree 2 rows, from 3 000 expanded algebras), but this is forced, because the walk only keeps algebras whose key equals the c0 LNA key. It is not a base rate for E-123. The informative comparison is inside E-123's own layered family (i -> {a_k} -> v -> {t1,t2}, n = 8): 0 of 2 704 members have an LNA/dual-LNA key, including 0 of 1 003 circuit-free members (J = 0) and 0 of 108 pure-W members (J != 0, the shape the walk actually realises). So the 962 / 900 key absences of E-123 are what every member of that family does; the test cannot separate non-W circuits from W or from circuit-free members there. It does NOT claim that non-W circuits could occur on walks (still open), nor that the family is representative of non-LNA-derived algebras.

## Evidence

Walk, n = 8 class 0, 3 000 expanded algebras (BFS, key-preserving mutations, as in E-116/skeptic_loose), every gate-admitted vertex; out-degree, parent in LNA keys, child in LNA keys, child in class c0:

| out-degree | rows | parent has LNA key | child key = c0 | child has LNA key but not c0 |
|---|---|---|---|---|
| 1 | 8 621 | 8 621 | 8 621 | 0 |
| 2 | 3 585 | 3 585 | 3 543 | 42 |
| 3 | 559 | 559 | 559 | 0 |
| 4 | 5 | 5 | 5 | 0 |

(Run on 600 expansions: 598 out-degree 2 rows, 2 not in c0; same pattern.) Parent keys: 100% by construction. Children: every legal child has some LNA key (18 of 18 classes hit is not claimed; only "in the union of the LNA and dual-LNA key sets of length 8", which has 11 keys).

Layered family, all set-partition members (circuit-free ones included, which E-123's script skips):

| m / n | category | J != 0 | members | LNA key |
|---|---|---|---|---|
| 3 / 7 | circuit-free | no | 107 | 0 |
| 3 / 7 | pure W | yes | 18 | 0 |
| 3 / 7 | non-W circuit | yes | 100 | 0 |
| 4 / 8 | circuit-free | no | 1 003 | 0 |
| 4 / 8 | pure W | yes | 108 | 0 |
| 4 / 8 | non-W circuit | yes | 1 593 | 0 |

Non-W circuit count 1 593 exceeds E-123's 900 because this table counts every non-W shape including those with a single path in J (E-123 excluded them); not reconciled further. Note E-123 already says its control was empty for m >= 3; the circuit-free row now shows the empty control is structural: the family has no LNA-keyed member at all (probably because of the two sinks t1, t2 with a shared source, a shape that LNA keys never take; not checked).

Observations a reader should not over-read: (a) the 42 out-degree 2 children with an LNA key outside c0 are presumably the rejected/non-tilting rows (the 61 rejects of E-116 came from a larger 5 145-expansion run); not matched row by row. (b) The key test is necessary only, as E-123 says.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/033/toolsmith_baserate.py layers 3   # 3 s
timeout 10m .venv/bin/python workshop/rounds/033/toolsmith_baserate.py layers 4   # 18 s
timeout 10m .venv/bin/python workshop/rounds/033/toolsmith_baserate.py walk 8 0 3000   # 190 s; output in toolsmith_baserate_walk_n8c0.txt
```

## Prior record

E-123 reports the 0 hits and admits "no base rate"; its pure-W control being empty at m >= 3 is recorded. The circuit-free zero and the 100% walk figure are not in `research/`. Nothing in RETRACTIONS touched. This result weakens E-123's last clause ("none ... has an LNA Coxeter key") from weak evidence to no evidence in that family; the working conclusion (exclusion must come from derived equivalence to an LNA) is unchanged but the key test cannot support it.

## Code changed

New `workshop/rounds/033/toolsmith_baserate.py` only (it contains one dead first `exec` line in the layers branch, harmless). No library change, no tests needed.

## Next

- Experimentalist/theorist: a family where the key test can fire, i.e. one that contains at least circuit-free members with an LNA key (e.g. LNA-derived algebras with a second source added, or the E-123 shape with the sinks attached to different vertices), then compare W vs non-W circuit hit rates. Without a positive control of that kind, key absence should not appear in STATE as evidence.
- Skeptic: match the 42 non-c0 out-degree 2 children to E-116's rejects.
