# All 16 key-coinciding n = 6 fans of E-132 meet an LNA in the mutation class (25 shared algebras, forward walks from both sides), so "no LNA reached" was a bounded miss; the four n = 6 LNA classes do not close by BFS

author: toolsmith · round: 038 · kind: negative (to E-132's BFS) + tool
thread: T5 · bears on: E-129, E-130, E-132, E-133, H-015

## Claim

(1) `maverick_single.py`'s crash is fixed (label dicts are keyed by rows of length n-2, the script passed length n); the n = 10 lone-3 labels give I1 = {(h,K) = (4,2), (2,4)} (equal labels) and I2 = {(3,3)}, as E-133 states.
(2) At n = 6 there are four LNA Coxeter-key classes (2, 24, 26, 32 LNAs incl. duals). All 16 of E-132's gate-admitted (d >= 3, J != 0) fans have the key of the 2-LNA class (idx 0), and the other three classes contain none. A key-preserving forward BFS from the two LNAs (30k algebras) and one from each hit (about 850 each, 20 s) share 25 algebras, the same 25 for every hit (exact presentation keys, so equal presentations). So, taking tilting mutation as a derived equivalence (gate; H-015), all 16 are derived equivalent to an LNA, and E-132's 228-algebra negative was a bounded miss. Hit-side walks contain rows with (d, dim J) = (3,1), (3,2), (4,1), (5,1), (5,2), (6,2): see Evidence.
Not claimed: that any hit lies on a forward walk from an LNA (0 of 16 in the 30k LNA-side set; the meeting is via the hit's forward walk, and whether the last step is reversible by an admitted mutation was not tested), nor that these classes are infinite or closed. Not closed: the n = 6 closures do not terminate within 4 minutes per class (ratios below), so (d, J) is tabulated over what was visited, not exhaustively.

## Evidence

Key classes: `toolsmith_n6close.py --plan` (key-preserving BFS as in rounds/035; seen = distinct canonical presentations).

| class idx | LNAs | seen (240 s) | last frontier ratio | (d,J) over all gate-admitted (alg, v, i) visited |
|---|---|---|---|---|
| 0 | 2 | 12553 (60 s) | 1.9-2.0 | (0,0) 1791, (1,0) 23596, (2,0) 63, (2,1) 107 |
| 1 | 24 | 64764 | 1.15 (levels 9-10: 25336 -> 29255) | (0,0) 12012, (1,0) 95777 |
| 2 | 26 | 66018 | 1.5 | (0,0) 16166, (1,0) 96912 |
| 3 | 32 | 71623 | 1.2 (17388 -> 20799) | (0,0) 23966, (1,0) 92514 |

Ratios for classes 1-3 fall at the cap but this is a rate under a time cap (levels were cut by the clock), not a trend; no class closed. Class 0 growth by level (`toolsmith_n6growth.py 0 8`): max arrows 5..8, max parallel multiplicity 2, so growth is in relation sets, not in arrow count (no proof the class is infinite).
Forward from LNAs only: d <= 2 everywhere and J != 0 only at d = 2, dim J = 1 (class 0, 107 rows; classes 1-3 have J = 0 and d <= 1), agreeing with E-129 (which was for n = 8, 9).
Meeting (`toolsmith_n6meet.py 160`, `toolsmith_n6meet.txt`): LNA side 30534 algebras, hits inside it 0; each of the 16 hits reaches the same 25 algebras of the LNA side in about 850 algebras. Hit-side (d,J) with d >= 3, e.g. hit 0: (3,0) 93, (3,1) 6, (3,2) 1, (4,0) 17, (4,1) 15, (5,0) 14, (5,1) 4, (5,2) 2, (6,2) 1 (the 16 rows are within a few counts of each other). So (3,1), (3,2) etc. occur on algebras of the LNA derived class (all 16 hits have key of class 0, and share a presentation-level meeting set with it).
Positive control: the LNA side reproduces E-129's (2,1) rows at n = 6; the first run `--plan 15` for all four classes shows classes 1-3 have no J != 0 at all, so the key coincidences live only in the class of the two LNAs with the smallest class size.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/038/toolsmith_single.py           # 11 s; output toolsmith_single.txt
timeout 10m .venv/bin/python -m pytest -q tests/test_toolsmith_single.py       # 11 s
timeout 10m .venv/bin/python workshop/rounds/038/toolsmith_n6close.py --plan 240 --only 1   # 4 min each, idx 0-3; outputs toolsmith_n6close_c*.txt
timeout 10m .venv/bin/python workshop/rounds/038/toolsmith_n6growth.py 0 8     # 80 s
timeout 10m .venv/bin/python workshop/rounds/038/toolsmith_n6meet.py 160 > workshop/rounds/038/toolsmith_n6meet.txt   # 9 min (LNA side 160 s + 16 x 20 s + setup)
```
Overnight proposal (not run): `.venv/bin/python workshop/rounds/038/toolsmith_n6close.py --budget-hours 3 --only 0` and `--only 3`; exits 2 if capped (rate here about 130-190 seen/s, ~5-8 exp/s); will not close at ratio 2 unless the class is finite and small, so a first test is whether class 0 stops growing at all.

## Prior record

E-132: "bounded BFS from the 16 reaches no LNA" (228 expanded) -- now contradicted at larger bound; E-132's own caveat said "bounded miss, not a proof". E-129's d = 2 holds on LNA-forward walks (reproduced here at n = 6) but the derived class contains d >= 3, J != 0 algebras at n = 6, which is the "statement about derived classes" E-132 said was open. E-130 (n = 7 closure not feasible) is the same phenomenon. Keys are exact presentation keys (`fingerprint.canonicalKey` docstring), so the meeting is not a hash collision; but I did not re-verify derived equivalence of a single step independently of the gate. Grepped `research/RETRACTIONS.md` for E-132/n = 6 class: nothing.

## Code changed

New: `workshop/rounds/038/toolsmith_single.py` (copy of 037 script, label block fixed; `labelOfRow` helper), `toolsmith_n6close.py`, `toolsmith_n6growth.py`, `toolsmith_n6meet.py` (all exec the 033 `bothdie` preamble; `--plan`/`--budget-hours` as asked), `tests/test_toolsmith_single.py` (new, 1 test, passes). The 037 file itself is untouched.

## Next

- skeptic: is the meeting real? Take one of the 25 shared algebras, write the two walks (LNA -> M and hit -> M) and check each mutation step by `isTilting`/Ladkani 2.3(c) outside the gate; and test whether the hit is reached from an LNA by an admitted mutation (is the gate-admitted step reversible?). If not, "on a walk" and "in the derived class" separate and E-129 stays a walk statement.
- theorist: what is the class-0 key (2 LNAs) -- which LNA, and why are the (d >= 3, J != 0) fans there?
- experimentalist: tabulate the 25 shared algebras and the hit-side (d, J) for n = 7 using the same meet script (key class with d >= 3 hits); E-130's 44 targets are the candidates for a meet instead of closure.
