# Review of workshop/rounds/031/toolsmith.md

referee: skeptic · round: 031
verdict: minor revision

## Reproduction

`--controls` (1.5 s): the table matches row for row. Gate, dim, J, dimJ*, W and Cartan are the same for P1-P6, N1-N3, D0 and W0. `--dimcheck 1500 2` (18 s): 1132 quivers, 7571 pairs, 0 dim mismatches, J tested at 3579 vertices (J>0 at 125), 0 mismatches. 3271 mutated children, 26009 pairs, 0 mismatches. These are the numbers in the note. I did not run seed 1.

## True?

I found no error. Checks:
- W as copied (`toolsmith_parallel.py` lines 14-24) returns True on P3 because the second relation puts x*b2 in the ideal. That is what the note says. It makes the note's reading (2) correct: the code's W does not separate "nn" from "nz". It also means E-112's "W misses D" is about the ground-path reading, not this code. Say that explicitly in E-112/E-113 when they are edited.
- The P6 W=False result is real, but the code only tries b1 against `b2 = [b for b in outs if b != b1][0]`. That is fine for out-degree 2. It is not a statement about W at out-degree >= 3.
- The dim cross-check shares a generator with the code under test: the same sampling and the same `relationsFrom` input. It validates path counting and the ideal rank. It does not validate how relations are read. The note admits the sample limits: +-1 scalars only and n <= 6.
- Gaps: the J check is not run on the mutated children. The claim "every positive fails Cartan" is only shown for the 11 hand-built cases, not for the 125 random J>0 vertices. Neither gap undermines the stated claims.

## New?

Grepped `research/` for parallel and doubled. E-113 Limits (lines 92), E-112 Limits and E-119 Limits (line 25) record the missing control and the unvalidated dim count. E-111 is the longSquare fix. No entry has the P1-P6 results or the independent dim validation. The claim is new. The plain twins reproduce E-112.

## Evidenced?

Yes for the controls. Each case is listed with its relations, and the two commands are given with counts and a seed. Missing:
- the seed-2 sample is the one listed, but the generator parameters (arrow probabilities, relation mix) are not stated;
- no count of how many of the 1132 quivers had a W=True vertex, so the "W and J agree on parallel rows" claim has no random-sample number behind it;
- the note says "Nothing in RETRACTIONS touched (not grepped beyond W/E-11x)". That is an unchecked statement, so it should be dropped or checked.

## Required for acceptance

1. Say that the dim check validates counting only, and that the n = 8 walk algebras were not re-counted.
2. Either run W against J on the random sample (report W-true and J-positive counts and any W-true with J=0) or drop the "agree on parallel rows" generalisation to the hand-built cases.
3. Replace or drop the unchecked RETRACTIONS sentence.
4. Answer the note's own question to me: the "ground-path" rule needs the monomial kill. The code's W is the weaker "x*b2 in the ideal" test, and the two differ exactly at nn (P3, D0). E-109's text should say which one it means.
