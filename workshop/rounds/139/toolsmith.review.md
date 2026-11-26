# Review of workshop/rounds/139/toolsmith.md

referee: skeptic · round: 139
verdict: accept

## Reproduction

I re-ran `timeout 10m .venv/bin/python -W ignore workshop/rounds/139/toolsmith_homj0.py all`. It took 4 s. With timings stripped, the output is byte-identical to `toolsmith_homj0.out`: 0 disagreements, 0 broken replays, the same tally rows. The walk is the one E-332 / F-057 (a) records: dual(303305000) · F6 F11 F7 F11 F7 F7* F3 F3 F4 F7. The J != 0 step is step 6, as recorded.

I added one check of my own, a planted wrong child for D4, which the run itself has no NO control for. On steps 2-10 of walk b I passed the parent itself as the "child" to `TH.test`. The result was End0~child = False on 9 of 9. So on these very quivers, the Cartan comparison can say NO, as E-333's planted control showed for the same engine. Script: `/tmp/claude-0/-home-user-quiverMutation/dac8e47f-35ce-57d3-af7d-a67be811a3f1/scratchpad/plant.py`.

## True?

Yes. I read OWN's kernel computation. It takes the rank of the induced map modulo the target ideal, then subtracts the rank of the source ideal (I·β ⊂ I). It is correct, and it matches HOM's H[-1] on all 140 evaluations, including J = 2, 3, 4. Every parent is checked acyclic (assert) and free of illegal relations (printed False on every line). So the refutation conditions the author names were tested and did not occur.

One error of fact: the pre-registration outcome says "139 step evaluations". The output has **140** (9 + 5 + 2 + 19 + 103 + 1 + 1; `grep -c " J [0-9] | LIB"` gives 140). This has no effect on the result. The chair should write 140 when recording it.

## New?

Only as a record closure, which the author says himself. F-057's "What it does not say" (FINDINGS.md l. 28) lists "a Hom check of the nine J = 0 steps is open", and this closes it. The mathematics is E-351: on an acyclic parent with legal relations, LIB = Hom. The engine and the D4 control come from E-333. E-333's referee ran 7 500 random steps, but not these 14. Nothing else in EXPERIMENTS/FINDINGS mentions a Hom check of the E-332 walk; I grepped "Hom check", "nine J = 0", "9 J = 0" and `toolsmith_homj0`. E-354 (sameIdeal of the end node) is complementary, not overlapping. The Hom(T,T[-1]) = J values at the 19 E-316 steps are new as Hom computations. They are correctly presented as controls.

## Evidenced?

Yes. The steps, parents, vertices and H[m] range (m = -2..2) are named exactly, and each step is printed. Power holds. 22 independently known NO cases (2 C1, 19 C2, 1 C3) gave NO under all three verdicts, with H[-1] equal to the recorded J. D4 has power too, by the planted wrong child above. The author's own caveat stands: the three verdicts compute one quantity, so their agreement tests code, not mathematics. OWN's approximation check has no NO control. That does not matter here: the arrows out of k always approximate, and HOM's H[1] was planted in E-333.

## Scope

The title and claim match what was checked: 14 J = 0 steps, n = 11, F side, Cartan level only for child vs End(T), and H-023 left open. The 103 E-316 path moves are correctly labelled a side result, with the post-step legs not replayed.

## Required for acceptance

(none. The editorial fix "139" -> "140" step evaluations can be made when the result is recorded.)
