# Toolsmith, round 108 (conference, phase 2 research)

author: toolsmith · round: 108 · kind: conference note (nothing run, no other file edited)

**Most promising question (T10 / H-027).** Does the E-287 tilting path B_n -> A_n pass through T whose End is NOT an LNA (or whose step the gate refuses), and does that happen at n = 7, 8, 9 with the same depth law (3, 4, 5?)? Why: it turns the existential Fosse Cor. 4.6 into an explicit certificate outside the move set, which is exactly what T10 says it lacks (join searches cannot test the J = 0 premise). It is cheap: the machinery (objects = D_n roots, Hom table, exchange BFS) already runs n = 6 depth 5 in 15 s, and the missing step is only "classify End of each intermediate T" (Cartan + relations check) plus n = 7 (42 objects).

**One question outside the main thread (T6/T7 neighbourhood: a library-side one).** A per-edge tilting certificate in the library: given a J != 0 key-kept edge, `checkJZero` already refuses J = 0 failures; the same Hom(T,T[m]) test (rounds/050 stepTest, 2 s per edge) could return tilting / silting-only for every one of the 70 failing n = 7 edges and the 7 rebuilt n = 9 edges (E-282). That decides, edge by edge, which "joins" are derived-equivalence proofs and which are not. Dormant since E-161/E-184 only because it was never made a function.

**Weakest claim the workshop relies on.** "25 of 25 joined" rests on the J = 0 premise (every gate-admitted J = 0 step is a derived equivalence), which E-282 supports only by showing the J != 0 edges are silting not tilting; the premise itself is untested by any search in this repo. Second weakest: my own E-287 provenance: generation (Bruning) is from memory, and "exchange" is not shown to be Aihara-Iyama mutation.

**Closest to a finding.** E-287: shortest tilting-exchange route B_n -> A_n is 3 (n = 5) / 4 (n = 6), End(T) = A_n verified for 7/7 and 9/9 hits, depth-1 control equals E-282's 8 vertex-local exchanges, and the object list is complete (count = D_n roots, 20/30). The one check it still needs: confirm generation, i.e. that each hit T is actually a tilting generator, by an independent route (e.g. check End(T) = A_n implies via Rickard that D^b(A_n) ~ D^b(B_n); equivalently compute the Cartan/Coxeter transformation of T-induced map and compare to the known Phi of B_n, or have scholar read Bruning). Then also n = 7 to see whether depth 4 is a law or an n = 6 accident.

**What I need.**
- theorist: End of each intermediate T on the depth-4 path (LNA or not, which step leaves the gate), and Fosse's explicit T to compare with mine.
- scholar: the Bruning classification statement (rigid + n summands + K0 unimodular => generates, D-type) and the Aihara-Iyama silting-order criterion for "exchange = mutation", with theorem numbers.
- skeptic: an independent re-derivation of the 30-object list at n = 6 (not by my cone closure), and a referee of the n = 7 run when it exists.
- experimentalist: none beyond the n = 7..9 job slots if I need > 3 processes.
