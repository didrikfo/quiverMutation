# Theorist notebook (rewritten round 006)

## What I believe now
- k(33x) = 2x is explained (round 006): in the interior `33x` is moved by the double mutation of F-032, NOT by the rule table (floating table is silent for
  x >= 5; it has 334<->333 only). L on the middle relation gives 33x@o -> 33(x-1)@(o+1), so c = x+o is conserved; the self-dual seed 333 joins chain c to the
  mirror of chain n-c; hence o <-> n-2x-o, and offsets o > s = n-2x (x-3 of them) have no room for a partner. Exact membership test 135/135 at n = 14,15,16.
- General formula for a drift family with self-dual seed x0, footprint x+w0: k = 2x + w0 - x0.
- Unproved: the end link (mirror chain D(y) -> 33 chain top, anchored/edge moves) and the upper bound (orbit has nothing else). Only computed.
- 44x has the same drift and seed but the prediction k = 2x-1 is FALSE: all offsets merge into the 333@0 orbit at n = 14..16. So the chain is a lower bound;
  exactness of 33x is a fact about what else attaches, not about the argument.
- 34x: k = x+3 = w, d = 0 (o <-> hi-o) for x = 5,7,8 and 344 at even n; 346 merges. Only 345 is explained (334@o -> 4444@o <- 345@o, so it rides the 33x chain).
- 45x: no drift; 456, 457 merge, 455 parity {0,2,4,6,8}{1,3,5,7} at 16 (4046-type, unexplained).
- Older (round 004): H-021' stands; interior/end-touch split explains none of the 7 failures; "s in cons" is tautological, do not sell it.

## What I tried
- Hand computation of L_t on 33x; theorist_{rule,step,local,path,chain,nbrs,words,link}.py in rounds/006; p-value-free exact-set test, not a fit.
- Did not try: x >= 9 (script cap: "3310" parses as 3,3,1,0), n >= 17, proving the end link, 34x at x = 7 by drift (no 33y row in its orbit, so another chain).

## Next
- Prove the end link and closure for one (n,c) by enumerating all moves out of the 310-row orbit (n = 16, x = 5, o = 4), then state it as a lemma.
- Find the drift family for 347/348/344 (orbit words: 3,y,y+1 at fixed offset is a second conserved direction: 345,356,367,378@3): same method, conserved quantity = offset.
- Ask why 44x and 346 fall into the c = 3 orbit (333@0 touches the source): is it always the source-touching chain top?
- Blind spots: the double mutation is interior-uniform by construction, so "H-020 theorem" (T4) is really a statement about doubleMutation, plus anchors at the ends; do not call it the rule table.
