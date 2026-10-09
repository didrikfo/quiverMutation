"""Relation-bearing control for the HH code (round 055).  Uses hh(n, rl, ok) of workshop/rounds/054/maverick_pq.py, where `ok` is the set of
pairs (a,b) joined by a nonzero path.  Algebras here are schurian (one path per pair, so a path is its endpoints).
 R: quiver 1,2 -> 3,4 -> 5 (crown plus a sink), all length-2 paths zero (rad^2 = 0, 6 arrows, 5 vertices, b1(Q) = 2).
    Independent value: no parallel pair for paths of length >= 2, so HH^m = 0 for m >= 2; HH^1 = outer derivations = arrow scalings
    mod inner = |Q1| - |Q0| + 1 = 2 (Cibils, rad^2 = 0); HH^0 = 1.  Expect (1,2).
 I: same quiver, incidence algebra (path 1->3->5 nonzero, both routes equal): poset 1,2<3,4<5 is a cone, contractible: expect (1,).
 so the only difference between R and I is the relation, and the code must see it.
 T: crown with a tail 1,2 -> 3,4 -> 5 -> 6, rad^2 = 0, b1 = 7 - 6 + 1 = 2... (arrows 4+2+1 = 7, vertices 6): expect (1,2).
 L: the same R relations on a LINEAR quiver 1->..->5 (an LNA, rl = 2222): expect (1,)  (2312.14699).
usage: python workshop/rounds/055/maverick_control.py"""
import sys
sys.path.insert(0, "workshop/rounds/054"); sys.path.insert(0, "workshop/rounds/029")
import maverick_pq as mp
def rad2(arrows, n):
    return {(a, a) for a in range(1, n + 1)} | set(arrows)
def incl(arrows, n):
    ok = rad2(arrows, n); ch = True
    while ch:
        ch = False
        for (a, b) in list(ok):
            for (c, d) in list(ok):
                if b == c and a != b and c != d and (a, d) not in ok: ok.add((a, d)); ch = True
    return ok
A = [(1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 5)]
print("R rad^2=0 crown+sink  ", mp.hh(5, None, rad2(A, 5)), "expect (1, 2)")
print("I incidence, same quiver", mp.hh(5, None, incl(A, 5)), "expect (1,)")
print("T rad^2=0 crown+tail   ", mp.hh(6, None, rad2(A + [(5, 6)], 6)), "expect (1, 2)")
print("L linear rad^2=0, n=5  ", mp.hh(5, (2, 2, 2, 2, 0)), "expect (1,)")
