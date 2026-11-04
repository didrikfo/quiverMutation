# Prediction (round 115, written before any run)
Method: round-114 script (key level, head ends, K = free head stretch, class = rows with the lone-3 key; K0 = last K with 2 image keys on K >= K1 rows).
P1: n=19 h=5 (class (5,10)/(10,5)): K0 = 5.   P2: n=19 h=6: K0 = 6.   P3: n=15 h=3 (class (3,8)): K0 = 3.
Falsified by K0 != min(h, n-4-h) in any of the three (a split at K = K0+1, or no split at K0). Caps: K1 chosen so K1 <= K0-1 if affordable, else stated.
P4 (non-lone, if time): n=13, class of the row with two 3s: r = lone-3 at h=2 plus a second 3 at h=7 (entries at positions 2 and 7, rest 0; check it is a valid row). Define class = all rows with the same Coxeter key. Prediction: image keys split for K <= some K0' <= 2 (the shorter free head stretch of the first 3 is 2), 1 key above; I predict K0' = 2, with counts not matching the lone-3 sequence 780/265/129/65. Falsifier: K0' != 2, or the lone-3 counts reappearing.
