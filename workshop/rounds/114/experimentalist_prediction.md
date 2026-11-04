# Prediction (written before any tail-end run, round 114)
Tail-end deletion (vertex n removed directly; image of a tail-free row r is r[:-1], checked against pwh.removeVertex(n,r,n)).
Class = rows whose Coxeter key equals the lone-3 key; K = number of end vertices no relation ends at.
P1: n=17 (K1=5): per-K image-key counts identical to head table: K=5 {412,265}, K=6 {136,129}, K=7 {65}; class rows 1007; 2 keys for K>=5,6; 1 key for K>=7.
P2: n=19 (K1=7): class rows 651, K=7 {265,257}, K=8 {129}; threshold K0=7.
Refuted by any different count (esp. a split at K=K0+1 or no split at K0).
P3 (if time): n=15, other lone-3 class h=4 ((4,7)/(7,4)): tail table equals the head table cell for cell (no prediction of the numbers themselves).
