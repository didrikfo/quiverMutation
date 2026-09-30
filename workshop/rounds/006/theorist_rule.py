"""Which table rule moves 33x@o to 33(x-1)@(o+1)?  Prints every rule of ALL_MOVES whose LHS is the window 3 3 x and whose RHS is 3 3 x-1 shifted.
  .venv/bin/python workshop/rounds/006/theorist_rule.py
"""
from quivermutation import lnaMoves
for d in lnaMoves.ALL_MOVES:
    w, lhs, rhs = d[0], d[1], d[2]
    if len(lhs) == 3 and [a for _, a in lhs][:2] == [3, 3] and len(rhs) == 3 and [a for _, a in rhs][:2] == [3, 3]:
        print(d)
