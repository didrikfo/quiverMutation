"""Positive controls: run-of-three bystander (7,3) or (10,3) plus a second bystander. Usage: no args."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import theorist_twoby as t
if __name__ == "__main__":
    cs = [((7,3),(12,2)), ((7,3),(13,2)), ((10,3),(5,2)), ((10,3),(6,2))]
    with Pool(4) as p:
        for r in p.imap_unordered(t.run, [(a, b, 3, 2) for a, b in cs]): print(*r, flush=True)
