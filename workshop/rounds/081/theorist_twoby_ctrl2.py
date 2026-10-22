"""Touching-run controls: run-of-three bystander plus a second bystander sharing an arrow with it. No args."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import theorist_twoby as t
if __name__ == "__main__":
    cs = [((7,3),(5,3)), ((10,3),(12,2))]
    with Pool(2) as p:
        for r in p.imap_unordered(t.run, [(a, b, 3, 2) for a, b in cs]): print(*r, flush=True)
