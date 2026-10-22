"""The three m=3,m=3 cells of theorist_twoby.py 3 2 6 12 that did not finish in 10 min with 4 jobs.
Usage: theorist_twoby_rest.py"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
import theorist_twoby as t
if __name__ == "__main__":
    done = {((6,2),(7,2)),((7,2),(11,2)),((6,2),(11,2)),((11,2),(12,2)),((6,2),(12,2)),((7,2),(12,2)),((7,2),(11,3)),((6,2),(11,3)),((6,3),(7,2)),((7,2),(12,3)),((11,2),(12,3)),((11,3),(12,2)),((6,3),(11,2)),((6,2),(12,3)),((6,3),(12,2))}
    todo = [c for c in t.cells(6, 12) if c not in done]
    print(todo, flush=True)
    with Pool(len(todo)) as p:
        for r in p.imap_unordered(t.run, [(a, b, 3, 2) for a, b in todo]): print(*r, flush=True)
