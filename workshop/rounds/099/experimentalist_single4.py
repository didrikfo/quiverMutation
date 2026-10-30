"""Controls for pair (8:4)(10:4) at k=4 or 5, chosen single-bystander cells. Usage: STEPS MARGIN "t,m;t,m" [--jobs J]"""
import sys
sys.path.insert(0,"workshop/rounds/099")
from multiprocessing import Pool
from experimentalist_single import run
if __name__=="__main__":
    st,mg=int(sys.argv[1]),int(sys.argv[2]); cs=[tuple(map(int,c.split(","))) for c in sys.argv[3].split(";")]
    jobs=int(sys.argv[sys.argv.index("--jobs")+1]) if "--jobs" in sys.argv else 4
    with Pool(jobs) as p:
        for r in p.imap_unordered(run,[(c,st,mg) for c in cs]): print(*r,flush=True)
