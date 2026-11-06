"""T7/H-010 round 119: gap-2 pair shapes at k>=4, per-cell wall cap, exceptions recorded (not fatal).
Usage: experimentalist_gap2.py STEPS MARGIN CAPSECS "s,l:s,l;..."  (wraps rounds/107 experimentalist_chain.run)"""
import sys, signal, importlib.util, time
spec = importlib.util.spec_from_file_location("chain", "workshop/rounds/107/experimentalist_chain.py")
ch = importlib.util.module_from_spec(spec); spec.loader.exec_module(ch)
def h(*a): raise TimeoutError
signal.signal(signal.SIGALRM, h)
steps, margin, cap = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
for c in sys.argv[4].split(";"):
    pat = tuple(tuple(map(int, r.split(","))) for r in c.split(":"))
    signal.alarm(cap)
    try: print(*ch.run((pat, steps, margin)), flush=True)
    except TimeoutError: print(pat, "CAP", cap, "s: undecided (not a verdict)", flush=True)
    except Exception as e: print(pat, "ERROR", type(e).__name__, str(e)[:90], flush=True)
    signal.alarm(0)
