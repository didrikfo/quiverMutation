"""Round 119 response: traceback of the k=5 ValueError, and an uncapped timing of one k=5 cell.
Usage: experimentalist_trace.py STEPS MARGIN "s,l:s,l"  (no cap; run under timeout)"""
import sys, time, traceback, importlib.util
spec = importlib.util.spec_from_file_location("chain", "workshop/rounds/107/experimentalist_chain.py")
ch = importlib.util.module_from_spec(spec); spec.loader.exec_module(ch)
steps, margin = int(sys.argv[1]), int(sys.argv[2])
pat = tuple(tuple(map(int, r.split(","))) for r in sys.argv[3].split(":"))
t0 = time.time()
try: print(*ch.run((pat, steps, margin)), flush=True)
except Exception: traceback.print_exc(file=sys.stdout); print("secs", round(time.time()-t0,1))
