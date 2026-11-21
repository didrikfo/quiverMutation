# Run a round-134 scholar script at another prime: scholar_prime.py p SCRIPT args...  (after the referee's scratch prime.py)
import sys, os, runpy
R = "/home/user/quiverMutation/workshop/rounds"
sys.path.insert(0, R + "/134"); sys.path.insert(0, R + "/089"); sys.path.insert(0, R + "/087")
p = int(sys.argv.pop(1))
import scholar_rankgrowth; scholar_rankgrowth.P = p
import maverick_bimod; maverick_bimod.P = p
import scholar_simple_nu; scholar_simple_nu.P = p
script = sys.argv.pop(1)
sys.argv[0] = script
g = runpy.run_path(script, run_name="__main__")
