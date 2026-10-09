"""Summarise experimentalist_orbscan4_n*.txt: merged words, the big orbit per n.  .venv/bin/python workshop/rounds/014/experimentalist_orbstats4.py"""
from collections import defaultdict
for n in (12,13,14,15):
    rows=[l.split() for l in open("workshop/rounds/014/experimentalist_orbscan4_n%d.txt"%n) if not l.startswith("done")]
    rows=[r for r in rows if r[2]]
    allclosed=all(r[6]=="1" for r in rows)
    mer=[r for r in rows if r[2]=="merged" and int(r[1])>1]
    orb=defaultdict(list)
    for r in mer: orb[(r[4],r[5])].append(r[0])
    print("n=%d words=%d merged=%d orbits_with_merged=%d capped=%s"%(n,len(rows),len(mer),len(orb),not allclosed))
    for k,v in sorted(orb.items(),key=lambda kv:-len(kv[1])): print("   orbit size %s: %d words %s"%(k[1],len(v)," ".join(v)))
    print("   34-containing merged:",[r[0] for r in mer if "34" in r[0]], " of ",[r[0] for r in rows if "34" in r[0]])
