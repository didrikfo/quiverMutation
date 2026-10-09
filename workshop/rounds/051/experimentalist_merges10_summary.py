"""Summarise workshop/rounds/051/experimentalist_merges10.jsonl (merges.py 10 --witness run)."""
import json, collections, sys
p = 'workshop/rounds/051/experimentalist_merges10.jsonl'
by = collections.defaultdict(lambda: [0, 0, 0, 0, 0.0])
members = collections.defaultdict(set)
for l in open(p):
    try: r = json.loads(l)
    except ValueError: continue
    b = by[r['depth']]; b[0] += 1
    b[1] += bool(r['reachedOrbits']); b[2] += bool(r['coveredOrbits']); b[3] += bool(r['alarms'])
    b[4] += r['seconds']
    members[r['depth']].add(tuple(r['member']))
for d in sorted(by):
    print('depth', d, 'searches', by[d][0], 'distinct members', len(members[d]), 'with link', by[d][1],
          'covered', by[d][2], 'alarms', by[d][3], 'cpu-s', round(by[d][4]))
