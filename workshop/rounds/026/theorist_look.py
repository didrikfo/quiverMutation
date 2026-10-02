import sys, json
D = json.load(open(sys.argv[1])); want = sys.argv[2] == '1'; lim = int(sys.argv[3])
k = 0
for r in D['recs']:
    v = r['v']; rels = r['rels']; outs = r['outs']
    per = [[rel for rel in rels if any(len(q) >= 2 and q[-1] == b and q[-2] == v for q in rel)] for b in outs]
    if not all(per) or bool(r['kerdim']) != want: continue
    if sorted(tuple(sorted(len(x) for x in p)) for p in per) != [(1, 2), (2,)]: continue
    s = lambda rel: ' = '.join('-'.join(map(str, q)) for q in rel) + ('' if len(rel) > 1 else ' = 0')
    print('v', v, 'outs', outs, 'kd', r['kerdim'], '|', ' ; '.join(s(x) for x in rels), '| arrows', r['arrows'])
    k += 1
    if k >= lim: break
