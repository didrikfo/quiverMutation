"""Round 085 (toolsmith): pairwise isomorphism (quiver+relations, vertex-relabelling) of the children of the key-kept perI steps.
usage: toolsmith_distinct.py PKL1 PKL2 ...   (pickles from toolsmith_d5pkl.py; first listed are compared pairwise)  Run from the repository root."""
import sys, pickle
files = sys.argv[1:]; sys.argv = ['x', 'none']
src = open('workshop/rounds/082/skeptic_join.py').read().split("def expand")[0]
sys.path.insert(0, '.'); exec(compile(src, 'j', 'exec'))
ch = []
for f in files:
    x = pickle.load(open(f, 'rb'))[0]; c = child(x['parent'], x['kind'], x['v']); ch.append(c)
sg = [sig(c) for c in ch]
print('files', [f.split('/')[-1] for f in files])
print('WL hashes distinct:', len(set(sg)), 'of', len(ch))
cls = []
for i, c in enumerate(ch):
    for cl in cls:
        if sg[cl[0]] == sg[i] and iso(ch[cl[0]], c): cl.append(i); break
    else: cls.append([i])
print('iso classes of children (up to vertex relabelling):', [[files[i].split('/')[-1] for i in cl] for cl in cls])
print('number of classes', len(cls))
print('quiver arrow counts', [c.quiver.number_of_edges() for c in ch], 'relations', [len(procedure.relationsFrom(c)) for c in ch])
print('parents iso classes too:')
ps = [x for x in (pickle.load(open(f, 'rb'))[0]['parent'] for f in files)]
pc = []
for i, c in enumerate(ps):
    for cl in pc:
        if iso(ps[cl[0]], c): cl.append(i); break
    else: pc.append([i])
print('parent classes', len(pc))
