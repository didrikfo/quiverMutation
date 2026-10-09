import sys,pickle
sys.argv=['x','/tmp/tsm/c1.pkl','1','5','6','400000','none']
src=open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("t0 = time.time()")[0]
exec(compile(src,'tp','exec'))
recs=[r for r in pickle.load(open('/tmp/tsm/c1.pkl','rb')) if r['kind']=='fail']
k=[fingerprint.canonicalKey(recs[i]['childObj']) for i in (14,15)]
print(len(recs),k[0]==k[1],k[0] is not None, [ (recs[i]['depth'],recs[i]['v']) for i in (14,15)])
