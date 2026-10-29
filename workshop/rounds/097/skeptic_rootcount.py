import sys, random, warnings
sys.path.insert(0,'.'); sys.path.insert(0,'workshop/rounds/097'); warnings.simplefilter('ignore')
from skeptic_omega import *
for (a,b,pat2) in ((4,6,False),(5,6,True)):
    rng=random.Random(97); tot=0; got=0
    for n in range(a,b+1):
        for _ in range(300):
            tot+=1; got += gen_parent(rng,n,pat2=pat2) is not None
    print(a,b,pat2,'attempts',tot,'roots',got)
