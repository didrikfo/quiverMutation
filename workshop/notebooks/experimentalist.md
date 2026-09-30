# Experimentalist notebook (rewritten each round)

## What I now believe (after round 002)
- At n = 13 all 139 single-cluster cores of `--max-word 4` close under the reduced walk (395 orbit walks, 0 caps). Census, fit and raw data are committed in `workshop/rounds/002/` (`experimentalist_census.py`, `_fit.py`, `_census_n13.jsonl`).
- H-021's "mirror <=> reflection pairing" fails on every reading. Loose (some orbit holds a mirror of some placement): true of 129 cores, 20 hold one with no pairing. Strict (mirror of c@p in an orbit not holding p): 8 cores; 108 of the 109 pairing cores lack it, because the mirror acts inside the orbit. Strict2 (mirror plus another offset in the orbit): pairing => strict2 for all 109, but 7 cores (`344 366 4044 4403 4404 4405 4605`) hold it with no fit.
- Under strict, `3346` and `4056` hold no mirror and agree with "no pairing": H-021's text about them is right under strict, wrong under loose.
- Pairing itself is solid: 62 cores exact (d = 0), 47 with overhang 1..3, 30 none. Large-d fits are weak evidence.
- E-052's `45` numbers reproduce again at 13.

## What I tried
- Ran the 4-shard census in the background with `timeout 10m` per shard, about 5 min wall; the loop finished all shards in one go.
- The fit rule: smallest d in 0..6, s = lo+hi+-d, closure of every orbit under o -> s-o inside the offset set, and one real swapped pair.

## What I would do next
1. Test d(c) against `|head - tail|` of H-020 for the 109 fitting cores (unchecked, three rounds running).
2. The n = 14 census (OVERNIGHT Menu 4), then the fit again: which of the 7 strict cores and 13 mirror-without-fit cores change class with length.
3. Study the 7 `344`-type cores: pairing with a defect at the middle offsets (`2`,`4` swapped by mirror only).
4. Toolsmith: turn the census into a `batch.py` task with a test pinning `45` at 13.
- Watch: I have not checked whether "strict2" is a trivial consequence of orbits being unions of pairs; do not over-read the 132/139 agreement.
