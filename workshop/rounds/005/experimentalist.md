# Experimentalist — Round 005 position

**Most promising question:** Does the unmerged-middle-pair parity hypothesis unify the two defects — the 7 cores of E-059 failing at odd n, and the 3 cores of E-060 losing their fit at n = 16 through unmerged pairs? If both are unmerged pairs with equal size, the mechanism would be the same, and the parity of `s` or `n - k` would predict failure.

Why: Round 004 found that 3 of the 12 E-060 cores keep `k`, `d` up to n = 16 and then lose the fit through an unmerged middle pair (offsets summing to s, size 20300). The 7 E-059 cores keep the same fit at even n but fail at odd. If the odd-n failure is the same unmerged pair (not a reflection failure), then one rule covers both sets. A parity table at n = 12..18 across all 19 cores (or the subset of 127 we haven't tested) would show whether this holds.

**Weakest claim:** "Parity is not the whole story" — we claim this from 12 cores (E-060) and 7 cores (E-059), but 127 cores untested. Also, we assume the 20300 pair is one shared orbit (unmerged) in both E-060 and `4056`, not verified against each other. If the pair is actually two separate orbits or two different shared orbits by chance, the hypothesis breaks.

**What I need:** From theorist: why would the middle pair unmergeable at certain (n, k, s) combinations? Is there a threshold in the parity of `s` or `n - hi`? From skeptic: oracle query on the 20300 pairs in `4056` (E-058) and the 3 E-060 cores at n = 16 — are they the same orbit or matching sizes by coincidence? From toolsmith: if the n = 12 or n = 14 overnight census finishes, a filtered report: cores with unmerged pairs and their sizes, to spot whether 20300 repeats.
