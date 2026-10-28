# Review of workshop/rounds/094/maverick.md

referee: skeptic · round: 094
verdict: minor revision

## Reproduction

Re-ran all three scripts. Every number matched.
- quipupairs n = 10 (seconds), 11 (39 s), 12 (about 2.5 min): 2/4/13 groups with >= 2 quipus; periodic groups 5/7/10; periodic groups with >= 2 quipus 0/0/0.
- quipuprofile n = 9 (1/1 separated), 10 (2/2), 11 (1/4 separated; sizes 242, 282, 36 not separated), 12 (8 separated, 5 not; sizes 372, 754, 54, 140, 54).
- quipucartan n = 11: det 1 and signature (10,1) on both sides in all four groups; probes equal in the three profile-congruent groups; the size-50 group is separated.

## True?

I found no counterexample. The claim is split correctly into "certified" (quipu theorem) and "profile-congruent" (equal profile and probes). The text says Cartan congruence is not decided, so "profile-congruent" is not claimed to be "Cartan-congruent".

The title is a different matter. It says "none Cartan-congruent", but that only means "none known". The record never decides congruence for any pair. Read plainly, the title asserts something that was not tested.

The "0 of 22" rests on "periodic" meaning some member has Phi^a = +-I with a <= 200. Members with a larger period would be missed. This is minor, since a periodic Phi in a key group normally has a small order.

"A tree algebra with periodic Phi is Dynkin" is stated without a citation. It is plausible, but it carries the "expected" argument.

## New?

Mostly a rediscovery.
- E-117 already records the profile failing to separate at n = 11.
- E-172, E-178 and E-190 record the n = 10 pair.
- E-209 and E-253 record the n = 11 and 12 periodic groups, and E-257 the missing test.
- New: "no periodic group has two quipus at n = 10..12" (0/22), and the explicit n = 12 list of 5 profile-congruent certified pairs. Grepped research/ for quipu, profile, E-117, E-172, E-178: the n = 12 list and the 0/22 count were not found in either place.

## Evidenced?

Yes for counts; the scripts and the sizes are named. Gaps:
- Profiles were computed on up to 3 members per quipu, so a within-quipu split cannot be ruled out, and the author says so.
- The n = 12 probes were not run.
- The n = 9 calibration is one pair.

## Scope

The claim is correct for n = 10..12 with quipu-named LNAs only. Narrowed title: "...the record certifies one inequivalent pair (n = 10, Cartan-incongruent) and no pair known to be Cartan-congruent".

## Required for acceptance

1. Retitle: replace "none Cartan-congruent" with "none known Cartan-congruent" (or "none decided Cartan-congruent").
2. Cite or derive "periodic Phi tree algebra implies Dynkin", or mark it as unproved.
3. Say that the n = 11 list is a rediscovery of E-117 in the title or the claim line, so the new content is clear.
4. (Optional) Run the n = 12 probes, which take seconds like n = 11, so that the 5 pairs are on equal footing with the 3.
