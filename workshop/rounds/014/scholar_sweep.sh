#!/bin/sh
# usage: scholar_sweep.sh "n:class n:class ..."  -- stop-on-reject walks, 540 s budget each, 4 at a time
for nc in $1; do echo $nc; done | xargs -P4 -I{} sh -c 'n=${1%%:*}; c=${1##*:}; timeout 10m .venv/bin/python workshop/rounds/014/scholar_walk.py $n --class $c --stop-on-reject --budget-sec 540 > workshop/rounds/014/scholar_walk_n${n}_c${c}.txt 2>&1' _ {}
