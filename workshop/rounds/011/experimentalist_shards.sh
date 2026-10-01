#!/bin/bash
# usage: experimentalist_shards.sh I ; one K=4 depth-6 shard, output saved
I=$1
s=$(date +%s)
timeout 10m .venv/bin/python workshop/rounds/010/toolsmith_verify.py 9 6 4 -1 --cand $I > workshop/rounds/011/experimentalist_cand$I.txt 2>&1
rc=$?
echo "cand $I rc $rc seconds $(( $(date +%s)-s ))" >> workshop/rounds/011/experimentalist_shard_times.txt
