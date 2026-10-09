#!/bin/bash
# usage: toolsmith_shards.sh tag size total  : run all shards of a frontier, 3 at a time, each under timeout 10m (repo root)
tag=$1; sz=$2; tot=$3; lo=0
while [ $lo -lt $tot ]; do
  hi=$((lo+sz)); [ $hi -gt $tot ] && hi=$tot
  while [ $(jobs -r | wc -l) -ge ${JOBS:-3} ]; do sleep 2; done
  [ -f /tmp/tsm/sh_${tag}_$lo.pkl ] || timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_depth7.py shard $tag $lo:$hi > /tmp/tsm/sh_${tag}_$lo.log 2>&1 &
  lo=$hi
done; wait
