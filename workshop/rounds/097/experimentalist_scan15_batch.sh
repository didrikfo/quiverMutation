#!/bin/bash
# usage (repo root): experimentalist_scan15_batch.sh  -- runs undone shards of 96, 4 workers, stops after 9 min wall; rerun until all 96 done
cd "$(dirname "$0")/../../.."
export D=workshop/rounds/097
for w in 0 1 2 3; do export w
 timeout 540 bash -c 'for s in $(seq $w 4 95); do
     [ -f $D/experimentalist_scan15_hits_7_$s.txt ] && continue
     timeout 400 .venv/bin/python $D/experimentalist_scan15.py scan $s 96 >> $D/experimentalist_scan15_log.txt
   done' &
done
wait
ls $D/experimentalist_scan15_hits_7_*.txt | wc -l
