#!/bin/bash
# round 050 toolsmith: resumable collector driver; usage: toolsmith_run.sh cls  (rerun until the pickle exists; each call < 10 min)
cls=$1; [ -f /tmp/tsm/c$cls.pkl ] && { echo done; exit 0; }
DEADLINE=${DEADLINE:-480} timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 $cls 20000 /tmp/tsm/c$cls.pkl 100
