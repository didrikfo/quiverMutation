"""Write compact (relation string, class) maps from the cached classifyLength CSVs (/tmp/s1del_A_<n>.csv).
usage: python workshop/rounds/065/experimentalist_dumpclasses.py   -> experimentalist_classes_n8.csv, _n9.csv
Regenerate the /tmp CSVs by running experimentalist_delcompat.py 9 (repeatedly under timeout 10m until it prints results)."""
import csv
from quivermutation import mutationClassTable as mct
for n in (8, 9):
    t = mct.MutationClassTable.fromCSV('/tmp/s1del_A_%d.csv' % n, n)
    with open('workshop/rounds/065/experimentalist_classes_n%d.csv' % n, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['relations', 'class'])
        for r in t.rows(): w.writerow([r[0], r[1]])
