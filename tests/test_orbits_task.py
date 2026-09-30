"""The orbit report of E-052: the orbits of a core's placements, reduced walk.

`45` at n = 13 is recorded in workshop round 002 (experimentalist): orbits
{0,5} 2386, {1,4} 1127, {2,3} 4217, {6} 447, each holding exactly its own
offsets' mirrors.  The smaller lengths are cheap companions.
"""

import batch


def test_45_at_13_is_the_recorded_report():
    record = batch.orbitCensus(13, '45', 1500000)
    assert record['offsets'] == list(range(7))
    got = [(o['held'], o['size'], o['closed']) for o in record['orbits']]
    assert got == [([0, 5], 2386, True), ([1, 4], 1127, True),
                   ([2, 3], 4217, True), ([6], 447, True)]
    assert all(o['mirrors'] == o['held'] for o in record['orbits'])


def test_344_at_13_has_cross_mirrors():
    record = batch.orbitCensus(13, '344', 1500000)
    got = [(o['held'], o['size'], o['closed']) for o in record['orbits']]
    assert got == [([0, 6], 42, True), ([1, 5], 34, True), ([2], 50, True),
                   ([3], 19, True), ([4], 50, True)]
    by_held = {tuple(o['held']): o['mirrors'] for o in record['orbits']}
    assert by_held[(2,)] == [4] and by_held[(4,)] == [2]


def test_a_core_that_does_not_fit_has_no_report():
    assert batch.orbitCensus(4, '45', 1000) is None


def test_the_task_resumes_from_its_ledger(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = ["orbits", "9", "--cores", "45,54", "--orbit-limit", "100000"]
    assert batch.main(args) == 0
    ledger = tmp_path / "logs" / "orbits-n9-w4a6-o100000.jsonl"
    first = ledger.read_text()
    assert batch.main(args) == 0
    assert ledger.read_text() == first
    assert batch.main(args + ["--summary"]) == 0


def test_a_partial_ledger_is_completed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    base = ["orbits", "9", "--orbit-limit", "100000"]
    assert batch.main(base + ["--cores", "45"]) == 0
    ledger = tmp_path / "logs" / "orbits-n9-w4a6-o100000.jsonl"
    first = ledger.read_text().splitlines()
    assert len(first) == 1
    assert batch.main(base + ["--cores", "45,54"]) == 0
    second = ledger.read_text().splitlines()
    assert len(second) == 2 and second[0] == first[0]
