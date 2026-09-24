"""The shape atlas: the census walk, the task, and the analysis over a ledger.

Checked at n = 5 and 6, where every LNA is in a quipu class and a depth-2 walk
takes a fraction of a second, so the instrument can be pinned against things
already known before it is run where nothing is.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import io
import json

import pytest

import batch
from quivermutation import jobs
from quivermutation import nakayama as nk
from quivermutation import search
from quivermutation import shapeAtlas as sa
from quivermutation import shapeKeys as sk


def test_a_walk_records_everything_the_label_exact_search_reaches():
    record = sa.walkStart(5, "300", 2)
    got = {search.quiverKey(sk.deserialise(node['quiver'])) for node in record['nodes']}
    expected = set(search.quiversReachedFrom(nk.LinearNakayamaAlgebra(5, "300"), 2,
                                             alsoDual = True))
    assert expected - {None} <= got


def test_the_first_node_is_the_start_and_every_edge_joins_recorded_nodes():
    record = sa.walkStart(5, "300", 2)
    first = record['nodes'][0]
    assert first['path'] == [] and first['depth'] == 0 and first['features']['isLine']
    size = len(record['nodes'])
    assert record['edges']
    for parent, vertex, child in record['edges']:
        assert 0 <= parent < size and 0 <= child < size and vertex != 0
    assert len({node['id'] for node in record['nodes']}) == size


def test_a_walk_is_json():
    record = sa.walkStart(5, "300", 2)
    assert json.loads(json.dumps(record)) == record


def test_every_n6_lna_is_tagged_with_a_quipu_class_and_mirrors_share_an_orbit():
    tags = sa.startTags(6)
    assert len(tags) == 42                                   # Catalan(5)
    assert all(cls.startswith('quipu:P^(') for _orbit, cls in tags.values())
    for row in tags:
        mirror = sa.rowString(batch.fm.mirrorRow(6, tuple(int(c) for c in row)))
        assert tags[row][0] == tags[mirror][0]


def test_sampling_keeps_every_leftover_and_one_per_quipu_class():
    tags = sa.startTags(9)
    chosen = sa.startsFor(9, sample = 1, seed = 0)
    leftovers = sorted(row for row, (_orbit, cls) in tags.items() if cls.startswith('orbit:'))
    assert len(leftovers) == 9                               # F-032, n = 9
    assert set(leftovers) <= set(chosen)
    quipuClasses = {cls for _orbit, cls in tags.values() if cls.startswith('quipu:')}
    assert len(chosen) == len(leftovers) + len(quipuClasses)
    assert sa.startsFor(9, sample = 1, seed = 0) == chosen   # deterministic


def test_the_task_writes_a_ledger_and_resumes_from_it(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    task = batch.TASKS['atlas']
    args = argparse.Namespace(length = 5, depth = 1, sample = 0, seed = 0)
    assert jobs.runTask(task, args, out = io.StringIO()) == 0
    records = jobs.Ledger(task.ledgerPath(args)).records()
    assert len(records) == 14                                # Catalan(4)
    assert records[0]['result']['nodes'][0]['depth'] == 0
    again = io.StringIO()
    assert jobs.runTask(task, args, out = again) == 0
    assert "nothing to do" in again.getvalue()
