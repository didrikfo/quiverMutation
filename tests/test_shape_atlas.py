"""The shape atlas: the census walk, the task, and the analysis over a ledger.

Checked at n = 5 and 6, where every LNA is in a quipu class and a depth-2 walk
takes a fraction of a second, so the instrument can be pinned against things
already known before it is run where nothing is.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import io
import json

import polars as pl
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


@pytest.fixture(scope = "module")
def atlas5():
    """Every LNA of length 5, walked to depth 2, resolved."""
    records = [{'unit': row, 'result': sa.walkStart(5, row, 2)} for row in sa.startsFor(5)]
    return {'records': records, 'tables': sa.resolve(records, 5)}


def test_resolve_gives_one_row_per_id_and_keys_at_every_level(atlas5):
    nodes = atlas5['tables']['nodes']
    assert nodes['id'].n_unique() == nodes.height
    for level in (0, 1, 2, 3):
        assert nodes['key{0}'.format(level)].null_count() == 0
    visits = atlas5['tables']['visits']
    assert set(visits['id'].to_list()) <= set(nodes['id'].to_list())


def test_coarser_levels_never_split_what_finer_levels_join(atlas5):
    nodes = atlas5['tables']['nodes']
    for finer, coarser in ((3, 2), (2, 1), (1, 0)):
        grouped = nodes.group_by('key{0}'.format(finer)).agg(
            pl.col('key{0}'.format(coarser)).n_unique().alias('n'))
        assert grouped['n'].max() == 1


def test_every_start_has_exactly_one_depth_zero_visit_and_it_is_a_line(atlas5):
    tables = atlas5['tables']
    roots = tables['visits'].filter(pl.col('depth') == 0).join(
        tables['nodes'].select('id', 'isLine'), on = 'id')
    assert roots['start'].n_unique() == roots.height == 14
    assert roots['isLine'].all()


def test_a_line_is_reached_from_its_own_class_only(atlas5):
    measures = sa.shapeMeasures(atlas5['tables'], 3)
    lines = measures.filter(pl.col('isLine'))
    assert lines.height > 0
    assert lines['classes'].max() == 1


def test_measures_are_in_range(atlas5):
    measures = sa.shapeMeasures(atlas5['tables'], 2)
    assert measures['returnRate'].min() >= 0 and measures['returnRate'].max() <= 1
    assert measures['starts'].min() >= 1
    assert (measures['classes'] <= measures['orbits']).all()


def test_transitions_start_and_end_at_a_line(atlas5):
    cycles = sa.transitions(atlas5['tables'], level = 2)
    assert cycles
    for count, shapes in cycles:
        assert count >= 1 and 1 <= len(shapes) <= 3 and 'LINE' not in shapes


def test_report_prints_the_sections(atlas5):
    out = io.StringIO()
    sa.report(atlas5['tables'], 2, 5, out)
    text = out.getvalue()
    assert 'hubs' in text and 'bridges' in text and 'line -> ' in text
