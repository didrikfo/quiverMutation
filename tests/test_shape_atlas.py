"""The shape atlas: the census walk, the task, and the analysis over a ledger.

Checked at n = 5 and 6, where every LNA is in a quipu class and a depth-2 walk
takes a fraction of a second, so the instrument can be pinned against things
already known before it is run where nothing is.  Spec:
docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import collections
import io
import json

import polars as pl
import pytest

import atlas
import batch
from quivermutation import atlasPage
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


def test_no_shape_is_shared_by_two_classes_at_n5(atlas5):
    """Every class at n = 5 is known and distinct; a shared L3 key would be a bug."""
    assert sa.candidateMerges(atlas5['tables']) == []


def _sharedNonLineKey(tables):
    visits = tables['visits'].filter(pl.col('depth') > 0).join(
        tables['nodes'].select('id', 'key3', 'isLine'), on = 'id').filter(~pl.col('isLine'))
    shared = visits.group_by('key3').agg(pl.col('start').n_unique().alias('n')).filter(
        pl.col('n') > 1).sort('key3')
    key = shared['key3'][0]
    rows = (visits.filter(pl.col('key3') == key).sort('depth', 'start')
            .unique('start', keep = 'first', maintain_order = True).head(2))
    return key, [sa._side(row) for row in rows.iter_rows(named = True)]


def test_replay_confirms_a_real_shared_key(atlas5):
    key, (first, second) = _sharedNonLineKey(atlas5['tables'])
    result = sa.replay(5, {'key': key, 'first': first, 'second': second})
    assert result['ok'], result


def test_replay_refuses_a_false_one(atlas5):
    key, (first, second) = _sharedNonLineKey(atlas5['tables'])
    result = sa.replay(5, {'key': key, 'first': first, 'second': dict(second, path = [])})
    assert not result['ok']


def test_validation_at_n5(atlas5):
    result = sa.validate(atlas5['tables'], 5, atlas5['records'], coverageSample = 5)
    assert result['certificates']['ok'], result['certificates']
    assert result['coverage']['ok'], result['coverage']
    assert result['quipuHub'] == 'skipped'
    assert 'shortSides' in result['squares']


def test_squares_verdict_is_ok_when_short_side_two_merely_dominates():
    # E-053 (n = 8): 957 returning squares with a short side of 2 against 36
    # with a short side of 3 -- H-022's "never a 3" is refuted, but F-027's
    # dominance of 2 still holds, so this must be ok.
    ok, majority = sa._squaresVerdict(collections.Counter({2: 957, 3: 36}), '2x2')
    assert ok and majority == 2


def test_squares_verdict_rejects_a_near_tie():
    ok, majority = sa._squaresVerdict(collections.Counter({2: 3, 3: 5}), '2x2')
    assert not ok and majority == 3


def test_squares_verdict_rejects_a_top_square_with_the_wrong_short_side():
    ok, majority = sa._squaresVerdict(collections.Counter({2: 5}), '3x3')
    assert not ok and majority == 2


def test_squares_verdict_rejects_no_squares_at_all():
    ok, majority = sa._squaresVerdict(collections.Counter(), None)
    assert not ok and majority is None


def test_a_drawing_has_a_node_per_vertex_and_a_path_per_arrow():
    data = sk.serialise(sk.deserialise(json.loads(json.dumps(
        sk.serialise(nk.LinearNakayamaAlgebra(5, "300"))))))
    svg = atlasPage.drawQuiver(data)
    assert svg.startswith('<svg') and svg.count('<circle') == 5 and svg.count('<path') >= 4


def test_the_page_carries_every_section(atlas5):
    page = atlasPage.render('Shape atlas n = 5', atlasPage.sectionsFrom(atlas5['tables'], 2, 3))
    assert page.startswith('<!doctype html>') and '<title>Shape atlas n = 5</title>' in page
    assert page.count('<svg') >= 3 and 'prefers-color-scheme: dark' in page


def test_the_command_line_reads_a_ledger(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(length = 5, depth = 2, sample = 0, seed = 0)
    jobs.runTask(batch.TASKS['atlas'], args, out = io.StringIO())
    assert atlas.main(['5', '--depth', '2', '--validate', '--page', 'page.html', '--top', '3']) == 0
    printed = capsys.readouterr().out
    assert 'hubs' in printed and 'candidate merges: 0' in printed and 'validation' in printed
    assert (tmp_path / 'page.html').exists()
    assert (tmp_path / 'logs' / 'atlas-n5-d2-s0-r0.nodes.parquet').exists()
