"""classifyLength: seed, search, name, resolve -- the whole hand-merge, automated.

Each of these checks the result against the published table in
arXiv:2305.06642: the number of classes must equal the number of quipus of that
order, each published class must land in one computed class, and nothing may be
left as a candidate.
"""

import collections

import pytest

import quipuForms as qf
import quiverMutation as qm
from helpers import quiet
from paper_classification import PAPER_CLASSES, relation_string


def classify(tmp_path, monkeypatch, length, depth=6):
    monkeypatch.chdir(tmp_path)
    return quiet(qm.classifyLength, length, depth, depth, None, False)


def assert_matches_the_paper(table, report, length, expectedSizes):
    classOf = {row[0]: row[1] for row in table.rows()}

    assert not table.unassignedRelationStrings()
    assert report["candidate"] == {}, "classes left unsettled"
    assert report["separated"] == {}, "classes proved distinct but sharing a polynomial"
    assert len(table.classNames()) == len(PAPER_CLASSES[length])

    # Every published class lands in exactly one computed class.
    for label, entries in PAPER_CLASSES[length].items():
        names = {classOf[relation_string(length, *entry)] for entry in entries}
        assert len(names) == 1, f"{label} split across {names}"

    # Distinct published classes land in distinct computed classes.
    names = [
        classOf[relation_string(length, *entries[0])]
        for entries in PAPER_CLASSES[length].values()
    ]
    assert len(set(names)) == len(names)

    sizes = collections.Counter(row[1] for row in table.rows())
    assert sorted(sizes.values()) == sorted(expectedSizes)
    assert sum(sizes.values()) == len(table)


@pytest.mark.slow
def test_length_6(tmp_path, monkeypatch):
    table, report = classify(tmp_path, monkeypatch, 6)
    assert_matches_the_paper(table, report, 6, [16, 13, 12, 1])


@pytest.mark.slow
def test_length_7(tmp_path, monkeypatch):
    table, report = classify(tmp_path, monkeypatch, 7)
    assert_matches_the_paper(table, report, 7, [54, 32, 29, 7, 6, 4])


@pytest.mark.slow
def test_length_8(tmp_path, monkeypatch):
    """n = 8 needs the resolve step.

    Seeding leaves class 340030 -- A_{8,(1,2,5)}^{(3,4,3)}, whose relations
    overlap too much for the theorem -- reaching nothing classified, because the
    search only walks right mutations and nothing seeded lies downstream of it.
    Searching from its relation dual finds the link.
    """
    table, report = classify(tmp_path, monkeypatch, 8)
    assert_matches_the_paper(
        table, report, 8, [133, 65, 64, 64, 40, 26, 13, 10, 9, 4, 1])


@pytest.mark.slow
def test_class_names_are_the_quipus_they_belong_to(tmp_path, monkeypatch):
    """Seeding names each class by its quipu rather than by an arbitrary LNA."""
    table, _report = classify(tmp_path, monkeypatch, 7)
    for className in table.classNames():
        assert className.startswith("P^("), className
        # The name must parse back to a quipu of the right order.
        parameters = next(
            (k, m) for k, m in [_parse(className)]
        )
        graph = qf.graphFromQuipuParameters(*parameters)
        assert graph.number_of_nodes() == 7


def _parse(name):
    """'P^(1,2)_(1,0,1)' -> ((1,0,1), (1,2))."""
    m, k = name[len("P^("):].split(")_(")
    return tuple(int(v) for v in k.rstrip(")").split(",")), tuple(int(v) for v in m.split(","))
