"""End-to-end runs of the classification pipeline for small lengths.

mutationSearch writes its CSV and per-class text files into the current working
directory, so each test runs in its own tmp_path.
"""

import collections
import csv

import pytest

import quiverMutation as qm
from helpers import quiet
from paper_classification import PAPER_CLASSES, relation_string


def run_search(tmp_path, monkeypatch, length, depth):
    monkeypatch.chdir(tmp_path)
    quiet(qm.mutationSearch, length, depth, 0, createNewCSVfile=True)
    with open(tmp_path / f"A_{length}_mutation_classes.csv", newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    assert header[:2] == ["Relations", "Mutation class"]
    return data


def classes_by_coxeter_polynomial(rows):
    """Group the search's classes by Coxeter polynomial.

    The search finds a lower bound on each class: it only merges two LNAs when
    its depth-first search actually walks a mutation path between them.  Classes
    that share a Coxeter polynomial are the candidates for merging, which is the
    step that was done by hand.
    """
    by_poly = collections.defaultdict(set)
    members = collections.defaultdict(list)
    for relations, class_name, _path, poly, _numbering in rows:
        by_poly[poly].add(class_name)
        members[poly].append(relations)
    return by_poly, members


def test_length_5_classification(tmp_path, monkeypatch):
    rows = run_search(tmp_path, monkeypatch, 5, 4)
    assert len(rows) == 14  # Catalan(4)

    by_poly, members = classes_by_coxeter_polynomial(rows)
    # The search finds both classes of A_5 outright, with no merging needed.
    assert len(by_poly) == 2
    assert all(len(names) == 1 for names in by_poly.values())
    assert sorted(len(m) for m in members.values()) == [6, 8]


@pytest.mark.slow
def test_length_6_classification(tmp_path, monkeypatch):
    rows = run_search(tmp_path, monkeypatch, 6, 6)
    assert len(rows) == 42  # Catalan(5)

    by_poly, members = classes_by_coxeter_polynomial(rows)
    # Four classes for n = 6 -- A_6, D_6, E_6 and the extended Dynkin D~_5 --
    # matching the four quipus of order 6 in the paper's table.
    assert len(by_poly) == len(PAPER_CLASSES[6]) == 4
    # The search leaves D_6 split in two: it does not find a mutation path
    # between A_{6,(1)}^{(3)} and A_{6,(3)}^{(3)} at this depth.
    assert sorted(len(names) for names in by_poly.values()) == [1, 1, 1, 2]
    assert sorted(len(m) for m in members.values()) == [1, 12, 13, 16]


@pytest.mark.slow
def test_length_7_classification(tmp_path, monkeypatch):
    """n = 7 is the first length where merging by hand does real work.

    The search leaves 11 classes, which fall into 6 groups by Coxeter
    polynomial -- exactly the 6 quipus of order 7 in the paper's table.
    """
    rows = run_search(tmp_path, monkeypatch, 7, 6)
    assert len(rows) == 132  # Catalan(6)

    by_poly, members = classes_by_coxeter_polynomial(rows)
    assert sum(len(names) for names in by_poly.values()) == 11
    assert len(by_poly) == len(PAPER_CLASSES[7]) == 6
    assert sorted(len(m) for m in members.values()) == [4, 6, 7, 29, 32, 54]


@pytest.mark.slow
def test_length_8_classification(tmp_path, monkeypatch):
    """The full n = 8 classification, the largest one the paper prints.

    The search leaves 28 classes, which fall into 11 groups by Coxeter
    polynomial -- exactly the 11 quipus of order 8 in the paper's table, with
    every published class landing inside a single group.
    """
    rows = run_search(tmp_path, monkeypatch, 8, 6)
    assert len(rows) == 429  # Catalan(7)

    by_poly, members = classes_by_coxeter_polynomial(rows)
    assert sum(len(names) for names in by_poly.values()) == 28
    assert len(by_poly) == len(PAPER_CLASSES[8]) == 11
    assert sorted(len(m) for m in members.values()) == [
        1, 4, 9, 10, 13, 26, 40, 64, 64, 65, 133,
    ]
    assert sum(len(m) for m in members.values()) == 429

    poly_of = {relations: poly for relations, _c, _p, poly, _n in rows}
    for label, entries in PAPER_CLASSES[8].items():
        polys = {poly_of[relation_string(8, *entry)] for entry in entries}
        assert len(polys) == 1, f"class {label} split across {polys}"


@pytest.mark.parametrize("length", [4, 5, 6])
def test_search_partition_refines_the_published_one(length, tmp_path, monkeypatch):
    """Every LNA the paper puts in one class must land in one Coxeter group.

    The search may split a published class (an unfound mutation path), but it
    must never merge two of them, and it must never disagree with the paper
    about which LNAs belong together.
    """
    rows = run_search(tmp_path, monkeypatch, length, 5)
    poly_of = {relations: poly for relations, _c, _p, poly, _n in rows}

    for label, entries in PAPER_CLASSES[length].items():
        polys = set()
        for starts, lengths in entries:
            key = relation_string(length, starts, lengths)
            assert key in poly_of, f"n={length}: {label} member {key!r} missing from CSV"
            polys.add(poly_of[key])
        assert len(polys) == 1, f"n={length} class {label} split across {polys}"
