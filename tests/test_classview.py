"""Reading a classification back: the query layer and the page it renders.

`classview` answers questions about a finished classification -- one row per
class, not one per LNA -- and `classpage` renders the same thing as a page.
What is worth pinning is the shape of the answers and the two places the
reading can go wrong: the *kind* a class gets from its hereditary form, and the
notation the page has to parse back out of a class name to draw its quipu.
"""

import json

import pytest

from quivermutation import classpage
from quivermutation import classview as cv
from quivermutation import mutationClassTable as mct
from quivermutation import nakayama as nk
from quivermutation import quipuForms as qf


def tableFrom(rows, length):
    """A table from (relations, class, path, coxeter, numbering, form) tuples."""
    return mct.MutationClassTable([list(row) for row in rows], length)


def classificationFrom(rows, length = 6):
    return cv.Classification.fromTable(tableFrom(rows, length))


NUMBERING = "1;2;3;4;5;6"

# Three classes of the three kinds that actually occur, plus an unplaced row.
SAMPLE = [
    ("1;2;3", "P^(2)_(1,1)", "", "poly-D", NUMBERING, "P^(2)_(1,1)"),
    ("1;2;3|3;4;5", "P^(2)_(1,1)", "4;1", "poly-D", NUMBERING, "P^(2)_(1,1)"),
    ("2;3;4", "P^(2)_(1,2)", "", "poly-E", NUMBERING, "P^(2)_(1,2)"),
    ("3345000", "3345000", "", "poly-T", NUMBERING, "C(2,4,4)"),
    ("3033030", "3033030", "", "poly-N", NUMBERING, mct.NOT_PIECEWISE_HEREDITARY),
    ("1;2;3;4;5", "", "", "", "", ""),
]


# -- what kind a class is -------------------------------------------------

@pytest.mark.parametrize("form, kind", [
    ("P^(1,4)_(1,0,1)", cv.QUIPU),
    ("P^(0)_(0,8)", cv.QUIPU),
    ("C(2,4,4)", cv.CANONICAL),
    (mct.NOT_PIECEWISE_HEREDITARY, cv.NOT_PIECEWISE_HEREDITARY),
    ("(((()())())())", cv.NON_QUIPU_TREE),
    ("wl:9fd0", cv.NOT_A_TREE),
    ("", cv.UNNAMED),
    ("   ", cv.UNNAMED),
    ("P^(2)_(1,1)|((()())())", cv.CONTRADICTORY),
])
def test_the_kind_a_hereditary_form_names(form, kind):
    assert cv.kindOfForm(form) == kind


def test_every_kind_is_listed():
    """`classes.py --kind` offers exactly these, so the list must be complete."""
    forms = ["P^(2)_(1,1)", "C(2,4)", mct.NOT_PIECEWISE_HEREDITARY, "(())", "wl:x", "",
             "a|b"]
    assert {cv.kindOfForm(f) for f in forms} == set(cv.KINDS)


# -- the summary ----------------------------------------------------------

def test_one_row_per_class_largest_first():
    classification = classificationFrom(SAMPLE)
    found = classification.classes()
    assert [c.name for c in found] == [
        "P^(2)_(1,1)", "3033030", "3345000", "P^(2)_(1,2)"]
    assert [c.size for c in found] == [2, 1, 1, 1]
    assert [c.kind for c in found] == [
        cv.QUIPU, cv.NOT_PIECEWISE_HEREDITARY, cv.CANONICAL, cv.QUIPU]


def test_the_unplaced_rows_are_not_a_class():
    """An empty class column means no search has reached the LNA yet."""
    classification = classificationFrom(SAMPLE)
    assert classification.unsettled() == ["1;2;3;4;5"]
    assert sum(c.size for c in classification.classes()) == 5
    assert classification.rowCount() == 6


def test_the_representative_is_the_member_that_needs_no_mutations():
    """The class' seed: shortest mutation path, then shortest relation string."""
    classification = classificationFrom(SAMPLE)
    assert classification.classNamed("P^(2)_(1,1)").representative == "1;2;3"


def test_members_come_back_with_their_mutation_paths_seed_first():
    classification = classificationFrom(SAMPLE)
    assert classification.members("P^(2)_(1,1)") == [
        ("1;2;3", "", NUMBERING),
        ("1;2;3|3;4;5", "4;1", NUMBERING),
    ]


def test_the_counts_by_kind_cover_every_class_and_every_lna():
    classification = classificationFrom(SAMPLE)
    counts = classification.kindCounts()
    assert counts[cv.QUIPU] == (2, 3)
    assert counts[cv.CANONICAL] == (1, 1)
    assert counts[cv.NOT_PIECEWISE_HEREDITARY] == (1, 1)
    assert sum(classes for classes, _ in counts.values()) == len(classification.classes())
    assert sum(members for _, members in counts.values()) == 5


# -- the filters the questions take ---------------------------------------

def test_the_filters():
    classification = classificationFrom(SAMPLE)
    assert [c.name for c in classification.classes(kind = cv.QUIPU)] == [
        "P^(2)_(1,1)", "P^(2)_(1,2)"]
    assert [c.name for c in classification.classes(quipu = "P^(2)_(1,2)")] == ["P^(2)_(1,2)"]
    assert [c.name for c in classification.classes(coxeterPolynomial = "poly-D")] == [
        "P^(2)_(1,1)"]
    assert [c.name for c in classification.classes(minSize = 2)] == ["P^(2)_(1,1)"]
    assert [c.name for c in classification.classes(maxSize = 1)] == [
        "3033030", "3345000", "P^(2)_(1,2)"]
    assert classification.classes(kind = cv.QUIPU, minSize = 5) == []


def test_a_shared_coxeter_polynomial_is_reported_and_a_unique_one_is_not():
    """Where the polynomial stops separating classes is the interesting question.

    Two classes carrying one polynomial is exactly the F-010 situation, and the
    view has to surface it rather than let it hide among the rows.
    """
    rows = [
        ("1;2;3", "A", "", "shared", NUMBERING, "P^(1,4)_(1,0,1)"),
        ("2;3;4", "B", "", "shared", NUMBERING, "P^(1,2)_(1,1,2)"),
        ("3;4;5", "C", "", "alone", NUMBERING, "P^(2)_(1,1)"),
    ]
    collisions = classificationFrom(rows).coxeterCollisions()
    assert list(collisions) == ["shared"]
    assert sorted(c.name for c in collisions["shared"]) == ["A", "B"]


# -- opening one off disk -------------------------------------------------

def test_a_table_round_trips_through_parquet_and_csv(tmp_path):
    table = tableFrom(SAMPLE, 6)
    table.writeCSV(str(tmp_path / "A_6_mutation_classes.csv"), header = True)
    table.writeParquet(str(tmp_path / "A_6_mutation_classes.parquet"))

    fromParquet = cv.Classification.forLength(6, str(tmp_path))
    fromCSV = cv.Classification.fromCSV(str(tmp_path / "A_6_mutation_classes.csv"), 6)
    assert [(c.name, c.size) for c in fromParquet.classes()] == \
        [(c.name, c.size) for c in fromCSV.classes()]

    (tmp_path / "A_6_mutation_classes.parquet").unlink()
    assert [c.name for c in cv.Classification.forLength(6, str(tmp_path)).classes()]


def test_a_length_that_was_never_classified_says_so(tmp_path):
    with pytest.raises(FileNotFoundError) as raised:
        cv.Classification.forLength(11, str(tmp_path))
    assert "classify.py 11" in str(raised.value)


# -- the page -------------------------------------------------------------

@pytest.mark.parametrize("name, expected", [
    ("P^(1,4)_(1,0,1)", ([1, 0, 1], [1, 4])),
    ("P^(0)_(0,8)", ([0, 8], [0])),
    ("P^(1,1,1)_(1,0,0,1)", ([1, 0, 0, 1], [1, 1, 1])),
    ("C(2,4,4)", None),
    ("3033030", None),
    ("P^(1,2)_(1,1)", None),          # one k too few to be a quipu name
])
def test_the_page_reads_the_quipu_back_out_of_a_class_name(name, expected):
    """The drawing comes from the name, so the parse has to be exact.

    `quipuForms` is the authority on the notation; this is the page reading it
    back, and the two must agree on what the parameters are.
    """
    assert classpage._quipuParameters(name) == expected
    if expected is not None:
        k, m = expected
        assert qf.formatQuipu((tuple(k), tuple(m))) == name


@pytest.mark.parametrize("relations, spans", [
    ("", []),
    ("1;2;3", [(1, 3)]),
    ("1;2;3|3;4;5;6", [(1, 3), (3, 6)]),
    ("1;2;3;4|4;5;6;7|6;7;8;9", [(1, 4), (4, 7), (6, 9)]),
])
def test_the_page_draws_a_relation_over_the_vertices_it_spans(relations, spans):
    """The arc has to sit over the relation's own span, which is its whole point:
    two arcs overlapping is two relations sharing arrows."""
    assert classpage._relationSpans(relations) == spans


def test_the_page_carries_the_whole_classification_and_renders_it():
    classification = classificationFrom(SAMPLE)
    data = classpage.collect(classification)
    assert data["length"] == 6
    assert data["lnaCount"] == 5
    assert data["unsettled"] == 1
    assert len(data["classes"]) == 4
    assert all(c["membersComplete"] for c in data["classes"])

    page = classpage.render(classification)
    assert "<title>A_6 Mutation Classes</title>" in page
    # the data is embedded, not fetched: no network for a page to work
    assert "://" not in page.split("<script id=\"data\"")[1].split("</script>")[0]
    assert json.loads(page.split('type="application/json">')[1].split("</script>")[0])
    for mutationClass in data["classes"]:
        assert mutationClass["name"] in page or mutationClass["name"] in json.dumps(data)

    standalone = classpage.renderStandalone(classification)
    assert standalone.startswith("<!doctype html>")
    assert page in standalone


def test_the_page_says_so_rather_than_silently_truncating(monkeypatch):
    """A length past the member limit still gets complete class summaries."""
    monkeypatch.setattr(classpage, "MEMBER_LIMIT", 2)
    data = classpage.collect(classificationFrom(SAMPLE))
    assert sum(c["size"] for c in data["classes"]) == 5      # every class counted
    assert sum(len(c["members"]) for c in data["classes"]) == 2
    assert not all(c["membersComplete"] for c in data["classes"])
    assert "at most" in classpage.SCRIPT


# -- against a real classification ----------------------------------------

@pytest.mark.slow
def test_against_a_real_classification():
    """The n = 6 classification, computed rather than written by hand.

    Four classes, 42 LNAs, every one a quipu class, and no two sharing a Coxeter
    polynomial -- which is the published answer for n = 6.
    """
    from helpers import quiet
    import quivermutation as qm

    table, _report = quiet(qm.classifyLength, 6, printOutput = False)
    classification = cv.Classification.fromTable(table)
    found = classification.classes()
    assert sum(c.size for c in found) == 42
    assert sorted(c.size for c in found) == [1, 12, 13, 16]
    assert {c.kind for c in found} == {cv.QUIPU}
    assert classification.coxeterCollisions() == {}
    assert classification.unsettled() == []

    for mutationClass in found:
        assert nk.LinearNakayamaAlgebra.fromRelationString(
            6, mutationClass.representative).quipuName() == mutationClass.name
        assert classpage._quipuParameters(mutationClass.name) is not None
