"""The table of mutation classes."""

import csv

import pytest

from quivermutation import mutationClassTable as mct


def sample_table():
    return mct.MutationClassTable(
        [
            ["", "000", "", "p_A", "1;2;3;4;5"],
            ["1;2;3", "000", "4;3;2", "p_A", "1;5;2;3;4"],
            ["1;2;3;4", "230", "1;1;1", "p_D", "3;2;5;1;4"],
            ["2;3;4;5", "", "", "", ""],
        ],
        lineLength=5,
    )


def test_lookup_by_relation_string():
    table = sample_table()
    assert len(table) == 4
    assert table.rowFor("1;2;3;4")[1] == "230"
    assert table.rowFor("")[1] == "000"           # A_5 itself, no relations
    assert table.rowFor("9;9;9") is None


def test_unassigned_rows_and_class_membership():
    table = sample_table()
    assert table.unassignedRelationStrings() == ["2;3;4;5"]
    assert table.classNames() == {"000", "230"}
    assert table.membersOfClass("000") == ["", "1;2;3"]


def test_assign_fills_a_row():
    table = sample_table()
    table.assign("2;3;4;5", "230", "1", "p_D", "2;1;3;4;5", "P^(2)_(1,1)")
    assert table.rowFor("2;3;4;5") == [
        "2;3;4;5", "230", "1", "p_D", "2;1;3;4;5", "P^(2)_(1,1)",
    ]
    assert table.unassignedRelationStrings() == []


def test_rename_class_moves_every_member():
    table = sample_table()
    table.renameClass("230", "300")
    assert table.classNames() == {"000", "300"}
    assert table.membersOfClass("300") == ["1;2;3;4"]


def test_classes_by_coxeter_polynomial_flags_merge_candidates():
    """A polynomial carrying two class names is a merge candidate."""
    table = mct.MutationClassTable([
        ["a", "X", "", "same_poly", ""],
        ["b", "Y", "", "same_poly", ""],
        ["c", "Z", "", "other_poly", ""],
    ])
    grouped = table.classesByCoxeterPolynomial()
    assert grouped == {"same_poly": {"X", "Y"}, "other_poly": {"Z"}}
    candidates = {p: names for p, names in grouped.items() if len(names) > 1}
    assert candidates == {"same_poly": {"X", "Y"}}


def test_csv_round_trip_with_and_without_a_header(tmp_path):
    table = sample_table()

    with_header = tmp_path / "with.csv"
    table.writeCSV(with_header, header=True)
    assert mct.MutationClassTable.fromCSV(with_header).rows() == table.rows()
    with open(with_header, newline="") as f:
        assert list(csv.reader(f))[0] == mct.COLUMNS

    without_header = tmp_path / "without.csv"
    table.writeCSV(without_header)
    assert mct.MutationClassTable.fromCSV(without_header).rows() == table.rows()


def test_parquet_round_trip(tmp_path):
    table = sample_table()
    path = tmp_path / "t.parquet"
    table.writeParquet(path)
    assert mct.MutationClassTable.fromParquet(path).rows() == table.rows()


def test_a_five_column_table_reads_back_with_an_empty_hereditary_column():
    """Tables written before the hereditary column existed must still load."""
    table = mct.MutationClassTable([["1;2;3", "000", "4;3", "p_A", "1;5;2;3;4"]])
    assert table.rowFor("1;2;3") == ["1;2;3", "000", "4;3", "p_A", "1;5;2;3;4", ""]


def test_hereditary_form_groups_classes_that_are_provably_the_same():
    table = mct.MutationClassTable([
        ["a", "X", "", "same_poly", "", "P^(3)_(1,1)"],
        ["b", "Y", "", "same_poly", "", "P^(3)_(1,1)"],
        ["c", "Z", "", "same_poly", "", "P^(2)_(1,2)"],
        ["d", "W", "", "other_poly", "", ""],
    ])
    assert table.classesByHereditaryForm() == {
        "P^(3)_(1,1)": {"X", "Y"},
        "P^(2)_(1,2)": {"Z"},
    }


def test_set_hereditary_form_for_class_touches_every_member():
    table = sample_table()
    table.setHereditaryFormForClass("000", "P^(0)_(0,4)")
    assert [row[5] for row in table.rows()] == ["P^(0)_(0,4)", "P^(0)_(0,4)", "", ""]


def test_dataframe_columns_are_all_strings():
    frame = sample_table().toDataFrame()
    assert frame.columns == mct.COLUMNS
    assert frame.height == 4
    assert all(str(dtype) == "String" for dtype in frame.dtypes)


def test_for_length_starts_empty():
    table = mct.MutationClassTable.forLength(4, ["", "1;2;3", "2;3;4"])
    assert len(table) == 3
    assert table.unassignedRelationStrings() == ["", "1;2;3", "2;3;4"]
    assert table.classNames() == set()
