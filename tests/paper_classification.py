"""The complete derived-equivalence classification of linear Nakayama algebras
of length n <= 8 with almost separate relations, transcribed from the table in

    D. Fosse, "Quipu quivers and Nakayama algebras with almost separate
    relations", arXiv:2305.06642, section "Complete classification for n <= 8".

Each entry maps a quipu to the list of Nakayama algebras derived equivalent to
its path algebra.  An algebra is written the way the paper writes it: the pair
(starts, lengths) of A_{n,(n_0,...,n_r)}^{(l_0,...,l_r)}, where the relation
with index i runs over l_i arrows starting at vertex n_i.  Only relations of
length >= 3 appear, since the paper shows a relation of length 2 never changes
the class of such an algebra.
"""

# quipu label -> list of (relation starts, relation lengths)
PAPER_CLASSES = {
    1: {"A_1": [((), ())]},
    2: {"A_2": [((), ())]},
    3: {"A_3": [((), ())]},
    4: {
        "A_4": [((), ())],
        "D_4": [((1,), (3,))],
    },
    5: {
        "A_5": [((), ())],
        "D_5": [((1,), (3,)), ((1,), (4,)), ((2,), (3,))],
    },
    6: {
        "A_6": [((), ())],
        "D_6": [((1,), (3,)), ((1,), (5,)), ((3,), (3,))],
        "E_6": [((1,), (4,)), ((2,), (3,)), ((2,), (4,))],
        "D~_5": [((1, 3), (3, 3))],
    },
    7: {
        "A_7": [((), ())],
        "D_7": [((1,), (3,)), ((1,), (6,)), ((4,), (3,))],
        "E_7": [
            ((1,), (4,)), ((1,), (5,)), ((2,), (3,)),
            ((2,), (5,)), ((3,), (3,)), ((3,), (4,)),
        ],
        "D~_6": [((1, 4), (3, 3))],
        "E~_6": [((2,), (4,))],
        "P_(1,0,2)^(1,1)": [
            ((1, 3), (3, 3)), ((1, 3), (3, 4)),
            ((1, 4), (4, 3)), ((2, 4), (3, 3)),
        ],
    },
    8: {
        "A_8": [((), ())],
        "D_8": [((1,), (3,)), ((1,), (7,)), ((5,), (3,))],
        "E_8": [
            ((1,), (4,)), ((1,), (6,)), ((2,), (3,)),
            ((2,), (6,)), ((4,), (3,)), ((4,), (4,)),
        ],
        "D~_7": [((1, 5), (3, 3))],
        "E~_7": [((1,), (5,)), ((3,), (3,)), ((3,), (5,))],
        "P_(1,1,2)^(1,1)": [
            ((1, 4), (3, 3)), ((1, 4), (3, 4)),
            ((1, 5), (4, 3)), ((2, 5), (3, 3)),
        ],
        "P_(1,0,3)^(1,1)": [
            ((1, 3), (3, 3)), ((1, 3), (3, 5)),
            ((1, 5), (5, 3)), ((3, 5), (3, 3)),
        ],
        "P_(1,0,0,1)^(1,1,1)": [((1, 3, 5), (3, 3, 3))],
        "P_(2,0,2)^(1,1)": [
            ((1, 4), (4, 3)), ((1, 4), (4, 4)),
            ((2, 4), (3, 3)), ((2, 4), (3, 4)),
        ],
        "P_(1,0,2)^(1,2)": [((1, 3), (3, 4)), ((2, 5), (4, 3))],
        "P_(2,3)^(2)": [((2,), (4,)), ((2,), (5,)), ((3,), (4,))],
    },
}


def rel_lengths(length, starts, lengths):
    """Paper notation -> the per-vertex relation-length list this repo uses.

    Entry i (0-based) of the result is the number of arrows in the relation
    starting at vertex i+1, or 0 if no relation starts there.  The list has
    length - 2 entries, since no relation can start at the last two vertices.
    """
    out = [0] * (length - 2)
    for start, rel_length in zip(starts, lengths):
        out[start - 1] = rel_length
    return out


def relation_string(length, starts, lengths):
    """Paper notation -> the 'a;b;c|d;e;f' relation string used in the CSVs."""
    paths = [
        ";".join(str(v) for v in range(start, start + rel_length + 1))
        for start, rel_length in zip(starts, lengths)
    ]
    return "|".join(sorted(paths))
