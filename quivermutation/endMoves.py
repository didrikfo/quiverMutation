"""The verified rewrites that hold against one end of the quiver.

A rule here is a rewrite of a window of arrows *anchored* to an end: it is
stated as holding with its window flush against the source or the sink of the
line, and `lnaMoves.windowStartsFor` offers it no other position.  That is not a
technicality.  The source has no arrow into it and the sink none out of it, so a
mutation there does something a mutation in the interior cannot, and the rewrite
built on it is true at the end and false everywhere else -- which is why these
could not be found by, or stated in, the search for rules that hold at every
position.

They matter out of all proportion to their two positions.  An isolated pair of
relations sharing two or more arrows is the commonest thing a classification
search still has to place (F-021) and cannot be pulled apart anywhere in the
interior (F-022); against an end it comes apart.  Adding these takes the share
of LNAs placed with no search at all from 51% to 60% at n = 9, and to 100% at
n = 6.

Where they came from, and what is not here.  `discover.py --anchor both
--max-arrows 5 --max-width 6` plants each pattern flush against each end of
A_11 and A_12, mutates within three vertices of it, and keeps the rewrites that
recur at both lengths; `verifyMove` then checks each at the four lengths its
window fits in, at the one position it claims.  That run verified **630** rules
that are not a floating rule restricted to an end.  Listed below are the **229**
of them that change the orbit partition at n <= 9 -- the rest reach nothing the
229 do not.  Re-run the command to get them all back.  Note what that curation
does and does not justify: a rule dropped here is redundant *at the lengths
measured*, and could in principle be the one that matters at n >= 10.

Most of them carry a **spectator** -- a relation inside the window that the
rewrite leaves exactly where it is.  188 of the 229 do, and that is the point of
them: the rules stated on a clean window are true but rarely match, because a
real LNA usually has something else nearby (research H-011).  They do not
compress into a handful of families the way the slide rules do: 190 distinct
rewrites once the spectators are set aside, so they are listed rather than
generated.

Each entry is (window width in arrows, relations before, relations after,
mutation sequence, which end), with relation starts and mutation vertices given
relative to the window's first arrow and a negative vertex meaning a left
mutation.  The count in each comment is how many applications `verifyMove`
confirmed, a confirmation being an admissible sequence landing on the predicted
LNA with the Coxeter polynomial kept.  E-023.
"""

DISCOVERED_END_MOVES = [
    (2, ((0, 2),), (), (1,), 'left'),   # 9 confirmed: window 2 arrows at the left end: (0:2)  ->  -   via [1]

    (2, ((0, 2),), (), (-3,), 'right'),   # 9 confirmed: window 2 arrows at the right end: (0:2)  ->  -   via [-3]

    (4, ((0, 2), (1, 2)), ((0, 2), (2, 2)), (-4,), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:2) (1:2)  ->  (0:2) (2:2)   via [-4]
    (4, ((0, 2), (2, 2)), ((0, 2), (1, 2)), (3,), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:2) (2:2)  ->  (0:2) (1:2)   via [3]
    (4, ((0, 2), (2, 2)), ((1, 2),), (1, 3), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:2) (2:2)  ->  (1:2)   via [1, 3]
    (4, ((0, 3), (1, 3)), ((0, 3), (2, 2)), (-4, -2), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:3) (1:3)  ->  (0:3) (2:2)   via [-4, -2]
    (4, ((0, 3), (2, 2)), ((0, 3), (1, 3)), (1, 1), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:3) (2:2)  ->  (0:3) (1:3)   via [1, 1]
    (4, ((0, 4),), ((1, 3),), (1, -3, -4), 'left'),   # 9 confirmed: window 4 arrows at the left end: (0:4)  ->  (1:3)   via [1, -3, -4]

    (4, ((0, 2), (1, 3)), ((0, 3), (1, 3)), (-5, -3), 'right'),   # 9 confirmed: window 4 arrows at the right end: (0:2) (1:3)  ->  (0:3) (1:3)   via [-5, -3]
    (4, ((0, 2), (2, 2)), ((1, 2),), (-3, -5), 'right'),   # 9 confirmed: window 4 arrows at the right end: (0:2) (2:2)  ->  (1:2)   via [-3, -5]
    (4, ((0, 2), (2, 2)), ((1, 2), (2, 2)), (-3,), 'right'),   # 9 confirmed: window 4 arrows at the right end: (0:2) (2:2)  ->  (1:2) (2:2)   via [-3]
    (4, ((0, 3), (1, 3)), ((0, 2), (1, 3)), (2, 2), 'right'),   # 9 confirmed: window 4 arrows at the right end: (0:3) (1:3)  ->  (0:2) (1:3)   via [2, 2]
    (4, ((0, 4),), ((0, 3),), (3, 2, -5), 'right'),   # 9 confirmed: window 4 arrows at the right end: (0:4)  ->  (0:3)   via [3, 2, -5]
    (4, ((1, 2), (2, 2)), ((0, 2), (2, 2)), (2,), 'right'),   # 9 confirmed: window 4 arrows at the right end: (1:2) (2:2)  ->  (0:2) (2:2)   via [2]

    (5, ((0, 2), (1, 2), (2, 2)), ((0, 2), (1, 2), (3, 2)), (-5,), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:2) (1:2) (2:2)  ->  (0:2) (1:2) (3:2)   via [-5]
    (5, ((0, 2), (1, 2), (2, 2)), ((1, 2), (3, 2)), (1, -5), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:2) (1:2) (2:2)  ->  (1:2) (3:2)   via [1, -5]
    (5, ((0, 2), (1, 2), (3, 2)), ((0, 2), (1, 2), (2, 2)), (4,), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:2) (1:2) (3:2)  ->  (0:2) (1:2) (2:2)   via [4]
    (5, ((0, 3), (1, 3)), ((0, 3), (3, 2)), (-4, -2, -5), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3)  ->  (0:3) (3:2)   via [-4, -2, -5]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 3), (1, 4)), (1, 1), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3) (2:3)  ->  (0:3) (1:4)   via [1, 1]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 4), (2, 3)), (1, 3, -5), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3) (2:3)  ->  (0:4) (2:3)   via [1, 3, -5]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 4), (3, 2)), (1, -5, -3), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3) (2:3)  ->  (0:4) (3:2)   via [1, -5, -3]
    (5, ((0, 3), (1, 3), (2, 3)), ((1, 3),), (2, 1, 2), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3) (2:3)  ->  (1:3)   via [2, 1, 2]
    (5, ((0, 3), (1, 3), (3, 2)), ((1, 3), (2, 3)), (2, 1, 2), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:3) (3:2)  ->  (1:3) (2:3)   via [2, 1, 2]
    (5, ((0, 3), (1, 4)), ((0, 3), (1, 3), (2, 3)), (-4, -4), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:4)  ->  (0:3) (1:3) (2:3)   via [-4, -4]
    (5, ((0, 3), (1, 4)), ((0, 4), (3, 2)), (-4, -5, -2), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:4)  ->  (0:4) (3:2)   via [-4, -5, -2]
    (5, ((0, 3), (1, 4)), ((1, 3),), (1, -4, 1), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (1:4)  ->  (1:3)   via [1, -4, 1]
    (5, ((0, 3), (2, 2)), ((0, 3), (3, 2)), (-5,), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (2:2)  ->  (0:3) (3:2)   via [-5]
    (5, ((0, 3), (3, 2)), ((0, 3), (1, 3)), (1, 4, 1), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (3:2)  ->  (0:3) (1:3)   via [1, 4, 1]
    (5, ((0, 3), (3, 2)), ((0, 3), (2, 2)), (4,), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:3) (3:2)  ->  (0:3) (2:2)   via [4]
    (5, ((0, 4), (1, 4)), ((0, 4), (2, 3)), (-5, -5), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (1:4)  ->  (0:4) (2:3)   via [-5, -5]
    (5, ((0, 4), (2, 3)), ((0, 3), (1, 3), (2, 3)), (3, -5, -2), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (2:3)  ->  (0:3) (1:3) (2:3)   via [3, -5, -2]
    (5, ((0, 4), (2, 3)), ((0, 4), (1, 4)), (1, 1), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (2:3)  ->  (0:4) (1:4)   via [1, 1]
    (5, ((0, 4), (2, 3)), ((0, 4), (3, 2)), (-5, -5), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (2:3)  ->  (0:4) (3:2)   via [-5, -5]
    (5, ((0, 4), (2, 3)), ((1, 3), (3, 2)), (-5, -2, -3), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (2:3)  ->  (1:3) (3:2)   via [-5, -2, -3]
    (5, ((0, 4), (3, 2)), ((0, 3), (1, 3), (2, 3)), (1, -3, 4), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (3:2)  ->  (0:3) (1:3) (2:3)   via [1, -3, 4]
    (5, ((0, 4), (3, 2)), ((0, 3), (1, 4)), (1, 4, 3), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (3:2)  ->  (0:3) (1:4)   via [1, 4, 3]
    (5, ((0, 4), (3, 2)), ((0, 4), (2, 3)), (1, 1), 'left'),   # 9 confirmed: window 5 arrows at the left end: (0:4) (3:2)  ->  (0:4) (2:3)   via [1, 1]

    (5, ((0, 2), (1, 3), (2, 3)), ((0, 3), (1, 3)), (-5, -3, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (1:3) (2:3)  ->  (0:3) (1:3)   via [-5, -3, -6]
    (5, ((0, 2), (1, 4)), ((0, 3), (1, 3), (2, 3)), (4, -6, -3), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (1:4)  ->  (0:3) (1:3) (2:3)   via [4, -6, -3]
    (5, ((0, 2), (1, 4)), ((0, 3), (1, 4)), (-6, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (1:4)  ->  (0:3) (1:4)   via [-6, -6]
    (5, ((0, 2), (1, 4)), ((0, 4), (2, 3)), (-6, -3, -4), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (1:4)  ->  (0:4) (2:3)   via [-6, -3, -4]
    (5, ((0, 2), (2, 2), (3, 2)), ((1, 2), (2, 2), (3, 2)), (-3,), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (2:2) (3:2)  ->  (1:2) (2:2) (3:2)   via [-3]
    (5, ((0, 2), (2, 3)), ((0, 3), (1, 3)), (-6, 3, -4), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (2:3)  ->  (0:3) (1:3)   via [-6, 3, -4]
    (5, ((0, 2), (2, 3)), ((1, 2), (2, 3)), (-3,), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (2:3)  ->  (1:2) (2:3)   via [-3]
    (5, ((0, 2), (2, 3)), ((1, 3), (2, 3)), (-3, -6, -4), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:2) (2:3)  ->  (1:3) (2:3)   via [-3, -6, -4]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 2), (1, 4)), (2, 4, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:3) (2:3)  ->  (0:2) (1:4)   via [2, 4, -6]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 3), (1, 4)), (2, -6, -4), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:3) (2:3)  ->  (0:3) (1:4)   via [2, -6, -4]
    (5, ((0, 3), (1, 3), (2, 3)), ((0, 4), (2, 3)), (-6, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:3) (2:3)  ->  (0:4) (2:3)   via [-6, -6]
    (5, ((0, 3), (1, 3), (2, 3)), ((1, 3),), (-5, -5, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:3) (2:3)  ->  (1:3)   via [-5, -5, -6]
    (5, ((0, 3), (1, 4)), ((0, 2), (1, 3)), (2, 5, 4), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:4)  ->  (0:2) (1:3)   via [2, 5, 4]
    (5, ((0, 3), (1, 4)), ((0, 2), (1, 4)), (2, 2), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:4)  ->  (0:2) (1:4)   via [2, 2]
    (5, ((0, 3), (1, 4)), ((0, 3), (1, 3), (2, 3)), (2, -4, 5), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:4)  ->  (0:3) (1:3) (2:3)   via [2, -4, 5]
    (5, ((0, 3), (1, 4)), ((0, 4), (1, 4)), (-6, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:3) (1:4)  ->  (0:4) (1:4)   via [-6, -6]
    (5, ((0, 4), (1, 4)), ((0, 3), (1, 4)), (2, 2), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:4) (1:4)  ->  (0:3) (1:4)   via [2, 2]
    (5, ((0, 4), (2, 3)), ((0, 2), (1, 4)), (3, 2, 5), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:4) (2:3)  ->  (0:2) (1:4)   via [3, 2, 5]
    (5, ((0, 4), (2, 3)), ((0, 3), (1, 3), (2, 3)), (3, 3), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:4) (2:3)  ->  (0:3) (1:3) (2:3)   via [3, 3]
    (5, ((0, 4), (2, 3)), ((1, 3),), (-6, 3, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (0:4) (2:3)  ->  (1:3)   via [-6, 3, -6]
    (5, ((1, 2), (2, 2), (3, 2)), ((0, 2), (2, 2)), (2, -6), 'right'),   # 9 confirmed: window 5 arrows at the right end: (1:2) (2:2) (3:2)  ->  (0:2) (2:2)   via [2, -6]
    (5, ((1, 2), (2, 2), (3, 2)), ((0, 2), (2, 2), (3, 2)), (2,), 'right'),   # 9 confirmed: window 5 arrows at the right end: (1:2) (2:2) (3:2)  ->  (0:2) (2:2) (3:2)   via [2]
    (5, ((1, 2), (2, 3)), ((0, 2), (2, 3)), (2,), 'right'),   # 9 confirmed: window 5 arrows at the right end: (1:2) (2:3)  ->  (0:2) (2:3)   via [2]

    (6, ((0, 2), (1, 2), (2, 2)), ((0, 2), (1, 2), (4, 2)), (-5, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:2) (2:2)  ->  (0:2) (1:2) (4:2)   via [-5, -6]
    (6, ((0, 2), (1, 2), (2, 3)), ((0, 2), (1, 3), (2, 3), (3, 3)), (-6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:2) (2:3)  ->  (0:2) (1:3) (2:3) (3:3)   via [-6, -6]
    (6, ((0, 2), (1, 2), (2, 3)), ((1, 3), (2, 3), (3, 3)), (1, -6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:2) (2:3)  ->  (1:3) (2:3) (3:3)   via [1, -6, -6]
    (6, ((0, 2), (1, 2), (4, 2)), ((0, 2), (1, 2), (2, 2)), (5, 4), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:2) (4:2)  ->  (0:2) (1:2) (2:2)   via [5, 4]
    (6, ((0, 2), (1, 3), (2, 3)), ((0, 2), (2, 3), (3, 3)), (-6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (2:3)  ->  (0:2) (2:3) (3:3)   via [-6, -6]
    (6, ((0, 2), (1, 3), (2, 4)), ((0, 2), (1, 4), (3, 3)), (3, -5, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (2:4)  ->  (0:2) (1:4) (3:3)   via [3, -5, -6]
    (6, ((0, 2), (1, 3), (2, 4)), ((0, 3), (1, 4), (4, 2)), (-5, -6, -3), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (2:4)  ->  (0:3) (1:4) (4:2)   via [-5, -6, -3]
    (6, ((0, 2), (1, 3), (3, 2)), ((0, 2), (1, 3), (4, 2)), (-6,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (3:2)  ->  (0:2) (1:3) (4:2)   via [-6]
    (6, ((0, 2), (1, 3), (3, 2)), ((1, 3), (4, 2)), (1, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (3:2)  ->  (1:3) (4:2)   via [1, -6]
    (6, ((0, 2), (1, 3), (4, 2)), ((0, 2), (1, 3), (3, 2)), (5,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (4:2)  ->  (0:2) (1:3) (3:2)   via [5]
    (6, ((0, 2), (1, 3), (4, 2)), ((1, 3), (3, 2)), (1, 5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:3) (4:2)  ->  (1:3) (3:2)   via [1, 5]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 2), (1, 3), (2, 4)), (4, 3, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:4) (3:3)  ->  (0:2) (1:3) (2:4)   via [4, 3, -6]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 3), (2, 3), (3, 3)), (4, -6, -3), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (1:4) (3:3)  ->  (0:3) (1:3) (2:3) (3:3)   via [4, -6, -3]
    (6, ((0, 2), (2, 2), (3, 2)), ((0, 2), (2, 2), (4, 2)), (-6,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (2:2) (3:2)  ->  (0:2) (2:2) (4:2)   via [-6]
    (6, ((0, 2), (2, 3), (4, 2)), ((0, 2), (1, 3), (2, 3), (3, 3)), (3, 3), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (2:3) (4:2)  ->  (0:2) (1:3) (2:3) (3:3)   via [3, 3]
    (6, ((0, 2), (2, 3), (4, 2)), ((1, 3), (2, 3), (3, 3)), (1, 3, 3), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:2) (2:3) (4:2)  ->  (1:3) (2:3) (3:3)   via [1, 3, 3]
    (6, ((0, 3), (1, 3), (2, 4)), ((0, 3), (1, 5)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:3) (2:4)  ->  (0:3) (1:5)   via [1, 1]
    (6, ((0, 3), (1, 3), (3, 2)), ((0, 3), (2, 2), (4, 2)), (-4, -2, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:3) (3:2)  ->  (0:3) (2:2) (4:2)   via [-4, -2, -6]
    (6, ((0, 3), (1, 3), (4, 2)), ((0, 3), (2, 2), (3, 2)), (-4, -2, 5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:3) (4:2)  ->  (0:3) (2:2) (3:2)   via [-4, -2, 5]
    (6, ((0, 3), (1, 4), (2, 4)), ((1, 4),), (2, 1, 2), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (2:4)  ->  (1:4)   via [2, 1, 2]
    (6, ((0, 3), (1, 4), (3, 3)), ((0, 4), (2, 3)), (1, 1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (3:3)  ->  (0:4) (2:3)   via [1, 1, 1]
    (6, ((0, 3), (1, 4), (3, 3)), ((1, 3), (2, 4)), (1, -4, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (3:3)  ->  (1:3) (2:4)   via [1, -4, 1]
    (6, ((0, 3), (1, 4), (3, 3)), ((1, 4), (2, 4)), (2, 1, 2), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (3:3)  ->  (1:4) (2:4)   via [2, 1, 2]
    (6, ((0, 3), (1, 4), (4, 2)), ((0, 4), (2, 3), (3, 3)), (1, 1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (4:2)  ->  (0:4) (2:3) (3:3)   via [1, 1, 1]
    (6, ((0, 3), (1, 4), (4, 2)), ((1, 4), (3, 3)), (2, 1, 2), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:4) (4:2)  ->  (1:4) (3:3)   via [2, 1, 2]
    (6, ((0, 3), (1, 5)), ((0, 3), (1, 3), (2, 4)), (-4, -4), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (1:5)  ->  (0:3) (1:3) (2:4)   via [-4, -4]
    (6, ((0, 3), (2, 2)), ((0, 3), (4, 2)), (-5, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (2:2)  ->  (0:3) (4:2)   via [-5, -6]
    (6, ((0, 3), (2, 2), (3, 2)), ((0, 3), (1, 3), (4, 2)), (1, 1, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (2:2) (3:2)  ->  (0:3) (1:3) (4:2)   via [1, 1, -6]
    (6, ((0, 3), (2, 2), (3, 2)), ((0, 3), (2, 2), (4, 2)), (-6,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (2:2) (3:2)  ->  (0:3) (2:2) (4:2)   via [-6]
    (6, ((0, 3), (2, 2), (4, 2)), ((0, 3), (1, 3), (3, 2)), (1, 1, 5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (2:2) (4:2)  ->  (0:3) (1:3) (3:2)   via [1, 1, 5]
    (6, ((0, 3), (2, 2), (4, 2)), ((0, 3), (2, 2), (3, 2)), (5,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (2:2) (4:2)  ->  (0:3) (2:2) (3:2)   via [5]
    (6, ((0, 3), (4, 2)), ((0, 3), (2, 2)), (5, 4), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:3) (4:2)  ->  (0:3) (2:2)   via [5, 4]
    (6, ((0, 4), (1, 4), (2, 4)), ((0, 4), (1, 5)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (1:4) (2:4)  ->  (0:4) (1:5)   via [1, 1]
    (6, ((0, 4), (1, 4), (2, 4)), ((0, 4), (2, 3), (3, 3)), (-5, -5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (1:4) (2:4)  ->  (0:4) (2:3) (3:3)   via [-5, -5]
    (6, ((0, 4), (1, 4), (3, 3)), ((0, 4), (2, 4)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (1:4) (3:3)  ->  (0:4) (2:4)   via [1, 1]
    (6, ((0, 4), (1, 5)), ((0, 4), (1, 4), (2, 4)), (-5, -5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (1:5)  ->  (0:4) (1:4) (2:4)   via [-5, -5]
    (6, ((0, 4), (2, 3)), ((0, 3), (1, 4), (3, 3)), (-6, -6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (2:3)  ->  (0:3) (1:4) (3:3)   via [-6, -6, -6]
    (6, ((0, 4), (2, 3)), ((0, 4), (4, 2)), (-5, -5, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (2:3)  ->  (0:4) (4:2)   via [-5, -5, -6]
    (6, ((0, 4), (2, 3), (3, 3)), ((0, 3), (1, 4), (4, 2)), (-6, -6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (2:3) (3:3)  ->  (0:3) (1:4) (4:2)   via [-6, -6, -6]
    (6, ((0, 4), (2, 4)), ((0, 4), (1, 4), (3, 3)), (-5, -5), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (2:4)  ->  (0:4) (1:4) (3:3)   via [-5, -5]
    (6, ((0, 4), (3, 2)), ((0, 4), (4, 2)), (-6,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (3:2)  ->  (0:4) (4:2)   via [-6]
    (6, ((0, 4), (4, 2)), ((0, 4), (3, 2)), (5,), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:4) (4:2)  ->  (0:4) (3:2)   via [5]
    (6, ((0, 5), (1, 5)), ((0, 5), (2, 4)), (-6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (1:5)  ->  (0:5) (2:4)   via [-6, -6]
    (6, ((0, 5), (2, 4)), ((0, 5), (1, 5)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (2:4)  ->  (0:5) (1:5)   via [1, 1]
    (6, ((0, 5), (2, 4)), ((0, 5), (3, 3)), (-6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (2:4)  ->  (0:5) (3:3)   via [-6, -6]
    (6, ((0, 5), (3, 3)), ((0, 5), (2, 4)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (3:3)  ->  (0:5) (2:4)   via [1, 1]
    (6, ((0, 5), (3, 3)), ((0, 5), (4, 2)), (-6, -6), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (3:3)  ->  (0:5) (4:2)   via [-6, -6]
    (6, ((0, 5), (4, 2)), ((0, 5), (3, 3)), (1, 1), 'left'),   # 9 confirmed: window 6 arrows at the left end: (0:5) (4:2)  ->  (0:5) (3:3)   via [1, 1]

    (6, ((0, 2), (1, 3), (2, 4)), ((0, 3), (1, 4), (4, 2)), (-5, -6, -3), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:3) (2:4)  ->  (0:3) (1:4) (4:2)   via [-5, -6, -3]
    (6, ((0, 2), (1, 3), (4, 2)), ((0, 3), (1, 3), (2, 3)), (-5, -5, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:3) (4:2)  ->  (0:3) (1:3) (2:3)   via [-5, -5, -7]
    (6, ((0, 2), (1, 3), (4, 2)), ((0, 3), (1, 3), (2, 3), (4, 2)), (-5, -5), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:3) (4:2)  ->  (0:3) (1:3) (2:3) (4:2)   via [-5, -5]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 3), (2, 3), (3, 3)), (4, -6, -3), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:4) (3:3)  ->  (0:3) (1:3) (2:3) (3:3)   via [4, -6, -3]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 3), (2, 4)), (-7, -7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:4) (3:3)  ->  (0:3) (1:3) (2:4)   via [-7, -7, -7]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 4)), (-6, -6, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:4) (3:3)  ->  (0:3) (1:4)   via [-6, -6, -7]
    (6, ((0, 2), (1, 4), (3, 3)), ((0, 4), (2, 3), (4, 2)), (-6, -3, -4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:4) (3:3)  ->  (0:4) (2:3) (4:2)   via [-6, -3, -4]
    (6, ((0, 2), (1, 5)), ((0, 3), (1, 5)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (1:5)  ->  (0:3) (1:5)   via [-7, -7]
    (6, ((0, 2), (2, 2), (3, 3)), ((1, 2), (2, 2), (3, 3)), (-3,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:2) (3:3)  ->  (1:2) (2:2) (3:3)   via [-3]
    (6, ((0, 2), (2, 2), (3, 3)), ((1, 2), (2, 3), (3, 3)), (-3, -7, -5), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:2) (3:3)  ->  (1:2) (2:3) (3:3)   via [-3, -7, -5]
    (6, ((0, 2), (2, 3), (3, 3)), ((1, 2), (2, 2), (3, 3)), (-3, 4, 4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:3) (3:3)  ->  (1:2) (2:2) (3:3)   via [-3, 4, 4]
    (6, ((0, 2), (2, 3), (4, 2)), ((1, 2), (2, 3)), (-3, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:3) (4:2)  ->  (1:2) (2:3)   via [-3, -7]
    (6, ((0, 2), (2, 3), (4, 2)), ((1, 2), (2, 3), (4, 2)), (-3,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:3) (4:2)  ->  (1:2) (2:3) (4:2)   via [-3]
    (6, ((0, 2), (2, 4)), ((1, 2), (2, 4)), (-3,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:4)  ->  (1:2) (2:4)   via [-3]
    (6, ((0, 2), (2, 4)), ((1, 3), (2, 4)), (-3, -7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (2:4)  ->  (1:3) (2:4)   via [-3, -7, -7]
    (6, ((0, 2), (3, 2), (4, 2)), ((2, 2), (3, 2), (4, 2)), (-3, -4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (3:2) (4:2)  ->  (2:2) (3:2) (4:2)   via [-3, -4]
    (6, ((0, 2), (3, 3)), ((2, 2), (3, 3)), (-3, -4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:2) (3:3)  ->  (2:2) (3:3)   via [-3, -4]
    (6, ((0, 3), (1, 3), (2, 4)), ((0, 2), (1, 4), (3, 3)), (2, 2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:3) (2:4)  ->  (0:2) (1:4) (3:3)   via [2, 2, 2]
    (6, ((0, 3), (1, 4), (2, 4)), ((0, 4), (2, 4)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:4) (2:4)  ->  (0:4) (2:4)   via [-7, -7]
    (6, ((0, 3), (1, 4), (3, 3)), ((0, 4), (1, 4)), (-6, -6, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:4) (3:3)  ->  (0:4) (1:4)   via [-6, -6, -7]
    (6, ((0, 3), (1, 4), (3, 3)), ((0, 4), (2, 3)), (-7, 4, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:4) (3:3)  ->  (0:4) (2:3)   via [-7, 4, -7]
    (6, ((0, 3), (1, 4), (3, 3)), ((1, 3), (2, 4)), (-7, -7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:4) (3:3)  ->  (1:3) (2:4)   via [-7, -7, -7]
    (6, ((0, 3), (1, 5)), ((0, 2), (1, 5)), (2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:5)  ->  (0:2) (1:5)   via [2, 2]
    (6, ((0, 3), (1, 5)), ((0, 4), (1, 5)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:3) (1:5)  ->  (0:4) (1:5)   via [-7, -7]
    (6, ((0, 4), (1, 4), (2, 4)), ((0, 5), (2, 4)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (1:4) (2:4)  ->  (0:5) (2:4)   via [-7, -7]
    (6, ((0, 4), (1, 4), (3, 3)), ((1, 4),), (-6, -6, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (1:4) (3:3)  ->  (1:4)   via [-6, -6, -7]
    (6, ((0, 4), (1, 5)), ((0, 3), (1, 5)), (2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (1:5)  ->  (0:3) (1:5)   via [2, 2]
    (6, ((0, 4), (1, 5)), ((0, 5), (1, 5)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (1:5)  ->  (0:5) (1:5)   via [-7, -7]
    (6, ((0, 4), (2, 3), (3, 3)), ((0, 5), (3, 3)), (-7, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (2:3) (3:3)  ->  (0:5) (3:3)   via [-7, -7]
    (6, ((0, 4), (2, 4)), ((0, 3), (1, 4), (2, 4)), (3, 3), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:4) (2:4)  ->  (0:3) (1:4) (2:4)   via [3, 3]
    (6, ((0, 5), (1, 5)), ((0, 4), (1, 5)), (2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:5) (1:5)  ->  (0:4) (1:5)   via [2, 2]
    (6, ((0, 5), (2, 4)), ((0, 4), (1, 4), (2, 4)), (3, 3), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:5) (2:4)  ->  (0:4) (1:4) (2:4)   via [3, 3]
    (6, ((0, 5), (3, 3)), ((0, 4), (2, 3), (3, 3)), (4, 4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (0:5) (3:3)  ->  (0:4) (2:3) (3:3)   via [4, 4]
    (6, ((1, 2), (2, 2), (3, 3)), ((0, 2), (2, 2), (3, 3)), (2,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:2) (3:3)  ->  (0:2) (2:2) (3:3)   via [2]
    (6, ((1, 2), (2, 2), (3, 3)), ((0, 2), (2, 3), (3, 3)), (2, -7, -5), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:2) (3:3)  ->  (0:2) (2:3) (3:3)   via [2, -7, -5]
    (6, ((1, 2), (2, 3), (3, 3)), ((0, 2), (2, 2), (3, 3)), (2, 4, 4), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:3) (3:3)  ->  (0:2) (2:2) (3:3)   via [2, 4, 4]
    (6, ((1, 2), (2, 3), (4, 2)), ((0, 2), (2, 3)), (2, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:3) (4:2)  ->  (0:2) (2:3)   via [2, -7]
    (6, ((1, 2), (2, 3), (4, 2)), ((0, 2), (2, 3), (4, 2)), (2,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:3) (4:2)  ->  (0:2) (2:3) (4:2)   via [2]
    (6, ((1, 2), (2, 4)), ((0, 2), (2, 4)), (2,), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:2) (2:4)  ->  (0:2) (2:4)   via [2]
    (6, ((1, 3), (2, 4)), ((0, 2), (2, 4)), (3, 2, 3), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:3) (2:4)  ->  (0:2) (2:4)   via [3, 2, 3]
    (6, ((1, 3), (2, 4)), ((0, 3), (1, 4), (3, 3)), (2, 2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:3) (2:4)  ->  (0:3) (1:4) (3:3)   via [2, 2, 2]
    (6, ((1, 3), (3, 2), (4, 2)), ((0, 3), (1, 3), (2, 3)), (2, 2, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:3) (3:2) (4:2)  ->  (0:3) (1:3) (2:3)   via [2, 2, -7]
    (6, ((1, 3), (3, 2), (4, 2)), ((0, 3), (1, 3), (2, 3), (4, 2)), (2, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (1:3) (3:2) (4:2)  ->  (0:3) (1:3) (2:3) (4:2)   via [2, 2]
    (6, ((2, 2), (3, 2), (4, 2)), ((0, 2), (3, 2)), (3, 2, -7), 'right'),   # 9 confirmed: window 6 arrows at the right end: (2:2) (3:2) (4:2)  ->  (0:2) (3:2)   via [3, 2, -7]
    (6, ((2, 2), (3, 2), (4, 2)), ((0, 2), (3, 2), (4, 2)), (3, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (2:2) (3:2) (4:2)  ->  (0:2) (3:2) (4:2)   via [3, 2]
    (6, ((2, 2), (3, 3)), ((0, 2), (3, 3)), (3, 2), 'right'),   # 9 confirmed: window 6 arrows at the right end: (2:2) (3:3)  ->  (0:2) (3:3)   via [3, 2]

    (7, ((0, 2), (1, 2), (2, 4)), ((0, 2), (1, 3), (2, 4), (3, 4)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:2) (2:4)  ->  (0:2) (1:3) (2:4) (3:4)   via [-7, -7]
    (7, ((0, 2), (1, 2), (2, 4)), ((1, 3), (2, 4), (3, 4)), (1, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:2) (2:4)  ->  (1:3) (2:4) (3:4)   via [1, -7, -7]
    (7, ((0, 2), (1, 2), (3, 3)), ((0, 2), (1, 3), (2, 3), (3, 4)), (-7, 4, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:2) (3:3)  ->  (0:2) (1:3) (2:3) (3:4)   via [-7, 4, -7]
    (7, ((0, 2), (1, 3), (2, 4)), ((0, 2), (1, 4), (2, 4), (3, 4)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:3) (2:4)  ->  (0:2) (1:4) (2:4) (3:4)   via [-7, -7]
    (7, ((0, 2), (1, 3), (3, 2)), ((0, 2), (1, 3), (5, 2)), (-6, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:3) (3:2)  ->  (0:2) (1:3) (5:2)   via [-6, -7]
    (7, ((0, 2), (1, 4), (3, 3)), ((0, 3), (1, 3), (2, 4), (4, 3)), (-7, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:4) (3:3)  ->  (0:3) (1:3) (2:4) (4:3)   via [-7, -7, -7]
    (7, ((0, 2), (1, 4), (4, 2)), ((0, 2), (1, 4), (5, 2)), (-7,), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:4) (4:2)  ->  (0:2) (1:4) (5:2)   via [-7]
    (7, ((0, 2), (1, 4), (4, 2)), ((1, 4), (5, 2)), (1, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:4) (4:2)  ->  (1:4) (5:2)   via [1, -7]
    (7, ((0, 2), (1, 5)), ((0, 3), (1, 5), (2, 5)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (1:5)  ->  (0:3) (1:5) (2:5)   via [-7, -7]
    (7, ((0, 2), (2, 2), (3, 3)), ((1, 2), (2, 3), (3, 3), (4, 3)), (-3, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (2:2) (3:3)  ->  (1:2) (2:3) (3:3) (4:3)   via [-3, -7, -7]
    (7, ((0, 2), (2, 3), (4, 2)), ((0, 2), (2, 3), (5, 2)), (-7,), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (2:3) (4:2)  ->  (0:2) (2:3) (5:2)   via [-7]
    (7, ((0, 2), (2, 3), (4, 2)), ((1, 2), (2, 3), (5, 2)), (-3, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (2:3) (4:2)  ->  (1:2) (2:3) (5:2)   via [-3, -7]
    (7, ((0, 2), (2, 3), (4, 2)), ((2, 3), (5, 2)), (1, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (2:3) (4:2)  ->  (2:3) (5:2)   via [1, -7]
    (7, ((0, 2), (2, 4)), ((1, 3), (2, 4), (3, 4)), (-3, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:2) (2:4)  ->  (1:3) (2:4) (3:4)   via [-3, -7, -7]
    (7, ((0, 3), (1, 3), (2, 4)), ((0, 4), (1, 4), (2, 4), (3, 4)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (1:3) (2:4)  ->  (0:4) (1:4) (2:4) (3:4)   via [-7, -7]
    (7, ((0, 3), (1, 4), (2, 4)), ((0, 4), (2, 4), (3, 4)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (1:4) (2:4)  ->  (0:4) (2:4) (3:4)   via [-7, -7]
    (7, ((0, 3), (1, 4), (3, 3)), ((0, 4), (2, 3), (3, 4)), (-7, 4, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (1:4) (3:3)  ->  (0:4) (2:3) (3:4)   via [-7, 4, -7]
    (7, ((0, 3), (1, 4), (3, 3)), ((1, 3), (2, 4), (4, 3)), (-7, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (1:4) (3:3)  ->  (1:3) (2:4) (4:3)   via [-7, -7, -7]
    (7, ((0, 3), (1, 5)), ((0, 4), (1, 5), (2, 5)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (1:5)  ->  (0:4) (1:5) (2:5)   via [-7, -7]
    (7, ((0, 3), (2, 2), (3, 3)), ((0, 3), (2, 3), (3, 3), (4, 3)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (2:2) (3:3)  ->  (0:3) (2:3) (3:3) (4:3)   via [-7, -7]
    (7, ((0, 3), (2, 3), (3, 3)), ((0, 3), (3, 3), (4, 3)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (2:3) (3:3)  ->  (0:3) (3:3) (4:3)   via [-7, -7]
    (7, ((0, 3), (2, 3), (4, 2)), ((0, 3), (2, 3), (5, 2)), (-7,), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:3) (2:3) (4:2)  ->  (0:3) (2:3) (5:2)   via [-7]
    (7, ((0, 4), (1, 4), (2, 4)), ((0, 5), (2, 4), (3, 4)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (1:4) (2:4)  ->  (0:5) (2:4) (3:4)   via [-7, -7]
    (7, ((0, 4), (1, 4), (3, 3)), ((1, 4), (5, 2)), (-6, -6, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (1:4) (3:3)  ->  (1:4) (5:2)   via [-6, -6, -7]
    (7, ((0, 4), (1, 5)), ((0, 5), (1, 5), (2, 5)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (1:5)  ->  (0:5) (1:5) (2:5)   via [-7, -7]
    (7, ((0, 4), (2, 3), (3, 3)), ((0, 5), (3, 3), (4, 3)), (-7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (2:3) (3:3)  ->  (0:5) (3:3) (4:3)   via [-7, -7]
    (7, ((0, 4), (3, 2), (4, 2)), ((0, 4), (2, 3), (5, 2)), (1, 1, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (3:2) (4:2)  ->  (0:4) (2:3) (5:2)   via [1, 1, -7]
    (7, ((0, 4), (3, 2), (4, 2)), ((0, 4), (3, 2), (5, 2)), (-7,), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:4) (3:2) (4:2)  ->  (0:4) (3:2) (5:2)   via [-7]
    (7, ((0, 5), (2, 4)), ((0, 3), (1, 5), (3, 4)), (-7, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:5) (2:4)  ->  (0:3) (1:5) (3:4)   via [-7, -7, -7]
    (7, ((0, 5), (3, 3)), ((0, 4), (1, 5), (4, 3)), (-7, -7, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:5) (3:3)  ->  (0:4) (1:5) (4:3)   via [-7, -7, -7]
    (7, ((0, 5), (3, 3)), ((0, 5), (5, 2)), (-6, -6, -7), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:5) (3:3)  ->  (0:5) (5:2)   via [-6, -6, -7]
    (7, ((0, 5), (4, 2)), ((0, 5), (5, 2)), (-7,), 'left'),   # 9 confirmed: window 7 arrows at the left end: (0:5) (4:2)  ->  (0:5) (5:2)   via [-7]

    (7, ((1, 2), (2, 2), (3, 4)), ((0, 2), (2, 2), (3, 4)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:2) (3:4)  ->  (0:2) (2:2) (3:4)   via [2]
    (7, ((1, 2), (2, 2), (3, 4)), ((0, 2), (2, 3), (3, 4)), (2, -8, -8), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:2) (3:4)  ->  (0:2) (2:3) (3:4)   via [2, -8, -8]
    (7, ((1, 2), (2, 2), (4, 3)), ((0, 2), (2, 2), (3, 3), (4, 3)), (2, 5, 5), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:2) (4:3)  ->  (0:2) (2:2) (3:3) (4:3)   via [2, 5, 5]
    (7, ((1, 2), (2, 3), (3, 4)), ((0, 2), (2, 3), (3, 4)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:3) (3:4)  ->  (0:2) (2:3) (3:4)   via [2]
    (7, ((1, 2), (2, 3), (4, 3)), ((0, 2), (2, 3), (4, 3)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:3) (4:3)  ->  (0:2) (2:3) (4:3)   via [2]
    (7, ((1, 2), (2, 3), (5, 2)), ((0, 2), (2, 3)), (2, -8), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:3) (5:2)  ->  (0:2) (2:3)   via [2, -8]
    (7, ((1, 2), (2, 3), (5, 2)), ((0, 2), (2, 3), (5, 2)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:3) (5:2)  ->  (0:2) (2:3) (5:2)   via [2]
    (7, ((1, 2), (2, 4), (5, 2)), ((0, 2), (2, 4)), (2, -8), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:4) (5:2)  ->  (0:2) (2:4)   via [2, -8]
    (7, ((1, 2), (2, 4), (5, 2)), ((0, 2), (2, 4), (5, 2)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:4) (5:2)  ->  (0:2) (2:4) (5:2)   via [2]
    (7, ((1, 2), (2, 5)), ((0, 2), (2, 5)), (2,), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (2:5)  ->  (0:2) (2:5)   via [2]
    (7, ((1, 2), (3, 3), (5, 2)), ((0, 2), (3, 3)), (2, -8), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:2) (3:3) (5:2)  ->  (0:2) (3:3)   via [2, -8]
    (7, ((1, 3), (2, 3), (3, 4)), ((0, 3), (1, 3), (2, 5)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:3) (3:4)  ->  (0:3) (1:3) (2:5)   via [2, 2]
    (7, ((1, 3), (2, 3), (4, 3)), ((0, 3), (1, 3), (4, 3)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:3) (4:3)  ->  (0:3) (1:3) (4:3)   via [2, 2]
    (7, ((1, 3), (2, 4), (3, 4)), ((0, 2), (2, 4)), (3, 2, 3), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:4) (3:4)  ->  (0:2) (2:4)   via [3, 2, 3]
    (7, ((1, 3), (2, 4), (4, 3)), ((0, 4), (2, 3), (3, 4)), (2, -5, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:4) (4:3)  ->  (0:4) (2:3) (3:4)   via [2, -5, 2]
    (7, ((1, 3), (2, 4), (5, 2)), ((0, 3), (1, 4), (3, 3), (4, 3)), (2, 2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:4) (5:2)  ->  (0:3) (1:4) (3:3) (4:3)   via [2, 2, 2]
    (7, ((1, 3), (2, 5)), ((0, 2), (2, 5)), (3, 2, 3), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:5)  ->  (0:2) (2:5)   via [3, 2, 3]
    (7, ((1, 3), (2, 5)), ((0, 3), (1, 5), (3, 4)), (2, 2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (2:5)  ->  (0:3) (1:5) (3:4)   via [2, 2, 2]
    (7, ((1, 3), (3, 2), (4, 3)), ((0, 3), (1, 3), (2, 3), (4, 3)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (3:2) (4:3)  ->  (0:3) (1:3) (2:3) (4:3)   via [2, 2]
    (7, ((1, 3), (3, 2), (5, 2)), ((0, 3), (1, 3), (2, 3), (4, 2)), (2, 2, 6), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (3:2) (5:2)  ->  (0:3) (1:3) (2:3) (4:2)   via [2, 2, 6]
    (7, ((1, 3), (4, 2), (5, 2)), ((0, 4), (2, 3), (3, 3), (5, 2)), (2, -5, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:3) (4:2) (5:2)  ->  (0:4) (2:3) (3:3) (5:2)   via [2, -5, 2]
    (7, ((1, 4), (2, 4), (3, 4)), ((0, 4), (1, 4), (2, 5)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (2:4) (3:4)  ->  (0:4) (1:4) (2:5)   via [2, 2]
    (7, ((1, 4), (2, 4), (4, 3)), ((0, 4), (1, 4), (3, 4)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (2:4) (4:3)  ->  (0:4) (1:4) (3:4)   via [2, 2]
    (7, ((1, 4), (2, 5)), ((0, 4), (1, 5), (4, 3)), (2, 2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (2:5)  ->  (0:4) (1:5) (4:3)   via [2, 2, 2]
    (7, ((1, 4), (3, 3), (4, 3)), ((0, 4), (1, 4), (2, 4), (3, 4)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (3:3) (4:3)  ->  (0:4) (1:4) (2:4) (3:4)   via [2, 2]
    (7, ((1, 4), (3, 3), (5, 2)), ((0, 4), (1, 4), (2, 4), (5, 2)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (3:3) (5:2)  ->  (0:4) (1:4) (2:4) (5:2)   via [2, 2]
    (7, ((1, 4), (4, 2), (5, 2)), ((0, 4), (1, 4), (3, 3)), (2, 2, -8), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (4:2) (5:2)  ->  (0:4) (1:4) (3:3)   via [2, 2, -8]
    (7, ((1, 4), (4, 2), (5, 2)), ((0, 4), (1, 4), (3, 3), (5, 2)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (4:2) (5:2)  ->  (0:4) (1:4) (3:3) (5:2)   via [2, 2]
    (7, ((1, 4), (5, 2)), ((0, 4), (1, 4), (3, 3)), (2, 6, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:4) (5:2)  ->  (0:4) (1:4) (3:3)   via [2, 6, 2]
    (7, ((1, 5), (3, 4)), ((0, 5), (1, 5), (2, 5)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:5) (3:4)  ->  (0:5) (1:5) (2:5)   via [2, 2]
    (7, ((1, 5), (4, 3)), ((0, 5), (1, 5), (3, 4)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:5) (4:3)  ->  (0:5) (1:5) (3:4)   via [2, 2]
    (7, ((1, 5), (5, 2)), ((0, 5), (1, 5), (4, 3)), (2, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (1:5) (5:2)  ->  (0:5) (1:5) (4:3)   via [2, 2]
    (7, ((2, 2), (3, 3), (5, 2)), ((0, 2), (3, 3), (5, 2)), (3, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (2:2) (3:3) (5:2)  ->  (0:2) (3:3) (5:2)   via [3, 2]
    (7, ((2, 2), (3, 4)), ((0, 2), (3, 4)), (3, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (2:2) (3:4)  ->  (0:2) (3:4)   via [3, 2]
    (7, ((3, 2), (4, 3)), ((0, 2), (4, 3)), (4, 3, 2), 'right'),   # 9 confirmed: window 7 arrows at the right end: (3:2) (4:3)  ->  (0:2) (4:3)   via [4, 3, 2]

    (8, ((0, 2), (1, 4), (4, 2)), ((0, 2), (1, 4), (6, 2)), (-7, -8), 'left'),   # 9 confirmed: window 8 arrows at the left end: (0:2) (1:4) (4:2)  ->  (0:2) (1:4) (6:2)   via [-7, -8]
    (8, ((0, 2), (1, 4), (4, 2)), ((1, 4), (6, 2)), (1, -7, -8), 'left'),   # 9 confirmed: window 8 arrows at the left end: (0:2) (1:4) (4:2)  ->  (1:4) (6:2)   via [1, -7, -8]
    (8, ((0, 3), (2, 3), (4, 2)), ((0, 3), (2, 3), (6, 2)), (-7, -8), 'left'),   # 9 confirmed: window 8 arrows at the left end: (0:3) (2:3) (4:2)  ->  (0:3) (2:3) (6:2)   via [-7, -8]
    (8, ((0, 5), (4, 2)), ((0, 5), (6, 2)), (-7, -8), 'left'),   # 9 confirmed: window 8 arrows at the left end: (0:5) (4:2)  ->  (0:5) (6:2)   via [-7, -8]

    (8, ((2, 2), (3, 3), (5, 3)), ((0, 2), (3, 3), (5, 3)), (3, 2), 'right'),   # 9 confirmed: window 8 arrows at the right end: (2:2) (3:3) (5:3)  ->  (0:2) (3:3) (5:3)   via [3, 2]
    (8, ((2, 2), (3, 4), (6, 2)), ((0, 2), (3, 4)), (3, 2, -9), 'right'),   # 9 confirmed: window 8 arrows at the right end: (2:2) (3:4) (6:2)  ->  (0:2) (3:4)   via [3, 2, -9]
    (8, ((2, 2), (3, 4), (6, 2)), ((0, 2), (3, 4), (6, 2)), (3, 2), 'right'),   # 9 confirmed: window 8 arrows at the right end: (2:2) (3:4) (6:2)  ->  (0:2) (3:4) (6:2)   via [3, 2]
    (8, ((2, 2), (3, 5)), ((0, 2), (3, 5)), (3, 2), 'right'),   # 9 confirmed: window 8 arrows at the right end: (2:2) (3:5)  ->  (0:2) (3:5)   via [3, 2]
]
