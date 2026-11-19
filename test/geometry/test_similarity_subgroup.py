r"""
Checks of the finitely generated subgroups of the group of similarities of the
plane, against hard-coded generators.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Vincent Delecroix
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************

import pickle
import re

import pytest

import sage.misc.sage_unittest

from sage.all import AA, QQ, matrix

from flatsurf.geometry.similarity import SimilarityGroup, SimilaritySubgroup

MEMBERSHIP = re.escape("membership in a finitely generated subgroup is not decided")
SKIP = ["_test_an_element", "_test_one", "_test_some_elements"]


@pytest.fixture
def F():
    return SimilarityGroup(QQ)


@pytest.fixture
def rot(F):
    return F(matrix(QQ, 2, [0, -1, 1, 0]))


@pytest.fixture
def tr(F):
    return F(1, 0, 1, 0)


@pytest.fixture
def r(F):
    return F(matrix(QQ, 2, [1, 0, 0, -1]))


@pytest.fixture
def A(F, rot, tr):
    return F.subgroup([rot, tr])


@pytest.fixture
def B(F, rot, r):
    return F.subgroup([rot, r])


def test_generators(F, A, B, rot, tr, r):
    assert isinstance(A, SimilaritySubgroup)
    assert A.gens() == (rot, tr)
    assert B.gens() == (rot, r)
    assert A.ngens() == B.ngens() == 2
    assert A.ambient() is F
    assert B.ambient() is F
    assert A.base_ring() is QQ
    assert repr(B) == f"Subgroup with 2 generators of {F!r}"


def test_orientable(F, A, B, rot, tr, r):
    assert A.is_orientable()
    assert not B.is_orientable()
    assert not F.subgroup([r]).is_orientable()
    assert F.subgroup([]).is_orientable()


def test_orientation_preserving_subgroup(F, A, B, rot, r):
    assert A.orientation_preserving_subgroup() is A

    P = B.orientation_preserving_subgroup()
    assert P.ngens() == 1
    assert P.gens() == (rot,)
    assert P.is_orientable()
    assert P.orientation_preserving_subgroup() is P

    # The reflection in the line y = -x, and the generators t r and r t
    # left once r r and r ~r are dropped, with r t the inverse of t r.
    t = F(matrix(QQ, 2, [0, -1, -1, 0]))
    Q = F.subgroup([r, t]).orientation_preserving_subgroup()
    assert Q.gens() == (t * r,)


def test_trivial_subgroup(F):
    T = F.subgroup([])
    assert T.gens() == ()
    assert T.ngens() == 0
    assert T.is_orientable()
    assert T.orientation_preserving_subgroup() is T
    assert T.one().is_one()
    assert list(T.matrix_group().gens())[0].matrix() == matrix(QQ, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])
    sage.misc.sage_unittest.TestSuite(T).run(skip=SKIP, raise_on_failure=True)


def test_change_ring(F, B):
    G = F.change_ring(AA)
    assert G is SimilarityGroup(AA)
    assert G.change_ring(QQ) is F
    BA = B.change_ring(AA)
    assert BA.ambient() is G
    assert BA.change_ring(QQ) is B


def test_unique_representation_and_pickling(F, A, B, rot, tr, r):
    assert F.subgroup([rot, tr]) is A
    assert F.subgroup((rot, r)) is B
    assert F.subgroup([r, rot]) is not B
    assert pickle.loads(pickle.dumps(B)) is B
    assert pickle.loads(pickle.dumps(A)) is A


def test_membership_refused(A, rot):
    with pytest.raises(NotImplementedError, match=MEMBERSHIP):
        rot in A
    with pytest.raises(NotImplementedError, match=MEMBERSHIP):
        A(rot)


def test_matrix_group(B):
    M = B.matrix_group()
    assert M.base_ring() is QQ
    assert [g.matrix() for g in M.gens()] == [
        matrix(QQ, 3, [0, -1, 0, 1, 0, 0, 0, 0, 1]),
        matrix(QQ, 3, [1, 0, 0, 0, -1, 0, 0, 0, 1]),
    ]


def test_test_suite(A, B):
    sage.misc.sage_unittest.TestSuite(A).run(skip=SKIP, raise_on_failure=True)
    sage.misc.sage_unittest.TestSuite(B).run(skip=SKIP, raise_on_failure=True)
