r"""
Tests for similarities of the plane: matrix input, reflections and
orientation-reversing similarities from vectors.
"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2025 Vincent Delecroix
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
# ********************************************************************

import pytest

from sage.all import QQ, RR, ZZ, NumberField, matrix, vector, polygen, VectorSpace

from flatsurf.geometry.similarity import (
    SimilarityGroup,
    reflection,
    similarity_from_vectors,
)


@pytest.mark.parametrize(
    "p, q, rows",
    [
        ((0, 0), (1, 1), [0, 1, 0, 1, 0, 0, 0, 0, 1]),
        ((1, 0), (1, 1), [-1, 0, 2, 0, 1, 0, 0, 0, 1]),
        (
            (0, 0),
            (2, 1),
            [QQ(3) / 5, QQ(4) / 5, 0, QQ(4) / 5, -QQ(3) / 5, 0, 0, 0, 1],
        ),
        (
            (1, 2),
            (3, 3),
            [
                QQ(3) / 5,
                QQ(4) / 5,
                -QQ(6) / 5,
                QQ(4) / 5,
                -QQ(3) / 5,
                QQ(12) / 5,
                0,
                0,
                1,
            ],
        ),
    ],
)
def test_reflection(p, q, rows):
    R = reflection(p, q)
    assert R == matrix(QQ, 3, rows)
    assert R.is_immutable()
    assert R.det() == -1
    assert R * R == 1
    assert R * vector(QQ, (p[0], p[1], 1)) == vector(QQ, (p[0], p[1], 1))
    assert R * vector(QQ, (q[0], q[1], 1)) == vector(QQ, (q[0], q[1], 1))
    assert SimilarityGroup(QQ)(R).is_reflection()


def test_reflection_number_field():
    x = polygen(QQ)
    K = NumberField(x**2 - 3, "a", embedding=RR(1.7))
    a = K.gen()
    R = reflection((0, 0), (1, a))
    assert R == matrix(K, 3, [-QQ(1) / 2, a / 2, 0, a / 2, QQ(1) / 2, 0, 0, 0, 1])
    assert R.base_ring() is K
    assert SimilarityGroup(K)(R).is_reflection()


def test_reflection_ring():
    assert reflection((0, 0), (1, 1)).base_ring() is QQ
    assert reflection((0, 0), (1, 1), ring=QQ).base_ring() is QQ


def test_reflection_equal_points():
    with pytest.raises(ValueError):
        reflection((1, 1), (1, 1))


@pytest.mark.parametrize(
    "u, v, rows",
    [
        ((1, 0), (0, 1), [0, 1, 1, 0]),
        ((1, 0), (1, 0), [1, 0, 0, -1]),
        ((5, 0), (3, 4), [QQ(3) / 5, QQ(4) / 5, QQ(4) / 5, -QQ(3) / 5]),
        ((2, 1), (1, -2), [QQ(4) / 5, -QQ(3) / 5, -QQ(3) / 5, -QQ(4) / 5]),
    ],
)
def test_similarity_from_vectors_reversing(u, v, rows):
    V = VectorSpace(QQ, 2)
    u, v = V(u), V(v)
    m = similarity_from_vectors(u, v, orientation_preserving=False)
    assert m == matrix(QQ, 2, rows)
    assert m.det() == -1
    assert m * u == v


def test_similarity_from_vectors_preserving_unchanged():
    V = VectorSpace(QQ, 2)
    u = V((1, 0))
    assert similarity_from_vectors(u, u) == 1
    assert similarity_from_vectors(u, V((0, 1))) == matrix(QQ, 2, [0, -1, 1, 0])


def test_matrix_input_round_trip():
    G = SimilarityGroup(QQ)
    g = G((1, -QQ(2) / 3, 1, 1, -1))
    assert g.matrix() == matrix(QQ, 3, [1, -QQ(2) / 3, 1, -QQ(2) / 3, -1, 1, 0, 0, 1])
    assert G(g.matrix()) == g


def test_matrix_input_integer_entries():
    G = SimilarityGroup(QQ)
    assert G(matrix(ZZ, 2, [0, -1, 1, 0])) == G((0, 1, 0, 0))
    assert G(matrix(ZZ, 3, [0, 1, 0, 1, 0, 0, 0, 0, 1])) == G((0, 1, 0, 0, -1))
    assert G(matrix(ZZ, 3, [0, 1, 0, 1, 0, 0, 0, 0, 1])).parent() is G


def test_matrix_input_refused():
    G = SimilarityGroup(QQ)
    with pytest.raises(ValueError, match="not a similarity matrix"):
        G(matrix(QQ, 3, [1, 0, 0, 0, 1, 0, 1, 0, 1]))
    with pytest.raises(ValueError, match="not a similarity matrix"):
        G(matrix(QQ, 3, [1, 1, 0, 0, 1, 0, 0, 0, 1]))
