r"""
Checks of the orientation reversing isometries of the hyperbolic plane
(reflections and glide reflections) against hard-coded matrices, traces and
axes.
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

import re

import pytest

from sage.all import QQ, I, diagonal_matrix, dumps, loads, matrix

from flatsurf import HyperbolicPlane
from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup


FORM3 = [1, 2, -3]
REFLECTION_LIFT = re.escape("an orientation reversing isometry has no lift to SL(2, R)")


@pytest.fixture
def G():
    return HyperbolicIsometryGroup(QQ)


@pytest.fixture
def G3():
    return HyperbolicIsometryGroup(QQ, form=FORM3)


@pytest.fixture
def H():
    return HyperbolicPlane(QQ)


@pytest.fixture
def H3():
    return HyperbolicPlane(QQ, form=FORM3)


def rational(rows):
    return matrix(QQ, [[QQ(x) for x in row] for row in rows])


B = [[QQ(5) / 4, 0, QQ(3) / 4], [0, 1, 0], [QQ(3) / 4, 0, QQ(5) / 4]]
FLIP = diagonal_matrix(QQ, [1, -1, 1])


def standard_rows():
    """
    The orientation reversing isometries of the standard group: the input,
    the matrix of the isometry, its trace and whether it is a reflection.
    """
    return [
        (
            FLIP,
            rational([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),
            1,
            True,
        ),
        (
            matrix(QQ, B) * FLIP,
            matrix(
                QQ,
                [[QQ(5) / 4, 0, QQ(3) / 4], [0, -1, 0], [QQ(3) / 4, 0, QQ(5) / 4]],
            ),
            QQ(3) / 2,
            False,
        ),
        (
            matrix(QQ, [[1, 0], [0, -1]]),
            rational([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            1,
            True,
        ),
        (
            matrix(QQ, [[0, 1], [1, 0]]),
            rational([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),
            1,
            True,
        ),
        (
            matrix(QQ, [[2, 0], [0, -1]]),
            matrix(
                QQ,
                [[-1, 0, 0], [0, QQ(5) / 4, QQ(3) / 4], [0, QQ(3) / 4, QQ(5) / 4]],
            ),
            QQ(3) / 2,
            False,
        ),
        (
            matrix(QQ, [[2, 3], [1, -2]]),
            matrix(
                QQ,
                [
                    [-QQ(1) / 7, QQ(8) / 7, -QQ(4) / 7],
                    [QQ(8) / 7, -QQ(1) / 7, QQ(4) / 7],
                    [QQ(4) / 7, -QQ(4) / 7, QQ(9) / 7],
                ],
            ),
            1,
            True,
        ),
    ]


@pytest.mark.parametrize("index", range(6))
def test_standard_group(G, index):
    x, M, trace, is_reflection = standard_rows()[index]
    g = G(x)
    assert g.matrix() == M
    assert g.det() == -1
    assert g.trace() == trace
    assert not g.is_orientable()
    assert g.is_reflection() == is_reflection
    assert g.is_glide_reflection() == (not is_reflection)
    assert (g * g).is_one() == is_reflection
    assert (~g * g).is_one()


def test_other_form(G3):
    mirrors = [
        diagonal_matrix(QQ, [-1, 1, 1]),
        matrix(
            QQ,
            [
                [QQ(1) / 3, -QQ(4) / 3, 0],
                [-QQ(2) / 3, -QQ(1) / 3, 0],
                [0, 0, 1],
            ],
        ),
    ]
    for M in mirrors:
        g = G3(M)
        assert g.matrix() == M
        assert g.det() == -1
        assert g.trace() == 1
        assert g.is_reflection()
        assert not g.is_glide_reflection()
        assert not g.is_orientable()
        assert (g * g).is_one()

    glide = G3.some_elements()[3] * G3(FLIP)
    assert glide.matrix() == matrix(
        QQ,
        [
            [QQ(19) / 13, 0, -QQ(24) / 13],
            [0, -1, 0],
            [-QQ(8) / 13, 0, QQ(19) / 13],
        ],
    )
    assert glide.det() == -1
    assert glide.trace() == QQ(25) / 13
    assert glide.is_glide_reflection()
    assert not glide.is_reflection()
    assert not (glide * glide).is_one()
    assert G3(glide.matrix()) == glide


def test_orientation_preserving(G, G3):
    for g in (G.one(), G.an_element(), G(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))):
        assert g.det() == 1
        assert g.is_orientable()
        assert not g.is_reflection()
        assert not g.is_glide_reflection()
    assert G3.an_element().det() == 1


def test_axes_standard(G, H):
    point = H.point

    assert G(FLIP).axis() == H.geodesic(
        point(0, 0, model="klein"), point(QQ(1) / 2, 0, model="klein")
    ).unoriented()

    along = H.geodesic(H(I), H(QQ(3) / 5 + QQ(4) / 5 * I))
    assert G(matrix(QQ, B)).axis() == along
    assert G(matrix(QQ, B) * FLIP).axis() == along
    assert (~G(matrix(QQ, B))).axis() == H.geodesic(H(I), H(-QQ(3) / 5 + QQ(4) / 5 * I))

    assert G(matrix(QQ, [[2, 1], [1, 1]])).axis() == H.geodesic(
        H(I), H(QQ(3) / 2 + QQ(1) / 2 * I)
    )
    assert G(matrix(QQ, [[2, 0], [0, -1]])).axis() == H.geodesic(H(I), H(2 * I))


def test_axes_other_form(G3, H3):
    point = H3.point

    assert G3(diagonal_matrix(QQ, [-1, 1, 1])).axis() == H3.geodesic(
        point(0, 0, model="klein"), point(0, QQ(1) / 2, model="klein")
    ).unoriented()
    mirror = G3(
        matrix(
            QQ,
            [[QQ(1) / 3, -QQ(4) / 3, 0], [-QQ(2) / 3, -QQ(1) / 3, 0], [0, 0, 1]],
        )
    )
    assert mirror.axis() == H3.geodesic(
        point(0, 0, model="klein"), point(1, -QQ(1) / 2, model="klein")
    ).unoriented()

    along = H3.geodesic(
        point(0, 0, model="klein"), point(-QQ(24) / 19, 0, model="klein")
    )
    assert G3.some_elements()[3].axis() == along
    assert (G3.some_elements()[3] * G3(FLIP)).axis() == along


def test_refusals(G, G3):
    reflection = G(FLIP)
    glide = G(matrix(QQ, [[2, 0], [0, -1]]))
    for g in (reflection, glide):
        with pytest.raises(ValueError, match="fixes no single point; see axis"):
            g.fixed_point()
        with pytest.raises(ValueError, match=REFLECTION_LIFT):
            g._half_lift()
        with pytest.raises(ValueError, match=REFLECTION_LIFT):
            g.psl2_matrix()
    with pytest.raises(ValueError, match=REFLECTION_LIFT):
        G3(diagonal_matrix(QQ, [-1, 1, 1])).psl2_matrix()

    with pytest.raises(ValueError, match="every point is fixed by the identity"):
        G.one().axis()
    with pytest.raises(ValueError, match="an elliptic isometry has no axis"):
        G([[0, -1, 0], [1, 0, 0], [0, 0, 1]]).axis()
    with pytest.raises(ValueError, match="a parabolic isometry has no axis"):
        G(matrix(QQ, [[1, 1], [0, 1]])).axis()

    with pytest.raises(ValueError, match=re.escape("matrix is not in O(2, 1)+")):
        G(diagonal_matrix(QQ, [1, 1, -1]))
    with pytest.raises(ValueError, match="matrix does not define an isometry"):
        G(matrix(QQ, [[1, 1], [1, 1]]))


def test_pickling(G, G3):
    for g in (G(FLIP), G3(diagonal_matrix(QQ, [-1, 1, 1]))):
        assert loads(dumps(g)) == g
