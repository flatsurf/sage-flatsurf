r"""
Checks of the hyperbolic polygons built from their angles and the cosh of
some of their sides, against hard-coded values.

Sources of the formulas: the second cosine rule, S. Katok, Fuchsian groups
(1992), Theorem 1.5.2, p. 16, for triangles and regular polygons; the
formula for right-angled hexagons, P. Buser, Geometry and spectra of
compact Riemann surfaces (1992), Chapter 2, and their existence, Lemma
1.7.1, for hexagons, pentagons and quadrilaterals. The fields, angles and
areas below are exact values computed with sage-flatsurf in AA by the walk
of the polygon, written over their minimal field.
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

import pytest

from sage.all import AA, QQ, NumberField, polygen

from flatsurf import HyperbolicPlane, hyperbolic_polygons

from hyperbolic_surface_oracles import (
    OCTAGON,
    OCTAGON_SIDE_COSH2,
    TRIANGLE_DOUBLES,
    octagon_vertices,
    triangle_double_vertices,
)

y = polygen(QQ, "y")

# Rows: name, constructor, arguments, defining polynomial of the field of the
# coordinates (None for QQ, an integer when only the degree is known), angles,
# area, computed with sage-flatsurf from the formulas cited in the module
# docstring (Katok, Thm 1.5.2; Buser, Chapter 2).
PRIMITIVES = [
    (
        "rational hexagon",
        "right_angled_hexagon",
        (QQ(5) / 4, QQ(5) / 4, QQ(145) / 24),
        None,
        (QQ(1) / 4,) * 6,
        QQ(1) / 2,
    ),
    (
        "equilateral hexagon",
        "right_angled_hexagon",
        (QQ(5) / 4, QQ(5) / 4, QQ(5) / 4),
        y**2 - 6,
        (QQ(1) / 4,) * 6,
        QQ(1) / 2,
    ),
    (
        "triangle (2, 8, 8)",
        "triangle",
        (QQ(1) / 4, QQ(1) / 16, QQ(1) / 16),
        y**4 - 2 * y**2 - 1,
        (QQ(1) / 4, QQ(1) / 16, QQ(1) / 16),
        QQ(1) / 8,
    ),
    (
        "octagon",
        "regular_ngon",
        (8, QQ(1) / 8),
        y**4 - 2 * y**2 - 1,
        (QQ(1) / 8,) * 8,
        QQ(2),
    ),
    (
        "pentagon of angle 1/14",
        "right_angled_pentagon",
        (QQ(1) / 14, QQ(5) / 4, QQ(5) / 4),
        6,
        (QQ(1) / 4, QQ(1) / 14, QQ(1) / 4, QQ(1) / 4, QQ(1) / 4),
        QQ(3) / 7,
    ),
    (
        "pentagon of angle 3/8",
        "right_angled_pentagon",
        (QQ(3) / 8, QQ(5) / 4, QQ(5) / 4),
        y**4 - 56 * y**2 - 72 * y + 414,
        (QQ(1) / 4, QQ(3) / 8, QQ(1) / 4, QQ(1) / 4, QQ(1) / 4),
        QQ(1) / 8,
    ),
    (
        "quadrilateral (1/4, 1/8)",
        "right_angled_quadrilateral",
        (QQ(1) / 4, QQ(1) / 8, QQ(5) / 4),
        y**2 - 34,
        (QQ(1) / 4, QQ(1) / 4, QQ(1) / 8, QQ(1) / 4),
        QQ(1) / 8,
    ),
    (
        "quadrilateral (1/6, 1/6)",
        "right_angled_quadrilateral",
        (QQ(1) / 6, QQ(1) / 6, QQ(5) / 4),
        y**2 - 3,
        (QQ(1) / 4, QQ(1) / 6, QQ(1) / 6, QQ(1) / 4),
        QQ(1) / 6,
    ),
]


@pytest.mark.parametrize(
    "name, constructor, arguments, field, angles, area",
    PRIMITIVES,
    ids=[row[0] for row in PRIMITIVES],
)
def test_primitive(name, constructor, arguments, field, angles, area):
    P = getattr(hyperbolic_polygons, constructor)(*arguments)
    assert P.angles() == angles
    assert P.area() == area
    assert P.vertices()[0].coordinates(model="klein") == (0, 0)
    assert P.vertices()[1].coordinates(model="klein")[1] == 0
    K = P.base_ring()
    if field is None:
        assert K is QQ
    elif isinstance(field, int):
        assert K.degree() == field
    else:
        assert K.is_isomorphic(NumberField(field, "b"))


def test_rational_hexagon_vertices():
    P = hyperbolic_polygons.right_angled_hexagon(
        QQ(5) / 4, QQ(5) / 4, QQ(145) / 24
    )
    assert P.base_ring() is QQ
    assert [v.coordinates(model="klein") for v in P.vertices()] == [
        (0, 0),
        (QQ(56) / 65, 0),
        (QQ(56) / 65, QQ(99) / 325),
        (QQ(840) / 4121, QQ(4015) / 4121),
        (QQ(1344) / 9425, QQ(143) / 145),
        (0, QQ(143) / 145),
    ]


def test_equilateral_hexagon_seams():
    # The seams have cosh = 5: (5/4 + 25/16) / (9/16) = 5.
    P = hyperbolic_polygons.right_angled_hexagon(QQ(5) / 4, QQ(5) / 4, QQ(5) / 4)
    assert P._cosh2_sides() == (25, QQ(25) / 16) * 3


def test_octagon_side_of_oracle():
    K, vertices = octagon_vertices()
    r = K.gen()
    side = sum(c * r**i for i, c in enumerate(OCTAGON_SIDE_COSH2))
    P = hyperbolic_polygons.regular_ngon(8, OCTAGON["interior_angle"])
    assert {AA(c) for c in P._cosh2_sides()} == {AA(side)}
    assert P.area() == OCTAGON["area_over_2pi"]


@pytest.mark.parametrize(
    "row",
    [row for row in TRIANGLE_DOUBLES if 7 not in row[0]],
    ids=lambda row: str(row[0]),
)
def test_triangle_of_triangle_double(row):
    K, front, _ = triangle_double_vertices(row)
    angles = [QQ(1) / (2 * m) for m in row[0]]
    T = hyperbolic_polygons.triangle(*angles)
    assert T.base_ring().is_isomorphic(K)
    assert [tuple(AA(c) for c in v.coordinates(model="klein")) for v in T.vertices()] == [
        tuple(AA(c) for c in v) for v in front
    ]


def test_ring():
    T = hyperbolic_polygons.triangle(QQ(1) / 4, QQ(1) / 16, QQ(1) / 16, ring=AA)
    assert T.base_ring() is AA
    assert T.angles() == (QQ(1) / 4, QQ(1) / 16, QQ(1) / 16)
    P = hyperbolic_polygons.right_angled_hexagon(5, 5, 5, ring=AA)
    assert P.base_ring() is AA


@pytest.mark.parametrize(
    "constructor, arguments",
    [
        ("triangle", (QQ(1) / 4, QQ(1) / 8, QQ(1) / 8)),
        ("triangle", (0, QQ(1) / 8, QQ(1) / 8)),
        ("triangle", (QQ(1) / 4, QQ(1) / 16, QQ(1) / 16, QQ)),
        ("regular_ngon", (4, QQ(1) / 4)),
        ("regular_ngon", (2, QQ(1) / 8)),
        ("right_angled_hexagon", (1, 2, 2)),
        ("right_angled_pentagon", (QQ(3) / 8, QQ(145) / 24, QQ(5) / 4)),
        ("right_angled_pentagon", (QQ(1) / 3, 3, QQ(5) / 4)),
        ("right_angled_pentagon", (QQ(1) / 2, 2, 2)),
        ("right_angled_quadrilateral", (QQ(3) / 8, QQ(1) / 16, QQ(145) / 24)),
        ("right_angled_quadrilateral", (QQ(1) / 4, QQ(1) / 4, 2)),
    ],
)
def test_refused(constructor, arguments):
    with pytest.raises(ValueError):
        getattr(hyperbolic_polygons, constructor)(*arguments)


def test_existing_with_obtuse_angle():
    P = hyperbolic_polygons.right_angled_pentagon(QQ(1) / 3, 2, QQ(5) / 4)
    assert P.angles() == (QQ(1) / 4, QQ(1) / 3, QQ(1) / 4, QQ(1) / 4, QQ(1) / 4)
    Q = hyperbolic_polygons.right_angled_quadrilateral(
        QQ(3) / 8, QQ(1) / 16, QQ(5) / 4
    )
    assert Q.angles() == (QQ(1) / 4, QQ(3) / 8, QQ(1) / 16, QQ(1) / 4)


def test_plane_polygon_from_angles():
    T = HyperbolicPlane(AA).polygon(angles=[QQ(1) / 4, QQ(1) / 16, QQ(1) / 16])
    assert T.angles() == (QQ(1) / 4, QQ(1) / 16, QQ(1) / 16)
    with pytest.raises(ValueError):
        HyperbolicPlane(QQ).polygon(angles=[QQ(1) / 4, QQ(1) / 6, QQ(1) / 14])
    with pytest.raises(NotImplementedError):
        HyperbolicPlane(AA).polygon(angles=[QQ(1) / 8] * 4)
    with pytest.raises(ValueError):
        HyperbolicPlane(AA).polygon(angles=[QQ(1) / 8] * 2)
