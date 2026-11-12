r"""
Checks of the hyperbolic polygons built from their angles and the cosh of
some of their sides, against hard-coded values.

Sources of the formulas: the second cosine rule, S. Katok, Fuchsian groups
(1992), Theorem 1.5.2, p. 16, for triangles and regular polygons; the
formula for right-angled hexagons, P. Buser, Geometry and spectra of
compact Riemann surfaces (1992), Theorem 2.4.1 (i) and their existence,
Theorem 2.4.2; the relation of the pentagon, Example 2.2.7 (v), and Theorem
2.3.4; the relations of the quadrilateral, Theorem 2.3.1 and Theorem 2.2.6.
The fields, angles and areas below are exact values computed with
sage-flatsurf in AA by the walk of the polygon, written over their minimal
field.
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

from sage.all import (
    AA,
    QQ,
    RDF,
    RR,
    NumberField,
    QuadraticField,
    RealField,
    polygen,
)

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
# docstring (Katok, Thm 1.5.2; Buser, Chapter 2, Theorems 2.3.1 and 2.4.1,
# Example 2.2.7).
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
    assert [P.edge_cosh_squared(i) for i in range(6)] == [25, QQ(25) / 16] * 3
    assert P._cosh2_sides() == (25, QQ(25) / 16) * 3


def test_octagon_side_of_oracle():
    K, vertices = octagon_vertices()
    r = K.gen()
    side = sum(c * r**i for i, c in enumerate(OCTAGON_SIDE_COSH2))
    P = hyperbolic_polygons.regular_ngon(8, OCTAGON["interior_angle"])
    assert {AA(P.edge_cosh_squared(i)) for i in range(8)} == {AA(side)}
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


def test_ring_conversion():
    # The triangle (1/4, 1/12, 1/12) is over y^2 - 6 by default; vertex 1 is
    # (a/3, 0) with a the square root of 6, computed with sage-flatsurf.
    K = QuadraticField(6, "a")
    a = K.gen()
    T = hyperbolic_polygons.triangle(QQ(1) / 4, QQ(1) / 12, QQ(1) / 12, ring=K)
    assert T.base_ring() is K
    assert T.vertices()[1].coordinates(model="klein") == (a / 3, 0)
    assert T.angles() == (QQ(1) / 4, QQ(1) / 12, QQ(1) / 12)


@pytest.mark.parametrize("ring", [RR, RDF, RealField(100)])
@pytest.mark.parametrize(
    "constructor, arguments",
    [
        ("triangle", (QQ(1) / 8, QQ(1) / 8, QQ(1) / 8)),
        ("regular_ngon", (8, QQ(1) / 16)),
        ("right_angled_hexagon", (QQ(5) / 4, QQ(5) / 4, QQ(5) / 4)),
        ("right_angled_pentagon", (QQ(1) / 6, QQ(5) / 4, QQ(5) / 4)),
        ("right_angled_quadrilateral", (QQ(1) / 4, QQ(1) / 8, QQ(5) / 4)),
    ],
)
def test_inexact_ring_refused(constructor, arguments, ring):
    with pytest.raises(ValueError, match="the ring must be exact, not"):
        getattr(hyperbolic_polygons, constructor)(*arguments, ring=ring)


def test_inexact_ring_message():
    with pytest.raises(
        ValueError,
        match="the ring must be exact, not Real Field with 53 bits of precision",
    ):
        hyperbolic_polygons.triangle(QQ(1) / 8, QQ(1) / 8, QQ(1) / 8, ring=RR)


@pytest.mark.parametrize("n", [QQ(11) / 2, 5.5])
def test_regular_ngon_non_integral(n):
    with pytest.raises(
        ValueError, match="the number of vertices must be an integer"
    ) as error:
        hyperbolic_polygons.regular_ngon(n, QQ(1) / 8)
    assert error.value.__cause__ is None
    assert error.value.__suppress_context__


def test_polygon_coordinates_error_suppresses_context():
    with pytest.raises(ValueError, match="do not lie in") as error:
        hyperbolic_polygons.regular_ngon(8, QQ(1) / 16, ring=QuadraticField(2, "a"))
    assert error.value.__suppress_context__


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
    assert T.parent() is HyperbolicPlane(AA)
    assert T == hyperbolic_polygons.triangle(
        QQ(1) / 4, QQ(1) / 16, QQ(1) / 16, ring=AA
    )
    with pytest.raises(ValueError):
        HyperbolicPlane(QQ).polygon(angles=[QQ(1) / 4, QQ(1) / 6, QQ(1) / 14])
    with pytest.raises(NotImplementedError):
        HyperbolicPlane(AA).polygon(angles=[QQ(1) / 8] * 4)
    with pytest.raises(ValueError):
        HyperbolicPlane(AA).polygon(angles=[QQ(1) / 8] * 2)
