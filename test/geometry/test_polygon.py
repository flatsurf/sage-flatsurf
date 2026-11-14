r"""
Test properties of polygons.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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


@pytest.mark.parametrize("n", [4, 6, 8, 10])
def test_get_point_position(n):
    from sage.all import vector

    from flatsurf import polygons

    inner = polygons.regular_ngon(n)
    inner = inner.translate(-inner.centroid())
    outer = inner * 2

    vertices = []
    for i in range(0, n, 2):
        vertices.append(inner.vertex(i))
        vertices.append(outer.vertex(i))
        vertices.append(outer.vertex(i + 1))
        vertices.append(inner.vertex(i + 1))

    from flatsurf import Polygon

    P = Polygon(vertices=vertices)

    for v in P.vertices():
        assert P.get_point_position(v).is_vertex()

    for i in range(len(P.vertices())):
        p = (P.vertex(i) + P.vertex(i + 1)) / 2
        assert P.get_point_position(p).is_in_edge_interior()

    for v in P.vertices():
        p = v + vector((1 / 1024, 0))
        p_position = P.get_point_position(p)

        q = v + vector((1 / 2048, 0))
        q_position = P.get_point_position(q)

        assert str(p_position) == str(q_position)

    for v in P.vertices():
        p = v - vector((1 / 1024, 0))
        p_position = P.get_point_position(p)

        q = v - vector((1 / 2048, 0))
        q_position = P.get_point_position(q)

        assert str(p_position) == str(q_position)


def _edge_polygons():
    from flatsurf import Polygon, polygons

    return [
        polygons.square(),
        Polygon(vertices=[(0, 0), (2, 0), (1, 1), (0, 2)]),
    ]


@pytest.mark.parametrize("P", _edge_polygons())
def test_edge_segment(P):
    n = len(P.vertices())

    for i in [-1] + list(range(n)):
        segment = P.edge_segment(i)
        assert segment.vector() == P.edge_vector(i)
        assert segment.start() == P.vertex(i)
        assert segment.end() == P.vertex(i + 1)

        reverse = -segment
        assert reverse.start() == P.vertex(i + 1)
        assert reverse.end() == P.vertex(i)

    assert P.edge_segments() == [P.edge_segment(i) for i in range(n)]
    assert [s.vector() for s in P.edge_segments()] == P.edge_vectors()


def test_edge_vectors():
    from flatsurf import Polygon

    P = Polygon(vertices=[(0, 0), (2, 0), (1, 1), (0, 2)])
    assert [tuple(v) for v in P.edge_vectors()] == [
        (2, 0),
        (-1, 1),
        (-1, 1),
        (0, -2),
    ]


@pytest.mark.parametrize("P", _edge_polygons())
def test_edge_deprecated(P):
    with pytest.warns(UserWarning, match="edge_vector"):
        assert P.edge(0) == P.edge_vector(0)

    with pytest.warns(UserWarning, match="edge_vector"):
        assert P.edges() == P.edge_vectors()


def _inexact_polygons():
    from sage.all import RDF, RR

    from flatsurf import Polygon, polygons

    return {
        "vertices-RDF": lambda: Polygon(
            vertices=[(0, 0), (1, 0), (0, 1)], base_ring=RDF
        ),
        "vertices-RDF-unchecked": lambda: Polygon(
            vertices=[(0, 0), (1, 0), (0, 1)], base_ring=RDF, check=False
        ),
        "edges-RR": lambda: Polygon(edges=[(1, 0), (-1, 1), (0, -1)], base_ring=RR),
        "vertices-RR-computed": lambda: Polygon(vertices=[(0, 0), (1.0, 0), (0, 1)]),
        "vertices-RDF-computed": lambda: Polygon(
            vertices=[(0, 0), (float(1), 0), (0, 1)]
        ),
        "angles-RR": lambda: Polygon(angles=(1, 1, 1), base_ring=RR),
        "change_ring-RDF": lambda: polygons.square().change_ring(RDF),
        "square-RDF": lambda: polygons.square(base_ring=RDF),
        "regular_ngon-RR": lambda: polygons.regular_ngon(8, field=RR),
    }


@pytest.mark.parametrize(
    "name",
    [
        "vertices-RDF",
        "vertices-RDF-unchecked",
        "edges-RR",
        "vertices-RR-computed",
        "vertices-RDF-computed",
        "angles-RR",
        "change_ring-RDF",
        "square-RDF",
        pytest.param(
            "regular_ngon-RR",
            marks=pytest.mark.filterwarnings("ignore::UserWarning"),
        ),
    ],
)
def test_inexact_ring_refused(name):
    with pytest.raises(ValueError, match="the ring must be exact, not"):
        _inexact_polygons()[name]()


def test_exact_rings_accepted():
    from sage.all import AA, QQ

    from flatsurf import Polygon, polygons

    assert polygons.square().change_ring(AA).base_ring() is AA
    assert (
        Polygon(vertices=[(0, 0), (1.0, 0), (0, 1)], base_ring=QQ).base_ring() is QQ
    )


def _euclidean_pairs():
    from sage.all import AA, QuadraticField

    from flatsurf import euclidean_polygons, polygons

    K = QuadraticField(2, "sqrt2")
    sqrt2 = K.gen()

    return {
        "square": (lambda: euclidean_polygons.square(), lambda: polygons.square()),
        "square-AA": (
            lambda: euclidean_polygons.square(base_ring=AA),
            lambda: polygons.square(base_ring=AA),
        ),
        "rectangle": (
            lambda: euclidean_polygons.rectangle(1, sqrt2),
            lambda: polygons.rectangle(1, sqrt2),
        ),
        "triangle": (
            lambda: euclidean_polygons.triangle(3, 4, 5),
            lambda: polygons.triangle(3, 4, 5),
        ),
        "regular_ngon": (
            lambda: euclidean_polygons.regular_ngon(8),
            lambda: polygons.regular_ngon(8),
        ),
        "regular_ngon-AA": (
            lambda: euclidean_polygons.regular_ngon(3, base_ring=AA),
            lambda: polygons.regular_ngon(3, field=AA),
        ),
    }


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize(
    "name",
    [
        "square",
        "square-AA",
        "rectangle",
        "triangle",
        "regular_ngon",
        "regular_ngon-AA",
    ],
)
def test_euclidean_polygons_agree_with_polygons(name):
    new, old = _euclidean_pairs()[name]
    new, old = new(), old()
    assert new == old
    assert new.category() == old.category()


def test_euclidean_triangle_angles():
    from sage.all import QQ

    from flatsurf import euclidean_polygons

    assert euclidean_polygons.triangle(1, 2, 3).angles() == (
        QQ(1) / 12,
        QQ(1) / 6,
        QQ(1) / 4,
    )
    assert euclidean_polygons.triangle(
        QQ(1) / 12, QQ(1) / 6, QQ(1) / 4
    ) == euclidean_polygons.triangle(1, 2, 3)
    assert euclidean_polygons.triangle(2, 5, 7).angles() == (
        QQ(1) / 14,
        QQ(5) / 28,
        QQ(1) / 4,
    )
    assert euclidean_polygons.triangle(QQ(1) / 6, QQ(1) / 8, QQ(1) / 8).angles() == (
        QQ(1) / 5,
        QQ(3) / 20,
        QQ(3) / 20,
    )


def _euclidean_errors():
    from sage.all import QQ, RDF, RR, QuadraticField

    from flatsurf import euclidean_polygons as E

    sqrt2 = QuadraticField(2, "sqrt2").gen()
    not_in_QQ = "the coordinates of the polygon do not lie in Rational Field"
    inexact = "the ring must be exact, not"

    return {
        "square-0": (lambda: E.square(0), ValueError, "side must be positive"),
        "square--1": (lambda: E.square(-1), ValueError, "side must be positive"),
        "square-RDF": (
            lambda: E.square(base_ring=RDF),
            ValueError,
            "the ring must be exact, not Real Double Field",
        ),
        "square-QQ": (lambda: E.square(sqrt2, base_ring=QQ), ValueError, not_in_QQ),
        "rectangle-width": (
            lambda: E.rectangle(0, 1),
            ValueError,
            "width must be positive",
        ),
        "rectangle-height": (
            lambda: E.rectangle(1, 0),
            ValueError,
            "height must be positive",
        ),
        "triangle-0": (
            lambda: E.triangle(0, 1, 1),
            ValueError,
            "angles must be positive rationals",
        ),
        "triangle-float": (lambda: E.triangle(1.0, 2.0, 3.0), TypeError, None),
        "triangle-QQ": (lambda: E.triangle(1, 2, 3, base_ring=QQ), ValueError, not_in_QQ),
        "triangle-RR": (
            lambda: E.triangle(1, 2, 3, base_ring=RR),
            ValueError,
            inexact,
        ),
        "regular_ngon-QQ": (
            lambda: E.regular_ngon(8, base_ring=QQ),
            ValueError,
            not_in_QQ,
        ),
        "regular_ngon-RR": (
            lambda: E.regular_ngon(8, base_ring=RR),
            ValueError,
            inexact,
        ),
        "regular_ngon-2": (
            lambda: E.regular_ngon(2),
            ValueError,
            "a polygon has at least 3 vertices",
        ),
        "regular_ngon-float": (
            lambda: E.regular_ngon(5.5),
            ValueError,
            "the number of vertices must be an integer",
        ),
        "regular_ngon-rational": (
            lambda: E.regular_ngon(QQ(7) / 2),
            ValueError,
            "the number of vertices must be an integer",
        ),
    }


@pytest.mark.parametrize("name", list(_euclidean_errors()))
def test_euclidean_polygons_errors(name):
    call, error, message = _euclidean_errors()[name]
    with pytest.raises(error, match=message):
        call()


def test_euclidean_polygons_base_ring():
    from sage.all import AA, QQ, QuadraticField

    from flatsurf import euclidean_polygons as E

    K = QuadraticField(2, "sqrt2")
    sqrt2 = K.gen()

    assert E.triangle(1, 2, 3, base_ring=AA) == E.triangle(1, 2, 3)
    assert E.triangle(1, 2, 3, base_ring=AA).base_ring() is AA
    assert E.rectangle(1, sqrt2, base_ring=AA).base_ring() is AA

    P = E.regular_ngon(8, base_ring=K)
    assert P == E.regular_ngon(8)
    assert P.base_ring() is K

    for P in [E.regular_ngon(4), E.regular_ngon(4, base_ring=QQ)]:
        assert P.base_ring() is QQ
        assert [tuple(v) for v in P.vertices()] == [(0, 0), (1, 0), (1, 1), (0, 1)]
