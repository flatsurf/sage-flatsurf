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


def _both_geometries():
    r"""
    Return the rings, hyperbolic planes, vertices, angles and points used by
    the tests of ``Polygon`` in both geometries, by name.
    """
    from sage.all import AA, QQ, RR, QuadraticField

    from flatsurf import HyperbolicPlane

    half = QQ(1) / 2
    K = QuadraticField(2, "s")
    s = K.gen()

    H1 = HyperbolicPlane(QQ)
    HK = HyperbolicPlane(K)
    HA = HyperbolicPlane(AA)
    HR = HyperbolicPlane(RR)

    V = [(0, 0), (half, 0), (0, half)]
    W = [(0, 0), (s / 4, 0), (0, s / 4)]
    A = (QQ(1) / 6, QQ(1) / 8, QQ(1) / 8)

    p0, p1, p2 = (H1.point(x, y, model="klein") for (x, y) in V)
    pK = HK.point(s / 4, 0, model="klein")
    qK = HK.point(0, s / 4, model="klein")
    r0, r1, r2 = (HR.point(x, y, model="klein") for (x, y) in V)

    return {
        "QQ": QQ,
        "AA": AA,
        "K": K,
        "s": s,
        "H1": H1,
        "HK": HK,
        "HA": HA,
        "HR": HR,
        "V": V,
        "W": W,
        "A": A,
        "p": (p0, p1, p2),
        "pK": pK,
        "qK": qK,
        "r": (r0, r1, r2),
    }


def _hyperbolic_equal():
    r"""
    Return the calls of ``Polygon`` that give the same polygon as
    ``HyperbolicPlane.polygon`` or ``hyperbolic_polygons.triangle``.
    """
    from flatsurf import Polygon, hyperbolic_polygons

    d = _both_geometries()
    H1, V, A, QQ = d["H1"], d["V"], d["A"], d["QQ"]
    p0, p1, p2 = d["p"]
    half = QQ(1) / 2

    def expected():
        return H1.polygon(vertices=V)

    return {
        "geometry": (lambda: Polygon(vertices=V, geometry="hyperbolic"), expected),
        "Geometry": (lambda: Polygon(vertices=V, geometry="Hyperbolic"), expected),
        "angles": (
            lambda: Polygon(angles=A, geometry="hyperbolic"),
            lambda: hyperbolic_polygons.triangle(*A),
        ),
        "parent": (lambda: Polygon(vertices=V, parent=H1), expected),
        "parent-geometry": (
            lambda: Polygon(vertices=V, parent=H1, geometry="hyperbolic"),
            expected,
        ),
        "parent-base_ring": (
            lambda: Polygon(vertices=V, parent=H1, base_ring=QQ),
            expected,
        ),
        "points": (lambda: Polygon(vertices=[p0, p1, p2]), expected),
        "points-numbers": (lambda: Polygon(vertices=[p0, (half, 0), p2]), expected),
    }


@pytest.mark.parametrize("name", list(_hyperbolic_equal()))
def test_polygon_hyperbolic_equal(name):
    from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

    call, expected = _hyperbolic_equal()[name]
    P, Q = call(), expected()
    assert isinstance(P, HyperbolicPolygon)
    assert P == Q
    assert P.parent() is Q.parent()


def _hyperbolic_in():
    r"""
    Return the calls of ``Polygon`` with the name of the plane of the result.
    """
    from flatsurf import Polygon

    d = _both_geometries()
    V, W, A, AA, QQ, K, s = d["V"], d["W"], d["A"], d["AA"], d["QQ"], d["K"], d["s"]
    HK, HA, HR = d["HK"], d["HA"], d["HR"]
    p0, p1, p2 = d["p"]
    r0, r1, r2 = d["r"]
    pK, qK = d["pK"], d["qK"]
    half = QQ(1) / 2

    return {
        "W": (lambda: Polygon(vertices=W, geometry="hyperbolic"), "HK"),
        "V-AA": (lambda: Polygon(vertices=V, geometry="hyperbolic", base_ring=AA), "HA"),
        "angles-AA": (
            lambda: Polygon(angles=A, geometry="hyperbolic", base_ring=AA),
            "HA",
        ),
        "half_plane": (
            lambda: Polygon(
                vertices=[(0, 1), (1, 1), (0, 2)],
                geometry="hyperbolic",
                model="half_plane",
            ),
            "H1",
        ),
        "non-convex": (
            lambda: Polygon(
                vertices=[(0, -half), (half, half), (0, 0), (-half, half)],
                geometry="hyperbolic",
            ),
            "H1",
        ),
        "parent-K": (lambda: Polygon(vertices=V, parent=HK), "HK"),
        "parent-angles-AA": (lambda: Polygon(angles=A, parent=HA), "HA"),
        "parent-K-points": (lambda: Polygon(vertices=[p0, p1, p2], parent=HK), "HK"),
        "parent-RR": (lambda: Polygon(vertices=V, parent=HR), "HR"),
        "parent-RR-points": (lambda: Polygon(vertices=[p0, p1, p2], parent=HR), "HR"),
        "points-K": (lambda: Polygon(vertices=[p0, pK, qK]), "HK"),
        "points-base_ring-K": (
            lambda: Polygon(vertices=[p0, p1, p2], base_ring=K),
            "HK",
        ),
        "points-RR": (lambda: Polygon(vertices=[r0, r1, r2]), "HR"),
        "points-mixed-RR": (lambda: Polygon(vertices=[p0, r1, r2]), "HR"),
        "points-numbers-K": (lambda: Polygon(vertices=[p0, (s / 4, 0), p2]), "HK"),
    }


@pytest.mark.parametrize("name", list(_hyperbolic_in()))
def test_polygon_hyperbolic_in(name):
    from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

    call, plane = _hyperbolic_in()[name]
    P = call()
    assert isinstance(P, HyperbolicPolygon)
    assert P.parent() is _both_geometries()[plane]


def test_polygon_hyperbolic_half_plane_vertices():
    from sage.all import I

    from flatsurf import Polygon

    H1 = _both_geometries()["H1"]
    P = Polygon(
        vertices=[(0, 1), (1, 1), (0, 2)], geometry="hyperbolic", model="half_plane"
    )
    assert list(P.vertices()) == [H1(I), H1(1 + I), H1(2 * I)]


def test_polygon_hyperbolic_angles():
    from flatsurf import Polygon

    d = _both_geometries()
    A, AA, HA = d["A"], d["AA"], d["HA"]
    assert Polygon(angles=A, geometry="hyperbolic", base_ring=AA).angles() == A
    assert Polygon(angles=A, parent=HA).angles() == A


def test_polygon_euclidean_unchanged():
    from sage.all import QQ

    from flatsurf import Polygon
    from flatsurf.geometry.polygon import EuclideanPolygon

    d = _both_geometries()
    V, A = d["V"], d["A"]
    P = Polygon(vertices=V)
    assert isinstance(P, EuclideanPolygon)
    assert Polygon(vertices=V, geometry="euclidean") == P
    assert Polygon(vertices=V, geometry="Euclidean") == P
    Q = Polygon(angles=A)
    assert isinstance(Q, EuclideanPolygon)
    assert Q.angles() == (QQ(1) / 5, QQ(3) / 20, QQ(3) / 20)


def _hyperbolic_errors():
    import re

    from sage.all import QQ, RR

    from flatsurf import Polygon
    from flatsurf.geometry.categories import EuclideanPolygons

    d = _both_geometries()
    V, W, A, AA = d["V"], d["W"], d["A"], d["AA"]
    H1, HR = d["H1"], d["HR"]
    p0, p1, p2 = d["p"]
    pK, qK = d["pK"], d["qK"]
    half = QQ(1) / 2
    edges = [(1, 0), (-1, 1), (0, -1)]
    not_in_QQ = "do not lie in Rational Field"
    inexact = "the ring must be exact, not Real"

    def hyp(**kwds):
        return lambda: Polygon(geometry="hyperbolic", **kwds)

    return {
        "model-euclidean": (
            lambda: Polygon(vertices=V, model="klein"),
            ValueError,
            "model can only be used with a hyperbolic polygon",
        ),
        "W-QQ": (
            hyp(vertices=W, base_ring=QQ),
            TypeError,
            re.escape("Unable to coerce 1/4*s to a rational"),
        ),
        "angles-QQ": (hyp(angles=A, base_ring=QQ), ValueError, not_in_QQ),
        "angles-sum": (
            hyp(angles=(1, 2, 3)),
            ValueError,
            "must sum to less than 1/2",
        ),
        "angles-4": (
            hyp(angles=(QQ(1) / 8,) * 4),
            NotImplementedError,
            "only triangles",
        ),
        "infinite": (
            hyp(vertices=[(0, 0), (2, 0), (0, 2)]),
            ValueError,
            re.escape("vertex (2, 0) is not a finite point"),
        ),
        "float": (
            hyp(vertices=[(0, 0), (0.5, 0), (0, 0.5)]),
            ValueError,
            "the ring must be exact, not Real Double Field",
        ),
        "base_ring-RR": (hyp(vertices=V, base_ring=RR), ValueError, inexact),
        "negative": (
            hyp(vertices=[(0, 0), (0, half), (half, 0)]),
            ValueError,
            "polygon has negative area",
        ),
        "poincare": (
            hyp(vertices=V, model="poincare"),
            NotImplementedError,
            "unsupported model 'poincare'",
        ),
        "parent-angles-QQ": (
            lambda: Polygon(angles=A, parent=H1),
            ValueError,
            not_in_QQ,
        ),
        "parent-angles-RR": (
            lambda: Polygon(angles=A, parent=HR),
            ValueError,
            inexact,
        ),
        "parent-base_ring": (
            lambda: Polygon(vertices=V, parent=H1, base_ring=AA),
            ValueError,
            "base_ring is Algebraic Real Field but parent is Hyperbolic Plane over Rational Field",
        ),
        "points-infinite": (
            lambda: Polygon(vertices=[H1(0), p1, p2]),
            ValueError,
            "vertex 0 is not a finite point",
        ),
        "a": (
            lambda: Polygon(vertices=V, geometry="spherical"),
            ValueError,
            re.escape("geometry must be \"euclidean\" or \"hyperbolic\", not 'spherical'"),
        ),
        "a-object": (
            lambda: Polygon(vertices=V, geometry=H1.geometry),
            ValueError,
            "not Exact geometry over Rational Field",
        ),
        "b": (
            lambda: Polygon(vertices=V, geometry="euclidean", parent=H1),
            ValueError,
            re.escape('parent is a hyperbolic plane but geometry is "euclidean"'),
        ),
        "c": (
            lambda: Polygon(vertices=V, parent=QQ),
            ValueError,
            "parent must be a HyperbolicPlane, not Rational Field",
        ),
        "d": (
            lambda: Polygon(vertices=[p0, p1, p2], geometry="euclidean"),
            ValueError,
            re.escape('the vertices are hyperbolic points but geometry is "euclidean"'),
        ),
        "e": (
            lambda: Polygon(vertices=[p0, pK, qK], parent=H1),
            ValueError,
            "the vertices do not convert into Hyperbolic Plane over Rational Field",
        ),
        "f": (
            hyp(edges=edges),
            ValueError,
            re.escape("a hyperbolic polygon has no edge vectors; give its vertices"),
        ),
        "f-parent": (
            lambda: Polygon(edges=edges, parent=H1),
            ValueError,
            re.escape("a hyperbolic polygon has no edge vectors; give its vertices"),
        ),
        "lengths": (
            hyp(angles=A, lengths=[1]),
            NotImplementedError,
            "a hyperbolic polygon cannot be built from lengths",
        ),
        "category": (
            hyp(vertices=V, category=EuclideanPolygons(QQ)),
            ValueError,
            "the category of a hyperbolic polygon cannot be chosen",
        ),
        "nothing": (
            hyp(),
            NotImplementedError,
            "one of vertices or angles must be set",
        ),
        "both": (
            hyp(vertices=V, angles=A),
            ValueError,
            "a hyperbolic polygon is built from vertices or from angles, not both",
        ),
        "angles-model": (
            hyp(angles=A, model="half_plane"),
            ValueError,
            "model can only be used with vertices",
        ),
    }


@pytest.mark.parametrize("name", list(_hyperbolic_errors()))
def test_polygon_hyperbolic_errors(name):
    call, error, message = _hyperbolic_errors()[name]
    with pytest.raises(error, match=message):
        call()


def test_polygon_hyperbolic_no_common_base_ring():
    from sage.all import AA, QQ, NumberField, QuadraticField, polygen

    from flatsurf import HyperbolicPlane, Polygon

    p0 = HyperbolicPlane(QQ).point(0, 0, model="klein")

    K = QuadraticField(2, "s")
    L = QuadraticField(3, "t")
    pK = HyperbolicPlane(K).point(K.gen() / 4, 0, model="klein")
    qL = HyperbolicPlane(L).point(0, L.gen() / 4, model="klein")
    with pytest.raises(
        ValueError,
        match="the vertices lie in hyperbolic planes with no common base ring",
    ):
        Polygon(vertices=[p0, pK, qL])

    x = polygen(QQ)
    K = NumberField(x**2 - 2, "s", embedding=AA(2).sqrt())
    L = NumberField(x**2 - 3, "t", embedding=AA(3).sqrt())
    pK = HyperbolicPlane(K).point(K.gen() / 4, 0, model="klein")
    qL = HyperbolicPlane(L).point(0, L.gen() / 4, model="klein")
    assert Polygon(vertices=[p0, pK, qL]).parent() is HyperbolicPlane(AA)
