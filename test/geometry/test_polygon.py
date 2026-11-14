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


@pytest.mark.filterwarnings("ignore::UserWarning")
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
        "regular_ngon-RR",
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
