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
