r"""
Test properties of polygons.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023-2025 Julian Rüth
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
    outer = inner.base_ring()(2) * inner

    vertices = []
    for i in range(0, n, 2):
        vertices.append(inner.corner(i).vector())
        vertices.append(outer.corner(i).vector())
        vertices.append(outer.corner(i + 1).vector())
        vertices.append(inner.corner(i + 1).vector())

    from flatsurf import Polygon

    P = Polygon(vertices=vertices)

    for v in P.corners():
        assert P.get_point_position(v.vector()).is_vertex()

    for i in range(len(P.corners())):
        p = (P.corner(i).vector() + P.corner(i + 1).vector()) / 2
        assert P.get_point_position(p).is_in_edge_interior()

    for v in P.corners():
        p = v.vector() + vector((1 / 1024, 0))
        p_position = P.get_point_position(p)

        q = v.vector() + vector((1 / 2048, 0))
        q_position = P.get_point_position(q)

        assert str(p_position) == str(q_position)

    for v in P.corners():
        p = v.vector() - vector((1 / 1024, 0))
        p_position = P.get_point_position(p)

        q = v.vector() - vector((1 / 2048, 0))
        q_position = P.get_point_position(q)

        assert str(p_position) == str(q_position)
