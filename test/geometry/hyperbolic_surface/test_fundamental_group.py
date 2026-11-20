r"""
Generators, vertex loops, free reduction and equality of the fundamental
groups of surfaces built from polygons.
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

from sage.all import QQ

from flatsurf import (
    MutableOrientedSimilaritySurface,
    Polygon,
    similarity_surfaces,
    translation_surfaces,
)


def unit_square_surface(gluings):
    r"""
    Return the unit square with the given pairs of edges glued.
    """
    S = MutableOrientedSimilaritySurface(QQ)
    S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]), label=0)
    for e, f in gluings:
        S.glue((0, e), (0, f))
    S.set_immutable()
    return S


def test_vertex_loops_square_torus():
    S = translation_surfaces.square_torus()
    F = S.fundamental_group(punctured_vertices=True)
    a, b = F.gens()
    loops = F.vertex_loops()
    assert list(loops) == [S(0, 0)]
    assert loops[S(0, 0)] == ~a * ~b * a * b


def test_vertex_loops_square_glued_by_rotations():
    S = unit_square_surface([(0, 1), (2, 3)])
    F = S.fundamental_group(punctured_vertices=True)
    x, y = F.gens()
    loops = F.vertex_loops()
    assert len(loops) == 3
    assert loops[S(0, 1)] == x
    assert loops[S(0, 3)] == y
    assert loops[S(0, 0)] == ~x * ~y
    assert (loops[S(0, 3)] * loops[S(0, 1)] * loops[S(0, 0)]).is_one()


@pytest.mark.parametrize(
    "surface, ngens",
    [
        (translation_surfaces.regular_octagon, 4),
        (translation_surfaces.octagon_and_squares, 6),
    ],
)
def test_number_of_generators(surface, ngens):
    F = surface().fundamental_group(punctured_vertices=True)
    assert len(F.gens()) == ngens


def test_cylinder():
    S = unit_square_surface([(1, 3)])
    F = S.fundamental_group(punctured_vertices=True)
    assert len(F.gens()) == 1
    assert F.vertex_loops() == {}


def test_free_reduction():
    S = translation_surfaces.square_torus()
    F = S.fundamental_group(punctured_vertices=True)
    assert F([0, 2]).is_one()


def test_self_glued_involution():
    P = Polygon(vertices=[(0, 0), (2, 0), (1, 4), (0, 5)])
    S = similarity_surfaces.self_glued_polygon(P)
    F = S.fundamental_group(punctured_vertices=True)
    h = F.gens()[0]
    assert not h.is_one()
    assert (h * h).is_one()


def test_default_group_equality():
    S = translation_surfaces.square_torus()
    F = S.fundamental_group()
    a, b = F.gens()
    assert a * b == a * b
    with pytest.raises(NotImplementedError):
        a * b == b * a
