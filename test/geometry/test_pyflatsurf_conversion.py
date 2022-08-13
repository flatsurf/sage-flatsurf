# -*- coding: utf-8 -*-
r"""
Discriminant loci in H(1,1)
"""
######################################################################
# This file is part of sage-flatsurf.
#
#       Copyright (C) 2020 Vincent Delecroix
#
# sage-flatsurf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# sage-flatsurf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
######################################################################

import sys
import pytest

pytest.importorskip("pyflatsurf")

from sage.all import QQ, AA
from flatsurf import translation_surfaces
from flatsurf.geometry.pyflatsurf_conversion import from_pyflatsurf, to_pyflatsurf


@pytest.fixture
def sage_flatsurf_surface_sample():

    surfaces = []

    from sage.all import SymmetricGroup
    G = SymmetricGroup(3)
    r = G('(1,2,3)')
    u = G('(1,2)')
    O = translation_surfaces.origami(r, u)
    surfaces.append(O)
    surfaces.append(AA(2).sqrt() * O)

    for n in [3, 5, 7]:
        surfaces.append(translation_surfaces.veech_double_n_gon(n))

    for g in [3, 4]:
        surfaces.append(translation_surfaces.arnoux_yoccoz(g))

    for n in [3, 17]:
        surfaces.append(translation_surfaces.ward(n))

    return surfaces

@pytest.fixture
def pyflatsurf_surface_sample():
    from pyflatsurf.vector import Vectors
    from pyflatsurf.factory import make_surface

    surfaces = []
    V2 = Vectors(QQ)
    vertices = [(1, -3, 2, -1, 6, -5), (-2, 4, -6, 5, -4, 3)]
    vectors = [V2((1, 1)), V2((-1, 0)), V2((0, -1)), V2((1, 1)), V2((-1, 0)), V2((0, -1))]
    surfaces.append(make_surface(vertices, [v.vector for v in vectors]))

    return surfaces


def test_origami1():
    from sage.all import SymmetricGroup
    G = SymmetricGroup(2)
    r = u = G('(1,2)')
    O = translation_surfaces.origami(r, u)
    S = to_pyflatsurf(O)
    assert str(S) == "FlatTriangulationCombinatorial(vertices = (1, 5, -4, 3, -5, -2)(-1, -6, 4, -3, 6, 2), faces = (1, 2, -5)(-1, -2, 6)(3, 4, 5)(-3, -4, -6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (1, 1), 6: (1, 1)}"


def test_origami2():
    from sage.all import SymmetricGroup
    G = SymmetricGroup(3)
    r = G('(1,2,3)')
    u = G('(1,2)')
    O = translation_surfaces.origami(r, u)
    S = to_pyflatsurf(O)
    assert str(S) == "FlatTriangulationCombinatorial(vertices = (1, 7, -4, -6, -9, 4, -3, 8, 2, -1, -8, -5, 6, 9, 5, 3, -7, -2), faces = (1, 2, -7)(-1, -2, 8)(3, 4, 7)(-3, 5, -8)(-4, -9, 6)(-5, 9, -6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (1, 1), 8: (1, 1), 9: (1, 1)}"


def test_to_pyflatsurf(sage_flatsurf_surface_sample):
    for S in sage_flatsurf_surface_sample:
        f = S._pyflatsurf
        T = f.codomain()

        assert f.domain() is S

        for e in S.edge_iterator():
            assert f.half_edge_from_pyflatsurf(f.half_edge_to_pyflatsurf(e)) == e, S

        for h in T.halfEdges():
            e = f.half_edge_from_pyflatsurf(h)
            assert e is None or f.half_edge_to_pyflatsurf(e) == h, S


def test_from_pyflatsurf(pyflatsurf_surface_sample):
    from flatsurf.geometry.pyflatsurf_conversion import FlatsurfConverter
    for T in pyflatsurf_surface_sample:
        f = FlatsurfConverter.from_pyflatsurf(T)
        S = f.domain()

        assert not S.is_mutable()
        assert f.codomain() == T
        
        for e in S.edge_iterator():
            assert f.half_edge_from_pyflatsurf(f.half_edge_to_pyflatsurf(e)) == e, T

        for h in T.halfEdges():
            e = f.half_edge_from_pyflatsurf(h)
            assert f.half_edge_to_pyflatsurf(e) == h, T
