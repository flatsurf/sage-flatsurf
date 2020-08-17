# -*- coding: utf-8 -*-
r"""
The gothic locus

From the article McMullen-Mukamel-Wright (2017).
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

pytest.importorskip('pyflatsurf')

from sage.all import QQ, AA, NumberField, polygen
from flatsurf import translation_surfaces, GL2ROrbitClosure
from surface_dynamics import AbelianStratum

def test_gothic_generic():
    x = polygen(QQ)
    K = NumberField(x**3 - 2, 'a', embedding=AA(2)**QQ((1,3)))
    a = K.gen()
    S = translation_surfaces.cathedral(a, a**2)
    O = GL2ROrbitClosure(S)
    assert O.ambient_stratum() == AbelianStratum(2, 2, 2)
    for d in O.decompositions(4, 50):
        O.update_tangent_space_from_flow_decomposition(d)
    assert O.dimension() == O.absolute_dimension() == 4
    assert O.field_of_definition() == QQ

def test_gothic_veech():
    x = polygen(QQ)
    K = NumberField(x**2 - 2, 'sqrt2', embedding=AA(2)**QQ((1,2)))
    sqrt2 = K.gen()
    x = QQ((1,2))
    y = 1
    a = x + y * sqrt2
    b = -3*x -QQ((3,2)) + 3*y*sqrt2
    S = translation_surfaces.cathedral(a,b)
    O = GL2ROrbitClosure(S)
    assert O.ambient_stratum() == AbelianStratum(2, 2, 2)
    for d in O.decompositions(4, 50):
        assert d.parabolic()
#        assert d.decomposition.cylinder_diagram()[0].stratum() == O.surface.stratum()
        O.update_tangent_space_from_flow_decomposition(d)
    assert O.dimension() == O.absolute_dimension() == 2
    assert O.field_of_definition() == O.base_ring()
