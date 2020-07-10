#!/usr/bin/env python3
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

pytest.importorskip('pyflatsurf')

from sage.all import polygen, NumberField, AA, QQ
from flatsurf import translation_surfaces, GL2ROrbitClosure

# TODO: Does not work because of sage-flatsurf
#def test_D9():
#    x = polygen(QQ)
#    K = NumberField(x**3 - 2, 'a', embedding=AA(2)**QQ((1,3)))
#    a = K.gen()
#    S = translation_surfaces.mcmullen_genus2_prototype(2,1,0,-1,a/4)
#    O = GL2ROrbitClosure(S)
#    for d in O.decompositions(5):
#        ncyl, nmin, nund = d.num_cylinders_minimals_undetermined()
#        assert (nund == 0)
#        assert ((nmin == 0) or (ncyl == 0 and 1 <= nmin <= 2))
#        O.update_tangent_space_from_flow_decomposition(d)
#    assert O.U.dimension() == 3

def test_D33():
    S = translation_surfaces.mcmullen_genus2_prototype(4,2,1,1,QQ((1,4)))
    O = GL2ROrbitClosure(S)
    for d in O.decompositions(5, 100):
        ncyl, nmin, nund = d.num_cylinders_minimals_undetermined()
        assert (nund == 0)
        assert ((nmin == 0) or (ncyl == 0 and nmin == 2))
        O.update_tangent_space_from_flow_decomposition(d)
    assert O.dimension() == 3
