# -*- coding: utf-8 -*-
r"""
Thurston-Veech constructions
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

from sage.all import QQ

from flatsurf import GL2ROrbitClosure
from flatsurf.geometry.thurston_veech import ThurstonVeech

def test_H2():
    # Calta-Mcmullen theorem: TV constructions in H(2) are Veech surfaces
    TV = ThurstonVeech('(1,2)', '(1,3)')
    
    for hm, vm, bound, discriminant in [((1,1),(1,1),5,5),
                    ((1,3),(2,5),5,1), ((1,1),(2,1),5,17)]:
        S = TV(hm, vm)
        O = GL2ROrbitClosure(S)
        dh = O.decomposition((1,0))
        assert dh.parabolic() is True
        dv = O.decomposition((0,1))
        assert dv.parabolic() is True
        if discriminant == 1:
            assert O.field_of_definition() is QQ
        else:
            assert O.field_of_definition().discriminant().squarefree_part() == discriminant
        assert O.is_teichmueller_curve(bound) is not False

def test_H11():
    # McMullen theorem: TV constructions in H(1,1) are either Veech (discriminant 5) or
    # eigenform loci
    TV = ThurstonVeech('(1,2)(3,4)', '(2,3)')

    for hm, vm, bound, discriminant in [((1, 1), (1, 1, 1), 5, 1),
            ((1,2), (1,1,1), 5, 3), ((1,1),(2,1,1), 5, 5),
            ((1,1),(1,2,1),5, 1), ((1,2),(1,2,1), 5, 41)]:
        S = TV(hm, vm)
        O = GL2ROrbitClosure(S)
        dh = O.decomposition((1,0))
        assert dh.parabolic() is True
        dv = O.decomposition((0,1))
        assert dv.parabolic() is True

        if discriminant == 1:
            assert O.field_of_definition() is QQ
        else:
            assert O.field_of_definition().degree() == 2
            assert O.field_of_definition().discriminant().squarefree_part() == discriminant

        for d in O.decompositions(bound):
            O.update_tangent_space_from_flow_decomposition(d)

        assert O.dimension() == 2 or O.dimension() == 3
