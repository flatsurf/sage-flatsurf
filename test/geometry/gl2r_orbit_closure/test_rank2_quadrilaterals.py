r"""
Eskin-McMullen-Mukamel-Write quadrilaterals from the article
"Billiards, quadrilaterals and moduli spaces"

Note that in some cases, the billiards might give a Veech surface and not
the ambient rank 2 locus.
"""
# ****************************************************************************
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
# ****************************************************************************

import pytest

import itertools

pytest.importorskip("pyflatsurf")  # noqa

from sage.all import AA, QQ
from flatsurf import EuclideanPolygonsWithAngles, similarity_surfaces, GL2ROrbitClosure


# TODO: the test for field of definition with is_isomorphic() does not check
# for embeddings... though for quadratic fields it does not matter much.
@pytest.mark.parametrize(
    "a,b,c,d,l1,l2,veech,discriminant",
    [
        (1, 1, 1, 7, 1, 1, True, 5),
        (1, 1, 1, 7, AA(2).sqrt(), 1, False, 5),
        (1, 1, 1, 9, 1, 1, True, 3),
        (1, 1, 1, 9, AA(2).sqrt(), 1, False, 1),
        (2, 2, 3, 13, 1, 1, False, 5),
        (1, 1, 2, 8, 1, 1, False, 1),
        (1, 1, 2, 12, 1, 1, False, 2),
        (1, 2, 2, 11, 1, 1, False, 2),
        (1, 2, 2, 15, 1, 1, False, 5),
    ],
)
def test_rank2_quadrilateral(a, b, c, d, l1, l2, veech, discriminant):
    E = EuclideanPolygonsWithAngles(a, b, c, d)
    P = E([l1, l2], normalized=True)
    B = similarity_surfaces.billiard(P, rational=True)
    S = B.minimal_cover(cover_type="translation")
    S = S.erase_marked_points()
    S, _ = S.normalized_coordinates()
    orbit_closure = GL2ROrbitClosure(S)
    assert orbit_closure.ambient_stratum() == E.billiard_unfolding_stratum(
        cover_type="translation"
    )
    D = itertools.islice(orbit_closure.decompositions(9, 40), 50)
    if veech:
        assert S.base_ring().degree() <= 2
        for dec in D:
            orbit_closure.update_tangent_space_from_flow_decomposition(dec)
            assert dec.parabolic()
        assert orbit_closure.dimension() == orbit_closure.absolute_dimension() == 2, (
            orbit_closure.dimension(),
            orbit_closure.absolute_dimension(),
        )
    else:
        for dec in D:
            orbit_closure.update_tangent_space_from_flow_decomposition(dec)
        assert orbit_closure.absolute_dimension() == orbit_closure.dimension() == 4, (
            orbit_closure.dimension(),
            orbit_closure.absolute_dimension(),
        )

    if discriminant == 1:
        assert orbit_closure.field_of_definition() == QQ
    else:
        assert orbit_closure.field_of_definition().degree() == 2
        assert (
            orbit_closure.field_of_definition().discriminant().squarefree_part()
            == discriminant
        )
