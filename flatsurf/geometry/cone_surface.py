#*****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from flatsurf.geometry.similarity_surface import SimilaritySurface

class ConeSurface(SimilaritySurface):
    r"""
    A Euclidean cone surface.
    """

    def angles(self):
        r"""
        Return the set of angles around the vertices of the surface.

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.regular_octagon().angles()
            [3]
        """
        if not self.is_finite():
            raise NotImplementedError("the set of edges is infinite!")

        edges = [pair for pair in self.edge_iterator()]
        edges = set(edges)
        angles = []
        while edges:
            p,e = edges.pop()
            angle = self.polygon(p).angle(e)
            pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
            while pp != p or ee != e:
                edges.remove((pp,ee))
                angle += self.polygon(pp).angle(ee)
                pp,ee = self.opposite_edge(pp,(ee-1)%self.polygon(pp).num_edges())
            angles.append(angle)
        return angles

    def genus(self):
        """
        Return the genus of the surface.

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.octagon_and_squares().genus()
            3

            sage: from flatsurf import *
            sage: T = polygons.triangle(3,4,5)
            sage: B = RationalConeSurface(similarity_surfaces.billiard(T))
            sage: B.genus()
            0
            sage: B.minimal_cover("translation").genus()
            3
        """
        return sum(a - 1 for a in self.angles()) // 2 + 1

    def area(self):
        r"""
        Return the area of this surface.
        """
        return self._s.area()

