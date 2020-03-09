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

from .similarity_surface import SimilaritySurface

class ConeSurface(SimilaritySurface):
    r"""
    A Euclidean cone surface.
    """
    def angles(self, numerical=False, return_adjacent_edges=False):
        r"""
        Return the set of angles around the vertices of the surface.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: T = polygons.triangle(3, 4, 5)
            sage: S = similarity_surfaces.billiard(T)
            sage: S.angles()
            [1/3, 1/4, 5/12]
            sage: S.angles(numerical=True)   # abs tol 1e-14
            [0.333333333333333, 0.250000000000000, 0.416666666666667]

            sage: S.angles(return_adjacent_edges=True)
            [(1/3, [(0, 1), (1, 2)]), (1/4, [(0, 0), (1, 0)]), (5/12, [(1, 1), (0, 2)])]
        """
        if not self.is_finite():
            raise NotImplementedError("the set of edges is infinite!")

        edges = [pair for pair in self.edge_iterator()]
        edges = set(edges)
        angles = []

        if return_adjacent_edges:
            while edges:
                p,e = edges.pop()
                adjacent_edges = [(p,e)]
                angle = self.polygon(p).angle(e, numerical=numerical)
                pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
                while pp != p or ee != e:
                    edges.remove((pp,ee))
                    adjacent_edges.append((pp,ee))
                    angle += self.polygon(pp).angle(ee, numerical=numerical)
                    pp,ee = self.opposite_edge(pp,(ee-1)%self.polygon(pp).num_edges())
                angles.append((angle, adjacent_edges))
        else:
            while edges:
                p,e = edges.pop()
                angle = self.polygon(p).angle(e, numerical=numerical)
                pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
                while pp != p or ee != e:
                    edges.remove((pp,ee))
                    angle += self.polygon(pp).angle(ee, numerical=numerical)
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
