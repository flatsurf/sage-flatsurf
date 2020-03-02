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

import itertools

from .surface import Surface
from .half_dilation_surface import HalfDilationSurface
from .rational_cone_surface import RationalConeSurface

from sage.rings.all import QQ

class HalfTranslationSurface(HalfDilationSurface, RationalConeSurface):
    r"""
    A half translation surface has gluings between polygons whose monodromy is +I or -I.
    """
    def angles(self):
        r"""
        Return the set of angles around the vertices of the surface.

        EXAMPLES::

            sage: import flatsurf.geometry.similarity_surface_generators as sfg
            sage: sfg.translation_surfaces.regular_octagon().angles()
            [3]
            sage: sfg.translation_surfaces.veech_2n_gon(5).angles()
            [2, 2]
            sage: sfg.translation_surfaces.veech_2n_gon(6).angles()
            [5]
            sage: sfg.translation_surfaces.veech_double_n_gon(5).angles()
            [3]
            sage: sfg.translation_surfaces.cathedral(1, 1).angles()
            [3, 3, 3]
        """
        edges = set(self.edge_iterator())
        angles = []
        while edges:
            pair = p,e = next(iter(edges))
            ve = self.polygon(p).edge(e)
            angle = 0
            while pair in edges:
                edges.remove(pair)
                ppair = pp,ee = self.opposite_edge(p,(e-1)%self.polygon(p).num_edges())
                vee = self.polygon(pp).edge(ee)
                angle += (ve[0] > 0 and vee[0] <= 0) or (ve[0] < 0 and vee[0] >= 0) or (ve[0] == vee[0] == 0)
                pair, p, e, ve = ppair, pp, ee, vee
            angles.append(QQ((angle,2)))
        return angles

    def _test_edge_matrix(self, **options):
        r"""
        Check the compatibility condition
        """
        tester = self._tester(**options)
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if self.is_finite():
            it = self.label_iterator()
        else:
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                # Warning: check the matrices computed from the edges,
                # rather the ones overriden by TranslationSurface.
                m = SimilaritySurface.edge_matrix(self,lab,e)
                tester.assertTrue(m.is_one() or (-m).is_one(),
                    "edge_matrix between edge e={} and e'={} has matrix\n{}\nwhich is neither a translation nor a rotation by pi".format((lab,e), self.opposite_edge((lab,e)), m))

