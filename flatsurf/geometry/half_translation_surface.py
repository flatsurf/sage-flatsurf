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

from .surface import Surface
from .half_dilation_surface import HalfDilationSurface
from .rational_cone_surface import RationalConeSurface

class HalfTranslationSurface(HalfDilationSurface, RationalConeSurface):
    r"""
    A half translation surface has gluings between polygons whose monodromy is +I or -I.
    """
    
    def _test_edge_matrix(self, **options):
        r"""
        Check the compatibility condition
        """
        tester = self._tester(**options)
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                # Warning: check the matrices computed from the edges,
                # rather the ones overriden by TranslationSurface.
                tester.assertTrue(SimilaritySurface.edge_matrix(self,lab,e).is_one() or \
                    (-SimilaritySurface.edge_matrix(self,lab,e)).is_one(), \
                    "edge_matrix of edge "+str((lab,e))+" is not a translation or rotation by pi.")

