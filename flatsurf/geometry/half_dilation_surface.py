# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
#                          2023 Julian Rüth
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

# Move to lazy.py?
from flatsurf.geometry.surface import OrientedSimilaritySurface

from sage.misc.cachefunc import cached_method


class GL2RImageSurface(OrientedSimilaritySurface):
    @cached_method
    def polygon(self, lab):

    def opposite_edge(self, p, e):
        if self._det_sign == 1:
            return self._s.opposite_edge(p, e)
        else:
            polygon = self._s.polygon(p)
            pp, ee = self._s.opposite_edge(p, len(polygon.vertices()) - 1 - e)
            polygon2 = self._s.polygon(pp)
            return pp, len(polygon2.vertices()) - 1 - ee
