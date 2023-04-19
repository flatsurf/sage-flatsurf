# ****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#                          2023 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .polygon import wedge_product
from .half_dilation_surface import HalfDilationSurface

from sage.rings.all import QQ, AA
from sage.matrix.constructor import matrix


class HalfTranslationSurface(HalfDilationSurface):
    r"""
    A half translation surface has gluings between polygons whose monodromy is +I or -I.
    """

    def __init__(self, surface, category=None):
        from flatsurf.geometry.categories import HalfTranslationSurfaces
        super().__init__(surface, category or surface.category() & HalfTranslationSurfaces().Oriented())

