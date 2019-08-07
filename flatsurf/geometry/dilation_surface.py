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

from flatsurf.geometry.half_dilation_surface import HalfDilationSurface

class DilationSurface(HalfDilationSurface):
    r"""
    Dilation surface.

    A dilation surface is a (G,X) structure on a surface for the group
    of positive dilatations `G = \RR_+` acting on the plane `X = \RR^2`.
    """
    pass
