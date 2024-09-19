r"""
Do not use this module, use :mod:`flatsurf.geometry.pyflatsurf.conversion` instead.

EXAMPLES::

    sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
    doctest:warning
    ...
    UserWarning: the flatsurf.geometry.pyflatsurf_conversion module has been deprecated and will be removed in a future version of sage-flatsurf; use flatsurf.geometry.pyflatsurf.conversion instead.
    Note that you might also want to call .pyflatsurf().codomain().flat_triangulation() on your surface now to convert a surface to the pyflatsurf FlatTriangulation

"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2024 Julian Rüth
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
# ********************************************************************
from flatsurf.geometry.pyflatsurf.conversion import *

import warnings

warnings.warn(
    "the flatsurf.geometry.pyflatsurf_conversion module has been deprecated and will be removed in a future version of sage-flatsurf; use flatsurf.geometry.pyflatsurf.conversion instead. "
    "Note that you might also want to call .pyflatsurf().codomain().flat_triangulation() on your surface now to convert a surface to the pyflatsurf FlatTriangulation"
)
