r"""
Similarity surfaces.
"""
# *********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
# *********************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.sage_unittest import TestSuite

from sage.structure.parent import Parent

from sage.rings.infinity import Infinity

from sage.rings.all import ZZ, QQ, AA, NumberField

from sage.modules.free_module_element import vector

from .similarity import SimilarityGroup
from .polygon import ConvexPolygons, wedge_product

from .surface import Surface, Surface_dict, Surface_list, LabelComparator
from .surface_objects import Singularity, SaddleConnection, SurfacePoint
from .circle import Circle


