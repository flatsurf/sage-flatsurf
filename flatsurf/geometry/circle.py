#TODO: Delete this file.
r"""
This class contains methods useful for working with circles.

This will be used to build a LazyDelaunayTriangulation class which will compute the
Delaunay decomposition for infinite surfaces.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
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

from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector


# TODO: This shoul be a method of EuclideanPlane.
def circle_from_three_points(p, q, r, base_ring=None):
    r"""
    Construct a circle from three points on the circle.
    """
    if base_ring is None:
        base_ring = p.base_ring()
    V2 = VectorSpace(base_ring.fraction_field(), 2)
    V3 = VectorSpace(base_ring.fraction_field(), 3)

    v1 = V3((p[0] + q[0], p[1] + q[1], 2))
    v2 = V3((p[1] - q[1], q[0] - p[0], 0))
    line1 = v1.cross_product(v2)
    v1 = V3((p[0] + r[0], p[1] + r[1], 2))
    v2 = V3((p[1] - r[1], r[0] - p[0], 0))
    line2 = v1.cross_product(v2)
    center_3 = line1.cross_product(line2)
    if center_3[2].is_zero():
        raise ValueError("The three points lie on a line.")
    center = V2((center_3[0] / center_3[2], center_3[1] / center_3[2]))
    return Circle(center, (p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)


