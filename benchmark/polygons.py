# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
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
from flatsurf import polygons


def time_triangle(angles):
    r"""
    Time how quickly we can form a triangle from a given triple of angles.

    The runtime of this should be mostly bound by the time it takes to
    construct the underlying number field.
    """
    polygons.triangle(*angles)


time_triangle.params = [
    # We pass the parameters as sets since otherwise asv tries to be smart and
    # runs with the columns of the following matrix instead of the rows.
    {3, 4, 5},
    {22, 23, 24},
    {23, 24, 25},
    {24, 25, 26},
    {26, 48, 75},
]
