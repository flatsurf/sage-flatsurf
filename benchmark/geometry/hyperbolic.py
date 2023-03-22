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
from flatsurf.geometry.hyperbolic import HyperbolicPlane


H = HyperbolicPlane()


def time_HyperbolicPlane_polygon(half_spaces):
    r"""
    Time how quickly we can form a convex polygon from a set of half spaces.

    The bottleneck here should be sorting the half spaces by angle, i.e., this
    should run in O(nlogn).
    """
    H.polygon(half_spaces)


class SilentSet(set):
    r"""
    A set that prints more compactly, i.e., without listing all its elements so
    we do not break asv's console output.
    """

    def __repr__(self):
        if len(self) <= 1:
            return repr(set(self))
        return f"[{next(iter(self))}, … {len(self) - 1} more]"


time_HyperbolicPlane_polygon.params = [
    SilentSet([H.half_circle(0, 1).right_half_space()]),
    SilentSet([H.half_circle(i, 1).right_half_space() for i in range(32)]),
    SilentSet([H.half_circle(i, 1).right_half_space() for i in range(1024)]),
    SilentSet([H.half_circle(i, 1).right_half_space() for i in range(2048)]),
    SilentSet([H.half_circle(i, 1).right_half_space() for i in range(4096)]),
]
