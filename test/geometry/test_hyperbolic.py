r"""
Test functionality of the HyperbolicPlane.
"""
# ****************************************************************************
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
# ****************************************************************************

import pytest


def random_half_spaces(H, count, predicate=lambda half_space: True):
    r"""
    Return ``count`` random half spaces that satisfy ``predicate``.
    """
    half_spaces = []

    while len(half_spaces) < count:
        half_space = H.random_element("half_space")

        if predicate(half_space):
            half_spaces.append(half_space)

    return half_spaces


@pytest.mark.repeat(1024)
def test_intersection_point():
    r"""
    Intersect random half spaces that contain a certain point and do not
    contain some other point.
    """
    from flatsurf.geometry.hyperbolic import HyperbolicPlane

    H = HyperbolicPlane()

    inside = H.random_element("point")

    while True:
        outside = H.random_element("point")
        if inside != outside:
            break

    from sage.all import randint

    half_spaces = random_half_spaces(
        H,
        randint(1, 8),
        lambda half_space: inside in half_space and outside not in half_space,
    )

    intersection = H.polygon(half_spaces)

    # By construction, inside must be in the polygon and outside must not be.
    assert (
        inside in intersection
    ), f"{inside} not in {intersection} which was found to be the intersection of {half_spaces}"
    assert (
        outside not in intersection
    ), f"{outside} in {intersection} which was found to be the intersection of {half_spaces}"

    # Check that polygon is idempotent.
    assert H.polygon(intersection.half_spaces()) == intersection

    # Check that the vertices of the intersection are correct.
    for i, e in enumerate(half_spaces):
        for f in half_spaces[i + 1 :]:
            p = e.boundary().intersection(f.boundary())
            if p.is_point():
                assert all(p in half_space for half_space in half_spaces) == (
                    p in intersection
                )
