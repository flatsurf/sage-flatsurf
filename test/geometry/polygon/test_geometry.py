r"""
Test basic geometry methods used in polygon construction.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2020-2023 Julian Rüth
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

from sage.all import QQ, randint


@pytest.mark.repeat(1024)
def test_is_parallel():
    from flatsurf.geometry.euclidean import is_parallel

    V = QQ**2

    while True:
        v = V.random_element()
        if v:
            break
    assert is_parallel(v, 2 * v)
    assert not is_parallel(v, -v)


@pytest.mark.repeat(100)
def test_is_anti_parallel():
    from flatsurf.geometry.euclidean import is_anti_parallel

    V = QQ**2

    while True:
        v = V.random_element()
        if v:
            break

    assert not is_anti_parallel(v, v)
    assert not is_anti_parallel(v, 2 * v)
    assert is_anti_parallel(v, -v)


@pytest.mark.repeat(4096)
def test_segment_intersect():
    from flatsurf.geometry.euclidean import is_segment_intersecting

    while True:
        us = (randint(-4, 4), randint(-4, 4))
        ut = (randint(-4, 4), randint(-4, 4))
        vs = (randint(-4, 4), randint(-4, 4))
        vt = (randint(-4, 4), randint(-4, 4))
        if us != ut and vs != vt:
            break

    ans1 = is_segment_intersecting((us, ut), (vs, vt))
    ans2 = is_segment_intersecting((ut, us), (vs, vt))
    ans3 = is_segment_intersecting((us, ut), (vt, vs))
    ans4 = is_segment_intersecting((ut, us), (vt, vs))
    ans5 = is_segment_intersecting((vs, vt), (us, ut))
    ans6 = is_segment_intersecting((vt, vs), (us, ut))
    ans7 = is_segment_intersecting((vs, vt), (ut, us))
    ans8 = is_segment_intersecting((vt, vs), (ut, us))
    assert ans1 == ans2 == ans3 == ans4 == ans5 == ans6 == ans7 == ans8, (
        us,
        ut,
        vs,
        vt,
        ans1,
        ans2,
        ans3,
        ans4,
        ans5,
        ans6,
        ans7,
        ans8,
    )


def test_is_between():
    from flatsurf.geometry.euclidean import is_between

    V = QQ**2

    vecs = [
        V((1, 0)),
        V((2, 1)),
        V((1, 1)),
        V((1, 2)),
        V((0, 1)),
        V((-1, 2)),
        V((-1, 1)),
        V((-2, 1)),
        V((-1, 0)),
        V((-2, -1)),
        V((-1, -1)),
        V((-1, -2)),
        V((0, -1)),
        V((1, -2)),
        V((1, -1)),
        V((2, -1)),
    ]
    for i, a in enumerate(vecs):
        for j, b in enumerate(vecs):
            for k, c in enumerate(vecs):
                if a == b or a == c or b == c:
                    continue
                assert is_between(a, b, c) == (i < k < j or k < j < i or j < i < k)
