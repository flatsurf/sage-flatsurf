# -*- coding: utf-8 -*-
r"""
Interface with pyflatsurf
"""
#*********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2019      Vincent Delecroix
#                      2019-2020 Julian Rüth
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
#*********************************************************************

from sage.all import QQ

from pyflatsurf import flatsurf
from pyflatsurf.factory import make_surface
from pyflatsurf.vector import Vectors

from .translation_surface import TranslationSurface

def _check_data(vp, fp, vec):
    r"""
    Check consistency of data

    vp - vector permutation
    fp - face permutation
    vec - vectors of the flat structure
    """
    assert isinstance(vp, list)
    assert isinstance(fp, list)
    assert isinstance(vec, list)

    n = len(vp) - 1

    assert n%2 == 0, n
    assert len(fp) == n+1
    assert len(vec) == n//2

    assert vp[0] is None
    assert fp[0] is None

    for i in range(1, n//2+1):
        # check fp/vp consistency
        assert fp[-vp[i]] == i, i

        # check that each face is a triangle and that vectors sum up to zero
        j = fp[i]
        k = fp[j]
        assert i != j and i != k and fp[k] == i, (i, j, k)
        vi = vec[i-1] if i >= 1 else -vec[-i-1]
        vj = vec[j-1] if j >= 1 else -vec[-j-1]
        vk = vec[k-1] if k >= 1 else -vec[-k-1]
        v = vi + vj + vk
        assert v.x() == 0, v.x()
        assert v.y() == 0, v.y()

def _cycle_decomposition(p):
    n = len(p) - 1
    assert n%2 == 0
    cycles = []
    unseen = [True] * (n+1)
    for i in list(range(-n//2+1, 0)) + list(range(1, n//2)):
        if unseen[i]:
            j = i
            cycle = []
            while unseen[j]:
                unseen[j] = False
                cycle.append(j)
                j = p[j]
            cycles.append(cycle)
    return cycles

def to_pyflatsurf(S):
    r"""
    Given S a translation surface from sage-flatsurf builds a
    flatsurf::FlatTriangulation from libflatsurf/pyflatsurf.
    """
    if not isinstance(S, TranslationSurface):
        raise TypeError("S must be a translation surface")
    if not S.is_finite():
        raise ValueError("the surface S must be finite")

    S = S.triangulate()

    # populate half edges and vectors
    n = sum(S.polygon(lab).num_edges() for lab in S.label_iterator())
    half_edge_labels = {}   # map: (face lab, edge num) in faltsurf -> integer
    vec = []                # vectors
    k = 1                   # half edge label in {1, ..., n}
    for t0, t1 in S.edge_iterator(gluings=True):
        if t0 in half_edge_labels:
            continue

        half_edge_labels[t0] = k
        half_edge_labels[t1] = -k

        f0, e0 = t0
        p = S.polygon(f0)
        vec.append(p.edge(e0))

        k += 1

    # compute vertex and face permutations
    vp = [None] * (n+1)  # vertex permutation
    fp = [None] * (n+1)  # face permutation
    for t in S.edge_iterator(gluings=False):
        e = half_edge_labels[t]
        j = (t[1] + 1) % S.polygon(t[0]).num_edges()
        fp[e] = half_edge_labels[(t[0], j)]
        vp[fp[e]] = -e


    # convert the vp permutation into cycles
    verts = _cycle_decomposition(vp)
    V = Vectors(S.base_ring())
    vec = [V(v).vector for v in vec]
    _check_data(vp, fp, vec)

    return make_surface(verts, vec)
