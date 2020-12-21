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

from sage.all import QQ, AA

from pyflatsurf import flatsurf
from pyflatsurf.factory import make_surface
from pyflatsurf.vector import Vectors

from .translation_surface import TranslationSurface
from .surface import Surface_list
from .polygon import ConvexPolygons

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
    Given S a translation surface from sage-flatsurf return a
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

    # find a finite SageMath base ring that contains all the coordinates
    base_ring = S.base_ring()
    if base_ring is AA:
        from sage.rings.qqbar import number_field_elements_from_algebraics
        from itertools import chain
        base_ring = number_field_elements_from_algebraics(list(chain(*[list(v) for v in vec])), embedded=True)[0]

    V = Vectors(base_ring)
    vec = [V(v).vector for v in vec]

    _check_data(vp, fp, vec)

    return make_surface(verts, vec)

def sage_base_ring(T):
    r"""
    Given T a flatsurf::FlatTriangulation from libflatsurf/pyflatsurf, return a
    SageMath base ring over which from_pyflatsurf() can construct a translation
    surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, sage_base_ring # optional: pyflatsurf
        sage: S = to_pyflatsurf(translation_surfaces.veech_double_n_gon(5)) # optional: pyflatsurf
        sage: ring, to_ring = sage_base_ring(S) # optional: pyflatsurf
        sage: ring # optional: pyflatsurf
        Number Field in a with defining polynomial x^4 - 5*x^2 + 5 with a = 1.902113032590308?

    """
    import cppyy

    def maybe_type(t):
        try:
            return t()
        except AttributeError:
            # The type constructed by t might not exist because the required C++ library has not been loaded.
            return None

    coordinate = type(T.fromHalfEdge(1).x())
    if coordinate is maybe_type(lambda: cppyy.gbl.mpz_class):
        to_ring = ring = ZZ
    elif coordinate is maybe_type(lambda: cppyy.gbl.mpq_class):
        to_ring = ring = QQ
    elif coordinate is maybe_type(lambda: cppyy.gbl.eantic.renf_elem_class):
        from pyeantic import RealEmbeddedNumberField
        rng = RealEmbeddedNumberField(T.fromHalfEdge(1).x().parent())
        ring = rng.number_field
        to_ring = lambda x: ring(rng(x))
    elif coordinate is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.IntegerRing]):
        from pyexactreal import ExactReals
        to_ring = ring = ExactReals(ZZ)
    elif coordinate is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.RationalField]):
        from pyexactreal import ExactReals
        to_ring = ring = ExactReals(QQ)
    elif coordinate is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.NumberField]):
        from pyexactreal import ExactReals
        to_ring = ring = ExactReals(T.fromHalfEdge(1).x().module().ring().parameters)
    else:
        raise NotImplementedError("unknown coordinate ring %s"%(coordinate,))

    return ring, to_ring

def from_pyflatsurf(T):
    r"""
    Given T a flatsurf::FlatTriangulation from libflatsurf/pyflatsurf, return a
    sage-flatsurf translation surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, from_pyflatsurf # optional: pyflatsurf
        sage: S = translation_surfaces.veech_double_n_gon(5) # optional: pyflatsurf
        sage: from_pyflatsurf(to_pyflatsurf(S)) # optional: pyflatsurf
        TranslationSurface built from 6 polygons

    """
    import cppyy

    ring, to_ring = sage_base_ring(T)

    S = Surface_list(ring)
    P = ConvexPolygons(ring)
    V = P.module()

    half_edges = {}

    for face in T.faces():
        a, b, c = map(cppyy.gbl.flatsurf.HalfEdge, face)

        vectors = [T.fromHalfEdge(he) for he in face]
        vectors = [V([to_ring(v.x()), to_ring(v.y())]) for v in vectors]
        triangle = P(vectors)
        face_id = S.add_polygon(triangle)

        assert(a not in half_edges)
        half_edges[a] = (face_id, 0)
        assert(b not in half_edges)
        half_edges[b] = (face_id, 1)
        assert(c not in half_edges)
        half_edges[c] = (face_id, 2)

    for half_edge, (face, id) in half_edges.items():
        _face, _id = half_edges[-half_edge]
        S.change_edge_gluing(face, id, _face, _id)

    S.set_immutable()

    return TranslationSurface(S)
