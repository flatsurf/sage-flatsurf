# -*- coding: utf-8 -*-
r"""
Interface with pyflatsurf
"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2019      Vincent Delecroix
#                      2019-2022 Julian Rüth
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

from sage.all import QQ, AA, ZZ

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

    assert n % 2 == 0, n
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
    assert n % 2 == 0
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


class FlatsurfConverter:
    r"""
    Conversion between sage-flatsurf and libflatsurf/pyflatsurf.
    """
    def __init__(self, S):
        from flatsurf.features import pyflatsurf_feature
        pyflatsurf_feature.require()
        from flatsurf.geometry.similarity_surface import SimilaritySurface

        if isinstance(S, SimilaritySurface):
            self._init_from_sage_flatsurf(S)
        else:
            self._init_from_pyflatsurf(S)

    def _init_from_sage_flatsurf(self, S):
        from flatsurf.geometry.translation_surface import TranslationSurface
        if not isinstance(S, TranslationSurface):
            raise TypeError("S must be a translation surface")
        if not S.is_finite():
            raise ValueError("the surface S must be finite")

        # number of half-edges in the triangulation
        n = 3 * sum(S.polygon(label).num_edges() - 2 for label in S.label_iterator())

        self._half_edge_map = {} # (face label, edge number) -> id
        self._half_edge_section = [] # index -> (face label, edge number)
        self._triangulation_edges = [] # (face label, vertex start, vertex end)

        # populate half edges in the surface
        k = 1
        vectors = []
        for t0, t1 in S.edge_iterator(gluings=True):
            if t0 in self._half_edge_map:
                continue

            self._half_edge_map[t0] = k
            self._half_edge_map[t1] = -k
            self._half_edge_section.append(t0)
            self._half_edge_section.append(t1)
            f0, e0 = t0
            p = S.polygon(f0)
            vectors.append(p.edge(e0))
            k += 1

        # possibly triangulate and compute face permutation
        from .polygon import triangulate, build_faces
        fp = [None] * (n + 1)
        for label in S.label_iterator():
            poly = S.polygon(label)
            ne = poly.num_edges()
            if ne == 3:
                # already a triangle
                for i in range(3):
                    fp[self._half_edge_map[(label, i)]] = self._half_edge_map[(label, (i+1) % 3)]
            else:
                local_edges = {}
                for i in range(ne):
                    local_edges[(i, (i+1)%ne)] = self._half_edge_map[(label, i)]
                new_edges = triangulate(poly.vertices())
                assert len(new_edges) == ne - 3
                for a, b in new_edges:
                    local_edges[(a, b)] = k
                    local_edges[(b, a)] = -k
                    self._half_edge_section.append(None)
                    self._half_edge_section.append(None)
                    k += 1
                    self._triangulation_edges.append((label, a, b))
                    vectors.append(poly.vertex(b) - poly.vertex(a))
                for a, b, c in build_faces(ne, new_edges):
                    fp[local_edges[(a, b)]] = local_edges[(b, c)]
                    fp[local_edges[(b, c)]] = local_edges[(c, a)]
                    fp[local_edges[(c, a)]] = local_edges[(a, b)]

        assert k == 1 + n // 2, (k, n)
        assert all(fp[i] is not None and fp[-i] is not None for i in range(1, n//2 + 1))

        # compute vertex permutation from face permutation
        vp = [None] * (n + 1)
        for i in range(1, n // 2 + 1):
            vp[fp[-i]] = i
            vp[fp[i]] = -i

        # find a finite SageMath base ring that contains all the coordinates
        base_ring = S.base_ring()
        if base_ring is AA:
            from sage.rings.qqbar import number_field_elements_from_algebraics
            from itertools import chain
            base_ring = number_field_elements_from_algebraics(list(chain(*[list(v) for v in vectors])), embedded=True)[0]

        # build flatsurf data
        import pyflatsurf
        from pyflatsurf.vector import Vectors
        from pyflatsurf.factory import make_surface

        # A model of the vector space R² in libflatsurf, e.g., to represent the
        # vector associated to a saddle connection.
        self._V2 = Vectors(base_ring)
        vectors = [self._V2(v).vector for v in vectors]
        _check_data(vp, fp, vectors)
        self._pyflatsurf_surface = make_surface(_cycle_decomposition(vp), vectors)
        self._surface = S

    def _init_from_pyflatsurf(self, T):
        n = len(T.halfEdges())

        self._half_edge_map = {} # (face label, edge number) -> id
        self._half_edge_section = [None] * n # index -> (face label, edge number)
        self._triangulation_edges = [] # (face label, vertex start, vertex end)

        ring = sage_ring(T)

        from flatsurf.geometry.surface import Surface_list
        S = Surface_list(ring)

        from flatsurf.geometry.polygon import ConvexPolygons
        P = ConvexPolygons(ring)

        V = P.module()

        import pyflatsurf
        for face in T.faces():
            a, b, c = map(pyflatsurf.flatsurf.HalfEdge, face)

            vectors = [T.fromHalfEdge(he) for he in face]
            vectors = [V([ring(to_sage_ring(v.x())), ring(to_sage_ring(v.y()))]) for v in vectors]
            triangle = P(vectors)
            face_id = S.add_polygon(triangle)

            ea = (face_id, 0)
            eb = (face_id, 1)
            ec = (face_id, 2)
            self._half_edge_map[ea] = a.id()
            self._half_edge_map[eb] = b.id()
            self._half_edge_map[ec] = c.id()

            assert self._half_edge_section[a.index()] is None
            self._half_edge_section[a.index()] = ea
            assert self._half_edge_section[b.index()] is None
            self._half_edge_section[b.index()] = eb
            assert self._half_edge_section[c.index()] is None
            self._half_edge_section[c.index()] = ec

        assert all(e is not None for e in self._half_edge_section[1:]), self._half_edge_section

        for (face, id), k in self._half_edge_map.items():
            h = pyflatsurf.flatsurf.HalfEdge(k)
            _face, _id = self._half_edge_section[(-h).index()]
            S.change_edge_gluing(face, id, _face, _id)

        S.set_immutable()

        from flatsurf.geometry.translation_surface import TranslationSurface

        self._surface = TranslationSurface(S)
        self._pyflatsurf_surface = T

    def sage_flatsurf_surface(self):
        return self._surface

    def pyflatsurf_surface(self):
        return self._pyflatsurf_surface

    def half_edge_to_pyflatsurf(self, e):
        r"""
        INPUT:

        - ``e`` - a pair ``(face_label, edge_number)``
        """
        import pyflatsurf
        return pyflatsurf.flatsurf.HalfEdge(self._half_edge_map[e])

    def half_edge_from_pyflatsurf(self, h):
        return self._half_edge_section[h.index()]

    def point_to_pyflatsurf(self, p):
        raise NotImplementedError

    def point_from_pyflatsurf(self, p):
        raise NotImplementedError

    def flow_decomposition(self, v, limit=-1):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatsurfConverter
            sage: S = translation_surfaces.cathedral(1, 1)
            sage: f = S.pyflatsurf_converter() # optional: pyflatsurf
            sage: f.flow_decomposition((1, 0)) # optional: pyflatsurf
            FlowDecomposition with 5 cylinders, 0 minimal components and 0 undetermined components
        """
        v = self._V2(v)

        import pyflatsurf

        decomposition = pyflatsurf.flatsurf.makeFlowDecomposition(self._pyflatsurf_surface, v.vector)

        u = self._V2._isomorphic_vector_space(v)
        if limit != 0:
            decomposition.decompose(int(limit))
        return decomposition

    def flow_decompositions(self, bound, limit=-1, bfs=False):
        limit = int(limit)

        connections = self._pyflatsurf_surface.connections().bound(int(bound))
        if bfs:
            connections = connections.byLength()

        slopes = None

        from flatsurf.features import cppyy_feature
        cppyy_feature.require()
        import cppyy

        for connection in connections:
            direction = connection.vector()
            if slopes is None:
                slopes = cppyy.gbl.std.set[type(direction), type(direction).CompareSlope]()
            if slopes.find(direction) != slopes.end():
                continue
            slopes.insert(direction)
            yield self.flow_decomposition(direction, limit)

    def flow_decompositions_depth_first(self, bound, limit=-1):
        return self.flow_decompositions(bound, bfs=False, limit=limit)

    def flow_decompositions_breadth_first(self, bound, limit=-1):
        return self.flow_decompositions(bound, bfs=True, limit=limit)


def to_pyflatsurf(S):
    r"""
    Given S a translation surface from sage-flatsurf return a
    flatsurf::FlatTriangulation from libflatsurf/pyflatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
        sage: T = translation_surfaces.cathedral(1, 1)
        sage: to_pyflatsurf(T) # optional: pyflatsurf
        FlatTriangulationCombinatorial...
    """
    return S.pyflatsurf_converter().pyflatsurf_surface()


def sage_ring(surface):
    r"""
    Return the SageMath ring over which the pyflatsurf surface ``surface`` can
    be constructed in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, sage_ring # optional: pyflatsurf
        sage: S = to_pyflatsurf(translation_surfaces.veech_double_n_gon(5)) # optional: pyflatsurf
        sage: sage_ring(S) # optional: pyflatsurf
        Number Field in a with defining polynomial x^4 - 5*x^2 + 5 with a = 1.902113032590308?

    """
    from sage.all import Sequence
    vectors = [surface.fromHalfEdge(e.positive()) for e in surface.edges()]
    return Sequence([to_sage_ring(v.x()) for v in vectors] + [to_sage_ring(v.y()) for v in vectors]).universe()


def to_sage_ring(x):
    r"""
    Given a coordinate of a flatsurf::Vector, return a SageMath element from
    which :meth:`from_pyflatsurf` can eventually construct a translation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.pyflatsurf_conversion import to_sage_ring  # optional: pyflatsurf
        sage: to_sage_ring(1R).parent()  # optional: pyflatsurf
        Integer Ring

    GMP coordinate types::

        sage: import cppyy  # optional: pyflatsurf
        sage: import pyeantic  # optional: pyflatsurf
        sage: to_sage_ring(cppyy.gbl.mpz_class(1)).parent()  # optional: pyflatsurf
        Integer Ring
        sage: to_sage_ring(cppyy.gbl.mpq_class(1, 2)).parent()  # optional: pyflatsurf
        Rational Field

    e-antic coordinate types::

        sage: import pyeantic  # optional: pyflatsurf
        sage: K = pyeantic.eantic.renf_class.make("a^3 - 3*a + 1", "a", "0.34 +/- 0.01", 64R)  # optional: pyflatsurf
        sage: to_sage_ring(K.gen()).parent()  # optional: pyflatsurf
        Number Field in a with defining polynomial x^3 - 3*x + 1 with a = 0.3472963553338607?

    exact-real coordinate types::

        sage: from pyexactreal import QQModule, RealNumber  # optional: pyflatsurf
        sage: M = QQModule(RealNumber.random())   # optional: pyflatsurf
        sage: to_sage_ring(M.gen(0R)).parent()  # optional: pyflatsurf
        Real Numbers as (Rational Field)-Module

    """
    from flatsurf.features import cppyy_feature
    cppyy_feature.require()
    import cppyy

    def maybe_type(t):
        try:
            return t()
        except AttributeError:
            # The type constructed by t might not exist because the required C++ library has not been loaded.
            return None

    if type(x) is int:
        return ZZ(x)
    elif type(x) is maybe_type(lambda: cppyy.gbl.mpz_class):
        return ZZ(str(x))
    elif type(x) is maybe_type(lambda: cppyy.gbl.mpq_class):
        return QQ(str(x))
    elif type(x) is maybe_type(lambda: cppyy.gbl.eantic.renf_elem_class):
        from pyeantic import RealEmbeddedNumberField
        real_embedded_number_field = RealEmbeddedNumberField(x.parent())
        return real_embedded_number_field.number_field(real_embedded_number_field(x))
    elif type(x) is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.IntegerRing]):
        from pyexactreal import ExactReals
        return ExactReals(ZZ)(x)
    elif type(x) is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.RationalField]):
        from pyexactreal import ExactReals
        return ExactReals(QQ)(x)
    elif type(x) is maybe_type(lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.NumberField]):
        from pyexactreal import ExactReals
        return ExactReals(x.module().ring().parameters)(x)
    else:
        raise NotImplementedError(f"unknown coordinate ring for element {x} which is a {type(x)}")


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

    TESTS:

    Verify that #137 has been resolved::

        sage: from flatsurf import polygons
        sage: from flatsurf.geometry.surface import Surface_list
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: from flatsurf.geometry.gl2r_orbit_closure import GL2ROrbitClosure
        sage: from flatsurf.geometry.pyflatsurf_conversion import from_pyflatsurf
        sage: P = polygons.regular_ngon(10)
        sage: S = Surface_list(P.base_ring())
        sage: S.add_polygon(P)
        0
        sage: for i in range(5): S.set_edge_pairing(0, i, 0, 5+i)
        sage: M = TranslationSurface(S)
        sage: M.set_immutable()
        sage: f = M.pyflatsurf_converter() # optional: pyflatsurf
        sage: D0 = list(f.flow_decompositions(2))[2]  # optional: pyflatsurf
        sage: T0 = D0.triangulation()  # optional: pyflatsurf
        sage: S0 = from_pyflatsurf(T0)  # optional: pyflatsurf
        sage: S0 # optional: pyflatsurf
        TranslationSurface built from 8 polygons
        sage: to_pyflatsurf(S0) is T0 # optional: pyflatsurf
        True
    """
    C = FlatsurfConverter(T)
    C._surface._converter = C
    return C._surface
