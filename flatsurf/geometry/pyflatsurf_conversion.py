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

from sage.misc.cachefunc import cached_method
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
    Conversion between a sage-flatsurf and an associated libflatsurf/pyflatsurf.

    It holds the following attributes

    - ``_domain`` : the source sage-flatsurf surface
    - ``_n``: number of half edges in a triangulation of ``_domain``
    - ``_half_edge_map``: dictionary ``(face_label, edge_number)`` -> ``pyflatsurf_edge_id``
    - ``_half_edge_section``: list ``pyflatsurf_edge_index`` -> ``(face_label, start_vertex, end_vertex)``
    - ``_fp``, ``_vp``: lists of respectively face permutation and vertex permutation
    """
    def __init__(self, surface, partial_triangulation=None):
        r"""
        INPUT:

        - ``surface`` -- a sage-flatsurf similarity surface

        - ``partial_triangulation`` -- an optional dictionary whose
          keys are polygon labels and values are list of pairs ``(i, j)``
          corresponding to edges that must appear in the resulting
          ``pyflatsurf`` triangulation

        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatsurfConverter
            sage: S = SymmetricGroup(3)
            sage: o = translation_surfaces.origami(S('(1,2)'), S('(1,3)'))
            sage: FlatsurfConverter(o).codomain()
            FlatTriangulationCombinatorial(vertices = (1, 7, -4, -6, -9, 4, -3, 8, 5, 3, -7, -2, 6, 9, 2, -1, -8, -5), faces = (1, 2, -7)(-1, -5, 8)(-2, 9, -6)(3, 4, 7)(-3, 5, -8)(-4, -9, 6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (1, 1), 8: (1, 1), 9: (1, 1)}
            sage: FlatsurfConverter(o, {1: [(0, 2)]}).codomain()
            FlatTriangulationCombinatorial(vertices = (1, 7, -4, -6, -9, 4, -3, 8, 5, 3, -7, -2, 6, 9, 2, -1, -8, -5), faces = (1, 2, -7)(-1, -5, 8)(-2, 9, -6)(3, 4, 7)(-3, 5, -8)(-4, -9, 6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (1, 1), 8: (1, 1), 9: (1, 1)}
            sage: FlatsurfConverter(o, {1: [(1, 3)]}).codomain()
            FlatTriangulationCombinatorial(vertices = (1, -4, -6, -9, 4, -7, -3, 8, 5, 3, -2, 6, 9, 2, 7, -1, -8, -5), faces = (1, 7, 4)(-1, -5, 8)(2, 3, -7)(-2, 9, -6)(-3, 5, -8)(-4, -9, 6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (-1, 1), 8: (1, 1), 9: (1, 1)}
            sage: FlatsurfConverter(o, {1: [(2, 0)]}).codomain()
            FlatTriangulationCombinatorial(vertices = (1, -7, -4, -6, -9, 4, -3, 8, 5, 3, 7, -2, 6, 9, 2, -1, -8, -5), faces = (1, 2, 7)(-1, -5, 8)(-2, 9, -6)(3, 4, -7)(-3, 5, -8)(-4, -9, 6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (-1, -1), 8: (1, 1), 9: (1, 1)}
            sage: FlatsurfConverter(o, {1: [(3, 1)]}).codomain()
            FlatTriangulationCombinatorial(vertices = (1, -4, -6, -9, 4, 7, -3, 8, 5, 3, -2, 6, 9, 2, -7, -1, -8, -5), faces = (1, -7, 4)(-1, -5, 8)(2, 3, 7)(-2, 9, -6)(-3, 5, -8)(-4, -9, 6)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, 0), 4: (0, -1), 5: (0, 1), 6: (1, 0), 7: (1, -1), 8: (1, 1), 9: (1, 1)}
        """
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if not isinstance(surface, SimilaritySurface):
            raise ValueError("surface must be a similarity surface")
        from flatsurf.geometry.translation_surface import TranslationSurface
        if not isinstance(surface, TranslationSurface):
            raise NotImplementedError("only implemented for translation surfaces")
        if not surface.is_finite():
            raise ValueError("the surface must be finite")
        if surface.is_mutable():
            raise ValueError("the surface must be immutable")

        self._domain = surface

        # number of half-edges in the triangulation
        self._n = 3 * sum(surface.polygon(label).num_edges() - 2 for label in surface.label_iterator())

        self._half_edge_map = {} # (face label, edge number) -> id
        self._half_edge_section = [] # index -> (face label, vertex start, vertex end)
        self._fp = [None] * (self._n + 1)

        # populate half edges in the surface
        m = 1
        vectors = []
        for t0, t1 in surface.edge_iterator(gluings=True):
            if t0 in self._half_edge_map:
                continue
            p0 = surface.polygon(t0[0])
            p1 = surface.polygon(t1[0])
            self._half_edge_map[t0] = m
            self._half_edge_map[t1] = -m
            self._half_edge_section.append((t0[0], t0[1], (t0[1] + 1) % p0.num_edges()))
            self._half_edge_section.append((t1[0], t1[1], (t1[1] + 1) % p1.num_edges()))
            m += 1

        # possibly triangulate and compute face permutation
        from .polygon import triangulate, build_faces
        for label in surface.label_iterator():
            poly = surface.polygon(label)
            ne = poly.num_edges()
            if ne == 3:
                # already a triangle
                for i in range(3):
                    self._fp[self._half_edge_map[(label, i)]] = self._half_edge_map[(label, (i+1) % 3)]
            else:
                local_edges = {}
                for i in range(ne):
                    local_edges[(i, (i+1)%ne)] = self._half_edge_map[(label, i)]
                forced_edges = [] if partial_triangulation is None else partial_triangulation.get(label, [])
                for a, b in forced_edges:
                    local_edges[(a, b)] = m
                    local_edges[(b, a)] = -m
                    self._half_edge_section.append((label, a, b))
                    self._half_edge_section.append((label, b, a))
                    m += 1
                for face in build_faces(ne, forced_edges):
                    additional_edges = triangulate([poly.vertex(a) for a in face])
                    for i, j in additional_edges:
                        a = face[i]
                        b = face[j]
                        local_edges[(a, b)] = m
                        local_edges[(b, a)] = -m
                        self._half_edge_section.append((label, a, b))
                        self._half_edge_section.append((label, b, a))
                        m += 1
                    for (i, j, k) in build_faces(len(face), additional_edges):
                        a = face[i]
                        b = face[j]
                        c = face[k]
                        self._fp[local_edges[(a, b)]] = local_edges[(b, c)]
                        self._fp[local_edges[(b, c)]] = local_edges[(c, a)]
                        self._fp[local_edges[(c, a)]] = local_edges[(a, b)]

        assert m == 1 + self._n // 2, (m, self._n)
        assert all(self._fp[i] is not None and self._fp[-i] is not None for i in range(1, self._n//2 + 1))
        assert len(self._half_edge_section) == self._n

        # compute vertex permutation from face permutation
        self._vp = [None] * (self._n + 1)
        for i in range(1, self._n // 2 + 1):
            self._vp[self._fp[-i]] = i
            self._vp[self._fp[i]] = -i

    @staticmethod
    def from_pyflatsurf(T):
        r"""
        Return a :class:`FlatsurfConverter` with a given pyflatsurf/libflatsurf FlatTriangulation as target.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatsurfConverter
            sage: from pyflatsurf.vector import Vectors # optional: pyflatsurf
            sage: from pyflatsurf.factory import make_surface # optional: pyflatsurf
            sage: surfaces = []
            sage: V2 = Vectors(QQ) # optional: pyflatsurf
            sage: vertices = [(1r, -3r, 2r, -1r, 6r, -5r), (-2r, 4r, -6r, 5r, -4r, 3r)]
            sage: vectors = [V2((1, 1)), V2((-1, 0)), V2((0, -1)), V2((1, 1)), V2((-1, 0)), V2((0, -1))] # optional: pyflatsurf
            sage: T = make_surface(vertices, [v.vector for v in vectors]) # optional: pyflatsurf
            sage: f = FlatsurfConverter.from_pyflatsurf(T) # optional: pyflatsurf
            sage: f.domain() # optional: pyflatsurf
            TranslationSurface built from 4 polygons
            sage: f.codomain() == T # optional: pyflatsurf
            True
        """
        n = len(T.halfEdges())

        ring = sage_ring(T)

        from flatsurf.geometry.surface import Surface_list
        S = Surface_list(ring)

        from flatsurf.geometry.polygon import ConvexPolygons
        P = ConvexPolygons(ring)
        V = P.module()

        half_edge_map = {}
        half_edge_section = [None] * n

        from flatsurf.features import pyflatsurf_feature
        pyflatsurf_feature.require()
        import pyflatsurf
        for face in T.faces():
            a, b, c = map(pyflatsurf.flatsurf.HalfEdge, face)

            vectors = [T.fromHalfEdge(he) for he in face]
            vectors = [V([ring(to_sage_ring(v.x())), ring(to_sage_ring(v.y()))]) for v in vectors]
            triangle = P(vectors)
            face_id = S.add_polygon(triangle)

            half_edge_map[(face_id, 0)] = a.id()
            half_edge_map[(face_id, 1)] = b.id()
            half_edge_map[(face_id, 2)] = c.id()

            assert half_edge_section[a.index()] is None
            half_edge_section[a.index()] = (face_id, 0, 1)
            assert half_edge_section[b.index()] is None
            half_edge_section[b.index()] = (face_id, 1, 2)
            assert half_edge_section[c.index()] is None
            half_edge_section[c.index()] = (face_id, 2, 0)

        assert all(e is not None for e in half_edge_section), half_edge_section

        for (face, id), k in half_edge_map.items():
            h = pyflatsurf.flatsurf.HalfEdge(k)
            _face, _id, _id_next = half_edge_section[(-h).index()]
            S.change_edge_gluing(face, id, _face, _id)
        S.set_immutable()
        from flatsurf.geometry.translation_surface import TranslationSurface
        surface = TranslationSurface(S)
        surface.set_immutable()

        vp = [None] * (n + 1)
        fp = [None] * (n + 1)
        for h in T.halfEdges():
            vp[h.id()] = T.nextAtVertex(h).id()
            fp[h.id()] = T.nextInFace(h).id()

        f = FlatsurfConverter.__new__(FlatsurfConverter)
        f._domain = surface
        f._n = n
        f._half_edge_map = half_edge_map
        f._half_edge_section = half_edge_section
        f._vp = vp
        f._fp = fp

        return f

    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError('incomparable objects {} and {}'.format(type(self), type(other)))
        return self._domain == other._domain and self._half_edge_section == other._half_edge_section

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError('incomparable objects {} and {}'.format(type(self), type(other)))
        return self._domain != other._domain or self._half_edge_section != other._half_edge_section

    def domain(self):
        r"""
        Return the ``sage-flatsurf`` surface which is the source of this converter.
        """
        return self._domain

    @cached_method
    def base_ring(self):
        r"""
        Return a SageMath base ring that contains all the edge coordinates
        """
        S = self._domain
        base_ring = S.base_ring()
        if base_ring is AA:
            from sage.rings.qqbar import number_field_elements_from_algebraics
            from itertools import chain
            coords = []
            for label in S.label_iterator():
                for e in S.polygon(label).edges():
                    coords.append(e[0])
                    coords.append(e[1])
            base_ring = number_field_elements_from_algebraics(coords, embedded=True)[0]
        return base_ring

    @cached_method
    def vector_space(self):
        r"""
        A model of the vector space R² in libflatsurf, e.g., to represent the
        vector associated to a saddle connection.
        """
        from flatsurf.features import pyflatsurf_feature
        pyflatsurf_feature.require()
        from pyflatsurf.vector import Vectors
        return Vectors(self.base_ring())

    @cached_method
    def vectors(self):
        r"""
        Return the list of holonomy vectors associated to the (positively oriented) edges.
        """
        V = self.vector_space()
        vectors = []
        for m in range(0, len(self._half_edge_section), 2):
            label, a, b = self._half_edge_section[m]
            poly = self._domain.polygon(label)
            vectors.append(V(poly.vertex(b) - poly.vertex(a)))
        return vectors

    # TODO: this should only be weakly cached
    @cached_method
    def codomain(self):
        r"""
        Return the ``pyflatsurf`` surface corresponding to the target of this converter.
        """
        from pyflatsurf.factory import make_surface
        vectors = [v.vector for v in self.vectors()]
        _check_data(self._vp, self._fp, vectors)
        return make_surface(_cycle_decomposition(self._vp), vectors)

    def half_edge_to_pyflatsurf(self, e):
        r"""
        INPUT:

        - ``e`` - a pair ``(face_label, edge_number)``
        """
        import pyflatsurf
        return pyflatsurf.flatsurf.HalfEdge(self._half_edge_map[e])

    def half_edge_from_pyflatsurf(self, h):
        label, a, b = self._half_edge_section[h.index()]
        poly = self.domain().polygon(label)
        return (label, a) if b == (a + 1) % poly.num_edges() else None

    def saddle_connection_to_pyflatsurf(self, sc):
        raise NotImplementedError

    def saddle_connection_from_pyflatsurf(self, sc, check=True):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.square_torus()
            sage: for sc in s._pyflatsurf.codomain().connections().bound(3r):
            ....:     print(s._pyflatsurf.saddle_connection_from_pyflatsurf(sc))
            Saddle connection ...
            ...

            sage: set(s.saddle_connections(3*3)) == set(s._pyflatsurf.saddle_connection_from_pyflatsurf(sc) for sc in s._pyflatsurf.codomain().connections().bound(3r))
            True
        """
        h_start = sc.source()
        f_start, e_start, _ = self._half_edge_section[h_start.index()]
        h_end = sc.target()
        f_end, e_end, _ = self._half_edge_section[h_end.index()]
        v = self.vector_space()(sc.vector())
        v = self.vector_space()._isomorphic_vector_space(v)
        S = self.domain()
        poly = S.polygon(f_start)
        from .surface_objects import SaddleConnection
        return SaddleConnection(surface=S, start_data=(f_start, e_start),
                direction=v, end_data=(f_end, e_end), end_direction=-v,
                holonomy=v, end_holonomy=-v, check=check)

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
            sage: f = S._pyflatsurf # optional: pyflatsurf
            sage: f.flow_decomposition((1, 0)) # optional: pyflatsurf
            FlowDecomposition with 5 cylinders, 0 minimal components and 0 undetermined components
        """
        V2 = self.vector_space()
        v = V2(v)

        import pyflatsurf

        decomposition = pyflatsurf.flatsurf.makeFlowDecomposition(self.codomain(), v.vector)

        u = V2._isomorphic_vector_space(v)
        if limit != 0:
            decomposition.decompose(int(limit))
        return decomposition

    def flow_decompositions(self, bound, limit=-1, bfs=False):
        limit = int(limit)

        connections = self.codomain().connections().bound(int(bound))
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
    return S._pyflatsurf.codomain()


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
        sage: f = M._pyflatsurf # optional: pyflatsurf
        sage: D0 = list(f.flow_decompositions(2))[2]  # optional: pyflatsurf
        sage: T0 = D0.triangulation()  # optional: pyflatsurf
        sage: from_pyflatsurf(T0)  # optional: pyflatsurf
        TranslationSurface built from 8 polygons
    """
    return FlatsurfConverter.from_pyflatsurf(T).domain()
