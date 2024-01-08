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
    assert len(fp) == n + 1
    assert len(vec) == n // 2

    assert vp[0] is None
    assert fp[0] is None

    for i in range(1, n // 2 + 1):
        # check fp/vp consistency
        assert fp[-vp[i]] == i, i

        # check that each face is a triangle and that vectors sum up to zero
        j = fp[i]
        k = fp[j]
        assert i != j and i != k and fp[k] == i, (i, j, k)
        vi = vec[i - 1] if i >= 1 else -vec[-i - 1]
        vj = vec[j - 1] if j >= 1 else -vec[-j - 1]
        vk = vec[k - 1] if k >= 1 else -vec[-k - 1]
        v = vi + vj + vk
        assert v.x() == 0, v.x()
        assert v.y() == 0, v.y()


def _cycle_decomposition(p):
    n = len(p) - 1
    assert n % 2 == 0
    cycles = []
    unseen = [True] * (n + 1)
    for i in list(range(-n // 2 + 1, 0)) + list(range(1, n // 2)):
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
    from flatsurf.geometry.categories import TranslationSurfaces

    if S not in TranslationSurfaces():
        raise TypeError("S must be a translation surface")
    if not S.is_finite_type():
        raise ValueError("the surface S must be finite")

    S = S.triangulate()

    # populate half edges and vectors
    n = sum(len(S.polygon(lab).vertices()) for lab in S.labels())
    half_edge_labels = {}  # map: (face lab, edge num) in faltsurf -> integer
    vec = []  # vectors
    k = 1  # half edge label in {1, ..., n}
    for t0, t1 in S.gluings():
        if t0 in half_edge_labels:
            continue

        half_edge_labels[t0] = k
        half_edge_labels[t1] = -k

        f0, e0 = t0
        p = S.polygon(f0)
        vec.append(p.edge(e0))

        k += 1

    # compute vertex and face permutations
    vp = [None] * (n + 1)  # vertex permutation
    fp = [None] * (n + 1)  # face permutation
    for t in S.edges():
        e = half_edge_labels[t]
        j = (t[1] + 1) % len(S.polygon(t[0]).vertices())
        fp[e] = half_edge_labels[(t[0], j)]
        vp[fp[e]] = -e

    # convert the vp permutation into cycles
    verts = _cycle_decomposition(vp)

    # find a finite SageMath base ring that contains all the coordinates
    base_ring = S.base_ring()
    if base_ring is AA:
        from sage.rings.qqbar import number_field_elements_from_algebraics
        from itertools import chain

        base_ring = number_field_elements_from_algebraics(
            list(chain(*[list(v) for v in vec])), embedded=True
        )[0]

    from flatsurf.features import pyflatsurf_feature

    pyflatsurf_feature.require()
    from pyflatsurf.vector import Vectors

    V = Vectors(base_ring)
    vec = [V(v).vector for v in vec]

    _check_data(vp, fp, vec)

    from pyflatsurf.factory import make_surface

    return make_surface(verts, vec)


def sage_ring(surface):
    r"""
    Return the SageMath ring over which the pyflatsurf surface ``surface`` can
    be constructed in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, sage_ring # optional: pyflatsurf
        sage: S = to_pyflatsurf(translation_surfaces.veech_double_n_gon(5)) # optional: pyflatsurf  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: sage_ring(S) # optional: pyflatsurf
        Number Field in a with defining polynomial x^4 - 5*x^2 + 5 with a = 1.902113032590308?

    """
    from sage.all import Sequence

    vectors = [surface.fromHalfEdge(e.positive()) for e in surface.edges()]
    return Sequence(
        [to_sage_ring(v.x()) for v in vectors] + [to_sage_ring(v.y()) for v in vectors]
    ).universe()


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
    elif type(x) is maybe_type(
        lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.IntegerRing]
    ):
        from pyexactreal import ExactReals

        return ExactReals(ZZ)(x)
    elif type(x) is maybe_type(
        lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.RationalField]
    ):
        from pyexactreal import ExactReals

        return ExactReals(QQ)(x)
    elif type(x) is maybe_type(
        lambda: cppyy.gbl.exactreal.Element[cppyy.gbl.exactreal.NumberField]
    ):
        from pyexactreal import ExactReals

        return ExactReals(x.module().ring().parameters)(x)
    else:
        raise NotImplementedError(
            f"unknown coordinate ring for element {x} which is a {type(x)}"
        )


def from_pyflatsurf(T):
    r"""
    Given T a flatsurf::FlatTriangulation from libflatsurf/pyflatsurf, return a
    sage-flatsurf translation surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, from_pyflatsurf # optional: pyflatsurf
        sage: S = translation_surfaces.veech_double_n_gon(5) # optional: pyflatsurf
        sage: T = from_pyflatsurf(to_pyflatsurf(S)) # optional: pyflatsurf
        sage: T  # optional: pyflatsurf
        Translation Surface in H_2(2) built from 6 isosceles triangles

    TESTS::

        sage: from flatsurf.geometry.categories import TranslationSurfaces
        sage: T in TranslationSurfaces()  # optional: pyflatsurf
        True

    Verify that #137 has been resolved::

        sage: from flatsurf import polygons, MutableOrientedSimilaritySurface
        sage: from flatsurf.geometry.gl2r_orbit_closure import GL2ROrbitClosure
        sage: from flatsurf.geometry.pyflatsurf_conversion import from_pyflatsurf
        sage: P = polygons.regular_ngon(10)
        sage: S = MutableOrientedSimilaritySurface(P.base_ring())
        sage: S.add_polygon(P)
        0
        sage: for i in range(5): S.glue((0, i), (0, 5+i))
        sage: S.set_immutable()
        sage: M = S
        sage: X = GL2ROrbitClosure(M)  # optional: pyflatsurf
        sage: D0 = list(X.decompositions(2))[2]  # optional: pyflatsurf
        sage: T0 = D0.triangulation()  # optional: pyflatsurf
        sage: from_pyflatsurf(T0)  # optional: pyflatsurf
        Translation Surface in H_2(1^2) built from 2 isosceles triangles and 6 triangles

    """
    from flatsurf.features import pyflatsurf_feature

    pyflatsurf_feature.require()
    import pyflatsurf

    ring = sage_ring(T)

    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

    S = MutableOrientedSimilaritySurface(ring)

    from flatsurf.geometry.polygon import Polygon

    V = ring**2

    half_edges = {}

    for face in T.faces():
        a, b, c = map(pyflatsurf.flatsurf.HalfEdge, face)

        vectors = [T.fromHalfEdge(he) for he in face]
        vectors = [
            V([ring(to_sage_ring(v.x())), ring(to_sage_ring(v.y()))]) for v in vectors
        ]
        triangle = Polygon(edges=vectors)
        face_id = S.add_polygon(triangle)

        assert a not in half_edges
        half_edges[a] = (face_id, 0)
        assert b not in half_edges
        half_edges[b] = (face_id, 1)
        assert c not in half_edges
        half_edges[c] = (face_id, 2)

    for half_edge, (face, id) in half_edges.items():
        _face, _id = half_edges[-half_edge]
        S.glue((face, id), (_face, _id))

    S.set_immutable()
    return S
