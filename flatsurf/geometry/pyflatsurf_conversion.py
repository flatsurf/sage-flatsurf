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


def to_pyflatsurf(S):
    r"""
    Given S a translation surface from sage-flatsurf return a
    flatsurf::FlatTriangulation from libflatsurf/pyflatsurf.
    """
    from flatsurf.geometry.translation_surface import TranslationSurface

    if not isinstance(S, TranslationSurface):
        raise TypeError("S must be a translation surface")
    if not S.is_finite():
        raise ValueError("the surface S must be finite")

    S = S.triangulate()

    from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

    return Surface_pyflatsurf._from_flatsurf(S.underlying_surface())[
        0
    ]._flat_triangulation


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
        sage: X = GL2ROrbitClosure(M)  # optional: pyflatsurf
        sage: D0 = list(X.decompositions(2))[2]  # optional: pyflatsurf
        sage: T0 = D0.triangulation()  # optional: pyflatsurf
        sage: from_pyflatsurf(T0)  # optional: pyflatsurf
        TranslationSurface built from 8 polygons

    """
    from flatsurf.features import pyflatsurf_feature

    pyflatsurf_feature.require()
    import pyflatsurf

    ring = sage_ring(T)

    from flatsurf.geometry.surface import Surface_list

    S = Surface_list(ring)

    from flatsurf.geometry.polygon import ConvexPolygons

    P = ConvexPolygons(ring)

    V = P.module()

    half_edges = {}

    for face in T.faces():
        a, b, c = map(pyflatsurf.flatsurf.HalfEdge, face)

        vectors = [T.fromHalfEdge(he) for he in face]
        vectors = [
            V([ring(to_sage_ring(v.x())), ring(to_sage_ring(v.y()))]) for v in vectors
        ]
        triangle = P(vectors)
        face_id = S.add_polygon(triangle)

        assert a not in half_edges
        half_edges[a] = (face_id, 0)
        assert b not in half_edges
        half_edges[b] = (face_id, 1)
        assert c not in half_edges
        half_edges[c] = (face_id, 2)

    for half_edge, (face, id) in half_edges.items():
        _face, _id = half_edges[-half_edge]
        S.change_edge_gluing(face, id, _face, _id)

    S.set_immutable()

    from flatsurf.geometry.translation_surface import TranslationSurface

    return TranslationSurface(S)
