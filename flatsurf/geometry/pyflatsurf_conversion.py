r"""
Conversion of sage-flatsurf objects to libflatsurf/pyflatsurf and vice versa.

Ideally, there should be no need to call the functionality in this module
directly. Interaction with libflatsurf/pyflatsurf should be handled
transparently by the sage-flatsurf objects. Even for authors of sage-flatsurf
there should essentially never be a need to use this module directly since
objects should provide a ``_pyflatsurf()`` method that returns a
:class:`Conversion` to libflatsurf/pyflatsurf.

EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
    sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
    sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)  # random output due to deprecation warnings

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C)      2019 Vincent Delecroix
#                      2019-2023 Julian Rüth
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

from flatsurf.features import pyflatsurf_feature, pyeantic_feature


class Conversion:
    r"""
    Generic base class for a conversion from sage-flatsurf to pyflatsurf.

    INPUT:

    - ``domain`` -- the domain of this conversion, can be ``None``, when the
      domain cannot be represented that easily with a sage-flatsurf object.

    - ``codomain`` -- the codomain of this conversion, can be ``None``, when
      the codomain cannot be represented that easily with
      libflatsurf/pyflatsurf object.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
        sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf_conversion import Conversion
        sage: isinstance(conversion, Conversion)
        True

    """

    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    def domain(self):
        r"""
        Return the domain of this conversion, a sage-flatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)
            sage: conversion.domain()
            <flatsurf.geometry.surface.Surface_dict object at 0x...>

        """
        if self._domain is not None:
            return self._domain

        raise NotImplementedError("this conversion does not implement domain() yet")

    def codomain(self):
        r"""
        Return the codomain of this conversion, a pyflatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)
            sage: conversion.codomain()
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, -8, 9, 7, -6, -9, 4, 3, -2, -4, 8, 5, -7, 6, -5), faces = (1, 2, 3)(-1, -5, 8)(-2, -3, 4)(-4, -9, -8)(5, 6, 7)(-6, -7, 9)) with vectors {...}

        """
        if self._codomain is not None:
            return self._codomain

        raise NotImplementedError("this conversion does not implement codomain() yet")

    def __call__(self, x):
        r"""
        Return the conversion at an element of :meth:`domain` and return the
        corresponding pyflatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: conversion(p)

        """
        raise NotImplementedError("this conversion does not implement a mapping of elements yet")

    def section(self, y):
        r"""
        Return the conversion of an element of :meth:`codomain` and return the
        corresponding sage-flatsurf object.

        This is the inverse of :meth:`__call__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion.section(q)

        """
        raise NotImplementedError("this conversion does not implement a section yet")

    def __repr__(self):
        r"""
        Return a printable representation of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConverter.to_pyflatsurf(S)
            Conversion from <flatsurf.geometry.surface.Surface_dict object at 0x...> to
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, -8, 9, 7, -6, -9, 4, 3, -2, -4, 8, 5, -7, 6, -5),
                                           faces = (1, 2, 3)(-1, -5, 8)(-2, -3, 4)(-4, -9, -8)(5, 6, 7)(-6, -7, 9))
                                           with vectors {1: ((1/2 ~ 0.50000000), (1/2*a^3 - 1*a ~ 1.5388418)),
                                                         2: ((-1/2*a^2 + 1 ~ -0.80901699), (-1/2*a^3 + 3/2*a ~ -0.58778525)),
                                                         3: ((1/2*a^2 - 3/2 ~ 0.30901699), (-1/2*a ~ -0.95105652)),
                                                         4: ((-1/2 ~ -0.50000000), (-1/2*a^3 + 1*a ~ -1.5388418)),
                                                         5: ((-1/2*a^2 + 1/2 ~ -1.3090170), (-1/2*a ~ -0.95105652)),
                                                         6: (1, 0),
                                                         7: ((1/2*a^2 - 3/2 ~ 0.30901699), (1/2*a ~ 0.95105652)),
                                                         8: ((-1/2*a^2 + 1 ~ -0.80901699), (1/2*a^3 - 3/2*a ~ 0.58778525)),
                                                         9: ((1/2*a^2 - 1/2 ~ 1.3090170), (1/2*a ~ 0.95105652))}

        """
        codomain = self.codomain()

        if hasattr(codomain, "__cpp_name__"):
            codomain = codomain.__cpp_name__

        return f"Conversion from {self.domain()} to {codomain}"


class RingConversion(Conversion):
    r"""
    A conversion between a SageMath ring and a C/C++ ring.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConverter
        sage: conversion = RingConverter.to_pyflatsurf_from_elements([1, 2, 3])

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: isinstance(conversion, RingConversion)
        True

    """


class FlatTriangulationConversion(Conversion):
    def __init__(self, domain, codomain, label_to_half_edges):
        super().__init__(domain, codomain)

        self._label_to_half_edges = label_to_half_edges


class Converter:
    r"""
    Abstract base class for converters from sage-flatsurf to
    pyflatsurf/libflatsurf.
    """
    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to the ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

        """
        raise NotImplementedError("this converter does not implement conversion to pyflatsurf yet")

    @classmethod
    def from_pyflatsurf(cls, codomain, domain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)
            sage: FlatTriangulationConverter.from_pyflatsurf(conversion.codomain())

        """
        raise NotImplementedError("this converter does not implement conversion from pyflatsurf yet")

    @classmethod
    def to_pyflatsurf_from_elements(cls, elements, codomain=None):
        r"""
        Return a :class:`Conversion` that converts the sage-flatsurf
        ``elements`` to  ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConverter
            sage: RingConverter.to_pyflatsurf_from_elements([1, 2, 3])
            Conversion from Integer Ring to __gmp_expr<__mpz_struct[1],__mpz_struct[1]>

        """
        from sage.all import Sequence

        return cls.to_pyflatsurf(domain=Sequence(elements).universe(), codomain=codomain)


class RingConverter(Converter):
    r"""
    Converts elements of a SageMath ring to something that
    libflatsurf/pyflatsurf can understand and vice-versa.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConverter
        sage: RingConverter.to_pyflatsurf(ZZ)
        Conversion from Integer Ring to __gmp_expr<__mpz_struct[1],__mpz_struct[1]>

    """
    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`Conversion` that converts the SageMath ring ``domain``
        to something that libflatsurf/pyflatsurf can understand.

        INPUT:

        - ``domain`` -- a ring

        - ``codomain`` -- a C/C++ type or ``None`` (default: ``None``); if
          ``None``, the corresponding type is constructed.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConverter
            sage: RingConverter.to_pyflatsurf(QQ)
            Conversion from Rational Field to __gmp_expr<__mpq_struct[1],__mpq_struct[1]>

        """
        if codomain is None:
            from sage.all import ZZ, QQ, NumberFields, RR
            if domain is ZZ:
                import cppyy
                codomain = cppyy.gbl.mpz_class
            elif domain is QQ:
                import cppyy
                codomain = cppyy.gbl.mpq_class
            elif domain in NumberFields():
                if not domain.embeddings(RR):
                    raise NotImplementedError("cannot determine pyflatsurf ring for not real-embedded number fields yet")

                pyeantic_feature.require()

                from pyeantic import RealEmbeddedNumberField
                codomain = RealEmbeddedNumberField(domain)
            else:
                raise NotImplementedError(f"cannot determine pyflatsurf ring corresponding to {domain} yet")

        return RingConversion(domain, codomain)

    @classmethod
    def from_pyflatsurf(cls, codomain, domain=None):
        raise NotImplementedError


class FlatTriangulationConverter(Converter):
    r"""
    Converts :class:`Surface` objects to ``FlatTriangulation`` objects and
    vice-versa.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
        sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

    """
    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to the ``codomain``.

        INPUT:

        - ``domain`` -- a :class:`Surface`

        - ``codomain`` -- a ``FlatTriangulation`` or ``None`` (default:
          ``None``); if ``None``, the corresponding ``FlatTriangulation`` is
          constructed.

        .. NOTE::

            The ``codomain``, if given, must be indistinguishable from the
            codomain that this method would construct automatically.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)

        """
        pyflatsurf_feature.require()

        from flatsurf.geometry.surface import Surface

        if not isinstance(domain, Surface):
            raise TypeError("domain must be a Surface")
        if not domain.is_finite():
            raise ValueError("domain must be finite")
        if not domain.is_triangulated():
            raise ValueError("domain must be triangulated")

        if codomain is None:
            vertex_permutation = cls._pyflatsurf_vertex_permutation(domain)
            vectors = cls._pyflatsurf_vectors(domain)

            from pyflatsurf.factory import make_surface
            codomain = make_surface(vertex_permutation, vectors)

        return FlatTriangulationConversion(domain, codomain, cls._pyflatsurf_labels(domain))

    @classmethod
    def _pyflatsurf_labels(cls, domain):
        r"""
        Return a mapping of the edges of the polygons of ``domain`` to half
        edges numbered compatibly with libflatsurf/pyflatsurf, i.e., by
        consecutive integers such that the opposite of an edge has the negative
        of that edge's label.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConverter._pyflatsurf_labels(S)
            {(0, 0): 1,
             (0, 1): 2,
             (0, 2): 3,
             (1, 0): 4,
             (1, 1): -2,
             (1, 2): -3,
             (ExtraLabel(...), 0): 5,
             (ExtraLabel(...), 1): 6,
             (ExtraLabel(...), 2): 7,
             (ExtraLabel(...), 0): -1,
             (ExtraLabel(...), 1): -5,
             (ExtraLabel(...), 2): 8,
             (ExtraLabel(...), 0): 9,
             (ExtraLabel(...), 1): -6,
             (ExtraLabel(...), 2): -7,
             (ExtraLabel(...), 0): -4,
             (ExtraLabel(...), 1): -9,
             (ExtraLabel(...), 2): -8}

        """
        labels = {}
        for half_edge, opposite_half_edge in domain.edge_gluing_iterator():
            if half_edge in labels:
                continue

            labels[half_edge] = len(labels) // 2 + 1
            labels[opposite_half_edge] = -labels[half_edge]

        return labels

    @classmethod
    def _pyflatsurf_vectors(cls, domain):
        r"""
        Return the vectors of the positive half edges in ``domain`` in order.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConverter._pyflatsurf_vectors(S)
            [((1/2 ~ 0.50000000), (1/2*a^3 - 1*a ~ 1.5388418)),
             ((-1/2*a^2 + 1 ~ -0.80901699), (-1/2*a^3 + 3/2*a ~ -0.58778525)),
             ((1/2*a^2 - 3/2 ~ 0.30901699), (-1/2*a ~ -0.95105652)),
             ((-1/2 ~ -0.50000000), (-1/2*a^3 + 1*a ~ -1.5388418)),
             ((-1/2*a^2 + 1/2 ~ -1.3090170), (-1/2*a ~ -0.95105652)),
             (1, 0),
             ((1/2*a^2 - 3/2 ~ 0.30901699), (1/2*a ~ 0.95105652)),
             ((-1/2*a^2 + 1 ~ -0.80901699), (1/2*a^3 - 3/2*a ~ 0.58778525)),
             ((1/2*a^2 - 1/2 ~ 1.3090170), (1/2*a ~ 0.95105652))]

        """
        labels = cls._pyflatsurf_labels(domain)

        vectors = [None] * (len(labels) // 2)

        for (polygon, edge), half_edge in labels.items():
            if half_edge < 0:
                continue

            vectors[half_edge - 1] = domain.polygon(polygon).edge(edge)

        ring_conversion = RingConverter.to_pyflatsurf_from_elements([vector[0] for vector in vectors] + [vector[1] for vector in vectors])

        from pyflatsurf.vector import Vectors
        vector_space = Vectors(ring_conversion.codomain())

        return [vector_space(vector).vector for vector in vectors]

    @classmethod
    def _pyflatsurf_vertex_permutation(cls, domain):
        r"""
        Return the permutation of half edges around vertices of ``domain`` in
        cycle notation.

        The permutation uses integers, as provided by
        :meth:`_pyflatsurf_labels`.

        INPUT:

        - ``domain`` -- a :class:`Surface`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConverter._pyflatsurf_vertex_permutation(S)
            [[1, -3, 2, -1, -8, 9, 7, -6, -9, 4, 3, -2, -4, 8, 5, -7, 6, -5]]

        """
        pyflatsurf_labels = cls._pyflatsurf_labels(domain)

        vertex_permutation = {}

        for polygon, edge in domain.edge_iterator():
            pyflatsurf_edge = pyflatsurf_labels[(polygon, edge)]

            next_edge = (edge + 1) % domain.polygon(polygon).num_edges()
            pyflatsurf_next_edge = pyflatsurf_labels[(polygon, next_edge)]

            vertex_permutation[pyflatsurf_next_edge] = -pyflatsurf_edge

        return cls._cycle_decomposition(vertex_permutation)

    @classmethod
    def _cycle_decomposition(self, permutation):
        r"""
        Return a permutation in cycle notation.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: FlatTriangulationConverter._cycle_decomposition({})
            []
            sage: FlatTriangulationConverter._cycle_decomposition({1: 1, -1: -1})
            [[1], [-1]]
            sage: FlatTriangulationConverter._cycle_decomposition({1: 2, 2: 1, 3: 4, 4: 5, 5: 3})
            [[1, 2], [3, 4, 5]]

        """
        cycles = []

        elements = set(permutation.keys())

        while elements:
            cycle = [elements.pop()]
            while True:
                cycle.append(permutation[cycle[-1]])

                if cycle[-1] == cycle[0]:
                    cycle.pop()
                    cycles.append(cycle)
                    break

                elements.remove(cycle[-1])

        return cycles

    @classmethod
    def from_pyflatsurf(cls, codomain, domain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to ``codomain``.

        INPUT:

        - ``codomain`` -- a ``FlatTriangulation``

        - ``domain`` -- a :class:`Surface` or ``None`` (default: ``None``); if
          ``None``, the corresponding :class:`Surface` is constructed.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConverter
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConverter.to_pyflatsurf(S)
            sage: FlatTriangulationConverter.from_pyflatsurf(conversion.codomain())

        """
        pyflatsurf_feature.require()

        if domain is not None:
            return cls.to_pyflatsurf(domain=domain, codomain=codomain)

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




def to_pyflatsurf(S):
    r"""
    Given S a translation surface from sage-flatsurf return a
    flatsurf::FlatTriangulation from libflatsurf/pyflatsurf.
    """
    import warnings
    warnings.warn("to_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use FlatTriangulationConverter.to_pyflatsurf(surface.triangulate().underlying_surface()).codomain() instead.")

    return FlatTriangulationConverter.to_pyflatsurf(S).codomain()


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
    # TODO: Remove this method without deprecation. We are almost certain that nobody was using it.
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
    # TODO: Remove this method without deprecation. We are almost certain that nobody was using it.
    import warnings
    warnings.warn("to_sage_ring() is deprecated and will be removed in a future version of sage-flaturf. Use RingConverter.from_pyflatsurf_element(x).section(x) instead.")

    return RingConverter.from_pyflatsurf_element(x).section(x)

    from flatsurf.features import cppyy_feature

    cppyy_feature.require()
    import cppyy

    def maybe_type(t):
        try:
            return t()
        except AttributeError:
            # The type constructed by t might not exist because the required C++ library has not been loaded.
            return None

    from sage.all import ZZ
    if type(x) is int:
        return ZZ(x)
    elif type(x) is maybe_type(lambda: cppyy.gbl.mpz_class):
        return ZZ(str(x))
    elif type(x) is maybe_type(lambda: cppyy.gbl.mpq_class):
        from sage.all import QQ
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
    import warnings
    warnings.warn("from_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use TranslationSurface(FlatTriangulationConverter.from_pyflatsurf(surface).domain()) instead.")

    return TranslationSurface(FlatTriangulationConverter.from_pyflatsurf(T).domain())
