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
    sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
    sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
    sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)  # random output due to deprecation warnings

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

from sage.misc.cachefunc import cached_method

from flatsurf.features import pyflatsurf_feature


# TODO: Move this fix to pyeantic
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pyeantic

    def unwrap_intrusive_ptr(K):
        import cppyy
        if isinstance(K, pyeantic.eantic.renf_class):
            K = cppyy.gbl.boost.intrusive_ptr['const eantic::renf_class'](K)
        if isinstance(K, cppyy.gbl.boost.intrusive_ptr['const eantic::renf_class']):
            ptr = K.get()
            ptr.__intrusive__ = K
            K = ptr
        return K

    import pyeantic.cppyy_eantic
    pyeantic.cppyy_eantic.unwrap_intrusive_ptr = unwrap_intrusive_ptr

    pyeantic.eantic.renf = lambda *args: unwrap_intrusive_ptr(pyeantic.eantic.renf_class.make(*args))


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
        sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
        sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf_conversion import Conversion
        sage: isinstance(conversion, Conversion)
        True

    """

    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to the ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

        """
        raise NotImplementedError("this converter does not implement conversion to pyflatsurf yet")

    @classmethod
    def from_pyflatsurf(cls, codomain, domain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: FlatTriangulationConversion.from_pyflatsurf(conversion.codomain())
            Conversion from ...

        """
        raise NotImplementedError("this converter does not implement conversion from pyflatsurf yet")

    @classmethod
    def to_pyflatsurf_from_elements(cls, elements, codomain=None):
        r"""
        Return a :class:`Conversion` that converts the sage-flatsurf
        ``elements`` to  ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: RingConversion.to_pyflatsurf_from_elements([1, 2, 3])
            Conversion from Integer Ring to __gmp_expr<__mpz_struct[1],__mpz_struct[1]>

        """
        from sage.all import Sequence

        return cls.to_pyflatsurf(domain=Sequence(elements).universe(), codomain=codomain)

    def domain(self):
        r"""
        Return the domain of this conversion, a sage-flatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.domain()
            Translation Surface in H_2(2) built from 6 isosceles triangles

        """
        if self._domain is not None:
            return self._domain

        raise NotImplementedError(f"{type(self).__name} does not implement domain() yet")

    def codomain(self):
        r"""
        Return the codomain of this conversion, a pyflatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.codomain()
            FlatTriangulationCombinatorial(vertices = ...) with vectors {...}

        """
        if self._codomain is not None:
            return self._codomain

        raise NotImplementedError(f"{type(self).__name__} does not implement codomain() yet")

    def __call__(self, x):
        r"""
        Return the conversion at an element of :meth:`domain` and return the
        corresponding pyflatsurf object.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: conversion(p)
            ((-1/4*a^2 + 1/2*a + 1/2 ~ 0.54654802), (1/4*a^2 - 3/4 ~ 0.15450850), (1/4 ~ 0.25000000)) in (1, 2, 3)

        """
        raise NotImplementedError(f"{type(self).__name__} does not implement a mapping of elements yet")

    def section(self, y):
        r"""
        Return the conversion of an element of :meth:`codomain` and return the
        corresponding sage-flatsurf object.

        This is the inverse of :meth:`__call__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion.section(q)
            Point (0, 1/2) of polygon 0

        """
        raise NotImplementedError(f"{type(self).__name__} does not implement a section yet")

    def __repr__(self):
        r"""
        Return a printable representation of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: FlatTriangulationConversion.to_pyflatsurf(S)
            Conversion from Translation Surface in H_2(2) built from 6 isosceles triangles to FlatTriangulationCombinatorial(vertices = ...) with vectors ...

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
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: conversion = RingConversion.to_pyflatsurf_from_elements([1, 2, 3])

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: isinstance(conversion, RingConversion)
        True

    """
    @staticmethod
    def _ring_conversions():
        r"""
        Return the available ring conversion types.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: conversions = RingConversion._ring_conversions()
            sage: list(conversions)  # random output, depends on the installed packages
            [<class 'flatsurf.geometry.pyflatsurf_conversion.RingConversion_eantic'>]

        """
        from flatsurf.features import pyeantic_feature

        yield RingConversion_gmp

        yield RingConversion_int

        if pyeantic_feature.is_present():
            yield RingConversion_eantic

    @classmethod
    def _create_conversion(cls, domain=None, codomain=None):
        r"""
        Return a conversion from ``domain`` to ``codomain``.

        Return ``None`` if this conversion cannot handle this
        ``domain``/``codomain`` combination.

        At least one of ``domain`` and ``codomain`` must not be ``None``.

        INPUT:

        - ``domain`` -- a SageMath ring (default: ``None``)

        - ``codomain`` -- a ring that pyflatsurf understands (default: ``None``)

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_eantic
            sage: RingConversion_eantic._create_conversion(QuadraticField(2))
            Conversion from Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? to NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])

        """
        raise NotImplementedError(f"{cls.__name__} does not implement _create_conversion() yet")

    @classmethod
    def _deduce_codomain(cls, elements):
        r"""
        Given elements from pyflatsurf, deduce which pyflatsurf ring it lives
        in.

        Return ``None`` if no (single) conversion can handle these elements.

        INPUT:

        - ``elements`` -- any sequence of objects that pyflatsurf understands

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_eantic
            sage: conversion = RingConversion.to_pyflatsurf(QuadraticField(2))
            sage: element = conversion.codomain().gen()

            sage: RingConversion_eantic._deduce_codomain([element])
            NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])

        """
        raise NotImplementedError(f"{cls.__name__} does not implement _deduce_codomain() yet")

    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`RingConversion` that converts the SageMath ring ``domain``
        to something that libflatsurf/pyflatsurf can understand.

        INPUT:

        - ``domain`` -- a ring

        - ``codomain`` -- a C/C++ type or ``None`` (default: ``None``); if
          ``None``, the corresponding type is constructed.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: RingConversion.to_pyflatsurf(QQ)
            Conversion from Rational Field to __gmp_expr<__mpq_struct[1],__mpq_struct[1]>

        """
        for conversion_type in RingConversion._ring_conversions():
            conversion = conversion_type._create_conversion(domain=domain, codomain=codomain)
            if conversion is not None:
                return conversion

        raise NotImplementedError(f"cannot determine pyflatsurf ring corresponding to {domain} yet")

    @classmethod
    def from_pyflatsurf_from_flat_triangulation(cls, flat_triangulation, domain=None):
        r"""
        Return a :class:`RingConversion` that can map ``domain`` to the ring
        over which ``flat_triangulation`` is defined.

        INPUT:

        - ``flat_triangulation`` -- a libflatsurf ``FlatTriangulation``

        - ``domain`` -- a SageMath ring, or ``None`` (default: ``None``); if
          ``None``, the ring is determined automatically.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion, RingConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: flat_triangulation = FlatTriangulationConversion.to_pyflatsurf(S).codomain()
            sage: conversion = RingConversion.from_pyflatsurf_from_flat_triangulation(flat_triangulation)

        Note that this conversion does not roundtrip back to the same SageMath
        ring. An e-antic ring only has a single variable name but a SageMath
        number field has a variable name and a (potentially different) variable
        name in the defining polynomial::

            sage: conversion.domain() is S.base_ring()
            False
            sage: conversion.domain()
            Number Field in a with defining polynomial x^4 - 5*x^2 + 5 with a = 1.902113032590308?
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.902113032590308?

        We can explicitly specify the domain to create the same conversion again::

            sage: conversion = RingConversion.from_pyflatsurf_from_flat_triangulation(flat_triangulation, domain=S.base_ring())
            sage: conversion.domain() is S.base_ring()
            True

        """
        return cls.from_pyflatsurf_from_elements([flat_triangulation.fromHalfEdge(he).x() for he in flat_triangulation.halfEdges()] + [flat_triangulation.fromHalfEdge(he).y() for he in flat_triangulation.halfEdges()], domain=domain)

    @classmethod
    def from_pyflatsurf_from_elements(cls, elements, domain=None):
        r"""
        Return a :class:`RingConversion` that can map ``domain`` to the ring of
        ``elements``.

        INPUT:

        - ``elements`` -- a sequence of elements that libflatsurf/pyflatsurf understands

        - ``domain`` -- a SageMath ring, or ``None`` (default: ``None``); if
          ``None``, the ring is determined automatically.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion

            sage: import gmpxxyy
            sage: conversion = RingConversion.from_pyflatsurf_from_elements([gmpxxyy.mpz()])

        """
        for conversion_type in RingConversion._ring_conversions():
            codomain = conversion_type._deduce_codomain(elements)
            if codomain is None:
                continue
            conversion = conversion_type._create_conversion(domain=domain, codomain=codomain)
            if conversion is not None:
                return conversion

        raise NotImplementedError(f"cannot determine a SageMath ring for {elements} yet")

    @classmethod
    def from_pyflatsurf(cls, codomain, domain=None):
        r"""
        Return a :class:`RingConversion` that maps ``domain`` to ``codomain``.

        INPUT:

        - ``codomain`` -- a libflatsurf/pyflatsurf type or ring

        - ``domain`` -- a SageMath ring, or ``None`` (default: ``None``); if
          ``None``, the ring is determined automatically.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: conversion = RingConversion.to_pyflatsurf(QQ)

            sage: RingConversion.from_pyflatsurf(conversion.codomain())
            Conversion from Rational Field to __gmp_expr<__mpq_struct[1],__mpq_struct[1]>

        """
        for conversion_type in RingConversion._ring_conversions():
            conversion = conversion_type._create_conversion(domain=domain, codomain=codomain)
            if conversion is not None:
                return conversion

        raise NotImplementedError(f"cannot determine pyflatsurf ring corresponding to {codomain} yet")


class RingConversion_eantic(RingConversion):
    r"""
    A conversion from a SageMath number field to an e-antic real embedded number field.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: conversion = RingConversion.to_pyflatsurf(domain=QuadraticField(2))

        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_eantic
        sage: isinstance(conversion, RingConversion_eantic)
        True

    """

    @classmethod
    def _create_conversion(cls, domain=None, codomain=None):
        r"""
        Implements :meth:`RingConversion._create_conversion`.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_eantic
            sage: RingConversion_eantic._create_conversion(domain=QuadraticField(3))
            Conversion from Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878? to NumberField(a^2 - 3, [1.732050807568877293527446341505872367 +/- 2.90e-37])

        """
        if domain is None and codomain is None:
            raise ValueError("at least one of domain and codomain must be set")

        if domain is None:
            from pyeantic import RealEmbeddedNumberField
            renf = RealEmbeddedNumberField(codomain)
            domain = renf.number_field

        from sage.all import NumberFields, QQ, RR

        if domain not in NumberFields():
            return None
        if domain is QQ:
            # GMP should handle the rationals
            return None

        if not domain.embeddings(RR):
            raise NotImplementedError("cannot determine pyflatsurf ring for not real-embedded number fields yet")
        if not domain.is_absolute():
            raise NotImplementedError("cannot determine pyflatsurf ring for a relative number field since there are no relative fields in e-antic yet")

        if codomain is None:
            from pyeantic import RealEmbeddedNumberField
            renf = RealEmbeddedNumberField(domain)
            codomain = unwrap_intrusive_ptr(renf.renf)

        return RingConversion_eantic(domain, codomain)

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this conversion.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = QuadraticField(2)
            sage: conversion = RingConversion.to_pyflatsurf(domain)
            sage: conversion(domain.gen())
            (a ~ 1.4142136)

        """
        import sage.structure.element

        parent = sage.structure.element.parent(x)

        if parent is not self.domain():
            raise ValueError(f"argument must be in the domain of this conversion but {x} is in {parent} and not in {self.domain()}")

        return self._pyrenf()(list(x)).renf_elem

    @cached_method
    def _pyrenf(self):
        r"""
        Return the pyeantic ``RealEmbeddedNumberField`` that wraps the codomain
        of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = QuadraticField(2)
            sage: conversion = RingConversion.to_pyflatsurf(domain)
            sage: conversion._pyrenf()
            Real Embedded Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095?

        """
        from pyeantic import RealEmbeddedNumberField

        return RealEmbeddedNumberField(self.codomain())

    def section(self, y):
        r"""
        Return the preimage of ``y`` under this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = QuadraticField(2)
            sage: conversion = RingConversion.to_pyflatsurf(domain)
            sage: gen = conversion(domain.gen())
            sage: conversion.section(gen)
            a

        """
        return self.domain()(self._pyrenf()(y))

    @classmethod
    def _deduce_codomain(cls, elements):
        r"""
        Given elements from e-antic, deduce the e-antic ring they live in.

        This implements :meth:`RingConversion._deduce_codomain`.

        Return ``None`` if this is not an e-antic element or if no single
        e-antic ring can handle them.

        INPUT:

        - ``elements`` -- a sequence of e-antic elements or something else

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_eantic
            sage: conversion = RingConversion.to_pyflatsurf(QuadraticField(2))
            sage: element = conversion.codomain().gen()

            sage: RingConversion_eantic._deduce_codomain([element])
            NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])
            sage: RingConversion_eantic._deduce_codomain([element, 2*element])
            NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])

            sage: import cppyy
            sage: RingConversion_eantic._deduce_codomain([element, cppyy.gbl.eantic.renf_elem_class()])
            NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])

        """
        import pyeantic

        ring = None

        for element in elements:
            if not isinstance(element, pyeantic.eantic.renf_elem_class):
                return None

            element_ring = unwrap_intrusive_ptr(element.parent())
            if ring is None or ring.degree() == 1:
                ring = element_ring
            elif element_ring != ring and element_ring.degree() != 1:
                return None

        return ring


class RingConversion_int(RingConversion):
    r"""
    Conversion between SageMath integers and machine long long integers.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: conversion = RingConversion.to_pyflatsurf(domain=int)

        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_int
        sage: isinstance(conversion, RingConversion_int)
        True

    """

    @classmethod
    def _create_conversion(cls, domain=None, codomain=None):
        r"""
        Implements :meth:`RingConversion._create_conversion`.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_int
            sage: RingConversion_int._create_conversion(domain=int)
            Conversion from <class 'int'> to long long

        """
        if domain is None and codomain is None:
            raise ValueError("at least one of domain and codomain must be set")

        import cppyy

        longlong = getattr(cppyy.gbl, 'long long')

        if domain in [None, int] and codomain in [None, longlong]:
            return RingConversion_int(int, longlong)

        return None

    @classmethod
    def _deduce_codomain(cls, elements):
        r"""
        Given long longs, return the long long type.

        This implements :meth:`RingConversion._deduce_codomain`.

        Return ``None`` if these are not long long element.

        INPUT:

        - ``element`` -- a long long or something else

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_int
            sage: conversion = RingConversion.to_pyflatsurf(int)
            sage: element = conversion.codomain()(1R)

            sage: RingConversion_int._deduce_codomain([element])
            <class 'cppyy.gbl.long long'>

        """
        import cppyy

        if not elements:
            return None

        longlong = getattr(cppyy.gbl, 'long long')

        for element in elements:
            if not isinstance(element, longlong):
                return None

        return longlong

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this conversion.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = int
            sage: conversion = RingConversion.to_pyflatsurf(int)
            sage: conversion(1R)
            1

        """
        if not isinstance(x, int):
            raise ValueError("argument must be an int")

        return self.codomain()(x)

    def section(self, y):
        r"""
        Return the preimage of ``y`` under this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = int
            sage: conversion = RingConversion.to_pyflatsurf(domain)
            sage: y = conversion(3R)
            sage: conversion.section(y)
            3

        """
        return self.domain()(str(y))


class RingConversion_gmp(RingConversion):
    r"""
    Conversion between SageMath integers and rationals and GMP mpz and mpq.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: conversion = RingConversion.to_pyflatsurf(domain=ZZ)

        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_gmp
        sage: isinstance(conversion, RingConversion_gmp)
        True

    """

    @classmethod
    def _create_conversion(cls, domain=None, codomain=None):
        r"""
        Implements :meth:`RingConversion._create_conversion`.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion_gmp
            sage: RingConversion_gmp._create_conversion(domain=ZZ)
            Conversion from Integer Ring to __gmp_expr<__mpz_struct[1],__mpz_struct[1]>

        """
        if domain is None and codomain is None:
            raise ValueError("at least one of domain and codomain must be set")

        from sage.all import ZZ, QQ
        import cppyy

        if domain is None:
            if codomain == cppyy.gbl.mpz_class:
                domain = ZZ
            elif codomain is cppyy.gbl.mpq_class:
                domain = QQ
            else:
                return None

        if codomain is None:
            if domain is ZZ:
                codomain = cppyy.gbl.mpz_class
            elif domain is QQ:
                codomain = cppyy.gbl.mpq_class
            else:
                return None

        if domain is ZZ and codomain is cppyy.gbl.mpz_class:
            return RingConversion_gmp(domain, codomain)
        if domain is QQ and codomain is cppyy.gbl.mpq_class:
            return RingConversion_gmp(domain, codomain)

        return None

    @classmethod
    def _deduce_codomain(cls, elements):
        r"""
        Given elements from GMP, return the GMP ring they live in.

        This implements :meth:`RingConversion._deduce_codomain`.

        Return ``None`` if there not (compatible) GMP elements.

        INPUT:

        - ``elements`` -- a sequence of GMP elements or something else

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_gmp
            sage: conversion = RingConversion.to_pyflatsurf(ZZ)
            sage: element = conversion.codomain()(1)

            sage: RingConversion_gmp._deduce_codomain([element])
            <class cppyy.gbl.__gmp_expr<__mpz_struct[1],__mpz_struct[1]> at 0x...>

        """
        import cppyy

        ring = None
        for element in elements:
            if isinstance(element, cppyy.gbl.mpz_class):
                if ring is None or ring is cppyy.gbl.mpz_class:
                    ring = cppyy.gbl.mpz_class
                else:
                    return None

            if isinstance(element, cppyy.gbl.mpq_class):
                if ring is None or ring is cppyy.gbl.mpq_class:
                    ring = cppyy.gbl.mpq_class
                else:
                    return None

        return ring

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this conversion.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = QQ
            sage: conversion = RingConversion.to_pyflatsurf(QQ)
            sage: conversion(1 / 3)
            1/3

        """
        import sage.structure.element

        parent = sage.structure.element.parent(x)

        if parent is not self.domain():
            raise ValueError(f"argument must be in the domain of this conversion but {x} is in {parent} and not in {self.domain()}")

        return self.codomain()(str(x))

    def section(self, y):
        r"""
        Return the preimage of ``y`` under this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
            sage: domain = QQ
            sage: conversion = RingConversion.to_pyflatsurf(domain)
            sage: y = conversion(1/3)
            sage: conversion.section(y)
            1/3

        """
        return self.domain()(str(y))


class VectorSpaceConversion(Conversion):
    r"""
    Converts vectors in a SageMath vector space into libflatsurf ``Vector<T>``s.

    EXAMPLES::

        sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
        sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)

    """
    @classmethod
    def to_pyflatsurf(cls, domain, codomain=None):
        r"""
        Return a :class:`Conversion` from ``domain`` to ``codomain``.

        INPUT:

        - ``domain`` -- a SageMath ``VectorSpace``

        - ``codomain`` -- a libflatsurf ``Vector<T>`` type or ``None``
          (default: ``None``); if ``None``, the type is determined
          automatically.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)

        """
        pyflatsurf_feature.require()

        ring_conversion = RingConversion.to_pyflatsurf(domain.base_ring())

        if codomain is None:
            import pyflatsurf

            T = ring_conversion.codomain()
            if not isinstance(T, type):
                T = type(T)

            codomain = pyflatsurf.flatsurf.Vector[T]

        return VectorSpaceConversion(domain, codomain)

    @classmethod
    def from_pyflatsurf_from_elements(cls, elements, domain=None):
        r"""
        Return a :class:`Conversion` that converts the pyflatsurf ``elements``
        into vectors in ``domain``.

        INPUT:

        - ``domain`` -- a SageMath ``VectorSpace`` or ``None``; if ``None``, it
          is determined automatically.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: codomain = VectorSpaceConversion.to_pyflatsurf(QQ^2).codomain()

            sage: VectorSpaceConversion.from_pyflatsurf_from_elements([codomain()])
            Conversion from Vector space of dimension 2 over Rational Field to flatsurf::cppyy::Vector<__gmp_expr<__mpq_struct[1],__mpq_struct[1]> >

        """
        if domain is None:
            ring_conversion = RingConversion.from_pyflatsurf_from_elements([element.x() for element in elements] + [element.y() for element in elements])

            from sage.all import VectorSpace
            domain = VectorSpace(ring_conversion.domain(), 2)

        return VectorSpaceConversion.to_pyflatsurf(domain=domain)

    @cached_method
    def _vectors(self):
        r"""
        Return the pyflatsurf ``Vectors`` helper for the codomain of this conversion.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)
            sage: conversion._vectors()
            Flatsurf Vectors over Rational Field

        """
        from pyflatsurf.vector import Vectors

        return Vectors(self._ring_conversion().codomain())

    @cached_method
    def _ring_conversion(self):
        r"""
        Return the :class:`RingConversion` that maps the coordinates of the domain to the coordinates of the codomain.

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)
            sage: conversion._ring_conversion()
            Conversion from Rational Field to __gmp_expr<__mpq_struct[1],__mpq_struct[1]>

        """
        return RingConversion.to_pyflatsurf(self.domain().base_ring())

    def __call__(self, vector):
        r"""
        Return a ``Vector<T>`` corresponding to the SageMath ``vector``.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)
            sage: conversion(vector(QQ, [1, 2]))
            (1, 2)

        """
        ring_conversion = self._ring_conversion()
        return self._vectors()([ring_conversion(coordinate) for coordinate in vector]).vector

    def section(self, vector):
        r"""
        Return the SageMath vector corresponding to the pyflatsurf ``vector``.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
            sage: conversion = VectorSpaceConversion.to_pyflatsurf(QQ^2)
            sage: v = vector(QQ, [1, 2])
            sage: conversion.section(conversion(v)) == v
            True

        """
        ring_conversion = self._ring_conversion()
        return self.domain()([ring_conversion.section(vector.x()), ring_conversion.section(vector.y())])


class FlatTriangulationConversion(Conversion):
    r"""
    Converts a :class:`Surface` object to a ``FlatTriangulation`` object and
    vice-versa.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
        sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

    """

    def __init__(self, domain, codomain, label_to_half_edge):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.square_torus().triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: isinstance(conversion, FlatTriangulationConversion)
            True

        """
        super().__init__(domain=domain, codomain=codomain)

        self._label_to_half_edge = label_to_half_edge
        self._half_edge_to_label = {half_edge: label for (label, half_edge) in label_to_half_edge.items()}

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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

        """
        pyflatsurf_feature.require()

        from flatsurf.geometry.categories import TranslationSurfaces

        if domain not in TranslationSurfaces():
            raise TypeError("domain must be a translation surface")
        if not domain.is_finite_type():
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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: FlatTriangulationConversion._pyflatsurf_labels(S)
            {(0, 0): 1,
             (0, 1): 2,
             (0, 2): 3,
             (1, 0): 6,
             (1, 1): -2,
             (1, 2): -3,
             (2, 0): -4,
             (2, 1): 7,
             (2, 2): 8,
             (3, 0): -1,
             (3, 1): 4,
             (3, 2): 5,
             (4, 0): -9,
             (4, 1): -7,
             (4, 2): -8,
             (5, 0): -6,
             (5, 1): 9,
             (5, 2): -5}

        """
        labels = {}
        for half_edge, opposite_half_edge in domain.gluings():
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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: FlatTriangulationConversion._pyflatsurf_vectors(S)
            [((1/2 ~ 0.50000000), (1/2*a^3 - 1*a ~ 1.5388418)),
             ((-1/2*a^2 + 1 ~ -0.80901699), (-1/2*a^3 + 3/2*a ~ -0.58778525)),
             ((1/2*a^2 - 3/2 ~ 0.30901699), (-1/2*a ~ -0.95105652)),
             ((1/2*a^2 - 1/2 ~ 1.3090170), (1/2*a ~ 0.95105652)),
             ((-1/2*a^2 + 1 ~ -0.80901699), (1/2*a^3 - 3/2*a ~ 0.58778525)),
             ((-1/2 ~ -0.50000000), (-1/2*a^3 + 1*a ~ -1.5388418)),
             (1, 0),
             ((1/2*a^2 - 3/2 ~ 0.30901699), (1/2*a ~ 0.95105652)),
             ((-1/2*a^2 + 1/2 ~ -1.3090170), (-1/2*a ~ -0.95105652))]

        """
        labels = cls._pyflatsurf_labels(domain)

        vectors = [None] * (len(labels) // 2)

        for (polygon, edge), half_edge in labels.items():
            if half_edge < 0:
                continue

            vectors[half_edge - 1] = domain.polygon(polygon).edge(edge)

        vector_conversion = VectorSpaceConversion.to_pyflatsurf_from_elements(vectors)
        return [vector_conversion(vector) for vector in vectors]

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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: FlatTriangulationConversion._pyflatsurf_vertex_permutation(S)
            [[1, -3, 2, -1, -5, -9, 8, -7, 9, 6, 3, -2, -6, 5, -4, -8, 7, 4]]

        """
        pyflatsurf_labels = cls._pyflatsurf_labels(domain)

        vertex_permutation = {}

        for polygon, edge in domain.edges():
            pyflatsurf_edge = pyflatsurf_labels[(polygon, edge)]

            next_edge = (edge + 1) % len(domain.polygon(polygon).vertices())
            pyflatsurf_next_edge = pyflatsurf_labels[(polygon, next_edge)]

            vertex_permutation[pyflatsurf_next_edge] = -pyflatsurf_edge

        return cls._cycle_decomposition(vertex_permutation)

    @classmethod
    def _cycle_decomposition(self, permutation):
        r"""
        Return a permutation in cycle notation.

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: FlatTriangulationConversion._cycle_decomposition({})
            []
            sage: FlatTriangulationConversion._cycle_decomposition({1: 1, -1: -1})
            [[1], [-1]]
            sage: FlatTriangulationConversion._cycle_decomposition({1: 2, 2: 1, 3: 4, 4: 5, 5: 3})
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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: FlatTriangulationConversion.from_pyflatsurf(conversion.codomain())
            Conversion from ...

        """
        pyflatsurf_feature.require()

        if domain is None:
            ring_conversion = RingConversion.from_pyflatsurf_from_flat_triangulation(codomain)

            from flatsurf import MutableOrientedSimilaritySurface, Polygon
            domain = MutableOrientedSimilaritySurface(ring_conversion.domain())

            half_edge_to_polygon_edge = {}

            from sage.all import VectorSpace
            vector_conversion = VectorSpaceConversion.to_pyflatsurf(VectorSpace(ring_conversion.domain(), 2))

            for (a, b, c) in codomain.faces():
                vectors = [codomain.fromHalfEdge(a), codomain.fromHalfEdge(b), codomain.fromHalfEdge(c)]
                vectors = [vector_conversion.section(vector) for vector in vectors]
                triangle = Polygon(edges=vectors)

                label = (a.id(), b.id(), c.id())

                domain.add_polygon(triangle, label=label)

                half_edge_to_polygon_edge[a] = (label, 0)
                half_edge_to_polygon_edge[b] = (label, 1)
                half_edge_to_polygon_edge[c] = (label, 2)

            for half_edge, (polygon, edge) in half_edge_to_polygon_edge.items():
                opposite = half_edge_to_polygon_edge[-half_edge]
                domain.glue((polygon, edge), opposite)

            domain.set_immutable()

        return cls.to_pyflatsurf(domain=domain, codomain=codomain)

    @cached_method
    def ring_conversion(self):
        r"""
        Return the conversion that maps the base ring of the domain of this
        conversion to the base ring of the codomain of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.ring_conversion()
            Conversion from Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.902113032590308? to NumberField(a^4 - 5*a^2 + 5, [1.902113032590307144232878666758764287 +/- 6.87e-37])

        """
        return RingConversion.to_pyflatsurf(domain=self.domain().base_ring())

    @cached_method
    def vector_space_conversion(self):
        r"""
        Return the conversion maps two-dimensional vectors over the base ring
        of the domain to ``Vector<T>`` for the codomain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.vector_space_conversion()
            Conversion from Vector space of dimension 2 over Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.902113032590308? to flatsurf::cppyy::Vector<eantic::renf_class>

        """
        from sage.all import VectorSpace

        return VectorSpaceConversion.to_pyflatsurf(VectorSpace(self.ring_conversion().domain(), 2))

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this conversion.

        INPUT:

        - ``x`` -- an object defined on the domain, e.g., a
          :class:`SurfacePoint`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

        We map a point::

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: conversion(p)
            ((-1/4*a^2 + 1/2*a + 1/2 ~ 0.54654802), (1/4*a^2 - 3/4 ~ 0.15450850), (1/4 ~ 0.25000000)) in (1, 2, 3)

        We map a half edge::

            sage: conversion((0, 0))
            1

        """
        from flatsurf.geometry.surface_objects import SurfacePoint

        if isinstance(x, SurfacePoint):
            return self._image_point(x)
        if isinstance(x, tuple) and len(x) == 2:
            return self._image_half_edge(*x)

        raise NotImplementedError(f"cannot map {type(x)} from sage-flatsurf to pyflatsurf yet")

    def section(self, y):
        r"""
        Return the preimage of ``y`` under this conversion.

        INPUT:

        - ``y`` -- an object defined in the codomain, e.g., a pyflatsurf
          ``Point``

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

        We roundtrip a point::

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion.section(q) == p
            True

        We roundtrip a half edge::

            sage: half_edge = conversion((0, 0))
            sage: conversion.section(half_edge)
            (0, 0)

        """
        import pyflatsurf

        if isinstance(y, pyflatsurf.flatsurf.Point[type(self.codomain())]):
            return self._preimage_point(y)
        if isinstance(y, pyflatsurf.flatsurf.HalfEdge):
            return self._preimage_half_edge(y)

        raise NotImplementedError(f"cannot compute the preimage of a {type(y)} in sage-flatsurf yet")

    def _image_point(self, p):
        r"""
        Return the image of the :class:`SurfacePoint` ``p`` under this conversion.

        This is a helper method for :meth:`__call__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: conversion._image_point(p)
            ((-1/4*a^2 + 1/2*a + 1/2 ~ 0.54654802), (1/4*a^2 - 3/4 ~ 0.15450850), (1/4 ~ 0.25000000)) in (1, 2, 3)

        """
        if p.surface() is not self.domain():
            raise ValueError("point is not a point in the domain of this conversion")

        label = next(iter(p.labels()))
        coordinates = next(iter(p.coordinates(label)))

        import pyflatsurf
        return pyflatsurf.flatsurf.Point[type(self.codomain())](self.codomain(), self((label, 0)), self.vector_space_conversion()(coordinates))

    def _preimage_point(self, q):
        r"""
        Return the preimage of the point ``q`` in the domain of this conversion.

        This is a helper method for :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion._preimage_point(q)
            Point (0, 1/2) of polygon 0

        """
        face = q.face()
        label, edge = self.section(face)
        coordinates = self.vector_space_conversion().section(q.vector(face))
        while edge != 0:
            coordinates += self.domain().polygon(label).edge(edge - 1)
            edge -= 1

        coordinates += self.domain().polygon(label).vertex(0)

        from flatsurf.geometry.surface_objects import SurfacePoint
        return SurfacePoint(self.domain(), label, coordinates)

    def _image_half_edge(self, label, edge):
        r"""
        Return the half edge that ``edge`` of polygon ``label`` maps to under this conversion.

        This is a helper method for :meth:`__call__`.

        INPUT:

        - ``label`` -- an arbitrary polygon label in the :meth:`domain`

        - ``edge`` -- an integer, the identifier of an edge in the polygon ``label``

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: conversion._image_half_edge(0, 0)
            1

        """
        import pyflatsurf
        return pyflatsurf.flatsurf.HalfEdge(self._label_to_half_edge[(label, edge)])

    def _preimage_half_edge(self, half_edge):
        r"""
        Return the preimage of the ``half_edge`` in the domain of this
        conversion as a pair ``(label, edge)``.

        This is a helper method for :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().codomain()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: import pyflatsurf
            sage: conversion._preimage_half_edge(pyflatsurf.flatsurf.HalfEdge(1R))
            (0, 0)

        """
        return self._half_edge_to_label[half_edge.id()]


def to_pyflatsurf(S):
    r"""
    Given S a translation surface from sage-flatsurf return a
    flatsurf::FlatTriangulation from libflatsurf/pyflatsurf.
    """
    import warnings
    warnings.warn("to_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use FlatTriangulationConversion.to_pyflatsurf(surface.triangulate()).codomain() instead.")

    return FlatTriangulationConversion.to_pyflatsurf(S.triangulate()).codomain()


def sage_ring(surface):
    r"""
    Return the SageMath ring over which the pyflatsurf surface ``surface`` can
    be constructed in sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf, sage_ring # optional: pyflatsurf
        sage: S = to_pyflatsurf(translation_surfaces.veech_double_n_gon(5)) # optional: pyflatsurf  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: sage_ring(S) # optional: pyflatsurf
        doctest:warning
        ...
        UserWarning: to_sage_ring() is deprecated and will be removed in a future version of sage-flaturf. Use RingConversion.from_pyflatsurf_from_elements([x]).section(x) instead.
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
        sage: import cppyy  # optional: pyflatsurf
        sage: from sage.structure.element import parent  # optional: pyflatsurf
        sage: parent(to_sage_ring(getattr(cppyy.gbl, 'long long')(1R)))  # optional: pyflatsurf
        <class 'int'>

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
    warnings.warn("to_sage_ring() is deprecated and will be removed in a future version of sage-flaturf. Use RingConversion.from_pyflatsurf_from_elements([x]).section(x) instead.")

    return RingConversion.from_pyflatsurf_from_elements([x]).section(x)

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
        raise NotImplementedError
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
        doctest:warning
        ...
        UserWarning: from_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use TranslationSurface(FlatTriangulationConversion.from_pyflatsurf(surface).domain()) instead.
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
    import warnings
    warnings.warn("from_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use TranslationSurface(FlatTriangulationConversion.from_pyflatsurf(surface).domain()) instead.")

    return FlatTriangulationConversion.from_pyflatsurf(T).domain()
