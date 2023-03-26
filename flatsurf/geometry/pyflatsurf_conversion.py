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
    sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: FlatTriangulationConversion.from_pyflatsurf(conversion.codomain())
            <flatsurf.geometry.surface.Surface_dict object at 0x...>

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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.domain()
            <flatsurf.geometry.surface.Surface_dict object at 0x...>

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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.codomain()
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, -8, 9, 7, -6, -9, 4, 3, -2, -4, 8, 5, -7, 6, -5), faces = (1, 2, 3)(-1, -5, 8)(-2, -3, 4)(-4, -9, -8)(5, 6, 7)(-6, -7, 9)) with vectors {...}

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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion.section(q)
            Surface point located at (0, 1/2) in polygon 0

        """
        raise NotImplementedError(f"{type(self).__name__} does not implement a section yet")

    def __repr__(self):
        r"""
        Return a printable representation of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConversion.to_pyflatsurf(S)
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
    def _deduce_codomain(cls, element):
        r"""
        Given an element from pyflatsurf, deduce which pyflatsurf ring it lives
        in.

        Return ``None`` if this conversion cannot handle such an element.

        INPUT:

        - ``element`` -- any object that pyflatsurf understands

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_eantic
            sage: conversion = RingConversion.to_pyflatsurf(QuadraticField(2))
            sage: element = conversion.codomain().gen()

            sage: RingConversion_eantic._deduce_codomain(element)
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
        return cls.from_pyflatsurf_from_element(flat_triangulation.fromHalfEdge(1).x(), domain=domain)

    @classmethod
    def from_pyflatsurf_from_element(cls, element, domain=None):
        r"""
        Return a :class:`RingConversion` that can map ``domain`` to the ring of
        ``element``.

        INPUT:

        - ``element`` -- an element that libflatsurf/pyflatsurf understands

        - ``domain`` -- a SageMath ring, or ``None`` (default: ``None``); if
          ``None``, the ring is determined automatically.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion

            sage: import gmpxxyy
            sage: conversion = RingConversion.from_pyflatsurf_from_element(gmpxxyy.mpz())

        """
        for conversion_type in RingConversion._ring_conversions():
            codomain = conversion_type._deduce_codomain(element)
            if codomain is None:
                continue
            conversion = conversion_type._create_conversion(domain=domain, codomain=codomain)
            if conversion is not None:
                return conversion

        raise NotImplementedError(f"cannot determine a SageMath ring for a {type(element)} yet")

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
        sage: conversion = RingConversion.to_pyflatsurf(domain = QuadraticField(2))

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

            sage: from flatsurf import translation_surfaces
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

        from pyeantic import RealEmbeddedNumberField

        # TODO: Is this very slow? I guess we should cache this one.
        pyrenf = RealEmbeddedNumberField(self.codomain())
        return pyrenf(list(x)).renf_elem

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
        from pyeantic import RealEmbeddedNumberField

        pyrenf = RealEmbeddedNumberField(self.codomain())
        x = self.domain()(pyrenf(y))
        return x

    @classmethod
    def _deduce_codomain(cls, element):
        r"""
        Given an element from e-antic, deduce the e-antic ring it lives in.

        This implements :meth:`RingConversion._deduce_codomain`.

        Return ``None`` if this is not an e-antic element.

        INPUT:

        - ``element`` -- an e-antic element or something else

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_eantic
            sage: conversion = RingConversion.to_pyflatsurf(QuadraticField(2))
            sage: element = conversion.codomain().gen()

            sage: RingConversion_eantic._deduce_codomain(element)
            NumberField(a^2 - 2, [1.4142135623730950488016887242096980786 +/- 7.96e-38])

        """
        import pyeantic

        if not isinstance(element, pyeantic.eantic.renf_elem_class):
            return None

        return unwrap_intrusive_ptr(element.parent())


class RingConversion_gmp(RingConversion):
    r"""
    Conversion between SageMath integers and rationals and GMP mpz and mpq.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion
        sage: conversion = RingConversion.to_pyflatsurf(domain = ZZ)

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
    def _deduce_codomain(cls, element):
        r"""
        Given an element from GMP, return the GMP ring it lives in.

        This implements :meth:`RingConversion._deduce_codomain`.

        Return ``None`` if this is not a GMP element.

        INPUT:

        - ``element`` -- a GMP element or something else

        EXAMPLES::

            sage: from flatsurf.geometry.pyflatsurf_conversion import RingConversion, RingConversion_gmp
            sage: conversion = RingConversion.to_pyflatsurf(ZZ)
            sage: element = conversion.codomain()(1)

            sage: RingConversion_gmp._deduce_codomain(element)
            <class cppyy.gbl.__gmp_expr<__mpz_struct[1],__mpz_struct[1]> at 0x...>

        """
        import cppyy

        if isinstance(element, cppyy.gbl.mpz_class):
            return cppyy.gbl.mpz_class

        if isinstance(element, cppyy.gbl.mpq_class):
            return cppyy.gbl.mpq_class

        return None


class FlatTriangulationConversion(Conversion):
    r"""
    Converts a :class:`Surface` object to a ``FlatTriangulation`` object and
    vice-versa.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
        sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
        sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

    """

    def __init__(self, domain, codomain, label_to_half_edge):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.square_torus().triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConversion._pyflatsurf_labels(S)
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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConversion._pyflatsurf_vectors(S)
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

        ring_conversion = RingConversion.to_pyflatsurf_from_elements([vector[0] for vector in vectors] + [vector[1] for vector in vectors])

        from pyflatsurf.vector import Vectors
        vector_space = Vectors(ring_conversion.codomain())

        # TODO: Use _image_vector()
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
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: FlatTriangulationConversion._pyflatsurf_vertex_permutation(S)
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: FlatTriangulationConversion.from_pyflatsurf(conversion.codomain())
            <flatsurf.geometry.surface.Surface_dict object at 0x...>

        """
        pyflatsurf_feature.require()

        if domain is not None:
            return cls.to_pyflatsurf(domain=domain, codomain=codomain)

        ring_conversion = RingConversion.from_pyflatsurf_from_flat_triangulation(codomain)

        from flatsurf.geometry.surface import Surface_dict
        domain = Surface_dict(ring_conversion.domain())

        from flatsurf.geometry.polygon import ConvexPolygons
        polygons = ConvexPolygons(ring_conversion.domain())

        half_edge_to_polygon_edge = {}

        for (a, b, c) in codomain.faces():
            # TODO: Use _preimage_vector()
            vectors = [codomain.fromHalfEdge(a), codomain.fromHalfEdge(b), codomain.fromHalfEdge(c)]
            vectors = [(ring_conversion.section(vector.x()), ring_conversion.section(vector.y())) for vector in vectors]
            triangle = polygons(vectors)

            label = (a.id(), b.id(), c.id())

            domain.add_polygon(triangle, label=label)

            half_edge_to_polygon_edge[a] = (label, 0)
            half_edge_to_polygon_edge[b] = (label, 1)
            half_edge_to_polygon_edge[c] = (label, 2)

        for half_edge, (polygon, edge) in half_edge_to_polygon_edge.items():
            opposite = half_edge_to_polygon_edge[-half_edge]
            domain.change_edge_gluing(polygon, edge, *opposite)

        domain.set_immutable()

        return domain

    @cached_method
    def ring_conversion(self):
        r"""
        Return the conversion that maps the base ring of the domain of this
        conversion to the base ring of the codomain of this conversion.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)
            sage: conversion.ring_conversion()
            Conversion from Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.902113032590308? to NumberField(a^4 - 5*a^2 + 5, [1.902113032590307144232878666758764287 +/- 6.87e-37])

        """
        return RingConversion.to_pyflatsurf(domain=self.domain().base_ring())

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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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

    def _image_vector(self, vector):
        r"""
        Return the SageMath ``vector`` over the base ring of ``domain`` as a
        libflatsurf ``Vector`` over the base ring of ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: v = vector(conversion.domain().base_ring(), (1, 2))
            sage: conversion._image_vector(v)
            (1, 2)

        """
        ring_conversion = self.ring_conversion()

        from pyflatsurf.vector import Vectors

        vector_space = Vectors(ring_conversion.codomain())
        return vector_space(list(vector)).vector

    def _preimage_vector(self, vector):
        r"""
        Return the SageMath vector corresponding to the pyflatsurf ``vector``.

        This is inverse to :meth:`_image_vector`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: v = vector(conversion.domain().base_ring(), (1, 2))
            sage: conversion._preimage_vector(conversion._image_vector(v)) == v
            True

        """
        ring_conversion = self.ring_conversion()

        import sage.all

        return sage.all.vector(ring_conversion.domain(), (ring_conversion.section(vector.x()), ring_conversion.section(vector.y())))

    def _image_point(self, p):
        r"""
        Return the image of the :class:`SurfacePoint` ``p`` under this conversion.

        This is a helper method for :meth:`__call__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
        return pyflatsurf.flatsurf.Point[type(self.codomain())](self.codomain(), self((label, 0)), self._image_vector(coordinates))

    def _preimage_point(self, q):
        r"""
        Return the preimage of the point ``q`` in the domain of this conversion.

        This is a helper method for :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
            sage: conversion = FlatTriangulationConversion.to_pyflatsurf(S)

            sage: p = SurfacePoint(S, 0, (0, 1/2))
            sage: q = conversion(p)
            sage: conversion._preimage_point(q)
            Surface point located at (0, 1/2) in polygon 0

        """
        face = q.face()
        label, edge = self.section(face)
        coordinates = self._preimage_vector(q.vector(face))
        while edge != 0:
            coordinates -= self.domain().edge(label, edge - 1)
            edge -= 1

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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
            sage: S = translation_surfaces.veech_double_n_gon(5).triangulate().underlying_surface()
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
    warnings.warn("to_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use FlatTriangulationConversion.to_pyflatsurf(surface.triangulate().underlying_surface()).codomain() instead.")

    return FlatTriangulationConversion.to_pyflatsurf(S).codomain()


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
    warnings.warn("to_sage_ring() is deprecated and will be removed in a future version of sage-flaturf. Use RingConversion.from_pyflatsurf_element(x).section(x) instead.")

    return RingConversion.from_pyflatsurf_element(x).section(x)

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
    warnings.warn("from_pyflatsurf() is deprecated and will be removed in a future version of sage-flatsurf. Use TranslationSurface(FlatTriangulationConversion.from_pyflatsurf(surface).domain()) instead.")

    from flatsurf.geometry.translation_surface import TranslationSurface
    return TranslationSurface(FlatTriangulationConversion.from_pyflatsurf(T).domain())
