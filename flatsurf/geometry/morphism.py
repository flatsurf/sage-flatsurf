r"""
Morphisms between Surfaces

.. NOTE::

    One should think of these as maps between surfaces. However, they might
    often not be maps on the points of the surface but just on homology, or in
    some cases, they might not be meaningful as maps in any sense.

    Technically, this means that these are "morphisms in the category of
    objects" where being a morphism does not really have any mathematical
    meaning.

.. NOTE::

    Our morphism infrastructure contains quite a few workarounds to make the
    SageMath machinery work. The fundamental problem that we are facing is that
    our parents (surfaces) are not unique representations. However, surfaces do
    implement equality if they are indistinguishable (and this is a good idea
    to make pickling work.) SageMath has the assumption that if S == T (which
    in SageMath normally implies S is T) that then Hom(S) is Hom(T). We could
    implement this, but then Hom(T).domain() is not T but S. Instead, we opted
    for tricking the coercion machinery into allowing our non-unique homsets.
    Namely, when comparing morphisms, we do not coerce them into a common
    parent first but compare them directly.

EXAMPLES:

We can use morphisms to follow a surface through a retriangulation process::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon()
    sage: morphism = S.subdivide_edges(3)
    sage: morphism = morphism.codomain().subdivide() * morphism
    sage: T = morphism.codomain()

    sage: morphism
    Composite morphism:
      From: Translation Surface in H_2(2) built from a regular octagon
      To:   Translation Surface in H_2(2, 0^9) built from 8 isosceles triangles and 16 triangles
      Defn:   Edge-Subdivision morphism:
              From: Translation Surface in H_2(2) built from a regular octagon
              To:   Translation Surface in H_2(2, 0^8) built from a regular octagon with 16 marked vertices
            then
              Marked-Point-Insertion morphism:
              From: Translation Surface in H_2(2, 0^8) built from a regular octagon with 16 marked vertices
              To:   Translation Surface in H_2(2, 0^9) built from 8 isosceles triangles and 16 triangles

We can then map points through the morphism::

    sage: p = S(0, (0, 0))
    sage: p
    Vertex 0 of polygon 0

    sage: q = morphism(p)
    sage: q
    Vertex 0 of polygon (0, 0)

A non-singular point::

    sage: p = S(0, (1, 1))
    sage: p
    Point (1, 1) of polygon 0

    sage: q = morphism(p)
    sage: q
    Point (1, 1) of polygon (0, 5)

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023-2024 Julian Rüth
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

from sage.categories.homset import Homset
from sage.misc.cachefunc import cached_method
from sage.categories.morphism import Morphism, IdentityMorphism as IdentityMorphism_sage
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import Ring
from flatsurf.geometry.surface import OrientedSimilaritySurface


class UnknownRing(UniqueRepresentation, Ring):
    r"""
    A placeholder for a SageMath ring that has been lost in the process of
    creating a morphism.

    Ideally, every morphism has a domain and a codomain. However, when
    migrating code that did not produce a morphism originally, it can be
    complicated to get a hold of the actual domain/codomain of a morphism.

    Instead, it is often convenient to set the domain/codomain to ``None``,
    i.e., stating that the domain/codomain is unknown. When this happens, also
    the ring over which the domain/codomain is defined is technically unknown.
    We use this placeholder ring in these situations since a surface requires a
    ring it is defined over.

    This ring provides no functionality other than being in the category of
    rings.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
        sage: S = translation_surfaces.square_torus()
        sage: S = MutableOrientedSimilaritySurface.from_surface(S)

        sage: triangulation = S.triangulate(in_place=True)
        sage: triangulation.domain()
        Unknown Surface
        sage: triangulation.domain().base_ring()
        The Unknown Ring

    TESTS::

        sage: from flatsurf.geometry.morphism import UnknownRing
        sage: isinstance(triangulation.domain().base_ring(), UnknownRing)
        True
        sage: isinstance(triangulation.codomain().base_ring(), UnknownRing)
        False

        sage: TestSuite(triangulation.codomain().base_ring()).run()

    """

    def __init__(self):
        from sage.all import ZZ

        super().__init__(ZZ)

    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)

            sage: triangulation = S.triangulate(in_place=True)
            sage: triangulation.domain().base_ring()
            The Unknown Ring

        """
        return "The Unknown Ring"


class UnknownSurface(UniqueRepresentation, OrientedSimilaritySurface):
    r"""
    A placeholder surface for a morphism's domain or codomain when that
    codomain is unknown or mutable.

    In SageMath a morphism must have an explicit domain and codomain. However,
    the domain of a morphism might not actually be known, for example when
    deforming a mutable surface, or exposing it might break things when it is
    mutable.

    In such cases, we replace the domain with this unknown surface which has no
    functionality other than being a surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
        sage: S = translation_surfaces.square_torus()
        sage: S = MutableOrientedSimilaritySurface.from_surface(S)
        sage: S
        Translation Surface built from a square

        sage: triangulation = S.triangulate(in_place=True)
        sage: S
        Translation Surface built from 2 isosceles triangles

        sage: triangulation.domain()
        Unknown Surface
        sage: triangulation.codomain()
        Unknown Surface

    TESTS::

        sage: from flatsurf.geometry.morphism import UnknownSurface
        sage: isinstance(triangulation.domain(), UnknownSurface)
        True

        sage: TestSuite(triangulation.domain()).run()

    For pragmatic reasons, all unknown surfaces are equal. Mathematically, this
    is not correct but otherwise pickling breaks and we need a lot of special
    casing for the unknown surfaces everywhere::

        sage: triangulation.domain() is triangulation.codomain()
        True

    """

    def is_mutable(self):
        r"""
        Return whether this surface can be mutated.

        Since nothing about this surface can be changed, we return that it is immutable.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.is_mutable()
            False

        """
        return False

    def _an_element_(self):
        r"""
        Return a point on this surface for testing.

        Since this surface cannot represent any points, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S._an_element_()

        """
        raise NotImplementedError("cannot produce points in an unknown surface")

    def roots(self):
        r"""
        Return the root labels of the connected components of this surface.

        Since this surface doesn't really have any labels, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.roots()

        """
        raise NotImplementedError("cannot determine root labels in an unknown surface")

    def is_finite_type(self):
        r"""
        Return whether this surface is built from a finite number of polygons.

        Since we do not record whether this unknown surface stands in for a
        surface of finite type, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.is_finite_type()

        """
        raise NotImplementedError(
            "cannot determine whether an unknown surface is of finite type"
        )

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        Since we do not record whether this unknown surface stands in for a
        surface that is compact, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.is_compact()

        """
        raise NotImplementedError(
            "cannot determine whether an unknown surface is compact"
        )

    def is_with_boundary(self):
        r"""
        Return whether this surface has a boundary of unglued edges.

        Since we do not record whether this unknown surface stands in for a
        surface that has a boundary, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.is_with_boundary()

        """
        raise NotImplementedError(
            "cannot determine whether an unknown surface has boundary"
        )

    def opposite_edge(self, label, edge):
        r"""
        Return the edge that is glued to the ``edge`` of the polygon with
        ``label``.

        Since we do not know anything about the gluings of this surface, we
        throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.opposite_edge(0, 0)

        """
        raise NotImplementedError("cannot determine how the unknown surface is glued")

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        Since we do not know the polygons that make up this surface, we throw an error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S.polygon(0)

        """
        raise NotImplementedError("cannot determine polygons of the unknown surface")

    def _test_an_element(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_components(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_elements(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_elements_eq_reflexive(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_elements_eq_symmetric(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_elements_eq_transitive(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_elements_neq(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_gluings(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_labels_polygons(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_refined_category(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _test_some_elements(self, **options):
        # This generic tests does not make sense on the unknown surface and is herefore disabled.
        pass

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: S = S.triangulate(in_place=False).domain()
            sage: S

        """
        return "Unknown Surface"


class MorphismSpace(Homset):
    r"""
    A set of morphisms between structures attached to surfaces.

    This is a base class for surface morphisms but also for other structures
    such as morphisms in homology.

    .. NOTE::

        Since surfaces are not unique parents, we need to override some
        functionallity of the SageMath Homset here to make pickling work
        correctly.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
        sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
        sage: triangulation = S.triangulate()

        sage: homset = triangulate.parent()
        sage: homset

    TESTS::

        sage: from flatsurf.geometry.morphism import SurfaceMorphismSpace
        sage: isinstance(homset, SurfaceMorphismSpace)
        True

        sage: TestSuite(homset).run()

    """

    def __init__(self, domain, codomain, category=None, base=None):
        from sage.categories.all import Objects

        self._category = category or Objects()
        super().__init__(domain, codomain, category=self._category, base=base)

    def base_ring(self):
        r"""
        Return the base ring of this morphism space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: S.triangulate().base_ring()

        """
        base_ring = super().base_ring()

        if base_ring is None:
            import warnings

            warnings.warn(
                f"This morphism set has no base ring. Are you trying to get the base ring of a surface? Use .codomain().base_ring() instead."
            )
            return self.codomain().base_ring()

        return base_ring

    def _an_element_(self):
        r"""
        Return a morphism in this set (for testing.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: End(S).an_element()

        """
        if self.is_endomorphism_set():
            return self.identity()
        raise NotImplementedError(
            f"cannot create a morphism from {self.domain()} to {self.codomain()} yet"
        )

    def identity(self):
        r"""
        Return the identical morphism from the domain to the codomain of this
        space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: End(S).identity()

        """
        raise NotImplementedError("this morphism space does not implement an identity morphism yet")

    def _test_associativity(self, **options):
        # We would fail tests for associativity because we do not actually
        # decide whether two morphisms are "equal" but just whether they are
        # indistinguishable. Consequently, identity * identity != identity. We
        # disable tests that fails because of this.
        pass

    def _test_one(self, **options):
        # We would fail this test because we do not actually decide whether two
        # morphisms are "equal" but just whether they are indistinguishable.
        # Consequently, identity * identity != identity. We disable tests that
        # fails because of this.
        pass

    # TODO: Explain why this fails.
    # def _test_prod(self, **options):
    #     pass

    def __reduce__(self):
        r"""
        Return a serializable version of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: loads(dumps(End(S))) == End(S)
            True

        """
        raise NotImplementedError("this morphism space does not implement pickling yet")

    def __eq__(self, other):
        r"""
        Return whether this space is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: End(S) == End(S)
            True

        """
        if type(self) is not type(other):
            return False
        return (
            self.domain() == other.domain()
            and self.codomain() == other.codomain()
            and self.category() == other.category()
        )

    def __hash__(self):
        r"""
        Return a hash value for this space that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: hash(End(S)) == hash(End(S))
            True

        """
        return hash((self.domain(), self.codomain()))


class SurfaceMorphismSpace(MorphismSpace):
    r"""
    The set of morphisms from surface ``domain`` to surface ``codomain``.

    .. NOTE::

        Since surfaces are not unique parents, we need to override some
        functionallity of the SageMath Homset here to make pickling work
        correctly.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
        sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
        sage: identity = S.erase_marked_points()

        sage: homset = identity.parent()
        sage: homset
        Surface Endomorphisms of Translation Surface in H_2(2) built from 3 squares

    TESTS::

        sage: from flatsurf.geometry.morphism import SurfaceMorphismSpace
        sage: isinstance(homset, SurfaceMorphismSpace)
        True

        sage: TestSuite(homset).run()

    """

    def identity(self):
        r"""
        Return the identical morphism from the domain to the codomain of this
        space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: End(S).identity()

        """
        if self.is_endomorphism_set():
            return IdentityMorphism._create_morphism(self.domain())
        return super().identity()

    def __repr__(self):
        r"""
        Return a printable representation of this space.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: End(S)

        """
        if self.domain() is self.codomain():
            return f"Surface Endomorphisms of {self.domain()}"
        return f"Surface Morphisms from {self.domain()} to {self.codomain()}"

    def __reduce__(self):
        r"""
        Return a serializable version of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: loads(dumps(End(S))) == End(S)
            True

        """
        return SurfaceMorphismSpace, (self.domain(), self.codomain(), self._category)


class SurfaceMorphism(Morphism):
    r"""
    Abstract base class for all morphisms that map from a ``domain`` surface to
    a ``codomain`` surface.

    INPUT:

    - ``domain`` -- a surface or ``None``; if ``None`` (or if the domain is
      mutable), the domain is replaced with the :class:`UnknownSurface` which
      means that almost all queries related to the domain are going to fail.

    - ``codomain`` -- a surface or ``None``; if ``None`` (or if the codomain is
      mutable), the codomain is replaced with the :class:`UnknownSurface` which
      means that almost all queries related to the codomain are going to fail.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
        sage: morphism.domain()
        Translation Surface in H_1(0) built from a square
        sage: morphism.codomain()
        Translation Surface in H_1(0) built from a rectangle

    TESTS::

        sage: from flatsurf.geometry.morphism import SurfaceMorphism
        sage: isinstance(morphism, SurfaceMorphism)
        True

    """

    def __init__(self, parent, category=None):
        if category is None:
            from sage.categories.all import Objects

            category = Objects()

        super().__init__(parent)

    @classmethod
    def _parent(cls, domain, codomain):
        r"""
        Return the homset containing the surface morphisms from ``domain`` to
        ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: from flatsurf.geometry.morphism import SurfaceMorphism
            sage: SurfaceMorphism._parent(S, S)
            Surface Endomorphisms of Translation Surface in H_1(0) built from a square

        """
        if domain is None:
            domain = UnknownSurface(UnknownRing())
        elif domain.is_mutable():
            domain = UnknownSurface(domain.base_ring())

        if codomain is None:
            codomain = UnknownSurface(UnknownRing())
        elif codomain.is_mutable():
            codomain = UnknownSurface(codomain.base_ring())

        return SurfaceMorphismSpace(domain, codomain)

    @classmethod
    def _create_morphism(cls, domain, codomain, *args, **kwargs):
        r"""
        Return a morphism of this type from ``domain`` to ``codomain``.

        Any additional parameters are passed on the to the constructor of this
        type.

        .. NOTE::

            All morphisms must be created through this method so that they have
            their parent and category set correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: from flatsurf.geometry.morphism import IdentityMorphism
            sage: IdentityMorphism._create_morphism(S)
            Identity endomorphism of Translation Surface in H_1(0) built from a square

        """
        parent = cls._parent(domain, codomain)
        return parent.__make_element_class__(cls)(parent, *args, **kwargs)

    def __getattr__(self, name):
        r"""
        Redirect attribute lookup to the codomain.

        EXAMPLES::

        A lot of methods that used to return a surface now return a morphism.
        To make transition of existing code easier, we look up attributes that
        cannot be found on the morphism up on the codomain and issue a
        deprecation warning::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism.is_triangulated()
            doctest:warning
            ...
            UserWarning: This methods returns a morphism instead of a surface. Use .codomain().is_triangulated to access the surface instead of the morphism.
            False

            sage: morphism.codomain().is_triangulated()
            False

            sage: morphism.is_foo()
            Traceback (most recent call last):
            ...
            AttributeError: ... has no attribute 'is_foo'

        """
        if name in ["__cached_methods", "_cached_methods"]:
            # Do not redirect __cached_methods to the surface since we do not
            # want to get the morphism and the surface cache mixed up.
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        try:
            attr = getattr(self.codomain(), name)
        except AttributeError:
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        import warnings

        warnings.warn(
            f"This methods returns a morphism instead of a surface. Use .codomain().{name} to access the surface instead of the morphism."
        )

        return attr

    def section(self):
        r"""
        Return a section of this morphism from its codomain to its domain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism.section()
            Linear morphism:
              From: Translation Surface in H_1(0) built from a rectangle
              To:   Translation Surface in H_1(0) built from a square
              Defn: [1/2   0]
                    [  0   1]

        """
        return SectionMorphism._create_morphism(self)

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.
        """
        return type(self).__name__

    def change(self, domain=None, codomain=None, check=True):
        r"""
        Return a copy of this morphism with the domain or codomain replaced
        with ``domain`` and ``codomain``, respectively.

        For this to work, the ``domain`` must be trivially a replacement for
        the original domain and the ``codomain`` must be trivially a
        replacement for the original codomain. This method is sometimes useful
        to implement new morphisms. It should not be necessary to call this
        method otherwise. This method is usually used when the domain or
        codomain was originally mutable or to replace the domain or codomain
        with another indistinguishable domain or codomain.

        INPUT:

        - ``domain`` -- a surface (default: ``None``); if set, the surfaces
          replaces the domain of this morphism

        - ``codomain`` -- a surface (default: ``None``); if set, the surfaces
          replaces the codomain of this morphism

        - ``check`` -- a boolean (default: ``True``); whether to check
          compatibility of the ``domain`` and ``codomain`` with the data
          defining the original morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism.domain()
            Unknown Surface

            sage: S.set_immutable()
            sage: morphism = morphism.change(domain=S)
            sage: morphism.domain()
            Translation Surface in H_1(0) built from a square

        ::

            sage: from flatsurf import translation_surfaces, MutableOrientedSimilaritySurface
            sage: S = translation_surfaces.square_torus()
            sage: S = MutableOrientedSimilaritySurface.from_surface(S)
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=True)
            sage: morphism.domain()
            Unknown Surface
            sage: morphism.codomain()
            Unknown Surface

            sage: S.set_immutable()
            sage: morphism = morphism.change(codomain=S)
            sage: morphism.domain()
            Unknown Surface
            sage: morphism.codomain()
            Translation Surface in H_1(0) built from a rectangle

        """
        if domain is not None:
            raise NotImplementedError(
                f"a {type(self).__name__} cannot swap out its domain yet"
            )
        if codomain is not None:
            raise NotImplementedError(
                f"a {type(self).__name__} cannot swap out its codomain yet"
            )

        return self

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- a point of the domain of this morphism or an object defined
          on the domain such as a homology class, a saddle connection or a
          cohomology class. (Mapping of most of these inputs might not be
          implemented, either because it makes no real mathematical sense or
          because it has simply not been implemented yet.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

        The image of a point::

            sage: morphism(S(0, (1/2, 0)))
            Point (1, 0) of polygon 0

            sage: morphism(_)
            Traceback (most recent call last):
            ...
            ValueError: point must be in the domain of this morphism

        The image of a saddle connection::

            sage: saddle_connection = next(iter(S.saddle_connections(1)))
            sage: saddle_connection
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 2 of polygon 0
            sage: morphism(saddle_connection)
            Saddle connection (2, 0) from vertex 0 of polygon 0 to vertex 2 of polygon 0

            sage: morphism(_)
            Traceback (most recent call last):
            ...
            ValueError: saddle connection must be in the domain of this morphism

        The image of a homology class::

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: a, b = H.gens()
            sage: a
            B[(0, 1)]
            sage: morphism(a)
            B[(0, 1)]

            sage: morphism(_)
            Traceback (most recent call last):
            ...
            ValueError: homology class must be defined over the domain of this morphism

        The image of a tangent vector::

            sage: t = S.tangent_vector(0, (0, 0), (1, 1))
            sage: morphism(t)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, 1)

        TESTS::

            sage: morphism(42)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot map Integer yet through ...

        """
        from flatsurf.geometry.surface_objects import SurfacePoint

        if isinstance(x, SurfacePoint):
            if x.parent() is not self.domain():
                raise ValueError("point must be in the domain of this morphism")
            image = self._image_point(x)
            assert image.parent() is self.codomain()
            return image

        raise NotImplementedError(f"cannot map {type(x).__name__} through {self} yet")

    def _image_point(self, p):
        r"""
        Return the image of the surface point ``p`` under this morphism.

        This is a helper method for :meth:`__call__`.

        Subclasses should implement this method if the morphism is meaningful
        on the level of points.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

        The image of a point::

            sage: morphism(S(0, (1/2, 0)))  # indirect doctest
            Point (1, 0) of polygon 0

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute the image of a point yet"
        )

    def _section_point(self, q):
        r"""
        Return a preimage of the surface point ``q`` under this morphism.

        This is a helper method for :meth:`__call__` of the :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: T = morphism.codomain()

            sage: q = T(0, (1, 0)); q
            Point (1, 0) of polygon 0
            sage: (morphism * morphism.section())(q)
            Point (1, 0) of polygon 0

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute a preimage of a point yet"
        )

    def _test_section_point(self, **options):
        r"""
        Verify that :meth:`_section_point` actually produces a section.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism._test_section_point()

        """
        tester = self._tester(**options)

        section = self.section()
        identity = self * section

        for q in tester.some_elements(self.codomain().some_elements()):
            tester.assertEqual(identity(q), q)

    def _image_homology(self, g, codomain=None):
        r"""
        Return the image of the homology class ``g`` under this morphism.

        This is a helper method for :meth:`__call__`.

        Subclasses can override this method if the morphism is meaningful on
        the level of homology.

        However, it's usually easier to override :meth:`_image_homology_edge`,
        :meth:`_image_homology_gen`, or :meth:`_image_homology_matrix` to
        support mapping homology classes.

        INPUT:

        - ``codomain`` -- a simplicial homology or ``None`` (default:
          ``None``); if set, the homology where the result should live,
          otherwise, the result will live in the :meth:`homology` of the
          :meth:`codomain`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

        The image of a homology class::

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: a, b = H.gens()

            sage: morphism(a)
            B[(0, 1)]

        """
        if g.parent().surface() is not self.domain():
            raise ValueError("g must be a homology class over this morphism's domain")

        if codomain is None:
            codomain = self.codomain().homology()

        assert codomain.surface() is self.codomain(), "codomain must be a homology of the codomain() of this morphism"

        return g.parent().hom(self._image_homology_matrix(domain=g.parent(), codomain=codomain), codomain=codomain)(g)

    def _section_homology(self, h, codomain=None):
        r"""
        Return a preimage of the homology class ``h`` under this morphism.

        This is a helper method for :meth:`__call__` of the :meth:`section`.

        Subclasses can override this method if the morphism is meaningful on
        the level of homology.

        However, it's usually easier to override :meth:`_section_homology_edge`,
        :meth:`_section_homology_gen`, or :meth:`_section_homology_matrix` to
        support mapping homology classes.

        INPUT:

        - ``codomain`` -- a simplicial homology or ``None`` (default: ``None``);
          if set, the homology where the result should live; otherwise, the
          result will live in the :meth:`homology` of the :meth:`domain`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: T = morphism.codomain()

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: a
            B[(0, 1)]
            sage: (morphism * morphism.section())(a)
            B[(0, 1)]

        """
        return SurfaceMorphism._image_homology(self.section(), h, codomain=codomain)

    def _test_section_homology(self, **options):
        r"""
        Verify that :meth:`_section_homology` actually produces a
        section.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism._test_section_homology()

        """
        tester = self._tester(**options)

        section = self.section()
        identity = self * section

        homology = self.codomain().homology()

        for q in tester.some_elements(homology.some_elements()):
            tester.assertEqual(q.parent().hom(identity)(q), q)

    @cached_method
    def _image_homology_matrix(self, domain, codomain):
        r"""
        Return the matrix `M` describing how this morphism acts on homology,
        i.e., for a homology class given by a vector `c` with respect to a
        basis of homology of the domain, the image is `M c` with respect to the
        basis of homology of the codomain.

        This is a helper method for :meth:`__call__` and
        :meth:`_image_homology`.

        Subclasses can override this method if the morphism is meaningful (and
        linear) on the level of homology.

        However, it is often easier to override :meth:`_image_homology_edge` or
        :meth:`_image_homology_gen`.

        INPUT:

        - ``domain`` -- a simplicial homology over the :meth:`domain`

        - ``codomain`` -- a simplicial homology over the :meth:`codomain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._image_homology_matrix()
            [1 0]
            [0 1]

        """
        assert domain.surface() is self.domain()
        assert codomain.surface() is self.codomain()

        domain_gens = domain.gens()
        codomain_gens = codomain.gens()

        from sage.all import matrix, ZZ

        M = matrix(ZZ, len(codomain_gens), len(domain_gens), sparse=True)

        for x, domain_gen in enumerate(domain_gens):
            image = self._image_homology_gen(domain_gen, codomain=codomain)
            for y, codomain_gen in enumerate(codomain_gens):
                M[y, x] = image.coefficient(codomain_gen)

        M.set_immutable()
        return M

    def _section_homology_matrix(self, domain, codomain):
        r"""
        Return the matrix describing a section of this morphism on the level of
        homology, see :meth:`_image_homology_matrix`.

        INPUT:

        - ``domain`` -- a simplicial homology over the :meth:`codomain`

        - ``codomain`` -- a simplicial homology over the :meth:`domain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._section_homology_matrix()
            [1 0]
            [0 1]

        """
        M = self._image_homology_matrix(domain=codomain, codomain=domain)
        if M.rank() != M.nrows():
            raise NotImplementedError(
                "cannot compute section of homology matrix because map is not onto in homology"
            )
        return M.pseudoinverse().change_ring(M.base_ring())

    def _image_homology_gen(self, gen, codomain):
        r"""
        Return the image of a generator of homology ``gen``.

        This is a helper method for :meth:`__call__` and
        :meth:`_image_homology_matrix`.

        Subclasses can override this method if the morphism is meaningful (and
        linear) on the level of homology.

        However, it is often easier to override :meth:`_image_homology_edge` instead.

        INPUT:

        - ``codomain`` -- a simplicial homology over the :meth:`codomain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(S)
            sage: a, b = H.gens()

            sage: morphism._image_homology_gen(a)
            B[(0, 1)]
            sage: morphism._image_homology_gen(b)
            B[(0, 0)]

        """
        assert codomain.surface() is self.codomain()

        chain = gen._chain
        image = codomain.zero()
        for label, edge in chain.support():
            coefficient = chain[(label, edge)]
            assert coefficient
            image += coefficient * self._image_homology_edge(label, edge, codomain=codomain)

        return codomain(image)

    def _section_homology_gen(self, gen, codomain):
        r"""
        Return a preimage of the homology generator ``gen``.

        This is a helper method for :meth:`_image_homology_matrix` of
        :meth:`section`. But usually this is not invoked since we compute the
        section with linear algebra in :meth:`_section_homology_matrix`.

        INPUT:

        - ``codomain`` -- a simplicial homology over the :meth:`domain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: from flatsurf import SimplicialHomology
            sage: H = SimplicialHomology(morphism.codomain())
            sage: a, b = H.gens()

            sage: morphism._section_homology_gen(a)
            B[(0, 1)]

        """
        return SurfaceMorphism._image_homology_gen(self.section(), gen, codomain=codomain)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Return the image of the homology class generated by ``edge`` in the
        polygon ``label`` under this morphism.

        This is a helper method for :meth:`__call__` and
        :meth:`_image_homolyg_gen`.

        Subclasses can override this method if the morphism is meaningful (and
        linear) on the level of homology.

        INPUT:

        - ``codomain`` -- a simplicial homology over the :meth:`codomain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._image_homology_edge(0, 0)
            B[(0, 0)]
            sage: morphism._image_homology_edge(0, 1)
            B[(0, 1)]

        """
        assert codomain.surface() is self.codomain()

        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute the image of an edge yet"
        )

    def _section_homology_edge(self, label, edge, codomain):
        r"""
        Return a preimage of an edge in homology.

        This is a helper method for :meth:`_image_homology_matrix` of
        :meth:`section`. But usually this is not invoked since we compute the
        section with linear algebra in :meth:`_section_homology_matrix`.

        INPUT:

        - ``codomain`` -- a simplicial homology over the :meth:`domain`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._section_homology_edge(0, 0)

        """
        return SurfaceMorphism._image_homology_edge(self.section(), label, edge, codomain=codomain)

    def _composition(self, other):
        r"""
        Return the composition of this morphism and ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: g = f.codomain().apply_matrix(matrix([[1/2, 0], [0, 1]]), in_place=False)
            sage: g * f
            Composite morphism:
              From: Translation Surface in H_1(0) built from a square
              To:   Translation Surface in H_1(0) built from a square
              Defn:   Linear morphism:
                      From: Translation Surface in H_1(0) built from a square
                      To:   Translation Surface in H_1(0) built from a rectangle
                      Defn: [2 0]
                            [0 1]
                    then
                      Linear morphism:
                      From: Translation Surface in H_1(0) built from a rectangle
                      To:   Translation Surface in H_1(0) built from a square
                      Defn: [1/2   0]
                            [  0   1]

        """
        if other.codomain() is not self.domain():
            raise ValueError(
                f"morphisms cannot be composed because domain of\n{self}\nis not compatible with codomain of\n{other}"
            )

        if isinstance(other, IdentityMorphism):
            return self

        if isinstance(self, IdentityMorphism):
            return other

        return CompositionMorphism._create_morphism(self, other)

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: H = End(S)

            sage: H.identity() == H.identity()
            True

        """
        raise NotImplementedError(f"{type(self).__name__} does not implement __eq__ yet")

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: H = End(S)

            sage: hash(H.identity()) == hash(H.identity())
            True

        """
        raise NotImplementedError(f"{type(self).__name__} does not implement __hash__ yet")


class SectionMorphism(SurfaceMorphism):
    r"""
    The formal section of a morphism.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.regular_octagon()
        sage: morphism = S.subdivide_edges(2)
        sage: section = morphism.section()
        sage: section
        Section morphism:
          From: Translation Surface in H_2(2, 0^4) built from a regular octagon with 8 marked vertices
          To:   Translation Surface in H_2(2) built from a regular octagon
          Defn: Section of Edge-Subdivision morphism:
                  From: Translation Surface in H_2(2) built from a regular octagon
                  To:   Translation Surface in H_2(2, 0^4) built from a regular octagon with 8 marked vertices

    TESTS::

        sage: from flatsurf.geometry.morphism import SectionMorphism
        sage: isinstance(section, SectionMorphism)
        True

        sage: TestSuite(section).run()

    """

    def __init__(self, parent, morphism, category=None):
        self._morphism = morphism
        super().__init__(parent, category=category)

    def section(self):
        r"""
        Return a section of this section.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: morphism = S.subdivide_edges(2)
            sage: section = morphism.section()
            sage: section.section() == morphism
            True

        """
        return self._morphism

    @classmethod
    def _create_morphism(cls, morphism):
        r"""
        Return a formal section of ``morphism``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: morphism = S.subdivide_edges(2)

            sage: from flatsurf.geometry.morphism import SectionMorphism
            sage: SectionMorphism._create_morphism(morphism)
            Section morphism:
              From: Translation Surface in H_2(2, 0^4) built from a regular octagon with 8 marked vertices
              To:   Translation Surface in H_2(2) built from a regular octagon
              Defn: Section of Edge-Subdivision morphism:
                      From: Translation Surface in H_2(2) built from a regular octagon
                      To:   Translation Surface in H_2(2, 0^4) built from a regular octagon with 8 marked vertices

        """
        return super()._create_morphism(
            morphism.codomain(), morphism.domain(), morphism
        )

    def _image_point(self, p):
        r"""
        Implements :meth:`SurfaceMorphism._image_point`, see
        :meth:`SurfaceMorphism._section_point`.
        """
        return self._morphism._section_point(p)

    def _image_homology(self, g, codomain=None):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology`, see
        :meth:`SurfaceMorphism._section_homology`.
        """
        return self._morphism._section_homology(g, codomain=codomain)

    def _image_homology_gen(self, gen, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_gen`, see
        :meth:`SurfaceMorphism._section_homology_gen`.
        """
        assert codomain.surface() is self.codomain()

        return self._morphism._section_homology_gen(gen, codomain=codomain)

    def _image_homology_matrix(self, domain, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_matrix`, see
        :meth:`SurfaceMorphism._section_homology_matrix`.
        """
        assert domain.surface() is self.domain()
        assert codomain.surface() is self.codomain()

        return self._morphism._section_homology_matrix(domain=domain, codomain=codomain)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: g = f.section()

            sage: g._image_homology_edge((0, 0), 0, S.homology())
            B[(0, 0)]

        """
        assert codomain.surface() is self.codomain()

        domain = self.domain().homology()

        for (l, e) in self.codomain().edges():
            if self._morphism._image_homology_edge(l, e, codomain=domain) == domain((label, edge)):
                return codomain((l, e))

        raise NotImplementedError("not a single edge maps to this edge, cannot implement preimage of this edge yet")

    def __eq__(self, other):
        r"""
        Return whether this section is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: f.section() == f.section()
            True

        """
        if not isinstance(other, SectionMorphism):
            return False

        return self._morphism == other._morphism

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: hash(f.section()) == hash(f.section())
            True

        """
        return -hash(self._morphism)

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: f.section()

        """
        return "Section"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: f.section()

        """
        return f"Section of {self._morphism}"


class CompositionMorphism(SurfaceMorphism):
    r"""
    The formal composition of two morphisms between surfaces.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.octagon()
        sage: f = S.triangulate()
        sage: T = f.codomain()
        sage: g = T.delaunay_triangulate()

        sage: composition = g * f
        sage: composition

    TESTS::

        sage: from flatsurf.geometry.morphism import CompositionMorphism
        sage: isinstance(composition, CompositionMorphism)
        sage: TestSuite(composition).run()

    """
    def __init__(self, parent, lhs, rhs, category=None):
        super().__init__(parent, category=category)

        # The morphisms as they are executed (in chronological order.)
        self._morphisms = []
        for morphism in [rhs, lhs]:
            if isinstance(morphism, CompositionMorphism):
                self._morphisms.extend(morphism._morphisms)
            else:
                self._morphisms.append(morphism)

    @classmethod
    def _create_morphism(cls, lhs, rhs):
        r"""
        Create a composition of ``lhs`` to be executed after ``rhs``.

        There should be no need to use this method directly, use ``*`` instead.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()

            sage: from flatsurf.geometry.morphism import CompositionMorphism
            sage: CompositionMorphism._create_morphism(g, f) == g * f
            True

        """
        return super()._create_morphism(rhs.domain(), lhs.codomain(), lhs, rhs)

    def __call__(self, x):
        r"""
        Evaluate this morphism on the point ``x``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()
            sage: h = g * f

            sage: p = S(0, 0)

            sage: q = h(p)

            sage: q == g(f(p))
            True

        """
        for morphism in self._morphisms:
            x = morphism(x)
        return x

    def _image_homology_matrix(self, domain, codomain):
        r"""
        Helper method to compute the image of homology under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()

            sage: H = S.homology()
            sage: h.hom(g * f)

            sage: h.matrix()

            sage: (g * f)._image_homology_matrix(H, g.codomain().homology())

        """
        assert domain.surface() is self.domain()
        assert codomain.surface() is self.codomain()

        matrix = IdentityMorphism._create_morphism(self.codomain())._image_homology_matrix(domain=codomain, codomain=codomain)

        for morphism in self._morphisms[:0:-1]:
            matrix *= morphism._image_homology_matrix(domain=morphism.domain().homology(), codomain=codomain)
            codomain = morphism.domain().homology()

        matrix *= self._morphisms[0]._image_homology_matrix(domain=domain, codomain=codomain)

        return matrix

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()

            sage: g * f == g * f
            True

        """
        if not isinstance(other, CompositionMorphism):
            return False

        return self._parent == other._parent and self._morphisms == other._morphisms

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible wit :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()

            sage: hash(g * f) == hash(g * f)

        """
        return hash(tuple(self._morphisms))

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()
            sage: g * f

        """
        return "Composite"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()
            sage: g * f

        """
        return "\nthen ".join(str(morphism) for morphism in self._morphisms)

    @cached_method
    def section(self):
        r"""
        Return a section of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon()
            sage: f = S.triangulate()
            sage: T = f.codomain()
            sage: g = T.delaunay_triangulate()
            sage: h = g * f
            sage: h

            sage: h.section()

        """
        from sage.all import prod

        return prod(morphism.section() for morphism in self._morphisms)


class IdentityMorphism(SurfaceMorphism, IdentityMorphism_sage):
    r"""
    The identity morphism from a surface to itself.

    EXAMPLES::

       sage: from flatsurf import translation_surfaces
       sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
       sage: identity = S.erase_marked_points()
       sage: identity
       Identity endomorphism of Translation Surface in H_2(2) built from 3 squares

    TESTS::

        sage: from flatsurf.geometry.morphism import IdentityMorphism
        sage: isinstance(identity, IdentityMorphism)
        True

        sage: TestSuite(identity).run()

    The morphism can be evaluated at points (inherited from the base class)::

        sage: identity(S(0, 0))
        Vertex 0 of polygon 0

    This is the neutral element for composition::

        sage: triangulation = S.triangulate()
        sage: triangulation * identity == triangulation
        True

        sage: End(triangulation.codomain()).identity() * triangulation == triangulation
        True

    """

    def __init__(self, parent, category=None):
        if parent.domain() is not parent.codomain():
            raise ValueError("domain and codomain of identity must be identical")

        super().__init__(parent, category=category)

    @classmethod
    def _create_morphism(cls, domain):
        r"""
        Return the identity morphism on ``domain``.

        EXAMPLES::

            sage: from flatsurf.geometry.morphism import translation_surfaces, IdentityMorphism
            sage: S = translation_surfaces.square_torus()
            sage: f = IdentiyMorphism._create_morphism(S)
            sage: f == End(S).identity()

        """
        return super()._create_morphism(domain, domain)

    def section(self):
        r"""
        Return an inverse of this morphism, i.e., this morphism itself.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: identity = S.erase_marked_points()
            sage: identity.section() == identity
            True

        """
        return self

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Helper method to compute the image of homology under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = End(S).identity()
            sage: f._image_homology_edge(0, 0, codomain=S.homology())
            B[(0, 0)]

        """
        assert codomain.surface() is self.codomain()

        return codomain((label, edge))

    def _repr_type(self):
        r"""
        Return a printable representation of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: identity = S.erase_marked_points()
            sage: identity
            Identity endomorphism of Translation Surface in H_2(2) built from 3 squares

        """
        return "Identity"

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: identity = S.erase_marked_points()
            sage: identitz = S.erase_marked_points()
            sage: identity == identitz
            True

        A morphism that is the identity but still distinguishable from the
        plain identity::

            sage: identitz = S.apply_matrix(matrix([[1, 0], [0, 1]]), in_place=False)
            sage: identity == identitz
            False

        """
        if not isinstance(other, IdentityMorphism):
            return False
        return self.domain() == other.domain()


class SurfaceMorphism_factorization(SurfaceMorphism):
    r"""
    A morphism between surfaces that is backed by another morphism, usually a
    :class:`CompositionMorphism`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.pyflatsurf()
        sage: f

        sage: f._factorization()

    TESTS::

        sage: from flatsurf.geometry.morphism import SurfaceMorphism_factorization
        sage: isinstance(f, SurfaceMorphism_factorization)
        True

        sage: TestSuite(f).run()

    """
    def _image_homology(self, g, codomain=None):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: a, b = S.homology().gens()
            sage: f._image_homology(a)

        """
        return self._factorization()._image_homology(g, codomain=codomain)

    def _section_homology(self, h, codomain=None):
        r"""
        Implements :meth:`SurfaceMorphism._section_homology` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: a, b = f.codomain().homology().gens()
            sage: f._section_homology(a)

        """
        return self._factorization()._section_homology(h, codomain=codomain)

    def _image_homology_matrix(self, domain, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_matrix` by invoking
        the underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._image_homology_matrix(S, f.codomain())

        """
        assert domain.surface() is self.domain()
        assert codomain.surface() is self.codomain()

        return self._factorization()._image_homology_matrix(domain=domain, codomain=codomain)

    def _section_homology_matrix(self, domain, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_matrix` by invoking
        the underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._section_homology_matrix(f.codomain(), S)

        """
        assert domain.surface() is self.codomain()
        assert codomain.surface() is self.domain()

        return self._factorization()._section_homology_matrix(domain=domain, codomain=codomain)

    def _image_homology_gen(self, gen, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_gen` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: a, b = S.homology().gens()
            sage: f._image_homology_gen(a)

        """
        assert codomain.surface() is self.codomain()

        return self._factorization()._image_homology_gen(gen, codomain=codomain)

    def _section_homology_gen(self, gen, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._section_homology_gen` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: a, b = f.codomain().homology().gens()
            sage: f._section_homology_gen(a)

        """
        assert codomain.surface() is self.domain()

        return self._factorization()._section_homology_gen(gen, codomain=codomain)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implement :meth:`SurfaceMorphism._image_homology_edge` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._image_homology_edge(0, 0, codomain=f.codomain().homology())

        """
        assert codomain.surface() is self.codomain()

        return self._factorization()._image_homology_edge(label, edge, codomain)

    def _section_homology_edge(self, label, edge, codomain):
        r"""
        Implement :meth:`SurfaceMorphism._section_homology_edge` by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._section_homology_edge((1, 2, 3), 0, codomain=f.domain().homology())

        """
        assert codomain.surface() is self.domain()

        return self._factorization()._image_homology_edge(label, edge, codomain)

    def _image_point(self, p):
        r"""
        Return the image of a point under this morphism by invoking the
        underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()
            sage: f(S(0, 0))

        """
        return self._factorization()._image_point(p)

    def _section_point(self, q):
        r"""
        Return a preimage of a point under this morphism by determining a
        preimage of the underlying morphism of this factorization.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()
            sage: f._section_point(f.codomain()((1, 2, 3), 0))

        """
        return self._factorization()._section_point(q)

    def _factorization(self):
        r"""
        Return the morphism underlying this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._factorization()

        """
        raise NotImplementedError("this factorization morphism is not functional because it does not implement factorization()")

    def _test_factorization(self, **options):
        r"""
        Verify that :meth:`_factorization` has been implemented correctly.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.pyflatsurf()

            sage: f._test_factorization()

        """
        tester = self._tester(**options)

        tester.assertTrue(self.domain() is self._factorization().domain())
        tester.assertTrue(self.codomain() is self._factorization().codomain())


class NamedFactorizationMorphism(SurfaceMorphism_factorization):
    r"""
    A morphism that is implement by another morphism but has a more appealing
    name when printing.

    This is commonly used when implementing new morphisms that rely on a
    sequence of existing morphisms while providing a nice interface to the user
    of sage-flatsurf.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.triangulate()
        sage: g = f.codomain().apply_matrix(matrix([[1, 2], [0, 1]]), in_place=False)

        sage: from flatsurf.geometry.morphism import NamedFactorizationMorphism
        sage: h = NamedFactorizationMorphism._create_morphism(S, g.codomain(), "Foo", g * f)
        sage: h

    TESTS::

        sage: from flatsurf.geometry.
        sage: isinstance(h, NamedFactorizationMorphism)
        True

        sage: TestSuite(h).run()

    """
    def __init__(self, parent, name, factorization):
        super().__init__(parent)

        self._name = name
        self.__factorization = factorization

        if self.__factorization.domain() is not self.domain():
            raise ValueError
        if self.__factorization.codomain() is not self.codomain():
            raise ValueError

    def _repr_type(self):
        r"""
        Helper method for :meth;`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: g = f.codomain().apply_matrix(matrix([[1, 2], [0, 1]]), in_place=False)

            sage: from flatsurf.geometry.morphism import NamedFactorizationMorphism
            sage: h = NamedFactorizationMorphism._create_morphism(S, g.codomain(), "Foo", g * f)
            sage: h

        """
        return self._name

    def _factorization(self):
        r"""
        Return the morphism underlying this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: g = f.codomain().apply_matrix(matrix([[1, 2], [0, 1]]), in_place=False)

            sage: from flatsurf.geometry.morphism import NamedFactorizationMorphism
            sage: h = NamedFactorizationMorphism._create_morphism(S, g.codomain(), "Foo", g * f)

            sage: h._factorization()

        """
        return self.__factorization


class TriangulationMorphism_base(SurfaceMorphism):
    r"""
    Abstract base class for morphisms from a surface to its triangulation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces

        sage: G = SymmetricGroup(4)
        sage: S = translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
        sage: f = S.triangulate()
        sage: f
        Triangulation morphism:
          From: Origami defined by r=(1,2,3,4) and u=(1,4,2,3)
          To:   Triangulation of Origami defined by r=(1,2,3,4) and u=(1,4,2,3)

    TESTS::

        sage: from flatsurf.geometry.morphism import TriangulationMorphism_base
        sage: isinstance(f, TriangulationMorphism_base)
        True

        sage: TestSuite(f).run()

    """

    def _image_edge(self, label, edge):
        r"""
        Return the image of the ``edge`` of the polygon with ``label`` in the
        codomain of this triangulation.

        Subclasses must implement this method to make this morphism functional.

        Since this is a triangulation, all edges of the domain are preserved in
        the codomain.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: S = translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: f = S.triangulate()
            sage: f._image_edge(0, 0)

        """
        raise NotImplementedError("this triangulation morphism is not functional because it does not implement _image_edge() yet")

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :class:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: G = SymmetricGroup(4)
            sage: S = translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: triangulation = S.triangulate()
            sage: triangulation._image_homology_edge(1, 0)
            [(1, (1, 0), 0)]

        """
        assert codomain.surface() is self.codomain()

        return codomain(self._image_edge(label, edge))

    def _image_point(self, p):
        r"""
        Implements :class:`SurfaceMorphism._image_point`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: G = SymmetricGroup(4)
            sage: S = translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: triangulation = S.triangulate()
            sage: triangulation._image_point(S(1, 0))
            Vertex 0 of polygon (1, 0)
            sage: triangulation._image_point(S(1, (1/2, 1/2)))
            Point (1/2, 1/2) of polygon (1, 0)

        """
        preimage_label, preimage_coordinates = p.representative()
        preimage_polygon = self.domain().polygon(preimage_label)

        for preimage_edge in range(len(preimage_polygon.edges())):
            relative_coordinates = preimage_coordinates - preimage_polygon.vertex(
                preimage_edge
            )
            image_label, image_edge = self._image_edge(preimage_label, preimage_edge)
            image_polygon = self.codomain().polygon(image_label)
            image_coordinates = image_polygon.vertex(image_edge) + relative_coordinates
            if image_polygon.contains_point((image_coordinates)):
                return self.codomain()(image_label, image_coordinates)

        assert (
            False
        ), "point must be in one of the polygons that split up the original polygon"

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: S = translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: f = S.triangulate()
            sage: f

        """
        return "Triangulation"


class TriangulationMorphism_LazyTriangulatedSurface(TriangulationMorphism_base):
    r"""
    A triangulation that maps a surface into a class :class:`LazyTriangulatedSurface`.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.triangulate()
        sage: f

    TESTS::

        sage: from flatsurf.geometry.morphism import TriangulationMorphism_lazyTriangulatedSurface
        sage: isinstance(f, LazyTriangulatedSurface)
        True
        sage: TestSuite(f).run()

    """
    def _image_edge(self, label, edge):
        r"""
        Implements :meth:`TriangulationMorphism_base._image_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.triangulate()
            sage: f._image_edge(0, 0)
            ((0, 0), 0)

        """
        return self.codomain()._image(label)[1][edge]

    def __eq__(self, other):
        r"""
        Return whether this triangulation is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.triangulate() == S.triangulate()
            True

        """
        if not isinstance(other, TriangulationMorphism_LazyTriangulatedSurface):
            return False

        return self.parent() == other.parent()

    def __hash__(self):
        r"""
        Return a hash value for this triangulation that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: hash(S.triangulate()) == hash(S.triangulate())
            True

        """
        return hash(self.parent())


class DelaunayTriangulationMorphism_delaunay_decomposition(TriangulationMorphism_base):
    r"""
    A triangulation that maps a a Delaunay decomposition into the Delaunay
    triangulation from which it was created.

    This morphism is used internally in :meth:`delaunay_decompose`, the actual
    delaunay decomposition, is the inverse of this morphism.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.delaunay_decompose()._factorization()._morphisms[1].section()

    TESTS::

        sage: from flatsurf.geometry.morphism import DelaunayTriangulationMorphism_delaunay_decomposition
        sage: isinstance(f, DelaunayTriangulationMorphism_delaunay_decomposition)
        True

        sage: TestSuite(f).run()

    """
    def _image_edge(self, label, edge):
        r"""
        Implements :meth:`TriangulationMorphism_base._image_edge.`

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1].section()
            sage: f._image_edge((0, 0), 0)
            ((0, 0), 0)

        """
        _, edges = self.domain()._cell(label)
        return edges[edge]

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1].section()
            sage: g = S.delaunay_decompose()._factorization()._morphisms[1].section()
            sage: f == g
            True

        """
        if not isinstance(other, DelaunayTriangulationMorphism_delaunay_decomposition):
            return False

        return self.parent() == other.parent()

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1].section()
            sage: g = S.delaunay_decompose()._factorization()._morphisms[1].section()
            sage: hash(f) == hash(g)
            True

        """
        return hash(self.parent())


class DelaunayTriangulationMorphism(SurfaceMorphism):
    r"""
    A morphism from a triangulated surface to its Delaunay triangulation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.delaunay_triangulate()._morphisms[1]
        sage: f

    TESTS::

        sage: isinstance(f, DelaunayTriangulationMorphism)
        True

        sage: TestSuite(f).run()

    """
    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_triangulate()._morphisms[1]

            sage: f._image_homology_edge((0, 0), 0, codomain=f.codomain().homology())

        """
        gens, chain_matrix = self._image_homology_edge_chain_matrix()

        image = chain_matrix.column(gens.index((label, edge)))

        return sum(c * codomain(gens[i]) for (i, c) in zip(image.support(), image.coefficients()))

    @cached_method
    def _image_homology_edge_chain_matrix(self):
        r"""
        Return the matrix that describes the action of this morphism in
        homology on the generators of the chain module.

        Helper method for :meth:`_image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_triangulate()._morphisms[1]
            sage: f._image_homology_edge_chain_matrix()

        """
        for label in self.codomain().labels():
            self.codomain()._certify(label)

        gens = tuple(self.domain().edges())

        from sage.all import identity_matrix, ZZ
        chain_matrix = identity_matrix(ZZ, len(gens), sparse=True)

        for flip in self.codomain()._flips[::-1]:
            chain_matrix *= self._image_homology_edge_chain_matrix_flip(flip, gens)

        return gens, chain_matrix

    def _image_homology_edge_chain_matrix_flip(self, flip, gens):
        r"""
        Return the matrix that describes the action of a single ``flip`` on
        homology on the generators ``gens`` of the chain module.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_triangulate()._morphisms[1]
            sage: f._image_homology_edge_chain_matrix_flip((0, 0), tuple(S.edges()))

        TESTS:

        This method should be migrated into the morphism describing a triangle
        flip but currently there is no such morphism::

            sage: isinstance(S.triangle_flip(0, 0), SurfaceMorphism)
            False

        """
        from sage.all import zero_matrix, ZZ

        (label, edge), (opposite_label, opposite_edge) = flip


        matrix = zero_matrix(ZZ, len(gens), sparse=True)

        for i, (lbl, e) in enumerate(gens):
            if lbl != label and lbl != opposite_label:
                matrix[i, i] = 1
                continue

        def image(preimage, *images):
            def denormalize(lbl, e):
                if lbl == label:
                    return lbl, (e + edge) % 3
                if lbl == opposite_label:
                    return lbl, (e + opposite_edge) % 3
                assert False

            preimage = denormalize(*preimage)
            for image in images:
                image = denormalize(*image)
                matrix[gens.index(image), gens.index(preimage)] += 1

        image((label, 0), (opposite_label, 1), (label, 2))
        image((label, 1), (opposite_label, 2))
        image((label, 2), (label, 1))

        image((opposite_label, 0), (label, 1), (opposite_label, 2))
        image((opposite_label, 1), (label, 2))
        image((opposite_label, 2), (opposite_label, 1))

        return matrix

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_triangulate()._morphisms[1]
            sage: g = S.delaunay_triangulate()._morphisms[1]
            sage: f == g
            True

        """
        if not isinstance(other, DelaunayTriangulationMorphism):
            return False

        return self.parent() == other.parent()

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_triangulate()._morphisms[1]
            sage: g = S.delaunay_triangulate()._morphisms[1]
            sage: hash(f) == hash(g)
            True

        """
        return hash(self.parent())


class DelaunayDecompositionMorphism(SectionMorphism):
    r"""
    A morphism describing the Delaunay decomposition of a Delaunay triangulation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.delaunay_decompose()._factorization()._morphisms[1]
        sage: f

    This morphism is implemented as the formal section of a triangulation
    morphism::

        sage: f.section()

    TESTS::

        sage: from flatsurf.geometry.morphism import DelaunayDecompositionMorphism
        sage: isinstance(DelaunayDecompositionMorphism)
        True

        sage: TestSuite(f).run()

    """
    def __init__(self, parent, category=None):
        super().__init__(parent, DelaunayTriangulationMorphism_delaunay_decomposition._create_morphism(parent.codomain(), parent.domain()), category=category)

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1]
            sage: f

        """
        return "Delaunay decomposition"

    @classmethod
    def _create_morphism(cls, domain, codomain):
        r"""
        Return a morphism from ``domain`` to ``codomain``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().delaunay_triangulate().codomain()
            sage: T = S.delaunay_decompose().codomain()

            sage: from flatsurf.geometry.morphism import DelaunayDecompositionMorphism
            sage: DelaunayDecompositionMorphism._create_morphism(S, T)

        """
        return super()._create_morphism(DelaunayTriangulationMorphism_delaunay_decomposition._create_morphism(codomain, domain))

    def __eq__(self, other):
        r"""
        Return whether this morphism is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1]
            sage: g = S.delaunay_decompose()._factorization()._morphisms[1]

            sage: f == g
            True

        """
        if not isinstance(other, DelaunayDecompositionMorphism):
            return False

        return self.parent() == other.parent()

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.delaunay_decompose()._factorization()._morphisms[1]
            sage: g = S.delaunay_decompose()._factorization()._morphisms[1]

            sage: hash(f) == hash(g)
            True

        """
        return hash(self.parent())


class DelaunayDecompositionIsomorphism(SurfaceMorphism_factorization):
    r"""
    An isomorphism between Delaunay decompositions.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: T = S.relabel({0: 1})
        sage: isomorphism = S.delaunay_decompose(codomain=T)
        sage: isomorphism

    TESTS::

        sage: from flatsurf.geometry.morphism import DelaunayDecompositionIsomorphism
        sage: isinstance(isomorphism, DelaunayDecompositonIsomorphism)
        True

        sage: TestSuite(isomorphism).run()

    """
    @cached_method
    def _factorization(self):
        r"""
        Return the sequence of morphism that implement this isomorphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: T = S.relabel({0: 1})
            sage: isomorphism = S.delaunay_decompose(codomain=T)
            sage: isomorphism._factorization()

        """
        from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

        domain_to_pyflatsurf = Surface_pyflatsurf._from_flatsurf(self.domain())
        codomain_to_pyflatsurf = Surface_pyflatsurf._from_flatsurf(self.codomain())

        from pyflatsurf import flatsurf
        deformation = domain_to_pyflatsurf.codomain().flat_triangulation().isomorphism(codomain_to_pyflatsurf.codomain().flat_triangulation(), flatsurf.ISOMORPHISM.DELAUNAY_CELLS).value()

        from flatsurf.geometry.pyflatsurf.morphism import Morphism_Deformation
        return codomain_to_pyflatsurf.section() * Morphism_Deformation._create_morphism(domain_to_pyflatsurf.codomain(), codomain_to_pyflatsurf.codomain(), deformation) * domain_to_pyflatsurf


class GL2RMorphism(SurfaceMorphism):
    r"""
    A morphism that describes the application of a matrix in `GL_2(\mathbb{R})`
    on each polygon of a surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
        sage: f

    TESTS::

        sage: from flatsurf.geometry.morphism import GL2RMorphism
        sage: isinstance(f, GL2RMorphism)
        True

        sage: TestSuite(f).run()

    """
    def __init__(self, parent, m, category=None):
        super().__init__(parent, category=category)

        from sage.all import matrix
        self._matrix = matrix(self.domain().base_ring(), m, immutable=True)

    def section(self):
        r"""
        Return an inverse of this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f.section()

        """
        return self._create_morphism(self.codomain(), self.domain(), ~self._matrix)

    def _image_point(self, point):
        r"""
        Return the image of the ``point`` under this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f(S(0, 0))
            Vertex 0 of polygon 0

        """
        label, coordinates = point.representative()
        return self.codomain()(label, self._matrix * coordinates)

    def _image_homology_edge(self, label, edge, codomain):
        r"""
        Implements :meth:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f._image_homology_edge(0, 0, codomain=f.codomain().homology())

        """
        assert codomain.surface() is self.codomain()

        return codomain((label, edge))

    def _repr_type(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f

        """
        return "Linear"

    def _repr_defn(self):
        r"""
        Helper method for :meth:`__repr__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f

        """
        return repr(self._matrix)

    def __eq__(self, other):
        r"""
        Return whether this morphis is indistinguishable from ``other``

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: g = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: f == g
            True

        """
        if not isinstance(other, GL2RMorphism):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    def __hash__(self):
        r"""
        Return a hash value for this morphism that is compatible with :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: f = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: g = S.apply_matrix(matrix([[1, 2], [0, 1]], in_place=False)
            sage: hash(f) == hash(g)
            True

        """
        return hash(self._matrix)
