r"""
Morphisms between Surfaces

.. NOTE::

    One should think of these as maps between surfaces. However, they might
    often not be maps on the points of the surface but just on homology, or in
    some cases, they might not be meaningful as maps anywhere.

    Technically, these are at worst morphisms in the category of objects where
    being a morphism does not really have any mathematical meaning.

.. NOTE::

    Our morphism infrastructure contains quite a few workarounds to make the
    SageMath machinery work. The fundemantel problem that we are facing is that
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
from sage.categories.morphism import Morphism


class SurfaceMorphismSpace(Homset):
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

    def __init__(self, domain, codomain, category=None):
        from sage.categories.all import Objects

        self.__category = category or Objects()
        super().__init__(domain, codomain, category=self.__category)

    def _an_element_(self):
        if self.is_endomorphism_set():
            return self.identity()
        raise NotImplementedError(
            f"cannot create a morphism from {self.domain()} to {self.codomain()} yet"
        )

    def identity(self):
        if self.is_endomorphism_set():
            return IdentityMorphism._create_morphism(self.domain())
        return super().identity()

    def __repr__(self):
        if self.domain() is self.codomain():
            return f"Surface Endomorphisms of {self.domain()}"
        return f"Surface Morphisms from {self.domain()} to {self.codomain()}"

    # We fail tests for associativity because we do not actually decide whether
    # two morphisms are "equal" but just whether they are indistinguishable.
    # Consequently, identity * identity != identity.
    # We disable tests that fails because of this.
    def _test_associativity(self, **options):
        pass

    def _test_one(self, **options):
        pass

    def _test_prod(self, **options):
        pass

    def __reduce__(self):
        return SurfaceMorphismSpace, (self.domain(), self.codomain(), self.__category)

    def __eq__(self, other):
        if not isinstance(other, SurfaceMorphismSpace):
            return False
        return (
            self.domain() == other.domain()
            and self.codomain() == other.codomain()
            and self.category() == other.category()
        )

    def __hash__(self):
        return hash((self.domain(), self.codomain()))


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

        A lot of methods that used to return a surface now return a morphism.
        To make transition of existing code easier, we look up attributes that
        cannot be found on the morphism up on the codomain and issue a
        deprecation warning.

        EXAMPLES::

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
        from flatsurf.geometry.surface_objects import SurfacePoint_base

        if isinstance(x, SurfacePoint_base):
            if x.parent() is not self.domain():
                raise ValueError("point must be in the domain of this morphism")
            image = self._image_point(x)
            assert image.parent() is self.codomain()
            return image

        from flatsurf.geometry.homology import SimplicialHomologyClass

        if isinstance(x, SimplicialHomologyClass):
            if x.parent().surface() is not self.domain():
                raise ValueError(
                    "homology class must be defined over the domain of this morphism"
                )
            image = self._image_homology(x)
            assert image.parent().surface() is self.codomain()
            return image

        from flatsurf.geometry.saddle_connection import SaddleConnection_base

        if isinstance(x, SaddleConnection_base):
            if x.surface() is not self.domain():
                raise ValueError(
                    "saddle connection must be in the domain of this morphism"
                )
            image = self._image_saddle_connection(x)
            assert image.surface() is self.codomain()
            return image

        from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector

        if isinstance(x, SimilaritySurfaceTangentVector):
            if x.surface() is not self.domain():
                raise ValueError(
                    "tangent vector must be on the domain of this morphism"
                )
            image = self._image_tangent_vector(x)
            assert image.surface() is self.codomain()
            return image

        raise NotImplementedError(f"cannot map {type(x).__name__} yet through {self}")

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

        Not all morphisms are meaningful on the level of points::

            # TODO: Add an example of such a morphism.

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

    def _image_saddle_connection(self, c):
        r"""
        Return the image of saddle connection ``c`` under this morphism.

        This is a helper method for :meth:`__call__`.

        Subclasses should implement this method if the morphism is meaningful
        on the level of saddle connections.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

        The image of a saddle connection::

            sage: from flatsurf.geometry.saddle_connection import SaddleConnection
            sage: c = SaddleConnection.from_vertex(S, 0, 0, (1, 1))
            sage: c
            Saddle connection (1, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0

            sage: morphism(c)
            Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0

        Not all morphisms are meaningful on the level of saddle connections::

            sage: morphism = S.subdivide()
            sage: morphism(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: a InsertMarkedPointsInFaceMorphism_with_category cannot compute the image of a saddle connection yet

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute the image of a saddle connection yet"
        )

    def _section_saddle_connection(self, t):
        r"""
        Return a preimage of a saddle connection ``t`` under this morphism.

        This is a helper method for :meth:`__call__` of the :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: T = morphism.codomain()

            sage: from flatsurf.geometry.saddle_connection import SaddleConnection
            sage: t = SaddleConnection.from_vertex(T, 0, 0, (2, 1)); t
            Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0
            sage: (morphism * morphism.section())(t)
            Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute a preimage of a saddle connection yet"
        )

    def _test_section_saddle_connection(self, **options):
        r"""
        Verify that :meth:`_section_saddle_connection` actually produces a
        section.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism._test_section_saddle_connection()

        """
        tester = self._tester(**options)

        section = self.section()
        identity = self * section

        from itertools import islice

        for q in tester.some_elements(islice(self.codomain().saddle_connections(), 16)):
            tester.assertEqual(identity(q), q)

    def _image_tangent_vector(self, t):
        r"""
        Return the image of the tangent vector ``v`` under this morphism.

        This is a helper method for :meth:`__call__`.

        Subclasses should implement this method if the morphism is meaningful
        on the level of tangent vectors.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

        The image of a tangent vector::

            sage: t = S.tangent_vector(0, (0, 0), (1, 1))
            sage: morphism(t)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, 1)

        Not all morphisms are meaningful on the level of tangent vectors::

            # TODO: Add an example of such a morphism

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute the image of a tangent vector yet"
        )

    def _section_tangent_vector(self, q):
        r"""
        Return a preimage of the tangent vector ``q`` under this morphism.

        This is a helper method for :meth:`__call__` of the :meth:`section`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: T = morphism.codomain()

            sage: q = T.tangent_vector(0, (0, 0), (2, 1)); q
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, 1)
            sage: (morphism * morphism.section())(q)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, 1)

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute a preimage of a tangent vector yet"
        )

    def _test_section_tangent_vector(self, **options):
        r"""
        Verify that :meth:`_section_tangent_vector` actually produces a
        section.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism._test_section_tangent_vector()

        """
        tester = self._tester(**options)

        section = self.section()
        identity = self * section

        for q in tester.some_elements(self.codomain().tangent_bundle().some_elements()):
            tester.assertEqual(identity(q), q)

    def _image_homology(self, g):
        r"""
        Return the image of the homology class ``g`` under this morphism.

        This is a helper method for :meth:`__call__`.

        Subclasses can override this method if the morphism is meaningful on
        the level of homology.

        However, it's usually easier to override :meth:`_image_homology_edge`,
        :meth:`_image_homology_gen`, or :meth:`_image_homology_matrix` to
        support mapping homology classes.

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

        Not all morphisms are meaningful on the level of homology::

            # TODO: Add an example of such a morphism

        """
        from flatsurf.geometry.homology import SimplicialHomology

        codomain_homology = SimplicialHomology(self.codomain())

        from sage.all import vector

        image = self._image_homology_matrix() * vector(g._homology())

        homology, to_chain, to_homology = codomain_homology._homology()

        image = sum(
            coefficient * gen for (coefficient, gen) in zip(image, homology.gens())
        )

        return codomain_homology(to_chain(image))

    def _section_homology(self, h):
        r"""
        Return a preimage of the homology class ``h`` under this morphism.

        This is a helper method for :meth:`__call__` of the :meth:`section`.

        Subclasses can override this method if the morphism is meaningful on
        the level of homology.

        However, it's usually easier to override :meth:`_section_homology_edge`,
        :meth:`_section_homology_gen`, or :meth:`_section_homology_matrix` to
        support mapping homology classes.

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
        # Invoke the generic machinery by computing a homology matrix instead
        # TODO; No! Instead solve_right() with the image matrix.
        return SurfaceMorphism._image_homology(self.section(), h)

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
            tester.assertEqual(identity(q), q)

    @cached_method
    def _image_homology_matrix(self):
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

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._image_homology_matrix()
            [1 0]
            [0 1]

        """
        from flatsurf.geometry.homology import SimplicialHomology

        domain_homology = SimplicialHomology(self.domain())
        codomain_homology = SimplicialHomology(self.codomain())

        domain_gens = domain_homology.gens()
        codomain_gens = codomain_homology.gens()

        from sage.all import matrix, ZZ

        M = matrix(ZZ, len(codomain_gens), len(domain_gens), sparse=True)

        for x, domain_gen in enumerate(domain_gens):
            image = self._image_homology_gen(domain_gen)
            for y, codomain_gen in enumerate(codomain_gens):
                M[y, x] = image.coefficient(codomain_gen)

        M.set_immutable()
        return M

    def _section_homology_matrix(self):
        r"""
        Return the matrix describing a section of this morphism on the level of
        homology, see :meth:`_image_homology_matrix`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._section_homology_matrix()
            [1 0]
            [0 1]

        """
        M = self._image_homology_matrix()
        if M.rank() != M.nrows():
            raise NotImplementedError(
                "cannot compute section of homology matrix because map is not onto in homology"
            )
        return M.pseudoinverse()

    def _image_homology_gen(self, gen):
        r"""
        Return the image of a generator of homology ``gen``.

        This is a helper method for :meth:`__call__` and
        :meth:`_image_homology_matrix`.

        Subclasses can override this method if the morphism is meaningful (and
        linear) on the level of homology.

        However, it is often easier to override :meth:`_image_homology_edge` instead.

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
        from flatsurf.geometry.homology import SimplicialHomology

        codomain_homology = SimplicialHomology(self.codomain())

        chain = gen._chain
        image = codomain_homology.zero()
        for label, edge in chain.support():
            coefficient = chain[(label, edge)]
            assert coefficient
            for (multiplicity, label, edge) in self._image_homology_edge(label, edge):
                image += multiplicity * coefficient * image.parent()((label, edge))

        return codomain_homology(image)

    def _section_homology_gen(self, gen):
        r"""
        Return a preimage of the homology generator ``gen``.

        This is a helper method for :meth:`_image_homology_matrix` of
        :meth:`section`. But usually this is not invoked since we compute the
        section with linear algebra in :meth:`_section_homology_matrix`.
        """
        return SurfaceMorphism._image_homology_gen(self.section(), gen)

    def _image_homology_edge(self, label, edge):
        r"""
        Return the image of the homology class generated by ``edge`` in the
        polygon ``label`` under this morphism.

        Returns the image as a sequence of triples (multiplicity, label, edge).

        This is a helper method for :meth:`__call__` and
        :meth:`_image_homolyg_gen`.

        Subclasses can override this method if the morphism is meaningful (and
        linear) on the level of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)

            sage: morphism._image_homology_edge(0, 0)
            [(1, 0, 0)]
            sage: morphism._image_homology_edge(0, 1)
            [(1, 0, 1)]

        """
        raise NotImplementedError(
            f"a {type(self).__name__} cannot compute the image of an edge yet"
        )

    def _section_homology_edge(self, label, edge):
        r"""
        Return a preimage of an edge in homology.

        This is a helper method for :meth:`_image_homology_matrix` of
        :meth:`section`. But usually this is not invoked since we compute the
        section with linear algebra in :meth:`_section_homology_matrix`.
        """
        return SurfaceMorphism._image_homology_edge(self.section(), label, edge)

    def __mul__(self, other):
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
                f"morphisms cannot be composed because domain of {self} is not compatible with codomain of {other}"
            )

        if isinstance(other, IdentityMorphism):
            return self

        if isinstance(self, IdentityMorphism):
            return other

        return CompositionMorphism._create_morphism(self, other)

    def push_vector_forward(self, tangent_vector):
        import warnings

        warnings.warn(
            "push_vector_forward() has been deprecated and will be removed in a future version of sage-flatsurf; call the morphism with the tangent vector instead, i.e., instead of morphism.push_vector_forward(t) use morphism(t)"
        )

        return self(tangent_vector)

    def pull_vector_back(self, tangent_vector):
        import warnings

        warnings.warn(
            "pull_vector_back() has been deprecated and will be removed in a future version of sage-flatsurf; call a section of morphism with the tangent vector instead, i.e., instead of morphism.pull_vector_back(t) use morphism.section()(t)"
        )

        return self.section()(tangent_vector)

    def __eq__(self, other):
        raise NotImplementedError("morphism does not implement __eq__ yet")

    def __hash__(self):
        return hash((self.domain(), self.codomain()))


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

    # Relay evaluating this morphism to asking the original morphism for a section.
    def _image_point(self, x):
        return self._morphism._section_point(x)

    def _image_homology(self, x):
        return self._morphism._section_homology(x)

    def _image_homology_gen(self, x):
        return self._morphism._section_homology_gen(x)

    def _image_saddle_connection(self, x):
        return self._morphism._section_saddle_connection(x)

    def _image_tangent_vector(self, x):
        return self._morphism._section_tangent_vector(x)

    def _image_homology_matrix(self):
        M = self._morphism._image_homology_matrix()
        return M.parent()(M.inverse())

    def _image_homology_edge(self, label, edge):
        for (l, e) in self.codomain().edges():
            if self._morphism._image_homology_edge(l, e) == [(1, label, edge)]:
                return [(1, l, e)]

        raise NotImplementedError

    def __eq__(self, other):
        if not isinstance(other, SectionMorphism):
            return False

        return self._morphism == other._morphism

    def _repr_type(self):
        return "Section"

    def _repr_defn(self):
        return f"Section of {self._morphism}"


class CompositionMorphism(SurfaceMorphism):
    # TODO: docstring
    def __init__(self, parent, lhs, rhs, category=None):
        # TODO: docstring
        super().__init__(parent, category=category)

        self._morphisms = []
        for morphism in [rhs, lhs]:
            if isinstance(morphism, CompositionMorphism):
                self._morphisms.extend(morphism._morphisms)
            else:
                self._morphisms.append(morphism)

    @classmethod
    def _create_morphism(cls, lhs, rhs):
        return super()._create_morphism(rhs.domain(), lhs.codomain(), lhs, rhs)

    def __call__(self, x):
        for morphism in self._morphisms:
            x = morphism(x)
        return x

    def _image_homology_matrix(self):
        from sage.all import prod

        return prod(
            morphism._image_homology_matrix() for morphism in self._morphisms[::-1]
        )

    def _image_homology_edge(self, label, edge):
        image = [(1, label, edge)]
        for morphism in self._morphisms:
            image = [
                (c * d, ll, ee)
                for (c, l, e) in image
                for (d, ll, ee) in morphism._image_homology_edge(l, e)
            ]

        return image

    def __eq__(self, other):
        if not isinstance(other, CompositionMorphism):
            return False

        return self._parent == other._parent and self._morphisms == other._morphisms

    def __hash__(self):
        return hash(tuple(self._morphisms))

    def _repr_type(self):
        return "Composite"

    def _repr_defn(self):
        return "  " + "\nthen\n  ".join(str(morphism) for morphism in self._morphisms)

    @cached_method
    def section(self):
        from sage.all import prod

        return prod(morphism.section() for morphism in self._morphisms)


class IdentityMorphism(SurfaceMorphism):
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

    """

    def __init__(self, parent, category=None):
        if parent.domain() is not parent.codomain():
            raise ValueError("domain and codomain of identity must be identical")

        super().__init__(parent, category=category)

    @classmethod
    def _create_morphism(cls, domain):
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

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this morphism, i.e., ``x`` itself.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: identity = S.erase_marked_points()
            sage: identity(S(0, 0))
            Vertex 0 of polygon 0

        """
        return x

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

    def __rmul__(self, other):
        r"""
        Concatenate ``other`` and this morphism.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: identity = S.erase_marked_points()
            sage: identity * identity
            Identity endomorphism of Translation Surface in H_2(2) built from 3 squares

        """
        super().__rmul__(other)
        return other


class SurfaceMorphism_factorization(SurfaceMorphism):
    def _image_homology(self, x, codomain=None):
        return self._factorization()._image_homology(x, codomain=codomain)


class NamedFactorizationMorphism(SurfaceMorphism_factorization):
    def __init__(self, parent, name, factorization):
        super().__init__(parent)

        self._name = name

    def _repr_type(self):
        return self._name


class TriangulationMorphism(SurfaceMorphism):
    r"""
    Morphism from a surface to its triangulation.

    EXAMPLES::

        sage: import flatsurf

        sage: G = SymmetricGroup(4)
        sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
        sage: S.triangulate()
        Triangulation morphism:
          From: Origami defined by r=(1,2,3,4) and u=(1,4,2,3)
          To:   Triangulation of Origami defined by r=(1,2,3,4) and u=(1,4,2,3)

    """

    def _image_edge(self, label, edge):
        return self.codomain()._triangulation(label)[1][edge]

    def _image_homology_edge(self, label, edge):
        r"""
        Implements :class:`SurfaceMorphism._image_homology_edge`.

        EXAMPLES::

            sage: import flatsurf

            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: triangulation = S.triangulate()
            sage: triangulation._image_homology_edge(1, 0)
            [(1, (1, 0), 0)]

        """
        return [(1, *self._image_edge(label, edge))]

    def _image_saddle_connection(self, connection):
        (label, edge) = connection.start()

        label, edge = self._image_edge(label, edge)

        from flatsurf.geometry.euclidean import ccw

        while (
            ccw(
                connection.direction().vector(),
                -self.codomain()
                .polygon(label)
                .edges()[(edge - 1) % len(self.codomain().polygon(label).edges())],
            )
            <= 0
        ):
            (label, edge) = self.codomain().opposite_edge(
                label, (edge - 1) % len(self.codomain().polygon(label).edges())
            )

        # TODO: This is extremely slow.
        from flatsurf.geometry.saddle_connection import SaddleConnection

        return SaddleConnection.from_vertex(
            self.codomain(), label, edge, connection.direction().vector()
        )

    def _section_saddle_connection(self, connection):
        (label, edge) = connection.start()

        domain_label = self.codomain()._reference_label(label)
        triangulation = self.codomain()._triangulation(domain_label)[1].inverse

        while (label, edge) not in triangulation:
            label, edge = self.codomain().opposite_edge(label, edge)
            edge = (edge + 1) % len(self.codomain().polygon(label).edges())

        # TODO: This is extremely slow.
        from flatsurf.geometry.saddle_connection import SaddleConnection

        return SaddleConnection.from_vertex(
            self.domain(),
            domain_label,
            triangulation[(label, edge)],
            connection.direction(),
        )

    def _image_point(self, p):
        r"""
        Implements :class:`SurfaceMorphism._image_point`.

        EXAMPLES::

            sage: import flatsurf

            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
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
        return "Triangulation"

    def __eq__(self, other):
        if not isinstance(other, TriangulationMorphism):
            return False

        return self.parent() == other.parent()


class DelaunayTriangulationMorphism(SurfaceMorphism):
    pass


class DelaunayDecompositionMorphism(SurfaceMorphism):
    def __eq__(self, other):
        if not isinstance(other, DelaunayDecompositionMorphism):
            return False

        return self.parent() == other.parent()

    def _repr_type(self):
        return "Delaunay Decomposition"


class GL2RMorphism(SurfaceMorphism):
    def __init__(self, parent, m, category=None):
        super().__init__(parent, category=category)

        from sage.all import matrix

        self._matrix = matrix(m, immutable=True)

    def change(self, domain=None, codomain=None, check=True):
        return type(self)._create_morphism(
            domain=domain or self.domain(),
            codomain=codomain or self.codomain(),
            m=self._matrix,
            category=self.category_for(),
        )

    def section(self):
        return self._create_morphism(self.codomain(), self.domain(), ~self._matrix)

    def _image_point(self, point):
        label, coordinates = point.representative()
        return self.codomain()(label, self._matrix * coordinates)

    def _image_saddle_connection(self, connection):
        if self._matrix.det() <= 0:
            raise NotImplementedError(
                "cannot compute the image of a saddle connection for this matrix yet"
            )

        from flatsurf.geometry.saddle_connection import SaddleConnection

        return SaddleConnection(
            surface=self.codomain(),
            start=connection.start(),
            end=connection.end(),
            holonomy=self._matrix * connection.holonomy(),
            end_holonomy=self._matrix * connection.end_holonomy(),
            check=False,
        )

    def _image_homology_edge(self, label, edge):
        return [(1, label, edge)]

    def _image_tangent_vector(self, t):
        return self.codomain().tangent_vector(
            t.polygon_label(), self._matrix * t.point(), self._matrix * t.vector()
        )

    def __eq__(self, other):
        if not isinstance(other, GL2RMorphism):
            return False

        return self.parent() == other.parent() and self._matrix == other._matrix

    def _repr_type(self):
        return "Linear"

    def _repr_defn(self):
        return repr(self._matrix)


class DelaunayDecompositionIsomorphism(SurfaceMorphism_factorization):
    pass
