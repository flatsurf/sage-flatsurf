r"""
Morphisms between Surfaces

.. NOTE::

    One should think of these as maps between surfaces. However, they might
    often not be maps on the points of the surface but just on homology, or in
    some cases, they might not be meaningful as maps anywhere.

    Technically, these are at worst morphisms in the category of objects where
    being a morphism does not really have any mathematical meaning.

EXAMPLES:

We can use morphisms to follow a surface through a retriangulation process::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon()
    sage: morphism = S.subdivide_edges(3)
    sage: morphism = morphism.codomain().subdivide() * morphism
    sage: T = morphism.codomain()

    sage: morphism
    Generic morphism:
      From: Translation Surface in H_2(2) built from a regular octagon
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
#        Copyright (C) 2023 Julian Rüth
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
from sage.rings.ring import Ring
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.morphism import Morphism
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

    """

    def __init__(self):
        from sage.all import ZZ
        super().__init__(ZZ)

    def _repr_(self):
        return "The Unknown Ring"


class UnknownSurface(OrientedSimilaritySurface):
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

    """

    def _repr_(self):
        return "Unknown Surface"


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
    def __init__(self, domain, codomain, category=None):
        if domain is None:
            domain = UnknownSurface(UnknownRing())
        elif domain.is_mutable():
            domain = UnknownSurface(domain.base_ring())

        if codomain is None:
            codomain = UnknownSurface(UnknownRing())
        elif codomain.is_mutable():
            codomain = UnknownSurface(codomain.base_ring())

        if category is None:
            from sage.categories.all import Objects
            category = Objects()

        from sage.all import Hom
        parent = Hom(domain, codomain, category=category)
        super().__init__(parent)

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
        if name in ["__cached_methods"]:
            # Do not redirect __cached_methods to the surface since we do not
            # want to get the morphism and the surface cache mixed up.
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        try:
            attr = getattr(self.codomain(), name)
        except AttributeError:
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        import warnings
        warnings.warn(f"This methods returns a morphism instead of a surface. Use .codomain().{name} to access the surface instead of the morphism.")

        return attr

    def section(self):
        r"""
        Return a section of this morphism from its codomain to its domain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: morphism = S.apply_matrix(matrix([[2, 0], [0, 1]]), in_place=False)
            sage: morphism.section()
            Generic morphism:
              From: Translation Surface in H_1(0) built from a rectangle
              To:   Translation Surface in H_1(0) built from a square

        """
        return SectionMorphism(self)

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
            Saddle connection in direction (1, 0) with start data (0, 0) and end data (0, 2)
            sage: morphism(saddle_connection)
            Saddle connection in direction (1, 0) with start data (0, 0) and end data (0, 2)

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

        The image of a cohomology class::

            sage: from flatsurf import SimplicialCohomology
            sage: H = SimplicialCohomology(S)
            sage: a, b = H.gens()
            sage: a
            {B[(0, 1)]: 1}
            sage: morphism(a)
            {B[(0, 1)]: 1}

            sage: morphism(_)
            Traceback (most recent call last):
            ...
            ValueError: cohomology class must be defined over the domain of this morphism

        The image of a tangent vector::

            sage: t = S.tangent_vector(0, (0, 0), (1, 1))
            sage: morphism(t)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, 1)

        TESTS::

            sage: morphism(42)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot map Integer through this morphism yet

        """
        from flatsurf.geometry.surface_objects import SurfacePoint
        if isinstance(x, SurfacePoint):
            if x.parent() is not self.domain():
                raise ValueError("point must be in the domain of this morphism")
            image = self._image_point(x)
            assert image.parent() is self.codomain()
            return image

        from flatsurf.geometry.homology import SimplicialHomologyClass
        if isinstance(x, SimplicialHomologyClass):
            if x.parent().surface() is not self.domain():
                raise ValueError("homology class must be defined over the domain of this morphism")
            image = self._image_homology(x)
            assert image.parent().surface() is self.codomain()
            return image

        from flatsurf.geometry.cohomology import SimplicialCohomologyClass
        if isinstance(x, SimplicialCohomologyClass):
            if x.parent().surface() is not self.domain():
                raise ValueError("cohomology class must be defined over the domain of this morphism")
            image = self._image_cohomology(x)
            assert image.parent().surface() is self.codomain()
            return image

        from flatsurf.geometry.saddle_connection import SaddleConnection
        if isinstance(x, SaddleConnection):
            if x.surface() is not self.domain():
                raise ValueError("saddle connection must be in the domain of this morphism")
            image = self._image_saddle_connection(x)
            assert image.surface() is self.codomain()
            return image

        from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
        if isinstance(x, SimilaritySurfaceTangentVector):
            if x.surface() is not self.domain():
                raise ValueError("tangent vector must be on the domain of this morphism")
            image = self._image_tangent_vector(x)
            assert image.surface() is self.codomain()
            return image

        raise NotImplementedError(f"cannot map {type(x).__name__} through this morphism yet")

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

            TODO: Add an example of such a morphism.

        """
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of a point yet")

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
            sage: c = SaddleConnection(S, (0, 0), (1, 1))
            sage: c
            Saddle connection in direction (1, 1) with start data (0, 0) and end data (0, 2)

            sage: morphism(c)
            Saddle connection in direction (1, 1/2) with start data (0, 0) and end data (0, 2)

        Not all morphisms are meaningful on the level of saddle connections::

            sage: morphism = S.subdivide()
            sage: morphism(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: a SubdivideMorphism cannot compute the image of a saddle connection yet

        """
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of a saddle connection yet")

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

        Not all morphisms are meaningful on the level of saddle connections::

            TODO: Add an example of such a morphism

        """
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of a tangent vector yet")

    def _image_homology(self, γ):
        # TODO: docstring

        from flatsurf.geometry.homology import SimplicialHomology
        codomain_homology = SimplicialHomology(self.codomain())

        from sage.all import vector
        image = self._image_homology_matrix() * vector(γ._homology())

        homology, to_chain, to_homology = codomain_homology._homology()

        image = sum(coefficient * gen for (coefficient, gen) in zip(image, homology.gens()))

        return codomain_homology(to_chain(image))

    def _image_edge(self, label, edge):
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of an edge yet")

    @cached_method
    def _image_homology_matrix(self):
        # TODO: docstring
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

    def _image_homology_gen(self, gen):
        from flatsurf.geometry.homology import SimplicialHomology
        codomain_homology = SimplicialHomology(self.codomain())

        chain = gen._chain
        image = codomain_homology.zero()
        for label, edge in chain.support():
            coefficient = chain[(label, edge)]
            assert coefficient
            for step in self._image_edge(label, edge):
                image += coefficient * image.parent()(step)

        return codomain_homology(image)

    def _image_cohomology(self, f):
        # TODO: docstring
        # TODO: Is this really a meaningful operation or should we restrict to computing sections of cohomology? Here and in __call__.
        from sage.all import ZZ, matrix

        from flatsurf.geometry.homology import SimplicialHomology

        domain_homology = SimplicialHomology(self.domain())
        codomain_homology = SimplicialHomology(self.codomain())

        M = matrix(ZZ, [
            [domain_homology_gen_image.coefficient(codomain_homology_gen) for codomain_homology_gen in codomain_homology.gens()]
            for domain_homology_gen_image in [self(domain_homology_gen) for domain_homology_gen in domain_homology.gens()]
        ])

        from sage.all import vector
        values = M.solve_right(vector([f(gen) for gen in domain_homology.gens()]))

        from flatsurf.geometry.cohomology import SimplicialCohomology
        codomain_cohomology = SimplicialCohomology(self.codomain())
        return codomain_cohomology({gen: value for (gen, value) in zip(codomain_homology.gens(), values)})

    def __mul__(self, other):
        # TODO: docstring
        return CompositionMorphism(self, other)

    def push_vector_forward(self, tangent_vector):
        import warnings
        warnings.warn("push_vector_forward() has been deprecated and will be removed in a future version of sage-flatsurf; call the morphism with the tangent vector instead, i.e., instead of morphism.push_vector_forward(t) use morphism(t)")

        return self(tangent_vector)

    def pull_vector_back(self, tangent_vector):
        import warnings
        warnings.warn("pull_vector_back() has been deprecated and will be removed in a future version of sage-flatsurf; call a section of morphism with the tangent vector instead, i.e., instead of morphism.pull_vector_back(t) use morphism.section()(t)")

        return self.section()(tangent_vector)


class IdentityMorphism(SurfaceMorphism):
    # TODO: docstring
    def __init__(self, domain):
        super().__init__(domain, domain)

    def __call__(self, x):
        return x

    def _image_homology(self, γ):
        return γ

    def _image_edge(self, label, edge):
        return [(label, edge)]


class SectionMorphism(SurfaceMorphism):
    def __init__(self, morphism):
        self._morphism = morphism
        super().__init__(morphism.codomain(), morphism.domain())

    def _image_homology_matrix(self):
        M = self._morphism._image_homology_matrix()
        return M.parent()(M.inverse())

    def _image_edge(self, label, edge):
        for (l, e) in self.codomain().edges():
            if self._morphism._image_edge(l, e) == (label, edge):
                return (l, e)

        raise NotImplementedError


class CompositionMorphism(SurfaceMorphism):
    # TODO: docstring
    def __init__(self, lhs, rhs):
        # TODO: docstring
        super().__init__(rhs.domain(), lhs.codomain())

        self._lhs = lhs
        self._rhs = rhs

    def __call__(self, x):
        # TODO: docstring
        return self._lhs(self._rhs(x))

    def _image_homology_matrix(self):
        return self._lhs._image_homology_matrix() * self._rhs._image_homology_matrix()

    def _image_edge(self, label, edge):
        return self._lhs._image_edge(*self._rhs._image_edge(label, edge))


class SubdivideMorphism(SurfaceMorphism):
    # TODO: docstring
    def _image_point(self, p):
        # TODO: docstring
        from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
        to_pyflatsurf = FlatTriangulationConversion.to_pyflatsurf(domain=self.codomain())

        label = next(iter(p.labels()))
        coordinates = next(iter(p.coordinates(label))) - self.domain().polygon(label).vertex(0)

        from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
        to_pyflatsurf_vector = VectorSpaceConversion.to_pyflatsurf(coordinates.parent())

        def ccw(v, w):
            r"""
            Return whether v->w describe a non-clockwise turn.
            """
            return v[0] * w[1] >= w[0] * v[1]

        initial_edge = ((label, 0), 0)
        mid_edge = ((label, len(self.codomain().polygons()) // len(self.domain().polygons()) - 1), 1)

        face = mid_edge if ccw(self.codomain().polygon(mid_edge[0]).edge(mid_edge[1]), coordinates) else initial_edge
        face = to_pyflatsurf(face)
        coordinates = to_pyflatsurf_vector(coordinates)

        import pyflatsurf

        p = pyflatsurf.flatsurf.Point[type(to_pyflatsurf.codomain())](to_pyflatsurf.codomain(), face, coordinates)

        return to_pyflatsurf.section(p)

    def _image_homology(self, γ):
        # TODO: homology should allow to write this code more naturally somehow
        # TODO: docstring
        from flatsurf.geometry.homology import SimplicialHomology

        H = SimplicialHomology(self.codomain())
        C = H.chain_module(dimension=1)

        image = H()

        for gen in γ.parent().gens():
            for step in gen._path():
                codomain_gen = (step, 0)
                sgn = 1
                if codomain_gen not in H.simplices():
                    sgn = -1
                    codomain_gen = self.codomain().opposite_edge(*codomain_gen)
                image += sgn * γ.coefficient(gen) * H(C(codomain_gen))

        return image


class SubdivideEdgesMorphism(SurfaceMorphism):
    # TODO: docstring

    def __init__(self, domain, codomain, parts):
        super().__init__(domain, codomain)
        self._parts = parts

    def _image_point(self, p):
        r"""
        Return the image of ``p`` in the codomain of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: morphism = S.subdivide_edges(2)

            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: p = SurfacePoint(S, 0, (1, 1))

            sage: morphism._image_point(p)
            Point (1, 1) of polygon 0

        """
        from flatsurf.geometry.surface_objects import SurfacePoint

        if not p.surface() == self.domain():
            raise ValueError

        label = next(iter(p.labels()))
        return SurfacePoint(self.codomain(), label, next(iter(p.coordinates(label))))

    # def _image_homology(self, γ):
    #     # TODO: homology should allow to write this code more naturally somehow
    #     # TODO: docstring
    #     from flatsurf.geometry.homology import SimplicialHomology

    #     H = SimplicialHomology(self.codomain())
    #     C = H.chain_module(dimension=1)

    #     image = H()

    #     for gen in γ.parent().gens():
    #         for step in gen._path():
    #             for p in range(self._parts):
    #                 codomain_gen = (step[0], step[1] * self._parts + p)
    #                 sgn = 1
    #                 if codomain_gen not in H.simplices():
    #                     sgn = -1
    #                     codomain_gen = self.codomain().opposite_edge(*codomain_gen)
    #                 image += sgn * γ.coefficient(gen) * H(C(codomain_gen))

    #     return image


class TrianglesFlipMorphism(SurfaceMorphism):
    # TODO: docstring
    def __init__(self, domain, codomain, flip_sequence):
        # TODO: docstring
        super().__init__(domain, codomain)
        assert not self.domain().is_mutable()
        self._flip_sequence = flip_sequence

    @cached_method
    def _flip_morphism(self):
        domains = [self.domain()]

        for flip in self._flip_sequence:
            from flatsurf.geometry.surface import Surface_dict
            codomain = Surface_dict(surface=domains[-1], mutable=True)
            codomain.flip(*flip)
            codomain.set_immutable()
            domains.append(codomain)

        assert domains[-1] == self.codomain()
        domains.pop()
        domains.append(self.codomain())

        # TODO: Get rid of this trivial step.
        morphism = IdentityMorphism(self.domain())
        for i, flip in enumerate(self._flip_sequence):
            morphism = TriangleFlipMorphism(domains[i], domains[i + 1], flip) * morphism

        return morphism

    def _image_homology(self,  γ):
        return self._flip_morphism()._image_homology(γ)

    def _image_homology_matrix(self):
        return self._flip_morphism()._image_homology_matrix()

    def _image_point(self, p):
        return self._flip_morphism()(p)


class TriangleFlipMorphism(SurfaceMorphism):
    # TODO: docstring
    def __init__(self, domain, codomain, flip):
        # TODO: docstring
        super().__init__(domain, codomain)
        self._flip = flip

    def _image_edge(self, label, edge):
        opposite_flip = self.domain().opposite_edge(*self._flip)

        def find_in_codomain(vector):
            candidates = [(l, e) for l in [self._flip[0], opposite_flip[0]] for e in [0, 1, 2]]
            candidates = [(l, e) for (l, e) in candidates if self.codomain().polygon(l).edge(e) == vector]
            assert candidates
            if len(candidates) > 1:
                raise NotImplementedError("cannot break ties on such a non-translation surface yet")

            return candidates[0]

        if label == self._flip[0]:
            if edge == self._flip[1]:
                return self._image_edge(*self.domain().opposite_edge(label, (edge + 1) % 3)) + self._image_edge(*self.domain().opposite_edge(label, (edge + 2) % 3))

            return [find_in_codomain(self.domain().polygon(label).edge(edge))]
        if label == opposite_flip[0]:
            if edge == opposite_flip[1]:
                return self._image_edge(*self.domain().opposite_edge(label, (edge + 1) % 3)) + self._image_edge(*self.domain().opposite_edge(label, (edge + 2) % 3))

            return [find_in_codomain(self.domain().polygon(label).edge(edge))]

        return [(label, edge)]

    def _image_point(self, p):
        # TODO: This is a dumb way to do this (and wrong for non-translation surfaces?)

        from flatsurf.geometry.surface_objects import SurfacePoint
        affected = {self._flip[0], self.domain().opposite_edge(*self._flip)[0]}

        for label in p.labels():
            if label not in affected:
                return SurfacePoint(self.codomain(), label, next(iter(p.coordinates(label))))

        for label in p.labels():
            polygon = self.domain().polygon(label)
            for edge in range(polygon.num_edges()):
                cross_label, cross_edge = self.domain().opposite_edge(label, edge)
                if cross_label in affected:
                    continue

                Δ = next(iter(p.coordinates(label))) - polygon.vertex(edge)
                recross_label, recross_edge = self.codomain().opposite_edge(cross_label, cross_edge)
                recross_polygon = self.codomain().polygon(recross_label)
                recross_position = recross_polygon.vertex(recross_edge) + Δ
                if not recross_polygon.get_point_position(recross_position).is_inside():
                    continue

                return SurfacePoint(self.codomain(), recross_label, recross_position)

        raise NotImplementedError


class TriangulationMorphism(SurfaceMorphism):
    def _image_edge(self, label, edge):
        return [self.codomain()._triangulation(label)[1][edge]]

    def _image_saddle_connection(self, connection):
        (label, edge) = connection.start_data()

        (label, edge) = self._image_edge(label, edge)[0]

        from flatsurf.geometry.euclidean import ccw
        while ccw(connection.direction(), -self.codomain().polygon(label).edges()[(edge - 1) % len(self.codomain().polygon(label).edges())]) <= 0:
            (label, edge) = self.codomain().opposite_edge(label, (edge - 1) % len(self.codomain().polygon(label).edges()))

        # TODO: This is extremely slow.
        from flatsurf.geometry.saddle_connection import SaddleConnection
        return SaddleConnection(self.codomain(), (label, edge), direction=connection.direction())


class DelaunayDecompositionMorphism(SurfaceMorphism):
    def __repr__(self):
        # TODO: docstring
        return f"Delaunay cell decomposition of {self.domain()}"


class GL2RMorphism(SurfaceMorphism):
    def __init__(self, domain, codomain, m):
        super().__init__(domain, codomain)

        from sage.all import matrix
        self._matrix = matrix(m, immutable=True)

    def _image_point(self, point):
        label, coordinates = point.representative()
        return self.codomain()(label, self._matrix * coordinates)

    def _image_saddle_connection(self, connection):
        if self._matrix.det() <= 0:
            raise NotImplementedError("cannot compute the image of a saddle connection for this matrix yet")

        from flatsurf.geometry.saddle_connection import SaddleConnection
        return SaddleConnection(
            surface=self.codomain(),
            start_data=connection.start_data(),
            direction=self._matrix * connection.direction(),
            end_data=connection.end_data(),
            holonomy=self._matrix * connection.holonomy(),
            end_holonomy=self._matrix * connection.end_holonomy(),
            check=False)

    def _image_edge(self, label, edge):
        return [(label, edge)]

    def _image_tangent_vector(self, t):
        return self.codomain().tangent_vector(
            t.polygon_label(),
            self._matrix * t.point(),
            self._matrix * t.vector())
