r"""
EXAMPLES:

Usually, you do not interact with the types in this module directly but call
``minimal_cover()`` on a surface::

    sage: from flatsurf import polygons, similarity_surfaces
    sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5))

    sage: S.minimal_cover("translation")
    Minimal Translation Cover of Genus 0 Rational Cone Surface built from 2 right triangles
    sage: S.minimal_cover("half-translation")
    Minimal Half-Translation Cover of Genus 0 Rational Cone Surface built from 2 right triangles
    sage: S.minimal_cover("planar")
    Minimal Planar Cover of Genus 0 Rational Cone Surface built from 2 right triangles

"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2018-2019 W. Patrick Hooper
#                      2020-2022 Vincent Delecroix
#                      2021-2023 Julian Rüth
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
from sage.structure.element import parent
from sage.matrix.constructor import matrix
from sage.modules.free_module import VectorSpace

from flatsurf.geometry.surface import OrientedSimilaritySurface
from flatsurf.geometry.polygon import Polygon


def _is_finite(surface):
    r"""
    Return whether ``surface`` is a finite rational cone surface.
    """
    if not surface.is_finite_type():
        return False

    from flatsurf.geometry.categories import ConeSurfaces

    if surface in ConeSurfaces().Rational():
        return True

    surface = surface.reposition_polygons()

    for label in surface.labels():
        polygon = surface.polygon(label)

        for e in range(len(polygon.vertices())):
            m = surface.edge_matrix(label, e, projective=False)

            from flatsurf.geometry.euclidean import is_cosine_sine_of_rational

            if not is_cosine_sine_of_rational(m[0][0], m[0][1]):
                return False

    return True


# TODO: this would better be a morphism
class OrientedSimilaritySurfaceCover(OrientedSimilaritySurface):
    r"""
    A cover of surface.
    """
    def __init__(self, base_surface, base_ring=None, category=None):
        if base_ring is None:
            base_ring = base_surface.base_ring()

        if base_surface.is_mutable():
            if base_surface.is_finite_type():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                base_surface = MutableOrientedSimilaritySurface.from_surface(
                    base_surface
                )
                base_surface.set_immutable()
            else:
                raise ValueError(
                    "Can not construct MinimalPlanarCover of a surface that is mutable and infinite."
                )

        # TODO: we should drop this assumption and add tests
        if base_surface.is_with_boundary():
            raise TypeError("surface must be without boundary")

        self._base_surface = base_surface

        OrientedSimilaritySurface.__init__(
            self, base=base_ring, category=category
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces, translation_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.is_mutable()
            False

            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.is_mutable()
            False

            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S.is_mutable()
            False

        """
        return False

    def base_surface(self):
        r"""
        Return the surface of which ``self`` is a covering of.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces, translation_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.base_surface()
            Genus 0 Rational Cone Surface built from 2 right triangles

        """
        return self._base_surface

    def fiber_root(self, base_label):
        r"""
        Return a root label above the polygon with ``base_label`` in the base surface.

        Needs to be implemented in subclasses.
        """
        raise NotImplementedError

    def fiber_matrix(self, base_label, fiber, projective=True):
        r"""
        Return the matrix transformation (the section from the base surface to the cover).

        Needs to be implemented in subclasses.
        """
        raise NotImplementedError

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.roots()
            ((0, (1, 0)),)
            sage: len(S.labels())
            20

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.roots()
            ((0, (1, 0)),)
            sage: len(S.labels())
            10

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S.roots()
            ((0, (x, y) |-> (x, y)),)

        """
        return tuple(
            (base_label, self.fiber_root(base_label)) for base_label in self.base_surface().roots()
        )

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label`` in the cover.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        Warning: this method does not check that the input defines a valid polygon label.

        EXAMPLES::

            sage: from flatsurf import polygons, translation_surfaces, similarity_surfaces

            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.polygon((0, (1, 0)))
            Polygon(vertices=[(0, 0), (1, 0), (1/4*c^2 - 1/4, 1/4*c)])

            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: root = S.root()
            sage: S.polygon(root)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.polygon((0, (1, 0)))
            Polygon(vertices=[(0, 0), (1, 0), (1/4*c^2 - 1/4, 1/4*c)])

        """
        if not isinstance(label, tuple) or len(label) != 2:
            raise ValueError("invalid label")
        base_label, fiber = label
        m = self.fiber_matrix(base_label, fiber, projective=True)

        # TODO: change this when we have proper projective action of 3x3 matrices on polygons
        if m.det() < 0:
            raise NotImplementedError
        V3 = VectorSpace(self.base_ring(), 3)
        V2 = VectorSpace(self.base_ring(), 2)
        vertices_proj = [m * V3((x, y, 1)) for x, y in self.base_surface().polygon(base_label).vertices()]
        vertices_aff = [V2((x / z, y / z)) for x, y, z in vertices_proj]
        return Polygon(vertices=vertices_aff)

    def edge_matrix(self, p, e=None, projective=None):
        r"""
        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.edge_matrix((1, (1, 0)), 1, projective=True)
            [           -1             0 1/2*c^2 - 1/2]
            [            0            -1        -1/2*c]
            [            0             0             1]
            sage: S.edge_matrix((1, (1, 0)), 1, projective=False)
            [-1  0]
            [ 0 -1]
            sage: S.edge_matrix((1, (1, 0)), 1)
            doctest:warning
            ...
            UserWarning: the behavior of edge_matrix for similarity surfaces will change behavior in future version of sage-flatsurf; call with with projective=False to keep the old behavior or projective=True to switch to the forward compatible default version
            [-1  0]
            [ 0 -1]

        """
        if e is None:
            import warnings

            warnings.warn(
                "passing only a single tuple argument to edge_matrix() has been deprecated and will be deprecated in a future version of sage-flatsurf; pass the label and edge index as separate arguments instead"
            )
            p, e = p

        if projective is None:
            import warnings

            warnings.warn(
                "the behavior of edge_matrix for similarity surfaces will change behavior in future version of sage-flatsurf; call with with projective=False to keep the old behavior or projective=True to switch to the forward compatible default version"
            )

            projective = False

        if e < 0 or e >= len(self.polygon(p).vertices()):
            raise ValueError("invalid edge index for this polygon")

        base_label, fiber = p
        p1 = self.opposite_edge(p, e)
        if p1 is None:
            return None
        (base_label1, fiber1), _ = p1
        m = self.base_surface().edge_matrix(base_label, e, projective=projective)
        m0 = self.fiber_matrix(base_label, fiber, projective=projective)
        m1 = self.fiber_matrix(base_label1, fiber1, projective=projective)
        mm = m1 * m * m0.inverse()
        mm.set_immutable()
        return mm


class MinimalTranslationCover(OrientedSimilaritySurfaceCover):
    r"""
    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon, similarity_surfaces, polygons
        sage: from flatsurf.geometry.minimal_cover import MinimalTranslationCover
        sage: s = MutableOrientedSimilaritySurface(QQ)
        sage: s.add_polygon(Polygon(vertices=[(0,0),(5,0),(0,5)]))
        0
        sage: s.add_polygon(Polygon(vertices=[(0,0),(3,4),(-4,3)]))
        1
        sage: s.glue((0, 0), (1, 2))
        sage: s.glue((0, 1), (1, 1))
        sage: s.glue((0, 2), (1, 0))
        sage: s.set_immutable()
        sage: ss = s.minimal_cover("translation")
        sage: isinstance(ss, MinimalTranslationCover)
        True
        sage: ss.is_finite_type()
        True
        sage: len(ss.polygons())
        8
        sage: TestSuite(ss).run()

    The following is to test that unfolding is reasonably fast on the instances reported
    in https://github.com/flatsurf/sage-flatsurf/issues/47::

        sage: T = polygons.triangle(2, 13, 26)
        sage: S = similarity_surfaces.billiard(T)
        sage: S = S.minimal_cover("translation")
        sage: S
        Minimal Translation Cover of Genus 0 Rational Cone Surface built from 2 triangles

    TESTS::

        sage: from flatsurf.geometry.categories import TranslationSurfaces
        sage: S in TranslationSurfaces()
        True

    """

    def __init__(self, similarity_surface, category=None):
        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category &= TranslationSurfaces()

        category = category.WithoutBoundary()

        if _is_finite(similarity_surface):
            category = category.FiniteType()
            if similarity_surface.is_compact():
                category = category.Compact()
        else:
            category = category.InfiniteType()

        if similarity_surface.is_connected():
            category = category.Connected()

        OrientedSimilaritySurfaceCover.__init__(
            self, similarity_surface, category=category
        )


    def fiber_matrix(self, base_label, fiber, projective=True):
        r"""
        Return the matrix corresponding to the section of the covering from the
        polygon in the base surface with ``base_label`` to the polygon in the
        covering with ``fiber``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.fiber_matrix(0, (1, 0))
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: S.fiber_matrix(0, (1, -1))
            [ 1  1  0]
            [-1  1  0]
            [ 0  0  1]
            sage: S.fiber_matrix(0, (1, -1), projective=False)
            [ 1  1]
            [-1  1]

        """
        if not isinstance(fiber, tuple) or len(fiber) != 2:
            raise ValueError("invalid fiber {!r}".format(fiber))
        a, b = fiber
        a = self.base_ring().coerce(a)
        b = self.base_ring().coerce(b)
        if projective:
            return matrix(self.base_ring(), 3, [a, -b, 0, b, a, 0, 0, 0, 1])
        else:
            return matrix(self.base_ring(), 2, [a, -b, b, a])

    def fiber_root(self, base_label):
        r"""
        Return the root label of the fiber over ``base_label``.
        """
        return (self.base_ring().one(), self.base_ring().zero())

    @cached_method
    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.opposite_edge((0, (1, 0)), 0)
            ((1, (1, 0)), 2)

        """
        if len(label) != 2:
            raise ValueError("invalid label")
        pp, (a, b) = label  # this is the polygon m * base_surface.polygon(p)
        p2, e2 = self.base_surface().opposite_edge(pp, edge)
        m = self.base_surface().edge_matrix(p2, e2, projective=False)
        aa = a * m[0][0] - b * m[1][0]
        bb = b * m[0][0] + a * m[1][0]
        return ((p2, (aa, bb)), e2)

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S
            Minimal Translation Cover of Genus 0 Rational Cone Surface built from 2 right triangles

        """
        return f"Minimal Translation Cover of {repr(self.base_surface())}"

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.minimal_cover("translation")) == hash(S.minimal_cover("translation"))
            True

        """
        return hash(self.base_surface())

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: T = polygons.triangle(2, 13, 26)
            sage: S = similarity_surfaces.billiard(T)
            sage: S.minimal_cover("translation") == S.minimal_cover("translation")
            True

        ::

            sage: TT = polygons.triangle(2, 15, 26)
            sage: SS = similarity_surfaces.billiard(TT)
            sage: SS = SS.minimal_cover("translation")

            sage: S == SS
            False

        """
        if not isinstance(other, MinimalTranslationCover):
            return False

        return self.base_surface() == other.base_surface()


class MinimalHalfTranslationCover(OrientedSimilaritySurfaceCover):
    r"""
    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon, similarity_surfaces, polygons
        sage: from flatsurf.geometry.minimal_cover import MinimalHalfTranslationCover
        sage: s = MutableOrientedSimilaritySurface(QQ)
        sage: s.add_polygon(Polygon(vertices=[(0,0),(5,0),(0,5)]))
        0
        sage: s.add_polygon(Polygon(vertices=[(0,0),(3,4),(-4,3)]))
        1
        sage: s.glue((0, 0), (1, 2))
        sage: s.glue((0, 1), (1, 1))
        sage: s.glue((0, 2), (1, 0))
        sage: s.set_immutable()
        sage: ss = s.minimal_cover("half-translation")
        sage: isinstance(ss, MinimalHalfTranslationCover)
        True
        sage: ss.is_finite_type()
        True
        sage: len(ss.polygons())
        4
        sage: TestSuite(ss).run()

    The following is to test that unfolding is reasonably fast on the instances reported
    in https://github.com/flatsurf/sage-flatsurf/issues/47::

        sage: T = polygons.triangle(2, 13, 26)
        sage: S = similarity_surfaces.billiard(T)
        sage: S = S.minimal_cover("half-translation")
        sage: S
        Minimal Half-Translation Cover of Genus 0 Rational Cone Surface built from 2 triangles

    TESTS::

        sage: from flatsurf.geometry.categories import DilationSurfaces
        sage: S in DilationSurfaces()
        True

    """

    def __init__(self, similarity_surface, category=None):
        from flatsurf.geometry.categories import HalfTranslationSurfaces

        if category is None:
            category = HalfTranslationSurfaces()

        category &= HalfTranslationSurfaces()

        if not similarity_surface.is_compact():
            category = category.NotCompact()

        if _is_finite(similarity_surface):
            category = category.FiniteType()
            if similarity_surface.is_compact():
                category = category.Compact()
        else:
            category = category.InfiniteType()

        category = category.WithoutBoundary()

        if similarity_surface.is_connected():
            category = category.Connected()

        OrientedSimilaritySurfaceCover.__init__(
            self, similarity_surface, category=category
        )

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S
            Minimal Half-Translation Cover of Genus 0 Rational Cone Surface built from 2 right triangles

        """
        return f"Minimal Half-Translation Cover of {repr(self.base_surface())}"

    def fiber_matrix(self, base_label, fiber, projective=True):
        r"""
        Return the matrix corresponding to the section of the covering from the
        polygon in the base surface with ``base_label`` to the polygon in the
        covering with ``fiber``.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.fiber_matrix(0, (1, 1), projective=True)
            [ 1 -1  0]
            [ 1  1  0]
            [ 0  0  1]
            sage: S.fiber_matrix(0, (1, 1), projective=False)
            [ 1 -1]
            [ 1  1]

        """

        if not isinstance(fiber, tuple) or len(fiber) != 2:
            raise ValueError("invalid label {!r}".format(label))
        a, b = fiber
        if projective:
            return matrix(self.base_ring(), 3, [a, -b, 0, b, a, 0, 0, 0, 1])
        else:
            return matrix(self.base_ring(), 2, [a, -b, b, a])

    def fiber_root(self, base_label):
        return (self.base_ring().one(), self.base_ring().zero())

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.opposite_edge((0, (1, 0)), 0)
            ((1, (1, 0)), 2)

        """
        if len(label) != 2:
            raise ValueError("invalid label")
        pp, (a, b) = label  # this is the polygon m * ss.polygon(p)
        p2, e2 = self.base_surface().opposite_edge(pp, edge)
        m = self.base_surface().edge_matrix(pp, edge, projective=False)
        aa = a * m[0][0] + b * m[1][0]
        bb = b * m[0][0] - a * m[1][0]
        if aa > 0 or (aa == 0 and bb > 0):
            return ((p2, (aa, bb)), e2)
        else:
            return ((p2, (-aa, -bb)), e2)

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.minimal_cover("half-translation")) == hash(S.minimal_cover("half-translation"))
            True

        """
        return hash(self.base_surface())

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: T = polygons.triangle(2, 13, 26)
            sage: S = similarity_surfaces.billiard(T)
            sage: S.minimal_cover("half-translation") == S.minimal_cover("half-translation")
            True

        ::

            sage: TT = polygons.triangle(2, 15, 26)
            sage: SS = similarity_surfaces.billiard(TT)
            sage: SS = SS.minimal_cover("half-translation")

            sage: S == SS
            False

        """
        if not isinstance(other, MinimalHalfTranslationCover):
            return False

        return self.base_surface() == other.base_surface()


class MinimalPlanarCover(OrientedSimilaritySurfaceCover):
    r"""
    The minimal planar cover of a surface `S` is the smallest cover `C` so that
    the developing map from the universal cover `U` to the plane induces a well
    defined map from `C` to the plane. This is a translation surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s = translation_surfaces.square_torus()
        sage: from flatsurf.geometry.minimal_cover import MinimalPlanarCover
        sage: pc = s.minimal_cover("planar")
        sage: isinstance(pc, MinimalPlanarCover)
        True
        sage: pc.is_finite_type()
        False
        sage: sing = pc.singularity(pc.root(), 0, limit=4)
        doctest:warning
        ...
        UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
        doctest:warning
        ...
        UserWarning: limit has been deprecated as a keyword argument when creating points and will be removed without replacement in a future version of sage-flatsurf
        sage: len(sing.vertex_set())
        doctest:warning
        ...
        UserWarning: vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
        4
        sage: TestSuite(s).run()

    """

    def __init__(self, similarity_surface, base_label=None, category=None):
        if not similarity_surface.is_connected():
            raise NotImplementedError(
                "can only create a minimal planar cover of connected surfaces"
            )

        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category = category.WithoutBoundary()

        if similarity_surface.is_connected():
            category = category.Connected()

        if similarity_surface.is_compact():
            category = category.InfiniteType()
        elif not similarity_surface.is_finite_type():
            category = category.NotCompact()
        else:
            raise NotImplementedError("cannot determine category of planar cover of non-compact surface yet")

        OrientedSimilaritySurfaceCover.__init__(
            self, similarity_surface, category=category
        )

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S
            Minimal Planar Cover of Translation Surface in H_1(0) built from a square

        """
        return f"Minimal Planar Cover of {repr(self.base_surface())}"

    def fiber_matrix(self, base_label, fiber, projective=True):
        return fiber.matrix(projective)

    def fiber_root(self, base_label):
        from flatsurf.geometry.similarity import SimilarityGroup
        return SimilarityGroup(self.base_ring()).one()

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("planar")
            sage: root = S.root()
            sage: S.opposite_edge(root, 0)
            ((1, (x, y) |-> (x, y)), 2)

        """
        pp, m = label  # this is the polygon m * base_surface.polygon(pp)
        p2, e2 = self.base_surface().opposite_edge(pp, edge)
        me = self.base_surface().edge_transformation(pp, edge)
        mm = m * ~me
        return ((p2, mm), e2)

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.minimal_cover("planar")) == hash(S.minimal_cover("planar"))
            True

        """
        return hash((self.base_surface(), self.roots()))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: T = polygons.triangle(2, 13, 26)
            sage: S = similarity_surfaces.billiard(T)
            sage: S.minimal_cover("planar") == S.minimal_cover("planar")
            True

        """
        if not isinstance(other, MinimalPlanarCover):
            return False

        return self.base_surface() == other.base_surface() and self.roots() == other.roots()
