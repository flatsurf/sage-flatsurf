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
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.surface import OrientedSimilaritySurface


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
            m = surface.edge_matrix(label, e)

            from flatsurf.geometry.euclidean import is_cosine_sine_of_rational

            if not is_cosine_sine_of_rational(m[0][0], m[0][1]):
                return False

    return True


class MinimalTranslationCover(OrientedSimilaritySurface):
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
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite_type():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                similarity_surface = MutableOrientedSimilaritySurface.from_surface(
                    similarity_surface
                )
            else:
                raise NotImplementedError(
                    "can not construct MinimalTranslationCover of a surface that is mutable and infinite"
                )

        if similarity_surface.is_with_boundary():
            raise TypeError("surface must be without boundary")

        self._ss = similarity_surface

        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category &= TranslationSurfaces()

        category = category.WithoutBoundary()

        if _is_finite(self._ss):
            category = category.FiniteType()
        else:
            category = category.InfiniteType()

        if similarity_surface.is_connected():
            category = category.Connected()

        if self.is_compact():
            category = category.Compact()

        OrientedSimilaritySurface.__init__(
            self, self._ss.base_ring(), category=category
        )

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
            ((0, 1, 0),)

        """
        self._F = self._ss.base_ring()
        return tuple(
            (label, self._F.one(), self._F.zero()) for label in self._ss.roots()
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.is_mutable()
            False

        """
        return False

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().minimal_cover("translation")
            sage: S.is_compact()
            False

        ::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.is_compact()
            True

        """
        if not self._ss.is_compact():
            return False

        if not self._ss.is_rational_surface():
            return False

        return True

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("translation")
            sage: S.polygon((0, 1, 0))
            Polygon(vertices=[(0, 0), (1, 0), (1/4*c^2 - 1/4, 1/4*c)])

        """
        if not isinstance(label, tuple) or len(label) != 3:
            raise ValueError("invalid label {!r}".format(label))
        return matrix([[label[1], -label[2]], [label[2], label[1]]]) * self._ss.polygon(
            label[0]
        )

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
            sage: S.opposite_edge((0, 1, 0), 0)
            ((1, 1, 0), 2)

        """
        pp, a, b = label  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, edge)
        m = self._ss.edge_matrix(p2, e2)
        aa = a * m[0][0] - b * m[1][0]
        bb = b * m[0][0] + a * m[1][0]
        return ((p2, aa, bb), e2)

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
        return f"Minimal Translation Cover of {repr(self._ss)}"

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
        return hash(self._ss)

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

        return self._ss == other._ss


class MinimalHalfTranslationCover(OrientedSimilaritySurface):
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
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite_type():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                self._ss = MutableOrientedSimilaritySurface.from_surface(
                    similarity_surface
                )
                self._ss.set_immutable()
            else:
                raise ValueError(
                    "Can not construct MinimalTranslationCover of a surface that is mutable and infinite."
                )
        else:
            self._ss = similarity_surface

        if similarity_surface.is_with_boundary():
            raise TypeError(
                "can only build translation cover of surfaces without boundary"
            )

        from flatsurf.geometry.categories import HalfTranslationSurfaces

        if category is None:
            category = HalfTranslationSurfaces()

        category &= HalfTranslationSurfaces()

        if _is_finite(self._ss):
            category = category.FiniteType()
        else:
            category = category.InfiniteType()

        category = category.WithoutBoundary()

        if similarity_surface.is_connected():
            category = category.Connected()

        OrientedSimilaritySurface.__init__(
            self, self._ss.base_ring(), category=category
        )

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.roots()
            ((0, 1, 0),)

        """
        self._F = self._ss.base_ring()
        return tuple(
            (label, self._F.one(), self._F.zero()) for label in self._ss.roots()
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.is_mutable()
            False

        """
        return False

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
        return f"Minimal Half-Translation Cover of {repr(self._ss)}"

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.polygon((0, 1, 0))
            Polygon(vertices=[(0, 0), (1, 0), (1/4*c^2 - 1/4, 1/4*c)])

        """
        if not isinstance(label, tuple) or len(label) != 3:
            raise ValueError("invalid label {!r}".format(label))
        return matrix([[label[1], -label[2]], [label[2], label[1]]]) * self._ss.polygon(
            label[0]
        )

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf import polygons, similarity_surfaces
            sage: S = similarity_surfaces.billiard(polygons.triangle(2, 3, 5)).minimal_cover("half-translation")
            sage: S.opposite_edge((0, 1, 0), 0)
            ((1, 1, 0), 2)

        """
        pp, a, b = label  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, edge)
        m = self._ss.edge_matrix(pp, edge)
        aa = a * m[0][0] + b * m[1][0]
        bb = b * m[0][0] - a * m[1][0]
        if aa > 0 or (aa == 0 and bb > 0):
            return ((p2, aa, bb), e2)
        else:
            return ((p2, -aa, -bb), e2)

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
        return hash(self._ss)

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

        return self._ss == other._ss


class MinimalPlanarCover(OrientedSimilaritySurface):
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
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite_type():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

                self._ss = MutableOrientedSimilaritySurface.from_surface(
                    similarity_surface
                )
                self._ss.set_immutable()
            else:
                raise ValueError(
                    "Can not construct MinimalPlanarCover of a surface that is mutable and infinite."
                )
        else:
            self._ss = similarity_surface

        if base_label is not None:
            import warnings

            warnings.warn(
                "the keyword argument base_label of a minimal planar cover is ignored and will be removed in a future version of sage-flatsurf; it had no effect in previous versions of sage-flatsurf"
            )

        if not self._ss.is_connected():
            raise NotImplementedError(
                "can only create a minimal planar cover of connected surfaces"
            )

        # The similarity group containing edge identifications.
        self._sg = self._ss.edge_transformation(self._ss.root(), 0).parent()
        self._root = (self._ss.root(), self._sg.one())

        if similarity_surface.is_with_boundary():
            raise TypeError(
                "can only build translation cover of surfaces without boundary"
            )

        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category &= TranslationSurfaces().InfiniteType()

        category = category.WithoutBoundary()

        if similarity_surface.is_connected():
            category = category.Connected()

        OrientedSimilaritySurface.__init__(
            self, self._ss.base_ring(), category=category
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
        return f"Minimal Planar Cover of {repr(self._ss)}"

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space, i.e.,
        return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S.is_compact()
            False

        """
        return False

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S.roots()
            ((0, (x, y) |-> (x, y)),)

        """
        return (self._root,)

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: S.is_mutable()
            False

        """
        return False

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().minimal_cover("planar")
            sage: root = S.root()
            sage: S.polygon(root)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        if not isinstance(label, tuple) or len(label) != 2:
            raise ValueError("invalid label {!r}".format(label))
        return label[1](self._ss.polygon(label[0]))

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
        pp, m = label  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, edge)
        me = self._ss.edge_transformation(pp, edge)
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
        return hash((self._ss, self._root))

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

        return self._ss == other._ss and self._root == other._root
