r"""
Triangulations, Delaunay triangulations, and Delaunay decompositions of
infinite surfaces.

EXAMPLES:

Typically, you don't need to create these surfaces directly, they are created
by invoking methods on the underlying surfaces::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.infinite_staircase()
    sage: S.triangulate()
    Triangulation of The infinite staircase

    sage: S.delaunay_triangulation()
    Delaunay triangulation of The infinite staircase

    sage: S.delaunay_decomposition()
    Delaunay cell decomposition of The infinite staircase

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
#                          2023 Julian Rüth
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

from flatsurf.geometry.surface import (
    MutableOrientedSimilaritySurface_base,
    OrientedSimilaritySurface,
    Labels,
)


class LazyTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    A triangulated surface whose structure is computed on demand.

    Used to triangulate surfaces when in-place-triangulation is not requested.
    In particular, this is used to triangulate infinite type surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()
        sage: S = S.triangulate()

    TESTS::

        sage: from flatsurf.geometry.lazy import LazyTriangulatedSurface
        sage: isinstance(S, LazyTriangulatedSurface)
        True
        sage: TestSuite(S).run()  # long time (1s)

    """

    def __init__(self, surface, labels=None, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyTriangulatedSurface; use relabel() to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyTriangulatedSurface()"
                )

        if surface.is_mutable():
            raise ValueError("Surface must be immutable.")

        if labels is not None:
            labels = set(labels)

        if surface.is_finite_type():
            if labels == set(surface.labels()):
                labels = None

        self._reference = surface
        self._triangulated_reference_labels = labels

        OrientedSimilaritySurface.__init__(
            self,
            surface.base_ring(),
            category=category or self._reference.category(),
        )

    def _is_triangulated(self, reference_label):
        r"""
        Return whether the ``reference_label`` of the reference surface is
        explicitly triangulated in this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S = S.triangulate()
            sage: S._is_triangulated(0)
            True
            sage: S._is_triangulated(1)
            True

            sage: S = S.triangulate(label=0)
            sage: S._is_triangulated(0)
            True
            sage: S._is_triangulated(1)
            False

        """
        if self._triangulated_reference_labels is None:
            return True

        return reference_label in self._triangulated_reference_labels

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
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
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.roots()
            ((0, 0),)

        """
        return tuple(self._image(root)[0].root() for root in self._reference.roots())

    @cached_method
    def _image(self, reference_label):
        r"""
        Return a triangulation of the ``reference_label`` in the underlying
        (typically non-triangulated) reference surface.

        If the ``reference_label`` is not being triangulated, the return a
        surface just consisting of this polygon.

        INPUT:

        - ``reference_label`` -- a polygon label in the reference surface that
          we are triangulating.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S._image(0)
            (Translation Surface with boundary built from 2 isosceles triangles,
             bidict({0: ((0, 0), 0), 1: ((0, 0), 1), 2: ((0, 1), 1), 3: ((0, 1), 2)}))

            sage: S = translation_surfaces.infinite_staircase().triangulate(label=1)
            sage: S._image(0)
            (Translation Surface with boundary built from a square,
             bidict({0: (0, 0), 1: (0, 1), 2: (0, 2), 3: (0, 3)}))

        """
        reference_polygon = self._reference.polygon(reference_label)

        if not self._is_triangulated(reference_label):
            from flatsurf import MutableOrientedSimilaritySurface

            triangulation = MutableOrientedSimilaritySurface(
                self._reference.base_ring()
            )
            triangulation.add_polygon(reference_polygon)

            from bidict import bidict

            return triangulation, bidict(
                {e: (reference_label, e) for e in range(len(reference_polygon.edges()))}
            )

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        return MutableOrientedSimilaritySurface._triangulate(
            self._reference, reference_label
        )

    def _reference_label(self, label):
        r"""
        Return the label of the underlying (untriangulated) reference surface
        which led to the creation of the polygon with ``label`` in this
        triangulation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S._reference_label((0, 0))
            0

        """
        if label in self._reference.labels():
            if not self._is_triangulated(label):
                return label
            if len(self._reference.polygon(label).vertices()) == 3:
                return label

        if not isinstance(label, tuple):
            raise KeyError(label)

        if len(label) != 2:
            raise KeyError(label)

        if label[0] not in self._reference.labels():
            raise KeyError(label)

        if label not in self._image(label[0])[0].labels():
            raise KeyError(label)

        return label[0]

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.polygon((0, 0))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        reference_label = self._reference_label(label)

        if not self._is_triangulated(reference_label):
            return self._reference.polygon(reference_label)

        triangulation, _ = self._image(reference_label)

        return triangulation.polygon(label)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.opposite_edge((0, 0), 0)
            ((1, 1), 1)

        """
        reference_label = self._reference_label(label)

        triangulation, outer_edges = self._image(reference_label)

        if (label, edge) in outer_edges.values():
            # pylint does not understand the bidict return type, so we disable a failing check here
            # pylint: disable=unsubscriptable-object
            reference_edge = outer_edges.inverse[(label, edge)]
            (
                opposite_reference_label,
                opposite_reference_edge,
            ) = self._reference.opposite_edge(reference_label, reference_edge)
            opposite_triangulation, opposite_outer_edges = self._image(
                opposite_reference_label
            )
            return opposite_outer_edges[opposite_reference_edge]

        return triangulation.opposite_edge(label, edge)

    def is_triangulated(self, limit=None):
        r"""
        Return whether this surface is triangulated.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.is_triangulated()
            True

            sage: S = translation_surfaces.infinite_staircase().triangulate(label=0)
            sage: S.is_triangulated()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot decide whether this (potentially infinite type) surface is triangulated

        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        if self._triangulated_reference_labels is None:
            return True

        return super().is_triangulated(limit=limit)

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.triangulate()) == hash(S.triangulate())
            True

        """
        return hash(self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.triangulate() == S.triangulate()
            True

        """
        if not isinstance(other, LazyTriangulatedSurface):
            return False

        return (
            self._reference == other._reference
            and self._triangulated_reference_labels
            == other._triangulated_reference_labels
        )

    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S.labels()
            ((0, 0), (1, 1), (-1, 1), (0, 1), (1, 0), (2, 0), (-1, 0), (-2, 0), (2, 1), (3, 1), (-2, 1), (-3, 1), (3, 0), (4, 0), (-3, 0), (-4, 0), …)

        """
        return TriangulationLabels(self, finite=self._reference.is_finite_type())

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: S
            Triangulation of The infinite staircase

        """
        if self._triangulated_reference_labels is not None:
            return f"Partial Triangulation of {self._reference!r}"

        return f"Triangulation of {self._reference!r}"


class TriangulationLabels(Labels):
    r"""
    The labels of a triangulation of a (possibly infinite) surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase().triangulate()
        sage: labels = S.labels()

    TESTS::

        sage: from flatsurf.geometry.lazy import TriangulationLabels
        sage: isinstance(labels, TriangulationLabels)
        True

    """

    def __contains__(self, label):
        r"""
        Return whether ``label`` is present as a label in this triangulation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().triangulate()
            sage: labels = S.labels()
            sage: 0 in labels
            False
            sage: (0, 0) in labels
            True
            sage: (0, 1) in labels
            True
            sage: (0, 2) in labels
            False

        """
        try:
            self._surface._reference_label(label)
        except KeyError:
            return False

        return True


class LazyOrientedSimilaritySurface(OrientedSimilaritySurface):
    r"""
    A surface that forwards all queries to an underlying reference surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()
        sage: T = matrix([[2, 0], [0, 1]]) * S

        sage: from flatsurf.geometry.lazy import LazyOrientedSimilaritySurface
        sage: isinstance(T, LazyOrientedSimilaritySurface)
        True

    """

    def __init__(self, base_ring, reference, category=None):
        super().__init__(base_ring, category=category or reference.category())
        self._reference = reference

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S

            sage: S.is_compact()
            True

        """
        return self._reference.is_compact()

    def is_translation_surface(self, positive=True):
        r"""
        Return whether this surface is a translation surface, i.e., glued
        edges can be transformed into each other by translations.

        This implements
        :meth:`flatsurf.geometry.categories.similarity_surfaces.SimilaritySurfaces.ParentMethods.is_translation_surface`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S

            sage: S.is_translation_surface()
            True

        """
        return self._reference.is_translation_surface(positive=positive)

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.roots()
            (0,)

        """
        return self._reference.roots()

    def labels(self):
        r"""
        Return the labels of this surface which are just the labels of the
        underlying reference surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.labels()
            (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, …)

        """
        return self._reference.labels()

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: T = matrix([[2, 0], [0, 1]]) * S
            sage: T.polygon(0)
            Polygon(vertices=[(0, 0), (2, 0), (2, 1), (0, 1)])

        """
        return self._reference.polygon(label)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: T = matrix([[2, 0], [0, 1]]) * S
            sage: T.opposite_edge(0, 0)
            (1, 2)

        """
        return self._reference.opposite_edge(label, edge)

    def is_mutable(self):
        r"""
        Return whether this surface could be changing.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: T = matrix([[2, 0], [0, 1]]) * S
            sage: T.is_mutable()
            False

        """
        return self._reference.is_mutable()


class GL2RImageSurface(LazyOrientedSimilaritySurface):
    r"""
    The GL(2,R) image of an oriented similarity surface obtained by applying a
    matrix to each polygon while keeping the gluings intact.

    EXAMPLE::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.octagon_and_squares()
        sage: r = matrix(ZZ,[[0, 1], [1, 0]])
        sage: SS = r * S

        sage: S.canonicalize() == SS.canonicalize()
        True

    TESTS::

        sage: TestSuite(SS).run()

        sage: from flatsurf.geometry.lazy import GL2RImageSurface
        sage: isinstance(SS, GL2RImageSurface)
        True

    """

    def __init__(self, reference, m, category=None):
        if reference.is_mutable():
            if not reference.is_finite_type():
                raise NotImplementedError(
                    "cannot apply matrix to mutable surface of infinite type"
                )

            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

            reference = MutableOrientedSimilaritySurface.from_surface(reference)

        self._reference = reference

        from sage.structure.element import get_coercion_model

        cm = get_coercion_model()
        base_ring = cm.common_parent(m.base_ring(), self._reference.base_ring())

        from sage.all import matrix

        self._matrix = matrix(base_ring, m, immutable=True)

        super().__init__(
            base_ring, reference, category=category or self._reference.category()
        )

    @cached_method
    def _sgn(self):
        r"""
        Return the sign of the determinant of the matrix underlying this
        surface, i.e., whether the matrix reversed orientation or not.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S
            sage: S._sgn()
            -1

        """
        return self._matrix.det().sign()

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S

            sage: S.polygon(0)
            Polygon(vertices=[(0, 0), (a, -a), (a + 2, -a), (2*a + 2, 0), (2*a + 2, 2), (a + 2, a + 2), (a, a + 2), (0, 2)])

        """
        return self._matrix * self._reference.polygon(label)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S

            sage: S.opposite_edge(0, 0)
            (2, 0)

        """
        reference_edge = edge
        if self._sgn() == -1:
            reference_edge = len(self.polygon(label).edges()) - 1 - edge

        opposite_label, opposite_edge = self._reference.opposite_edge(
            label, reference_edge
        )

        if self._sgn() == -1:
            opposite_edge = (
                len(self._reference.polygon(opposite_label).edges()) - 1 - opposite_edge
            )

        return opposite_label, opposite_edge

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])
            sage: S = r * S

            sage: S.is_mutable()
            False

        """
        return False

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: matrix([[0, 1], [1, 0]]) * S
            Translation Surface in H_3(4) built from 2 squares and a regular octagon
            sage: matrix([[0, 2], [1, 0]]) * S
            Translation Surface in H_3(4) built from a rhombus, a rectangle and an octagon

        ::

            sage: m = matrix([[2, 1], [1, 1]])
            sage: m * translation_surfaces.infinite_staircase()
            GL2R image of The infinite staircase

        """
        if self.is_finite_type():
            from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

            S = MutableOrientedSimilaritySurface.from_surface(self)
            S.set_immutable()
            return repr(S)

        return f"GL2R image of {self._reference!r}"

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: r = matrix(ZZ,[[0, 1], [1, 0]])

            sage: hash(r * S) == hash(r * S)
            True

        """
        return hash((self._reference, self._matrix))

    def __eq__(self, other):
        r"""
        Return whether this image is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: m = matrix(ZZ,[[0, 1], [1, 0]])
            sage: m * S == m * S
            True

        """
        if not isinstance(other, GL2RImageSurface):
            return False

        return (
            self._reference == other._reference
            and self._matrix == other._matrix
            and self.base_ring() == other.base_ring()
        )

    def is_triangulated(self, limit=None):
        r"""
        Return whether this surface is triangulated.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.octagon_and_squares()
            sage: m = matrix(ZZ,[[0, 1], [1, 0]])
            sage: (m * S).is_triangulated()
            False

        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        return self._reference.is_triangulated(limit=limit)


class LazyMutableOrientedSimilaritySurface(
    LazyOrientedSimilaritySurface, MutableOrientedSimilaritySurface_base
):
    r"""
    A helper surface for :class:`LazyDelaunayTriangulatedSurface`.

    A mutable wrapper of an (infinite) reference surface. When a polygon is not
    present in this wrapper yet, it is taken from the reference surface and can
    then be modified.

    .. NOTE::

        This surface does not implement the entire surface interface correctly.
        It just supports the operations in the way that they are necessary to
        make :class:`LazyDelaunayTriangulatedSurface` work.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()

        sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
        sage: T = LazyMutableOrientedSimilaritySurface(S)
        sage: p = T.polygon(0)
        sage: p
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
        sage: q = p * 2

        sage: S.replace_polygon(0, q)
        Traceback (most recent call last):
        ...
        AttributeError: '_InfiniteStaircase_with_category' object has no attribute 'replace_polygon'...

        sage: T.replace_polygon(0, q)
        sage: T.polygon(0)
        Polygon(vertices=[(0, 0), (2, 0), (2, 2), (0, 2)])

    """

    def __init__(self, surface, category=None):
        from flatsurf.geometry.categories import SimilaritySurfaces

        if surface not in SimilaritySurfaces().Oriented().WithoutBoundary():
            raise NotImplementedError("cannot handle surfaces with boundary yet")

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        self._surface = MutableOrientedSimilaritySurface(surface.base_ring())

        super().__init__(
            surface.base_ring(), surface, category=category or surface.category()
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``True``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.is_mutable()
            True

        """
        return True

    def replace_polygon(self, label, polygon):
        r"""
        Swap out the polygon with the label ``label`` with ``polygon``.

        The polygons must have the same number of sides since gluings are kept.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.replace_polygon(0, T.polygon(0))

        """
        self._ensure_polygon(label)
        return self._surface.replace_polygon(label, polygon)

    def glue(self, x, y):
        r"""
        Glue the (label, edge) pair ``x`` with the pair ``y`` in this surface.

        This unglues any existing gluings of these edges.

        .. NOTE::

            After a sequence of such glue operations, no edges must be unglued.
            Otherwise, gluings get copied over from the underlying surface with
            confusing side effects.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.gluings()
            (((0, 0), (1, 2)), ((0, 1), (-1, 3)), ((0, 2), (1, 0)), ((0, 3), (-1, 1)), ((1, 0), (0, 2)), ((1, 1), (2, 3)), ((1, 2), (0, 0)), ((1, 3), (2, 1)), ((-1, 0), (-2, 2)),
             ((-1, 1), (0, 3)), ((-1, 2), (-2, 0)), ((-1, 3), (0, 1)), ((2, 0), (3, 2)), ((2, 1), (1, 3)), ((2, 2), (3, 0)), ((2, 3), (1, 1)), …)
            sage: T.glue((0, 0), (1, 0))
            sage: T.glue((1, 2), (0, 2))
            sage: T.gluings()
            (((0, 0), (1, 0)), ((0, 1), (-1, 3)), ((0, 2), (1, 2)), ((0, 3), (-1, 1)), ((1, 0), (0, 0)), ((1, 1), (2, 3)), ((1, 2), (0, 2)), ((1, 3), (2, 1)), ((-1, 0), (-2, 2)),
             ((-1, 1), (0, 3)), ((-1, 2), (-2, 0)), ((-1, 3), (0, 1)), ((2, 0), (3, 2)), ((2, 1), (1, 3)), ((2, 2), (3, 0)), ((2, 3), (1, 1)), …)

        """
        return self._surface.glue(x, y)

    def _ensure_gluings(self, label):
        r"""
        Make sure that the surface used to internally represent this surface
        has copied over all the gluings for ``label`` from the underlying
        surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T._ensure_polygon(0)
            sage: T._ensure_gluings(0)

        """
        self._ensure_polygon(label)
        for edge in range(len(self._surface.polygon(label).vertices())):
            cross = self._surface.opposite_edge(label, edge)
            if cross is None:
                cross_label, cross_edge = self._reference.opposite_edge(label, edge)
                self._ensure_polygon(cross_label)

                assert (
                    self._surface.opposite_edge(cross_label, cross_edge) is None
                ), "surface must not have a boundary"

                # Note that we cannot detect whether something has been
                # explicitly unglued. So we just reestablish any gluings of
                # this edge.
                self._surface.glue((label, edge), (cross_label, cross_edge))

    def _ensure_polygon(self, label):
        r"""
        Make sure that the surface used to internally represent this surface
        has copied over the polygon ``label`` from the underlying surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T._ensure_polygon(0)

        """
        if label not in self._surface.labels():
            self._surface.add_polygon(self._reference.polygon(label), label=label)

    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        self._ensure_polygon(label)
        return self._surface.polygon(label)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()

            sage: from flatsurf.geometry.lazy import LazyMutableOrientedSimilaritySurface
            sage: T = LazyMutableOrientedSimilaritySurface(S)
            sage: T.opposite_edge(0, 0)
            (1, 2)

        """
        self._ensure_polygon(label)
        self._ensure_gluings(label)
        cross_label, cross_edge = self._surface.opposite_edge(label, edge)
        self._ensure_polygon(cross_label)
        self._ensure_gluings(cross_label)
        return cross_label, cross_edge


class LazyDelaunayTriangulatedSurface(OrientedSimilaritySurface):
    r"""
    Delaunay triangulation of an (infinite type) surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
        sage: len(S.polygon(S.root()).vertices())
        3
        sage: TestSuite(S).run()  # long time (.8s)
        sage: S.is_delaunay_triangulated()
        True

        sage: from flatsurf.geometry.lazy import LazyDelaunayTriangulatedSurface
        sage: isinstance(S, LazyDelaunayTriangulatedSurface)
        True

    ::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: S = chamanara_surface(QQ(1/2))
        sage: m = matrix([[2,1],[1,1]])**4
        sage: S = (m*S).delaunay_triangulation()
        sage: TestSuite(S).run()  # long time (1s)
        sage: S.is_delaunay_triangulated()
        True
        sage: TestSuite(S).run()  # long time (.5s)

        sage: from flatsurf.geometry.lazy import LazyDelaunayTriangulatedSurface
        sage: isinstance(S, LazyDelaunayTriangulatedSurface)
        True

    """

    def __init__(self, similarity_surface, direction=None, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyDelaunayTriangulatedSurface; use relabel() to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyDelaunayTriangulatedSurface()"
                )

        if direction is not None:
            direction = None
            import warnings

            warnings.warn(
                "the direction keyword argument has been deprecated for LazyDelaunayTriangulatedSurface and will be removed in a future version of sage-flatsurf; "
                "its value is ignored in this version of sage-flatsurf; if you see this message when restoring a pickle, the object might not be fully functional"
            )

        if similarity_surface.is_mutable():
            raise ValueError("surface must be immutable")

        if not similarity_surface.is_connected():
            raise NotImplementedError("surface must be connected")

        if not similarity_surface.is_triangulated():
            raise ValueError("surface must be triangulated")

        self._reference = similarity_surface

        # This surface will converge to the Delaunay Triangulation
        self._surface = LazyMutableOrientedSimilaritySurface(similarity_surface)

        # Set of labels corresponding to known delaunay polygons
        self._certified_labels = set()

        # Triangle flips (as morphisms) that have been performed so far.
        self._flips = []

        # Triangulate the base polygon
        root = self._surface.root()

        # Certify the base polygon (or apply flips...)
        while not self._certify_or_improve(root):
            pass

        OrientedSimilaritySurface.__init__(
            self,
            self._surface.base_ring(),
            category=category or self._surface.category(),
        )

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
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
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.roots()
            ((0, 0),)

        """
        return self._surface.roots()

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.polygon((0, 0))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        if label not in self.labels():
            raise ValueError("no polygon with this label")

        if label not in self._certified_labels:
            # If the label is not final in this surface, we walk the surface
            # and thereby certify its polygons until we find that label.
            # Note that this is somewhat inefficient since we start the walk
            # from the start every time. However, the certification process is
            # what consumes time, so unless there are a lot of labels, this
            # should not impact performance too much.
            for certified_label in self._walk():
                if label == certified_label:
                    assert label in self._certified_labels
                    break

        return self._surface.polygon(label)

    def _walk(self):
        r"""
        Return an iterator that walks the labels of the surface in order.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()

            sage: from itertools import islice
            sage: list(islice(S._walk(), 3))
            [(0, 0), (1, 1), (-1, 1)]

        """
        visited = set()
        from collections import deque

        next = deque(
            [(self.root(), 0), (self.root(), 1), (self.root(), 2)],
        )
        while next:
            label, edge = next.popleft()
            if label in visited:
                continue

            yield label

            visited.add(label)

            for edge in range(3):
                next.append(self.opposite_edge(label, edge))

    @cached_method
    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.opposite_edge((0, 0), 0)
            ((1, 1), 1)

        """
        self.polygon(label)
        while True:
            cross_label, cross_edge = self._surface.opposite_edge(label, edge)
            if self._certify_or_improve(cross_label):
                break

        return self._surface.opposite_edge(label, edge)

    def _certify_or_improve(self, label):
        r"""
        This method attempts to develop the circumscribing disk about the polygon
        with label ``label`` into the surface.

        The method returns True if this is successful. In this case the label
        is added to the set _certified_labels. It returns False if it failed to
        develop the disk into the surface. (In this case the original polygon was
        not a Delaunay triangle.

        The algorithm divides any non-certified polygon in self._s it encounters
        into triangles. If it encounters a pair of triangles which need a diagonal
        flip then it does the flip.
        """
        if label in self._certified_labels:
            # Already certified.
            return True
        p = self._surface.polygon(label)
        assert len(p.vertices()) == 3

        c = p.circumscribing_circle()

        # Develop through each of the 3 edges:
        for e in range(3):
            edge_certified = False
            # This keeps track of a chain of polygons the disk develops through:
            edge_stack = []

            # We repeat this until we can verify that the portion of the circle
            # that passes through the edge e developes into the surface.
            while not edge_certified:
                if len(edge_stack) == 0:
                    # Start at the beginning with label l and edge e.
                    # The 3rd coordinate in the tuple represents what edge to develop
                    # through in the triangle opposite this edge.
                    edge_stack = [(label, e, 1, c)]
                ll, ee, step, cc = edge_stack[len(edge_stack) - 1]

                lll, eee = self._surface.opposite_edge(ll, ee)

                if lll not in self._certified_labels:
                    ppp = self._surface.polygon(lll)
                    assert len(ppp.vertices()) == 3

                    if self._surface._delaunay_edge_needs_flip(ll, ee):
                        # Perform the flip
                        self._surface.triangle_flip(ll, ee, in_place=True)

                        # If we touch the original polygon, then we return False.
                        if label == ll or label == lll:
                            return False
                        # We might have flipped a polygon from earlier in the chain
                        # In this case we need to trim the stack down so that we recheck
                        # that polygon.
                        for index, tup in enumerate(edge_stack):
                            if tup[0] == ll or tup[0] == lll:
                                edge_stack = edge_stack[:index]
                                break
                        # The following if statement makes sure that we check both subsequent edges of the
                        # polygon opposite the last edge listed in the stack.
                        if len(edge_stack) > 0:
                            ll, ee, step, cc = edge_stack.pop()
                            edge_stack.append((ll, ee, 1, cc))
                        continue

                    # If we reach here then we know that no flip was needed.
                    ccc = self._surface.edge_transformation(ll, ee) * cc

                    # Check if the disk passes through the next edge in the chain.
                    lp = ccc.line_segment_position(
                        ppp.vertex((eee + step) % 3), ppp.vertex((eee + step + 1) % 3)
                    )
                    if lp == 1:
                        # disk passes through edge and opposite polygon is not certified.
                        edge_stack.append((lll, (eee + step) % 3, 1, ccc))
                        continue

                    # We reach this point if the disk doesn't pass through the edge eee+step of polygon lll.

                # Either lll is already certified or the disk didn't pass
                # through edge (lll,eee+step)

                # Trim off unnecessary edges off the stack.
                # prune_count=1
                ll, ee, step, cc = edge_stack.pop()
                if step == 1:
                    # if we have just done step 1 (one edge), move on to checking
                    # the next edge.
                    edge_stack.append((ll, ee, 2, cc))
                # if we have pruned an edge, continue to look at pruning in the same way.
                while step == 2 and len(edge_stack) > 0:
                    ll, ee, step, cc = edge_stack.pop()
                    # prune_count= prune_count+1
                    if step == 1:
                        edge_stack.append((ll, ee, 2, cc))
                if len(edge_stack) == 0:
                    # We're done with this edge
                    edge_certified = True
        self._certified_labels.add(label)
        return True

    def is_triangulated(self, limit=None):
        r"""
        Return whether this surface is triangulated, which it naturally is.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.is_triangulated()
            True

        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_triangulated() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        return True

    def is_delaunay_triangulated(self, limit=None):
        r"""
        Return whether this surface is Delaunay triangulated, which it
        naturally is.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.is_delaunay_triangulated()
            True

        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_delaunay_triangulated() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        return True

    def labels(self):
        r"""
        Return the labels of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.labels`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()
            sage: S.labels()
            ((0, 0), (1, 1), (-1, 1), (0, 1), (1, 0), (2, 0), (-1, 0), (-2, 0), (2, 1), (3, 1), (-2, 1), (-3, 1), (3, 0), (4, 0), (-3, 0), (-4, 0), …)

        """
        return self._surface.labels()

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.delaunay_triangulation()) == hash(S.delaunay_triangulation())
            True

        """
        return hash((self._reference))

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.delaunay_triangulation() == S.delaunay_triangulation()
            True

        """
        if not isinstance(other, LazyDelaunayTriangulatedSurface):
            return False

        return self._reference == other._reference

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_triangulation()

        """
        reference = self._reference
        if isinstance(reference, LazyTriangulatedSurface):
            reference = reference._reference
        return f"Delaunay triangulation of {reference!r}"


class LazyDelaunaySurface(OrientedSimilaritySurface):
    r"""
    Delaunay cell decomposition of a (possibly infinite type) surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.infinite_staircase()
        sage: m = matrix([[2, 1], [1, 1]])
        sage: S = (m * S).delaunay_decomposition()

        sage: S.polygon(S.root())
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: S.is_delaunay_decomposed()
        True

        sage: TestSuite(S).run()  # long time (2s)

        sage: from flatsurf.geometry.lazy import LazyDelaunaySurface
        sage: isinstance(S, LazyDelaunaySurface)
        True

    ::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: S = chamanara_surface(QQ(1/2))
        sage: m = matrix([[3, 4], [-4, 3]]) * matrix([[4, 0],[0, 1/4]])
        sage: S = (m * S).delaunay_decomposition()
        sage: S.is_delaunay_decomposed()
        True

        sage: TestSuite(S).run()  # long time (1.5s)

        sage: from flatsurf.geometry.lazy import LazyDelaunaySurface
        sage: isinstance(S, LazyDelaunaySurface)
        True

    """

    def __init__(self, similarity_surface, direction=None, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyDelaunaySurface; use relabel() to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyDelaunaySurface()"
                )

        if direction is not None:
            direction = None
            import warnings

            warnings.warn(
                "the direction keyword argument has been deprecated for LazyDelaunayTriangulatedSurface and will be removed in a future version of sage-flatsurf; "
                "its value is ignored in this version of sage-flatsurf; if you see this message when restoring a pickle, the object might not be fully functional"
            )

        if similarity_surface.is_mutable():
            raise ValueError("surface must be immutable")

        if not similarity_surface.is_delaunay_triangulated():
            raise ValueError("surface must be triangulated")

        self._reference = similarity_surface

        super().__init__(
            similarity_surface.base_ring(),
            category=category or similarity_surface.category(),
        )

    @cached_method
    def polygon(self, label):
        r"""
        Return the polygon with ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.polygon`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.polygon((0, 0))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        if label not in self._reference.labels():
            raise ValueError("no polygon with this label")

        cell, edges = self._cell(label)

        if label != self._label(cell):
            raise ValueError("no polygon with this label")

        edges = [self._reference.polygon(edge[0]).edge(edge[1]) for edge in edges]

        from flatsurf import Polygon

        return Polygon(edges=edges)

    @cached_method
    def _label(self, cell):
        r"""
        Return a canonical label for the Delaunay cell that is made up by the
        Delaunay triangles ``cell``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S._label(frozenset({
            ....:     (0, 0),
            ....:     (0, 1)}))
            (0, 0)

        """
        for label in self._reference.labels():
            if label in cell:
                return label

    @cached_method
    def _normalize_label(self, label):
        r"""
        Return a canonical label for the Delaunay cell that contains the
        Delaunay triangle ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S._normalize_label((0, 0))
            (0, 0)
            sage: S._normalize_label((0, 1))
            (0, 0)

        """
        cell, _ = self._cell(label)
        return self._label(cell)

    @cached_method
    def _cell(self, label):
        r"""
        Return the labels of the Delaunay triangles that contain the Delaunay
        triangle ``label`` together with the interior edges in that cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()

        This cell (a square) is formed by two triangles that form a cylinder,
        i.e., the two triangles are glued at two of their edges::

            sage: S._cell((0, 0))
            (frozenset({(0, 0), (0, 1)}),
             [((0, 0), 0), ((0, 0), 1), ((0, 1), 1), ((0, 1), 2)])

        """
        edges = []
        cell = set()
        explore = [(label, 2), (label, 1), (label, 0)]

        while explore:
            triangle, edge = explore.pop()

            cell.add(triangle)

            delaunay = self._reference._delaunay_edge_needs_join(triangle, edge)

            if not delaunay:
                edges.append((triangle, edge))
                continue

            cross_triangle, cross_edge = self._reference.opposite_edge(triangle, edge)

            for shift in [2, 1]:
                next_triangle, next_edge = cross_triangle, (cross_edge + shift) % 3

                if (next_triangle, next_edge) in edges:
                    raise NotImplementedError
                if (next_triangle, next_edge) in explore:
                    raise NotImplementedError

                explore.append((next_triangle, next_edge))

        cell = frozenset(cell)
        normalized_label = self._label(cell)
        if normalized_label != label:
            return self._cell(normalized_label)

        return cell, edges

    @cached_method
    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and edge index when crossing over the ``edge``
        of the polygon ``label``.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.opposite_edge`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.opposite_edge((0, 0), 0)
            ((1, 1), 2)

        """
        if label not in self._reference.labels():
            raise ValueError

        cell, edges = self._cell(label)

        if label != self._label(cell):
            raise ValueError

        edge = edges[edge]

        cross_triangle, cross_edge = self._reference.opposite_edge(*edge)

        cross_cell, cross_edges = self._cell(cross_triangle)
        cross_label = self._label(cross_cell)
        return cross_label, cross_edges.index((cross_triangle, cross_edge))

    def roots(self):
        r"""
        Return root labels for the polygons forming the connected
        components of this surface.

        This implements
        :meth:`flatsurf.geometry.categories.polygonal_surfaces.PolygonalSurfaces.ParentMethods.roots`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.roots()
            ((0, 0),)

        """
        return self._reference.roots()

    def is_compact(self):
        r"""
        Return whether this surface is compact as a topological space.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_compact`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.is_compact()
            False

        """
        return self._reference.is_compact()

    def is_mutable(self):
        r"""
        Return whether this surface is mutable, i.e., return ``False``.

        This implements
        :meth:`flatsurf.geometry.categories.topological_surfaces.TopologicalSurfaces.ParentMethods.is_mutable`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase().delaunay_decomposition()
            sage: S.is_mutable()
            False

        """
        return False

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        :meth:`__eq__`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: hash(S.delaunay_decomposition()) == hash(S.delaunay_decomposition())
            True

        """
        return hash(self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        See :meth:`SimilaritySurfaces.FiniteType._test_eq_surface` for details
        on this notion of equality.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.lazy import LazyDelaunaySurface
            sage: S = translation_surfaces.infinite_staircase()
            sage: m = matrix([[2, 1], [1, 1]])
            sage: S = (m * S).delaunay_triangulation()
            sage: S == S
            True

        """
        if not isinstance(other, LazyDelaunaySurface):
            return False

        return self._reference == other._reference

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.delaunay_decomposition()
            Delaunay cell decomposition of The infinite staircase

        """
        reference = self._reference
        if isinstance(reference, LazyDelaunayTriangulatedSurface):
            reference = reference._reference
        if isinstance(reference, LazyTriangulatedSurface):
            reference = reference._reference
        return f"Delaunay cell decomposition of {reference!r}"

    def is_delaunay_decomposed(self, limit=None):
        r"""
        Return whether this surface is decomposed into Delaunay cells, which it
        naturally is.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.infinite_staircase()
            sage: S.delaunay_decomposition().is_delaunay_decomposed()
            True

        """
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for is_delaunay_decomposed() and will be removed from a future version of sage-flatsurf; "
                "if you rely on this check, you can try to run this method on MutableOrientedSimilaritySurface.from_surface(surface, labels=surface.labels()[:limit])"
            )

        return True


class LazyRelabeledSurface(LazyOrientedSimilaritySurface):
    r"""
    A relabeled surface which forwards all requests to an underlying reference
    surface after translation of labels.

    Subclasses may override ``_to_reference_label`` and
    ``_from_reference_label`` to establish a custom mapping of labels.
    Otherwise, labels are mapped to the non-negative integers in order.

    EXAMPLES::

        sage: from flatsurf.geometry.chamanara import chamanara_surface
        sage: S = chamanara_surface(1/2)

    TESTS::

        sage: from flatsurf.geometry.lazy import LazyRelabeledSurface
        sage: isinstance(S, LazyRelabeledSurface)
        True

    """

    def __init__(self, reference, category=None):
        super().__init__(
            base_ring=reference.base_ring(),
            reference=reference,
            category=category or reference.category(),
        )

    def _to_reference_label(self, label):
        r"""
        Return the image of ``label`` in the underlying reference surface.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S._to_reference_label(0)
            (0, 1, 0)

        """
        return self._reference.labels()[label]

    def _from_reference_label(self, reference_label):
        r"""
        Return the preimage of ``reference_label`` in the labels of this surface.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S._from_reference_label((0, 1, 0))
            0

        """
        label = 0
        for lbl in self._reference.labels():
            if lbl == reference_label:
                return label
            label += 1

        raise ValueError

    def _test_label_map(self, **options):
        r"""
        Verify that :meth:`_to_reference_label` and
        :meth:`_from_reference_label` are inverse to each other.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S._test_label_map()

        """
        tester = self._tester(**options)

        from itertools import islice

        for label in islice(self._reference.labels(), 30):
            tester.assertEqual(
                label, self._to_reference_label(self._from_reference_label(label))
            )

        for label in islice(
            range(30)
            if not self.is_finite_type()
            else range(len(self._reference.labels())),
            30,
        ):
            tester.assertEqual(
                label, self._from_reference_label(self._to_reference_label(label))
            )

    def labels(self):
        r"""
        Return the labels of this surface after renaming.

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.labels()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, …)

        """
        from flatsurf.geometry.surface import LabelsFromView

        if self._reference.is_finite_type():
            return LabelsFromView(self, range(len(self._reference.labels())))

        from sage.all import NN

        return LabelsFromView(self, NN)

    def polygon(self, label):
        r"""
        Return the polygon with ``label`` in this surface.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.polygon(0)
            Polygon(vertices=[(0, 0), (1, 0), (-1, 2), (-1, 1)])

        """
        return self._reference.polygon(self._to_reference_label(label))

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon label and its edge that is across from the polygon
        with ``label`` and its ``edge``.

        Return ``None`` if there is no polygon glued to that edge.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.opposite_edge(0, 0)
            (1, 0)
            sage: S.opposite_edge(1, 0)
            (0, 0)

        """
        label = self._to_reference_label(label)
        opposite_edge = self._reference.opposite_edge(label, edge)

        if opposite_edge is None:
            return None

        opposite_label, opposite_edge = opposite_edge
        opposite_label = self._from_reference_label(opposite_label)

        return opposite_label, opposite_edge

    def roots(self):
        r"""
        Return the labels of the roots of the components of this surface.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.roots()
            (0,)

        """
        return tuple(
            self._from_reference_label(label) for label in self._reference.roots()
        )

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with
        ``__eq__``.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: hash(S) == hash(S)
            True

        """
        return hash(self._reference)

    def __eq__(self, other):
        r"""
        Return whether this surface and ``other`` are indistinguishable.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: T = chamanara_surface(1/2)
            sage: S == T
            True

        """
        if type(self) != type(other):
            # Since we encourage subclassing this surface, we are very strict here.
            return False

        return self._reference == other._reference

    def _repr_(self):
        r"""
        Return a printable representation of this surface.

        Since the relabeling is often just done to make the labels a bit easier
        to work with, we do not mention it when printing this surface.

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S
            Minimal Translation Cover of Chamanara surface with parameter 1/2

        """
        return repr(self._reference)
