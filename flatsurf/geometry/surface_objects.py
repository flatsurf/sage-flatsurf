r"""
Geometric objects on surfaces.

This includes singularities, saddle connections and cylinders.

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks
"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2017-2020 W. Patrick Hooper
#                      2017-2020 Vincent Delecroix
#                           2023 Julian Rüth
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
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.plot.graphics import Graphics
from sage.plot.polygon import polygon2d
from sage.rings.qqbar import AA
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject
from sage.structure.element import Element

from flatsurf.geometry.similarity import SimilarityGroup


def Singularity(similarity_surface, label, v, limit=None):
    r"""
    Return the point of ``similarity_surface`` at the ``v``-th vertex of the
    polygon ``label``.

    If the surface is infinite, the ``limit`` can be set. In this case the
    construction of the singularity is successful if the sequence of vertices
    hit by passing through edges closes up in ``limit`` or less steps.

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import TranslationSurfaceGenerators
        sage: s=TranslationSurfaceGenerators.veech_2n_gon(5)
        sage: from flatsurf.geometry.surface_objects import Singularity
        sage: sing=Singularity(s, 0, 1)
        doctest:warning
        ...
        UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
        sage: print(sing)
        Vertex 1 of polygon 0
        sage: TestSuite(sing).run()

    """
    import warnings

    warnings.warn(
        "Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead."
    )
    return similarity_surface.point(
        label, similarity_surface.polygon(label).vertex(v), limit=limit
    )


class SurfacePoint(Element):
    r"""
    A point on ``surface``.

    INPUT:

    - ``surface`` -- a similarity surface

    - ``label`` -- a polygon label for the polygon with respect to which the
      ``point`` coordinates can be made sense of

    - ``point`` -- coordinates of a point in the polygon ``label`` or the index
      of the vertex of the polygon with ``label``

    - ``ring`` -- a SageMath ring or ``None`` (default: ``None``); the
      coordinate ring for ``point``

    - ``limit`` -- an integer or ``None`` (default: ``None`` for an unlimited
      number of steps); if this is a singularity of the surface, then this
      limits the number of edges that are crossed to determine all the edges
      adjacent to that singularity. An error is raised if the limit is
      insufficient.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces

        sage: permutation = SymmetricGroup(2)('(1, 2)')
        sage: S = translation_surfaces.origami(permutation, permutation)

    A point can have a single representation with coordinates when it is in
    interior of a polygon::

        sage: S.point(1, (1/2, 1/2))
        Point (1/2, 1/2) of polygon 1

    A point can have two representations when it is in interior of an edge::

        sage: p = S.point(1, (1/2, 0))
        sage: q = S.point(2, (1/2, 1))
        sage: p == q
        True
        sage: p.coordinates(2)
        ((1/2, 1),)

    A point can have even more representations when it is a vertex::

        sage: S.point(1, (0, 0))
        Vertex 0 of polygon 1

    TESTS:

    Verify that #275 has been resolved, i.e., points on the boundary can be
    created::

        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
        sage: S = MutableOrientedSimilaritySurface(QQ)
        sage: S.add_polygon(Polygon(vertices=[(0,0), (-1, -1), (1,0)]))
        0
        sage: S.add_polygon(Polygon(vertices=[(0,0), (0, 1), (-1,-1)]))
        1
        sage: S.glue((0, 0), (1, 2))
        sage: S.set_immutable()
        sage: S
        Translation Surface with boundary built from 2 triangles

        sage: S(0, (0, 0))
        Vertex 0 of polygon 0
        sage: S(0, (1/2, 0))
        Point (1/2, 0) of polygon 0
        sage: S(0, (1, 0))
        Vertex 2 of polygon 0

    """

    def __init__(self, surface, label, point, ring=None, limit=None):
        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument when creating points and will be removed without replacement in a future version of sage-flatsurf"
            )

        self._surface = surface

        if ring is None:
            ring = surface.base_ring()

        if ring is not surface.base_ring():
            import warnings

            warnings.warn(
                "the ring parameter is deprecated and will be removed in a future version of sage-flatsurf; define the surface over a larger ring instead so that this points' coordinates live in the base ring"
            )

        polygon = surface.polygon(label)

        from sage.all import ZZ

        if point in ZZ:
            point = surface.polygon(label).vertex(point)

        point = (ring**2)(point)
        point.set_immutable()

        position = polygon.get_point_position(point)

        if not position.is_inside():
            raise NotImplementedError(
                "point must be positioned within the polygon with the given label"
            )

        if position.is_in_interior():
            self._representatives = {(label, point)}
        elif position.is_in_edge_interior():
            self._representatives = {(label, point)}

            opposite_edge = surface.opposite_edge(label, position.get_edge())
            if opposite_edge is not None:
                cross_label, cross_edge = opposite_edge
                cross_point = surface.edge_transformation(label, position.get_edge())(
                    point
                )
                cross_point.set_immutable()

                self._representatives.add((cross_label, cross_point))
        elif position.is_vertex():
            self._representatives = set()

            source_edge = position.get_vertex()

            def collect_representatives(label, source_edge, direction, limit):
                def rotate(label, source_edge, direction):
                    if direction == -1:
                        source_edge = (source_edge - 1) % len(
                            surface.polygon(label).vertices()
                        )
                    opposite_edge = surface.opposite_edge(label, source_edge)
                    if opposite_edge is None:
                        return None
                    label, source_edge = opposite_edge
                    if direction == 1:
                        source_edge = (source_edge + 1) % len(
                            surface.polygon(label).vertices()
                        )
                    return label, source_edge

                while True:
                    self._representatives.add((label, source_edge))

                    # Rotate to the next edge that is leaving at the vertex
                    rotated = rotate(label, source_edge, direction)
                    if rotated is None:
                        # Surface is disconnected
                        return

                    label, source_edge = rotated

                    if limit is not None:
                        limit -= 1
                        if limit < 0:
                            raise ValueError(
                                "number of edges at singularity exceeds limit"
                            )

                    if (label, source_edge) in self._representatives:
                        break

            # Collect respresentatives of edge walking clockwise and counterclockwise
            collect_representatives(label, source_edge, 1, limit)
            collect_representatives(label, source_edge, -1, limit)

            self._representatives = {
                (label, surface.polygon(label).vertex(vertex))
                for (label, vertex) in self._representatives
            }
        else:
            raise NotImplementedError

        self._representatives = frozenset(self._representatives)

        super().__init__(surface)

    def surface(self):
        r"""
        Return the surface containing this point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (1/2, 1/2))
            sage: p.surface() is S
            True

        """
        return self._surface

    def is_vertex(self):
        r"""
        Return whether this point is a singularity of the surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))

            sage: p.is_vertex()
            True

        """
        label, coordinates = self.representative()
        position = self.surface().polygon(label).get_point_position(coordinates)
        return position.is_vertex()

    def one_vertex(self):
        r"""
        Return a pair (l, v) from the equivalence class of this singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.one_vertex()  # random output: depends on the Python version
            doctest:warning
            ...
            UserWarning: one_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use (label, coordinates) = point.representative(); vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
            (2, 1)

        """
        import warnings

        warnings.warn(
            "one_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use (label, coordinates) = point.representative(); vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead"
        )
        label, coordinates = self.representative()
        vertex = (
            self.surface().polygon(label).get_point_position(coordinates).get_vertex()
        )
        return label, vertex

    def representatives(self):
        r"""
        Return the representatives of this point as pairs of polygon labels and
        coordinates.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.representatives()
            frozenset({(1, (0, 0)), (1, (1, 1)), (2, (0, 1)), (2, (1, 0))})

        """
        return self._representatives

    def representative(self):
        r"""
        Return a representative of this point, i.e., the first of
        :meth:`representatives`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.representative()  # random output: depends on the Python version
            (2, (1, 0))

        """
        return next(iter(self.representatives()))

    def vertex_set(self):
        r"""
        Return the list of pairs (l, v) in the equivalence class of this singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: list(p.vertex_set())  # random output: ordering depends on the Python version
            doctest:warning
            ...
            UserWarning: vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
            [(2, 1), (1, 2), (1, 0), (2, 3)]

        """
        import warnings

        warnings.warn(
            "vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead"
        )

        return [
            (
                label,
                self.surface()
                .polygon(label)
                .get_point_position(coordinates)
                .get_vertex(),
            )
            for label, coordinates in self.representatives()
        ]

    def contains_vertex(self, label, v=None):
        r"""
        Checks if the pair ``(label, v)`` is in the equivalence class returning
        true or false. If ``v`` is ``None``, the both the pair ``(label, v)``
        is passed as a single parameter in ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.contains_vertex((1, 0))
            doctest:warning
            ...
            UserWarning: contains_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use the == operator instead
            doctest:warning
            ...
            UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
            True
            sage: p.contains_vertex(label=1, v=0)
            True

        """
        import warnings

        warnings.warn(
            "contains_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use the == operator instead"
        )

        if v is None:
            label, v = label

        return Singularity(self.surface(), label, v) == self

    def num_coordinates(self):
        r"""
        Return the number of different coordinate representations of the point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.num_coordinates()
            doctest:warning
            ...
            UserWarning: num_coordinates() is deprecated and will be removed in a future version of sage-flatsurf; use len(representatives()) instead.
            4

        """
        import warnings

        warnings.warn(
            "num_coordinates() is deprecated and will be removed in a future version of sage-flatsurf; use len(representatives()) instead."
        )

        return len(self._representatives)

    def labels(self):
        r"""
        Return the labels of polygons containing the point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

        For a point in the interior of polygon, there is exactly one label::

            sage: p = S.point(0, (1/2, 1/2))
            sage: p.labels()
            {0}

        For a point in the interior of an edge of a polygon, there can be up to
        two labels::

            sage: p = S.point(0, (0, 1/2))
            sage: p.labels()
            {0, 2}

        For a point at a vertex, there can be more labels::

            sage: p = S.point(0, (0, 0))
            sage: p.labels()
            {0, 1, 2}

        """
        return {label for (label, _) in self._representatives}

    def coordinates(self, label):
        r"""
        Return coordinates for the point in the in the polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: p = S.point(0, (0, 0))
            sage: p.coordinates(0)  # random output: order depends on the Python version
            ((0, 0), (1, 0), (0, 1), (1, 1))

        """
        return tuple(
            coordinates for (l, coordinates) in self._representatives if l == label
        )

    def graphical_surface_point(self, graphical_surface=None):
        r"""
        Return a
        :class:`flatsurf.graphical.surface_point.GraphicalSurfacePoint` to
        represent this point graphically.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: G = p.graphical_surface_point()

        """
        from flatsurf.graphical.surface_point import GraphicalSurfacePoint

        return GraphicalSurfacePoint(self, graphical_surface=graphical_surface)

    def plot(self, *args, **kwargs):
        r"""
        Return a plot of this point.

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (0, 0))
            sage: p.plot()
            ...Graphics object consisting of 1 graphics primitive

        .. jupyter-execute::

            sage: p = S.point(0, (0, 25/12))
            sage: p.plot()
            ...Graphics object consisting of 1 graphics primitive

        """
        graphical_surface = None
        if args:
            graphical_surface = args[0]
            args = args[1:]

        return self.graphical_surface_point(graphical_surface=graphical_surface).plot(
            *args, **kwargs
        )

    def __repr__(self):
        r"""
        Return a printable representation of this point.

        TESTS::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: S.point(0, (1/2, 1/2))
            Point (1/2, 1/2) of polygon 0

        """

        def render(label, coordinates):
            if self.is_vertex():
                vertex = (
                    self.surface()
                    .polygon(label)
                    .get_point_position(coordinates)
                    .get_vertex()
                )
                return "Vertex {} of polygon {}".format(vertex, label)

            return "Point {} of polygon {}".format(coordinates, label)

        # We pick a specific representative to make our lives easier when doctesting
        return min(
            render(label, coordinates)
            for (label, coordinates) in self.representatives()
        )

    def __eq__(self, other):
        r"""
        Return whether this point is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p == p
            True

            sage: q = S.point(0, (1/2, 1/3))
            sage: p == q
            False

        TESTS:

        Verify that points can be compared to non-points so they can be put into sets and dicts with other objects::

            sage: p == 42
            False

        """
        if self is other:
            return True
        if not isinstance(other, SurfacePoint):
            return False
        if not self._surface == other._surface:
            return False
        return self._representatives == other._representatives

    def _test_category(self, **options):
        r"""
        Check that this point inherits from the element class of its surface's
        category.

        Overridden to disable these tests when this is a point of a mutable
        surface since the category might then change as the surface becomes
        immutable.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p._test_category()

        """
        if self.surface().is_mutable():
            return

        super()._test_category(**options)

    def __hash__(self):
        r"""
        Return a hash value of this point.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (0, 0))
            sage: q = S.point(4, (0, 0))
            sage: hash(p) == hash(q)
            True

        """
        return hash(self._representatives)

    def __ne__(self, other):
        r"""
        Return whether this point is distinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p != p
            False

            sage: q = S.point(0, (1/2, 1/3))
            sage: p != q
            True

        """
        return not (self == other)


class SaddleConnection_base(SageObject):
    def __init__(self, surface):
        self._surface = surface


class SaddleConnection(SaddleConnection_base):
    r"""
    Represents a saddle connection on a SimilaritySurface.

    TESTS::

        sage: from flatsurf.geometry.saddle_connection import SaddleConnection
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.cathedral(1, 2)
        sage: SaddleConnection.from_vertex(S, 1, 8, (0, -1))
        Saddle connection (0, -2) from vertex 8 of polygon 1 to vertex 5 of polygon 1

    """

    # TODO: The constructor should not do so much work. Missing parameters
    # should be filled in by calling a static factory method.
    # TODO: direction and end_direction should not be part of the data. It's
    # just the tangent vector given by the holonomy.
    # TODO: start and end should probably be just a tangent vector? And the
    # naming of start() and end() should also reflect that? Then start() and
    # end() could actually just return the surface points instead of returning
    # a tuple.
    def __init__(
        self,
        surface,
        start,
        direction=None,
        end=None,
        end_direction=None,
        holonomy=None,
        end_holonomy=None,
        check=True,
        limit=None,
        start_data=None,
        end_data=None,
    ):
        r"""
        TODO: Cleanup documentation.

        Construct a saddle connection on a SimilaritySurface.

        The only necessary parameters are the surface, start_data, and direction
        (to start). If there is missing data that can not be inferred from the surface
        type, then a straight-line trajectory will be computed to confirm that this is
        indeed a saddle connection. The trajectory will pass through at most limit
        polygons before we give up.

        Details of the parameters are provided below.

        Parameters
        ----------
        surface : a SimilaritySurface
            which will contain the saddle connection being constructed.
        start : a pair
            consisting of the label of the polygon where the saddle connection starts
            and the starting vertex.
        direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction the saddle connection is moving in (in the
            coordinates of the initial polygon).
        end : a pair
            consisting of the label of the polygon where the saddle connection terminates
            and the terminating vertex.
        end_direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction to move backward from the end point (in the
            coordinates of the terminal polygon). If the surface is a DilationSurface
            or better this will be the negation of the direction vector. If the surface
            is a HalfDilation surface or better, then this will be either the direction
            vector or its negation. In either case the value can be inferred from the
            end.
        holonomy : 2-dimensional vector with entries in the base_ring of the surface
            the holonomy of the saddle connection measured from the start. To compute this
            you develop the saddle connection into the plane starting from the starting
            polygon.
        end_holonomy : 2-dimensional vector with entries in the base_ring of the surface
            the holonomy of the saddle connection measured from the end (with the opposite
            orientation). To compute this you develop the saddle connection into the plane
            starting from the terminating polygon. For a translation surface, this will be
            the negation of holonomy, and for a HalfTranslation surface it will be either
            equal to holonomy or equal to its negation. In both these cases the end_holonomy
            can be inferred and does not need to be passed to the constructor.
        check : boolean
            If all data above is provided or can be inferred, then when check=False this
            geometric data is not verified. With check=True the data is always verified
            by straight-line flow. Erroroneous data will result in a ValueError being thrown.
            Defaults to true.
        limit :
            The combinatorial limit (in terms of number of polygons crossed) to flow forward
            to check the saddle connection geometry.

        TESTS:

        Arguments are validated. If the direction points out of the polygon, no saddle connection can be created::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.cathedral(1, 2)
            sage: from flatsurf.geometry.saddle_connection import SaddleConnection
            sage: SaddleConnection(S, (1, 5), (1, 1))
            Traceback (most recent call last):
            ...
            ValueError: Singular point with vector pointing away from polygon

        """
        from flatsurf.geometry.categories import SimilaritySurfaces

        if surface not in SimilaritySurfaces():
            raise TypeError("surface must be a similarity surface")

        if start_data is not None:
            import warnings

            warnings.warn(
                "start_data has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use start instead"
            )
            start = start_data
            del start_data

        if start is None:
            raise ValueError("start must be specified to create a SaddleConnection")

        if end_data is not None:
            import warnings

            warnings.warn(
                "end_data has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use end instead"
            )
            end = end_data
            del end_data

        if direction is not None:
            import warnings

            if holonomy is None:
                warnings.warn(
                    "direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future of sage-flatsurf; if you want to create a SaddleConnection without specifying the holonomy, use SaddleConnection.from_vertex() instead."
                )
                c = SaddleConnection.from_vertex(
                    surface, *start, direction, limit=limit
                )
                start = c._start
                holonomy = c._holonomy
                end = c._end
                end_holonomy = c._end_holonomy
            else:
                warnings.warn(
                    "direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future of sage-flatsurf; there is no need to pass this argument anymore when the holonomy is specified."
                )
            del direction

        if holonomy is None:
            raise ValueError("holonomy must be specified to create a SaddleConnection")

        if end is None or end_holonomy is None:
            import warnings

            warnings.warn(
                "end and end_holonomy must be provided as a keyword argument for SaddleConnection() in future versions of sage-flatsurf; use SaddleConnection.from_vertex() instead to create a SaddleConnection without specifying these."
            )
            c = SaddleConnection.from_vertex(surface, *start, holonomy, limit=limit)
            start = c._start
            holonomy = c._holonomy
            end = c._end
            end_holonomy = c._end_holonomy

        if end_direction is not None:
            import warnings

            warnings.warn(
                "end_direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; the end direction is deduced from the end holonomy automatically"
            )
            del end_direction

        if limit is not None:
            import warnings

            warnings.warn(
                "limit has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use SaddleConnection.from_vertex() to search with a limit instead"
            )
            del limit

        super().__init__(surface)

        V = self._surface.base_ring() ** 2

        self._start = tuple(start)
        self._holonomy = V(holonomy)
        self._start, self._holonomy = self._normalize(self._start, self._holonomy)
        self._holonomy.set_immutable()

        self._end = tuple(end)
        self._end_holonomy = V(end_holonomy)
        self._end, self._end_holonomy = self._normalize(self._end, self._end_holonomy)
        self._end_holonomy.set_immutable()

    def surface(self):
        return self._surface

    def _normalize(self, start, holonomy):
        r"""
        Normalize the ``start`` and ``holonomy`` data describing this saddle
        connection.

        When the saddle connection is parallel to a polygon's edge, there can
        be two different descriptions of the same saddle connection.

        EXAMPLES::

            sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=((0, 0), (1, 0), (1, 1))))
            0

            sage: S.glue((0, 0), (0, 0))
            sage: S.glue((0, 1), (0, 1))
            sage: S.glue((0, 2), (0, 2))

            sage: S.set_immutable()

            sage: from flatsurf import SaddleConnection
            sage: SaddleConnection.from_vertex(surface=S, label=0, vertex=0, direction=(1, 0))
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 0 of polygon 0
            sage: SaddleConnection.from_vertex(surface=S, label=0, vertex=0, direction=(1, 1))
            Saddle connection (-1, -1) from vertex 2 of polygon 0 to vertex 2 of polygon 0

        ::

            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=((0, 0), (1, 0), (1, 1))))
            0

            sage: S.set_immutable()

            sage: from flatsurf import SaddleConnection
            sage: SaddleConnection(surface=S, start=(0, 0), direction=(1, 0))
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 1 of polygon 0
            sage: SaddleConnection(surface=S, start=(0, 0), direction=(1, 1))
            Saddle connection (1, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0

        """
        label = start[0]
        polygon = self._surface.polygon(label)
        previous_edge = (start[1] - 1) % len(polygon.vertices())
        if holonomy == -polygon.edge(previous_edge):
            opposite_edge = self._surface.opposite_edge(label, previous_edge)
            if opposite_edge is not None:
                return (
                    opposite_edge,
                    self._surface.edge_transformation(label, previous_edge).derivative()
                    * holonomy,
                )

        return start, holonomy

    def __neg__(self):
        return SaddleConnection(
            surface=self._surface,
            start=self._end,
            end=self._start,
            holonomy=self._end_holonomy,
            end_holonomy=self._holonomy,
            check=False,
        )

    @classmethod
    def from_half_edge(self, surface, label, edge):
        r"""
        Return a saddle connection along the ``edge`` in the polygon ``label``
        of ``surface``.

        INPUT:

        - ``surface`` -- a similarity surface

        - ``label`` -- a polygon label in ``surface``

        - ``edge`` -- the index of an edge in the polygon with ``label``

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SaddleConnection
            sage: S = translation_surfaces.square_torus()

            sage: SaddleConnection.from_half_edge(S, 0, 0)
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 2 of polygon 0

        Saddle connections in a surface with boundary::

            sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=((0, 0), (1, 0), (1, 1))))
            0
            sage: S.set_immutable()

            sage: c = SaddleConnection.from_half_edge(S, 0, 0); c
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 1 of polygon 0
            sage: -c
            Saddle connection (-1, 0) from vertex 1 of polygon 0 to vertex 0 of polygon 0

            sage: c == -c
            False

        Saddle connections in a surface with self-glued edges::

            sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=((0, 0), (1, 0), (1, 1))))
            0
            sage: S.glue((0, 0), (0, 0))
            sage: S.glue((0, 1), (0, 1))
            sage: S.glue((0, 2), (0, 2))

            sage: c = SaddleConnection.from_half_edge(S, 0, 0); c
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 0 of polygon 0
            sage: -c
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 0 of polygon 0

            sage: c == -c
            True

        """
        polygon = surface.polygon(label)
        holonomy = polygon.edge(edge)

        return SaddleConnection(
            surface=surface,
            start=(label, edge),
            end=(label, (edge + 1) % len(polygon.vertices())),
            holonomy=holonomy,
            end_holonomy=-holonomy,
            check=False,
        )

    @classmethod
    def from_vertex(cls, surface, label, vertex, direction, limit=None):
        r"""
        Return the saddle connection emanating from the ``vertex`` of the
        polygon with ``label`` following the ray ``direction``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SaddleConnection
            sage: S = translation_surfaces.square_torus()

            sage: SaddleConnection.from_vertex(S, 0, 0, (1, 0))
            Saddle connection (1, 0) from vertex 0 of polygon 0 to vertex 2 of polygon 0
            sage: SaddleConnection.from_vertex(S, 0, 0, (2, 1))
            Saddle connection (2, 1) from vertex 0 of polygon 0 to vertex 2 of polygon 0
            sage: SaddleConnection.from_vertex(S, 0, 0, (0, 1))
            Saddle connection (0, 1) from vertex 1 of polygon 0 to vertex 3 of polygon 0

        TESTS::

            sage: from flatsurf.geometry.saddle_connection import SaddleConnection
            sage: from flatsurf import translation_surfaces

            sage: S = translation_surfaces.cathedral(1, 2)
            sage: SaddleConnection.from_vertex(S, 1, 8, (0, -1))
            Saddle connection (0, -2) from vertex 8 of polygon 1 to vertex 5 of polygon 1

        """
        from flatsurf.geometry.ray import Rays

        R = Rays(surface.base_ring())
        direction = R(direction)

        tangent_vector = surface.tangent_vector(
            label, surface.polygon(label).vertex(vertex), direction.vector()
        )
        trajectory = tangent_vector.straight_line_trajectory()

        if limit is None:
            from sage.all import infinity

            limit = infinity

        trajectory.flow(steps=limit)

        if not trajectory.is_saddle_connection():
            raise ValueError(
                "no saddle connection in this direction within the specified limit"
            )

        end_tangent_vector = trajectory.terminal_tangent_vector()

        assert (
            not trajectory.segments()[0].is_edge() or len(trajectory.segments()) == 1
        ), "when the saddle connection is an edge it must not consist of more than that edge"

        if trajectory.segments()[0].is_edge():
            # When the saddle connection is just an edge, the similarity
            # computation below can be wrong when that edge is glued to the
            # same polygon. Namely, the similarity is missing the final factor
            # that comes from that gluing. E.g., in the Cathedral test case
            # above, the vertical saddle connection connecting the vertices of
            # polygon 1 at (1, 3) and (2, 1) by going in direction (0, -1) is
            # misinterpreted as the connection going in direction (1, -2).
            return SaddleConnection.from_half_edge(
                surface, label, trajectory.segments()[0].edge()
            )

        from flatsurf.geometry.similarity import SimilarityGroup

        one = SimilarityGroup(surface.base_ring()).one()
        segments = list(trajectory.segments())[1:]

        from sage.all import prod

        similarity = prod(
            [
                surface.edge_transformation(
                    segment.start().polygon_label(),
                    segment.start().position().get_edge(),
                )
                for segment in segments
            ],
            one,
        )

        holonomy = (
            similarity(trajectory.segments()[-1].end().point())
            - trajectory.initial_tangent_vector().point()
        )
        end_holonomy = (~similarity.derivative()) * holonomy

        return SaddleConnection(
            surface=surface,
            start=(label, vertex),
            holonomy=holonomy,
            end=(end_tangent_vector.polygon_label(), end_tangent_vector.vertex()),
            end_holonomy=end_holonomy,
        )

    @cached_method
    def direction(self):
        r"""
        Return a ray parallel to the :meth:`holonomy`.
        """
        from flatsurf.geometry.ray import Rays

        return Rays(self._holonomy.base_ring())(self._holonomy)

    @cached_method
    def end_direction(self):
        r"""
        Return a ray parallel to the :meth:`end_holonomy`.
        """
        from flatsurf.geometry.ray import Rays

        return Rays(self._holonomy.base_ring())(self._end_holonomy)

    def start_data(self):
        r"""
        Return the pair (l, v) representing the label and vertex of the corresponding polygon
        where the saddle connection originates.
        """
        import warnings

        warnings.warn(
            "start_data() has been deprecated and will be removed from a future version of sage-flatsurf; use start() instead."
        )

        return self.start()

    def start(self):
        # TODO: Document that surface()(*start()) produces the actual vertex.
        return self._start

    def end_data(self):
        r"""
        Return the pair (l, v) representing the label and vertex of the corresponding polygon
        where the saddle connection terminates.
        """
        import warnings

        warnings.warn(
            "end_data() has been deprecated and will be removed from a future version of sage-flatsurf; use end() instead."
        )

        return self.end()

    def end(self):
        # TODO: Document that surface()(*end()) produces the actual vertex.
        return self._end

    def holonomy(self):
        r"""
        Return the holonomy vector of the saddle connection (measured from the start).

        In a SimilaritySurface this notion corresponds to developing the saddle connection into the plane
        using the initial chart coming from the initial polygon.
        """
        return self._holonomy

    def length(self):
        r"""
        In a cone surface, return the length of this saddle connection. Since
        this may not lie in the field of definition of the surface, it is
        returned as an element of the Algebraic Real Field.
        """
        from flatsurf.geometry.categories import ConeSurfaces

        if self._surface not in ConeSurfaces():
            raise NotImplementedError(
                "length of a saddle connection only makes sense for cone surfaces"
            )

        from sage.all import vector, AA

        return vector(AA, self._holonomy).norm()

    def length_squared(self):
        holonomy = self.holonomy()
        return holonomy[0] ** 2 + holonomy[1] ** 2

    def end_holonomy(self):
        r"""
        Return the holonomy vector of the saddle connection (measured from the end).

        In a SimilaritySurface this notion corresponds to developing the saddle connection into the plane
        using the initial chart coming from the initial polygon.
        """
        return self._end_holonomy

    def start_tangent_vector(self):
        r"""
        Return a tangent vector to the saddle connection based at its
        :meth:`start`.
        """
        return self._surface.tangent_vector(
            self._start[0],
            self._surface.polygon(self._start[0]).vertex(self._start[1]),
            self.direction().vector(),
        )

    @cached_method(key=lambda self, limit, cache: None)
    def trajectory(self, limit=1000, cache=None):
        r"""
        Return a straight line trajectory representing this saddle connection.
        Fails if the trajectory passes through more than limit polygons.
        """
        if cache is not None:
            import warnings

            warnings.warn(
                "The cache keyword argument of trajectory() is ignored. Trajectories are always cached."
            )

        v = self.start_tangent_vector()
        traj = v.straight_line_trajectory()
        traj.flow(limit)
        if not traj.is_saddle_connection():
            raise ValueError(
                "Did not obtain saddle connection by flowing forward. Limit="
                + str(limit)
            )

        return traj

    def plot(self, *args, **options):
        r"""
        Equivalent to ``.trajectory().plot(*args, **options)``
        """
        return self.trajectory().plot(*args, **options)

    def end_tangent_vector(self):
        r"""
        Return a tangent vector to the saddle connection based at its
        :meth:`end`.
        """
        return self._surface.tangent_vector(
            self._end[0],
            self._surface.polygon(self._end[0]).vertex(self._end[1]),
            self._end_direction.vector(),
        )

    def invert(self):
        r"""
        Return this saddle connection but with opposite orientation.
        """
        return SaddleConnection(
            self._surface,
            self._end,
            self._end_direction,
            self._start,
            self._direction,
            self._end_holonomy,
            self._holonomy,
            check=False,
        )

    def intersections(self, traj, count_singularities=False, include_segments=False):
        r"""
        See documentation of :meth:`~.straight_line_trajectory.AbstractStraightLineTrajectory.intersections`
        """
        return self.trajectory().intersections(
            traj, count_singularities, include_segments
        )

    def intersects(self, traj, count_singularities=False):
        r"""
        See documentation of :meth:`~.straight_line_trajectory.AbstractStraightLineTrajectory.intersects`
        """
        return self.trajectory().intersects(
            traj, count_singularities=count_singularities
        )

    def __eq__(self, other):
        r"""
        Return whether this saddle connection is indistinguishable from
        ``other``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: connections = S.saddle_connections(13)  # random output due to deprecation warning from cppyy

            sage: connections[0] == connections[0]
            True
            sage: connections[0] == connections[1]
            False


        TESTS:

        Verify that saddle connections can be compared to arbitrary objects (so
        they can be put into dicts with other objects)::

            sage: connections[0] == 42
            False

        ::

            sage: len(connections)
            32
            sage: len(set(connections))
            32

        """
        if self is other:
            return True

        if not isinstance(other, SaddleConnection):
            return False

        if self._surface != other._surface:
            return False

        if self._start != other._start:
            return False

        if self._holonomy != other._holonomy:
            return False

        return True

    def __hash__(self):
        return hash((self._start, self._holonomy))

    def _test_geometry(self, **options):
        # Test that this saddle connection actually exists on the surface.
        SaddleConnection(
            self._surface,
            self._start,
            self._direction,
            self._end,
            self._end_direction,
            self._holonomy,
            self._end_holonomy,
            check=True,
        )

    def __repr__(self):
        return f"Saddle connection {self.holonomy()} from vertex {self.start()[1]} of polygon {self.start()[0]} to vertex {self.end()[1]} of polygon {self.end()[0]}"

    def _test_inverse(self, **options):
        # Test that inverting works properly.
        SaddleConnection(
            self._surface,
            self._end,
            self._end_direction,
            self._start,
            self._direction,
            self._end_holonomy,
            self._holonomy,
            check=True,
        )

    def is_closed(self):
        return self.surface()(*self.start()) == self.surface()(*self.end())


    def _homology_(self, H):
        r"""
        Return this saddle connection as a chain of edges in the homology group ``H``.

        EXAMPLES::

            sage: from flatsurf import *
            sage: S = translation_surfaces.mcmullen_L(1,1,1,1)
            sage: H = S.homology()
            sage: holonomy = lambda x: sum(coeff * S.polygon(label).edge(e) for (label, e), coeff in dict(x._chain).items())
            sage: for sc in S.saddle_connections(5):
            ....:     h = H(sc)
            ....:     hol1 = holonomy(h)
            ....:     hol2 = sc.holonomy()
            ....:     assert hol1 == hol2, (sc, h, hol1, hol2)
        """
        from .homology import SimplicialHomologyGroup

        if not isinstance(H, SimplicialHomologyGroup):
            raise TypeError("H must be a SimplicialHomologyGroup")

        surface = self._surface
        if H._surface != surface:
            raise ValueError("homology and surface do not match")

        traj = self.start_tangent_vector().straight_line_trajectory()
        traj.flow(Infinity)
        segments = traj.segments()

        if len(segments) == 1 and segments[0].is_edge():
            # NOTE: in the special case the saddle connection is an edge,
            # the behavior of the corresponding segment is not appropriate
            # to the generic code afterwards. Namely, the start and the end
            # belongs to two different polygons. See
            # https://github.com/flatsurf/sage-flatsurf/issues/309
            label = segments[0].polygon_label()
            e = segments[0].edge()
            return H((label, e))

        h = H.zero()
        for s in segments:
            label = s.polygon_label()
            n = len(surface.polygon(label).vertices())

            start = s.start()
            if start.position().is_in_edge_interior():
                # pick the next vertex clockwise
                i = (start.position().get_edge() + 1) % n
            else:
                assert start.position().is_vertex()
                i = start.vertex()

            end = s.end()
            if end.position().is_in_edge_interior():
                # pick the previous vertex clockwise
                j = end.position().get_edge()
            else:
                assert end.position().is_vertex()
                j = end.vertex()

            if i < j:
                h += sum(H((label, e)) for e in range(i, j))
            else:
                h += sum(H((label, e)) for e in range(i, n))
                h += sum(H((label, e)) for e in range(j))

        return h


class Cylinder(SageObject):
    r"""
    Represents a cylinder in a SimilaritySurface. A cylinder for these purposes is a
    topological annulus in a surface bounded by a finite collection of saddle connections
    meeting at 180 degree angles.

    To Do
    -----
    * Support cylinders whose monodromy is a dilation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s = translation_surfaces.octagon_and_squares()
        sage: from flatsurf.geometry.surface_objects import Cylinder
        sage: cyl = Cylinder(s, 0, [2, 3, 3, 3, 2, 0, 1, 3, 2, 0])
        sage: cyl.initial_label()
        0
        sage: cyl.edges()
        (2, 3, 3, 3, 2, 0, 1, 3, 2, 0)
        sage: # a = sqrt(2) below.
        sage: cyl.area()
        2*a + 4
        sage: cyl.circumference().minpoly()
        x^4 - 680*x^2 + 400
        sage: cyl.holonomy()
        (8*a + 12, 4*a + 6)
    """

    def __init__(self, s, label0, edges):
        r"""
        Construct a cylinder on the surface `s` from an initial label and a
        sequence of edges crossed.

        Parameters
        ----------
        s: A SimilaritySurface
            the surface containing the cylinder
        label0: An initial label
            representing a polygon the cylinder passes through.
        edges: a list
            giving the sequence of edges the cylinder crosses until it closes.
        """
        self._surface = s
        self._label0 = label0
        self._edges = tuple(edges)
        ss = s.minimal_cover(cover_type="planar")
        SG = SimilarityGroup(s.base_ring())
        labels = [(label0, SG.one())]  # labels of polygons on the cover ss.
        for e in edges:
            labels.append(ss.opposite_edge(labels[-1], e)[0])
        if labels[0][0] != labels[-1][0]:
            raise ValueError("Combinatorial path does not close.")
        trans = labels[-1][1]
        if not trans.is_translation():
            raise NotImplementedError(
                "Only cylinders with translational monodromy are currently supported"
            )
        m = trans.matrix()
        v = vector(s.base_ring(), (m[0][2], m[1][2]))  # translation vector
        from flatsurf.geometry.euclidean import ccw

        p = ss.polygon(labels[0])
        e = edges[0]
        min_y = ccw(v, p.vertex(e))
        max_y = ccw(v, p.vertex((e + 1) % len(p.vertices())))
        if min_y >= max_y:
            raise ValueError("Combinatorial data does not represent a cylinder")

        # Stores the vertices where saddle connections starts:
        min_list = [0]
        max_list = [0]

        for i in range(1, len(edges)):
            e = edges[i]
            p = ss.polygon(labels[i])
            y = ccw(v, p.vertex(e))
            if y == min_y:
                min_list.append(i)
            elif y > min_y:
                min_list = [i]
                min_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")
            y = ccw(v, p.vertex((e + 1) % len(p.vertices())))
            if y == max_y:
                max_list.append(i)
            elif y < max_y:
                max_list = [i]
                max_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")

        # Extract the saddle connections on the right side:
        from flatsurf.geometry.surface_objects import SaddleConnection

        sc_set_right = set()
        vertices = []
        for i in min_list:
            label = labels[i]
            p = ss.polygon(label)
            vertices.append((i, p.vertex(edges[i])))
        i, vert_i = vertices[-1]
        vert_i = vert_i - v
        j, vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v) * li[1])
            lio = ss.opposite_edge(li, edges[i])
            lj = labels[j]
            sc = SaddleConnection(
                s,
                (lio[0][0], (lio[1] + 1) % len(ss.polygon(lio[0]).vertices())),
                (~lio[0][1])(vert_j) - (~lio[0][1])(vert_i),
            )
            sc_set_right.add(sc)
        i = j
        vert_i = vert_j
        for j, vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                li = (li[0], SG(-v) * li[1])
                lio = ss.opposite_edge(li, edges[i])
                lj = labels[j]
                sc = SaddleConnection(
                    s,
                    (lio[0][0], (lio[1] + 1) % len(ss.polygon(lio[0]).vertices())),
                    (~lio[0][1])(vert_j) - (~lio[0][1])(vert_i),
                    limit=j - i,
                )
                sc_set_right.add(sc)
            i = j
            vert_i = vert_j

        # Extract the saddle connections on the left side:
        sc_set_left = set()
        vertices = []
        for i in max_list:
            label = labels[i]
            p = ss.polygon(label)
            vertices.append((i, p.vertex((edges[i] + 1) % len(p.vertices()))))
        i, vert_i = vertices[-1]
        vert_i = vert_i - v
        j, vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v) * li[1])
            lio = ss.opposite_edge(li, edges[i])
            lj = labels[j]
            sc = SaddleConnection(
                s,
                (lj[0], (edges[j] + 1) % len(ss.polygon(lj).vertices())),
                (~lj[1])(vert_i) - (~lj[1])(vert_j),
            )
            sc_set_left.add(sc)
        i = j
        vert_i = vert_j
        for j, vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                lio = ss.opposite_edge(li, edges[i])
                lj = labels[j]
                sc = SaddleConnection(
                    s,
                    (lj[0], (edges[j] + 1) % len(ss.polygon(lj).vertices())),
                    (~lj[1])(vert_i) - (~lj[1])(vert_j),
                )
                sc_set_left.add(sc)
            i = j
            vert_i = vert_j
        self._boundary1 = frozenset(sc_set_right)
        self._boundary2 = frozenset(sc_set_left)
        self._boundary = frozenset(self._boundary1.union(self._boundary2))

        edge_intersections = []
        i = min_list[0]
        label = labels[i]
        p = ss.polygon(label)
        right_point = p.vertex(edges[i])  # point on the right boundary
        i = max_list[0]
        label = labels[i]
        p = ss.polygon(label)
        left_point = p.vertex((edges[i] + 1) % len(p.vertices()))
        from flatsurf.geometry.euclidean import solve

        for i in range(len(edges)):
            label = labels[i]
            p = ss.polygon(label)
            e = edges[i]
            v1 = p.vertex(e)
            v2 = p.vertex((e + 1) % len(p.vertices()))
            a, b = solve(left_point, v, v1, v2 - v1)
            w1 = (~(label[1]))(v1 + b * (v2 - v1))
            a, b = solve(right_point, v, v1, v2 - v1)
            w2 = (~(label[1]))(v1 + b * (v2 - v1))
            edge_intersections.append((w1, w2))

        polygons = []
        pair1 = edge_intersections[-1]
        l1 = labels[-2][0]
        e1 = edges[-1]
        for i in range(len(edges)):
            l2 = labels[i][0]
            pair2 = edge_intersections[i]
            e2 = edges[i]
            trans = s.edge_transformation(l1, e1)
            pair1p = (trans(pair1[0]), trans(pair1[1]))
            polygon_verts = [pair1p[0], pair1p[1]]
            if pair2[1] != pair1p[1]:
                polygon_verts.append(pair2[1])
            if pair2[0] != pair1p[0]:
                polygon_verts.append(pair2[0])

            from flatsurf import Polygon

            polygons.append(
                (l2, Polygon(vertices=polygon_verts, base_ring=s.base_ring()))
            )
            l1 = l2
            pair1 = pair2
            e1 = e2
        self._polygons = tuple(polygons)

    def surface(self):
        return self._surface

    def initial_label(self):
        r"""
        Return one label on the surface that the cylinder passes through.
        """
        return self._label0

    def edges(self):
        r"""
        Return a tuple of edge numbers representing the edges crossed
        when the cylinder leaves the polygon with `initial_label` until
        it returns by closing.
        """
        return self._edges

    def boundary(self):
        r"""
        Return the set of saddle connections in the boundary, oriented so that
        the surface is on the left.
        """
        return self._boundary

    def polygons(self):
        r"""
        Return a list of pairs each consisting of a label and a polygon.
        Each polygon represents a sub-polygon of the polygon on the surface
        with the given label. The union of these sub-polygons form the
        cylinder. The subpolygons are listed in cyclic order.
        """
        return self._polygons

    @cached_method
    def area(self):
        r"""
        Return the area of this cylinder if it is contained in a ConeSurface.
        """
        from flatsurf.geometry.categories import ConeSurfaces

        if self._surface not in ConeSurfaces():
            raise NotImplementedError("area only makes sense for cone surfaces")

        area = 0
        for label, p in self.polygons():
            area += p.area()
        return area

    def plot(self, **options):
        r"""
        Plot this cylinder in coordinates used by a graphical surface. This
        plots this cylinder as a union of subpolygons. Only the intersections
        with polygons visible in the graphical surface are shown.

        Parameters other than `graphical_surface` are passed to `polygon2d`
        which is called to render the polygons.

        Parameters
        ----------
        graphical_surface : a GraphicalSurface
            If not provided or `None`, the plot method uses the default graphical
            surface for the surface.
        """
        if "graphical_surface" in options and options["graphical_surface"] is not None:
            gs = options["graphical_surface"]
            if gs.get_surface() != self._surface:
                raise ValueError("graphical surface for the wrong surface")
            del options["graphical_surface"]
        else:
            gs = self._surface.graphical_surface()
        plt = Graphics()
        for label, p in self.polygons():
            if gs.is_visible(label):
                gp = gs.graphical_polygon(label)
                t = gp.transformation()
                pp = t(p)
                poly = polygon2d(pp.vertices(), **options)
                plt += poly.plot()
        return plt

    @cached_method
    def labels(self):
        r"""
        Return the set of labels that this cylinder passes through.
        """
        polygons = self.polygons()
        return frozenset([label for label, p in polygons])

    def boundary_components(self):
        r"""
        Return a set of two elements: the set of saddle connections on
        the right and left sides. Saddle connections are oriented so that
        the surface is on the left.
        """
        return frozenset([self._boundary1, self._boundary2])

    def next(self, sc):
        r"""
        Return the next saddle connection as you move around the cylinder boundary
        moving from sc in the direction of its orientation.
        """
        if sc not in self._boundary:
            raise ValueError

        v = sc.end_tangent_vector()
        v = v.clockwise_to(-v.vector())
        from flatsurf.geometry.euclidean import is_parallel

        for sc2 in self._boundary:
            if sc2.start_data() == (
                v.polygon_label(),
                v.vertex(),
            ) and is_parallel(sc2.direction(), v.vector()):
                return sc2
        raise ValueError("Failed to find next saddle connection in boundary set.")

    def previous(self, sc):
        r"""
        Return the previous saddle connection as you move around the cylinder boundary
        moving from sc in the direction opposite its orientation.
        """
        if sc not in self._boundary:
            raise ValueError
        v = sc.start_tangent_vector()
        v = v.counterclockwise_to(-v.vector())
        from flatsurf.geometry.euclidean import is_parallel

        for sc2 in self._boundary:
            if sc2.end_data() == (v.polygon_label(), v.vertex()) and is_parallel(
                sc2.end_direction(), v.vector()
            ):
                return sc2
        raise ValueError("Failed to find previous saddle connection in boundary set.")

    @cached_method
    def holonomy(self):
        r"""
        In a translation surface, return one of the two holonomy vectors of the cylinder,
        which differ by a sign.
        """
        from flatsurf.geometry.categories import TranslationSurfaces

        if self._surface not in TranslationSurfaces():
            raise NotImplementedError(
                "holonomy currently only computable for translation surfaces"
            )

        V = self._surface.base_ring() ** 2
        total = V.zero()
        for sc in self._boundary1:
            total += sc.holonomy()

        # Debugging:
        total2 = V.zero()
        for sc in self._boundary2:
            total2 += sc.holonomy()
        assert (
            total + total2 == V.zero()
        ), "Holonomy of the two boundary components should sum to zero."

        return total

    @cached_method
    def circumference(self):
        r"""
        In a cone surface, return the circumference, i.e., the length
        of a geodesic loop running around the cylinder. Since this may
        not lie in the field of definition of the surface, it is returned
        as an element of the Algebraic Real Field.
        """
        from flatsurf.geometry.categories import ConeSurfaces

        if self._surface not in ConeSurfaces():
            raise NotImplementedError(
                "circumference only makes sense for cone surfaces"
            )

        total = 0
        for sc in self._boundary1:
            total += sc.length()
        return total
