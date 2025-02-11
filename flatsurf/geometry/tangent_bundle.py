r"""
.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks
"""

# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2022 W. Patrick Hooper
#                      2016-2022 Vincent Delecroix
#                      2022-2025 Julian Rüth
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

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.misc.cachefunc import cached_method

rotate_limit = 100


# TODO: Rename to make this inaccessible in a breaking change.
class TangentVector(Element):
    r"""
    Tangent vector to a similarity surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces

    Examples on edges in direction of edges::

        sage: s = translation_surfaces.square_torus()
        sage: s.tangent_vector(0, (1/2, 0), (1, 0))
        (1, 0) at (1/2, 0) of polygon 0
        sage: s.tangent_vector(0, (1/2, 0), (-1, 0))
        (-1, 0) at (1/2, 1) of polygon 0
        sage: s.tangent_vector(0, (1/2, 1), (1, 0))
        (1, 0) at (1/2, 0) of polygon 0
        sage: s.tangent_vector(0, (1/2, 1), (-1, 0))
        (-1, 0) at (1/2, 1) of polygon 0

        sage: s.tangent_vector(0, (0, 1/2), (0, 1))
        (0, 1) at (1, 1/2) of polygon 0
        sage: s.tangent_vector(0, (0, 1/2), (0, -1))
        (0, -1) at (0, 1/2) of polygon 0
        sage: s.tangent_vector(0, (1, 1/2), (0, 1))
        (0, 1) at (1, 1/2) of polygon 0
        sage: s.tangent_vector(0, (1, 1/2), (0, -1))
        (0, -1) at (0, 1/2) of polygon 0

    Examples on vertices in direction of edges::

        sage: s = translation_surfaces.square_torus()
        sage: s.tangent_vector(0, (0, 0), (1, 0))
        (1, 0) at vertex 0 of polygon 0
        sage: s.tangent_vector(0, (1, 0), (-1, 0))
        (-1, 0) at vertex 2 of polygon 0
        sage: s.tangent_vector(0, (0, 1), (1, 0))
        (1, 0) at vertex 0 of polygon 0
        sage: s.tangent_vector(0, (1, 1), (-1, 0))
        (-1, 0) at vertex 2 of polygon 0

        sage: s.tangent_vector(0, (0, 0), (0, 1))
        (0, 1) at vertex 1 of polygon 0
        sage: s.tangent_vector(0, (0, 1), (0, -1))
        (0, -1) at vertex 3 of polygon 0
        sage: s.tangent_vector(0, (1, 0), (0, 1))
        (0, 1) at vertex 1 of polygon 0
        sage: s.tangent_vector(0, (1, 1), (0, -1))
        (0, -1) at vertex 3 of polygon 0

    TESTS:

    Verify that tangent vectors can be formed in surfaces with non-convex
    polygons::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
        sage: L = Polygon(vertices=[(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)])
        sage: S = MutableOrientedSimilaritySurface(QQ)
        sage: S.add_polygon(L)
        0
        sage: S.set_immutable()

        sage: S.tangent_vector(0, (0, 0), (1, 1))
        (1, 1) at vertex 0 of polygon 0
        sage: S.tangent_vector(0, (1, 1), (-1, -1))
        (-1, -1) at vertex 3 of polygon 0
        sage: S.tangent_vector(0, (1, 1), (-1, 1))
        (-1, 1) at vertex 3 of polygon 0
        sage: S.tangent_vector(0, (1, 1), (1, -1))
        (1, -1) at vertex 3 of polygon 0

    """

    def __init__(self, tangent_bundle, polygon_label, point, vector):
        super().__init__(tangent_bundle)

        # TODO: Use parent instead
        self._bundle = tangent_bundle
        # TODO: Rename to label
        self._polygon_label = polygon_label

        self._point = (tangent_bundle.base_ring()**2)(point)
        self._point.set_immutable()

        self._vector = (tangent_bundle.base_ring()**2)(vector)
        self._vector.set_immutable()

    def __repr__(self):
        return f"{self.vector()!r} at {self.surface()(*self.base_point())._repr_representative(*self.base_point(), uppercase=False, shortened=True)}"

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, TangentVector):
            return False
        if self.parent() != other.parent():
            return False
        if self._polygon_label != other._polygon_label:
            return False
        if self._point != other._point:
            return False
        if self._vector != other._vector:
            return False

        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.square_torus()
            sage: for y in [0,1]:
            ....:     for d in [1,-1]:
            ....:         h = hash(s.tangent_vector(0, (1/2, y), (d, 0)))
        """
        return hash((self._bundle, self._polygon_label, self._point, self._vector))

    def surface(self):
        r"""Return the underlying surface."""
        return self._bundle.surface()

    # TODO: Is this a good interface?
    def is_based_at_singularity(self):
        r"""
        Return the truth value of the statement 'the base point for this vector is a singularity.'
        """
        return self.position().is_vertex()

    # TODO: Is this a good interface?
    def vertex(self):
        r"""Return the index of the vertex."""
        return self.position().get_vertex()

    # TODO: Is this a good interface?
    def is_in_boundary_of_polygon(self):
        r"""
        Return the truth value of the statement
        'the base point for this vector lies on the boundary of
        one of the polygons making up the surface.'
        """
        return self.position().is_in_boundary()

    # TODO: Hide as an internal method.
    @cached_method
    def position(self):
        r"""
        Return the PolygonPosition representing the location of
        the basepoint of the vector in the polygon that contains it.
        """
        return self.surface().polygon(self._polygon_label).get_point_position(self._point)

    # TODO: Deprecate to parent()
    def bundle(self):
        r"""Return the tangent bundle containing this vector."""
        return self._bundle

    # TODO: Deprecate to access this through the base point.
    def polygon_label(self):
        return self._polygon_label

    # TODO: Deprecate to do this manually.
    def polygon(self):
        return self.surface().polygon(self.polygon_label())

    # TODO: Deprecate. What's the good interface?
    def point(self):
        r"""
        Return the base point of this tangent vector as a vector.

        The coordinates of output are given with respect to the polygon it
        belongs to.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces

            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (1/2,0), (0,1))
            sage: v.point()
            (1/2, 0)
            sage: parent(_)
            Vector space of dimension 2 over Rational Field
        """
        return self._point

    # TODO: Maybe add a representative=False here?
    def base_point(self):
        r"""
        Return the base point of this tangent vector as a pair ``(label,
        coordinates)``.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces

            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (1/2, 0), (0, 1))
            sage: v.base_point()
            (1, (1/2, 0))
            sage: s(*v.base_point())
            Point (1/2, 0) of polygon 0
            sage: _.representatives()
            frozenset({(0, (1/2, 0)), (1, (1/2, 0))})

        """
        return (self._polygon_label, self._point)

    def vector(self):
        r"""
        Return the coordinates of this vector within the assigned polygon.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces

            sage: s = similarity_surfaces.example()
            sage: v = s.tangent_vector(0, (1/2,0), (0,1))
            sage: v.vector()
            (0, 1)
            sage: parent(_)
            Vector space of dimension 2 over Rational Field
        """
        return self._vector

    # TODO: Is this a good interface?
    def edge_pointing_along(self):
        r"""
        Returns the pair of (p,e) where p is the polygon label at the base point,
        and e is the edge this vector points along or none if it does not point
        along an edge. Here pointing along means that the vector is based at
        a vertex and represents the vector joining this edge to the next vertex."""
        if self.is_based_at_singularity():
            e = self.vertex()
            if self.vector() == self.polygon().edge(e):
                return (self.polygon_label(), e)
        return None

    # TODO: Is this a good interface?
    def differs_by_scaling(self, another_tangent_vector):
        r"""
        Returns true if the other vector just differs by scaling. This means they should lie
        in the same polygon, be based at the same point, and point in the same direction.
        """
        from flatsurf.geometry.euclidean import is_parallel

        return (
            self.polygon_label() == another_tangent_vector.polygon_label()
            and self.point() == another_tangent_vector.point()
            and is_parallel(self.vector(), another_tangent_vector.vector())
        )

    # TODO: Deprecate to operator-
    def invert(self):
        r"""
        Returns the negation of this tangent vector.
        Raises a ValueError if the vector is based at a singularity.'
        """
        if self.is_based_at_singularity():
            raise ValueError("Can't invert tangent vector based at a singularity.")
        return self._bundle(
            self.polygon_label(), self.point(), -self.vector()
        )

    # TODO: Is this a good interface?
    def forward_to_polygon_boundary(self):
        r"""
        Flows forward (in the direction of the tangent vector) until the end
        of the polygon is reached.
        Returns the tangent vector based at the endpoint which point backward along the trajectory.

        NOTES::

            We return the backward trajectory, because continuing forward does not make sense if a
            singularity is reached. You can obtain the forward vector by subsequently applying invert().

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()

            sage: v = s.tangent_vector(0, (0,0), (3,-1))
            sage: v
            (3, -1) at vertex 0 of polygon 0
            sage: v2 = v.forward_to_polygon_boundary()
            sage: v2
            (-3, 1) at (2, -2/3) of polygon 0
            sage: v2.invert()
            (4, -3) at (2/3, 2) of polygon 1
        """
        p = self.polygon()
        point2, pos2 = p.flow_to_exit(self.point(), self.vector())
        new_vector = self._bundle(
            self.polygon_label(), point2, -self.vector()
        )
        return new_vector

    def straight_line_trajectory(self):
        r"""
        Return the straight line trajectory associated to this vector.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: s = translation_surfaces.square_torus()
            sage: v = s.tangent_vector(0, (0, 0), (1, 1))
            sage: v.straight_line_trajectory()
            Straight line trajectory made of 1 segments from (0, 0) in polygon 0 to (1, 1) in polygon 0
            sage: l = v.straight_line_trajectory()
            sage: l
            Straight line trajectory made of 1 segments from (0, 0) in polygon 0 to (1, 1) in polygon 0
            sage: l.is_saddle_connection()
            True

            sage: v = s.tangent_vector(0, (0, 0), (1, 1 + AA(5).sqrt()), ring=AA)
            doctest:warning
            ...
            UserWarning: the ring parameter has been deprecated in tangent_bundle(); call change_ring() on the underlying surface instead
            sage: l = v.straight_line_trajectory()
            sage: l.flow(20)
            sage: l.segment(20)
            Segment in polygon 0 starting at (0.9442719099991588?, 0) and ending at (1, 0.1803398874989485?)
        """
        from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectory

        return StraightLineTrajectory(self)

    def rotate(self, direction, ccw=True, preserve_scaling=True):
        r"""
        Return this tangent vector rotated so that it is parallel with
        ``direction``.

        INPUT:

        - ``direction`` -- a non-zero vector

        - ``ccw`` -- a boolean (default: ``True``); whether the rotation should
          be counterclockwise or clockwise

        - ``preserve_scaling`` -- a boolean (default: ``True); whether to
          return a tangent vector with equal length or any positive multiple
          thereof

        EXAMPLES:

        When ``direction`` is parallel with the tangent vector, no rotation is
        performed::

            sage: from flatsurf import translation_surfaces

            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)
            sage: t = S.tangent_vector(0, (0, 0), (1, 1))
            sage: t
            (1, 1) at vertex 0 of polygon 0
            sage: t.rotate((1, 1))
            (1, 1) at vertex 0 of polygon 0

        You can force a rotation in such a case by rotating twice by an angle π::

            sage: t.rotate((-1, -1)).rotate((1, 1))
            (1, 1) at vertex 0 of polygon 1

        The result might not be representable over the base field without
        taking some square roots. In that case, ``preserve_scaling`` can be set
        to return a positive multiple of the rotated tangent vector::

            sage: t = S.tangent_vector(0, (0, 0), (1, 0))
            sage: t
            (1, 0) at vertex 0 of polygon 0
            sage: t.rotate((1, 1))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(1/2) to a rational
            sage: t.rotate((1, 1), preserve_scaling=False)
            (1, 1) at vertex 0 of polygon 0

        The scaling of ``direction`` has no effect on the rotation. The length
        of the tangent vector is preserved::

            sage: t.rotate((0, 1))
            (0, 1) at vertex 1 of polygon 2
            sage: t.rotate((0, 2))
            (0, 1) at vertex 1 of polygon 2

        However, in a dilation surface, apparent scaling can happen when
        rotating across an edge::

            sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 2)]))
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.glue((0, 1), (0, 3))
            sage: S.set_immutable()

            sage: t = S.tangent_vector(0, (0, 0), (1, 1))
            sage: t
            (1, 1) at vertex 0 of polygon 0
            sage: t.rotate((-1, 1))
            (-1/2, 1/2) at vertex 1 of polygon 0

        The ``direction`` is always made sense of in the polygon in which the
        tangent vector is currently represented. In non-translation surfaces
        this can be confusing, e.g., when trying to rotate a tangent vector by
        an angle π, make sure to pick ``direction`` relative to the correct
        polygon::

            sage: t = S.tangent_vector(0, (0, 2), (1, -1))
            sage: t.rotate((-1, 1), preserve_scaling=False)
            (-1/2, 1/2) at vertex 1 of polygon 0
            sage: t
            (1, 0) at vertex 0 of polygon 0

        In other words, the angle of rotation that is needed is determined by
        comparing ``direction`` to the :meth:`vector` of this tangent vector.
        That rotation is then performed and the resulting :meth:`vector` might
        not be parallel with ``direction`` if edges are glued with rotations::

            sage: t.rotate((-1, -1), preserve_scaling=False)
            (-1, 0) at vertex 2 of polygon 0

        """
        surface = self.surface()

        direction = (surface.base_ring()**2)(direction)
        if not direction:
            raise ValueError("direction must be non-zero")

        if preserve_scaling:
            from flatsurf.geometry.euclidean import rotate
        else:
            def rotate(v, direction): return direction

        label = self.polygon_label()
        polygon = surface.polygon(label)
        point = self.point()
        vector = self.vector()
        sector_start = vector

        if not self.position().is_vertex():
            return self._bundle(label, point, rotate(vector, direction))

        # TODO: Maybe use edge_ccw here instead of rolling our own.
        vertex = self.position().get_vertex()
        while True:
            point = polygon.vertex(vertex)
            edge = (vertex - 1 if ccw else vertex) % len(polygon.edges())
            sector_end = polygon.edge(edge)
            if ccw:
                sector_end = -sector_end

            # We determine whether the target direction is between the vector
            # and the following edge.
            # If it is, the rotated tangent vector is based in this polygon,
            # otherwise, we need to turn into the next polygon.

            from flatsurf.geometry.euclidean import is_between
            if is_between(sector_start if ccw else sector_end, sector_end if ccw else sector_start, direction, strict=False):
                return self._bundle(label, point, rotate(vector, direction))

            opposite_label, opposite_edge = surface.opposite_edge(label, edge)
            direction = surface.edge_matrix(label, edge) * direction
            vector = surface.edge_matrix(label, edge) * vector
            sector_start = surface.edge_matrix(label, edge) * sector_end
            label = opposite_label
            polygon = surface.polygon(label)
            vertex = (opposite_edge if ccw else opposite_edge + 1) % len(polygon.edges())

    def clockwise_to(self, w, code=False):
        r"""
        Return the new tangent vector obtained by rotating this one in the clockwise
        direction until the vector is parallel to w, and scaling so that the length matches
        that of w.

        Note that we always do some rotation so that if w is parallel to this vector, then a
        -360 degree rotation is performed.

        The vector w must be nonzero.

        On an infinite surface, this is potentially an infinite calculation
        so we impose a limit (representing the maximal number of polygons
        that must be rotated through). This is the variable rotate_limit
        in this package.

        If code is True, we compute the sequences of numbers associated to edges
        crossed as a list. We return a pair consisting of the newly computing
        tangent vector an this code. This is currently only implemented when
        based at a singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.regular_octagon()
            sage: v = s.tangent_vector(0, (0, 0), (1, 1))
            sage: v.clockwise_to((-1, -1))
            doctest:warning
            ...
            UserWarning: clockwise_to() has been deprecated and will be removed in a future version of sage-flatsurf; use rotate(w, ccw=False, preserve_scaling=False) instead
            (-1, -1) at vertex 5 of polygon 0
            sage: v.clockwise_to((1, 1))
            (1, 1) at vertex 7 of polygon 0
            sage: v.clockwise_to((1, 1), code=True)
            ((1, 1) at vertex 7 of polygon 0, [0, 5, 2])
        """
        import warnings
        warnings.warn("clockwise_to() has been deprecated and will be removed in a future version of sage-flatsurf; use rotate(w, ccw=False, preserve_scaling=False) instead")

        if not w:
            raise ValueError("w must be non-zero")

        if not self.surface().is_translation_surface():
            raise NotImplementedError("clockwise_to() only implemented on translation surfaces")

        if self.is_based_at_singularity():
            s = self.surface()
            v1 = self.vector()
            label = self.polygon_label()
            vertex = self.vertex()
            v2 = s.polygon(label).edge(vertex)
            from sage.matrix.constructor import Matrix

            der = Matrix(s.base_ring(), [[1, 0], [0, 1]])
            if code:
                codes = []

            from flatsurf.geometry.euclidean import ccw

            for count in range(rotate_limit):
                if ccw(v2, w) >= 0 and ccw(w, v1) > 0:
                    # We've found it!
                    break
                if code:
                    codes.append(vertex)
                label2, edge2 = s.opposite_edge(label, vertex)
                der = der * s.edge_matrix(label2, edge2)
                v1 = der * (-s.polygon(label2).edge(edge2))
                label = label2
                vertex = (edge2 + 1) % len(s.polygon(label2).vertices())
                v2 = der * (s.polygon(label2).edge(vertex))
            assert count < rotate_limit, "Reached limit!"
            if code:
                return (
                    self.surface().tangent_vector(
                        label, s.polygon(label).vertex(vertex), w
                    ),
                    codes,
                )
            else:
                return self.surface().tangent_vector(
                    label, s.polygon(label).vertex(vertex), w
                )
        else:
            raise NotImplementedError(
                "Rotating tangent vectors is only implemented when at a singularity"
            )

    def counterclockwise_to(self, w, code=False):
        r"""
        Return the new tangent vector obtained by rotating this one in the counterclockwise
        direction until the vector is parallel to w, and scaling so that the length matches
        that of w.

        Note that we always do some rotation so that if w is parallel to this vector, then a
        360 degree rotation is performed.

        The vector w must be nonzero.

        On an infinite surface, this is potentially an infinite calculation
        so we impose a limit (representing the maximal number of polygons
        that must be rotated through). This is the variable rotate_limit
        in this package.

        If code is True, we compute the sequences of numbers associated to edges
        crossed as a list. We return a pair consisting of the newly computing
        tangent vector an this code. This is currently only implemented when
        based at a singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s=translation_surfaces.regular_octagon()
            sage: v=s.tangent_vector(0, (0, 0), (1, 1))
            sage: v.counterclockwise_to((-1, -1))
            doctest:warning
            ...
            UserWarning: counterclockwise_to() has been deprecated and will be removed in a future version of sage-flatsurf; use rotate(w, preserve_scaling=False) instead
            (-1, -1) at vertex 3 of polygon 0
            sage: v.counterclockwise_to((1, 1))
            (1, 1) at vertex 1 of polygon 0
            sage: v.counterclockwise_to((1, 1), code=True)
            ((1, 1) at vertex 1 of polygon 0, [7, 2, 5])

        """
        import warnings
        warnings.warn("counterclockwise_to() has been deprecated and will be removed in a future version of sage-flatsurf; use rotate(w, preserve_scaling=False) instead")

        if not w:
            raise ValueError("w must be non-zero")

        if not self.surface().is_translation_surface():
            raise NotImplementedError("counterclockwise_to() only implemented on translation surfaces")

        if self.is_based_at_singularity():
            s = self.surface()
            v1 = self.vector()
            label = self.polygon_label()
            vertex = self.vertex()
            previous_vertex = (vertex - 1 + len(s.polygon(label).vertices())) % len(
                s.polygon(label).vertices()
            )
            v2 = -s.polygon(label).edge(previous_vertex)
            from sage.matrix.constructor import Matrix

            der = Matrix(s.base_ring(), [[1, 0], [0, 1]])
            if code:
                codes = []

            from flatsurf.geometry.euclidean import ccw

            if not (ccw(v1, w) > 0 and ccw(w, v2) > 0):
                for count in range(rotate_limit):
                    label2, edge2 = s.opposite_edge(label, previous_vertex)
                    if code:
                        codes.append(previous_vertex)
                    der = der * s.edge_matrix(label2, edge2)
                    label = label2
                    vertex = edge2
                    previous_vertex = (
                        vertex - 1 + len(s.polygon(label).vertices())
                    ) % len(s.polygon(label).vertices())
                    v1 = der * (s.polygon(label).edge(vertex))
                    v2 = der * (-s.polygon(label).edge(previous_vertex))
                    if ccw(v1, w) >= 0 and ccw(w, v2) > 0:
                        # We've found it!
                        break
                assert count < rotate_limit, "Reached limit!"
            if code:
                return (
                    self.surface().tangent_vector(
                        label, s.polygon(label).vertex(vertex), w
                    ),
                    codes,
                )
            else:
                return self.surface().tangent_vector(
                    label, s.polygon(label).vertex(vertex), w
                )
        else:
            raise NotImplementedError(
                "Rotating tangent vectors is only implemented when at a singularity"
            )

    def plot(self, **kwargs):
        r"""
        Return a plot of this tangent vector.

        EXAMPLES:

        .. jupyter-execute::

            sage: import flatsurf
            sage: S = flatsurf.translation_surfaces.veech_double_n_gon(5)
            sage: v = S.tangent_vector(0, (1/8, 1/4), (1/2, 1/4))
            sage: S.plot() + v.plot()
            Graphics object consisting of 22 graphics primitives

        Any keyword arguments are passed on to the underlying plot method from SageMath:

        .. jupyter-execute::

            sage: S.plot() + v.plot(color="red")
            Graphics object consisting of 22 graphics primitives

        """
        return self.vector().plot(
            **{"start": self.point(), "width": 1, "arrowsize": 2, **kwargs}
        )

    def translate(self, holonomy, reverse=False):
        r"""
        Return the tangent vector obtained by translating the base point of
        this tangent vector by ``holonomy``.

        INPUT:

        - ``holonomy`` -- a vector that is parallel to this tangent vector

        - ``reverse`` -- a boolean (default: ``False``); whether to return the
          tangent vector obtained by translating (if ``False``), or the tangent
          vector that points back to the start point (if ``True``)

        EXAMPLES::

            sage: import flatsurf
            sage: S = flatsurf.translation_surfaces.veech_double_n_gon(5)
            sage: v = S.tangent_vector(0, (1/8, 1/4), (1/2, 1/4))

            sage: v.translate((1/2, 1/4))
            (1/2, 1/4) at (5/8, 1/2) of polygon 0

        A case that crosses over a polygon boundary at an edge::

            sage: v.translate((2, 1))
            (1/2, 1/4) at (-1/2*a^2 + 21/8, -1/2*a + 5/4) of polygon 0

            sage: v.translate((2, 1), reverse=True)
            (-1/2, -1/4) at (-1/2*a^2 + 21/8, -1/2*a + 5/4) of polygon 0

        A case that hits a singularity::

            sage: v = S.tangent_vector(0, (0, 0), (1, 0))
            sage: v.translate((1, 0))
            Traceback (most recent call last):
            ...
            ValueError: vector must not point out of the polygon at singularity
            sage: v.translate((1, 0), reverse=True)
            (-1, 0) at vertex 0 of polygon 1

        """
        holonomy = (self.surface().base_ring()**2)(holonomy)

        if not holonomy:
            if reverse:
                return self.invert()
            return self

        from flatsurf.geometry.euclidean import is_parallel
        if not is_parallel(holonomy, self.vector()):
            raise ValueError("translation must be parallel to the direction of this tangent vector")

        label = self.polygon_label()
        point = self.point()
        direction = self.vector()

        while True:
            point, holonomy, _ = self.surface().polygon(label).flow(point, holonomy)

            if not holonomy:
                if reverse:
                    return self._bundle(label, point, -direction)
                return self._bundle(label, point, direction)

            if self.surface()(label, point).angle() != 1:
                raise ValueError("translation must not cross over a singularity")

            # Turn the remaining holonomy into the next polygon
            holonomy = self._bundle(label, point, -holonomy)
            holonomy = holonomy.rotate(-holonomy.vector())
            holonomy = holonomy.vector()

            # Turn the tangent vector into the next polygon.
            exit = self._bundle(label, point, -direction)
            exit = exit.rotate(-exit.vector())
            label = exit.polygon_label()
            point = exit.point()
            direction = exit.vector()


# TODO: Rename to TangentBundle
class SimilaritySurfaceTangentBundle(Parent):
    r"""
    Construct the tangent bundle of a given similarity surface.

    Needs work: We should check for coercion from the base_ring of the surface

    TESTS::

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.square_torus()
        sage: TS = S.tangent_bundle()
        sage: TestSuite(TS).run()

    """
    Element = TangentVector

    def __init__(self, surface, ring, category=None):
        assert not surface.is_mutable()

        self._surface = surface

        if category is None:
            from flatsurf.geometry.categories.tangent_bundles import TangentBundles
            category = TangentBundles(surface.category())

        super().__init__(ring, category=category)

    def _an_element_(self):
        return self.edge(next(iter(self.surface().labels())), 0)

    def _element_constructor_(self, label, point, vector):
        r"""
        Construct a tangent vector from a polygon label, a point in the polygon and a vector. The point and the vector should have coordinates
        in the base field."""
        vector = (self.base_ring()**2)(vector)

        if not vector:
            raise ValueError("tangent vector must be non-zero")

        if label is None:
            label, point = self._call_unpack_point(point, vector)

        point = (self.base_ring()**2)(point)

        label, point, vector = self._normalize(label, point, vector)

        return self.element_class(
            self, label, point, vector
        )

    def _call_unpack_point(self, point, vector):
        r"""
        Check that ``point`` describes the base point of a tangent vector
        independent of its representation in one of the polygons of the
        surface and return one such representative.

        This is a helper method for :meth:`__call__`.
        """
        candidates = { (normalized[0], normalized[1]) for (label, coordinates) in point.representatives() if (normalized := self._normalize(label, coordinates, vector)) is not None }

        if not candidates:
            raise ValueError("vector does not describe a tangent vector at this point")
        if len(candidates) > 1:
            raise ValueError("vector does not describe a unique tangent vector at this point")

        return next(iter(candidates))

    def _normalize(self, label, point, vector):
        from flatsurf.geometry.euclidean import ccw, is_anti_parallel

        pos = self.surface().polygon(label).get_point_position(point)

        if pos.is_outside():
            raise ValueError("point must be a point of the polygon")

        if pos.is_in_interior():
            pass
        elif pos.is_in_edge_interior():
            edge = pos.get_edge()
            edge_vector = self.surface().polygon(label).edge(edge)

            if ccw(edge_vector, vector) < 0 or is_anti_parallel(edge_vector, vector):
                # Move point and vector to opposite edge.
                opposite_label, opposite_edge = self.surface().opposite_edge(label, edge)
                similarity = self.surface().edge_transformation(label, edge)

                label = opposite_label
                point = similarity(point)
                vector = similarity.derivative() * vector
        else:
            assert pos.is_vertex()

            from flatsurf.geometry.euclidean import is_between, is_parallel

            edges = self.surface()(label, pos.get_vertex()).edges_ccw((label, pos.get_vertex()))

            # Find the polygon for which the vector points into the polygon.
            for (((sector_start_label, sector_start_edge), _), ((sector_end_label, sector_end_edge), _)) in zip(edges[::2], edges[1::2]):
                assert sector_start_label == sector_end_label

                label = sector_start_label
                polygon = self.surface().polygon(label)
                point = polygon.vertex(sector_start_edge)

                sector_start_vector = polygon.edge(sector_start_edge)
                sector_end_vector = -polygon.edge(sector_end_edge);

                if is_parallel(sector_start_vector, vector) or is_between(sector_start_vector, sector_end_vector, vector):
                    break

                if not is_parallel(sector_end_vector, vector):
                    if self.surface()(label, point).angle() != 1:
                        raise ValueError("vector must not point out of the polygon at singularity")

                vector = self.surface().edge_matrix(sector_end_label, sector_end_edge) * vector
            else:
                assert False, "could not normalize tangent vector at vertex"

        return label, point, vector

    def __repr__(self):
        return "Tangent bundle of {!r} defined over {!r}".format(
            self._surface, self.base_ring()
        )

    def base_ring(self):
        return self._base

    field = base_ring

    def vector_space(self):
        r"""
        Return the vector space over the field of the bundle.
        """
        # TODO: Deprecate
        from sage.modules.free_module import VectorSpace

        return VectorSpace(self.base_ring(), 2)

    def surface(self):
        r"""Return the surface this bundle is over."""
        return self._surface

    def __eq__(self, other):
        if not isinstance(other, SimilaritySurfaceTangentBundle):
            return False
        return self.surface() == other.surface()

    def __hash__(self):
        return hash(self.surface())

    def edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels
        counterclockwise around the boundary of the polygon along the edge with the provided index.
        The length of the vector matches the length of the indexed edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: tb = s.tangent_bundle()
            sage: s.polygon(0)
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
            sage: tb.edge(0,0)
            (2, -2) at vertex 0 of polygon 0
        """
        polygon = self.surface().polygon(polygon_label)
        point = polygon.vertex(edge_index)
        vector = polygon.edge(edge_index)
        return self(polygon_label, point, vector)

    def clockwise_edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels
        *clockwise* around the boundary of the polygon along the edge with the provided index.
        The length of the vector matches the length of the indexed edge.
        Note that the point will be based in the polygon opposite the provided edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: tb = s.tangent_bundle()
            sage: s.polygon(0)
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
            sage: s.polygon(1)
            Polygon(vertices=[(0, 0), (2, 0), (1, 3)])
            sage: s.opposite_edge(0, 0)
            (1, 1)
            sage: tb.clockwise_edge(0,0)
            (-1, 3) at vertex 1 of polygon 1

        """
        polygon = self.surface().polygon(polygon_label)
        point = polygon.vertex(edge_index + 1)
        vector = -polygon.edge(edge_index)
        return self(polygon_label, point, vector)

# TODO: TangentBundle is a functor from XYZ surfaces to tangent bundles of XYZ surfaces.
# The super cate
