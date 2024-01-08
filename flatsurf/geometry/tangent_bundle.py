# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2022 W. Patrick Hooper
#                      2016-2022 Vincent Delecroix
#                      2022-2023 Julian Rüth
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

# Limit for clockwise_to and counter_clockwise_to in SimilaritySurfaceTangentVector.
rotate_limit = 100


class SimilaritySurfaceTangentVector:
    r"""
    Tangent vector to a similarity surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces

    Examples on edges in direction of edges::

        sage: s = translation_surfaces.square_torus()
        sage: s.tangent_vector(0, (1/2, 0), (1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1/2, 0) with vector (1, 0)
        sage: s.tangent_vector(0, (1/2, 0), (-1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1/2, 1) with vector (-1, 0)
        sage: s.tangent_vector(0, (1/2, 1), (1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1/2, 0) with vector (1, 0)
        sage: s.tangent_vector(0, (1/2, 1), (-1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1/2, 1) with vector (-1, 0)

        sage: s.tangent_vector(0, (0, 1/2), (0, 1))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 1/2) with vector (0, 1)
        sage: s.tangent_vector(0, (0, 1/2), (0, -1))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 1/2) with vector (0, -1)
        sage: s.tangent_vector(0, (1, 1/2), (0, 1))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 1/2) with vector (0, 1)
        sage: s.tangent_vector(0, (1, 1/2), (0, -1))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 1/2) with vector (0, -1)

    Examples on vertices in direction of edges::

        sage: s = translation_surfaces.square_torus()
        sage: s.tangent_vector(0, (0, 0), (1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (1, 0)
        sage: s.tangent_vector(0, (1, 0), (-1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 1) with vector (-1, 0)
        sage: s.tangent_vector(0, (0, 1), (1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (1, 0)
        sage: s.tangent_vector(0, (1, 1), (-1, 0))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 1) with vector (-1, 0)

        sage: s.tangent_vector(0, (0, 0), (0, 1))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 0) with vector (0, 1)
        sage: s.tangent_vector(0, (0, 1), (0, -1))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 1) with vector (0, -1)
        sage: s.tangent_vector(0, (1, 0), (0, 1))
        SimilaritySurfaceTangentVector in polygon 0 based at (1, 0) with vector (0, 1)
        sage: s.tangent_vector(0, (1, 1), (0, -1))
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 1) with vector (0, -1)

    """

    def __init__(self, tangent_bundle, polygon_label, point, vector):
        from flatsurf.geometry.euclidean import ccw, is_anti_parallel

        self._bundle = tangent_bundle
        p = self.surface().polygon(polygon_label)
        pos = p.get_point_position(point)
        if not vector:
            raise NotImplementedError("vector must be non-zero")
        if pos.is_in_interior():
            self._polygon_label = polygon_label
            self._point = point
            self._vector = vector
            self._position = pos
        elif pos.is_in_edge_interior():
            e = pos.get_edge()
            edge_v = p.edge(e)

            if ccw(edge_v, vector) < 0 or is_anti_parallel(edge_v, vector):
                # Need to move point and vector to opposite edge.
                label2, e2 = self.surface().opposite_edge(polygon_label, e)
                similarity = self.surface().edge_transformation(polygon_label, e)
                point2 = similarity(point)
                vector2 = similarity.derivative() * vector
                self._polygon_label = label2
                self._point = point2
                self._vector = vector2
                self._position = (
                    self.surface().polygon(label2).get_point_position(point2)
                )
            else:
                self._polygon_label = polygon_label
                self._point = point
                self._vector = vector
                self._position = pos
        elif pos.is_vertex():
            v = pos.get_vertex()
            p = self.surface().polygon(polygon_label)
            # subsequent edge:
            edge1 = p.edge(v)
            # prior edge:
            edge0 = p.edge((v - 1) % len(p.vertices()))
            wp1 = ccw(edge1, vector)
            wp0 = ccw(edge0, vector)
            if wp1 < 0 or wp0 < 0:
                raise ValueError(
                    "Singular point with vector pointing away from polygon"
                )
            if wp0 == 0:
                # vector points backward along edge 0
                label2, e2 = self.surface().opposite_edge(
                    polygon_label, (v - 1) % len(p.vertices())
                )
                similarity = self.surface().edge_transformation(
                    polygon_label, (v - 1) % len(p.vertices())
                )
                point2 = similarity(point)
                vector2 = similarity.derivative() * vector
                self._polygon_label = label2
                self._point = point2
                self._vector = vector2
                self._position = (
                    self.surface().polygon(label2).get_point_position(point2)
                )
            else:
                # vector points along edge1 in that directior or points into polygons interior
                self._polygon_label = polygon_label
                self._point = point
                self._vector = vector
                self._position = pos
        else:
            raise ValueError("Provided point lies outside the indexed polygon")

        self._point.set_immutable()
        self._vector.set_immutable()

    def __repr__(self):
        return (
            "SimilaritySurfaceTangentVector in polygon "
            + repr(self._polygon_label)
            + " based at "
            + repr(self._point)
            + " with vector "
            + repr(self._vector)
        )

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.surface() == other.surface()
                and self.polygon_label() == other.polygon_label()
                and self.point() == other.point()
                and self.vector() == other.vector()
            )
        return NotImplemented

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

    def is_based_at_singularity(self):
        r"""
        Return the truth value of the statement 'the base point for this vector is a singularity.'
        """
        return self._position.is_vertex()

    def vertex(self):
        r"""Return the index of the vertex."""
        return self._position.get_vertex()

    def is_in_boundary_of_polygon(self):
        r"""
        Return the truth value of the statement
        'the base point for this vector lies on the boundary of
        one of the polygons making up the surface.'
        """
        return self._position.is_in_boundary()

    def position(self):
        r"""
        Return the PolygonPosition representing the location of
        the basepoint of the vector in the polygon that contains it.
        """
        return self._position

    def bundle(self):
        r"""Return the tangent bundle containing this vector."""
        return self._bundle

    def polygon_label(self):
        return self._polygon_label

    def polygon(self):
        return self.surface().polygon(self.polygon_label())

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

    def invert(self):
        r"""
        Returns the negation of this tangent vector.
        Raises a ValueError if the vector is based at a singularity.'
        """
        if self.is_based_at_singularity():
            raise ValueError("Can't invert tangent vector based at a singularity.")
        return SimilaritySurfaceTangentVector(
            self.bundle(), self.polygon_label(), self.point(), -self.vector()
        )

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
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: s.polygon(0)
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
            sage: s.polygon(1)
            Polygon(vertices=[(0, 0), (2, 0), (1, 3)])
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
            sage: V = tb.surface().base_ring()**2
            sage: v = SimilaritySurfaceTangentVector(tb, 0, V((0,0)), V((3,-1)))
            sage: v
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (3, -1)
            sage: v2 = v.forward_to_polygon_boundary()
            sage: v2
            SimilaritySurfaceTangentVector in polygon 0 based at (2, -2/3) with vector (-3, 1)
            sage: v2.invert()
            SimilaritySurfaceTangentVector in polygon 1 based at (2/3, 2) with vector (4, -3)
        """
        p = self.polygon()
        point2, pos2 = p.flow_to_exit(self.point(), self.vector())
        # diff=point2-point
        new_vector = SimilaritySurfaceTangentVector(
            self.bundle(), self.polygon_label(), point2, -self.vector()
        )
        return new_vector

    def straight_line_trajectory(self):
        r"""
        Return the straight line trajectory associated to this vector.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: s = translation_surfaces.square_torus()
            sage: v = s.tangent_vector(0, (0,0), (1,1))
            sage: v.straight_line_trajectory()
            Straight line trajectory made of 1 segments from (0, 0) in polygon 0 to (1, 1) in polygon 0
            sage: l = v.straight_line_trajectory()
            sage: l
            Straight line trajectory made of 1 segments from (0, 0) in polygon 0 to (1, 1) in polygon 0
            sage: l.is_saddle_connection()
            True

            sage: v = s.tangent_vector(0, (0,0), (1,1+AA(5).sqrt()), ring=AA)
            sage: l = v.straight_line_trajectory()
            sage: l.flow(20)
            sage: l.segment(20)
            Segment in polygon 0 starting at (0.9442719099991588?, 0) and ending at (1, 0.1803398874989485?)
        """
        from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectory

        return StraightLineTrajectory(self)

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
        that must be rotated through.) This is the variable rotate_limit
        in this package.

        If code is True, we compute the sequences of numbers associated to edges
        crossed as a list. We return a pair consisting of the newly computing
        tangent vector an this code. This is currently only implemented when
        based at a singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s=translation_surfaces.regular_octagon()
            sage: v=s.tangent_vector(0,(0,0),(1,1))
            sage: v.clockwise_to((-1,-1))
            SimilaritySurfaceTangentVector in polygon 0 based at (0, a + 1) with vector (-1, -1)
            sage: v.clockwise_to((1,1))
            SimilaritySurfaceTangentVector in polygon 0 based at (-1/2*a, 1/2*a) with vector (1, 1)
            sage: v.clockwise_to((1,1), code=True)
            (SimilaritySurfaceTangentVector in polygon 0 based at (-1/2*a, 1/2*a) with vector (1, 1), [0, 5, 2])
        """
        if not w:
            raise ValueError("w must be non-zero")

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
                "Rotating tangent vectors is only implemnted when at a singularity"
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
        that must be rotated through.) This is the variable rotate_limit
        in this package.

        If code is True, we compute the sequences of numbers associated to edges
        crossed as a list. We return a pair consisting of the newly computing
        tangent vector an this code. This is currently only implemented when
        based at a singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s=translation_surfaces.regular_octagon()
            sage: v=s.tangent_vector(0,(0,0),(1,1))
            sage: v.counterclockwise_to((-1,-1))
            SimilaritySurfaceTangentVector in polygon 0 based at (1/2*a + 1, 1/2*a + 1) with vector (-1, -1)
            sage: v.counterclockwise_to((1,1))
            SimilaritySurfaceTangentVector in polygon 0 based at (1, 0) with vector (1, 1)
            sage: v.counterclockwise_to((1,1), code=True)
            (SimilaritySurfaceTangentVector in polygon 0 based at (1, 0) with vector (1, 1), [7, 2, 5])
        """
        if not w:
            raise ValueError("w must be non-zero")

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
                "Rotating tangent vectors is only implemnted when at a singularity"
            )

    def plot(self, **kwargs):
        r"""
        Return a plot of this tangent vector.

        EXAMPLES::

            sage: import flatsurf
            sage: S = flatsurf.translation_surfaces.veech_double_n_gon(5)
            sage: v = S.tangent_vector(0, (1/8, 1/4), (1/2, 1/4))
            sage: S.plot() + v.plot()
            Graphics object consisting of 22 graphics primitives

        Any keyword arguments are passed on to the underlying plot method from SageMath::

            sage: S.plot() + v.plot(color="red")
            Graphics object consisting of 22 graphics primitives

        """
        return self.vector().plot(
            **{"start": self.point(), "width": 1, "arrowsize": 2, **kwargs}
        )


class SimilaritySurfaceTangentBundle:
    r"""
    Construct the tangent bundle of a given similarity surface.

    Needs work: We should check for coercion from the base_ring of the surface
    """

    def __init__(self, similarity_surface, ring=None):
        self._s = similarity_surface
        if ring is None:
            self._base_ring = self._s.base_ring()
        else:
            self._base_ring = ring
        from sage.modules.free_module import VectorSpace

        self._V = VectorSpace(self._base_ring, 2)

    def __call__(self, polygon_label, point, vector):
        r"""
        Construct a tangent vector from a polygon label, a point in the polygon and a vector. The point and the vector should have coordinates
        in the base field."""
        return SimilaritySurfaceTangentVector(
            self, polygon_label, self._V(point), self._V(vector)
        )

    def __repr__(self):
        return "Tangent bundle of {!r} defined over {!r}".format(
            self._s, self._base_ring
        )

    def base_ring(self):
        return self._base_ring

    field = base_ring

    def vector_space(self):
        r"""
        Return the vector space over the field of the bundle.
        """
        return self._V

    def surface(self):
        r"""Return the surface this bundle is over."""
        return self._s

    def edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels
        counterclockwise around the boundary of the polygon along the edge with the provided index.
        The length of the vector matches the length of the indexed edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: s.polygon(0)
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
            sage: tb.edge(0,0)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, -2)
        """
        polygon = self.surface().polygon(polygon_label)
        point = polygon.vertex(edge_index)
        vector = polygon.edge(edge_index)
        return SimilaritySurfaceTangentVector(self, polygon_label, point, vector)

    def clockwise_edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels
        *clockwise* around the boundary of the polygon along the edge with the provided index.
        The length of the vector matches the length of the indexed edge.
        Note that the point will be based in the polygon opposite the provided edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: s.polygon(0)
            Polygon(vertices=[(0, 0), (2, -2), (2, 0)])
            sage: s.polygon(1)
            Polygon(vertices=[(0, 0), (2, 0), (1, 3)])
            sage: s.opposite_edge(0, 0)
            (1, 1)
            sage: tb.clockwise_edge(0,0)
            SimilaritySurfaceTangentVector in polygon 1 based at (2, 0) with vector (-1, 3)
        """
        polygon = self.surface().polygon(polygon_label)
        point = polygon.vertex(edge_index + 1)
        vector = -polygon.edge(edge_index)
        return SimilaritySurfaceTangentVector(self, polygon_label, point, vector)
