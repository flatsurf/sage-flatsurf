from sage.all import SageObject
from sage.misc.cachefunc import cached_method


# TODO: SaddleConnection should be an element in the space of SaddleConnections or the space of Paths rather?

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
        end_data=None
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
            warnings.warn("start_data has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use start instead")
            start = start_data
            del start_data

        if start is None:
            raise ValueError("start must be specified to create a SaddleConnection")

        if end_data is not None:
            import warnings
            warnings.warn("end_data has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use end instead")
            end = end_data
            del end_data

        if direction is not None:
            import warnings
            if holonomy is None:
                warnings.warn("direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future of sage-flatsurf; if you want to create a SaddleConnection without specifying the holonomy, use SaddleConnection.from_vertex() instead.")
                c = SaddleConnection.from_vertex(surface, *start, direction, limit=limit)
                start = c._start
                holonomy = c._holonomy
                end = c._end
                end_holonomy = c._end_holonomy
            else:
                warnings.warn("direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future of sage-flatsurf; there is no need to pass this argument anymore when the holonomy is specified.")
            del direction

        if holonomy is None:
            raise ValueError("holonomy must be specified to create a SaddleConnection")

        if end is None or end_holonomy is None:
            import warnings
            warnings.warn("end and end_holonomy must be provided as a keyword argument for SaddleConnection() in future versions of sage-flatsurf; use SaddleConnection.from_vertex() instead to create a SaddleConnection without specifying these.")
            c = SaddleConnection.from_vertex(surface, *start, holonomy, limit=limit)
            start = c._start
            holonomy = c._holonomy
            end = c._end
            end_holonomy = c._end_holonomy

        if end_direction is not None:
            import warnings
            warnings.warn("end_direction has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; the end direction is deduced from the end holonomy automatically")
            del end_direction

        if limit is not None:
            import warnings
            warnings.warn("limit has been deprecated as a keyword argument for SaddleConnection() and will be removed in a future version of sage-flatsurf; use SaddleConnection.from_vertex() to search with a limit instead")
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
                return opposite_edge, self._surface.edge_transformation(label, previous_edge).derivative() * holonomy

        return start, holonomy

    def __neg__(self):
        return SaddleConnection(
            surface=self._surface,
            start=self._end,
            end=self._start,
            holonomy=self._end_holonomy,
            end_holonomy=self._holonomy,
            check=False)

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

        tangent_vector = surface.tangent_vector(label, surface.polygon(label).vertex(vertex), direction.vector())
        trajectory = tangent_vector.straight_line_trajectory()

        if limit is None:
            from sage.all import infinity
            limit = infinity

        trajectory.flow(steps=limit)

        if not trajectory.is_saddle_connection():
            raise ValueError("no saddle connection in this direction within the specified limit")

        end_tangent_vector = trajectory.terminal_tangent_vector()

        assert not trajectory.segments()[0].is_edge() or len(trajectory.segments()) == 1, "when the saddle connection is an edge it must not consist of more than that edge"

        if trajectory.segments()[0].is_edge():
            # When the saddle connection is just an edge, the similarity
            # computation below can be wrong when that edge is glued to the
            # same polygon. Namely, the similarity is missing the final factor
            # that comes from that gluing. E.g., in the Cathedral test case
            # above, the vertical saddle connection connecting the vertices of
            # polygon 1 at (1, 3) and (2, 1) by going in direction (0, -1) is
            # misinterpreted as the connection going in direction (1, -2).
            return SaddleConnection.from_half_edge(surface, label, trajectory.segments()[0].edge())

        from flatsurf.geometry.similarity import SimilarityGroup
        one = SimilarityGroup(surface.base_ring()).one()
        segments = list(trajectory.segments())[1:]

        from sage.all import prod
        similarity = prod([
            surface.edge_transformation(segment.start().polygon_label(), segment.start().position().get_edge()) for segment in segments], one)

        holonomy = similarity(trajectory.segments()[-1].end().point()) - trajectory.initial_tangent_vector().point()
        end_holonomy = (~similarity.derivative()) * holonomy

        return SaddleConnection(
            surface=surface,
            start=(label, vertex),
            holonomy=holonomy,
            end=(end_tangent_vector.polygon_label(), end_tangent_vector.vertex()),
            end_holonomy=end_holonomy)

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
        warnings.warn("start_data() has been deprecated and will be removed from a future version of sage-flatsurf; use start() instead.")

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
        warnings.warn("end_data() has been deprecated and will be removed from a future version of sage-flatsurf; use end() instead.")

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
        return holonomy[0]**2 + holonomy[1]**2

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
            self._surface.polygon(self._start[0]).vertex(
                self._start[1]
            ),
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

    def homology(self):
        r"""
        Return a homology class (generated by edges) that is homologous to this saddle connection.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: connections = list(S.saddle_connections(13))
            sage: connections[-1].homology()
            -2*B[(0, 0)] - 3*B[(0, 1)]

        ::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.cathedral(1, 2)
            sage: connections = [connection for connection in S.saddle_connections(13) if connection.is_closed()]
            sage: connections[-1].homology()
            -B[(1, 1)] - 2*B[(1, 2)] - B[(1, 6)] - B[(3, 1)] + B[(3, 7)]

        """
        to_pyflatsurf = self._surface.pyflatsurf()

        connection = to_pyflatsurf(self)

        # TODO: We should probably make pyflatsurf chains and saddle
        # connections proper objects in sage-flatsurf so that they can be
        # mapped through to_pyflatsurf.section()
        chain = connection._connection.chain()

        chain = {e.positive().id(): chain[e] for e in to_pyflatsurf.codomain()._flat_triangulation.edges()}

        from sage.all import ZZ
        chain = {
            ([label for label in to_pyflatsurf.codomain().labels() if e in label][0], [label for label in to_pyflatsurf.codomain().labels() if e in label][0].index(e)): ZZ(multiplicity) for (e, multiplicity) in chain.items()}

        from flatsurf.geometry.homology import SimplicialHomology
        homology = SimplicialHomology(to_pyflatsurf.codomain())

        chain = sum(multiplicity * homology(e) for (e, multiplicity) in chain.items())

        return to_pyflatsurf.section()(chain)
