from sage.all import SageObject
from sage.misc.cachefunc import cached_method


class SaddleConnection(SageObject):
    r"""
    Represents a saddle connection on a SimilaritySurface.

    TESTS::

        sage: from flatsurf.geometry.saddle_connection import SaddleConnection
        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.cathedral(1, 2)
        sage: SaddleConnection(S, (1, 8), (0, -1))
        Saddle connection in direction (0, -1) with start data (1, 8) and end data (1, 5)

    """

    def __init__(
        self,
        surface,
        start_data,
        direction,
        end_data=None,
        end_direction=None,
        holonomy=None,
        end_holonomy=None,
        check=True,
        limit=1000,
    ):
        r"""
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
        start_data : a pair
            consisting of the label of the polygon where the saddle connection starts
            and the starting vertex.
        direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction the saddle connection is moving in (in the
            coordinates of the initial polygon).
        end_data : a pair
            consisting of the label of the polygon where the saddle connection terminates
            and the terminating vertex.
        end_direction : 2-dimensional vector with entries in the base_ring of the surface
            representing the direction to move backward from the end point (in the
            coordinates of the terminal polygon). If the surface is a DilationSurface
            or better this will be the negation of the direction vector. If the surface
            is a HalfDilation surface or better, then this will be either the direction
            vector or its negation. In either case the value can be inferred from the
            end_data.
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
        """
        from flatsurf.geometry.categories import SimilaritySurfaces

        if surface not in SimilaritySurfaces():
            raise TypeError

        self._surface = surface

        # Sanitize the direction vector:
        V = self._surface.base_ring().fraction_field() ** 2
        self._direction = V(direction)
        if self._direction == V.zero():
            raise ValueError("Direction must be nonzero.")
        # To canonicalize the direction vector we ensure its endpoint lies in the boundary of the unit square.
        xabs = self._direction[0].abs()
        yabs = self._direction[1].abs()
        if xabs > yabs:
            self._direction = self._direction / xabs
        else:
            self._direction = self._direction / yabs

        # Fix end_direction if not standard.
        if end_direction is not None:
            xabs = end_direction[0].abs()
            yabs = end_direction[1].abs()
            if xabs > yabs:
                end_direction = end_direction / xabs
            else:
                end_direction = end_direction / yabs

        self._surfacetart_data = tuple(start_data)

        if end_direction is None:
            from flatsurf.geometry.categories import DilationSurfaces

            # Attempt to infer the end_direction.
            if self._surface in DilationSurfaces().Positive():
                end_direction = -self._direction
            elif self._surface in DilationSurfaces() and end_data is not None:
                p = self._surface.polygon(end_data[0])
                from flatsurf.geometry.euclidean import ccw

                if (
                    ccw(p.edge(end_data[1]), self._direction) >= 0
                    and ccw(
                        p.edge(
                            (len(p.vertices()) + end_data[1] - 1) % len(p.vertices())
                        ),
                        self._direction,
                    )
                    > 0
                ):
                    end_direction = self._direction
                else:
                    end_direction = -self._direction

        if end_holonomy is None and holonomy is not None:
            # Attempt to infer the end_holonomy:
            from flatsurf.geometry.categories import (
                HalfTranslationSurfaces,
                TranslationSurfaces,
            )

            if self._surface in TranslationSurfaces():
                end_holonomy = -holonomy
            if self._surface in HalfTranslationSurfaces():
                if direction == end_direction:
                    end_holonomy = holonomy
                else:
                    end_holonomy = -holonomy

        if (
            end_data is None
            or end_direction is None
            or holonomy is None
            or end_holonomy is None
            or check
        ):
            v = self.start_tangent_vector()
            traj = v.straight_line_trajectory()
            traj.flow(limit)
            if not traj.is_saddle_connection():
                raise ValueError(
                    "Did not obtain saddle connection by flowing forward. Limit="
                    + str(limit)
                )
            tv = traj.terminal_tangent_vector()
            self._end_data = (tv.polygon_label(), tv.vertex())
            if end_data is not None:
                if end_data != self._end_data:
                    raise ValueError(
                        "Provided or inferred end_data="
                        + str(end_data)
                        + " does not match actual end_data="
                        + str(self._end_data)
                    )

            self._end_direction = tv.vector()
            # Canonicalize again.
            xabs = self._end_direction[0].abs()
            yabs = self._end_direction[1].abs()
            if xabs > yabs:
                self._end_direction = self._end_direction / xabs
            else:
                self._end_direction = self._end_direction / yabs
            if end_direction is not None:
                if end_direction != self._end_direction:
                    raise ValueError(
                        "Provided or inferred end_direction="
                        + str(end_direction)
                        + " does not match actual end_direction="
                        + str(self._end_direction)
                    )

            if traj.segments()[0].is_edge():
                # Special case (The below method causes error if the trajectory is just an edge.)
                self._holonomy = self._surface.polygon(start_data[0]).edge(
                    start_data[1]
                )
                self._end_holonomy = self._surface.polygon(self._end_data[0]).edge(
                    self._end_data[1]
                )
            else:
                from .similarity import SimilarityGroup

                sim = SimilarityGroup(self._surface.base_ring()).one()
                itersegs = iter(traj.segments())
                next(itersegs)
                for seg in itersegs:
                    sim = sim * self._surface.edge_transformation(
                        seg.start().polygon_label(), seg.start().position().get_edge()
                    )
                self._holonomy = (
                    sim(traj.segments()[-1].end().point())
                    - traj.initial_tangent_vector().point()
                )
                self._end_holonomy = -((~sim.derivative()) * self._holonomy)

            if holonomy is not None:
                if holonomy != self._holonomy:
                    print("Combinatorial length: " + str(traj.combinatorial_length()))
                    print("Start: " + str(traj.initial_tangent_vector().point()))
                    print("End: " + str(traj.terminal_tangent_vector().point()))
                    print("Start data:" + str(start_data))
                    print("End data:" + str(end_data))
                    raise ValueError(
                        "Provided holonomy "
                        + str(holonomy)
                        + " does not match computed holonomy of "
                        + str(self._holonomy)
                    )
            if end_holonomy is not None:
                if end_holonomy != self._end_holonomy:
                    raise ValueError(
                        "Provided or inferred end_holonomy "
                        + str(end_holonomy)
                        + " does not match computed end_holonomy of "
                        + str(self._end_holonomy)
                    )
        else:
            self._end_data = tuple(end_data)
            self._end_direction = end_direction
            self._holonomy = holonomy
            self._end_holonomy = end_holonomy

        # Make vectors immutable
        self._direction.set_immutable()
        self._end_direction.set_immutable()
        self._holonomy.set_immutable()
        self._end_holonomy.set_immutable()

    def surface(self):
        return self._surface

    def direction(self):
        r"""
        Returns a vector parallel to the saddle connection pointing from the start point.

        The will be normalized so that its $l_\infty$ norm is 1.
        """
        return self._direction

    def end_direction(self):
        r"""
        Returns a vector parallel to the saddle connection pointing from the end point.

        The will be normalized so that its `l_\infty` norm is 1.
        """
        return self._end_direction

    def start_data(self):
        r"""
        Return the pair (l, v) representing the label and vertex of the corresponding polygon
        where the saddle connection originates.
        """
        return self._surfacetart_data

    def end_data(self):
        r"""
        Return the pair (l, v) representing the label and vertex of the corresponding polygon
        where the saddle connection terminates.
        """
        return self._end_data

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

    def end_holonomy(self):
        r"""
        Return the holonomy vector of the saddle connection (measured from the end).

        In a SimilaritySurface this notion corresponds to developing the saddle connection into the plane
        using the initial chart coming from the initial polygon.
        """
        return self._end_holonomy

    def start_tangent_vector(self):
        r"""
        Return a tangent vector to the saddle connection based at its start.
        """
        return self._surface.tangent_vector(
            self._surfacetart_data[0],
            self._surface.polygon(self._surfacetart_data[0]).vertex(
                self._surfacetart_data[1]
            ),
            self._direction,
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
        Return a tangent vector to the saddle connection based at its start.
        """
        return self._surface.tangent_vector(
            self._end_data[0],
            self._surface.polygon(self._end_data[0]).vertex(self._end_data[1]),
            self._end_direction,
        )

    def invert(self):
        r"""
        Return this saddle connection but with opposite orientation.
        """
        return SaddleConnection(
            self._surface,
            self._end_data,
            self._end_direction,
            self._surfacetart_data,
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
            sage: connections = S.saddle_connections(13)

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
        if not self._surface == other._surface:
            return False
        if not self._direction == other._direction:
            return False
        if not self._surfacetart_data == other._surfacetart_data:
            return False
        # Initial data should determine the saddle connection:
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return 41 * hash(self._direction) - 97 * hash(self._surfacetart_data)

    def _test_geometry(self, **options):
        # Test that this saddle connection actually exists on the surface.
        SaddleConnection(
            self._surface,
            self._surfacetart_data,
            self._direction,
            self._end_data,
            self._end_direction,
            self._holonomy,
            self._end_holonomy,
            check=True,
        )

    def __repr__(self):
        return "Saddle connection in direction {} with start data {} and end data {}".format(
            self._direction, self._surfacetart_data, self._end_data
        )

    def _test_inverse(self, **options):
        # Test that inverting works properly.
        SaddleConnection(
            self._surface,
            self._end_data,
            self._end_direction,
            self._surfacetart_data,
            self._direction,
            self._end_holonomy,
            self._holonomy,
            check=True,
        )

    def is_closed(self):
        return self.start() == self.end()

    def start(self):
        return self._surface(*self.start_data())

    def end(self):
        return self._surface(*self.end_data())

    def homology(self):
        r"""
        Return a homology class (generated by edges) that is homologous to this saddle connection.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: connection = S.saddle_connections(13)[-1]
            sage: connection.homology()
            3*B[(0, 0)] - B[(0, 1)]

        ::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.cathedral(1, 2)
            sage: connections = [connection for connection in S.saddle_connections(13) if connection.is_closed()]
            sage: connections[-1].homology()
            B[(3, 7)]

        """
        to_pyflatsurf = self._surface.pyflatsurf()

        connection = to_pyflatsurf(self)

        # TODO: We should probably make chains and saddle connections proper objects
        # in sage-flatsurf so that they can be mapped through
        # to_pyflatsurf.section()
        chain = connection._connection.chain()

        chain = {e.positive().id(): chain[e] for e in to_pyflatsurf.codomain()._flat_triangulation.edges()}

        from sage.all import ZZ
        chain = {
            ([label for label in to_pyflatsurf.codomain().labels() if e in label][0], [label for label in to_pyflatsurf.codomain().labels() if e in label][0].index(e)): ZZ(multiplicity) for (e, multiplicity) in chain.items()}

        from flatsurf.geometry.homology import SimplicialHomology
        homology = SimplicialHomology(to_pyflatsurf.codomain())

        chain = sum(multiplicity * homology(e) for (e, multiplicity) in chain.items())

        return to_pyflatsurf.section()(chain)
