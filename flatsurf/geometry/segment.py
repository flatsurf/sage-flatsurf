r"""
Straight line segments on surfaces

EXAMPLES:

A segment can be created from a tangent vector and a holonomy vector::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.square_torus()

    sage: S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1))
    Vertex 0 of polygon 0 to +(1, 1)

Optionally, we can provide the tangent vector at the end point (since computing
it might be costly)::

    sage: S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1), S.tangent_vector(0, (1, 1), (-1, -1)), check=False)
    Vertex 0 of polygon 0 to +(1, 1)

We can also provide a starting point instead of a tangent vector if that
uniquely singles out the tangent vector given the holonomy::

    sage: S.segment(S(0, 0), (1, 1))
    Vertex 0 of polygon 0 to +(1, 1)

Segments are oriented, so we can construct the opposite of a segment::

    sage: s = S.segment(S(0, 0), (1, 1))
    sage: -s
    Vertex 2 of polygon 0 to +(-1, -1)

Segments with compatible start and end points can be combined into a
:class:`flatsurf.geometry.path.Path`:

    sage: from flatsurf import Path
    sage: Path(S, [s, s, -s])
    [Vertex 0 of polygon 0 to +(1, 1), Vertex 0 of polygon 0 to +(1, 1), Vertex 2 of polygon 0 to +(-1, -1)]


.. NOTE::

    Currently, there is some overlap between this module and
    :mod:`straight_line_trajectory`. Eventually that module should be subsumed
    into this module.

.. NOTE::

    Currently, there is some overlap between a path and
    :class:`flatsurf.geometry.surface_objects.SaddleConnection`. Eventually, a
    saddle connection should probably just be a path thath knows that it is a saddle
    connection. Or the underlying implementation might even be quite different,
    but from a user's perspective, the two should have the same interface.

"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2025 Julian Rüth
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


class Segment:
    r"""
    An oriented straight line segment on a surface.
    """
    def __init__(self, surface, start, holonomy, end=None, check=True):
        holonomy = (surface.base_ring()**2)(holonomy)

        from flatsurf.geometry.tangent_bundle import TangentVector
        if not isinstance(start, TangentVector):
            if start not in surface:
                raise TypeError("start must be a tangent vector or a point in the surface")

            start = surface.tangent_vector(None, start, holonomy)

        if check:
            if start.surface() is not surface:
                raise ValueError("endpoint must be a tangent vector in the given surface")

            from flatsurf.geometry.euclidean import is_parallel
            if not is_parallel(start.vector(), holonomy):
                raise ValueError("tangent vector must be parallel to the holonomy vector")

        self._start = start
        self._holonomy = holonomy
        self._end = None

        if end is not None:
            if not isinstance(end, TangentVector):
                end = surface.tangent_vector(None, end, -holonomy)

            if check:
                if end.surface() is not surface:
                    raise ValueError("endpoint must be a tangent vector in the given surface")

                from flatsurf.geometry.euclidean import is_anti_parallel
                if not is_anti_parallel(end.vector(), holonomy):
                    raise ValueError("tangent vector must be parallel to the holonomy vector")

                if self.end() != end:
                    raise ValueError("end must be the endpoint of the segment")

            self._end = end

    def start(self):
        r"""
        Return the tangent vector at the start of this oriented segment.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: s =  S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1))
            sage: s.start()
            (1, 1) at vertex 0 of polygon 0

        """
        return self._start

    def end(self):
        r"""
        Return the tangent vector at the end of this oriented segment.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: s =  S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1))
            sage: s.end()
            (-1, -1) at vertex 2 of polygon 0

        """
        if self._end is None:
            self._end = self._start.translate(self._holonomy, reverse=True)
        
        return self._end

    def __neg__(self):
        r"""
        Return the reversed version of this oriented segment.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: s =  S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1))
            sage: s
            Vertex 0 of polygon 0 to +(1, 1)
            sage: -s
            Vertex 2 of polygon 0 to +(-1, -1)

        """
        return Segment(surface=self._start.surface(), start=self.end(), holonomy=-self._holonomy, end=self._start, check=False)

    def __repr__(self):
        r"""
        Return a printable representation of this segment.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()

            sage: S.segment(S.tangent_vector(0, (0, 0), (1, 1)), (1, 1))
            Vertex 0 of polygon 0 to +(1, 1)

        """
        surface = self.start().surface()
        representative = self.start().base_point()
        start_point = surface(*representative)

        return f"{start_point._repr_representative(*representative, uppercase=True, shortened=True)} to +{self._holonomy!r}"
