r"""
Paths on surfaces

A path on a surface is a sequence of compatible
:class:`flatsurf.geometry.segment.Segment straight line segments`, i.e., such
that the end point of one segment is the start point of the next segment.

EXAMPLES::

    sage: from flatsurf import translation_surfaces, Path

    sage: S = translation_surfaces.square_torus()
    sage: p = Path(S, [S.segment(S(0, 0), (1, 0)), S.segment(S(0, 1), (0, 1)), S.segment(S(0, 2), (-1, 0)), S.segment(S(0, 3), (0, -1))])
    sage: p
    [Vertex 0 of polygon 0 to +(1, 0), Vertex 1 of polygon 0 to +(0, 1), Vertex 2 of polygon 0 to +(-1, 0), Vertex 3 of polygon 0 to +(0, -1)]

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


class Path:
    r"""
    A sequence of compatible segments in a surface.
    """
    def __init__(self, surface, segments, check=True):
        # TODO: Would be good to cast segments into the space of segments.
        self._surface = surface

        if check:
            for s, t in zip(segments, segments[1:]):
                if surface(*s.end().base_point()) != surface(*t.start().base_point()):
                    raise ValueError("end points of segments are not compatible to form a path")

            self._segments = segments

    def __repr__(self):
        r"""
        Return a printable representation of this path.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, Path

            sage: S = translation_surfaces.square_torus()
            sage: s = S.segment(S(0, 0), (1, 0))
            sage: Path(S, [s, s])
            [Vertex 0 of polygon 0 to +(1, 0), Vertex 0 of polygon 0 to +(1, 0)]

        """
        return repr(self._segments)
