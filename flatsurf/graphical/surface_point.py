# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix
#                     2013-2019 W. Patrick Hooper
#                     2022-2023 Julian Rüth
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

from sage.plot.point import point2d


class GraphicalSurfacePoint:
    def __init__(self, surface_point, graphical_surface=None):
        r"""
        Create a graphical segment from SurfacePoint.

        If a graphical_surface is provided the point is created on the graphical surface.
        Otherwise, we create it on the default graphical surface.
        """
        if graphical_surface is None:
            self._gs = surface_point.surface().graphical_surface()
        else:
            if surface_point.surface() != graphical_surface.get_surface():
                raise ValueError
            self._gs = graphical_surface
        self._sp = surface_point

    def surface_point(self):
        r"""
        Return the underlying SurfacePoint.
        """
        return self._sp

    def points(self):
        r"""Return the list of points as RDF vectors."""
        point_list = []
        for label in self._sp.labels():
            if self._gs.is_visible(label):
                for coord in self._sp.coordinates(label):
                    point_list.append(
                        self._gs.graphical_polygon(label).transform(coord)
                    )
        return point_list

    def plot(self, **options):
        r"""
        Plot the point (which might involve drawing several dots).

        The options are passed to point2d.

        If no "zorder" option is provided then we set "zorder" to 50.
        """
        if "zorder" not in options:
            options["zorder"] = 50
        return point2d(points=self.points(), **options)
