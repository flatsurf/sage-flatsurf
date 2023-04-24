r"""
The category of surfaces built by gluing Euclidean polygons.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for surfaces.

EXAMPLES::

    sage: from flatsurf import Surface_dict
    sage: C = Surface_dict(QQ).category()

    sage: from flatsurf.geometry.categories import RealProjectivePolygonalSurfaces
    sage: C.is_subcategory(RealProjectivePolygonalSurfaces())
    True

"""
# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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
# ####################################################################

from flatsurf.geometry.categories.surface_category import SurfaceCategory


class RealProjectivePolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing (Euclidean) polygons (or more
    generally, polygons in two-dimensional real-projective space.)

    EXAMPLES::

        sage: from flatsurf.geometry.categories import RealProjectivePolygonalSurfaces
        sage: RealProjectivePolygonalSurfaces()
        Category of real projective polygonal surfaces

    """

    def super_categories(self):
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces
        return [PolygonalSurfaces()]

    class ParentMethods:
        def graphical_surface(self, *args, **kwargs):
            if "cached" in kwargs:
                kwargs.pop("cached")  # TODO: Warn that this is gone.
            from flatsurf.graphical.surface import GraphicalSurface
            return GraphicalSurface(self, *args, **kwargs)

        def plot(self, **kwargs):
            # TODO: Deprecate all this.
            graphical_surface_keywords = {
                key: kwargs.pop(key)
                for key in [
                    "cached",
                    "adjacencies",
                    "polygon_labels",
                    "edge_labels",
                    "default_position_function",
                ]
                if key in kwargs
            }
            return self.graphical_surface(**graphical_surface_keywords).plot(**kwargs)

        # TODO: What should be the fate of this function?
        def plot_polygon(
            self,
            label,
            graphical_surface=None,
            plot_polygon=True,
            plot_edges=True,
            plot_edge_labels=True,
            edge_labels=None,
            polygon_options={"axes": True},
            edge_options=None,
            edge_label_options=None,
        ):
            r"""
            Returns a plot of the polygon with the provided label.

            Note that this method plots the polygon in its coordinates as opposed to
            graphical coordinates that the :func:``plot`` method uses. This makes it useful
            for visualizing the natural coordinates of the polygon.

            INPUT:

                - ``graphical_surface`` -- (default ``None``) If provided this function pulls graphical options
                  from the graphical surface. If not provided, we use the default graphical surface.

                - ``plot_polygon`` -- (default ``True``) If True, we plot the solid polygon.

                - ``polygon_options`` -- (default ``{"axes":True}``) Options for the rendering of the polygon.
                  These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_polygon`.
                  This should be either None or a dictionary.

                - ``plot_edges`` -- (default ``True``) If True, we plot the edges of the polygon as segments.

                - ``edge_options`` -- (default ``None``) Options for the rendering of the polygon edges.
                  These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_edge`.
                  This should be either None or a dictionary.

                - ``plot_edge_labels`` -- (default ``True``) If True, we plot labels on the edges.

                - ``edge_label_options`` -- (default ``None``) Options for the rendering of the edge labels.
                  These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_edge_label`.
                  This should be either None or a dictionary.

                - ``edge_labels`` -- (default ``None``) If None and plot_edge_labels is True, we write the edge
                  number on each edge. Otherwise edge_labels should be a list of strings of length equal to the
                  number of edges of the polygon. The strings will be printed on each edge.

            EXAMPLES::

                sage: from flatsurf import *
                sage: s = similarity_surfaces.example()
                sage: s.plot()
                ...Graphics object consisting of 13 graphics primitives
                sage: s.plot_polygon(1)
                ...Graphics object consisting of 7 graphics primitives

                sage: labels = []
                sage: p = s.polygon(1)
                sage: for e in range(p.num_edges()): \
                    labels.append(str(p.edge(e)))
                sage: s.plot_polygon(1, polygon_options=None, plot_edges=False, \
                    edge_labels=labels, edge_label_options={"color":"red"})
                ...Graphics object consisting of 4 graphics primitives
            """
            if graphical_surface is None:
                graphical_surface = self.graphical_surface()
            p = self.polygon(label)
            from flatsurf.graphical.polygon import GraphicalPolygon

            gp = GraphicalPolygon(p)

            if plot_polygon:
                if polygon_options is None:
                    o = graphical_surface.polygon_options
                else:
                    o = graphical_surface.polygon_options.copy()
                    o.update(polygon_options)
                plt = gp.plot_polygon(**o)

            if plot_edges:
                if edge_options is None:
                    o = graphical_surface.non_adjacent_edge_options
                else:
                    o = graphical_surface.non_adjacent_edge_options.copy()
                    o.update(edge_options)
                for e in range(p.num_edges()):
                    plt += gp.plot_edge(e, **o)

            if plot_edge_labels:
                if edge_label_options is None:
                    o = graphical_surface.edge_label_options
                else:
                    o = graphical_surface.edge_label_options.copy()
                    o.update(edge_label_options)
                for e in range(p.num_edges()):
                    if edge_labels is None:
                        el = str(e)
                    else:
                        el = edge_labels[e]
                    plt += gp.plot_edge_label(e, el, **o)
            return plt