r"""
The category of surfaces built by gluing Euclidean polygons.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for immutable surfaces.

EXAMPLES::

    sage: from flatsurf import MutableOrientedSimilaritySurface
    sage: C = MutableOrientedSimilaritySurface(QQ).category()

    sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
    sage: C.is_subcategory(EuclideanPolygonalSurfaces())
    True

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

"""

# ####################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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


class EuclideanPolygonalSurfaces(SurfaceCategory):
    r"""
    The category of surfaces built by gluing Euclidean polygons (or more
    generally, polygons in two-dimensional real space).

    EXAMPLES::

        sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
        sage: EuclideanPolygonalSurfaces()
        Category of euclidean polygonal surfaces

    """

    def super_categories(self):
        r"""
        The categories such surfaces are also automatically contained in,
        namely the category of surfaces built from polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygonalSurfaces
            sage: C = EuclideanPolygonalSurfaces()
            sage: C.super_categories()
            [Category of polygonal surfaces]

        """
        from flatsurf.geometry.categories.polygonal_surfaces import PolygonalSurfaces

        return [PolygonalSurfaces()]

    class ParentMethods:
        r"""
        Provides methods available to all surfaces that are built from polygons
        in the real plane.

        If you want to add functionality for such surfaces you most likely
        want to put it here.
        """

        def graphical_surface(self, *args, **kwargs):
            r"""
            Return a graphical representation of this surface.

            This method can be used to further configure or augment a plot
            beyond the possibilities of :meth:`plot`.

            The documentation of sage-flatsurf contains a section of example
            plots or consult the :mod:`flatsurf.graphical.surface` reference for all the
            details.

            EXAMPLES::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.graphical_surface()
                Graphical representation of Translation Surface in H_1(0) built from a square

            """
            if "cached" in kwargs:
                import warnings

                warnings.warn(
                    "The cached keyword has been removed from graphical_surface(). The keyword is ignored in this version of sage-flatsurf and will be dropped completely in a future version of sage-flatsurf. "
                    "The result of graphical_surface() is never cached anymore."
                )

                kwargs.pop("cached")

            from flatsurf.graphical.surface import GraphicalSurface

            return GraphicalSurface(self, *args, **kwargs)

        def plot(self, **kwargs):
            r"""
            Return a plot of this surface.

            The documentation of sage-flatsurf contains a section of example
            plots or consult the :mod:`flatsurf.graphical.surface` reference
            for all the details.

            EXAMPLES:

            .. jupyter-execute::

                sage: from flatsurf import translation_surfaces
                sage: S = translation_surfaces.square_torus()
                sage: S.plot()
                Graphics object consisting of 10 graphics primitives

            """
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

            EXAMPLES:

            .. jupyter-execute::

                sage: from flatsurf import similarity_surfaces
                sage: s = similarity_surfaces.example()
                sage: s.plot()
                ...Graphics object consisting of 13 graphics primitives

            .. jupyter-execute::

                sage: s.plot_polygon(1)
                ...Graphics object consisting of 7 graphics primitives

            .. jupyter-execute::

                sage: labels = []
                sage: p = s.polygon(1)
                sage: for e in range(len(p.vertices())):
                ....:     labels.append(str(p.edge(e)))
                sage: s.plot_polygon(1, polygon_options=None, plot_edges=False, edge_labels=labels, edge_label_options={"color":"red"})
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
                for e in range(len(p.vertices())):
                    plt += gp.plot_edge(e, **o)

            if plot_edge_labels:
                if edge_label_options is None:
                    o = graphical_surface.edge_label_options
                else:
                    o = graphical_surface.edge_label_options.copy()
                    o.update(edge_label_options)
                for e in range(len(p.vertices())):
                    if edge_labels is None:
                        el = str(e)
                    else:
                        el = edge_labels[e]
                    plt += gp.plot_edge_label(e, el, **o)
            return plt
