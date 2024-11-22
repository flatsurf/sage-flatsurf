r"""
.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

EXAMPLES:

.. jupyter-execute::

    sage: import flatsurf
    sage: flatsurf.translation_surfaces.veech_2n_gon(4).plot()
    ...Graphics object consisting of 18 graphics primitives

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#                     2022-2023 Julian Rüth <julian.rueth@fsfe.org>
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

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector


class GraphicalSurface:
    r"""
    This class manages the rendering of a SimilaritySurface.

    This class essentially consists of a collection of GraphicalPolygons which
    control how individual polygons are positioned. In addition, this class
    stores options which are passed to the polygons when they are rendered.

    Some setup features set in the constructor and can be set again later via
    :meth:`process_options`.

    The basic tasks of the class are to render the polygons, edges and labels.
    To customize a rendering, it is useful to know something about how this
    class works. (Apologies!)

    There are attributes which control whether or not certain objects are
    rendered, namely:

    - ``will_plot_polygons`` -- Whether to plot polygons which are right-side up.

    - ``will_plot_upside_down_polygons`` -- Whether to plot polygons which are
        upside down. Defaults to False.

    - ``will_plot_polygon_labels`` -- Whether to plot polygon labels.

    - ``will_plot_edges`` -- If this is False then no edges will be plotted.

    - ``will_plot_non_adjacent_edges`` = Whether to plot polygon edges which are
        not adjacent to the edge it is glued to.

    - ``will_plot_adjacent_edges`` -- Whether to plot polygon edges which are
        adjacent to the polygon they are glued to.

    - ``will_plot_self_glued_edges`` -- Whether to plot polygon edges which are
        glued to themselves.

    - ``will_plot_edge_labels`` -- Whether to plot polygon edges labels.

    - ``will_plot_zero_flags`` -- Whether to plot line segments from the
        baricenter to the zero vertex of each polygon. Useful in working out
        vertex and edge labels. Defaults to False.

    The :meth:`plot` method calls some other built in methods: :meth:`plot_polygon`,
    :meth:`plot_polygon_label`, :meth:`plot_edge` and :meth:`plot_edge_label`. These in turn
    call methods in :class:`~.polygon.GraphicalPolygon`.

    - ``polygon_options`` -- Options passed to :meth:`.polygon.GraphicalPolygon.plot_polygon` when
        plotting a polygon right-side up.

    - ``upside_down_polygon_options`` -- Options passed to
        :meth:`.polygon.GraphicalPolygon.plot_polygon` when plotting a polygon upside-down.

    - ``polygon_label_options`` -- Options passed to :meth:`.polygon.GraphicalPolygon.plot_label`
        when plotting a polygon label.

    - ``non_adjacent_edge_options`` -- Options passed to
        :meth:`.polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge which is not
        adjacent to the edge it is glued to.

    - ``adjacent_edge_options`` -- Options passed to
        :meth:`.polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge which is
        adjacent to the edge it is glued to.

    - ``self_glued_edge_options`` -- Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge`
        when plotting a polygon edge which is glued to itself.

    - ``edge_label_options`` -- Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge_label`
        when plotting a edge label.

    - ``zero_flag_options`` -- Options passed to
        :meth:`.polygon.GraphicalPolygon.plot_zero_flag` when plotting a zero_flag.

    INPUT:

    - ``similarity_surface`` -- a similarity surface

    - ``polygon_labels`` -- a boolean (default ``True``) whether the label
      of polygons are displayed

    - ``edge_labels`` -- option to control the display of edge labels. It
      can be one of

        - ``False`` or ``None`` for no labels

        - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
          name of the polygon to which it is glued

        - ``'number'`` -- to put on each side of each edge the number of the
          edge

        - ``'gluings and number'`` -- full information

        - ``'letter'`` -- add matching letters to glued edges in an arbitrary way

    - ``adjacencies`` -- a list of pairs ``(p,e)`` to be used to set
      adjacencies of polygons.

    - ``default_position_function`` -- a function mapping polygon labels to
      similarities describing the position of the corresponding polygon.

    If adjacencies is not defined and the surface is finite, make_all_visible()
    is called to make all polygons visible.

    EXAMPLES:

    .. jupyter-execute::

        sage: from flatsurf import similarity_surfaces
        sage: from flatsurf.graphical.surface import GraphicalSurface

        sage: s = similarity_surfaces.example()
        sage: gs = GraphicalSurface(s)
        sage: gs.polygon_options["color"]="red"
        sage: gs.plot()
        ...Graphics object consisting of 13 graphics primitives

    TESTS:

    Verify that surfaces with boundary can be plotted::

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

        sage: S.plot()
        Graphics object consisting of 9 graphics primitives

    """

    def __init__(
        self,
        surface,
        adjacencies=None,
        polygon_labels=True,
        edge_labels="gluings",
        default_position_function=None,
    ):
        self._ss = surface
        self._default_position_function = default_position_function
        self._polygons = {}

        if not self._ss.is_connected():
            raise NotImplementedError("cannot plot disconnected surfaces yet")

        self._visible = set(self._ss.roots())

        if adjacencies is None:
            if self._ss.is_finite_type():
                self.make_all_visible()
        self._edge_labels = None

        self.will_plot_polygons = True
        r"""
        Whether to plot polygons which are right-side up.
        """

        self.polygon_options = {"color": "lightgray"}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_polygon` when plotting a polygon right-side up."""

        self.will_plot_upside_down_polygons = False
        r"""
        Whether to plot polygons which are upside down
        """

        self.upside_down_polygon_options = {"color": "lightgray", "zorder": -1}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_polygon` when plotting a polygon upside-down."""

        self.will_plot_polygon_labels = True
        r"""
        Whether to plot polygon labels.
        """

        self.polygon_label_options = {
            "color": "black",
            "vertical_alignment": "center",
            "horizontal_alignment": "center",
        }
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_label` when plotting a polygon label."""

        self.will_plot_edges = True
        r"""
        If this is False then no edges will be plotted.
        """

        self.will_plot_non_adjacent_edges = True
        r"""
        Whether to plot polygon edges which are not adjacent to the edge it is glued to.
        """

        self.non_adjacent_edge_options = {"color": "blue"}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge
        which is not adjacent to the edge it is glued to."""

        self.will_plot_adjacent_edges = True
        r"""
        Whether to plot polygon edges which are adjacent to the polygon they are glued to.
        """

        self.adjacent_edge_options = {"color": "blue", "linestyle": ":"}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge`
        when plotting a polygon edge which is adjacent to the edge it is glued to."""

        self.will_plot_self_glued_edges = True
        r"""
        Whether to plot polygon edges which are glued to themselves.
        """

        self.self_glued_edge_options = {"color": "red"}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge
        which is glued to itself."""

        self.will_plot_edge_labels = True
        r"""
        Whether to plot polygon edges which are not adjacent to the polygon they are glued to.
        """

        self.edge_label_options = {"color": "blue"}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_edge_label` when plotting a polygon label."""

        self.will_plot_zero_flags = False
        r"""
        Whether to plot line segments from the baricenter to the zero vertex of each polygon.
        """

        self.zero_flag_options = {"color": "green", "thickness": 0.5}
        r"""Options passed to :meth:`.polygon.GraphicalPolygon.plot_zero_flag` when plotting a zero_flag."""

        self.process_options(
            adjacencies=adjacencies,
            polygon_labels=polygon_labels,
            edge_labels=edge_labels,
        )

    def process_options(
        self,
        adjacencies=None,
        polygon_labels=None,
        edge_labels=None,
        default_position_function=None,
    ):
        r"""
        Process the options listed as if the graphical_surface was first
        created.

        INPUT:

        Consult :class:`GraphicalSurface` for the possible arguments.

        TESTS::

            sage: from flatsurf import translation_surfaces

            sage: c = translation_surfaces.chamanara(1/2)
            sage: gs = c.graphical_surface()
            sage: gs.process_options(edge_labels='hey')
            Traceback (most recent call last):
            ...
            ValueError: invalid value for edge_labels (='hey')
        """
        if adjacencies is not None:
            for p, e in adjacencies:
                self.make_adjacent(p, e)

        if polygon_labels is not None:
            if not isinstance(polygon_labels, bool):
                raise ValueError("polygon_labels must be True, False or None.")
            self.will_plot_polygon_labels = polygon_labels

        if edge_labels is not None:
            if edge_labels is True:
                edge_labels = "gluings"

            if edge_labels is False:
                self._edge_labels = None
                self.will_plot_edge_labels = False
            elif edge_labels in ["gluings", "number", "gluings and number", "letter"]:
                self._edge_labels = edge_labels
                self.will_plot_edge_labels = True
            else:
                raise ValueError(
                    "invalid value for edge_labels (={!r})".format(edge_labels)
                )

        if default_position_function is not None:
            self._default_position_function = default_position_function

    def copy(self):
        r"""
        Make a copy of this GraphicalSurface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.octagon_and_squares()
            sage: gs = s.graphical_surface()
            sage: gs.will_plot_zero_flags = True
            sage: gs.graphical_polygon(1).transformation()
            (x, y) |-> (x + 2, y)
            sage: gs.make_adjacent(0,2)
            sage: gs.graphical_polygon(1).transformation()
            (x, y) |-> (x + (a + 4), y + (a + 2))
            sage: gs.polygon_options["color"]="yellow"
            sage: gs2 = gs.copy()
            sage: gs2 == gs
            False
            sage: gs2.will_plot_zero_flags
            True
            sage: gs2.graphical_polygon(1).transformation()
            (x, y) |-> (x + (a + 4), y + (a + 2))
            sage: gs2.polygon_options
            {'color': 'yellow'}
        """
        gs = GraphicalSurface(
            self.get_surface(),
            default_position_function=self._default_position_function,
        )

        # Copy plot options
        gs.will_plot_polygons = self.will_plot_polygons
        gs.polygon_options = dict(self.polygon_options)
        gs.will_plot_upside_down_polygons = self.will_plot_upside_down_polygons
        gs.upside_down_polygon_options = dict(self.upside_down_polygon_options)
        gs.will_plot_polygon_labels = self.will_plot_polygon_labels
        gs.polygon_label_options = dict(self.polygon_label_options)
        gs.will_plot_edges = self.will_plot_edges
        gs.will_plot_non_adjacent_edges = self.will_plot_non_adjacent_edges
        gs.non_adjacent_edge_options = dict(self.non_adjacent_edge_options)
        gs.will_plot_adjacent_edges = self.will_plot_adjacent_edges
        gs.adjacent_edge_options = dict(self.adjacent_edge_options)
        gs.will_plot_self_glued_edges = self.will_plot_self_glued_edges
        gs.self_glued_edge_options = dict(self.self_glued_edge_options)
        gs.will_plot_edge_labels = self.will_plot_edge_labels
        gs.edge_label_options = dict(self.edge_label_options)
        gs.will_plot_zero_flags = self.will_plot_zero_flags
        gs.zero_flag_options = dict(self.zero_flag_options)

        # Copy polygons and visible set.
        gs._polygons = {label: gp.copy() for label, gp in self._polygons.items()}
        gs._visible = set(self._visible)
        gs._edge_labels = self._edge_labels

        return gs

    def __repr__(self):
        return "Graphical representation of {!r}".format(self._ss)

    def visible(self):
        r"""
        Return a copy of the set of visible labels.
        """
        return set(self._visible)

    def is_visible(self, label):
        r"""
        Return whether the polygon with the given label is marked as visible.
        """
        return label in self._visible

    def make_visible(self, label):
        r"""
        Mark the polygon with the given label as visible.
        """
        self._visible.add(label)

    def hide(self, label):
        r"""
        Mark the polygon with the given label as invisible.
        """
        self._visible.remove(label)

    def make_all_visible(self, adjacent=None, limit=None):
        r"""
        Attempt to show all invisible polygons by walking over the surface.

        INPUT:

        - ``adjacent`` -- (default ``None``) whether the newly added polygon are
          set to be adjacent or not. Defaults to true unless a default_position_function was
          provided.

        - ``limit`` -- (default ``None``) maximal number of additional polygons to make visible

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import similarity_surfaces

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface()
            sage: g.make_all_visible()
            sage: g.plot()
            ...Graphics object consisting of 13 graphics primitives

        .. jupyter-execute::

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(adjacencies=[])
            sage: g.make_all_visible(adjacent=False)
            sage: g.plot()
            ...Graphics object consisting of 16 graphics primitives

        """
        if adjacent is None:
            adjacent = self._default_position_function is None
        if limit is None:
            if not self._ss.is_finite_type():
                raise NotImplementedError
            if adjacent:
                for label, poly in zip(self._ss.labels(), self._ss.polygons()):
                    for e in range(len(poly.vertices())):
                        opposite_edge = self._ss.opposite_edge(label, e)
                        if opposite_edge is None:
                            continue
                        l2, _ = opposite_edge
                        if not self.is_visible(l2):
                            self.make_adjacent(label, e)
            else:
                from flatsurf.geometry.similarity import SimilarityGroup

                T = SimilarityGroup(self._ss.base_ring())
                for label in self._ss.labels():
                    if not self.is_visible(label):
                        if self._default_position_function is None:
                            # No reasonable way to display the polygon, so we do this hack:
                            g = self.graphical_polygon(label)
                            poly = self._ss.polygon(label)
                            t = T(
                                (
                                    QQ(self.xmax() - g.xmin() + 1),
                                    QQ(-(g.ymin() + g.ymax()) / ZZ(2)),
                                )
                            )
                            g.set_transformation(t)
                        self.make_visible(label)
        else:
            if limit <= 0:
                raise ValueError
            if adjacent:
                i = 0
                for label, poly in zip(self._ss.labels(), self._ss.polygons()):
                    for e in range(len(poly.vertices())):
                        l2, e2 = self._ss.opposite_edge(label, e)
                        if not self.is_visible(l2):
                            self.make_adjacent(label, e)
                            i = i + 1
                            if i >= limit:
                                return
            else:
                from flatsurf.geometry.similarity import SimilarityGroup

                T = SimilarityGroup(self._ss.base_ring())
                i = 0
                for label in self._ss.labels():
                    if not self.is_visible(label):
                        if self._default_position_function is None:
                            # No reasonable way to display the polygon, so we do this hack:
                            g = self.graphical_polygon(label)
                            poly = self._ss.polygon(label)
                            t = T(
                                (
                                    QQ(self.xmax() - g.xmin() + 1),
                                    QQ(-(g.ymin() + g.ymax()) / ZZ(2)),
                                )
                            )
                            g.set_transformation(t)
                        self.make_visible(label)
                        i = i + 1
                        if i >= limit:
                            return

    def get_surface(self):
        r"""
        Return the underlying similarity surface.
        """
        return self._ss

    def xmin(self):
        r"""
        Return the minimal x-coordinate of a vertex of a visible graphical polygon.

        .. TODO::

            this should be xmin
        """
        return min([self.graphical_polygon(label).xmin() for label in self.visible()])

    def ymin(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.
        """
        return min([self.graphical_polygon(label).ymin() for label in self.visible()])

    def xmax(self):
        r"""
        Return the maximal x-coordinate of a vertex of a visible graphical polygon.
        """
        return max([self.graphical_polygon(label).xmax() for label in self.visible()])

    def ymax(self):
        r"""
        Return the minimal y-coordinate of a vertex of a visible graphical polygon.
        """
        return max([self.graphical_polygon(label).ymax() for label in self.visible()])

    def bounding_box(self):
        r"""
        Return the quadruple (x1,y1,x2,y2) where x1 and y1 are the minimal
        x- and y-coordinates of a visible graphical polygon and x2 and y2 are the
        maximal x-and y- coordinates of a visible graphical polygon.
        """
        return self.xmin(), self.ymin(), self.xmax(), self.ymax()

    def graphical_polygon(self, label):
        r"""
        Return the graphical_polygon with the given label.
        """
        if label in self._polygons:
            return self._polygons[label]
        else:
            t = None
            if self._default_position_function is not None:
                t = self._default_position_function(label)
            from flatsurf.graphical.polygon import GraphicalPolygon

            p = GraphicalPolygon(self._ss.polygon(label), transformation=t)
            self._polygons[label] = p
            return p

    def make_adjacent(self, p, e, reverse=False, visible=True):
        r"""
        Move the polygon across the prescribed edge so that is adjacent.
        The polygon moved is obtained from opposite_edge(p,e).

        If reverse=True then the polygon is moved so that there is a fold
        at the edge.

        If visible is True (by default), we also make the moved polygon visible.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface(adjacencies=[])
            sage: gs.graphical_polygon(0)
            GraphicalPolygon with vertices [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: gs.graphical_polygon(1)
            GraphicalPolygon with vertices [(0.0, 0.0), (2.0, 0.0), (1.0, 3.0)]
            sage: print("Polygon 0, edge 0 is opposite "+str(gs.opposite_edge(0,0)))
            Polygon 0, edge 0 is opposite (1, 1)
            sage: gs.make_adjacent(0,0)
            sage: gs.graphical_polygon(0)
            GraphicalPolygon with vertices [(0.0, 0.0), (2.0, -2.0), (2.0, 0.0)]
            sage: gs.graphical_polygon(1)
            GraphicalPolygon with vertices [(0.4, -2.8), (2.0, -2.0), (0.0, 0.0)]
        """
        pp, ee = self._ss.opposite_edge(p, e)
        if reverse:
            from flatsurf.geometry.similarity import SimilarityGroup

            G = SimilarityGroup(self._ss.base_ring())

            q = self._ss.polygon(p)
            a = q.vertex(e)
            b = q.vertex(e + 1)
            # This is the similarity carrying the origin to a and (1,0) to b:
            g = G(b[0] - a[0], b[1] - a[1], a[0], a[1])

            qq = self._ss.polygon(pp)
            aa = qq.vertex(ee + 1)
            bb = qq.vertex(ee)
            # This is the similarity carrying the origin to aa and (1,0) to bb:
            gg = G(bb[0] - aa[0], bb[1] - aa[1], aa[0], aa[1])

            reflection = G(
                self._ss.base_ring().one(),
                self._ss.base_ring().zero(),
                self._ss.base_ring().zero(),
                self._ss.base_ring().zero(),
                -1,
            )

            # This is the similarity carrying (a,b) to (aa,bb):
            g = gg * reflection * (~g)
        else:
            g = self._ss.edge_transformation(pp, ee)
        h = self.graphical_polygon(p).transformation()
        self.graphical_polygon(pp).set_transformation(h * g)
        if visible:
            self.make_visible(pp)

    def is_adjacent(self, p, e):
        r"""
        Returns the truth value of the statement
        'The polygon opposite edge (p,e) is adjacent to that edge.'

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(adjacencies=[])
            sage: g.is_adjacent(0,0)
            False
            sage: g.is_adjacent(0,1)
            False
            sage: g.make_all_visible(adjacent=True)
            sage: g.is_adjacent(0,0)
            True
            sage: g.is_adjacent(0,1)
            False

        TESTS:

        Verify that this works correctly for boundary edges::

            sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface
            sage: S = MutableOrientedSimilaritySurface(QQ)
            sage: S.add_polygon(Polygon(vertices=[(0,0), (-1, -1), (1,0)]))
            0
            sage: G = S.graphical_surface()
            sage: G.is_adjacent(0, 0)
            False

        """
        opposite_edge = self.opposite_edge(p, e)
        if opposite_edge is None:
            return False
        pp, ee = opposite_edge
        if not self.is_visible(pp):
            return False
        g = self.graphical_polygon(p)
        gg = self.graphical_polygon(pp)
        return g.transformed_vertex(e) == gg.transformed_vertex(
            ee + 1
        ) and g.transformed_vertex(e + 1) == gg.transformed_vertex(ee)

    def to_surface(
        self,
        point,
        v=None,
        label=None,
        ring=None,
        return_all=False,
        singularity_limit=None,
        search_all=False,
        search_limit=None,
    ):
        r"""Converts from graphical coordinates to similarity surface coordinates.

        A point always must be provided. If a vector v is provided then a
        SimilaritySurfaceTangentVector will be returned. If v is not provided, then a
        SurfacePoint is returned.

        INPUT:

        - ``point`` -- Coordinates of a point in graphical coordinates to be
            converted to surface coordinates.

        - ``v`` -- (default ``None``) If provided a tangent vector in graphical
            coordinates based at the provided point.

        - ``label`` -- (default ``None``) If provided, then we only convert
            points and tangent vectors in the corresponding graphical polygon.

        - ``ring`` -- (default ``None``) If provided, then objects returned
            will be defined over the given ring, otherwise we use the base_ring
            of the surface.

        - ``return_all`` -- (default ``False``) By default we return the first
            point or vector we find. However if the graphical polygons overlap,
            then a point or vector might correspond to more than one point or
            vector on the surface. If ``return_all`` is set to ``True`` then we
            return a set of all points we find instead.

        - ``singularity_limit`` -- (default ``None``) This only has an effect
            if returning a singular point (i.e., ``v`` is ``None``) and the
            surface is infinite. In this case, the singularity should be
            returned but it could be infinite. Then singularity_limit controls
            how far we look for the singularity to close. This value is passed
            to the ``point()`` method on the surface.

        - ``search_all`` -- (default ``False``) By default we look just in
            polygons with visible label. If set to ``True``, then we instead
            look in all labels.

        - ``search_limit`` -- (default ``None``) If ``search_all`` is ``True``,
            then we look at the first ``search_limit`` polygons instead of all
            polygons. This must be set to an positive integer if ``search_all``
            is and the surface is infinite.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: gs.to_surface((1,-2))
            Point (1, 1/2) of polygon 1
            sage: gs.to_surface((1,-2), v=(1,0))
            SimilaritySurfaceTangentVector in polygon 1 based at (1, 1/2) with vector (1, -1/2)

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.infinite_staircase()
            sage: gs = s.graphical_surface()
            sage: gs.to_surface((4,4), (1,1), search_all=True, search_limit=20)
            SimilaritySurfaceTangentVector in polygon 8 based at (0, 0) with vector (1, 1)

            sage: s = translation_surfaces.square_torus()
            sage: pc = s.minimal_cover(cover_type="planar")
            sage: gs = pc.graphical_surface()
            sage: gs.to_surface((3,2), search_all=True, search_limit=20)
            Vertex 0 of polygon (0, (x, y) |-> (x + 3, y + 2))
            sage: p = gs.to_surface((sqrt(3),sqrt(2)), ring=AA, search_all=True, search_limit=20)
            doctest:warning
            ...
            UserWarning: the ring parameter is deprecated and will be removed in a future version of sage-flatsurf; define the surface over a larger ring instead so that this points' coordinates live in the base ring
            sage: next(iter(p.coordinates(next(iter(p.labels()))))).parent()
            Vector space of dimension 2 over Algebraic Real Field
            sage: v = gs.to_surface((3/2,3/2),(sqrt(3),sqrt(2)),ring=AA,search_all=True, search_limit=20)
            sage: v.bundle()
            Tangent bundle of Minimal Planar Cover of Translation Surface in H_1(0) built from a square defined over Algebraic Real Field

        """
        if singularity_limit is not None:
            import warnings

            warnings.warn(
                "singularity_limit has been deprecated as a keyword argument for to_surface() and will be removed from a future version of sage-flatsurf"
            )

        surface = self.get_surface()

        point = vector(point, ring or surface.base_ring())

        if label is None:
            if search_all:
                if search_limit is None:
                    if surface.is_finite_type():
                        labels = surface.labels()
                    else:
                        raise ValueError(
                            "If search_all=True and the surface is infinite, then a search_limit must be provided."
                        )
                else:
                    from itertools import islice

                    labels = islice(surface.labels(), search_limit)
            else:
                labels = self.visible()
        else:
            labels = [label]

        points = set()

        for label in labels:
            gp = self.graphical_polygon(label)
            coords = gp.transform_back(point)

            pos = surface.polygon(label).get_point_position(coords)
            if not pos.is_inside():
                continue

            if v is None:
                points.add(
                    surface.point(label, coords, ring=ring, limit=singularity_limit)
                )
            else:
                direction = (~(gp.transformation().derivative())) * vector(v)
                if pos.is_vertex():
                    vertex = pos.get_vertex()
                    polygon = surface.polygon(label)

                    from flatsurf.geometry.euclidean import ccw

                    if ccw(polygon.edge(vertex), direction) < 0:
                        continue
                    if (
                        ccw(
                            polygon.edge((vertex - 1) % len(polygon.vertices())),
                            direction,
                        )
                        < 0
                    ):
                        continue

                points.add(
                    surface.tangent_vector(
                        label,
                        coords,
                        direction,
                        ring=ring,
                    )
                )
            if not return_all:
                return next(iter(points))

        if return_all:
            return points
        else:
            raise ValueError("Point or vector is not in a visible graphical polygon.")

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._ss.opposite_edge(p, e)

    def reset_letters(self, p, e):
        r"""
        Resets the letter dictionary for storing letters used in
        edge labeling if edge_labels="letter" is used.
        """
        try:
            del self._letters
            del self._next_letter
        except AttributeError:
            pass

    def _get_letter_for_edge(self, p, e):
        if not hasattr(self, "_letters"):
            self._letters = {}
            self._next_letter = 1
        try:
            return self._letters[(p, e)]
        except KeyError:
            # convert number to string
            nl = self._next_letter
            self._next_letter = nl + 1
            letter = ""
            while nl != 0:
                val = nl % 52
                if val == 0:
                    val = 52
                    letter = "Z" + letter
                elif val < 27:
                    letter = chr(97 + val - 1) + letter
                else:
                    letter = chr(65 + val - 27) + letter
                nl = (nl - val) / 52
            self._letters[(p, e)] = letter
            self._letters[self._ss.opposite_edge(p, e)] = letter
            return letter

    def edge_labels(self, lab):
        r"""
        Return the list of edge labels to be used for the polygon with label ``lab``.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(adjacencies=[])
            sage: g.edge_labels(0)
            ['1', '1', '1']
            sage: g.make_all_visible(adjacent=True)
            sage: g.edge_labels(0)
            [None, '1', '1']

            sage: g.make_adjacent(0,0)
            sage: g.edge_labels(0)
            [None, '1', '1']
            sage: g.edge_labels(1)
            ['0', None, '0']

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(adjacencies=[], edge_labels='number')
            sage: g.edge_labels(0)
            ['0', '1', '2']

            sage: g = s.graphical_surface(adjacencies=[], edge_labels='gluings and number')
            sage: g.edge_labels(0)
            ['0 -> (1, 1)', '1 -> (1, 2)', '2 -> (1, 0)']
            sage: g.make_all_visible(adjacent=True)
            sage: g.edge_labels(0)
            ['0', '1 -> (1, 2)', '2 -> (1, 0)']
        """
        if not self._edge_labels:
            return None

        s = self._ss
        g = self.graphical_polygon(lab)
        p = g.base_polygon()

        if self._edge_labels == "gluings":
            labels = []
            for e in range(len(p.vertices())):
                if self.is_adjacent(lab, e):
                    labels.append(None)
                elif s.opposite_edge(lab, e) is None:
                    labels.append(None)
                else:
                    llab, _ = s.opposite_edge(lab, e)
                    labels.append(str(llab))
        elif self._edge_labels == "number":
            labels = list(map(str, range(len(p.vertices()))))
        elif self._edge_labels == "gluings and number":
            labels = []
            for e in range(len(p.vertices())):
                if self.is_adjacent(lab, e):
                    labels.append(str(e))
                else:
                    labels.append("{} -> {}".format(e, s.opposite_edge(lab, e)))
        elif self._edge_labels == "letter":
            labels = []
            for e in range(len(p.vertices())):
                llab, ee = s.opposite_edge(lab, e)
                if not self.is_visible(llab) or self.is_adjacent(lab, e):
                    labels.append(None)
                else:
                    labels.append(self._get_letter_for_edge(lab, e))
        else:
            raise RuntimeError("invalid option for edge_labels")

        return labels

    def plot_polygon(self, label, graphical_polygon, upside_down):
        r"""
        Internal method for plotting polygons returning a Graphics object.

        Calls :meth:`.polygon.GraphicalPolygon.plot_polygon` passing
        the attribute ``upside_down_polygon_options`` if the polygon is upside down
        and ``polygon_options`` otherwise.

        Override this method for fine control of how the polygons are drawn.

        INPUT:

        - ``label`` -- The label of the polygon being plotted.

        - ``graphical_polygon`` -- The associated graphical polygon.

        - ``upside_down`` -- True if and only if the polygon will be rendered upside down.
        """
        if upside_down:
            return graphical_polygon.plot_polygon(**self.upside_down_polygon_options)
        else:
            return graphical_polygon.plot_polygon(**self.polygon_options)

    def plot_polygon_label(self, label, graphical_polygon, upside_down):
        r"""
        Internal method for plotting polygon labels returning a Graphics2D.

        Calls :meth:`.polygon.GraphicalPolygon.plot_label` passing the
        attribute ``polygon_label_options``.

        Override this method for fine control of how the polygons are drawn.

        INPUT:

        - ``label`` -- The label of the polygon being plotted.

        - ``graphical_polygon`` -- The associated graphical polygon.

        - ``upside_down`` -- True if and only if the polygon will be rendered upside down.
        """
        return graphical_polygon.plot_label(label, **self.polygon_label_options)

    def plot_edge(self, label, edge, graphical_polygon, is_adjacent, is_self_glued):
        r"""
        Internal method for plotting a polygon's edge returning a Graphics2D.

        The method calls :meth:`.polygon.GraphicalPolygon.plot_edge`.
        Depending on the geometry of the edge pair, it passes one of the attributes
        ``adjacent_edge_options``, ``self_glued_edge_options`` or ``non_adjacent_edge_options``.

        Override this method for fine control of how the edge is drawn.

        INPUT:

        - ``label`` -- The label of the polygon.

        - ``edge`` -- Integer representing the edge of the polygon.

        - ``graphical_polygon`` -- The associated graphical polygon.

        - ``is_adjacent`` -- True if and only if the polygon opposite this edge is visible and adjacent to this edge.
            In this case, plot_edge is called only once for this edge.

        - ``is_self_glued`` -- True if and only if the edge is glued to itself by a 180 degree rotation.
            This is never True if is_adjacent is True.
        """
        if is_adjacent:
            return graphical_polygon.plot_edge(edge, **self.adjacent_edge_options)
        elif is_self_glued:
            return graphical_polygon.plot_edge(edge, **self.self_glued_edge_options)
        else:
            return graphical_polygon.plot_edge(edge, **self.non_adjacent_edge_options)

    def plot_edge_label(self, p, e, edge_label, graphical_polygon):
        r"""
        Internal method for plotting an edge label.
        Calls :meth:`.polygon.GraphicalPolygon.plot_edge_label` passing
        the attribute ``edge_label_options``.

        Override this method for fine control of how the edge is drawn.

        INPUT:

        - ``p`` -- The label of the polygon.

        - ``e`` -- Integer representing the edge of the polygon.

        - ``edge_label`` -- A string containing the label to be printed on the edge.

        - ``graphical_polygon`` -- The associated graphical polygon.
        """
        return graphical_polygon.plot_edge_label(
            e, edge_label, **self.edge_label_options
        )

    def plot_zero_flag(self, label, graphical_polygon):
        r"""
        Internal method for plotting a polygon's zero_flag and returning a Graphics2D.
        Simply calls :meth:`.polygon.GraphicalPolygon.plot_zero_flag` passing
        the attribute` `zero_flag_options``.

        Override this method for fine control of how the edge is drawn.

        INPUT:

        - ``label`` -- The label of the polygon.

        - ``graphical_polygon`` -- The associated graphical polygon.
        """
        return graphical_polygon.plot_zero_flag(**self.zero_flag_options)

    def plot(self, **kwargs):
        r"""
        Return a plot of this surface.

        INPUT:

        - ``kwargs`` -- arguments are normally forwarded to the polygon
          plotting. However, prefixed arguments, e.g., ``polygon_label_color``,
          are routed correspondingly. Also, a dictionary suffixed with
          ``_options`` is merged with the existing options of this surface. See
          examples below for details.

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: from flatsurf.graphical.surface import GraphicalSurface
            sage: gs = GraphicalSurface(s)
            sage: gs.plot()
            ...Graphics object consisting of 13 graphics primitives

        Keyword arguments that end in ``_options`` are merged into the
        corresponding attribute before plotting; see :class:`GraphicalSurface`
        for a list of all supported ``_options``:

        .. jupyter-execute::

            sage: gs.plot(polygon_label_options={"color": "red"})
            ...Graphics object consisting of 13 graphics primitives

        Keyword arguments that are prefixed with such an aspect of plotting,
        are also merged into the corresponding attribute before plotting; see
        :class:`GraphicalSurface` for a list of all supported prefixes, i.e.,
        ``_options``:

        .. jupyter-execute::

            sage: gs.plot(polygon_label_color="red")
            ...Graphics object consisting of 13 graphics primitives

        All other arguments are passed to the polygon plotting itself:

        .. jupyter-execute::

            sage: gs.plot(fill=None)
            ...Graphics object consisting of 13 graphics primitives

        TESTS:

        Check that label options are handled correctly::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus()
            sage: S.plot(polygon_labels=True, edge_labels=True)
            ...Graphics object consisting of 10 graphics primitives
            sage: S.plot(polygon_labels=False, edge_labels=True)
            ...Graphics object consisting of 9 graphics primitives
            sage: S.plot(polygon_labels=True, edge_labels=False)
            ...Graphics object consisting of 6 graphics primitives
            sage: S.plot(polygon_labels=False, edge_labels=False)
            ...Graphics object consisting of 5 graphics primitives
        """
        if kwargs:
            surface = self.copy()

            options = [
                option
                for option in surface.__dict__
                if option.endswith("_options") and option != "process_options"
            ]
            # Sort recognized options so we do not pass polygon_label_color to
            # polygon_options as label_color.
            options.sort(key=lambda option: -len(option))

            for key, value in kwargs.items():
                if key in options:
                    setattr(surface, key, {**getattr(surface, key), **value})
                    continue

                for option in options:
                    prefix = option[: -len("options")]
                    if key.startswith(prefix) and key != prefix:
                        getattr(surface, option)[key[len(prefix) :]] = value
                        break
                else:
                    surface.polygon_options[key] = value

            return surface.plot()

        from sage.plot.graphics import Graphics

        p = Graphics()

        # Make sure we don't plot adjacent edges more than once.
        plotted_adjacent_edges = set()

        for label in self._visible:
            polygon = self.graphical_polygon(label)
            upside_down = polygon.transformation().sign() == -1

            # Plot the polygons
            if upside_down and self.will_plot_upside_down_polygons:
                p += self.plot_polygon(label, polygon, upside_down)
            elif self.will_plot_polygons:
                p += self.plot_polygon(label, polygon, upside_down)

            if self.will_plot_zero_flags:
                p += self.plot_zero_flag(label, polygon)

            # Add the polygon label
            if self.will_plot_polygon_labels:
                p += self.plot_polygon_label(label, polygon, upside_down)

            # Plot the edges
            if self.will_plot_edges:
                for i in range(len(self._ss.polygon(label).vertices())):
                    if self.is_adjacent(label, i):
                        if (
                            self.will_plot_adjacent_edges
                            and (label, i) not in plotted_adjacent_edges
                        ):
                            plotted_adjacent_edges.add(self._ss.opposite_edge(label, i))
                            p += self.plot_edge(label, i, polygon, True, False)
                    elif (label, i) == self._ss.opposite_edge(label, i):
                        # Self-glued edge
                        if self.will_plot_self_glued_edges:
                            p += self.plot_edge(label, i, polygon, False, True)
                    else:
                        if self.will_plot_non_adjacent_edges:
                            p += self.plot_edge(label, i, polygon, False, False)

            # Plot the edge labels.
            if self.will_plot_edge_labels:
                # get the edge labels
                edge_labels = self.edge_labels(label)
                if edge_labels is not None:
                    for i in range(len(self._ss.polygon(label).vertices())):
                        if edge_labels[i] is not None:
                            p += self.plot_edge_label(label, i, edge_labels[i], polygon)
        return p
