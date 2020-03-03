r"""
EXAMPLES::

    sage: import flatsurf
    sage: flatsurf.translation_surfaces.veech_2n_gon(4).plot()
    Graphics object consisting of 18 graphics primitives
"""
#*****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from flatsurf.geometry.similarity_surface import SimilaritySurface
from .polygon import *

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
    `process_options()`.

    The basic tasks of the class are to render the polygons, edges and labels.
    To customize a rendering, it is useful to know something about how this
    class works. (Apologies!)

    There are attributes which control whether or not certain objects are
    rendered, namely:

    - `will_plot_polygons` -- Whether to plot polygons which are right-side up.

    - `will_plot_upside_down_polygons` -- Whether to plot polygons which are
        upside down. Defaults to False.

    - `will_plot_polygon_labels` -- Whether to plot polygon labels.

    - `will_plot_edges` -- If this is False then no edges will be plotted.

    - `will_plot_non_adjacent_edges` = Whether to plot polygon edges which are
        not adjacent to the edge it is glued to.

    - `will_plot_adjacent_edges` -- Whether to plot polygon edges which are
        adjacent to the polygon they are glued to.

    - `will_plot_self_glued_edges` -- Whether to plot polygon edges which are
        glued to themselves.

    - `will_plot_edge_labels` -- Whether to plot polygon edges labels.

    - `will_plot_zero_flags` -- Whether to plot line segments from the
        baricenter to the zero vertex of each polygon. Useful in working out
        vertex and edge labels. Defaults to False.

    The `plot` method calls some other built in methods: `plot_polygon`,
    `plot_polygon_label`, `plot_edge` and `plot_edge_label`. These in turn
    call methods in `GraphicalPolygon`.

    - `polygon_options` -- Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_polygon` when
        plotting a polygon right-side up.

    - `upside_down_polygon_options` -- Options passed to
        :func:`graphical_polygon.GraphicalPolygon.plot_polygon` when plotting a polygon upside-down.

    - `polygon_label_options` -- Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_label`
        when plotting a polygon label.

    - `non_adjacent_edge_options` -- Options passed to
        :func:`graphical_polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge which is not
        adjacent to the edge it is glued to.

    - `self.adjacent_edge_options` -- Options passed to
        :func:`graphical_polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge which is
        adjacent to the edge it is glued to.

    - `self_glued_edge_options` -- Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_edge`
        when plotting a polygon edge which is glued to itself.

    - `edge_label_options` -- Options passed to :func:`graphical_polygon.GraphicalPolygon.edge_label`
        when plotting a edge label.

    - `zero_flag_options` -- Options passed to
        :func:`graphical_polygon.GraphicalPolygon.plot_zero_flag` when plotting a zero_flag.

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.graphical.surface import GraphicalSurface

        sage: s = similarity_surfaces.example()
        sage: gs = GraphicalSurface(s)
        sage: gs.polygon_options["color"]="red"
        sage: gs.plot()     # not tested (problem with matplotlib font caches on Travis)
        Graphics object consisting of 13 graphics primitives
    """

    def __init__(self, similarity_surface, adjacencies=None, polygon_labels=True, \
                 edge_labels="gluings", default_position_function = None):
        r"""
        Construct a GraphicalSurface from a similarity surface.

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
        """
        assert isinstance(similarity_surface, SimilaritySurface)
        self._ss = similarity_surface
        self._default_position_function = default_position_function
        self._polygons = {}
        self._visible = set([self._ss.base_label()])

        if adjacencies is None:
            if self._ss.is_finite():
                self.make_all_visible()
        self._edge_labels = None

        self.will_plot_polygons = True
        r"""
        Whether to plot polygons which are right-side up.
        """

        self.polygon_options = {"color":"lightgray"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_polygon` when plotting a polygon right-side up."""

        self.will_plot_upside_down_polygons = False
        r"""
        Whether to plot polygons which are upside down
        """

        self.upside_down_polygon_options = {"color":"lightgray", "zorder":-1}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_polygon` when plotting a polygon upside-down."""

        self.will_plot_polygon_labels = True
        r"""
        Whether to plot polygon labels.
        """

        self.polygon_label_options = {"color":"black", "vertical_alignment":"center", "horizontal_alignment":"center"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_label` when plotting a polygon label."""

        self.will_plot_edges = True
        r"""
        If this is False then no edges will be plotted.
        """

        self.will_plot_non_adjacent_edges = True
        r"""
        Whether to plot polygon edges which are not adjacent to the edge it is glued to.
        """

        self.non_adjacent_edge_options = {"color":"blue"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge
        which is not adjacent to the edge it is glued to."""

        self.will_plot_adjacent_edges = True
        r"""
        Whether to plot polygon edges which are adjacent to the polygon they are glued to.
        """

        self.adjacent_edge_options = {"color":"blue", "linestyle":":"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_edge`
        when plotting a polygon edge which is adjacent to the edge it is glued to."""

        self.will_plot_self_glued_edges = True
        r"""
        Whether to plot polygon edges which are glued to themselves.
        """

        self.self_glued_edge_options = {"color":"red"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_edge` when plotting a polygon edge
        which is glued to itself."""

        self.will_plot_edge_labels = True
        r"""
        Whether to plot polygon edges which are not adjacent to the polygon they are glued to.
        """

        self.edge_label_options = {"color":"blue"}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.edge_label` when plotting a polygon label."""

        self.will_plot_zero_flags = False
        r"""
        Whether to plot line segments from the baricenter to the zero vertex of each polygon.
        """

        self.zero_flag_options = {"color":"green", "thickness":0.5}
        r"""Options passed to :func:`graphical_polygon.GraphicalPolygon.plot_zero_flag` when plotting a zero_flag."""

        self.process_options(adjacencies=adjacencies,
                polygon_labels=polygon_labels, edge_labels=edge_labels)

    def process_options(self, adjacencies=None, polygon_labels=None, edge_labels=None, default_position_function = None):
        r"""
        Process the options listed as if the graphical_surface was first
        created.

        INPUT:

        - ``adjacencies`` -- a list of pairs ``(p,e)`` to be used to set
          adjacencies of polygons.

        - ``polygon_labels`` -- a boolean (default ``True``) whether the label
          of polygons are displayed

        - ``edge_labels`` -- option to control the display of edge labels. It
          can be one of

            - ``None`` for no change

            - ``False`` for no labels

            - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
              name of the polygon to which it is glued

            - ``'number'`` -- to put on each side of each edge the number of the
              edge

            - ``'gluings and number'`` -- full information

            - ``'letter'`` -- add matching letters to glued edges in an arbitrary way

        - ``default_position_function`` -- a function mapping polygon labels to
          similarities describing the position of the corresponding polygon.
          Note that this will not affect polygons which have already been
          positioned.

        TESTS::

            sage: from flatsurf import *

            sage: c = translation_surfaces.chamanara(1/2)
            sage: gs = c.graphical_surface()
            sage: gs.process_options(edge_labels='hey')
            Traceback (most recent call last):
            ...
            ValueError: invalid value for edge_labels (='hey')
        """
        if adjacencies is not None:
            for p,e in adjacencies:
                self.make_adjacent(p, e)
        if polygon_labels is not None:
            if not isinstance(polygon_labels, bool):
                raise ValueError("polygon_labels must be True, False or None.")
            self.will_plot_polygon_labels = polygon_labels
        if edge_labels is not None:
            if edge_labels is True:
                self.will_plot_edge_labels = True
                edge_labels = 'gluings'
            elif edge_labels is False:
                self.will_plot_edge_labels = False
                edge_labels = None
            elif edge_labels in ['gluings', 'number', 'gluings and number', 'letter']:
                self._edge_labels = edge_labels
                self.will_plot_edge_labels = True
            else:
                raise ValueError("invalid value for edge_labels (={!r})".format(edge_labels))
        if default_position_function is not None:
            self._default_position_function = default_position_function

    def copy(self):
        r"""
        Make a copy of this GraphicalSurface.

        EXAMPLES::

            sage: from flatsurf import *
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
        gs = GraphicalSurface(self.get_surface(), default_position_function = self._default_position_function)

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
        gs._polygons = {label:gp.copy() for label,gp in iteritems(self._polygons)}
        gs._visible = set(self._visible)
        gs._edge_labels = self._edge_labels

        return gs

    def __repr__(self):
        return "Graphical version of Similarity Surface {!r}".format(self._ss)

    def visible(self):
        r"""
        Return a copy of the set of visible labels.
        """
        return set(self._visible)

    def is_visible(self,label):
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

        EXAMPLES::

            sage: from flatsurf import *

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface()
            sage: g.make_all_visible()
            sage: g.plot()      # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 13 graphics primitives

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(cached=False, adjacencies=[])
            sage: g.make_all_visible(adjacent=False)
            sage: g.plot()      # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 16 graphics primitives
        """
        if adjacent is None:
            adjacent = (self._default_position_function is None)
        if limit is None:
            assert self._ss.is_finite()
            if adjacent:
                for l,poly in self._ss.walker().label_polygon_iterator():
                    for e in range(poly.num_edges()):
                        l2,e2 = self._ss.opposite_edge(l,e)
                        if not self.is_visible(l2):
                            self.make_adjacent(l,e)
            else:
                from flatsurf.geometry.similarity import SimilarityGroup
                T = SimilarityGroup(self._ss.base_ring())
                for l in self._ss.label_iterator():
                    if not self.is_visible(l):
                        if self._default_position_function is None:
                            # No reasonable way to display the polygon, so we do this hack:
                            g = self.graphical_polygon(l)
                            poly = self._ss.polygon(l)
                            sxmax = self.xmax()
                            pxmin = g.xmin()
                            t = T((QQ(self.xmax() - g.xmin() + 1),
                                QQ(-(g.ymin()+g.ymax()) / ZZ(2) )))
                            g.set_transformation(t)
                        self.make_visible(l)
        else:
            assert limit>0
            if adjacent:
                i = 0
                for l,poly in self._ss.walker().label_polygon_iterator():
                    for e in range(poly.num_edges()):
                        l2,e2 = self._ss.opposite_edge(l,e)
                        if not self.is_visible(l2):
                            self.make_adjacent(l,e)
                            i=i+1
                            if i>=limit:
                                return
            else:
                from flatsurf.geometry.similarity import SimilarityGroup
                T = SimilarityGroup(self._ss.base_ring())
                i = 0
                for l in self._ss.label_iterator():
                    if not self.is_visible(l):
                        if self._default_position_function is None:
                            # No reasonable way to display the polygon, so we do this hack:
                            g = self.graphical_polygon(l)
                            poly = self._ss.polygon(l)
                            sxmax = self.xmax()
                            pxmin = g.xmin()
                            t = T((QQ(self.xmax() - g.xmin() + 1),
                                QQ(-(g.ymin()+g.ymax()) / ZZ(2) )))
                            g.set_transformation(t)
                        self.make_visible(l)
                        i=i+1
                        if i>=limit:
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
        maximal x-and y- cordinates  of a visible graphical polygon.
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
            if not self._default_position_function is None:
                t=self._default_position_function(label)
            p = GraphicalPolygon(self._ss.polygon(label), transformation=t)
            self._polygons[label] = p
            return p

    def make_adjacent(self, p, e, reverse = False, visible = True):
        r"""
        Move the polygon across the prescribed edge so that is adjacent.
        The polygon moved is obtained from opposite_edge(p,e).

        If reverse=True then the polygon is moved so that there is a fold
        at the edge.

        If visible is True (by default), we also make the moved polygon visible.

        EXAMPLES::

            sage: from flatsurf import *
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
        pp,ee = self._ss.opposite_edge(p,e)
        if reverse:
            from flatsurf.geometry.similarity import SimilarityGroup
            G = SimilarityGroup(self._ss.base_ring())

            q = self._ss.polygon(p)
            a = q.vertex(e)
            b = q.vertex(e+1)
            # This is the similarity carrying the origin to a and (1,0) to b:
            g = G(b[0]-a[0],b[1]-a[1],a[0],a[1])

            qq = self._ss.polygon(pp)
            aa = qq.vertex(ee+1)
            bb = qq.vertex(ee)
            # This is the similarity carrying the origin to aa and (1,0) to bb:
            gg = G(bb[0]-aa[0],bb[1]-aa[1],aa[0],aa[1])

            reflection = G(
                self._ss.base_ring().one(),
                self._ss.base_ring().zero(),
                self._ss.base_ring().zero(),
                self._ss.base_ring().zero(),
                -1)

            # This is the similarity carrying (a,b) to (aa,bb):
            g = gg*reflection*(~g)
        else:
            g = self._ss.edge_transformation(pp,ee)
        h = self.graphical_polygon(p).transformation()
        self.graphical_polygon(pp).set_transformation(h*g)
        if visible:
            self.make_visible(pp)

    def make_adjacent_and_visible(self, p, e, reverse=False):
        r"""
        Move the polygon across the prescribed edge so that is adjacent,
        and make the moved polygon visible.
        """
        from sage.misc.superseded import deprecation
        deprecation(42, "Do not use .make_adjacent_and_visible(). Use .make_adjacent() instead.")
        self.make_adjacent(p, e, reverse=reverse)

    def is_adjacent(self,p,e):
        r"""
        Returns the truth value of the statement
        'The polygon opposite edge (p,e) is adjacent to that edge.'

        EXAMPLES::

            sage: from flatsurf import *
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
        """
        pp,ee = self.opposite_edge(p,e)
        if not self.is_visible(pp):
            return False
        g = self.graphical_polygon(p)
        gg = self.graphical_polygon(pp)
        return g.transformed_vertex(e) == gg.transformed_vertex(ee+1) and \
               g.transformed_vertex(e+1) == gg.transformed_vertex(ee)

    def to_surface(self, point, v=None, label=None, ring=None, return_all=False, \
                   singularity_limit=None, search_all = False, search_limit=None):
        r""" Converts from graphical coordinates to similarity surface coordinates.

        A point always must be provided. If a vector v is provided then a
        SimilaritySurfaceTangentVector will be returned. If v is not provided, then a
        SurfacePoint is returned.

        INPUT:

        - ``point`` -- Coordinates of a point in graphical coordinates to be
            converted to graphical coordinates.

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
            to ``SimilaritySurface.surface_point``.

        - ``search_all`` -- (default ``False``) By default we look just in
            polygons with visible label. If set to `True``, then we instead
            look in all labels.

        - ``search_limit`` -- (default ``None``) If ``search_all`` is ``True``,
            then we look at the first ``search_limit`` polygons instead of all
            polygons. This must be set to an positive integer if ``search_all``
            is and the surface is infinite.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: gs = s.graphical_surface()
            sage: gs.to_surface((1,-2))
            Surface point located at (1, 1/2) in polygon 1
            sage: gs.to_surface((1,-2), v=(1,0))
            SimilaritySurfaceTangentVector in polygon 1 based at (1, 1/2) with vector (1, -1/2)

            sage: s = translation_surfaces.infinite_staircase()
            sage: gs = s.graphical_surface()
            sage: gs.to_surface((4,4), (1,1), search_all=True, search_limit=20)
            SimilaritySurfaceTangentVector in polygon 8 based at (0, 0) with vector (1, 1)

            sage: s = translation_surfaces.square_torus()
            sage: pc = s.minimal_cover(cover_type="planar")
            sage: gs = pc.graphical_surface()
            sage: gs.to_surface((3,2), search_all=True, search_limit=20)
            Traceback (most recent call last):
            ...
            ValueError: To obtain a singularity on an infinite surface, singularity_limit must be set.
            sage: gs.to_surface((3,2), search_all=True, search_limit=20, singularity_limit=4)
            Surface point with 4 coordinate representations
            sage: p = gs.to_surface((sqrt(3),sqrt(2)), ring=AA, search_all=True, search_limit=20)
            sage: next(iter(p.coordinates(p.labels()[0]))).parent()
            Vector space of dimension 2 over Algebraic Real Field
            sage: v = gs.to_surface((3/2,3/2),(sqrt(3),sqrt(2)),ring=AA,search_all=True, search_limit=20)
            sage: v.bundle()
            Tangent bundle of TranslationSurface built from infinitely many polygons defined over Algebraic Real Field
        """
        if label is None:
            if return_all:
                ret = set()
            s = self.get_surface()
            if search_all:
                if search_limit is None:
                    if s.is_finite():
                        it = s.label_iterator()
                    else:
                        raise ValueError("If search_all=True and the surface is infinite, then a search_limit must be provided.")
                else:
                    from itertools import islice
                    it = islice(s.label_iterator(), search_limit)
            else:
                it = self.visible()
            for label in it:
                try:
                    val = self.to_surface(point, v=v, label=label, ring=ring, singularity_limit=singularity_limit)
                    if return_all:
                        ret.add(val)
                    else:
                        return val
                except AssertionError:
                    # Not in the polygon
                    pass
                except ValueError as e:
                    if e.args[0] == 'need a limit when working with an infinite surface':
                        raise ValueError("To obtain a singularity on an infinite surface, " + \
                            "singularity_limit must be set.")
                    # Otherwise it is not in the polygon.
            if return_all:
                return ret
            else:
                raise ValueError("Point or vector is not in a visible graphical_polygon.")
        else:
            gp = self.graphical_polygon(label)
            coords = gp.transform_back(point)
            s = self.get_surface()
            if v is None:
                return s.surface_point(label, coords, ring=ring, limit=singularity_limit)
            else:
                return s.tangent_vector(label, coords, (~(gp.transformation().derivative()))*vector(v), ring=ring)

    def opposite_edge(self, p, e):
        r"""
        Given the label ``p`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``pp``, ``ee``) to which this edge is glued.
        """
        return self._ss.opposite_edge(p,e)

    def reset_letters(self,p,e):
        r"""
        Resets the letter dictionary for storing letters used in
        edge labeling if edge_labels="letter" is used.
        """
        try:
            del self._letters
            del self._next_letter
        except:
            pass

    def _get_letter_for_edge(self, p, e):
        if not hasattr(self,"_letters"):
            self._letters={}
            self._next_letter=1
        try:
            return self._letters[(p,e)]
        except KeyError:
            # convert number to string
            nl = self._next_letter
            self._next_letter = nl + 1
            letter = ""
            while nl!=0:
                val = nl % 52
                if val==0:
                    val=52
                    letter = "Z" + letter
                elif val<27:
                    letter = chr(97+val-1) + letter
                else:
                    letter = chr(65+val-27) + letter
                nl = (nl-val)/52
            self._letters[(p,e)] = letter
            self._letters[self._ss.opposite_edge(p,e)] = letter
            return letter

    def edge_labels(self, lab):
        r"""
        Return the list of edge labels to be used for the polygon with label ``lab``.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(cached=False, adjacencies=[])
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
            sage: g = s.graphical_surface(cached=False, adjacencies=[], edge_labels='number')
            sage: g.edge_labels(0)
            ['0', '1', '2']

            sage: g = s.graphical_surface(cached=False, adjacencies=[], edge_labels='gluings and number')
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

        if self._edge_labels == 'gluings':
            ans = []
            for e in range(p.num_edges()):
                if self.is_adjacent(lab, e):
                    ans.append(None)
                else:
                    llab,ee = s.opposite_edge(lab,e)
                    ans.append(str(llab))
        elif self._edge_labels == 'number':
            ans = list(map(str, range(p.num_edges())))
        elif self._edge_labels == 'gluings and number':
            ans = []
            for e in range(p.num_edges()):
                if self.is_adjacent(lab, e):
                    ans.append(str(e))
                else:
                    ans.append("{} -> {}".format(e, s.opposite_edge(lab,e)))
        elif self._edge_labels == "letter":
            ans = []
            for e in range(p.num_edges()):
                llab,ee = s.opposite_edge(lab,e)
                if not self.is_visible(llab) or self.is_adjacent(lab, e):
                    ans.append(None)
                else:
                    ans.append(self._get_letter_for_edge(lab,e))
        else:
            raise RuntimeError("invalid option for edge_labels")

        return ans

    def plot_polygon(self, label, graphical_polygon, upside_down):
        r"""
        Internal method for plotting polygons returning a Graphics object.

        Calls :func:`graphical_polygon.GraphicalPolygon.plot_polygon` passing
        the attribute `upside_down_polygon_options` if the polygon is upside down
        and `polygon_options` otherwise.

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

        Calls :func:`graphical_polygon.GraphicalPolygon.plot_polygon_label` passing
        the attribute `polygon_label_options`.

        Override this method for fine control of how the polygons are drawn.

        INPUT:

        - ``label`` -- The label of the polygon being plotted.

        - ``graphical_polygon`` -- The associated graphical polygon.

        - ``upside_down`` -- True if and only if the polygon will be rendered upside down.
        """
        return graphical_polygon.plot_label(label,**self.polygon_label_options)

    def plot_edge(self, label, edge, graphical_polygon, is_adjacent, is_self_glued):
        r"""
        Internal method for plotting a polygon's edge returning a Graphics2D.

        The method calls :func:`graphical_polygon.GraphicalPolygon.plot_edge`.
        Depending on the geometry of the edge pair, it passes one of the attributes
        `adjacent_edge_options`, `self_glued_edge_options` or `non_adjacent_edge_options`.

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
        Calls :func:`graphical_polygon.GraphicalPolygon.plot_edge_label` passing
        the attribute `edge_label_options`.

        Override this method for fine control of how the edge is drawn.

        INPUT:

        - ``p`` -- The label of the polygon.

        - ``e`` -- Integer representing the edge of the polygon.

        - ``edge_label`` -- A string containing the label to be printed on the edge.

        - ``graphical_polygon`` -- The associated graphical polygon.
        """
        return graphical_polygon.plot_edge_label(e, edge_label, **self.edge_label_options)

    def plot_zero_flag(self, label, graphical_polygon):
        r"""
        Internal method for plotting a polygon's zero_flag and returning a Graphics2D.
        Simply calls :func:`graphical_polygon.GraphicalPolygon.plot_zero_flag` passing
        the attribute `zero_flag_options`.

        Override this method for fine control of how the edge is drawn.

        INPUT:

        - ``label`` -- The label of the polygon.

        - ``graphical_polygon`` -- The associated graphical polygon.
        """
        return graphical_polygon.plot_zero_flag(**self.zero_flag_options)


    def plot(self):
        r"""
        Returns a plot of the GraphicalSurface

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: from flatsurf.graphical.surface import GraphicalSurface
            sage: gs = GraphicalSurface(s)
            sage: gs.make_visible(1)
            sage: gs.plot()      # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 13 graphics primitives


        Check that label options are handled correctly::

            sage: S = translation_surfaces.square_torus()
            sage: S.plot(polygon_labels=True, edge_labels=True)     # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 10 graphics primitives
            sage: S.plot(polygon_labels=False, edge_labels=True)    # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 9 graphics primitives
            sage: S.plot(polygon_labels=True, edge_labels=False)    # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 6 graphics primitives
            sage: S.plot(polygon_labels=False, edge_labels=False)   # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 5 graphics primitives
        """
        from sage.plot.graphics import Graphics
        p = Graphics()

        # Make sure we don't plot adjacent edges more than once.
        plotted_adjacent_edges = set()

        for label in self._visible:
            polygon = self.graphical_polygon(label)
            upside_down = polygon.transformation().sign()==-1

            # Plot the polygons
            if upside_down and self.will_plot_upside_down_polygons:
                    p += self.plot_polygon(label, polygon, upside_down)
            elif self.will_plot_polygons:
                p += self.plot_polygon(label, polygon, upside_down)

            if self.will_plot_zero_flags:
                p += self.plot_zero_flag(label,polygon)

            # Add the polygon label
            if self.will_plot_polygon_labels:
                p += self.plot_polygon_label(label, polygon, upside_down)

            # Plot the edges
            if self.will_plot_edges:
                for i in range(self._ss.polygon(label).num_edges()):
                    if self.is_adjacent(label,i):
                        if self.will_plot_adjacent_edges and (label,i) not in plotted_adjacent_edges:
                            plotted_adjacent_edges.add(self._ss.opposite_edge(label,i))
                            p += self.plot_edge(label, i, polygon, True, False)
                    elif (label,i) == self._ss.opposite_edge(label,i):
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
                    for i in range(self._ss.polygon(label).num_edges()):
                        if edge_labels[i] is not None:
                            p += self.plot_edge_label(label, i, edge_labels[i], polygon)
        return p


