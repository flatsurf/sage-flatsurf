r"""
Two dimensional plot of surfaces

EXAMPLES::

The construction of a plot of a surface is split into three distinct
components:
- a layout (disposition of polygons in the plane)
- a labeller (choice of labels to use for polygons and edges)
- options (styles)
At the moment there is a single class for each of them but nothing prevents you
from writing your own layout, labeller or option handling::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.graphical.surface_plot import BasicLayout, SurfaceLabeller, SurfacePlotOptions

In order to plot a surface one needs to insantiate each of the above::

    sage: surface = translation_surfaces.arnoux_yoccoz(3)
    sage: u, v = surface.vertices()
    sage: layout = BasicLayout(surface)
    sage: labeller = SurfaceLabeller(surface)
    sage: option_handler = SurfacePlotOptions(polygons={"color": "pink", 0: {"color": "blue"}}, polygons_alpha=0.4,
    ....:              polygons_labels={0: {"string": "left wing"}, 1: {"string": "nef", "color": "firebrick"}, 2: {"string": "right wing"}},
    ....:              adjacent_edges={"color": "green"}, non_adjacent_edges_color="chartreuse", adjacent_edges_labels=False, edges_labels={(0, 0): {"string": "a"}},
    ....:               vertices={u: {"marker": "s", "color": "hotpink"}, v: {"marker": "h", "color": "lime"}})

The graphic is then made of a superposition of four different elements
- polygons rendered with sage ``polygon2d`` function,
- polygon labels rendered with sage ``text`` function,
- edge rendered with sage ``line2d`` function,
- edge labels rendered with sage ``text`` function.

The arguments to be sent to each function are built by the option handler
(respectively by the functions ``polygon_options``, ``polygon_label_options``,
`edge_options``, ``edge_label_options``). We illustrate this on a polygon and
an edge::

    sage: option_handler.polygon_options(0, layout)
    {'alpha': 0.400000000000000,
     'color': 'blue',
     'points': ((0, 0),
      (-1/2*alpha^2, -1/2*alpha^2 + 1/2),
      (-1/2, -1/2*alpha^2 + alpha - 1/2))}
    sage: polygon2d(**option_handler.polygon_options(0, layout))
    Graphics object consisting of 1 graphics primitive

    sage: option_handler.polygon_label_options(0, layout, labeller)
    {'color': 'gray',
     'string': 'left wing',
     'xy': (-1/6*alpha^2 - 1/6, -1/3*alpha^2 + 1/3*alpha)}
    sage: text(**option_handler.polygon_label_options(0, layout, labeller))
    Graphics object consisting of 1 graphics primitive

    sage: option_handler.edge_options(0, 0, layout)
    {'color': 'green',
     'linestyle': ':',
     'points': [(0, 0), (-1/2*alpha^2, -1/2*alpha^2 + 1/2)],
     'thickness': 0.5}
    sage: line2d(**option_handler.edge_options(0, 0, layout))
    Graphics object consisting of 1 graphics primitive

    sage: option_handler.edge_label_options(1, 0, layout, labeller)
    {'color': 'black',
     'horizontal_alignment': 'left',
     'string': 'a',
     'vertical_alignment': 'bottom',
     'xy': (-0.5942382152164581, 0.14630046089495533)}
    sage: text(**option_handler.edge_label_options(1, 0, layout, labeller))
    Graphics object consisting of 1 graphics primitive

    sage: u, v = surface.vertices()
    sage: option_handler.vertex_options(v, layout)
    {'color': 'lime',
     'marker': 'h',
     'points': [(1/2*alpha, 1/2*alpha^2 + alpha - 1/2),
      (-1/2, -1/2*alpha^2 + alpha - 1/2),
      (1/2, -1/2*alpha^2 + alpha - 1/2),
      (-1/2*alpha^2, -1/2*alpha^2 + 1/2),
      (-1/2*alpha^2, 1/2*alpha^2 - 1/2)],
     'size': 50}
    sage: point2d(**option_handler.vertex_options(v, layout))
    Graphics object consisting of 1 graphics primitive
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2024 Vincent Delecroix <20100.delecroix@gmail.com>
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

import collections

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.infinity import Infinity
from sage.modules.free_module_element import vector

from flatsurf.geometry.similarity import SimilarityGroup


# TODO: to be changed as Surface.plot()
def surface_plot(self, layout=None, labeller=None, option_handler=None, **kwds):
    r"""
    Return a 2d graphic representation of this surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.graphical.surface_plot import surface_plot
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: surface_plot(S)
        Graphics object consisting of ... graphics primitives
        sage: surface_plot(S, polygons={"color": "pink", 0: {"color": "blue"}}, polygons_alpha=0.4,
        ....:              polygons_labels={0: {"string": "left wing"}, 1: {"string": "nef", "color": "firebrick"}, 2: {"string": "right wing"}},
        ....:              adjacent_edges={"color": "green"}, non_adjacent_edges_color="chartreuse", adjacent_edges_labels=False, edges_labels={(0, 0): {"string": "a"}})
        Graphics object consisting of ... graphics primitives
    """
    from sage.plot.graphics import Graphics
    from sage.plot.line import line2d
    from sage.plot.point import point2d
    from sage.plot.polygon import polygon2d
    from sage.plot.text import text

    # TODO: customization
    if layout is None:
        layout = BasicLayout(self)

    # TODO: customization
    if labeller is None:
        labeller = SurfaceLabeller(self)

    # TODO: customization
    # TODO: do not choose vertex colors here!!
    if option_handler is None:
        from matplotlib import colormaps
        colormap = colormaps["hsv"]
        vertex_colors = {}
        vertices = self.vertices()
        for i, v in enumerate(vertices):
            vertex_colors[v] = {"color": colormap(float(i) / float(len(vertices)))[:3]}
        kwds["vertices"] = vertex_colors
        option_handler = SurfacePlotOptions(**kwds)
    elif kwds:
        raise ValueError

    G = Graphics()

    for label in layout.labels():
        opts = option_handler.polygon_options(label, layout)
        if opts is None:
            continue
        G += polygon2d(zorder=0, **opts)

        opts = option_handler.polygon_label_options(label, layout, labeller)
        if opts is not None:
            G += text(zorder=3, **opts)

    for label1, edge1, label2, edge2 in layout.adjacent_edges():
        opts1 = option_handler.edge_options(label1, edge1, layout)
        opts2 = option_handler.edge_options(label2, edge2, layout)
        points1 = opts1.pop("points")
        points2 = opts2.pop("points")
        if opts1 != opts2:
            raise ValueError("inconsistent options for adjacent edge ({}, {}), ({}, {})".format(label1, edge1, label2, edge2))
        if opts1 is None:
            continue
        G += line2d(points=points1, zorder=1, **opts1)

        opts1 = option_handler.edge_label_options(label1, edge1, layout, labeller)
        if opts1 is not None:
            G += text(zorder=3, **opts1)

        opts2 = option_handler.edge_label_options(label2, edge2, layout, labeller)
        if opts2 is not None:
            G += text(zorder=3, **opts2)

    for label, edge in layout.non_adjacent_edges():
        # NOTE: this also contains boundary edges, self-glued edges and edges
        # whose opposite edge belong to a hidden polygon
        opts = option_handler.edge_options(label, edge, layout)
        if opts is None:
            continue
        G += line2d(zorder=1, **opts)

        opts = option_handler.edge_label_options(label, edge, layout, labeller)
        if opts is None:
            continue
        G += text(zorder=3, **opts)

    for v in layout.vertices():
        opts = option_handler.vertex_options(v, layout)
        if opts is None:
            continue
        G += point2d(zorder=2, **opts)

    G.axes(False)
    G.set_aspect_ratio(1)
    return G


class SurfaceLayout:
    r"""
    Abstract class for surface layout.
    """
    def labels(self):
        r"""
        Return an iterator of the visible polygon labels.
        """
        raise NotImplementedError

    def adjacent_edges(self):
        r"""
        Iterator through quadruples ``(label1, edge1, label2, edge2)`` of pair of
        edges that are adjacent in this layout.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.graphical.surface_plot import BasicLayout
            sage: surface = translation_surfaces.arnoux_yoccoz(3)
            sage: layout = BasicLayout(surface, adjacencies=[(0, 0), (0, 1), (0, 2), (1, 2), (4, 0), (4, 1), (10, 0), (8, 0), (8, 1), (11, 0), (2, 1)])
            sage: sorted(layout.adjacent_edges())
            [(0, 0, 4, 2),
             (0, 1, 1, 1),
             (0, 2, 10, 2),
             (1, 2, 11, 2),
             (2, 1, 3, 1),
             (4, 0, 8, 2),
             (4, 1, 5, 1),
             (8, 0, 2, 0),
             (8, 1, 9, 1),
             (10, 0, 6, 2),
             (11, 0, 7, 2)]
        """
        visited = set()
        for label in self.labels():
            for edge in range(len(self._surface.polygon(label).vertices())):
                if (label, edge) in visited:
                    continue

                if self.is_adjacent(label, edge):
                    label2, edge2 = self._surface.opposite_edge(label, edge)
                    yield (label, edge, label2, edge2)
                    visited.add((label2, edge2))

    def non_adjacent_edges(self, vertex=None):
        r"""
        Iterator through ``(label, edge)`` that are edges which are not part
        of a pair of adjacent edges in this layout.

        If ``vertex`` is provided, only yield edges adjacent to it.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.graphical.surface_plot import BasicLayout
            sage: surface = translation_surfaces.arnoux_yoccoz(3)
            sage: layout = BasicLayout(surface, adjacencies=[(0, 0), (0, 1), (0, 2), (1, 2), (4, 0), (4, 1), (10, 0), (8, 0), (8, 1), (11, 0), (2, 1)])
            sage: sorted(layout.non_adjacent_edges())
            [(1, 0),
             (2, 2),
             (3, 0),
             (3, 2),
             (5, 0),
             (5, 2),
             (6, 0),
             (6, 1),
             (7, 0),
             (7, 1),
             (9, 0),
             (9, 2),
             (10, 1),
             (11, 1)]
        """
        for label in self.labels():
            for edge in range(len(self._surface.polygon(label).vertices())):
                opposite_edge = self._surface.opposite_edge(label, edge)
                if opposite_edge is None:
                    continue
                label2, edge2 = opposite_edge
                if not self.is_adjacent(label, edge):
                    if vertex is None or self._surface.point(label, self._surface.polygon(label).vertex(edge)) == vertex:
                        yield (label, edge)

    def vertices(self):
        r"""
        Return the vertices in the surface that are visible.
        """
        verts = set()
        for label, edge in self.non_adjacent_edges():
            verts.add(self._surface.point(label, self._surface.polygon(label).vertex(edge)))
        return verts

    def layout(self, label):
        raise NotImplementedError

    def polygon(self, label):
        if not self.is_visible(label):
            raise ValueError("invisible polygon")
        return self.layout(label)(self._surface.polygon(label))

    def is_visible(self, label):
        raise NotImplementedError

    def make_visible(self, label):
        raise NotImplementedError

    def hide(self, label):
        raise NotImplementedError

    def is_adjacent(self, label, edge):
        r"""
        Test whether the edge with the given polygon label and edge number is
        adjacent to its paired edge in this layout.
        """
        surface = self._surface
        opposite_edge = surface.opposite_edge(label, edge)
        if opposite_edge is None:
            return False
        label2, edge2 = opposite_edge
        if self.is_visible(label) and self.is_visible(label2):
            p = self.polygon(label)
            p2 = self.polygon(label2)
            return (p.vertex(edge) == p2.vertex(edge2 + 1) and
                    p.vertex(edge + 1) == p2.vertex(edge2))
        return False


class BasicLayout(SurfaceLayout):
    r"""
    Layout for similarity surfaces.

    A layout is a partial map from the polygons defining a similarity surface
    and similarities in the plane. Eaach similarity corresponds to a choice of
    a given disposition of the polygon in the plane. Note that the map is
    partial so that some of the polygons may be set invisible (in particular
    when the surface is infinite).
    """
    def __init__(self, surface, adjacencies=None):
        self._surface = surface
        self._group = SimilarityGroup(surface.base_ring())
        self._transformations = {}

        # NOTE: the labels in _transformations are the visible polygons
        r = next(iter(self._surface.labels()))
        self._transformations[r] = self._group.one()

        if adjacencies is not None:
            for label, edge in adjacencies:
                self.make_adjacent(label, edge)
        else:
            self.discover()

    def layout(self, label):
        r"""
        Return the similarity for the polygon with the given ``label``.
        """
        return self._transformations[label]

    def is_visible(self, label):
        return label in self._transformations

    def labels(self):
        return self._transformations.keys()

    def make_visible(self, label):
        r"""
        Make the polygon with the given ``label`` visible.

        The function tries to make it adjacent to one of its neighbors. If
        none of the neighboring polygon has a layout set, then it is plot
        positionned according to its coordinates.
        """
        for edge in range(len(self._surface.polygon(label).vertices())):
            label2, edge2 = self._surface.opposite_edge(label, edge)
            if (label2, edge2) in self._transformations:
                self.make_adjacent(label2, edge2)
                return

        # default
        self._transformation[label] = self._group.one()

    def discover(self, limit=None):
        r"""
        Discover the surface and set the layout.
        """
        if limit is None:
            limit = Infinity

        if limit == Infinity:
            if not self._surface.is_finite_type():
                raise NotImplementedError
            labels = self._surface.labels()
        else:
            labels = itertools.islice(self._surface.labels(), limit)

        if not self._transformations:
            self._transformations[next(self._surface.labels())] = self._group.one()

        # simple BFS
        todo = collections.deque(self._transformations)
        i = 0
        while todo and i < limit:
            label = todo.popleft()
            poly = self._surface.polygon(label)
            for edge in range(len(poly.vertices())):
                opposite_edge = self._surface.opposite_edge(label, edge)
                if opposite_edge is None:
                    continue
                label2, edge2 = opposite_edge
                if label2 not in self._transformations:
                    self.make_adjacent(label, edge)
                    todo.append(label2)
                    i += 1

    def hide(self, label):
        self._transformations.pop(label)

    def make_adjacent(self, p, e, reverse=False):
        if p not in self._transformations:
            raise ValueError("no layout set for polygon with label={}".format(p))
        pp, ee = self._surface.opposite_edge(p, e)
        if reverse:
            q = self._surface.polygon(p)
            a = q.vertex(e)
            b = q.vertex(e + 1)
            # This is the similarity carrying the origin to a and (1,0) to b:
            g = self._group(b[0] - a[0], b[1] - a[1], a[0], a[1])

            qq = self._surface.polygon(pp)
            aa = qq.vertex(ee + 1)
            bb = qq.vertex(ee)
            # This is the similarity carrying the origin to aa and (1,0) to bb:
            gg = self._group(bb[0] - aa[0], bb[1] - aa[1], aa[0], aa[1])

            reflection = G(
                self._surface.base_ring().one(),
                self._surface.base_ring().zero(),
                self._surface.base_ring().zero(),
                self._surface.base_ring().zero(),
                -1,
            )

            # This is the similarity carrying (a,b) to (aa,bb):
            g = gg * reflection * (~g)
        else:
            g = self._surface.edge_transformation(pp, ee)

        self._transformations[pp] = self._transformations[p] * g


class SurfaceLabeller:
    def __init__(self, surface):
        self._surface = surface
        self._letters = {}
        self._next_letter = 1

    def polygon_label(self, label):
        return str(label)

    def edge_label(self, label, edge):
        try:
            return self._letters[(label, edge)]
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
            self._letters[(label, edge)] = letter
            self._letters[self._surface.opposite_edge(label, edge)] = letter
            return letter


class SurfacePlotOptions:
    r"""
    Proxy between the layout, labeller and the sage plot code.

    EXAMPLES::

        sage: from flatsurf.graphical.surface_plot import SurfacePlotOptions
        sage: SurfacePlotOptions(adjacent_edges_labels_color="red")
        ...
        adjacent_edges_labels               : True
        ...
        adjacent_edges_labels_options       : {'color': 'red'}
        ...

    Individual options do not bypass flags::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.graphical.surface_plot import BasicLayout, SurfaceLabeller, SurfacePlotOptions
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: layout = BasicLayout(S, adjacencies=[(0, 0)])
        sage: labeller = SurfaceLabeller(S)
        sage: SP = SurfacePlotOptions(adjacent_edges_labels=False, edges_labels={(0, 0): {"string": "my label"}})
        sage: SP.edge_label_options(0, 0, layout, labeller)
        sage: SP = SurfacePlotOptions(adjacent_edges_labels=True, edges_labels={(0, 0): {"string": "my label"}})
        sage: SP.edge_label_options(0, 0, layout, labeller)
        {'color': 'black',
         'horizontal_alignment': 'right',
         'string': 'my label',
         'vertical_alignment': 'top',
         'xy': (-0.05356065604458417, 0.10179080927503625)}
    """
    def __init__(self, **kwds):
        # polygon options
        self._polygons = True
        self._folded_polygons = False

        self._polygons_options = {
            "color": "lightgray"
        }
        self._folded_polygons_options = {
            "color": "whitesmoke"
        }

        self._individual_polygons_options = {}

        # polygon label options
        self._polygons_labels = True
        self._folded_polygons_labels = False

        self._polygons_labels_options = {
            "color": "gray"
        }
        self._folded_polygons_labels_options = {}

        self._individual_polygons_labels_options = {}

        # edge/ray options
        self._edges = True
        self._adjacent_edges = True
        self._non_adjacent_edges = True
        self._boundary_edges = True
        self._self_glued_edges = True

        self._edges_options = {
            "color": "blue"
        }
        self._adjacent_edges_options = {
            "linestyle": ":",
            "thickness": 0.5
        }
        self._non_adjacent_edges_options = {}
        self._self_glued_edges_options = {
            "color": "green"
        }
        self._boundary_edges_options = {
            "color": "cyan",
            "thickness": 2
        }

        self._individual_edges_options = {}

        # edge label options
        self._adjacent_edges_labels = False
        self._non_adjacent_edges_labels = True
        self._boundary_edges_labels = False
        self._self_glued_edges_labels = False

        self._edges_labels_options = {
            "color": "black"
        }
        self._adjacent_edges_labels_options = {}
        self._non_adjacent_edges_labels_options = {}
        self._boundary_edges_labels_options = {}
        self._self_glued_edges_labels_options = {}

        self._individual_edges_labels_options = {}

        # vertex options
        self._vertices = True
        self._vertices_options = {
            "marker": "o",
            "color": "black",
            "size": 50
        }

        self._individual_vertices_options = {}

        self.set_options(**kwds)

    def set_options(self, **kwds):
        special_keys = ('polygons', 'folded_polygons', 'edges', 'boundary_edges', 'adjacent_edges', 'non_adjacent_edges', 'self_glued_edges', 'vertices')

        for skey in special_keys:
            for key in [skey, skey + '_labels']:
                value = kwds.pop(key, None)
                if value is None:
                    continue
                elif isinstance(value, bool):
                    setattr(self, '_' + key, value)
                elif isinstance(value, dict):
                    for k, v in value.items():
                        if isinstance(v, dict):
                            getattr(self, '_individual_' + key + '_options')[k] = v
                        elif isinstance(k, str):
                            setattr(self, '_' + key, True)
                            getattr(self, '_' + key + '_options')[k] = v
                        else:
                            raise ValueError("invalid option for {} got {}: {}".format(key, k, v))

        for key in list(kwds):
            for skey in special_keys:
                skeylab = skey + '_labels'
                if key.startswith(skeylab + '_'):
                    setattr(self, '_' + skeylab, True)
                    getattr(self, '_' + skeylab + '_options')[key[len(skeylab) + 1:]] = kwds.pop(key)

                elif key.startswith(skey + '_'):
                    setattr(self, '_' + skey, True)
                    getattr(self, '_' + skey + '_options')[key[len(skey) + 1:]] = kwds.pop(key)

        if kwds:
            raise ValueError("unknown keyword arguments: {}".format(", ".join(map(str, kwds))))

    def __str__(self):
        s = ["SurfacePlotOptions"]
        attrs = [# boolean attributes
                 "polygons", "folded_polygons",
                 "edges", "adjacent_edges", "non_adjacent_edges", "boundary_edges", "self_glued_edges",
                 "polygons_labels", "folded_polygons_labels",
                 "adjacent_edges_labels", "non_adjacent_edges_labels", "boundary_edges_labels", "self_glued_edges_labels",
                 # dictionaries
                 "polygons_options", "folded_polygons_options",  "individual_polygons_options",
                 "polygons_labels_options", "folded_polygons_labels_options", "individual_polygons_labels_options",
                 "edges_options", "adjacent_edges_options",  "non_adjacent_edges_options", "self_glued_edges_options", "boundary_edges_options", "individual_edges_options",
                 "edges_labels_options", "adjacent_edges_labels_options", "non_adjacent_edges_labels_options", "boundary_edges_labels_options", "self_glued_edges_labels_options",
                 "individual_edges_labels_options"]
        n = max(map(len, attrs)) + 1
        formatter = "{{:{}}}".format(n)
        for attr in attrs:
            s.append("  " + formatter.format(attr) + " : " + str(getattr(self, '_' + attr)))
        return "\n".join(s)

    __repr__ = __str__

    def polygon_options(self, label, layout=None):
        r"""
        Return ``None`` when the polygon should not be plotted or a dictionary of options.

       If not ``None``, the dictionary in the output has a key ``"points"`` which is the list of
       positions of the vertices to be rendered (which is the required argument of the sage
       ``polygon2d`` function).
        """
        if not self._polygons or not layout.is_visible(label):
            return
        opts = self._polygons_options.copy()
        opts.update(self._individual_polygons_options.get(label, []))
        if "points" not in opts:
            if layout is None:
                raise ValueError("a layout is needed as \"points\" is not preset for polygon with label {}".format(label))
            opts["points"] = layout.polygon(label).vertices()
        return opts

    def polygon_label_options(self, label, layout=None, labeller=None):
        r"""
        Return ``None`` when the polygon label should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output always has a key ``"string"`` and a key ``"xy"``
        which are the string label and position to be rendered (which are the required
        arguments of the sage ``text`` function).
        """
        if not self._polygons_labels or not layout.is_visible(label):
            return
        opts = self._polygons_labels_options.copy()
        opts.update(self._individual_polygons_labels_options.get(label, []))
        if "string" not in opts:
            if labeller is None:
                raise ValueError("a labeller is needed as \"string\" is not preset for polygon with label {}".format(label))
            opts["string"] = labeller.polygon_label(label)
        if "xy" not in opts:
            if layout is None:
                raise ValueError("a layout is needed as \"xy\" is not preset for polygon with label {}".format(label))
            verts = layout.polygon(label).vertices()
            opts["xy"] = sum(verts) / len(verts)
        return opts

    def _edge_opt(self, layout, label, edge, boundary, self_glued, adjacent, non_adjacent, boundary_options, self_glued_options, adjacent_options, non_adjacent_options, individual_options):
        opposite_edge = layout._surface.opposite_edge(label, edge)
        if opposite_edge is None:
            # boundary
            return non_adjacent_options if boundary else None
        else:
            label2, edge2 = opposite_edge
            if label == label2 and edge == edge2:
                # self-glued
                return self_glued_edge_options if self_glued else None
            elif layout.is_adjacent(label, edge):
                # adjacent
                return adjacent_options if adjacent else None
            else:
                # non-adjacent
                return non_adjacent_options if non_adjacent else None

        raise RuntimeError

    def edge_options(self, label, edge, layout=None):
        r"""
        Return ``None`` when the edge should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output has a key ``"points"`` which is the 2-tuple of
        start and end to be rendered (which is the required argument of the sage ``line2d``
        function).
        """
        if not layout.is_visible(label):
            return
        opts = self._edge_opt(layout, label, edge,
                              self._boundary_edges, self._self_glued_edges, self._adjacent_edges, self._non_adjacent_edges,
                              self._boundary_edges_options, self._self_glued_edges_options, self._adjacent_edges_options, self._non_adjacent_edges_options, self._individual_edges_options)
        if opts is None:
            return None

        opts0 = self._edges_options.copy()
        opts0.update(opts)
        opts = opts0

        opts.update(self._individual_edges_options.get((label, edge), []))

        if "points" not in opts:
            if layout is None:
                raise ValueError("a layout is needed as \"points\" is not preset for edge ({}, {})".format(polygon, label))
            p = layout.polygon(label)
            opts["points"] = [p.vertex(edge), p.vertex(edge + 1)]
        return opts

    # TODO: the positionning of the edge_label is delicate
    # matplotlib does not offer much
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.text.html#matplotlib.pyplot.text
    def edge_label_options(self, label, edge, layout=None, labeller=None):
        r"""
        Return ``None`` when the edge label should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output has a key ``"string"`` and a key ``"xy"``
        which are the string label and position to be rendered (which are the required
        arguments of the sage ``text`` function).
        """
        if not layout.is_visible(label):
            return
        opts = self._edge_opt(layout, label, edge,
                              self._boundary_edges_labels, self._self_glued_edges_labels, self._adjacent_edges_labels, self._non_adjacent_edges_labels,
                              self._boundary_edges_labels_options, self._self_glued_edges_labels_options, self._adjacent_edges_labels_options, self._non_adjacent_edges_labels_options, self._individual_edges_labels_options)

        if opts is None:
            return opts

        opts0 = self._edges_labels_options.copy()
        opts0.update(opts)
        opts = opts0

        opts.update(self._individual_edges_labels_options.get((label, edge), []))

        if "xy" not in opts:
            if layout is None:
                raise ValueError("a layout is needed as \"xy\" is not preset for the label of edge ({}, {})".format(label, edge))
            p = layout.polygon(label)
            v0 = p.vertex(edge)
            v1 = p.vertex(edge + 1)
            edge_vector = v1 - v0
            t = 0.3
            push_off = 0.01
            normal_vector = vector(RDF, (-edge_vector[1], edge_vector[0]))
            normal_vector /= normal_vector.norm()
            pos = (1 - t) * p.vertex(edge) + t * p.vertex(edge + 1) + push_off * normal_vector
            opts["xy"] = pos
        if "horizontal_alignment" not in opts:
            if edge_vector[1] < 0:
                opts["horizontal_alignment"] = "left"
            elif edge_vector[1] > 0:
                opts["horizontal_alignment"] = "right"
            else:
                opts["horizontal_alignment"] = "center"
        if "vertical_alignment" not in opts:
            if edge_vector[0] < 0:
                opts["vertical_alignment"] = "top"
            elif edge_vector[0] > 0:
                opts["vertical_alignment"] = "bottom"
            else:
                opts["vertical_alignment"] = "center"
        if "string" not in opts:
            if labeller is None:
                raise ValueError("a \"labeller\" is needed as \"string\" is not preset for the label of edge ({}, {})".format(label, edge))
            opts["string"] = labeller.edge_label(label, edge)
            edge_label = labeller.edge_label(label, edge)
            if edge_label is None:
                return
        return opts

    def vertex_options(self, vertex, layout=None):
        r"""
        Return ``None`` if the vertex ``v`` should not be plotted otherwise a dictionary of options.
        """
        if not self._vertices:
            return

        opts = self._vertices_options.copy()
        opts.update(self._individual_vertices_options.get(vertex, []))
        if "points" not in opts:
            if layout is None:
                raise ValueError("a layout is needed as \"points\" is not preset for vertex {}".format(vertex))
            points = [layout.polygon(label).vertex(edge) for (label, edge) in layout.non_adjacent_edges(vertex)]
            opts["points"] = points
        return opts
