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

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector

from flatsurf.geometry.similarity import SimilarityGroup


class DynamicalLayout:
    r"""
    Note: the layout could just be a wrapper on top of a similarity surface
    """
    def __init__(self, surface):
        self._surface = surface
        self._group = SimilarityGroup(surface.base_ring())
        self._transformations = {}

        r = next(iter(self._surface.labels()))
        self._transformations[r] = self._group.one()
        self.make_adjacencies_if_unset()

    def is_adjacent(self, label, edge):
        r"""
        Test whether the edge with the given polygon label and edge number is
        adjacent to its paired edge.
        """
        opposite_edge = self._surface.opposite_edge(label, edge)
        if opposite_edge is None:
            return False
        label2, edge2 = opposite_edge
        g = self[label]
        g2 = self[label2]
        p = self._surface.polygon(label)
        p2 = self._surface.polygon(label2)
        return (g is not None and
                g2 is not None and
                g(p.vertex(edge)) == g2(p2.vertex(edge2 + 1)) and
                g(p.vertex(edge + 1)) == g2(p2.vertex(edge2)))

    def make_adjacent(self, p, e, reverse=False):
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

    def make_adjacencies_if_unset(self, adjacencies=None):
        if adjacencies is None:
            adjacencies = ((label, e) for label, p in zip(self._surface.labels(), self._surface.polygons()) for e in range(len(self._surface.vertices())))
        for p, e in adjacencies:
            opposite_edge = self._surface.opposite_edge(p, e)
            if opposite_edge is None:
                continue
            pp, ee = opposite_edge
            if p in self._transformations and pp not in self._transformations:
                self.make_adjacent(p, e)

    def __getitem__(self, p):
        try:
            return self._transformations[p]
        except KeyError:
            return self._group.one()


class LetterEdgeLabeller:
    def __init__(self, surface):
        self._surface = surface
        self._letters = {}
        self._next_letter = 1

    def __getitem__(self, key):
        try:
            return self._letters[key]
        except KeyError:
            # convert number to string
            p, e = key
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
            self._letters[self._surface.opposite_edge(p, e)] = letter
            return letter


class SurfacePlot:
    r"""
    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.graphical.surface_plot import SurfacePlot
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: SP = SurfacePlot(S, polygons={"color": "pink", 0: {"color": "blue"}}, polygons_alpha=0.4,
        ....:                  polygons_labels={0: {"string": "left wing"}, 1: {"string": "nef", "color": "firebrick"}, 2: {"string": "right wing"}},
        ....:                  adjacent_edges={"color": "green"}, non_adjacent_edges_color="chartreuse", adjacent_edges_labels=False, edges_labels={(0, 0): {"string": "a"}})
        sage: SP.polygon_options(0)
        {'alpha': 0.400000000000000,
         'color': 'blue',
         'points': [(0, 0),
          (-1/2*alpha^2, -1/2*alpha^2 + 1/2),
          (-1/2, -1/2*alpha^2 + alpha - 1/2)]}
        sage: polygon2d(**SP.polygon_options(0))
        Graphics object consisting of 1 graphics primitive
        sage: SP.polygon_label_options(0)
        {'color': 'black',
         'string': 'left wing',
         'xy': (-1/6*alpha^2 - 1/6, -1/3*alpha^2 + 1/3*alpha)}
        sage: text(**SP.polygon_label_options(0))
        Graphics object consisting of 1 graphics primitive
        sage: SP.edge_options(0, 0)
        {'color': 'green',
         'linestyle': ':',
         'points': [(0, 0), (-1/2*alpha^2, -1/2*alpha^2 + 1/2)],
         'thickness': 0.5}
        sage: line2d(**SP.edge_options(0, 0))
        Graphics object consisting of 1 graphics primitive
        sage: SP.edge_label_options(0, 2)
        {'color': 'black', 'string': 'a', 'xy': (-1/4, -1/4*alpha^2 + 1/2*alpha - 1/4)}
        sage: text(**SP.edge_label_options(0, 2))
        Graphics object consisting of 1 graphics primitive

        sage: SP = SurfacePlot(S, adjacent_edges_labels_color="red")
        sage: SP
        ...
        adjacent_edges_labels               : True
        ...
        adjacent_edges_labels_options       : {'color': 'red'}
        ...

    Individual options bypass flags::

        sage: SP = SurfacePlot(S, adjacent_edges_labels=False, non_adjacent_edges_labels=False, edges_labels={(0, 0): {"string": "my label"}})
        sage: SP.edge_label_options(0, 0)['string']
        'my label'
    """
    def __init__(self, surface, **kwds):
        self._surface = surface
        self._layout = DynamicalLayout(surface)
        self._edge_labels = LetterEdgeLabeller(surface)

        self._visible = set(surface.labels())

        # polygon options
        self._polygons = True
        self._folded_polygons = False

        self._polygons_options = {"color": "lightgray"}
        self._folded_polygons_options = {"color": "whitesmoke"}

        self._individual_polygons_options = {}

        # polygon label options
        self._polygons_labels = True
        self._folded_polygons_labels = False

        self._polygons_labels_options = {"color": "black"}
        self._folded_polygons_labels_options = {}

        self._individual_polygons_labels_options = {}

        # edge options
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

        special_keys = ('polygons', 'folded_polygons', 'edges', 'boundary_edges', 'adjacent_edges', 'non_adjacent_edges')

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
        s = ["SurfacePlot of {}".format(self._surface)]
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

    def visible_labels(self):
        r"""
        Iterable of labels of polygons to be displayed
        """
        return self._visible

    def is_visible(self, label):
        return label in self._visible

    def make_visible(self, label):
        return self._visible.add(label)

    def hide(self, label):
        self._visible.remove(label)

    def polygon_options(self, label):
        r"""
        Return ``None`` when the polygon should not be plotted or a dictionary of options.

       If not ``None``, the dictionary in the output has a key ``"points"`` which is the list of
       positions of the vertices to be rendered (which is the required argument of the sage
       ``polygon2d`` function).
        """
        if not self._polygons and label not in self._individual_polygons_options:
            return
        opts = self._polygons_options.copy()
        opts.update(self._individual_polygons_options.get(label, []))
        if "points" not in opts:
            similarity = self._layout[label]
            opts["points"] = [similarity(v) for v in self._surface.polygon(label).vertices()]
        return opts

    def polygon_label_options(self, label):
        r"""
        Return ``None`` when the polygon label should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output always has a key ``"string"`` and a key ``"xy"``
        which are the string label and position to be rendered (which are the required
        arguments of the sage ``text`` function).
        """
        if not self._polygons_labels and label not in self._individual_polygons_labels_options:
            return
        opts = self._polygons_labels_options.copy()
        opts.update(self._individual_polygons_labels_options.get(label, []))
        if "string" not in opts:
            opts["string"] = str(label)
        if "xy" not in opts:
            similarity = self._layout[label]
            verts = [similarity(v) for v in self._surface.polygon(label).vertices()]
            opts["xy"] = sum(verts) / len(verts)
        return opts

    def _edge_opt(self, label, edge, boundary, self_glued, adjacent, non_adjacent, boundary_options, self_glued_options, adjacent_options, non_adjacent_options, individual_options):
        opposite_edge = self._surface.opposite_edge(label, edge)
        if (label, edge) in individual_options:
            boundary = self_glued = adjacent = non_adjacent = True
        if opposite_edge is None:
            # boundary
            return non_adjacent_options if boundary else None
        else:
            label2, edge2 = opposite_edge
            if label == label2 and edge == edge2:
                # self-glued
                return self_glued_edge_options if self_glued else None
            elif self._layout.is_adjacent(label, edge):
                # adjacent
                return adjacent_options if adjacent else None
            else:
                # non-adjacent
                return non_adjacent_options if non_adjacent else None

        raise RuntimeError

    def edge_options(self, label, edge):
        r"""
        Return ``None`` when the edge should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output has a key ``"points"`` which is the 2-tuple of
        start and end to be rendered (which is the required argument of the sage ``line2d``
        function).
        """
        opts = self._edge_opt(label, edge,
                              self._boundary_edges, self._self_glued_edges, self._adjacent_edges, self._non_adjacent_edges,
                              self._boundary_edges_options, self._self_glued_edges_options, self._adjacent_edges_options, self._non_adjacent_edges_options, self._individual_edges_options)
        if opts is None:
            return None

        opts0 = self._edges_options.copy()
        opts0.update(opts)
        opts = opts0

        opts.update(self._individual_edges_options.get((label, edge), []))

        if "points" not in opts:
            similarity = self._layout[label]
            p = self._surface.polygon(label)
            opts["points"] = [similarity(p.vertex(edge)), similarity(p.vertex(edge + 1))]
        return opts

    def edge_label_options(self, label, edge):
        r"""
        Return ``None`` when the edge label should not be plotted or a dictionary of options

        If not ``None``, the dictionary in the output has a key ``"string"`` and a key ``"xy"``
        which are the string label and position to be rendered (which are the required
        arguments of the sage ``text`` function).
        """
        opts = self._edge_opt(label, edge,
                              self._boundary_edges_labels, self._self_glued_edges_labels, self._adjacent_edges_labels, self._non_adjacent_edges_labels,
                              self._boundary_edges_labels_options, self._self_glued_edges_labels_options, self._adjacent_edges_labels_options, self._non_adjacent_edges_labels_options, self._individual_edges_labels_options)

        if opts is None:
            return opts

        opts0 = self._edges_labels_options.copy()
        opts0.update(opts)
        opts = opts0

        opts.update(self._individual_edges_labels_options.get((label, edge), []))

        if "xy" not in opts:
            p = self._surface.polygon(label)
            similarity = self._layout[label]
            opts["xy"] = similarity((p.vertex(edge) + p.vertex(edge + 1)) / 2)
        if "string" not in opts:
            opts["string"] = self._edge_labels[label, edge]
        return opts

    def plot(self):
        from sage.plot.graphics import Graphics
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d
        from sage.plot.text import text

        G = Graphics()

        for label in self._visible:
            opts = self.polygon_options(label)
            if opts is not None:
                G += polygon2d(**opts)
                opts = self.polygon_label_options(label)
                if opts is not None:
                    G += text(**opts)

            for edge in range(len(self._surface.polygon(label).vertices())):
                opts = self.edge_options(label, edge)
                if opts is not None:
                    G += line2d(**opts)
                    opts = self.edge_label_options(label, edge)
                    if opts is not None:
                        G += text(**opts)

        G.axes(False)
        G.set_aspect_ratio(1)
        return G


