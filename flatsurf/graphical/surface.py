r"""
Rendering of surfaces

.. NOTE::

    The documentation of sage-flatsurf contains a section of example plots that
    showcase all the ways in which plots can be customized.

EXAMPLES::

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


class GraphicalSurface:
    r"""
    Renders a surface.

    This is a base class to customize rendering of surfaces. Subclasses handle
    the details of rendering surfaces defined by euclidean polygons and
    surfaces defined by hyperbolic polygons.

    Objects of this type should not be created manually, instead use
    ``graphical_surface()`` and ``plot()`` on a surface.

    EXAMPLES::

        sage: import flatsurf
        sage: flatsurf.translation_surfaces.veech_2n_gon(4).graphical_surface()
        Graphical representation of Translation Surface in H_2(2) built from a regular octagon

    """

    DEFAULTS = {
        "polygon": {
            "axes": False,
            "color": "lightgray",
            "edgecolor": None,
        },
        "point": {
            # We choose z-order so that points are above everything else,
            # including text, and axes, see
            # https://matplotlib.org/3.1.1/gallery/misc/zorder_demo.html
            "zorder": 4,
        },
        "zero_flag": {"color": "green", "thickness": 0.5},
        "polygon_label": {
            "color": "black",
            "vertical_alignment": "center",
            "horizontal_alignment": "center",
        },
        "adjacent_edge": {
            "color": "blue",
            "linestyle": ":",
        },
        "non_adjacent_edge": {
            "color": "blue",
        },
        "self_glued_edge": {
            "color": "red",
        },
        "boundary_edge": {
            "color": "grey",
            "linestyle": "--",
        },
        "adjacent_edge_label": {},
        "non_adjacent_edge_label": {
            "color": "blue",
            "label": "letter",
        },
        "self_glued_edge_label": {},
        "boundary_edge_label": {
            "color": "grey",
            "label": "number",
        },
    }

    def __init__(self, surface, adjacencies=None, default_position_function=None, **kwargs):
        if adjacencies is not None:
            import warnings
            warnings.warn("The adjacencies keyword has been removed from this version of sage-flatsurf. Loop over the labels and call layout() instead.")

        if default_position_function is not None:
            import warnings
            warnings.warn("The default_position_function keyword has been removed from this version of sage-flatsurf. Call layout() instead or override _layout_one in a subclass.")

        self._surface = surface
        self._transformation = {}

        self._options = {key: dict(value) for (key, value) in self.DEFAULTS.items()}
        self._options["zero_flag"] = None

        self.apply_options(**kwargs)

    def copy(self):
        copy = type(self)(self._surface)
        copy._transformation = dict(self._transformation)

        from copy import deepcopy
        copy._options = deepcopy(self._options)
        return copy

    def get_surface(self):
        r"""
        Return the underlying surface.
        """
        # TODO: Deprecate
        return self._surface

    def _ensure_layout(self):
        if not self._transformation:
            self.layout(limit=None if self._surface.is_finite_type() else 16)

    def visible(self):
        r"""
        Return a copy of the set of visible labels.
        """
        # TODO: This should be replaced with more specific predicates.
        self._ensure_layout()
        return set(self._transformation)

    def is_visible(self, label):
        r"""
        Return whether the polygon with the given label is marked as visible.
        """
        # TODO: This should be replaced with more specific predicates.
        self._ensure_layout()
        return label in self._transformation

    def is_adjacent(self, label, edge):
        self._ensure_layout()

        opposite = self._surface.opposite_edge(label, edge)
        if opposite is None:
            return False

        if opposite[0] not in self._transformation:
            return False

        if opposite[0] == label:
            return opposite[1] == edge

        return tuple(self._graphical_vertices(label, edge)) == tuple(self._graphical_vertices(*opposite)[::-1])


    def is_polygon_visible(self, label):
        self._ensure_layout()
        return label in self._transformation and self._polygon_options(label)

    def is_zero_flag_visible(self, label):
        self._ensure_layout()
        return label in self._transformation and self._zero_flag_options(label)

    def is_polygon_label_visible(self, label):
        self._ensure_layout()
        return label in self._transformation and self._polygon_label_options(label)

    def is_edge_visible(self, label, edge):
        self._ensure_layout()

        if label not in self._transformation:
            return False

        options = self._edge_options(label, edge)

        if options is None:
            return False

        opposite = self._surface.opposite_edge(label, edge)
        if opposite is None:
            return True

        if opposite[0] == label:
            return True

        if opposite[0] not in self._transformation:
            return True

        if self.is_adjacent(label, edge):
            for lbl in self._surface.labels():
                if lbl == label:
                    return True
                if lbl == opposite[0]:
                    return False

            assert False, "label not in surface"

        return True

    def is_edge_label_visible(self, label, edge):
        return self._edge_label_options(label, edge) and self._edge_label(label, edge)

    def _polygon_options(self, label):
        return self._options["polygon"]

    def _zero_flag_options(self, label):
        return self._options["zero_flag"]

    def _polygon_label_options(self, label):
        return self._options["polygon_label"]

    def _edge_options(self, label, edge):
        opposite = self._surface.opposite_edge(label, edge)
        if opposite is None:
            return self._options["boundary_edge"]
        if opposite == (label, edge):
            return self._options["self_glued_edge"]
        if self.is_adjacent(label, edge):
            return self._options["adjacent_edge"]

        return self._options["non_adjacent_edge"]

    def _edge_label_options(self, label, edge):
        opposite = self._surface.opposite_edge(label, edge)
        if opposite is None:
            return self._options["boundary_edge_label"]
        elif opposite == (label, edge):
            return self._options["self_glued_edge_label"]
        elif self.is_adjacent(label, edge):
            return self._options["adjacent_edge_label"]

        return self._options["non_adjacent_edge_label"]

    def _edge_label(self, label, edge):
        r"""
        Return the label for the ``edge`` of the polygon with ``label``.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface()
            sage: g.layout(limit=1)
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            ['A', 'B', 'C']

            sage: g.layout()
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            ['A', 'B', None]

            sage: g.make_adjacent(0, 0)
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            [None, 'A', 'B']

            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface(non_adjacent_edge_label_label="number")
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            ['0', '1', None]

            sage: g = s.graphical_surface(non_adjacent_edge_label_label="gluings and number")
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            ['0 => (1, 1)', '1 => (1, 2)', None]
            sage: g = s.graphical_surface(non_adjacent_edge_label_label="gluings and number", adjacent_edge_label_label="number")
            sage: [g._edge_label(0, edge) for edge in [0, 1, 2]]
            ['0 => (1, 1)', '1 => (1, 2)', '2']

        """
        style = self._edge_label_options(label, edge).get("label", None)

        if not style:
            return None

        if style == "gluings":
            opposite_label, _ = self._surface.opposite_edge(label, edge)
            return str(opposite_label)

        elif style == "number":
            return str(edge)

        elif style == "gluings and number":
            opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)
            return f"{edge} => ({opposite_label}, {opposite_edge})"

        elif style == "letter":
            return self._edge_label_letter(label, edge)

        else:
            raise ValueError("unsupported edge label style")

    def _edge_label_letter(self, label, edge):
        # TODO: This is very slow. We need a local (explicit) cache here.
        letter = 0
        seen = set()

        # TODO: Iterate over visible polygons only.
        for ((l, e), (ol, oe)) in self._surface.gluings():
            if (ol, oe) in seen:
                continue
            if self.is_adjacent(l, e):
                continue

            seen.add((l, e))
            if (label, edge) == (l, e) or (label, edge) == (ol, oe):
                break

            letter += 1
        else:
            assert False, f"{(label, edge)} not glued"

        alphabet = [chr(65 + i) for i in range(26)] + [chr(97 + i) for i in range(26)]
        word = ''
        while True:
            word += alphabet[letter % 52]
            letter //= 52
            if not letter:
                return word

    def _point_options(self, point):
        return self._options["point"]

    def __repr__(self):
        return "Graphical representation of {!r}".format(self._surface)

    def apply_options(self, **kwargs):
        options = list(self._options)
        # Sort recognized options so we do not pass polygon_label_color to
        # polygon_options as label_color.
        options.sort(key=lambda option: -len(option))

        if "edge_labels" in kwargs:
            import warnings
            warnings.warn("edge_labels has been deprecated as a parameter to plots and will be removed in a future version of sage-flatsurf. Use adjacent_edge_labels, non_adjacent_edge_labels, self_glued_edge_labels, and boundary_edge_labels instead")

            edge_labels = kwargs.pop("edge_labels")
            if not edge_labels:
                # disable all edge labels
                kwargs["adjacent_edge_label_options"] = {}
                kwargs["non_adjacent_edge_label_options"] = {}
                kwargs["self_glued_edge_label_options"] = {}
                kwargs["boundary_edge_label_options"] = {}
            elif edge_labels is True:
                # enable default edge labels
                kwargs["adjacent_edge_label_options"] = self.DEFAULTS["adjacent_edge_label"]
                kwargs["non_adjacent_edge_label_options"] = self.DEFAULTS["non_adjacent_edge_label"]
                kwargs["self_glued_edge_label_options"] = self.DEFAULTS["self_glued_edge_label"]
                kwargs["boundary_edge_label_options"] = self.DEFAULTS["boundary_edge_label"]
            else:
                # use a fixed kind of edge labels everywhere (legacy compatibility)
                kwargs["adjacent_edge_label_options"] = {"color": "blue", "label": edge_labels}
                kwargs["non_adjacent_edge_label_options"] = {"color": "blue", "label": edge_labels}
                kwargs["self_glued_edge_label_options"] = {"color": "blue", "label": edge_labels}
                kwargs["boundary_edge_label_options"] = {"color": "blue", "label": edge_labels}

        for key, value in kwargs.items():
            for option in options:
                # Handle flags such as plot_zero_flags, zero_flag_options, and zero_flags
                if key in [f"plot_{option}s", f"{option}s", f"{option}_options"]:
                    if isinstance(value, dict):
                        self._options[option] = value
                    elif value:
                        self._options[option] = dict(type(self).DEFAULTS[option])
                    else:
                        self._options[option] = {}
                    break
                # Handle flags such as polygon_label_color.
                if key.startswith(option + "_"):
                    self._options[option][key[len(option) + 1:]] = value
                    break
            else:
                # The argument matched no recognized pattern, pass it to polygon plotting if possible.
                if key in ["alpha", "edgecolor", "fill", "hue", "legend_color", "legend_label", "linestyle", "rgbcolor", "thickness", "zorder"]:
                    self._options["polygon"][key] = value
                else:
                    raise ValueError(f"unrecognized option {key}")

    def make_adjacent(self, label, edge):
        return self.layout(label=label, algorithm="adjacent", edge=edge)

    def layout(self, label=None, limit=None, algorithm=None, **kwargs):
        if label is None:
            self._layout_many(limit=limit, algorithm=algorithm, **kwargs)
        else:
            if limit is not None:
                raise ValueError("limit must not be set when a label has been specified")
            self._layout_one(label=label, algorithm=algorithm, **kwargs)

    def _layout_many(self, limit=None, algorithm=None, **kwargs):
        if limit is None and not self._surface.is_finite_type():
            raise ValueError("limit must be set for surfaces of infinite type")

        labels = self._surface.labels()
        if limit is not None:
            from itertools import islice
            labels = islice(labels, limit)

        for label in labels:
            self.layout(label=label, algorithm=algorithm, **kwargs)

    def _layout_one(self, label, algorithm=None, **kwargs):
        if algorithm is None and label in self._transformation:
            return

        if algorithm is None:
            algorithm = "adjacent"

        if algorithm == "adjacent":
            self._layout_one_adjacent(label, **kwargs)
        else:
            raise NotImplementedError("unknown ayout algorithm")

    def hide(self, label):
        raise NotImplementedError

    def make_all_visible(self, adjacent=None, limit=None):
        import warnings
        warnings.warn("The make_all_visible() function has been deprecated and will be removed in a future version of sage-flatsurf. Use layout() instead.")

        self.layout(limit=limit)

    def plot(self, **kwargs):
        r"""
        Return a plot of this surface.

        INPUT:

        - ``kwargs`` -- arguments are normally forwarded to the polygon
          plotting. However, prefixed arguments, e.g., ``polygon_label_color``,
          are routed correspondingly. Also, a dictionary suffixed with
          ``_options`` is merged with the existing options of this surface.
          Finally, arguments prefixed with ``plot_`` control whether certain
          bits are plotted at all. See examples below for details.

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: from flatsurf.graphical.surface import GraphicalSimilaritySurface
            sage: gs = GraphicalSimilaritySurface(s)
            sage: gs.plot()
            ...Graphics object consisting of 13 graphics primitives

        Keyword arguments that end in ``_options`` are merged into the
        corresponding attribute before plotting; see :meth:`options` for all
        supported ``_options``::

            sage: gs.plot(polygon_label_options={"color": "red"})
            ...Graphics object consisting of 13 graphics primitives

        Keyword arguments that are prefixed with such an aspect of plotting,
        are also merged into the corresponding attribute before plotting; see
        :class:`GraphicalSurface` for a list of all supported prefixes, i.e.,
        ``_options``::

            sage: gs.plot(polygon_label_color="red")
            ...Graphics object consisting of 13 graphics primitives

        Keyword arguments that start with ``plot_``, control whether a certain
        aspect of the plotting happens or not::

            sage: gs.plot(plot_zero_flags=True)
            ...Graphics object consisting of 15 graphics primitives

        All other arguments are passed to the polygon plotting itself::

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
            self = self.copy()
            self.apply_options(**kwargs)

        self._ensure_layout()

        from sage.plot.graphics import Graphics

        g = Graphics()

        g += self._plot_polygons()
        g += self._plot_zero_flags()
        g += self._plot_polygon_labels()
        g += self._plot_edges()
        g += self._plot_edge_labels()

        return g

        # TODO
        # # Make sure we don't plot adjacent edges more than once.
        # plotted_adjacent_edges = set()

        # for label in self._visible:
        #     polygon = self.graphical_polygon(label)
        #     upside_down = polygon.transformation().sign() == -1

        #     # Plot the polygons
        #     if upside_down and self.will_plot_upside_down_polygons:
        #         p += self.plot_polygon(label, polygon, upside_down)
        #     elif self.will_plot_polygons:
        #         p += self.plot_polygon(label, polygon, upside_down)

        #     if self.will_plot_zero_flags:
        #         p += self.plot_zero_flag(label, polygon)

        #     # Add the polygon label
        #     if self.will_plot_polygon_labels:
        #         p += self.plot_polygon_label(label, polygon, upside_down)

        #     # Plot the edges
        #     if self.will_plot_edges:
        #         for i in range(len(self._ss.polygon(label).vertices())):
        #             if self.is_adjacent(label, i):
        #                 if (
        #                     self.will_plot_adjacent_edges
        #                     and (label, i) not in plotted_adjacent_edges
        #                 ):
        #                     plotted_adjacent_edges.add(self._ss.opposite_edge(label, i))
        #                     p += self.plot_edge(label, i, polygon, True, False)
        #             elif (label, i) == self._ss.opposite_edge(label, i):
        #                 # Self-glued edge
        #                 if self.will_plot_self_glued_edges:
        #                     p += self.plot_edge(label, i, polygon, False, True)
        #             else:
        #                 if self.will_plot_non_adjacent_edges:
        #                     p += self.plot_edge(label, i, polygon, False, False)

        #     # Plot the edge labels.
        #     if self.will_plot_edge_labels:
        #         # get the edge labels
        #         edge_labels = self.edge_labels(label)
        #         if edge_labels is not None:
        #             for i in range(len(self._ss.polygon(label).vertices())):
        #                 if edge_labels[i] is not None:
        #                     p += self.plot_edge_label(label, i, edge_labels[i], polygon)
        # return p

    def _layout_one_adjacent(self, label, edge=None):
        if edge is None:
            for e in range(len(self._surface.polygon(label).edges())):
                opposite = self._surface.opposite_edge(label, e)
                if opposite is not None and opposite[0] in self._transformation:
                    edge = e
                    break

        reference = None
        if edge is not None:
            opposite = self._surface.opposite_edge(label, edge)
            if opposite is None:
                raise ValueError("cannot make polygons adjacent along edge that is not glued")
            reference = opposite[0]

        if reference is None:
            self._layout_one_adjacent_root(label)
            return

        self._transformation[label] = self._transformation[reference] * self._surface.edge_transformation(label, edge)

    def _plot_polygons(self):
        from sage.plot.graphics import Graphics

        g = Graphics()

        for label in self._transformation:
            if self.is_polygon_visible(label):
                g += self._plot_polygon(label)

        return g

    def _plot_zero_flags(self):
        from sage.plot.graphics import Graphics

        g = Graphics()

        for label in self._transformation:
            if self.is_zero_flag_visible(label):
                g += self._plot_zero_flag(label)

        return g

    def _plot_polygon_labels(self):
        from sage.plot.graphics import Graphics

        g = Graphics()

        for label in self._transformation:
            if self.is_polygon_label_visible(label):
                g += self._plot_polygon_label(label)

        return g

    def _plot_edges(self):
        from sage.plot.graphics import Graphics

        g = Graphics()

        for label in self._transformation:
            for edge in range(len(self._surface.polygon(label).edges())):
                if self.is_edge_visible(label, edge):
                    g += self._plot_edge(label, edge)

        # # Keep track of plotted edges so we do not plot glued edges twice.
        # plotted = set()

        # for label in self._surface.labels():
        #     for edge in range(len(self._surface.polygon(label).edges())):

        #         if not options:
        #             continue

        #         plotted.add((label, edge))

        #         g += self._plot_edge(label, edge, **options)

        return g

    def _plot_edge_labels(self):
        from sage.plot.graphics import Graphics

        g = Graphics()

        for label in self._transformation:
            for edge in range(len(self._surface.polygon(label).edges())):
                if self.is_edge_label_visible(label, edge):
                    g += self._plot_edge_label(label, edge)

        return g

    def graphical_points(self, point):
        self._ensure_layout()

        for label, p in point.representatives():
            if not self.is_visible(label):
                continue

            yield self._graphical_point(p, self._transformation[label])

    def graphical_point(self, label, point):
        self._ensure_layout()
        return self._graphical_point(point, self._transformation[label])

    def plot_point(self, point, **kwargs):
        self._ensure_layout()

        from sage.plot.graphics import Graphics

        g = Graphics()

        for p in self.graphical_points(point):
            g += self._plot_point(p, **{**(self._point_options(point) or {}), **kwargs})

        return g

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
        r"""
        Converts from graphical coordinates to similarity surface coordinates.

        If a vector v is provided then a
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
            Traceback (most recent call last):
            ...
            ValueError: not in a visible graphical polygon
            sage: gs.to_surface((1, 1/2))
            Point (1, 1/2) of polygon 1
            sage: gs.to_surface((1,-2), v=(1,0))
            Traceback (most recent call last):
            ...
            ValueError: not in a visible graphical polygon
            sage: gs.to_surface((1, 1/2), v=(1, 0))
            SimilaritySurfaceTangentVector in polygon 1 based at (1, 1/2) with vector (1, 0)

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.infinite_staircase()
            sage: gs = s.graphical_surface()
            sage: gs.to_surface((1, 1/2), (1, 1))
            SimilaritySurfaceTangentVector in polygon -1 based at (0, 1/2) with vector (1, 1)

            sage: s = translation_surfaces.square_torus()
            sage: pc = s.minimal_cover(cover_type="planar")
            sage: gs = pc.graphical_surface()
            sage: gs.layout(limit=20)
            sage: gs.to_surface((3,2))
            Vertex 0 of polygon (0, (x, y) |-> (x + 3, y + 2))

        """
        import warnings
        if v is not None:
            warnings.warn("to_surface() has been deprecated and will be removed in a future version of sage-flatsurf; use tangent_vectors() instead.")
        else:
            warnings.warn("to_surface() has been deprecated and will be removed in a future version of sage-flatsurf; use points() instead.")

        if search_all:
            raise NotImplementedError("search_all is not supported as a keyword argument to to_surface anymore")

        if ring is not None and ring is not self._surface.base_ring():
            raise NotImplementedError("ring is not supported as a keyword argument to to_surface() anymore; use change_ring() on the underlying surface instead.")

        labels = None
        if label is not None:
            labels = [label]

        if v is not None:
            items = self.tangent_vectors(point, v, labels=labels)
        else:
            items = self.points(point, labels=labels)

        if return_all:
            return items

        for item in items:
            return item

        raise ValueError("not in a visible graphical polygon")

    def points(self, coordinates, labels=None):
        self._ensure_layout()

        points = set()

        for label in labels or self._transformation:
            point = self._point(label, coordinates)
            if point is None:
                continue
            if point in points:
                continue
            points.add(point)
            yield point

    def tangent_vectors(self, coordinates, direction, labels=None):
        self._ensure_layout()

        tangent_vectors = set()

        for label in labels or self._transformation:
            tangent_vector = self._tangent_vector(label, coordinates, direction)
            if tangent_vector is None:
                continue
            if tangent_vector in tangent_vectors:
                continue
            tangent_vectors.add(tangent_vector)
            yield tangent_vector


class GraphicalSimilaritySurface(GraphicalSurface):
    r"""
    Renders a similarity surface.

    TODO
    """

    def _layout_one_adjacent_root(self, label):
        # We cannot place this polygon adjacent to an existing one since
        # none of its neighbors has been placed yet. Therefore, we just
        # leave it at its natural position (we should improve this
        # eventually so that connected components do not overlap.)
        from flatsurf.geometry.similarity import SimilarityGroup

        T = SimilarityGroup(self._surface.base_ring())
        self._transformation[label] = T((0, 0))

    def graphical_polygon(self, label):
        self._ensure_layout()
        return self._surface.polygon(label).apply_similarity(self._transformation[label])

    def _graphical_vertices(self, label, edge):
        vertices = self.graphical_polygon(label).vertices()
        return vertices[edge], vertices[(edge + 1) % len(vertices)]

    def _plot_polygon(self, label):
        from sage.all import polygon2d, RDF
        return polygon2d(self.graphical_polygon(label).change_ring(RDF).vertices(), **self._polygon_options(label))

    def _plot_zero_flag(self, label):
        line_options = self._zero_flag_options(label)

        t = .5
        if "t" in line_options:
            t = float(line_options.pop("t"))

        polygon = self.graphical_polygon(label)
        start = polygon.vertices()[0]
        centroid = polygon.centroid()

        from sage.all import line2d
        return line2d((start, start + t * (centroid - start)), **line_options)

    def _plot_polygon_label(self, label):
        from sage.all import text
        return text(str(label), self.graphical_polygon(label).centroid(), **self._polygon_label_options(label))

    def _plot_edge(self, label, edge):
        from sage.all import line2d
        return line2d(self._graphical_vertices(label, edge), **self._edge_options(label, edge))

    def _plot_edge_label(self, label, edge):
        from sage.all import RDF

        e = self.graphical_polygon(label).edges()[edge]
        options = dict(self._edge_label_options(label, edge))
        options.pop("label", None)

        if "position" in options:
            if options["position"] not in ["inside", "outside", "edge"]:
                raise ValueError(
                    "The 'position' parameter must take the value 'inside', 'outside', or 'edge'."
                )
            pos = options.pop("position")
        else:
            pos = "inside"

        if pos == "outside":
            # position outside polygon.
            if "horizontal_alignment" in options:
                pass
            elif e[1] > 0:
                options["horizontal_alignment"] = "left"
            elif e[1] < 0:
                options["horizontal_alignment"] = "right"
            else:
                options["horizontal_alignment"] = "center"

            if "vertical_alignment" in options:
                pass
            elif e[0] > 0:
                options["vertical_alignment"] = "top"
            elif e[0] < 0:
                options["vertical_alignment"] = "bottom"
            else:
                options["vertical_alignment"] = "center"

        elif pos == "inside":
            # position inside polygon.
            if "horizontal_alignment" in options:
                pass
            elif e[1] < 0:
                options["horizontal_alignment"] = "left"
            elif e[1] > 0:
                options["horizontal_alignment"] = "right"
            else:
                options["horizontal_alignment"] = "center"

            if "vertical_alignment" in options:
                pass
            elif e[0] < 0:
                options["vertical_alignment"] = "top"
            elif e[0] > 0:
                options["vertical_alignment"] = "bottom"
            else:
                options["vertical_alignment"] = "center"

        else:
            # centered on edge.
            if "horizontal_alignment" in options:
                pass
            else:
                options["horizontal_alignment"] = "center"
            if "vertical_alignment" in options:
                pass
            else:
                options["vertical_alignment"] = "center"

        if "t" in options:
            t = RDF(options.pop("t"))
        else:
            t = 0.3

        if "push_off" in options:
            push_off = RDF(options.pop("push_off"))
        else:
            push_off = 0.03
        if pos == "outside":
            push_off = -push_off
        # Now push_off stores the amount it should be pushed into the polygon

        from sage.all import text

        V = RDF**2
        no = V((-e[1], e[0]))
        return text(self._edge_label(label, edge), self.graphical_polygon(label).vertices()[edge] + t * e + push_off * no, **options)

    def _plot_point(self, p, **kwargs):
        from sage.all import point2d
        return point2d(points=[p], **kwargs)

    def _graphical_point(self, point, similarity):
        return similarity(point)

    def _point(self, label, graphical_coordinates):
        from sage.all import vector
        graphical_coordinates = vector(graphical_coordinates)
        coordinates = (~self._transformation[label])(graphical_coordinates)
        polygon = self._surface.polygon(label)
        if not polygon.contains_point(coordinates):
            return None
        return self._surface.point(label, coordinates)

    def _tangent_vector(self, label, graphical_base, graphical_direction):
        from sage.all import vector
        graphical_base = vector(graphical_base)
        base = (~self._transformation[label])(graphical_base)
        position = self._surface.polygon(label).get_point_position(base)
        if self._point(label, graphical_base) is None:
            return None

        graphical_direction = vector(graphical_direction)
        direction = (~self._transformation[label].derivative()) * graphical_direction

        if position.is_vertex():
            vertex = position.get_vertex()
            polygon = self._surface.polygon(label)

            from flatsurf.geometry.euclidean import ccw
            if ccw(polygon.edge(vertex), direction) < 0:
                return None
            if (
                ccw(
                    polygon.edge((vertex - 1) % len(polygon.vertices())),
                    direction,
                )
                < 0
            ):
                return None

        return self._surface.tangent_vector(
            label,
            base,
            direction)


class GraphicalHyperbolicIsometrySurface(GraphicalSurface):
    r"""
    Renders a surface built from hyperbolic polygons glued by isometries.

    TODO
    """

    def _layout_one_adjacent_root(self, label):
        # We cannot place this polygon adjacent to an existing one since
        # none of its neighbors has been placed yet. Therefore, we just
        # leave it at its natural position (we should improve this
        # eventually so that connected components do not overlap.)
        from sage.all import matrix
        self._transformation[label] = matrix([[1, 0], [0, 1]])

    def graphical_polygon(self, label):
        isometry = self._transformation[label]
        return self._surface.polygon(label).apply_isometry(isometry)

    def graphical_edge(self, label, edge):
        isometry = self._transformation[label]
        return self._surface.polygon(label).edges()[edge].apply_isometry(isometry)

    def _graphical_vertices(self, label, edge):
        return self.graphical_edge(label, edge).vertices()

    def _plot_polygon(self, label):
        return self.graphical_polygon(label).plot(**self._polygon_options(label))

    def _plot_zero_flag(self, label):
        # Instead of plotting a proper "flag" we just try to color edge 0. When
        # the midpoint cannot be computed, the coloring might be misleading.
        # Also, the hyperbolic midpoint is probably not ideal.
        polygon = self.graphical_polygon(label)

        edge = polygon.edges()[0]

        start = edge.start()
        end = edge.end()
        if edge.is_finite():
            try:
                end = edge.midpoint()
            except ValueError:
                # Cannot compute midpoint over this ring.
                pass

        return start.segment(end).plot(**self._zero_flag_options(label))

    def _plot_polygon_label(self, label):
        # TODO
        from sage.all import Graphics
        return Graphics()

    def _plot_edge(self, label, edge):
        return self.graphical_edge(label, edge).plot(**self._edge_options(label, edge))

    def _plot_edge_label(self, label, edge):
        # TODO
        from sage.all import Graphics
        return Graphics()

    def _graphical_point(self, point, similarity):
        return point.apply_isometry(similarity)


class DELETEME:
    """
    This class essentially consists of a collection of GraphicalPolygons which
    control how individual polygons are positioned. In addition, this class
    stores options which are passed to the polygons when they are rendered.

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
        barycenter to the zero vertex of each polygon. Useful in working out
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

    EXAMPLES::

        sage: from flatsurf import similarity_surfaces
        sage: from flatsurf.graphical.surface import GraphicalSurface

        sage: s = similarity_surfaces.example()
        sage: gs = s.graphical_surface(polygon_color="red")
        sage: gs.plot()
        ...Graphics object consisting of 13 graphics primitives

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
        self.polygon_options = {"color": "lightgray"}
        self.will_plot_upside_down_polygons = False
        self.upside_down_polygon_options = {"color": "lightgray", "zorder": -1}
        self.will_plot_polygon_labels = True
        self.polygon_label_options = {
            "color": "black",
            "vertical_alignment": "center",
            "horizontal_alignment": "center",
        }
        self.will_plot_edges = True
        self.will_plot_non_adjacent_edges = True
        self.non_adjacent_edge_options = {"color": "blue"}
        self.will_plot_adjacent_edges = True
        self.adjacent_edge_options = {"color": "blue", "linestyle": ":"}
        self.will_plot_self_glued_edges = True
        self.self_glued_edge_options = {"color": "red"}
        self.will_plot_edge_labels = True
        self.edge_label_options = {"color": "blue"}
        self.will_plot_zero_flags = False
        self.zero_flag_options = {"color": "green", "thickness": 0.5}

        self.process_options(
            adjacencies=adjacencies,
            polygon_labels=polygon_labels,
            edge_labels=edge_labels,
        )

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

    def is_adjacent(self, p, e):
        r"""
        Returns the truth value of the statement
        'The polygon opposite edge (p,e) is adjacent to that edge.'

        EXAMPLES::

            sage: from flatsurf import similarity_surfaces
            sage: s = similarity_surfaces.example()
            sage: g = s.graphical_surface()
            sage: g.layout(label=s.root())
            sage: g.is_adjacent(0,0)
            False
            sage: g.is_adjacent(0,1)
            False
            sage: g.layout()
            sage: g.is_adjacent(0, 0)
            False
            sage: g.is_adjacent(0,1)
            False
            sage: g.is_adjacent(0,2)
            True

        """
        pp, ee = self.opposite_edge(p, e)
        if not self.is_visible(pp):
            return False
        g = self.graphical_polygon(p)
        gg = self.graphical_polygon(pp)
        return g.transformed_vertex(e) == gg.transformed_vertex(
            ee + 1
        ) and g.transformed_vertex(e + 1) == gg.transformed_vertex(ee)

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
