r"""
Plotting primitives for subsets of the hyperbolic plane

EXAMPLES:

Usually, the primitives defined here should not be used directly. Instead the
:meth:`flatsurf.geometry.hyperbolic.HyperbolicConvexSet.plot` method of
hyperbolic sets internally uses these primitives::

    sage: from flatsurf import HyperbolicPlane

    sage: H = HyperbolicPlane()

    sage: geodesic = H.vertical(0)

    sage: plot = geodesic.plot()

    sage: list(plot)
    [CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)), CartesianPathPlotCommand(code='RAYTO', args=(0, 1))])]

.. NOTE::

    The need for these primitives arises because SageMath has no good
    facilities to plot infinite objects such as lines and rays. However, these are
    needed to plot subsets of the hyperbolic plane in the upper half plane model.

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2023 Julian Rüth
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

from dataclasses import dataclass

from sage.plot.primitive import GraphicPrimitive
from sage.misc.cachefunc import cached_method
from sage.misc.decorators import options, rename_keyword

from flatsurf.geometry.hyperbolic import HyperbolicPoint, HyperbolicPlane


class CartesianPathPlot(GraphicPrimitive):
    r"""
    A plotted path in the hyperbolic plane, i.e., a sequence of commands and
    associated control points in the hyperbolic plane.

    The ``plot`` methods of most hyperbolic convex sets rely on such a path.
    Usually, such a path should not be produced directly.

    This can be considered a more generic version of ``sage.plot.line.Line``
    and ``sage.plot.polygon.Polygon`` since this is not limited to finite line
    segments. At the same time this generalizes matplotlib's ``Path`` somewhat,
    again by allowing infinite rays and lines.

    INPUT:

    - ``commands`` -- a sequence of :class:`HyperbolicPathPlotCommand`
      describing the path.

    - ``options`` -- a dict or ``None`` (the default), options to affect the
      plotting of the path; the options accepted are the same that ``Polygon``
      of :mod:`sage.plot.polygon` accepts.

    EXAMPLES:

    A geodesic plot as a single such path (wrapped in a SageMath graphics
    object)::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand
        sage: H = HyperbolicPlane()
        sage: P = H.vertical(0).plot()
        sage: isinstance(P[0], CartesianPathPlot)
        True

    The sequence of commands should always start with a move command to
    establish the starting point of the plot; note that coordinates are always
    given in the Cartesian two dimension plot coordinate system::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (0, 0))
        ....: ])

    After the initial move, a sequence of arcs can be drawn to represent
    objects in the upper half plane model; the parameters is the center and the
    end point of the arc::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
        ....:     CartesianPathPlotCommand("ARCTO", ((0, 0), (0, 1))),
        ....: ])

    We can also draw segments to represent objects in the Klein disk model::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
        ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
        ....: ])

    Additionally, we can draw rays to represent verticals in the upper half
    plane model; the parameter is the direction of the ray, i.e., (0, 1) for
    a vertical::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (0, 0)),
        ....:     CartesianPathPlotCommand("RAYTO", (0, 1)),
        ....: ])

    Similarly, we can also move the cursor to an infinitely far point in a
    certain direction. This can be used to plot a half plane, e.g., the point
    with non-negative real part::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (0, 0)),
        ....:     CartesianPathPlotCommand("MOVETOINFINITY", (1, 0)),
        ....:     CartesianPathPlotCommand("MOVETOINFINITY", (0, 1)),
        ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
        ....: ])

    In a similar way, we can also draw an actual line, here the real axis::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETOINFINITY", (-1, 0)),
        ....:     CartesianPathPlotCommand("RAYTO", (1, 0)),
        ....: ])

    Finally, we can draw an arc in clockwise direction; here we plot the point
    in the upper half plane of norm between 1 and 2::

        sage: P = CartesianPathPlot([
        ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
        ....:     CartesianPathPlotCommand("ARCTO", ((0, 0), (1, 0))),
        ....:     CartesianPathPlotCommand("MOVETO", (2, 0)),
        ....:     CartesianPathPlotCommand("RARCTO", ((0, 0), (-2, 0))),
        ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
        ....: ])

    .. SEEALSO::

        :func:`hyperbolic_path` to create a ``Graphics`` containing a
        :class:`CartesianPathPlot`, most likely you want to use that
        function if you want to use this functionality in plots of your own.

    """

    def __init__(self, commands, options=None):
        options = options or {}

        valid_options = self._allowed_options()
        for option in options:
            if option not in valid_options:
                raise RuntimeError(f"option {option} not valid")

        # We don't validate the commands here. The consumers below are going to
        # do that implicitly.
        self._commands = commands

        super().__init__(options)

    def _allowed_options(self):
        r"""
        Return the options that are supported by a path.

        We support all the options that are understood by a SageMath polygon.

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot
            sage: P = CartesianPathPlot([])
            sage: P._allowed_options()  # random output depending on the version of SageMath
            {'alpha': 'How transparent the figure is.',
             'edgecolor': 'The color for the border of filled polygons.',
             'fill': 'Whether or not to fill the polygon.',
             'hue': 'The color given as a hue.',
             'legend_color': 'The color of the legend text.',
             'legend_label': 'The label for this item in the legend.',
             'linestyle': 'The style of the enclosing line.',
             'rgbcolor': 'The color as an RGB tuple.',
             'thickness': 'How thick the border line is.',
             'zorder': 'The layer level in which to draw'}
            sage: P = CartesianPathPlot([], options={"alpha": .1})
            sage: P = CartesianPathPlot([], options={"beta": .1})
            Traceback (most recent call last):
            ...
            RuntimeError: option beta not valid

        """
        from sage.plot.polygon import Polygon

        return Polygon([], [], {})._allowed_options()

    def __repr__(self):
        r"""
        Return a printable representation of this plot for debugging purposes.

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
            ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
            ....: ])
            sage: P
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(-1, 0)), CartesianPathPlotCommand(code='LINETO', args=(0, 0))])

        """
        return f"CartesianPathPlot({self._commands})"

    def _render_on_subplot(self, subplot):
        r"""
        Render this path on ``subplot``.

        Matplotlib was not really made to draw things that extend to infinity.
        The trick here is to register a callback that redraws whenever the
        viewbox of the plot changes, e.g., as more objects are added to the
        plot or as the plot is dragged around.

        This implements the interface required by SageMath's
        ``GraphicPrimitive``.

        INPUT:

        - ``subplot`` -- the axes of a subplot

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
            ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
            ....: ])

            sage: from matplotlib.figure import Figure
            sage: figure = Figure()
            sage: subplot = figure.add_subplot(111)

            sage: P._render_on_subplot(subplot)
            sage: subplot.patches
            <Axes.ArtistList of 1 patches>

        """
        matplotlib_options = {
            key: value
            for (key, value) in self.options().items()
            if key
            not in {
                "alpha",
                "legend_color",
                "legend_label",
                "linestyle",
                "rgbcolor",
                "thickness",
            }
        }

        from matplotlib.path import Path

        path = Path([(0, 0)])

        from matplotlib.patches import PathPatch

        patch = PathPatch(path, **matplotlib_options)

        subplot.axes.add_patch(patch)

        options = self.options()

        fill = options.pop("fill", None)
        if fill:
            patch.set_fill(True)

        # Translate SageMath options to matplotlib style.
        if "thickness" in options:
            patch.set_linewidth(float(options.pop("thickness")))

        if "linestyle" in options:
            patch.set_linestyle(options.pop("linestyle"))

        if "alpha" in options:
            patch.set_alpha(float(options["alpha"]))

        from sage.plot.colors import to_mpl_color

        color = None
        if "rgbcolor" in options:
            color = to_mpl_color(options.pop("rgbcolor"))

        edge_color = options.pop("edgecolor", None)
        if edge_color is not None:
            edge_color = to_mpl_color(edge_color)

        if edge_color is None:
            if color is None:
                pass
            else:
                patch.set_color(color)
        else:
            patch.set_edgecolor(edge_color)
            if color is None:
                pass
            else:
                patch.set_facecolor(color)

        if "legend_label" in options:
            patch.set_label(options.pop("legend_label"))

        def redraw(_=None):
            r"""
            Redraw after the viewport has been rescaled to make sure that
            infinite rays reach the end of the viewport.
            """
            # We use ._path directly since .set_path is not available in old
            # versions of matplotlib
            patch._path = self._create_path(
                subplot.axes.get_xlim(), subplot.axes.get_ylim(), fill=fill
            )

        subplot.axes.callbacks.connect("ylim_changed", redraw)
        subplot.axes.callbacks.connect("xlim_changed", redraw)
        redraw()

    def _create_path(self, xlim, ylim, fill):
        r"""
        Create a matplotlib path for this primitive in the bounding box given
        by ``xlim`` and ``ylim``.

        This is a helper method for :meth:`_rendor_on_subplot`.

        INPUT:

        - ``xlim`` -- a pair of floats, the lower and upper horizontal bound of
          the view box

        - ``ylim`` -- a pair of floats, the lower and upper vertical bound of
          the view box

        - ``fill`` -- a boolean; whether the area enclosed by this path should
          be filled

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
            ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
            ....: ])

            sage: P._create_path([-1, 1], [-1, 1], fill=False)
            Path(array([[-1.,  0.],
                   [ 0.,  0.]]), array([1, 2], dtype=uint8))

        """
        # Handle the first command. This controls how the path starts.
        command = self._commands[0]

        if command.code == "MOVETO":
            pos = command.args
            direction = None
            vertices = [pos]
        elif command.code == "MOVETOINFINITY":
            # We cannot draw lines yet. We don't seem to need them for hyperbolic plots at the moment.
            raise NotImplementedError(
                "cannot draw a path that starts at an infinite point yet"
            )
        else:
            raise RuntimeError(f"path must not start with a {command.code} command")

        from matplotlib.path import Path

        codes = [Path.MOVETO]

        for command in self._commands[1:]:
            pos, direction = self._extend_path(
                vertices, codes, pos, direction, command, fill, xlim, ylim
            )

        return Path(vertices, codes)

    @staticmethod
    def _infinity(pos, direction, xlim, ylim):
        r"""
        Return the finite coordinates of a point that ressembles infinity
        starting from ``pos`` and going in ``direction``.

        This is a helper method for :meth:`_create_path`.

        INPUT:

        - ``pos`` -- a coordinate in the plane as a tuple

        - ``direction`` -- a direction in the plane as a tuple

        - ``xlim`` -- a pair of floats, the lower and upper horizontal bound of
          the view box

        - ``ylim`` -- a pair of floats, the lower and upper vertical bound of
          the view box

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot

            sage: CartesianPathPlot._infinity((0., 0.), (0., 1.), (-1024., 1024.), (-1024., 1024.))
            (0.000000000000000, 5120.00000000000)

            sage: CartesianPathPlot._infinity((0., 0.), (1., 1.), (-1024., 1024.), (-1024., 1024.))
            (3920.30937574010, 3920.30937574010)

            sage: CartesianPathPlot._infinity((1., 0.), (1., 1.), (-1024., 1024.), (-1024., 1024.))
            (3920.30937574010, 3919.30937574010)

        """
        from sage.all import vector

        direction = vector(direction)
        pos = vector(pos)

        from sage.all import infinity

        if direction[0]:
            λx = max(
                (xlim[0] - pos[0]) / direction[0], (xlim[1] - pos[0]) / direction[0]
            )
        else:
            λx = infinity

        if direction[1]:
            λy = max(
                (ylim[0] - pos[1]) / direction[1], (ylim[1] - pos[1]) / direction[1]
            )
        else:
            λy = infinity

        λ = min(λx, λy)

        # Additionally, we now move out a full plot size so we are sure
        # that no artifacts coming from any sweeps (see below) show up in
        # the final plot.
        plot_size = (xlim[1] - xlim[0]) + (ylim[1] - ylim[0])
        λ += plot_size / direction.norm()

        return pos + λ * direction

    @staticmethod
    def _extend_path(vertices, codes, pos, direction, command, fill, xlim, ylim):
        r"""
        Extend the matplotlib Path ``vertices`` and ``codes`` by realizing ``command``.

        INPUT:

        - ``vertices`` -- a list of control points

        - ``codes`` -- a list of matplotlib Path plot command

        - ``pos`` -- the current location of the plotted path before executing
          ``command``

        - ``direction`` -- ``None`` or the direction if the current position is
          infinite

        - ``command`` -- a :class:`CartesianPathPlotCommand`

        - ``fill`` -- a boolean, whether the interior of the plotted path
          should be filled

        - ``xlim`` -- a pair of floats, the lower and upper horizontal bound of
          the view box

        - ``ylim`` -- a pair of floats, the lower and upper vertical bound of
          the view box

        OUTPUT:

        The new ``pos`` and ``direction`` after executing the ``command``.

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand
            sage: from matplotlib.path import Path

            sage: vertices = [(0., 0.)]
            sage: codes = [Path.MOVETO]
            sage: pos = (0., 0.)
            sage: direction = None
            sage: command = CartesianPathPlotCommand("LINETO", (1., 1.))

            sage: CartesianPathPlot._extend_path(vertices, codes, pos, direction, command, False, (-1024., 1024.), (-1024., 1024.))
            ((1.00000000000000, 1.00000000000000), None)

            sage: vertices
            [(0.000000000000000, 0.000000000000000), (1.00000000000000, 1.00000000000000)]
            sage: codes
            [1, 2]

        """
        from matplotlib.path import Path

        def extend(path):
            vertices.extend(path.vertices[1:])
            codes.extend(path.codes[1:])

        if command.code == "LINETO":
            target = command.args

            if direction is not None:
                vertices.append(
                    CartesianPathPlot._infinity(target, direction, xlim, ylim)
                )
                codes.append(Path.LINETO)
                direction = None

            pos = target

            vertices.append(pos)
            codes.append(Path.LINETO)
        elif command.code == "RAYTO":
            if direction is None:
                direction = command.args

                vertices.append(CartesianPathPlot._infinity(pos, direction, xlim, ylim))
                codes.append(Path.LINETO)
            else:
                start = CartesianPathPlot._infinity(pos, direction, xlim, ylim)

                direction = command.args
                end = CartesianPathPlot._infinity(pos, direction, xlim, ylim)

                # Sweep the bounding box counterclockwise from start to end
                from sage.all import vector

                center = vector(((start[0] + end[0]) / 2, (start[1] + end[1]) / 2))

                extend(CartesianPathPlot._arc_path(center, start, end))

                vertices.append(end)
                codes.append(Path.LINETO if fill else Path.MOVETO)
        elif command.code == "ARCTO":
            target, center = command.args

            assert direction is None

            extend(CartesianPathPlot._arc_path(center, pos, target))

            pos = target
        elif command.code == "RARCTO":
            target, center = command.args

            assert direction is None

            extend(CartesianPathPlot._arc_path(center, target, pos, reverse=True))

            pos = target
        elif command.code == "MOVETO":
            target = command.args
            pos = target
            direction = None
            vertices.append(pos)
            codes.append(Path.MOVETO)
        else:
            raise RuntimeError(f"cannot draw {command.code} yet")

        return pos, direction

    @cached_method
    def get_minmax_data(self):
        r"""
        Return a bounding box for this plot.

        This implements the required interface of SageMath's GraphicPrimitive.

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
            ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
            ....: ])

            sage: P.get_minmax_data()
            {'xmax': 0.0, 'xmin': -1.0, 'ymax': 0.0, 'ymin': 0.0}

        ::

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (0, 0)),
            ....:     CartesianPathPlotCommand("RAYTO", (0, 1)),
            ....: ])

            sage: P.get_minmax_data()['ymax']
            1.0

        """
        try:
            from matplotlib.transforms import Bbox

            bbox = Bbox.null()

            for command in self._commands:
                if command.code in ["MOVETO", "LINETO"]:
                    pos = command.args
                    bbox.update_from_data_xy([pos], ignore=False)
                elif command.code in ["ARCTO", "RARCTO"]:
                    target, center = command.args
                    # We simplify the computation of the bounding box here to
                    # speed things up. Since these are hyperbolic arcs, their
                    # center must be at y=0 and the endpoints at y≥0.
                    # If we want to use this for non-hyperbolic objects, we'd
                    # have to use this instead:
                    # bbox = bbox.union(
                    #     [
                    #         bbox,
                    #         self._arc_path(
                    #             center, target, pos, reverse=command.code == "RARCTO"
                    #         ).get_extents(),
                    #     ]
                    # )
                    bbox = bbox.union(
                        [
                            bbox,
                            Bbox.from_bounds(*pos, 0, 0),
                            Bbox.from_bounds(*target, 0, 0),
                        ]
                    )
                    from sage.all import sgn

                    if sgn(pos[0] - center[0]) != sgn(target[0] - center[0]):
                        # The bounding box includes a point that is higher
                        # than the endpoints of the arc.
                        from math import sqrt

                        bbox = bbox.union(
                            [
                                bbox,
                                Bbox.from_bounds(
                                    center[0],
                                    sqrt((pos[0] - center[0]) ** 2 + pos[1] ** 2),
                                    0,
                                    0,
                                ),
                            ]
                        )

                    pos = target
                elif command.code in ["RAYTO", "MOVETOINFINITY"]:
                    # We just move "1" in the direction of the ray so that
                    # infinite rays are visible at all. If everything plotted
                    # is extremely small, then this is not ideal. Maybe we
                    # should make this relative to the entire plot size.
                    from sage.all import vector

                    direction = vector(command.args)
                    moved = vector(pos) + direction
                    bbox = bbox.union(
                        [bbox, Bbox.from_bounds(pos[0], moved[0], pos[1], moved[1])]
                    )
                else:
                    raise NotImplementedError(
                        f"cannot determine bounding box for {command.code} command"
                    )

            from sage.plot.plot import minmax_data

            return minmax_data(bbox.intervalx, bbox.intervaly, dict=True)
        except Exception as e:
            raise RuntimeError(e)

    @staticmethod
    def _arc_path(center, start, end, reverse=False):
        r"""
        Return a matplotlib path approximating an circular arc.

        This is a helper method for :meth:`_extend_path`.

        INPUT:

        - ``center`` -- the center point of the circle completing the arc

        - ``start`` -- the coordinates of the starting point of the arc

        - ``end`` -- the coordinates of the end point of the arc

        - ``reverse`` -- a boolean (default: ``False``); whether the arc should
          go clockwise from ``end`` to ``start`` instead of going
          counterclockwise from ``start`` to ``end``.

        EXAMPLES::

            sage: from flatsurf.graphical.hyperbolic import CartesianPathPlot

            sage: CartesianPathPlot._arc_path((0, 0), (1, 0), (0, 1))
            Path(array([[1.00000000e+00, 0.00000000e+00], ...

        .. NOTE::

            Currently, we are using 32 Bezier segments to approximate an arc.
            This is probably too many for small arcs and not enough for very
            big arcs and should be optimized.

        """
        from matplotlib.path import Path
        from math import atan2, pi
        from sage.all import vector

        unit_arc = Path.arc(
            atan2(start[1] - center[1], start[0] - center[0]) / pi * 180,
            atan2(end[1] - center[1], end[0] - center[0]) / pi * 180,
            n=32,
        )

        # Scale and translate the arc
        arc_vertices = (
            unit_arc.vertices
            * vector((start[0] - center[0], start[1] - center[1])).norm()
            + center
        )

        if reverse:
            arc_vertices = arc_vertices[::-1]

        return Path(arc_vertices, unit_arc.codes)


@dataclass
class HyperbolicPathPlotCommand:
    r"""
    A step in a hyperbolic plot.

    Such a step is independent of the model chosen for the plot. It merely
    draws a segment (``LINETO``) in the hyperbolic plan or performs a movement
    without drawing a segment (``MOVETO``).

    Each command has a parameter, the point to which the command moves.

    A sequence of such commands cannot be plotted directly. It is first
    converted into a sequence of :class:`CartesianPathPlotCommand` which
    realizes the commands in a specific hyperbolic model.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand
        sage: H = HyperbolicPlane()
        sage: HyperbolicPathPlotCommand("MOVETO", H(0))
        HyperbolicPathPlotCommand(code='MOVETO', target=0)

    """

    code: str  # Literal["MOVETO", "LINETO"] requires Python 3.8
    target: HyperbolicPoint

    def cartesian(self, model, cursor=None, fill=True, stroke=True):
        r"""
        Return the sequence of commands that realizes this plot in the
        Cartesian plot coordinate system.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` or ``"klein"``

        - ``cursor`` -- a point in the hyperbolic plane or ``None`` (the
          default); assume that the cursor has been positioned at ``cursor``
          before this command.

        - ``fill`` -- a boolean; whether to return commands that produce the
          correct polygon to represent the area of the polygon.

        - ``stroke`` -- a boolean; whether to return commands that produce the
          correct polygon to represent the lines of the polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand
            sage: H = HyperbolicPlane()
            sage: command = HyperbolicPathPlotCommand("MOVETO", H(0))
            sage: command.cartesian("half_plane")
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000))]
            sage: command.cartesian("klein")
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, -1.00000000000000))]

            sage: command = HyperbolicPathPlotCommand("LINETO", H(1))
            sage: command.cartesian("half_plane", cursor=H(0))
            [CartesianPathPlotCommand(code='RARCTO', args=((1.00000000000000, 0.000000000000000), (0.500000000000000, 0)))]
            sage: command.cartesian("klein", cursor=H(0))
            [CartesianPathPlotCommand(code='LINETO', args=(1.00000000000000, 0.000000000000000))]

            sage: command = HyperbolicPathPlotCommand("LINETO", H(oo))
            sage: command.cartesian("half_plane", cursor=H(1))
            [CartesianPathPlotCommand(code='RAYTO', args=(0, 1))]
            sage: command.cartesian("klein", cursor=H(1))
            [CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000))]

        """
        if cursor is None:
            if self.code != "MOVETO":
                raise ValueError(
                    "when no previous cursor position is specified, command must be MOVETO"
                )

            if model == "half_plane" and self.target == self.target.parent().infinity():
                return [CartesianPathPlotCommand("MOVETOINFINITY", (0, 1))]

            from sage.all import RR

            return [
                CartesianPathPlotCommand(
                    "MOVETO", self.target.change_ring(RR).coordinates(model=model)
                )
            ]

        if self.code == "LINETO":
            return HyperbolicPathPlotCommand.create_segment_cartesian(
                cursor, self.target, model=model
            )

        if self.code == "MOVETO":
            return HyperbolicPathPlotCommand.create_move_cartesian(
                cursor, self.target, model=model, fill=fill, stroke=stroke
            )

        raise NotImplementedError(
            "cannot convert this command to a Cartesian plot command yet"
        )

    @staticmethod
    def make_cartesian(commands, model, fill=True, stroke=True):
        r"""
        Return the sequence of :class:`CartesianPathPlotCommand` that realizes
        the hyperbolic ``commands`` in the ``model``.

        INPUT:

        - ``commands`` -- a sequence of :class:`HyperbolicPathPlotCommand`.

        - ``model`` -- one of ``"half_plane"`` or ``"klein"``

        - ``fill`` -- a boolean; whether to return commands that produce the
          correct polygon to represent the area of the polygon.

        - ``stroke`` -- a boolean; whether to return commands that produce the
          correct polygon to represent the lines of the polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand
            sage: H = HyperbolicPlane()

        A finite closed triangle in the hyperbolic plane::

            sage: commands = [
            ....:     HyperbolicPathPlotCommand("MOVETO", H(I)),
            ....:     HyperbolicPathPlotCommand("LINETO", H(I + 1)),
            ....:     HyperbolicPathPlotCommand("LINETO", H(2 * I)),
            ....:     HyperbolicPathPlotCommand("LINETO", H(I)),
            ....: ]

        And its corresponding plot in different models::

            sage: HyperbolicPathPlotCommand.make_cartesian(commands, model="half_plane")
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 1.00000000000000)),
             CartesianPathPlotCommand(code='RARCTO', args=((1.00000000000000, 1.00000000000000), (0.500000000000000, 0))),
             CartesianPathPlotCommand(code='ARCTO', args=((0.000000000000000, 2.00000000000000), (-1.00000000000000, 0))),
             CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000))]

            sage: HyperbolicPathPlotCommand.make_cartesian(commands, model="klein")
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
             CartesianPathPlotCommand(code='LINETO', args=(0.666666666666667, 0.333333333333333)),
             CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.600000000000000)),
             CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))]

        Asking for a polygon that works for both fill and stroke is not always
        possible::

            sage: commands = [
            ....:     HyperbolicPathPlotCommand("MOVETO", H(0)),
            ....:     HyperbolicPathPlotCommand("MOVETO", H(1)),
            ....:     HyperbolicPathPlotCommand("LINETO", H(oo)),
            ....:     HyperbolicPathPlotCommand("LINETO", H(0)),
            ....: ]

            sage: HyperbolicPathPlotCommand.make_cartesian(commands, model="half_plane")
            Traceback (most recent call last):
            ...
            ValueError: exactly one of fill & stroke must be set

            sage: HyperbolicPathPlotCommand.make_cartesian(commands, model="half_plane", fill=False, stroke=True)
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
             CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
             CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
             CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))]

            sage: HyperbolicPathPlotCommand.make_cartesian(commands, model="half_plane", fill=True, stroke=False)
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
             CartesianPathPlotCommand(code='LINETO', args=(1.00000000000000, 0.000000000000000)),
             CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
             CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))]

        """
        cartesian_commands = []
        cursor = None

        for command in commands:
            cartesian_commands.extend(
                command.cartesian(model=model, cursor=cursor, stroke=stroke, fill=fill)
            )
            cursor = command.target

        while cartesian_commands and cartesian_commands[-1].code.startswith("MOVE"):
            cartesian_commands.pop()

        return cartesian_commands

    @staticmethod
    def create_segment_cartesian(start, end, model):
        r"""
        Return a sequence of :class:`CartesianPathPlotCommand` that represent
        the closed boundary of a
        :class:`~flatsurf.geometry.hyperbolic.HyperbolicConvexPolygon`, namely
        the segment to ``end`` (from the previous position ``start``).

        This is a helper function for :meth:`cartesian`.

        INPUT:

        - ``start`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPoint`

        - ``end`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPoint`

        - ``model`` -- one of ``"half_plane"`` or ``"klein"`` in which model to
          realize this segment

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand
            sage: H = HyperbolicPlane()

        A finite segment in the hyperbolic plane; note that we assume that
        "cursor" is at ``start``, so only the command that goes to ``end`` is
        returned::

            sage: HyperbolicPathPlotCommand.create_segment_cartesian(H(I), H(2*I), model="half_plane")
            [CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 2.00000000000000))]

        An infinite segment::

            sage: HyperbolicPathPlotCommand.create_segment_cartesian(H(I), H(oo), model="half_plane")
            [CartesianPathPlotCommand(code='RAYTO', args=(0, 1))]

        A segment that is infinite on both ends; it looks the same because the
        starting point is not rendered here::

            sage: HyperbolicPathPlotCommand.create_segment_cartesian(H(0), H(oo), model="half_plane")
            [CartesianPathPlotCommand(code='RAYTO', args=(0, 1))]

        Note that this is a "closed" boundary of the polygon that is left of
        that segment unlike the "open" version produced by
        :meth:`create_move_cartesian` which contains the entire positive real
        axis::

            sage: HyperbolicPathPlotCommand.create_move_cartesian(H(0), H(oo), model="half_plane", stroke=True, fill=False)
            [CartesianPathPlotCommand(code='MOVETOINFINITY', args=(0, 1))]
            sage: HyperbolicPathPlotCommand.create_move_cartesian(H(0), H(oo), model="half_plane", stroke=False, fill=True)
            [CartesianPathPlotCommand(code='RAYTO', args=(1, 0)),
             CartesianPathPlotCommand(code='RAYTO', args=(0, 1))]

        The corresponding difference in the Klein model::

            sage: HyperbolicPathPlotCommand.create_segment_cartesian(H(0), H(oo), model="klein")
            [CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000))]
            sage: HyperbolicPathPlotCommand.create_move_cartesian(H(0), H(oo), model="klein", stroke=True, fill=False)
            [CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 1.00000000000000))]
            sage: HyperbolicPathPlotCommand.create_move_cartesian(H(0), H(oo), model="klein", stroke=False, fill=True)
            [CartesianPathPlotCommand(code='ARCTO', args=((0.000000000000000, 1.00000000000000), (0, 0)))]

        .. NOTE::

            Sometimes there are problems on a very small scale due to our usage
            of RR internally. We should probably use ball arithmetic to make
            things more robust.

        """
        if start == end:
            raise ValueError(
                f"cannot draw segment from point {start} to itself ({end})"
            )

        if model == "half_plane":
            from sage.all import RR

            if end == end.parent().infinity():
                return [
                    CartesianPathPlotCommand(
                        "RAYTO",
                        (0, 1),
                    )
                ]

            end_x, end_y = end.change_ring(RR).coordinates(model="half_plane")

            if start == start.parent().infinity():
                return [
                    CartesianPathPlotCommand("LINETO", (end_x, end_y)),
                ]

            start_x, start_y = start.change_ring(RR).coordinates(model="half_plane")

            # We should probably be more careful here and not just use a random
            # epsilon.
            if (start_x - end_x).abs() < (start_y - end_y).abs() * 1e-6:
                # This segment is (almost) vertical. We plot it as if it were
                # vertical to avoid numeric issues.
                return [CartesianPathPlotCommand("LINETO", (end_x, end_y))]

            real_hyperbolic_plane = HyperbolicPlane(RR)
            geodesic = real_hyperbolic_plane.geodesic(
                real_hyperbolic_plane.point(start_x, start_y, model="half_plane"),
                real_hyperbolic_plane.point(end_x, end_y, model="half_plane"),
            )
            center = (
                (geodesic.start().coordinates()[0] + geodesic.end().coordinates()[0])
                / 2,
                0,
            )

            return [
                CartesianPathPlotCommand(
                    "RARCTO" if start_x < end_x else "ARCTO", ((end_x, end_y), center)
                )
            ]
        elif model == "klein":
            from sage.all import RR

            return [
                CartesianPathPlotCommand(
                    "LINETO", end.change_ring(RR).coordinates(model="klein")
                )
            ]
        else:
            raise NotImplementedError("cannot draw segment in this model")

    @staticmethod
    def create_move_cartesian(start, end, model, stroke=True, fill=True):
        r"""
        Return a list of :class:`CartesianPathPlotCommand` that represent the
        open "segment" on the boundary of a polygon connecting ``start`` and
        ``end``.

        This is a helper function for :meth:`make_cartesian`.

        INPUT:

        - ``start`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPoint`

        - ``end`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPoint`

        - ``model`` -- one of ``"half_plane"`` or ``"klein"`` in which model to
          realize this segment

        - ``stroke`` -- a boolean (default: ``True``); whether this is part of
          a stroke path that is not filled

        - ``fill`` -- a boolean (default: ``True``); whether this is part of a
          filled path that is not stroked

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand
            sage: H = HyperbolicPlane()

            sage: HyperbolicPathPlotCommand.create_move_cartesian(H(0), H(1), "half_plane", stroke=True, fill=False)
            [CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000))]

        """
        if start == end:
            raise ValueError("cannot move from point to itself")

        if start.is_finite():
            raise ValueError(f"starting point of move must be ideal but was {start}")

        if end.is_finite():
            raise ValueError(f"end of move must be ideal but was {end}")

        if fill == stroke:
            raise ValueError("exactly one of fill & stroke must be set")

        if model == "half_plane":
            from sage.all import RR

            if not fill:
                if end == end.parent().infinity():
                    return [CartesianPathPlotCommand("MOVETOINFINITY", (0, 1))]

                return [
                    CartesianPathPlotCommand(
                        "MOVETO", end.change_ring(RR).coordinates()
                    )
                ]

            if start == start.parent().infinity():
                return [
                    CartesianPathPlotCommand("RAYTO", (-1, 0)),
                    CartesianPathPlotCommand(
                        "LINETO", end.change_ring(RR).coordinates()
                    ),
                ]

            if end == end.parent().infinity():
                return [
                    CartesianPathPlotCommand("RAYTO", (1, 0)),
                    CartesianPathPlotCommand("RAYTO", (0, 1)),
                ]

            if (
                start.change_ring(RR).coordinates()[0]
                < end.change_ring(RR).coordinates()[0]
            ):
                return [
                    CartesianPathPlotCommand(
                        "LINETO", end.change_ring(RR).coordinates()
                    )
                ]
            else:
                return [
                    CartesianPathPlotCommand("RAYTO", (1, 0)),
                    CartesianPathPlotCommand("RAYTO", (0, 1)),
                    CartesianPathPlotCommand("RAYTO", (-1, 0)),
                    CartesianPathPlotCommand(
                        "LINETO", end.change_ring(RR).coordinates()
                    ),
                ]

        elif model == "klein":
            from sage.all import RR

            if fill:
                return [
                    CartesianPathPlotCommand(
                        "ARCTO",
                        (end.change_ring(RR).coordinates(model="klein"), (0, 0)),
                    )
                ]
            else:
                return [
                    CartesianPathPlotCommand(
                        "MOVETO", end.change_ring(RR).coordinates(model="klein")
                    )
                ]
        else:
            raise NotImplementedError("cannot move in this model")


@dataclass
class CartesianPathPlotCommand:
    r"""
    A plot command in the plot coordinate system.

    EXAMPLES:

    Move the cursor to the origin of the coordinate system::

        sage: from flatsurf.graphical.hyperbolic import CartesianPathPlotCommand
        sage: P = CartesianPathPlotCommand("MOVETO", (0, 0))

    Draw a line segment to another point::

        sage: P = CartesianPathPlotCommand("LINETO", (1, 1))

    Draw a ray from the current position in a specific direction::

        sage: P = CartesianPathPlotCommand("RAYTO", (1, 1))

    Move the cursor to a point at infinity in a specific direction::

        sage: P = CartesianPathPlotCommand("MOVETOINFINITY", (0, 1))

    When already at a point at infinity, then this draws a line::

        sage: P = CartesianPathPlotCommand("RAYTO", (0, -1))

    When at a point at infinity, we can also draw a ray to a finite point::

        sage: P = CartesianPathPlotCommand("LINETO", (0, 0))

    Finally, we can draw counterclockwise and clockwise sectors of the circle,
    i.e., arcs by specifying the other endpoint and the center of the circle::

        sage: P = CartesianPathPlotCommand("ARCTO", ((2, 0), (1, 0)))
        sage: P = CartesianPathPlotCommand("RARCTO", ((0, 0), (1, 0)))

    .. SEEALSO::

        :class:`CartesianPathPlot` which draws a sequence of such commands with
        matplotlib.

        :meth:`HyperbolicPathPlotCommand.make_cartesian` to generate a sequence
        of such commands from a sequence of plot commands in the hyperbolic plane.

    """

    code: str  # Literal["MOVETO", "MOVETOINFINITY", "LINETO", "RAYTO", "ARCTO", "RARCTO"] requires Python 3.8
    args: tuple


@rename_keyword(color="rgbcolor")
@options(
    alpha=1,
    rgbcolor=(0, 0, 1),
    edgecolor=None,
    thickness=1,
    legend_label=None,
    legend_color=None,
    aspect_ratio=1.0,
    fill=True,
)
def hyperbolic_path(commands, model="half_plane", **options):
    r"""
    Return a SageMath ``Graphics`` object that represents the hyperbolic path
    encoded by ``commands``.

    INPUT:

    - ``commands`` -- a sequence of :class:`HyperbolicPathPlotCommand`

    - ``model`` -- one of ``"half_plane"`` or ``"klein"``

    Many additional keyword arguments are understood, see
    :class:`CartesianPathPlot` for details.

    EXAMPLES:

    .. jupyter-execute::

        sage: from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand, hyperbolic_path
        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: hyperbolic_path([
        ....:     HyperbolicPathPlotCommand("MOVETO", H(0)),
        ....:     HyperbolicPathPlotCommand("LINETO", H(I + 1)),
        ....:     HyperbolicPathPlotCommand("LINETO", H(oo)),
        ....:     HyperbolicPathPlotCommand("LINETO", H(I - 1)),
        ....:     HyperbolicPathPlotCommand("LINETO", H(0)),
        ....: ])
        ...Graphics object consisting of 2 graphics primitives

    .. SEEALSO::

        :meth:`flatsurf.geometry.hyperbolic.HyperbolicConvexSet.plot`

    """
    if options["thickness"] is None:
        if options["fill"] and options["edgecolor"] is None:
            options["thickness"] = 0
        else:
            options["thickness"] = 1

    from sage.plot.all import Graphics

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))

    try:
        if options.get("fill", None):
            g.add_primitive(
                CartesianPathPlot(
                    HyperbolicPathPlotCommand.make_cartesian(
                        commands, model=model, fill=True, stroke=False
                    ),
                    {**options, "thickness": 0},
                )
            )
        g.add_primitive(
            CartesianPathPlot(
                HyperbolicPathPlotCommand.make_cartesian(
                    commands, model=model, fill=False, stroke=True
                ),
                {**options, "fill": False},
            )
        )
    except Exception as e:
        raise RuntimeError(f"Failed to render hyperbolic path {commands}", e)

    if options["legend_label"]:
        g.legend(True)
        g._legend_colors = [options["legend_color"]]
    return g
