r"""
Two dimensional hyperbolic geometry.

EXAMPLES::

    sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

    sage: H = HyperbolicPlane(QQ)

    TODO: More examples.

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
#                      2022 Sam Freedman
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
######################################################################

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.decorators import options, rename_keyword
from sage.plot.primitive import GraphicPrimitive


class HyperbolicPlane(Parent, UniqueRepresentation):
    r"""
    The hyperbolic plane over a base ring.

    All objects in the plane must be specified over the given base ring. Note
    that, in some representations, objects might appear to live in a larger
    ring. E.g., when specifying a line by giving a center and the square of
    its radius in the half plane model, then the ideal endpoints of this line
    might have coordinates in the ring after adjoining a square root.

    The implemented elements of the plane are convex subsets such as (finite
    and infinite) points, geodesics, closed half planes, and closed convex
    polygons.

    ALGORITHM:

    We do not use a fixed representation of the hyperbolic plane internally but
    switch between the Poincaré half plane and the Klein model freely.

    For the Klein model, we use a unit disc centered at (0, 0). The map from
    the Poincaré half plane sends the imaginary unit `i` to the center at the
    origin, and sends 0 to (0, -1), 1 to (1, 0), -1 to (-1, 0) and infinity to
    (0, 1). The Möbius transformation

    .. MATH::

        z \mapsto \frac{z-i}{1 - iz}

    maps from the half plane model to the Poincaré disc model. We then post-compose this with the map that goes from the Poincaré disc model to the Klein disc model, which in polar coordinates sends

    .. MATH::

        (\phi, r)\mapsto \left(\phi, \frac{2r}{1 + r^2}\right).

    When we write this map out explicitly in Euclidean coordinates, we get

    .. MATH::

        (x, y) \mapsto \frac{1}{1 + x^2 + y^2}\left(2x, -1 + x^2 + y^2\right)

    and

    .. MATH::

        (x, y) \mapsto \frac{1}{1 - y}\left(x, \sqrt{1 - x^2 - y^2}\right),

    for its inverse.

    A geodesic in the Poincaré half plane is then given by an equation of the form

    .. MATH::

        a(x^2 + y^2) + bx + c = 0

    which converts to an equation in the Klein disc as

    .. MATH::

        (a + c) + bx + (a - c)y = 0.

    Note that the intersection of two geodesics with coefficients in a field
    `K` therefore has coordinates in the same field `K` in either model.

    Conversely, a geodesic's equation in the Klein disc

    .. MATH::

        a + bx + cy = 0

    corresponds to the equation

    .. MATH::

        (a + c)(x^2 + y^2) + 2bx + (a - c) = 0

    in the Poincaré half plane model.

    INPUT:

    - ``base_ring`` -- a base ring for the coefficients defining the equations
      of geodesics in the plane; defaults to the rational field if not
      specified.

    - ``category`` -- the category for this object; if not specified, defaults
      to sets with partial maps. Note that we do not use metric spaces here
      since the elements of this space are convex subsets of the hyperbolic
      plane and not just points so the elements do not satisfy the assumptions
      of a metric space.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

        sage: HyperbolicPlane()
        Hyperbolic Plane over Rational Field

    """

    @staticmethod
    def __classcall__(cls, base_ring=None, category=None):
        r"""
        Create the hyperbolic plane with normalized arguments to make it a
        unique SageMath parent.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane() is HyperbolicPlane(QQ)
            True

        """
        from sage.all import QQ
        base_ring = base_ring or QQ

        from sage.categories.all import Sets
        category = category or Sets()

        return super(HyperbolicPlane, cls).__classcall__(cls, base_ring=base_ring, category=category)

    def __init__(self, base_ring, category):
        r"""
        Create the hyperbolic plane over ``base_ring``.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: TestSuite(HyperbolicPlane(QQ)).run()
            sage: TestSuite(HyperbolicPlane(AA)).run()

        """
        super().__init__(category=category)
        self._base_ring = base_ring

    def _an_element_(self):
        r"""
        Return an element of the hyperbolic plane (mostly for testing.)

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().an_element()
            0

        """
        return self.real(0)

    def some_elements(self):
        r"""
        Return some representative convex subsets for automated testing.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().some_elements()
            [{}, ∞, 0, 1, -1]

        """
        # TODO: Return more elements.
        return [self.empty_set(), self.infinity(), self.real(0), self.real(1), self.real(-1)]

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as an element of the hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane(QQ)
            sage: H(H.an_element()) in H
            True

        Base ring elements can be converted to ideal points::

            sage: H(1)
            1

        The point at infinity in the half plane model can be written directly::

            sage: H(oo)
            ∞

        Complex numbers in the upper half plane can be converted to points in
        the hyperbolic plane::

            sage: H(I)
            I

        Elements can be converted between hyperbolic planes with compatible base rings::

            sage: HyperbolicPlane(AA)(H(1))
            1

        """
        if x.parent() is self:
            return x

        from sage.all import Infinity
        if x is Infinity:
            return self.infinity()

        if x in self.base_ring():
            return self.real(x)

        if isinstance(x, HyperbolicConvexSubset):
            return x.change_ring(self.base_ring())

        from sage.categories.all import NumberFields
        if x.parent() in NumberFields():
            K = x.parent()

            from sage.all import I
            if I not in K:
                raise NotImplementedError("cannot create a hyperbolic point from an element in a number field that does not contain the imaginary unit")

            return self.point(x.real(), x.imag(), model="half_plane")

        raise NotImplementedError("cannot create a subset of the hyperbolic plane from this element yet.")

    def base_ring(self):
        r"""
        Return the base ring over which objects in the plane are defined.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().base_ring()
            Rational Field

        """
        return self._base_ring

    def infinity(self):
        r"""
        Return the point at infinity in the Poincaré half plane model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().infinity()
            ∞

        """
        return self.projective(1, 0)

    def real(self, r):
        r"""
        Return the ideal point ``r`` on the real axis in the Poincaré half
        plane model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().real(-2)
            -2

        """
        return self.projective(r, 1)

    def projective(self, p, q):
        r"""
        Return the ideal point with projective coordinates ``[p: q]`` in the
        Poincaré half plane model.

        EXMAPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: H.projective(0, 1)
            0

            sage: H.projective(-1, 0)
            ∞

            sage: H.projective(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: one of p and q must not be zero

        """
        if p == 0 and q == 0:
            raise ValueError("one of p and q must not be zero")

        if q == 0:
            return self.point(0, 1, model="klein")

        p = self.base_ring()(p)
        q = self.base_ring()(q)

        return self.point(p/q, 0, model="half_plane")

    def point(self, x, y, model):
        r"""
        Return the point with coordinates (x, y) in the given model.

        When ``model`` is ``"half_plane"``, return the point `x + iy` in the upper half plane.

        When ``model`` is ``"klein"``, return the point (x, y) in the Klein disc.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.point(0, 1, model="half_plane")
            I

            sage: H.point(1, 2, model="half_plane")
            1 + 2*I

            sage: H.point(0, 1, model="klein")
            ∞

        """
        x = self.base_ring()(x)
        y = self.base_ring()(y)

        if model == "klein":
            return self.__make_element_class__(HyperbolicPoint)(self, x, y)
        if model == "half_plane":
            denominator = 1 + x*x + y*y
            return self.point(
                x=2*x / denominator,
                y=(-1 + x*x + y*y) / denominator,
                model="klein")

        raise NotImplementedError("unsupported model")

    def half_circle(self, center, radius_squared):
        r"""
        Return the geodesic centered around the real ``center`` and with
        ``radius_squared`` in the Poincaré half plane model. The geodesic is
        oriented such that the point at infinity is to its left.

        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_circle(0, 1)
            {(x^2 + y^2) - 1 = 0}

            sage: H.half_circle(1, 3)
            {(x^2 + y^2) - 2*x - 2 = 0}

            sage: H.half_circle(1/3, 1/2)
            {18*(x^2 + y^2) - 12*x - 7 = 0}

        """
        center = self.base_ring()(center)
        radius_squared = self.base_ring()(radius_squared)

        if radius_squared <= 0:
            raise ValueError("radius must be positive")

        # Represent this geodesic as a(x^2 + y^2) + b*x + c = 0
        a = 1
        b = -2*center
        c = center*center - radius_squared

        # Convert to the Klein model.
        return self.__make_element_class__(HyperbolicGeodesic)(self, a + c, b, a - c)

    def vertical(self, real):
        r"""
        Return the vertical geodesic at the ``real`` ideal point in the
        Poincaré half plane model. The geodesic is oriented such that it goes
        from ``real`` to the point at infinity.

        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0)
            {-x = 0}

            sage: H.vertical(1)
            {-x + 1 = 0}

            sage: H.vertical(-1)
            {-x - 1 = 0}

        """
        real = self.base_ring()(real)

        # Convert the equation -x + real = 0 to the Klein model.
        return self.__make_element_class__(HyperbolicGeodesic)(self, real, -1, -real)

    def geodesic(self, a, b):
        r"""
        Return the geodesic from `a` to `b`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geodesic(-1, 1)
            {(x^2 + y^2) - 1 = 0}

            sage: H.geodesic(0, I)
            {-x = 0}

            sage: H.geodesic(-1, I + 1)
            {2*(x^2 + y^2) - x - 3 = 0}

        """
        a = self(a)
        b = self(b)

        C = b._x - a._x
        B = a._y - b._y
        A = -(B * a._x + C * a._y)

        return self.__make_element_class__(HyperbolicGeodesic)(self, A, B, C)

    def half_space(self, geodesic):
        r"""
        Return the closed half plane that is on the left of ``geodesic``.

        Use the ``-`` operator to pass to the half space on the right.
        """
        raise NotImplementedError

    def intersection(self, subsets):
        r"""
        Return the intersection of convex ``subsets``.
        """
        raise NotImplementedError

    def empty_set(self):
        r"""
        Return an empty subset of this space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().empty_set()
            {}

        """
        return self.__make_element_class__(HyperbolicEmptySet)(self)

    def _repr_(self):
        r"""
        Return a printable representation of this hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: HyperbolicPlane(AA)
            Hyperbolic Plane over Algebraic Real Field

        """
        return f"Hyperbolic Plane over {repr(self.base_ring())}"


class HyperbolicConvexSubset(Element):
    r"""
    Base class for convex subsets of :class:`HyperbolicPlane`.
    """

    def _equations(self, model):
        r"""
        Return a minimal set of equations defining a set of half spaces such
        that this set is the intersection of these half spaces.

        The equations are given as triples ``a``, ``b``, ``c`` such that

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is in the half space if `a(x^2 + y^2) + bx + c ≥ 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy ≥ 0`.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.
        """
        # TODO: Check that all subclasses implement this.
        raise NotImplementedError("Convex sets must implement this method.")

    def is_subset(self, other):
        r"""
        Return whether the convex set ``other`` is a subset of this set.
        """
        raise NotImplementedError

    def intersection(self, other):
        r"""
        Return the intersection with the ``other`` convex set.
        """
        return self.parent().intersection([self, other])

    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this set.
        """
        raise NotImplementedError

    def is_finite(self):
        r"""
        Return whether all points in this set are finite.
        """
        raise NotImplementedError

    def change_ring(self, ring):
        r"""
        Return this set as an element of the hyperbolic plane over ``ring``.
        """
        # TODO: Check that all subclasses implement this.
        raise NotImplementedError

    def plot(self, model="half_plane", *kwds):
        r"""
        Return a plot of this subset.
        """
        # TODO: Check that all subclasses implement this.
        raise NotImplementedError

    def apply_isometry(self, isometry):
        r"""
        Return the image of this set under the isometry.

        INPUT:

        - ``isometry`` -- a matrix in `PGL(2,\mathbb{R})`

        """
        # TODO: Understand how isometries transform geodesics so we can
        # transform inequalities in the Klein model.
        raise NotImplementedError

    def _neg_(self):
        r"""
        Return the convex subset obtained by taking the negatives of the half
        spaces whose intersection define this set.
        """
        raise NotImplementedError


class HyperbolicHalfSpace(HyperbolicConvexSubset):
    r"""
    A closed half space of the hyperbolic plane.

    Use :meth:`HyperbolicPlane.half_space` to create a half plane.
    """

    def __init__(self, geodesic):
        self._geodesic = geodesic


class HyperbolicGeodesic(HyperbolicConvexSubset):
    r"""
    An oriented geodesic in the hyperbolic plane.

    Internally, we represent geodesics as the chords satisfying the equation `a
    + bx + cy=0` in the unit disc of the Klein model.

    The geodesic is oriented such that the half space `a + bx + cy ≥ 0` is on
    its left.
    """

    def __init__(self, parent, a, b, c):
        super().__init__(parent)
        self._a = a
        self._b = b
        self._c = c

    def start(self):
        r"""
        Return the ideal starting point of this geodesic.

        Note that this is only possible if the radius of this geodesic is a
        square in the base ring of the :class:`HyperbolicPlane`.
        """
        raise NotImplementedError

    def end(self):
        r"""
        Return the ideal end point of this geodesic.

        Note that this is only possible if the radius of this geodesic is a
        square in the base ring of the :class:`HyperbolicPlane`.
        """
        return (-self).start()

    def _richcmp_(self, other, op):
        r"""
        Return how this geodesic compares to ``other``.

        Geodesics are partially ordered by their slope in the Klein model.
        """
        raise NotImplementedError

    def _neg_(self):
        raise NotImplementedError

    def equation(self, model):
        r"""
        Return an equation for this geodesic as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is on the geodesic if it satisfies `a(x^2 + y^2) + bx + c = 0`.
          The coefficients are such that the half plane `a(x^2 + y^2) + bx + c
          ≥ 0` is on the left of the geodesic.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy = 0`. The sign of the coefficients is such that the half
          plane `a + bx + cy ≥ 0` is on the left of the geodesic.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.
        """
        a, b, c = self._a, self._b, self._c

        if model == "klein":
            return a, b, c

        if model == "half_plane":
            return a + c, 2*b, a - c

        raise NotImplementedError("cannot determine equation for this model yet")

    def _repr_(self):
        # Convert to the Poincaré half plane model as a(x^2 + y^2) + bx + c = 0.
        a, b, c = self.equation(model="half_plane")

        try:
            if self.parent().base_ring().is_exact():
                from sage.all import gcd
                d = gcd((a, b, c))
                a /= d
                b /= d
                c /= d
        except Exception:
            raise
            pass

        from sage.all import PolynomialRing
        R = PolynomialRing(self.parent().base_ring(), names="x")
        if a != 0:
            return "{" + repr(R([0, a]))[:-1] + "(x^2 + y^2)" + repr(R([c, b, 1]))[3:] + " = 0}"
        else:
            return "{" + repr(R([c, b])) + " = 0}"

    def plot(self, model="half_plane", **kwds):
        r"""
        Create a plot of this geodesic in the hyperbolic ``model``.

        Additional arguments are passed on to the underlying SageMath plotting methods.

        EXAMPLES::
        """
        a, b, c = self.equation(model=model)

        if model == "half_plane":
            if a == 0:
                # This is a vertical in the half plane model.
                x = -c/b

                return vertical(x, **kwds)

            else:
                # This is a half-circle in the half plane model.
                center = -(b/a)/2
                radius_squared = center*center - (c/a)

                from sage.plot.all import arc
                from sage.all import RR, pi
                return arc((RR(center), 0), RR(radius_squared).sqrt(), sector=(0, pi), **kwds)

        if model == "klein":
            raise NotImplementedError

        raise NotImplementedError("plotting not supported in this hyperbolic model")


class HyperbolicPoint(HyperbolicConvexSubset):
    r"""
    A (possibly infinite) point in the :class:`HyperbolicPlane`.

    Internally, we represent a point as the Euclidean coordinates in the unit
    disc of the Klein model.
    """

    def __init__(self, parent, x, y):
        super().__init__(parent)
        self._x = x
        self._y = y

    def coordinates(self, model="half_plane", ring=None):
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return projective coordinates in the
        Poincaré half plane model.

        If ``model`` is ``"klein"``, return Euclidean coordinates in the Klein model.

        If no ``ring`` has been specified, an appropriate extension of the base
        ring of the :class:`HyperbolicPlane` is chosen where these coordinates
        live.
        """
        x, y = self._x, self._y

        if model == "half_plane":
            denominator = 1 - y
            return (x / denominator, (1 - x*x - y*y).sqrt()/denominator)

        raise NotImplementedError

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: H.infinity() == H.projective(1, 0)
            True

        """
        from sage.structure.richcmp import rich_to_bool
        return rich_to_bool(op, other is HyperbolicPoint and other._x == self._x and other._y == self._y)

    def _repr_(self):
        if self._x == 0 and self._y == 1:
            return "∞"

        x, y = self.coordinates()

        # We represent x + y*I in R[[I]] so we do not have to reimplement printing ourselves.
        from sage.all import PowerSeriesRing
        return repr(PowerSeriesRing(self.parent().base_ring(), names="I")([x, y]))

    def change_ring(self, ring):
        return HyperbolicPlane(ring).point(self._x, self._y, model="klein")


class HyperbolicConvexPolygon(HyperbolicConvexSubset):
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`HyperbolicHalfSpace`s.
    """

    def __init__(self, parent, half_planes, assume_normalized=False):
        raise NotImplementedError

    def _normalize(self):
        r"""
        Normalize the internal list of half planes so that they describe the
        :meth:`boundary`.
        """
        raise NotImplementedError

    def equations(self):
        r"""
        Return the equations describing the boundary of this polygon.

        The output is minimal and sorted by slope in the Klein model.
        """
        raise NotImplementedError

    def edges(self):
        r"""
        Return the :class:`HyperbolicEdge`s defining this polygon.
        """
        raise NotImplementedError

    def vertices(self):
        r"""
        Return the vertices of this polygon, i.e., the end points of the
        :meth:`edges`.
        """
        raise NotImplementedError


class HyperbolicEdge(HyperbolicConvexSubset):
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.
    """

    def __init__(self, geodesic, start=None, end=None):
        raise NotImplementedError


class HyperbolicEmptySet(HyperbolicConvexSubset):
    r"""
    The empty subset of the hyperbolic plane.
    """

    def __init__(self, parent):
        super().__init__(parent)

    def _richcmp_(self, other, op):
        r"""
        Return how this set compares to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: H.empty_set() == H.empty_set()
            True

        """
        from sage.structure.richcmp import rich_to_bool
        return rich_to_bool(op, other is HyperbolicEmptySet)

    def _repr_(self):
        return "{}"


class Vertical(GraphicPrimitive):
    r"""
    A graphical ray going vertically up from (x, y).

    Used internally to, e.g., plot a vertical geodesic in the upper half plane
    model.

    This object should not be created directly (even inside this module) but by
    calling :meth:`vertical` below.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import Vertical
        sage: Vertical(0, 0)
        Vertical at (0, 0)


    """

    def __init__(self, x, y=0, options={}):
        valid_options = self._allowed_options()
        for option in options:
            if option not in valid_options:
                raise RuntimeError("Error in line(): option '%s' not valid." % option)

        self._x = x
        self._y = y
        super().__init__(options)

    def _allowed_options(self):
        r"""
        Return the options that are supported by a vertical.

        We support all the options that are understood by a SageMath line.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import Vertical
            sage: Vertical(0, 0)._allowed_options()
            {'alpha': 'How transparent the line is.',
             'hue': 'The color given as a hue.',
             'legend_color': 'The color of the legend text.',
             'legend_label': 'The label for this item in the legend.',
             'linestyle': "The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted).",
             'marker': 'the marker symbol (see documentation for line2d for details)',
             'markeredgecolor': 'the color of the marker edge',
             'markeredgewidth': 'the size of the marker edge in points',
             'markerfacecolor': 'the color of the marker face',
             'markersize': 'the size of the marker in points',
             'rgbcolor': 'The color as an RGB tuple.',
             'thickness': 'How thick the line is.',
             'zorder': 'The layer level in which to draw'}

        """
        from sage.plot.line import Line
        return Line([], [], {})._allowed_options()

    def __repr__(self):
        r"""
        Return a printable representation of this graphical primitive.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import Vertical
            sage: Vertical(0, 0)
            Vertical at (0, 0)

        """
        return f"Vertical at ({self._x}, {self._y})"

    def _render_on_subplot(self, subplot):
        r"""
        Render this vertical on the subplot.

        Matplotlib was not really made to draw things that extend to infinity.
        The trick here is to register a callback that redraws the vertical
        whenever the viewbox of the plot changes, e.g., as more objects are
        added to the plot.
        """
        # Rewrite options to only contain matplotlib compatible entries
        matplotlib_options = {
            key: value for (key, value) in self.options().items()
            if key not in {'alpha', 'legend_color', 'legend_label', 'linestyle', 'rgbcolor', 'thickness'}
        }

        from matplotlib.lines import Line2D
        line = Line2D([self._x, self._x], [self._y, self._y], **matplotlib_options)
        subplot.add_line(line)

        # Translate SageMath options to matplotlib style.
        options = self.options()
        line.set_alpha(float(options['alpha']))
        line.set_linewidth(float(options['thickness']))
        from sage.plot.colors import to_mpl_color
        line.set_color(to_mpl_color(options['rgbcolor']))
        line.set_label(options['legend_label'])

        def redraw(_=None):
            r"""
            Redraw the vertical after the viewport has been rescaled to
            make sure it reaches the top of the viewport.
            """
            ylim = max(self._y, subplot.axes.get_ylim()[1])
            line.set_ydata((self._y, ylim))

        subplot.axes.callbacks.connect('ylim_changed', redraw)
        redraw()

    def get_minmax_data(self):
        r"""
        Return the bounding box of this vertical.

        This box is used to make sure that the viewbox of the plot is zoomed
        such that the vertical is visible.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import Vertical
            sage: Vertical(1, 2).get_minmax_data()
            {'xmax': 1, 'xmin': 1, 'ymax': 2, 'ymin': 2}

        """
        from sage.plot.plot import minmax_data
        return minmax_data([self._x, self._x], [self._y, self._y], dict=True)


@rename_keyword(color='rgbcolor')
@options(alpha=1, rgbcolor=(0, 0, 1), thickness=1, legend_label=None, legend_color=None, aspect_ratio='automatic')
def vertical(x, y=0, **options):
    r"""
    Create a SageMath graphics object that describe a vertical ray going up
    from the coordinates ``x``, ``y``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import vertical
        sage: vertical(1, 2)
        Graphics object contsisting of 1 graphics primitive


    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Vertical(x, y, options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    return g
