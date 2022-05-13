r"""
Two dimensional hyperbolic geometry.

EXAMPLES::

    sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

    sage: H = HyperbolicPlane()

Points in the hyperbolic plane can be specified directly with coordinates
in the upper (complex) half plane::

    sage: H(1 + I)
    1 + I

The hyperbolic plane is defined over a fixed base ring; the rationals if no
base has been specified explicitly::

    sage: H(sqrt(2) + I)
    Traceback (most recent call last):
    ...
    TypeError: unable to convert sqrt(2) to a rational

We can use a bigger field instead::

    sage: HAA = HyperbolicPlane(AA)
    sage: HAA(sqrt(2) + I)
    1.414213562373095? + 1.000000000000000?*I

Given two points in the hyperbolic plane, we can form the geodesic they lay on::

    sage: a = H(I)
    sage: b = H(2*I)
    sage: H.geodesic(a, b)
    {-x = 0}

Note that such a geodesic is oriented. The orientation is such that when we
replace the ``=`` in the above representation with a ``≥``, we obtain the half
space on its left::

    sage: H.geodesic(a, b).left_half_space()
    {x ≤ 0}

A vertical can also be specified directly::

    sage: H.vertical(0)
    {-x = 0}

We can also create ideal, i.e., infinite, points in the hyperbolic plane and
construct the geodesic that connects them::

    sage: H(1)
    1

    sage: H(oo)
    ∞

    sage: H.geodesic(1, oo)
    {-x + 1 = 0}

The geodesic that is given by a half circle in the upper half plane can be
created directly by providing its midpoint and the square of its radius::

    sage: H.half_circle(0, 1)
    {(x^2 + y^2) - 1 = 0}

Geodesics can be intersected::

    sage: H.half_circle(0, 1).intersection(H.vertical(0))
    I

    sage: H.half_circle(0, 1).intersection(H.half_circle(0, 2))
    {}

The intersection of two geodesics might be an ideal point::

    sage: H.vertical(-1).intersection(H.vertical(1))
    ∞

General convex subsets of the hyperbolic plane can be constructed by
intersecting half spaces; this way we can construct (possibly unbounded) convex
polygons::

    sage: P = H.intersection(
    ....:   H.vertical(-1).right_half_space(),
    ....:   H.vertical(1).left_half_space(),
    ....:   H.half_circle(0, 2).left_half_space())

    sage: P
    {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

We can also intersect objects that are not half spaces::

    sage: P.intersection(H.vertical(0))
    {x = 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
#                      2022 Sam Freedman
#                      2022 Vincent Delecroix
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

from dataclasses import dataclass

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.decorators import options, rename_keyword
from sage.plot.primitive import GraphicPrimitive


class HyperbolicPlane(Parent, UniqueRepresentation):
    r"""
    The hyperbolic plane.

    All objects in the plane must be specified over the given base ring. Note
    that, in some representations, objects might appear to live in a larger
    ring. E.g., when specifying a line by giving a center and the square of
    its radius in the half plane model, then the ideal endpoints of this line
    might have coordinates in the ring after adjoining a square root.

    The implemented elements of the plane are convex subsets such as (finite
    and infinite) points, geodesics, closed half planes, and closed convex
    polygons.

    ALGORITHM:

    We usually display objects as if they were defined in the Poincaré half
    plane model. However, internally, we store most objects in a representation
    in the Klein model. In that model it tends to be easier to perform
    computations without having to extend the base ring and we can also rely on
    standard algorithms for geometry in the Euclidean plane.

    For the Klein model, we use a unit disk centered at (0, 0). The map from
    the Poincaré half plane sends the imaginary unit `i` to the center at the
    origin, and sends 0 to (0, -1), 1 to (1, 0), -1 to (-1, 0) and infinity to
    (0, 1). The Möbius transformation

    .. MATH::

        z \mapsto \frac{z-i}{1 - iz}

    maps from the half plane model to the Poincaré disk model. We then
    post-compose this with the map that goes from the Poincaré disk model to
    the Klein model, which in polar coordinates sends

    .. MATH::

        (\phi, r)\mapsto \left(\phi, \frac{2r}{1 + r^2}\right).

    When we write out the full map explicitly in Euclidean coordinates, we get

    .. MATH::

        (x, y) \mapsto \frac{1}{1 + x^2 + y^2}\left(2x, -1 + x^2 + y^2\right)

    and

    .. MATH::

        (x, y) \mapsto \frac{1}{1 - y}\left(x, \sqrt{1 - x^2 - y^2}\right),

    for its inverse.

    A geodesic in the Poincaré half plane is given by an equation of the form

    .. MATH::

        a(x^2 + y^2) + bx + c = 0

    which converts to an equation in the Klein model as

    .. MATH::

        (a + c) + bx + (a - c)y = 0.

    Conversely, a geodesic's equation in the Klein model

    .. MATH::

        a + bx + cy = 0

    corresponds to the equation

    .. MATH::

        (a + c)(x^2 + y^2) + 2bx + (a - c) = 0

    in the Poincaré half plane model.

    Note that the intersection of two geodesics defined by coefficients in a
    field `K` in the Klein model has coordinates in `K` in the Klein model.
    The corresponding statement is not true for the Poincaré half plane model.

    INPUT:

    - ``base_ring`` -- a base ring for the coefficients defining the equations
      of geodesics in the plane; defaults to the rational field if not
      specified.

    - ``category`` -- the category for this object; if not specified, defaults
      to sets. Note that we do not use metric spaces here since the elements of
      this space are convex subsets of the hyperbolic plane and not just points
      so the elements do not satisfy the assumptions of a metric space.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

        sage: HyperbolicPlane()
        Hyperbolic Plane over Rational Field

    ::

        sage: HyperbolicPlane(AA)
        Hyperbolic Plane over Algebraic Real Field

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
        # TODO: What should we do about this? We need HyperbolicPlane(RR) internally, but we do not really trust that it's fully functional.
        # if not base_ring.is_exact():
        #     # Much of the implementation might work over inexact rings,
        #     # * we did not really worry about precision issues here so unit
        #     #   tests should be added to check that everything works.
        #     # * if +infinity is in the base ring, then there might be problems
        #     #   in the upper half plane model.
        #     # * if NaN can be represented in the base ring, then there might be
        #     #   problems in many places where we do not expect this to show up.
        #     raise NotImplementedError("hyperbolic plane only implemented over exact rings")

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
            [{}, ∞, 0, 1, -1, ...]

        """
        from sage.all import ZZ

        return [self.empty_set(),
                # Points
                self.infinity(),
                self.real(0),
                self.real(1),
                self.real(-1),
                self.geodesic(0, 2).start(),
                # Geodesics
                self.vertical(1),
                self.half_circle(0, 1),
                self.half_circle(1, 3),
                # Half spaces
                self.vertical(0).left_half_space(),
                self.half_circle(0, 2).left_half_space(),
                # An unbounded polygon
                self.vertical(1).left_half_space().intersection(self.vertical(-1).right_half_space()),
                # An unbounded polygon which is bounded in the Euclidean plane
                self.vertical(1).left_half_space().intersection(self.vertical(-1).right_half_space()).intersection(self.geodesic(0, 1).left_half_space()).intersection(self.geodesic(0, -1).right_half_space()),
                # A bounded polygon
                self.geodesic(-ZZ(1)/3, 2).left_half_space().intersection(self.geodesic(ZZ(1)/3, -2).right_half_space()).intersection(self.geodesic(-ZZ(2)/3, 3).right_half_space()).intersection(self.geodesic(ZZ(2)/3, -3).left_half_space()),
                # An unbounded segment
                self.vertical(0).intersection(self.geodesic(-1, 1).left_half_space()),
                # A bounded segment
                self.vertical(0).intersection(self.geodesic(-2, 2).right_half_space()).intersection(self.geodesic(-ZZ(1)/2, ZZ(1)/2).left_half_space()),
                ]

    def _test_some_subsets(self, tester=None, **options):
        r"""
        Run test suite on some representative convex subsets.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane()._test_some_subsets()

        """
        is_sub_testsuite = tester is not None
        tester = self._tester(tester=tester, **options)

        for x in self.some_elements():
            tester.info(f"\n  Running the test suite of {x}")

            from sage.all import TestSuite
            TestSuite(x).run(verbose=tester._verbose, prefix=tester._prefix + "  ", raise_on_failure=is_sub_testsuite)
            tester.info(tester._prefix + " ", newline=False)

    def random_element(self, kind=None):
        r"""
        Return a random convex subset of this hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane(QQ)

        Make the following randomized tests reproducible::

            sage: set_random_seed(0)

        ::

            sage: H.random_element()
            {}

        Specific types of random subsets can be requested::

            sage: H.random_element("point")
            -1/2 + 1/95*I

            sage: H.random_element("geodesic")
            {-12*(x^2 + y^2) + 1144*x + 1159 = 0}

            sage: H.random_element("half_space")
            {648*(x^2 + y^2) + 1654*x + 85 ≤ 0}

            sage: H.random_element("segment")
            {3*(x^2 + y^2) - 5*x - 1 = 0} ∩ {21*(x^2 + y^2) + 88*x - 89 ≥ 0} ∩ {(x^2 + y^2) + 12*x - 14 ≤ 0}

            sage: H.random_element("polygon")
            {545260*(x^2 + y^2) + 2579907*x + 65638 ≤ 0} ∩ {144*(x^2 + y^2) - 455*x - 2677 ≤ 0} ∩ {18*(x^2 + y^2) + 91*x + 109 ≥ 0}

        """
        kinds = ["empty_set", "point", "geodesic", "half_space", "segment", "polygon"]

        if kind is None:
            from sage.all import randint
            kind = kinds[randint(0, len(kinds) - 1)]

        if kind == "empty_set":
            return self.empty_set()

        if kind == "point":
            return self.point(self.base_ring().random_element(), self.base_ring().random_element().abs(), model="half_plane", check=False)

        if kind == "geodesic":
            a = self.random_element("point")
            b = self.random_element("point")
            while b == a:
                b = self.random_element("point")

            return self.geodesic(a, b)

        if kind == "half_space":
            return self.random_element("geodesic").left_half_space()

        if kind == "segment":
            a = self.random_element("point")
            b = self.random_element("point")
            while a == b or (not a.is_finite() and not b.is_finite()):
                b = self.random_element("point")

            return self.segment(self.geodesic(a, b), start=a, end=b)

        if kind == "polygon":
            from sage.all import ZZ
            interior_points = [self.random_element("point") for i in range(ZZ.random_element().abs() + 3)]

            half_spaces = []

            while len(half_spaces) < len(interior_points):
                half_space = self.random_element("half_space")

                for p in interior_points:
                    if p in half_space:
                        continue

                    a, b, c = half_space.equation(model="klein")

                    x, y = p.coordinates(model="klein")

                    a = -(b*x + c*y)

                    half_space = self.half_space(a, b, c, model="klein")

                    assert p in half_space

                half_spaces.append(half_space)

            return self.polygon(half_spaces)

        raise ValueError(f"kind must be one of {kinds}")

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
        from sage.all import parent
        parent = parent(x)

        if parent is self:
            return x

        from sage.all import Infinity
        if x is Infinity:
            return self.infinity()

        if x in self.base_ring():
            return self.real(x)

        if isinstance(x, HyperbolicConvexSet):
            return x.change_ring(self.base_ring())

        from sage.categories.all import NumberFields
        if parent in NumberFields():
            K = parent

            from sage.all import I
            if I not in K:
                raise NotImplementedError("cannot create a hyperbolic point from an element in a number field that does not contain the imaginary unit")

            return self.point(x.real(), x.imag(), model="half_plane")

        from sage.all import SR
        if parent is SR:
            return self.point(x.real(), x.imag(), model="half_plane")

        from sage.categories.all import Rings
        if parent in Rings():
            raise ValueError(f"cannot convert this element in {parent} to the hyperbolic plane over {self.base_ring()}")

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

        return self.point(p/q, 0, model="half_plane", check=False)

    def point(self, x, y, model, check=True):
        r"""
        Return the point with coordinates (x, y) in the given model.

        When ``model`` is ``"half_plane"``, return the point `x + iy` in the upper half plane.

        When ``model`` is ``"klein"``, return the point (x, y) in the Klein model.

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
        if model is None:
            if isinstance(x, HyperbolicGeodesic):
                if y is not None:
                    raise ValueError("y must be none when x is a geodesic")

                point = self.__make_element_class__(HyperbolicPoint)(self, self(x), None)
            else:
                raise NotImplementedError("unsupported coordinates for implicit model")
        else:
            x = self.base_ring()(x)
            y = self.base_ring()(y)

            if model == "klein":
                point = self.__make_element_class__(HyperbolicPoint)(self, x, y)
            elif model == "half_plane":
                denominator = 1 + x*x + y*y
                return self.point(
                    x=2*x / denominator,
                    y=(-1 + x*x + y*y) / denominator,
                    model="klein",
                    check=check)
            else:
                raise NotImplementedError("unsupported model")

        if check:
            point._check()

        return point

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

        TESTS::

            sage: H.half_circle(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: radius must be positive

            sage: H.half_circle(0, -1)
            Traceback (most recent call last):
            ...
            ValueError: radius must be positive

            sage: H.half_circle(oo, 1)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert +Infinity to a rational

        """
        center = self.base_ring()(center)
        radius_squared = self.base_ring()(radius_squared)

        if radius_squared <= 0:
            raise ValueError("radius must be positive")

        # Represent this geodesic as a(x^2 + y^2) + b*x + c = 0
        a = 1
        b = -2*center
        c = center*center - radius_squared

        return self.geodesic(a, b, c, model="half_plane")

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
        return self.geodesic(real, -1, -real, model="klein")

    def geodesic(self, a, b, c=None, model=None, check=True):
        r"""
        Return a geodesic in the hyperbolic plane.

        If only ``a`` and ``b`` are given, return the geodesic going through the points
        ``a`` and then ``b``.

        If ``c`` is specified and ``model`` is ``"half_plane"``, return the
        geodesic given by the half circle

        .. MATH::

            a(x^2 + y^2) + bx + c = 0

        oriented such that the half plane

        .. MATH::

            a(x^2 + y^2) + bx + c \ge 0

        is to its left.

        If ``c`` is specified and ``model`` is ``"klein"``, return the
        geodesic given by the chord with the equation

        .. MATH::

            a + bx + cy = 0

        oriented such that the half plane

        .. MATH::

            a + bx + cy \ge 0

        is to its left.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geodesic(-1, 1)
            {(x^2 + y^2) - 1 = 0}

            sage: H.geodesic(0, I)
            {-x = 0}

            sage: H.geodesic(-1, I + 1)
            {2*(x^2 + y^2) - x - 3 = 0}

            sage: H.geodesic(2, -1, -3, model="half_plane")
            {2*(x^2 + y^2) - x - 3 = 0}

            sage: H.geodesic(-1, -1, 5, model="klein")
            {2*(x^2 + y^2) - x - 3 = 0}

        TESTS::

            sage: H.geodesic(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: points specifying a geodesic must be distinct

        Geodesics cannot be defined from points whose coordinates are over a
        quadratic field extension::

            sage: H.geodesic(H.half_circle(0, 2).start(), H.half_circle(0, 2).end())
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not a rational number

        """
        if c is None:
            a = self(a)
            b = self(b)

            if a == b:
                raise ValueError("points specifying a geodesic must be distinct")

            ax, ay = a.coordinates(model="klein")
            bx, by = b.coordinates(model="klein")

            C = bx - ax
            B = ay - by
            A = -(B * ax + C * ay)

            return self.geodesic(A, B, C, model="klein", check=check)

        if model is None:
            raise ValueError("a model must be specified when specifying a geodesic with coefficients")

        if model == "half_plane":
            # Convert to the Klein model.
            return self.geodesic(a + c, b, a - c, model="klein", check=check)

        if model == "klein":
            a = self.base_ring()(a)
            b = self.base_ring()(b)
            c = self.base_ring()(c)
            geodesic = self.__make_element_class__(HyperbolicGeodesic)(self, a, b, c)

            if check:
                geodesic = geodesic._normalize()
                geodesic._check()

            return geodesic

        raise NotImplementedError("cannot create geodesic from coefficients in this model")

    def half_space(self, a, b, c, model, check=True):
        r"""
        Return a closed half space from its equation in ``model``.

        If ``model`` is ``"half_plane"``, return the half space

        .. MATH::

            a(x^2 + y^2) + bx + c \ge 0

        in the upper half plane.

        If ``model`` is ``"klein"``, return the half space

        .. MATH::

            a + bx + cy \ge 0

        in the Klein model.

        ..SEEALSO::

            :meth:`HperbolicGeodesic.left_half_space`
            :meth:`HperbolicGeodesic.right_half_space`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_space(0, -1, 0, model="half_plane")
            {x ≤ 0}

        It is often easier to construct a half space as the space bounded by a geodesic::

            sage: H.vertical(0).left_half_space()
            {x ≤ 0}

        """
        geodesic = self.geodesic(a, b, c, model=model, check=check)

        return self.__make_element_class__(HyperbolicHalfSpace)(self, geodesic)

    def segment(self, geodesic, start=None, end=None, check=True, assume_normalized=False):
        r"""
        Return the segment on the ``geodesic`` bounded by ``start`` and ``end``.

        INPUT:

        - ``geodesic`` -- a :meth:`geodesic` in this space.

        - ``start`` -- ``None`` or a :meth:`point` on the ``geodesic``, e.g.,
          obtained from the :meth:`HyperbolicGeodesic.intersection` of
          ``geodesic`` with another geodesic. If ``None``, the segment starts
          at the infinite :meth:`HyperbolicGeodesic.start` point of the geodesic.

        - ``end`` -- ``None`` or a :meth:`point` on the ``geodesic``, same as
          ``start``; must be later on ``geodesic`` than ``start``.

        - ``check`` -- boolean (default: ``True``), whether validation is
          performed on the arguments.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        When neither ``start`` nor ``end`` are given, a geodesic is returned::

            sage: H.segment(H.vertical(0), start=None, end=None)
            {-x = 0}

        When only one endpoint is provided, the segment is infinite on one end::

            sage: H.segment(H.vertical(0), start=I, end=None)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        When both endpoints are provided, a proper closed segment is returned::

            sage: H.segment(H.vertical(0), start=I, end=2*I)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}

        However, ideal endpoints on the geodesic are ignored::

            sage: H.segment(H.vertical(0), start=0, end=oo)
            {-x = 0}

        A segment can be reduced to a single point::

            sage: H.segment(H.vertical(0), start=I, end=I)
            I

        The endpoints must have coordinates over the base ring::

            sage: H.segment(H.half_circle(0, 2), H.half_circle(0, 2).start(), H.half_circle(0, 2).end())
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not a rational number

        """
        geodesic = self(geodesic)

        if not isinstance(geodesic, HyperbolicGeodesic):
            raise TypeError("geodesic must be a geodesic")

        if start is not None:
            start = self(start)
            if not isinstance(start, HyperbolicPoint):
                raise TypeError("start must be a point")

        if end is not None:
            end = self(end)
            if not isinstance(end, HyperbolicPoint):
                raise TypeError("end must be a point")

        segment = self.__make_element_class__(HyperbolicSegment)(self, geodesic, start, end)

        if check:
            segment._check(require_normalized=False)

        if not assume_normalized:
            segment = segment._normalize()

        if check:
            segment._check(require_normalized=True)

        return segment

    def polygon(self, half_spaces, check=True, assume_sorted=False, assume_minimal=False):
        r"""
        Return the convex polygon obtained by intersecting ``half_spaces``.

        See :meth:`intersection` for algorithmic details.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A finite convex polygon::

            sage: H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 2).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ])
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        Redundant half spaces are removed from the final representation::

            sage: H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 2).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....:   H.half_circle(0, 6).right_half_space(),
            ....: ])
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        The vertices of the polygon can be at ideal points; this polygon has
        vertices at -1 and 1::

            sage: H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ])
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        Any set of half spaces defines a polygon, even if the edges do not even
        meed at ideal points::

            sage: H.polygon([
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 2).right_half_space(),
            ....: ])
            {(x^2 + y^2) - 2 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        However, when the resulting set is point, the result is not represented
        as a polygon anymore::

            sage: H.polygon([
            ....:   H.vertical(-1).left_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....: ])
            ∞

        We can force the creation of this set as a polygon which might be
        beneficial in some algorithmic applications::

            sage: H.polygon([
            ....:   H.vertical(-1).left_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....: ], check=False, assume_minimal=True)
            {x + 1 ≤ 0} ∩ {x - 1 ≥ 0}

        Note that forcing this mode does not remove redundant half spaces from
        the representation; we usually assume that the representation is
        minimal, so such a polygon might not behave correctly::

            sage: H.polygon([
            ....:   H.vertical(-1).left_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....:   H.vertical(2).right_half_space(),
            ....: ], check=False, assume_minimal=True)
            {x + 1 ≤ 0} ∩ {x - 1 ≥ 0} ∩ {x - 2 ≥ 0}

        We could manually pass to a minimal representation by rewriting the
        point as half spaces again::

            sage: minimal = H.polygon([
            ....:   H.vertical(-1).left_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....:   H.vertical(2).right_half_space(),
            ....: ])
            sage: H.polygon(minimal._half_spaces(), check=False, assume_minimal=True)
            {x ≤ 0} ∩ {x - 1 ≥ 0}

        Note that this chose half spaces not in the original set; you might
        also want to have a look at :meth:`HyperbolicConvexPolygon._normalize`
        for some ideas how to manually reduce the half spaces that are used in
        a polygon.

        Note that the same applies if the intersection of half spaces is empty
        or just a single half space::

            sage: empty = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....:   H.half_circle(0, 2).left_half_space(),
            ....: ])
            sage: type(empty)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicEmptySet_with_category'>

        ::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ])
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicHalfSpace_with_category'>

        The intersection of the half spaces is computed in time quasi-linear in
        the number of half spaces. The limiting factor is sorting the half
        spaces by :meth:`HyperbolicHalfSpace._less_than`. If we know that the
        half spaces are already sorted like that, we can make the process run
        in linear time by setting ``assume_sorted``.

            sage: H.polygon(H.infinity()._half_spaces(), assume_sorted=True)
            ∞

        """
        half_spaces = [self.coerce(half_space) for half_space in half_spaces]

        if not assume_sorted:
            half_spaces = HyperbolicHalfSpace._merge_sorted(*[[half_space] for half_space in half_spaces])

        polygon = self.__make_element_class__(HyperbolicConvexPolygon)(self, half_spaces)

        if check:
            polygon._check(require_normalized=False)

        if check or not assume_minimal:
            polygon = polygon._normalize()

        if check:
            polygon._check()

        return polygon

    def intersection(self, *subsets):
        r"""
        Return the intersection of convex ``subsets``.

        ALGORITHM:

        We compute the intersection of the
        :meth:`HyperbolicConvexSet._half_spaces` that make up the ``subsets``.
        That intersection can be computed in the Klein model where we can
        essentially reduce this problem to the intersection of half spaces in
        the Euclidean plane.

        The Euclidean intersection problem can be solved in time linear in the
        number of half spaces assuming that the half spaces are already sorted
        in a certain way. In particular, this is the case if there is only a
        constant number of ``subsets``. Otherwise, the algorithm is
        quasi-linear in the number of half spaces due to the added complexity
        of sorting.

        See :meth:`HyperbolicConvexPolygon._normalize` for more algorithmic details.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.intersection(H.vertical(0).left_half_space())
            {x ≤ 0}

            sage: H.intersection(H.vertical(0).left_half_space(), H.vertical(0).right_half_space())
            {x = 0}

        """
        subsets = [self.coerce(subset) for subset in subsets]

        half_spaces = [subset._half_spaces() for subset in subsets]

        half_spaces = HyperbolicHalfSpace._merge_sorted(*half_spaces)

        return self.polygon(half_spaces, assume_sorted=True, assume_minimal=False, check=False)

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


# TODO: Change richcmp to compare according to the subset relation.
class HyperbolicConvexSet(Element):
    r"""
    Base class for convex subsets of :class:`HyperbolicPlane`.

    .. NOTE::

        Concrete subclasses should apply the following rules.

        There should only be a single type to describe a certain subset:
        normally, a certain subset, say a point, should only described by a
        single class, namely :class:`HyperbolicPoint`. Of course, one could
        describe a point as a polygon delimited by some edges that all
        intersect in that single point, such objects should be avoided. Namely,
        the methods that create a subset, say :meth:`HyperbolicPlane.polygon`
        take care of this by calling a sets
        :meth:`HyperbolicConvexSet._normalize` to rewrite a set in its most
        natural representation. To get the denormalized representation, we can
        always set `check=False` when creating the object. For this to work,
        the `__init__` should not take care of any such normalization and
        accept any input that can possibly be made sense of.

        Comparison with :meth:`_richcmp_` should compare by inclusion of sets:
        Care has to be taken when implementing hashing here since it must be
        consistent with that partial ordering. In particular, we cannot assume
        that all objects are normalized, i.e., an object must return the hash
        of its normalization.

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicConvexSet
        sage: H = HyperbolicPlane()

        sage: isinstance(H(0), HyperbolicConvexSet)
        True

    """

    def _half_spaces(self):
        r"""
        Return a minimal set of half spaces whose intersection is this convex set.

        The half spaces are ordered by :meth:`HyperbolicHalfSpace._less_than`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).left_half_space()._half_spaces()
            [{x ≤ 0}]

            sage: H.vertical(0)._half_spaces()
            [{x ≤ 0}, {x ≥ 0}]

            sage: H(0)._half_spaces()
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

        """
        # TODO: Check that all subclasses implement this.
        raise NotImplementedError("Convex sets must implement this method.")

    def _check(self, require_normalized=True):
        r"""
        Validate this convex subset.

        If ``require_normalized``, we also check that the object has the
        correct implementation class, e.g., that a point is a
        :class:`HyperbolicPoint` and not say a :class:`HyperbolicSegment` of
        length zero.
        """
        pass

    def _normalize(self):
        r"""
        Return this set possibly rewritten in a simpler form.

        This method is only relevant for sets created with ``check=False``.
        Such sets might have been created in a non-canonical way, e.g., when
        creating a :class:`HyperbolicSegment` whose start and end point are ideal,
        then this is actually a geodesic and it shuold be described as such.
        """
        return self

    def intersection(self, other):
        r"""
        Return the intersection with the ``other`` convex set.
        """
        return self.parent().intersection(self, other)

    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this set.
        """
        for half_space in self._half_spaces():
            if point not in half_space:
                return False

        return True

    def is_finite(self):
        r"""
        Return whether all points in this set are finite.
        """
        raise NotImplementedError(f"this {type(self)} does not support checking finiteness")

    def change_ring(self, ring):
        r"""
        Return this set as an element of the hyperbolic plane over ``ring``.
        """
        raise NotImplementedError(f"this {type(self)} does not support changing the base ring")

    def _test_change_ring(self, **options):
        r"""
        Verify that this set implements :meth:`change_ring`.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_change_ring()

        """
        tester = self._tester(**options)
        tester.assertEqual(self, self.change_ring(self.parent().base_ring()))

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this subset.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).plot()
            Graphics object consisting of 1 graphics primitive

        """
        raise NotImplementedError(f"this {type(self)} does not support plotting")

    def _test_plot(self, **options):
        r"""
        Verify that this set implements :meth:`plot`.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_plot()

        """
        from sage.all import Graphics

        tester = self._tester(**options)
        if self != self.parent().infinity():
            tester.assertIsInstance(self.plot(), Graphics)
            tester.assertIsInstance(self.plot(model="half_plane"), Graphics)
        tester.assertIsInstance(self.plot(model="klein"), Graphics)

    def _check_isometry_klein(self, isometry):
        from sage.matrix.special import diagonal_matrix

        # TODO: check that isometry is actually a matrix?
        if isometry.nrows() != 3 or isometry.ncols() != 3 or not self.parent().base_ring().has_coerce_map_from(isometry.base_ring()):
            raise ValueError('invalid isometry')
        D = isometry.transpose() * diagonal_matrix([1, 1, -1]) * isometry
        if D[0, 1] or D[0, 2] or D[1, 0] or D[1, 2] or D[2, 0] or D[2, 1]:
            raise ValueError('invalid isometry')
        if D[0, 0].is_zero() or D[1, 1].is_zero() or D[2, 2].is_zero():
            raise ValueError('invalid isometry')
        if D[0, 0] != D[1, 1] or D[0, 0] != - D[2, 2]:
            raise ValueError('invalid isometry')

    def _apply_isometry_klein(self, isometry):
        raise NotImplementedError

    def apply_isometry(self, isometry, model="half_plane"):
        r"""
        Return the image of this set under the isometry.

        INPUT:

        - ``isometry`` -- a 2x2 matrix in `PGL(2,\mathbb{R})` or a 3x3 matrix in `SO(1, 2)`

        - ``model`` -- either ``"half_plane"`` or ``"klein"``
        """
        if model == "half_plane":
            isometry = sl2_to_so12(isometry)
            model = "klein"

        if model == "klein":
            self._check_isometry_klein(isometry)
            return self._apply_isometry_klein(isometry)

        raise NotImplementedError("applying isometry not supported in this hyperbolic model")

    def _neg_(self):
        r"""
        Return the convex subset obtained by taking the (closed) complements of
        the half spaces whose intersection define this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: -H.vertical(0).left_half_space()
            {x ≥ 0}

        """
        return self.parent().intersection(*[-half_space for half_space in self._half_spaces()])

    # TODO: Test that _richcmp_ can compare all kinds of sets by inclusion.
    # TODO: Provide hashing.

    def an_element(self):
        # TODO: Test that everything implements an_element().
        raise NotImplementedError

    @classmethod
    def _enhance_plot(self, plot, model):
        if model == "klein":
            from sage.all import circle
            plot = circle([0, 0], 1, fill=False, color='#d1d1d1', zorder=-1) + plot

        return plot

    def is_empty(self):
        return self.dimension() < 0

    def __bool__(self):
        return not self.is_empty()

    def dimension(self):
        # TODO: Test that this is an ZZ integer.
        raise NotImplementedError(f"{type(self)} does not implement dimension() yet")

    def is_point(self):
        return self.dimension() == 0


class HyperbolicHalfSpace(HyperbolicConvexSet):
    r"""
    A closed half space of the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)

        sage: H.half_circle(0, 1).left_half_space()
        {(x^2 + y^2) - 1 ≥ 0}

    """

    def __init__(self, parent, geodesic):
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicGeodesic):
            raise TypeError("geodesic must be a geodesic")

        self._geodesic = geodesic

    def equation(self, model):
        r"""
        Return an inequality for this half space as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is in the half space if it satisfies `a(x^2 + y^2) + bx + c \ge 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy \ge 0`.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.
        """
        return self._geodesic.equation(model=model)

    def _repr_(self):
        r"""
        Return a printable representation of this half space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: S = H.half_circle(0, 1).right_half_space()

            sage: S
            {(x^2 + y^2) - 1 ≤ 0}

            sage: -S
            {(x^2 + y^2) - 1 ≥ 0}

        """
        # Convert to the Poincaré half plane model as a(x^2 + y^2) + bx + c ≥ 0.
        a, b, c = self.equation(model="half_plane")

        try:
            from sage.all import gcd
            d = gcd((a, b, c))
            assert d > 0
            a /= d
            b /= d
            c /= d
        except Exception:
            pass

        # Remove any trailing - signs in the output.
        cmp = "≥"
        if a < 0 or (a == 0 and b < 0):
            a *= -1
            b *= -1
            c *= -1
            cmp = "≤"

        from sage.all import PolynomialRing
        R = PolynomialRing(self.parent().base_ring(), names="x")
        if a != 0:
            return f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} {cmp} 0}}"
        else:
            return f"{{{repr(R([c, b]))} {cmp} 0}}"

    def _half_spaces(self):
        r"""
        Implements :meth:`HyperbolicConvexSet._half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: S = H.vertical(0).left_half_space()
            sage: [S] == S._half_spaces()
            True

        """
        return [self]

    def _neg_(self):
        r"""
        Return the closure of the complement of this half space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: S = H.half_circle(0, 1).left_half_space()
            sage: -S
            {(x^2 + y^2) - 1 ≤ 0}

        """
        return self._geodesic.right_half_space()

    @staticmethod
    def _merge_sorted(*half_spaces):
        r"""
        Return the merge of lists of ``half_spaces``.

        The lists are assumed to be sorted by
        :meth:`HyperbolicHalfSpace._less_than` and are merged into a single
        list with that sorting.

        Naturally, when there are a lot of short lists, such a merge sort takes
        quasi-linear time. However, when there are only a few lists, this runs
        in linear time.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicHalfSpace
            sage: H = HyperbolicPlane()

            sage: HyperbolicHalfSpace._merge_sorted()
            []

            sage: HyperbolicHalfSpace._merge_sorted(H.real(0)._half_spaces())
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpace._merge_sorted(H.real(0)._half_spaces(), H.real(0)._half_spaces())
            [{(x^2 + y^2) + x ≤ 0}, {(x^2 + y^2) + x ≤ 0}, {x ≥ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpace._merge_sorted(*[[half_space] for half_space in H.real(0)._half_spaces() * 2])
            [{(x^2 + y^2) + x ≤ 0}, {(x^2 + y^2) + x ≤ 0}, {x ≥ 0}, {x ≥ 0}]

        """

        # A standard merge-sort implementation.
        count = len(half_spaces)

        if count == 0:
            return []

        if count == 1:
            return half_spaces[0]

        # The non-trivial base case.
        if count == 2:
            A = half_spaces[0]
            B = half_spaces[1]
            merged = []

            while A and B:
                if HyperbolicHalfSpace._less_than(A[-1], B[-1]):
                    merged.append(B.pop())
                else:
                    merged.append(A.pop())

            merged.reverse()

            return A + B + merged

        # Divide & Conquer recursively.
        return HyperbolicHalfSpace._merge_sorted(*(
            HyperbolicHalfSpace._merge_sorted(*half_spaces[: count // 2]),
            HyperbolicHalfSpace._merge_sorted(*half_spaces[count // 2:])
            ))

    @staticmethod
    def _less_than(lhs, rhs):
        # TODO: This is essentially atan2.
        r"""
        Return whether the half space ``lhs`` is smaller than ``rhs`` in a cyclic
        ordering of normal vectors, i.e., in an ordering that half spaces
        whether their normal points to the left/right, the slope of the
        geodesic, and finally by containment.

        This ordering is such that :meth:`HyperbolicPlane.intersection` can be
        computed in linear time for two hyperbolic convex sets.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicHalfSpace
            sage: H = HyperbolicPlane(QQ)

        A half space is equal to itself::

            sage: HyperbolicHalfSpace._less_than(H.vertical(0).left_half_space(), H.vertical(0).left_half_space())
            False

        A half space whose normal in the Klein model points to the left is
        smaller than one whose normal points to the right::

            sage: HyperbolicHalfSpace._less_than(H.vertical(1).left_half_space(), H.half_circle(0, 1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.vertical(0).left_half_space(), -H.vertical(0).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(-H.half_circle(-1, 1).left_half_space(), -H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(-H.half_circle(-1, 1).left_half_space(), -H.vertical(1/2).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.vertical(1).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.vertical(1/2).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True

        Half spaces are ordered by the slope of their normal in the Klein model::

            sage: HyperbolicHalfSpace._less_than(H.vertical(-1).left_half_space(), H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(-H.half_circle(-1, 1).left_half_space(), H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.half_circle(-1, 1).left_half_space(), -H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.vertical(0).left_half_space(), H.vertical(1).left_half_space())
            True

        Parallel half spaces in the Klein model are ordered by inclusion::

            sage: HyperbolicHalfSpace._less_than(-H.half_circle(-1, 1).left_half_space(), H.vertical(1/2).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(-H.vertical(1/2).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.half_circle(0, 2).left_half_space(), H.half_circle(0, 1).left_half_space())
            True
            sage: HyperbolicHalfSpace._less_than(H.half_circle(0, 1).right_half_space(), H.half_circle(0, 2).right_half_space())
            True

        Verify that comparisons are projective::

            sage: HyperbolicHalfSpace._less_than(H.geodesic(5, -5, -1, model="half_plane").left_half_space(), H.geodesic(5/13, -5/13, -1/13, model="half_plane").left_half_space())
            False
            sage: HyperbolicHalfSpace._less_than(H.geodesic(5/13, -5/13, -1/13, model="half_plane").left_half_space(), H.geodesic(5, -5, -1, model="half_plane").left_half_space())
            False

        """
        a, b, c = lhs.equation(model="klein")
        aa, bb, cc = rhs.equation(model="klein")

        def normal_points_left(b, c):
            return b < 0 or (b == 0 and c < 0)

        if normal_points_left(b, c) != normal_points_left(bb, cc):
            # The normal vectors of the half spaces in the Klein model are in
            # different half planes, one is pointing left, one is pointing
            # right.
            return normal_points_left(b, c)

        # The normal vectors of the half spaces in the Klein model are in the
        # same half plane, so we order them by slope.
        if b * bb == 0:
            if b == bb:
                # The normals are vertical and in the same half plane, so
                # they must be equal. We will order the half spaces by
                # inclusion later.
                cmp = 0
            else:
                # Exactly one of the normals is vertical; we order half spaces
                # such that that one is bigger.
                return bb == 0
        else:
            # Order by the slope of the normal.
            cmp = (b * bb).sign() * (c * bb - cc * b).sign()

        if cmp == 0:
            # The half spaces are parallel in the Klein model. We order them by
            # inclusion, i.e., by the offset in direction of the normal.
            if c * cc:
                cmp = c.sign() * (a * cc - aa * c).sign()
            else:
                assert b * bb
                cmp = b.sign() * (a * bb - aa * b).sign()

        return cmp < 0

    def boundary(self):
        r"""
        Return a geodesic on the boundary of this half space, oriented such that the half space is on its left.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: S = H.vertical(0).left_half_space()
            sage: S.boundary()
            {-x = 0}

        """
        return self._geodesic

    def __contains__(self, point):
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        return a + b * x + c * y >= 0

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicHalfSpace):
                return False
            return self._geodesic._richcmp_(other._geodesic, op)

    def plot(self, model="half_plane", **kwds):
        return self.parent().polygon([self], check=False, assume_minimal=True).plot(model=model, **kwds)

    def change_ring(self, ring):
        return self._geodesic.change_ring(ring).left_half_space()

    def dimension(self):
        from sage.all import ZZ
        return ZZ(2)


class HyperbolicGeodesic(HyperbolicConvexSet):
    r"""
    An oriented geodesic in the hyperbolic plane.

    Internally, we represent geodesics as the chords satisfying the equation `a
    + bx + cy=0` in the unit disk of the Klein model.

    The geodesic is oriented such that the half space `a + bx + cy ≥ 0` is on
    its left.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)

        sage: H.vertical(0)
        {-x = 0}

        sage: H.half_circle(0, 1)
        {(x^2 + y^2) - 1 = 0}

        sage: H.geodesic(H(I), 0)
        {x = 0}

    """

    def __init__(self, parent, a, b, c):
        super().__init__(parent)

        if not isinstance(a, Element) or a.parent() is not parent.base_ring():
            raise TypeError("a must be an element of the base ring")
        if not isinstance(b, Element) or b.parent() is not parent.base_ring():
            raise TypeError("b must be an element of the base ring")
        if not isinstance(c, Element) or c.parent() is not parent.base_ring():
            raise TypeError("c must be an element of the base ring")

        self._a = a
        self._b = b
        self._c = c

    def _check(self, require_normalized=True):
        if self.is_ultra_ideal():
            raise ValueError(f"equation {self._a} + ({self._b})*x + ({self._c})*y = 0 does not define a chord in the Klein model")

    def is_ultra_ideal(self):
        return self._b*self._b + self._c*self._c <= self._a*self._a

    def _half_spaces(self):
        r"""
        Implements :meth:`HyperbolicConvexSet._half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0)._half_spaces()
            [{x ≤ 0}, {x ≥ 0}]

        """
        return HyperbolicHalfSpace._merge_sorted([self.left_half_space()], [self.right_half_space()])

    def start(self):
        r"""
        Return the ideal starting point of this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).start()
            0

        The coordinates of the end points of the half circle of radius
        `\sqrt{2}` around 0 can not be written down in the rationals::

            sage: p = H.half_circle(0, 2).start()
            sage: p
            -1.41421356237310

            sage: p.coordinates()
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not a rational number

        Passing to a bigger field, the coordinates can be represented::

            sage: K.<a> = QQ.extension(x^2 - 2, embedding=1.4)
            sage: H.half_circle(0, 2).change_ring(K).start()
            -a

        """
        return self.parent().point(x=self, y=None, model=None, check=False)

    def end(self):
        r"""
        Return the ideal end point of this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).end()
            ∞

        The coordinates of the end points of the half circle of radius
        `\sqrt{2}` around 0 can not be written down in the rationals::

            sage: p = H.half_circle(0, 2).end()
            sage: p
            1.41421356237310

            sage: p.coordinates()
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not a rational number

        Passing to a bigger field, the coordinates can be represented::

            sage: K.<a> = QQ.extension(x^2 - 2, embedding=1.4)
            sage: H.half_circle(0, 2).change_ring(K).end()
            a

        """
        return (-self).start()

    def _neg_(self):
        r"""
        Return the reversed geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: -H.vertical(0)
            {x = 0}

        """
        return self.parent().geodesic(-self._a, -self._b, -self._c, model="klein")

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
            from sage.all import gcd
            d = gcd((a, b, c))
            assert d > 0
            a /= d
            b /= d
            c /= d
        except Exception:
            pass

        from sage.all import PolynomialRing
        R = PolynomialRing(self.parent().base_ring(), names="x")
        if a != 0:
            return f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} = 0}}"
        else:
            return f"{{{repr(R([c, b]))} = 0}}"

    def plot(self, model="half_plane", **kwds):
        r"""
        Create a plot of this geodesic in the hyperbolic ``model``.

        Additional arguments are passed on to the underlying SageMath plotting methods.
        """
        return self.parent().segment(self, start=None, end=None, check=False, assume_normalized=True).plot(model=model, **kwds)

    def change_ring(self, ring):
        r"""
        Return this geodesic in the Hyperbolic plane over ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(1).change_ring(QQ)
            {-x + 1 = 0}

            sage: H.vertical(AA(2).sqrt()).change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real ... to Rational

        """
        return HyperbolicPlane(ring).geodesic(self._a, self._b, self._c, model="klein")

    def left_half_space(self):
        r"""
        Return the closed half space to the left of this (oriented) geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(0).left_half_space()
            {x ≤ 0}

        """
        return self.parent().half_space(self._a, self._b, self._c, model="klein")

    def right_half_space(self):
        r"""
        Return the closed half space to the right of this (oriented) geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(0).right_half_space()
            {x ≥ 0}

        """
        return (-self).left_half_space()

    def an_element(self):
        a, b, c = self.equation(model="klein")

        return self.parent().geodesic(0, -c, b, model="klein", check=False)._intersection(self)

    def _configuration(self, other):
        r"""
        Return a classification of the angle between this
        geodesic and ``other`` in the Klein model.

        """
        # TODO: Can we make this public somehow?
        intersection = self._intersection(other)

        if intersection is None:
            orientation = (self._b * other._b + self._c * other._c).sign()

            assert orientation != 0

            if self == other:
                assert orientation > 0
                return "equal"

            if self == -other:
                assert orientation < 0
                return "negative"

            if orientation > 0:
                return "parallel"

            return "anti-parallel"

        tangent = (self._c, -self._b)
        orientation = (-tangent[0] * other._b - tangent[1] * other._c).sign()

        assert orientation != 0

        # TODO: Is convex/concave the right term?
        if orientation > 0:
            return "convex"

        # TODO: Use the bounding-box trick to get mostly rid of concave cases.
        return "concave"

    def _intersection(self, other):
        r"""
        Return the intersection of this geodesic and ``other`` in the Klein
        model or in the Euclidean plane if the intersection point is ultra
        ideal, i.e., not in the unit disk.

        Returns ``None`` if the lines do not intersect in a point.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

        ::

            sage: A = -H.vertical(0)
            sage: B = H.vertical(-1)
            sage: C = H.vertical(0)
            sage: A._intersection(B)
            ∞
            sage: A._intersection(C)
            sage: B._intersection(A)
            ∞
            sage: B._intersection(C)
            ∞
            sage: C._intersection(A)
            sage: C._intersection(B)
            ∞

        """
        if not isinstance(other, HyperbolicGeodesic):
            raise TypeError("can only intersect with another geodesic")

        # TODO: Reference the trac ticket that says that solving for 2×2 matrices is very slow.
        det = self._b * other._c - self._c * other._b

        if det == 0:
            return None

        x = (-other._c * self._a + self._c * other._a) / det
        y = (other._b * self._a - self._b * other._a) / det

        return self.parent().point(x, y, model="klein", check=False)

    def __contains__(self, point):
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        return a + b * x + c * y == 0

    def parametrize(self, point, model, check=True):
        if isinstance(point, HyperbolicPoint):
            if check and point not in self:
                raise ValueError("point must be on geodesic to be parametrized")

        if model == "euclidean":
            base = self.an_element().coordinates(model="klein")
            tangent = (self._c, -self._b)

            if isinstance(point, HyperbolicPoint):
                coordinate = 0 if tangent[0] else 1
                return (point.coordinates(model="klein")[coordinate] - base[coordinate]) / tangent[coordinate]

            λ = self.parent().base_ring()(point)

            return self.parent().point(
                x=base[0] + λ * tangent[0],
                y=base[1] + λ * tangent[1],
                model="klein",
                check=check)

        raise NotImplementedError("cannot parametrize a geodesic over this model yet")

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicGeodesic):
                return False
            if self._b:
                return self._b.sign() == other._b.sign() and self._a * other._b == other._a * self._b and self._c * other._b == other._c * self._b
            else:
                return self._c.sign() == other._c.sign() and self._a * other._c == other._a * self._c and self._b * other._c == other._b * self._c

        super()._richcmp_(other, op)

    def _apply_isometry_klein(self, isometry):
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: p0 = H(0)
            sage: p1 = H(1)
            sage: p2 = H(oo)
            sage: for (a, b, c, d) in [(2, 1, 1, 1), (1, 1, 0, 1), (1, 0, 1, 1), (2, 0, 0 , 1)]:
            ....:     m = matrix(2, [a, b, c, d])
            ....:     q0 = p0.apply_isometry(m)
            ....:     q1 = p1.apply_isometry(m)
            ....:     q2 = p2.apply_isometry(m)
            ....:     assert H.geodesic(p0, p1).apply_isometry(m) == H.geodesic(q0, q1)
            ....:     assert H.geodesic(p1, p0).apply_isometry(m) == H.geodesic(q1, q0)
            ....:     assert H.geodesic(p1, p2).apply_isometry(m) == H.geodesic(q1, q2)
            ....:     assert H.geodesic(p2, p1).apply_isometry(m) == H.geodesic(q2, q1)
            ....:     assert H.geodesic(p2, p0).apply_isometry(m) == H.geodesic(q2, q0)
            ....:     assert H.geodesic(p0, p2).apply_isometry(m) == H.geodesic(q0, q2)
        """
        from sage.modules.free_module_element import vector

        b, c, a = vector(self.parent().base_ring(), [self._b, self._c, self._a]) * isometry.inverse()
        return self.parent().geodesic(a, b, c, model="klein")

    def dimension(self):
        from sage.all import ZZ
        return ZZ(1)


class HyperbolicPoint(HyperbolicConvexSet):
    r"""
    A (possibly infinite) point in the :class:`HyperbolicPlane`.

    Internally, we typically represent a point as the Euclidean coordinates in
    the unit disk of the Klein model.

    Additionally, we allow points that are the (ideal) endpoints of geodesics
    even if these only have coordinates over a quadratic extension.
    """

    def __init__(self, parent, x, y):
        super().__init__(parent)

        if y is None:
            if not isinstance(x, HyperbolicGeodesic):
                raise TypeError("x must be a geodesic")

            self._coordinates = x
        else:
            if x.parent() is not parent.base_ring():
                raise TypeError("x must be an element of the base ring")
            if y.parent() is not parent.base_ring():
                raise TypeError("y must be an element of the base ring")

            self._coordinates = (x, y)

    def _check(self, require_normalized=True):
        if self.is_ultra_ideal():
            raise ValueError(f"point {self} is not in the unit disk in the Klein model")

    def is_ultra_ideal(self):
        x, y = self.coordinates(model="klein")
        return x*x + y*y > 1

    def _half_spaces(self):
        r"""
        Implements :meth:`HyperbolicConvexSet._half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H(I)._half_spaces()
            [{(x^2 + y^2) + 2*x - 1 ≤ 0}, {x ≥ 0}, {(x^2 + y^2) - 1 ≥ 0}]

            sage: H(I + 1)._half_spaces()
            [{x - 1 ≤ 0}, {(x^2 + y^2) - 3*x + 1 ≤ 0}, {(x^2 + y^2) - 2 ≥ 0}]

            sage: H.infinity()._half_spaces()
            [{x ≤ 0}, {x - 1 ≥ 0}]

            sage: H(0)._half_spaces()
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

            sage: H(-1)._half_spaces()
            [{x + 1 ≤ 0}, {(x^2 + y^2) - 1 ≤ 0}]

            sage: H(1)._half_spaces()
            [{(x^2 + y^2) - x ≤ 0}, {(x^2 + y^2) - 1 ≥ 0}]

            sage: H(2)._half_spaces()
            [{2*(x^2 + y^2) - 3*x - 2 ≥ 0}, {3*(x^2 + y^2) - 7*x + 2 ≤ 0}]

            sage: H(-2)._half_spaces()
            [{(x^2 + y^2) - x - 6 ≥ 0}, {2*(x^2 + y^2) + 3*x - 2 ≤ 0}]

            sage: H(1/2)._half_spaces()
            [{6*(x^2 + y^2) - x - 1 ≤ 0}, {2*(x^2 + y^2) + 3*x - 2 ≥ 0}]

            sage: H(-1/2)._half_spaces()
            [{2*(x^2 + y^2) + 7*x + 3 ≤ 0}, {2*(x^2 + y^2) - 3*x - 2 ≤ 0}]

        For ideal endpoints of geodesics that do not have coordinates over the
        base ring, we cannot produce defining half spaces since these would
        require equations over a quadratic extension as well::

            sage: H.half_circle(0, 2).start()._half_spaces()
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not a rational number

        """
        x0, y0 = self.coordinates(model="klein")

        if self.is_finite():
            return HyperbolicHalfSpace._merge_sorted(
                # x ≥ x0
                [self.parent().half_space(-x0, 1, 0, model="klein")],
                # y ≥ y0
                [self.parent().half_space(-y0, 0, 1, model="klein")],
                # x + y ≤ x0 + y0
                [self.parent().half_space(x0 + y0, -1, -1, model="klein")])
        else:
            return HyperbolicHalfSpace._merge_sorted(
                # left of the line from (0, 0) to this point
                [self.parent().half_space(0, -y0, x0, model="klein")],
                # right of a line to this point with a starting point right of (0, 0)
                [self.parent().half_space(-x0*x0 - y0*y0, y0 + x0, y0 - x0, model="klein")])

    def is_finite(self):
        if isinstance(self._coordinates, HyperbolicGeodesic):
            return False

        x, y = self.coordinates(model="klein")
        return x*x + y*y < 1

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
        # TODO: Implement ring

        if model == "half_plane":
            x, y = self.coordinates(model="klein")

            if x == 0 and y == 1:
                raise ValueError(f"{self} has no coordinates in the upper half plane")

            denominator = 1 - y

            # TODO: Do something better here
            if (1 - x*x - y*y).abs() < 1e-6:
                return (x / denominator, self.parent().base_ring().zero())

            return (x / denominator, (1 - x*x - y*y).sqrt(extend=False)/denominator)

        if model == "klein":
            if isinstance(self._coordinates, HyperbolicGeodesic):
                a, b, c = self._coordinates.equation(model="half_plane")

                if a == 0:
                    if b > 0:
                        self._coordinates = self.parent().infinity()._coordinates
                    else:
                        self._coordinates = self.parent().real(-c/b)._coordinates
                else:
                    discriminant = b*b - 4*a*c
                    root = discriminant.sqrt(extend=False)

                    endpoints = ((-b - root) / (2*a), (-b + root) / (2*a))

                    if a > 0:
                        self._coordinates = self.parent().real(min(endpoints))._coordinates
                    else:
                        self._coordinates = self.parent().real(max(endpoints))._coordinates

            return self._coordinates

        raise NotImplementedError

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other``, see
        :meth:`HyperbolicConvexSet._richcmp_`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: H.infinity() == H.projective(1, 0)
            True

        TESTS:

        We can compare points even though their coordinates are only defined
        over a quadratic extension::

            sage: H.half_circle(0, 2).start() == H.half_circle(0, 2).start()
            True

            sage: H.half_circle(0, 2).start() == H.half_circle(0, 2).end()
            False

            sage: H.half_circle(0, 2).start() == H(0)
            False

            sage: H.half_circle(0, 1).end() == H(1)
            True

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicPoint):
                return False
            if isinstance(self._coordinates, HyperbolicGeodesic):
                if isinstance(other._coordinates, HyperbolicGeodesic):
                    if self._coordinates == other._coordinates:
                        return True
                    if self._coordinates == -other._coordinates:
                        return False
                    intersection = self._coordinates._intersection(other._coordinates)
                    if intersection is None or intersection.is_ultra_ideal():
                        return False
                    if intersection.is_finite():
                        return False
                    return self._coordinates.parametrize(intersection, model="euclidean") < 0 and other._coordinates.parametrize(intersection, model="euclidean") < 0
                if other.is_finite():
                    return False
                if other in self._coordinates and self._coordinates.parametrize(other, model="euclidean") < 0:
                    return True
                return False
            elif isinstance(other._coordinates, HyperbolicGeodesic):
                return other == self
            return self.coordinates(model="klein") == other.coordinates(model="klein")

        super()._richcmp_(other, op)

    def _repr_(self):
        if self == self.parent().infinity():
            return "∞"

        try:
            x, y = self.coordinates()
        except ValueError:
            # Implement better printing.
            from sage.all import RR
            return repr(self.change_ring(RR))

        # We represent x + y*I in R[[I]] so we do not have to reimplement printing ourselves.
        if x not in self.parent().base_ring() or y not in self.parent().base_ring():
            x, y = self.coordinates(model="klein")
            return f"({repr(x)}, {repr(y)})"

        # TODO: Map the ultra-ideal points to the negative half plane.

        # TODO: This does not work when the coordinates are not in the base_ring.
        from sage.all import PowerSeriesRing
        return repr(PowerSeriesRing(self.parent().base_ring(), names="I")([x, y]))

    def change_ring(self, ring):
        if isinstance(self._coordinates, HyperbolicGeodesic):
            return HyperbolicPlane(ring).point(self._coordinates, y=None, model=None, check=False)

        return HyperbolicPlane(ring).point(*self._coordinates, model="klein", check=False)

    def _apply_isometry_klein(self, isometry):
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: for (a, b, c, d) in [(2, 1, 1, 1), (1, 1, 0, 1), (1, 0, 1, 1), (2, 0, 0 , 1)]:
            ....:     m = matrix(2, [a, b, c, d])
            ....:     assert H(0).apply_isometry(m) == H(b / d if d else oo)
            ....:     assert H(1).apply_isometry(m) == H((a + b) / (c + d) if c+d else oo)
            ....:     assert H(oo).apply_isometry(m) == H(a / c if c else oo)
        """
        from sage.modules.free_module_element import vector

        x, y = self.coordinates(model="klein")

        x, y, z = isometry * vector(self.parent().base_ring(), [x, y, 1])
        return self.parent().point(x / z, y / z, model="klein")

    def plot(self, model="half_plane", **kwds):
        if model == "half_plane":
            from sage.all import point
            plot = point(self.coordinates(model="half_plane"), **kwds)
        elif model == "klein":
            from sage.all import point
            plot = point(self.coordinates(model="klein"), **kwds)
        else:
            raise NotImplementedError

        return self._enhance_plot(plot, model=model)

    def dimension(self):
        from sage.all import ZZ
        return ZZ.zero()


class HyperbolicConvexPolygon(HyperbolicConvexSet):
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`HyperbolicHalfSpace`s.
    """

    def __init__(self, parent, half_spaces):
        super().__init__(parent)

        if not isinstance(half_spaces, list):
            raise TypeError("half_spaces must be a list of half spaces")

        self._halfspaces = half_spaces

    def _check(self, require_normalized=True):
        # TODO
        pass

    # TODO: Add examples.
    def _normalize(self):
        # TODO: Check docstring.
        r"""
        Return a convex set describing the intersection of the half spaces underlying this polygon.

        The half spaces are assumed to be sorted respecting :meth:`HyperbolicHalfSpace._less_than`.

        ALGORITHM:

        We compute the intersection of the half spaces in the Klein model in several steps:

        * Drop trivially redundant half spaces, e.g., repeated ones.
        * Handle the case that the intersection is empty or a single point, see :meth:`_euclidean_boundary`.
        * Compute the intersection of the corresponding half spaces in the Euclidean plane, see :meth:`_normalize_drop_euclidean_redundant`.
        * Remove redundant half spaces that make no contribution for the unit disk of the Klein model, see :meth:`_normalize_drop_unit_disk_redundant`.
        * Determine of which nature (point, segment, line, polygon) the intersection of half spaces is and return the resulting set.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        An instance that caused problems at some point:

            sage: P = polygon(
            ....:   H.geodesic(7, -4, -3, model="half_plane").left_half_space(),
            ....:   H.geodesic(1, -1, 0, model="half_plane").left_half_space(),
            ....:   H.vertical(1/2).right_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....:   H.geodesic(1, 4, -5, model="half_plane").left_half_space(),
            ....:   H.geodesic(50, 57, -43, model="half_plane").left_half_space(),
            ....:   H.geodesic(3, 2, -3, model="half_plane").left_half_space()
            ....: )
            sage: P._normalize()
            {x - 1 ≥ 0}

        """
        self = self._normalize_drop_trivially_redundant()

        if not self._halfspaces:
            raise NotImplementedError("cannot model intersection of no half spaces yet")

        # Find a segment on the boundary of the intersection.
        boundary = self._euclidean_boundary()

        if not isinstance(boundary, HyperbolicHalfSpace):
            # When there was no such segment, i.e., the intersection is empty
            # or just a point, we are done.
            return boundary

        # Compute a minimal subset of the half spaces that defines the intersection in the Euclidean plane.
        self = self._normalize_drop_euclidean_redundant(boundary)

        # Remove half spaces that make no contribution when restricting to the unit disk of the Klein model.
        return self._normalize_drop_unit_disk_redundant()

    def _normalize_drop_trivially_redundant(self):
        # TODO: Check docstring.
        r"""
        Return a sublist of ``half_spaces`` without changing their intersection
        by removing some trivially redundant half spaces.

        The ``half_spaces`` are assumed to be sorted consistent with :meth:`HyperbolicHalfSpace._less_than`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        Repeated half spaces are removed::

            sage: polygon(H.vertical(0).left_half_space(), H.vertical(0).left_half_space())._normalize_drop_trivially_redundant()
            {x ≤ 0}

        Inclusions of half spaces are simplified::

            sage: polygon(H.vertical(0).left_half_space(), H.geodesic(1/2, 2).left_half_space())._normalize_drop_trivially_redundant()
            {x ≤ 0}

        But only if the inclusion is already present when extending the half
        space from the Klein model to the entire Euclidean plane::

            sage: polygon(H.vertical(0).left_half_space(), H.vertical(1).left_half_space())._normalize_drop_trivially_redundant()
            {x ≤ 0} ∩ {x - 1 ≤ 0}

        TESTS:

        The intersection of two half circles centered at 0::

            sage: polygon(*(H.half_circle(0, 1)._half_spaces() + H.half_circle(0, 2)._half_spaces()))._normalize_drop_trivially_redundant()
            {(x^2 + y^2) - 1 ≤ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        """
        reduced = []

        for half_space in self._halfspaces:
            if reduced:
                a, b, c = half_space.equation(model="klein")
                A, B, C = reduced[-1].equation(model="klein")

                if c * B == C * b and b.sign() == B.sign() and c.sign() == C.sign():
                    # The half spaces are parallel in the Euclidean plane. Since we
                    # assume spaces to be sorted by inclusion, we can drop this
                    # space.
                    continue

            reduced.append(half_space)

        return self.parent().polygon(reduced, check=False, assume_sorted=True, assume_minimal=True)

    def _normalize_drop_euclidean_redundant(self, boundary):
        # TODO: Check docstring.
        # TODO: Should we remove boundary from the interface?
        r"""
        Return a minimal sublist of ``half_spaces`` that describe their
        intersection as half spaces of the Euclidean plane.

        Consider the half spaces in the Klein model. Ignoring the unit disk,
        they also describe half spaces in the Euclidean plane.

        The half space ``boundary`` must be one of the ``half_spaces`` that
        defines a boundary edge of the intersection polygon in the Euclidean
        plane.

        ALGORITHM:

        We use an approach similar to gift-wrapping (but from the inside) to remove
        redundant half spaces from the input list. We start from the
        ``boundary`` which is one of the minimal half spaces and extend to the
        full intersection by walking the sorted half spaces.

        Since we visit each half space once, this reduction runs in linear time
        in the number of half spaces.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        An intersection which is a single point on the boundary of the unit
        disk::

            sage: polygon(*H.infinity()._half_spaces())._normalize_drop_euclidean_redundant(
            ....:     boundary=H.vertical(1).right_half_space())
            {x ≤ 0} ∩ {x - 1 ≥ 0}

        An intersection which is a segment outside of the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(0).left_half_space())
            {(x^2 + y^2) + 4*x + 3 ≤ 0} ∩ {x ≤ 0} ∩ {9*(x^2 + y^2) + 32*x + 25 ≥ 0} ∩ {x ≥ 0}

        An intersection which is a polygon outside of the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(1, -2, 0, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.half_space(17/8, 2, -1, model="klein"))
            {(x^2 + y^2) + 4*x + 3 ≤ 0} ∩ {(x^2 + y^2) - 4*x + 1 ≥ 0} ∩ {9*(x^2 + y^2) + 32*x + 25 ≥ 0} ∩ {x ≥ 0}

        An intersection which is an (unbounded) polygon touching the unit disk::

            sage: polygon(
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.vertical(1).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(1).right_half_space())
            {x + 1 ≤ 0} ∩ {x - 1 ≥ 0}

        An intersection which is a segment touching the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.geodesic(-1, -2).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(0).left_half_space())
            {x + 1 ≤ 0} ∩  {x ≤ 0} ∩ {(x^2 + y^2) + 3*x + 2 ≥ 0} ∩ {x ≥ 0}

        An intersection which is a polygon inside the unit disk::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.geodesic(0, 1).left_half_space())
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) + x ≥ 0}

        A polygon which has no vertices inside the unit disk but intersects the unit disk::

            sage: polygon(
            ....:     H.geodesic(2, 3).left_half_space(),
            ....:     H.geodesic(-3, -2).left_half_space(),
            ....:     H.geodesic(-1/2, -1/3).left_half_space(),
            ....:     H.geodesic(1/3, 1/2).left_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.geodesic(1/3, 1/2).left_half_space())
            {6*(x^2 + y^2) - 5*x + 1 ≥ 0} ∩ {(x^2 + y^2) - 5*x + 6 ≥ 0} ∩ {(x^2 + y^2) + 5*x + 6 ≥ 0} ∩ {6*(x^2 + y^2) + 5*x + 1 ≥ 0}

        A single half plane::

            sage: polygon(
            ....:     H.vertical(0).left_half_space()
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(0).left_half_space())
            {x ≤ 0}

        A pair of anti-parallel half planes::

            sage: polygon(
            ....:     H.geodesic(1/2, 2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.geodesic(-1/2, -2).right_half_space())
            {2*(x^2 + y^2) - 5*x + 2 ≥ 0} ∩ {2*(x^2 + y^2) + 5*x + 2 ≥ 0}

        A pair of anti-parallel half planes in the upper half plane::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(1).left_half_space())
            {x - 1 ≤ 0} ∩ {x + 1 ≥ 0}

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(-1).right_half_space())
            {x - 1 ≤ 0} ∩ {x + 1 ≥ 0}

        A segment in the unit disk with several superfluous half planes at infinity::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(1/2).left_half_space(),
            ....:     H.vertical(1/3).left_half_space(),
            ....:     H.vertical(1/4).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.vertical(-1/2).right_half_space(),
            ....:     H.vertical(-1/3).right_half_space(),
            ....:     H.vertical(-1/4).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(0).left_half_space())
            {x ≤ 0} ∩ {4*x + 1 ≥ 0} ∩ {x ≥ 0}

        A polygon in the unit disk with several superfluous half planes::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.vertical(2).left_half_space(),
            ....:     H.vertical(-2).right_half_space(),
            ....:     H.geodesic(0, 1/2).left_half_space(),
            ....:     H.geodesic(0, -1/2).right_half_space(),
            ....:     H.vertical(3).left_half_space(),
            ....:     H.vertical(-3).right_half_space(),
            ....:     H.geodesic(0, 1/3).left_half_space(),
            ....:     H.geodesic(0, -1/3).right_half_space(),
            ....: )._normalize_drop_euclidean_redundant(boundary=H.vertical(1).left_half_space())
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) + x ≥ 0}

        TESTS:

        An example that caused trouble at some point::

            sage: P = polygon(
            ....:   H.geodesic(7, -4, -3, model="half_plane").left_half_space(),
            ....:   H.geodesic(1, -1, 0, model="half_plane").left_half_space(),
            ....:   H.vertical(1/2).right_half_space(),
            ....:   H.vertical(1).right_half_space(),
            ....:   H.geodesic(1, 4, -5, model="half_plane").left_half_space(),
            ....:   H.geodesic(50, 57, -43, model="half_plane").left_half_space(),
            ....:   H.geodesic(3, 2, -3, model="half_plane").left_half_space()
            ....: )
            sage: P._normalize_drop_euclidean_redundant(boundary=P._euclidean_boundary())
            {(x^2 + y^2) - x ≥ 0} ∩ {2*x - 1 ≥ 0} ∩ {x - 1 ≥ 0}

        """
        # TODO: Make all other assumptions clear in the interface.

        half_spaces = self._halfspaces

        half_spaces = half_spaces[half_spaces.index(boundary):] + half_spaces[:half_spaces.index(boundary)]
        half_spaces.reverse()

        required_half_spaces = [half_spaces.pop()]

        while half_spaces:
            A = required_half_spaces[-1]
            B = half_spaces.pop()
            C = half_spaces[-1] if half_spaces else required_half_spaces[0]

            # Determine whether B is redundant, i.e., whether the intersection
            # A, B, C and A, C are the same.
            # Since we know that A is required and the space non-empty, the
            # question here is whether C blocks the line of sight from A to B.

            # We distinguish cases, depending of the nature of the intersection of A and B.
            AB = A.boundary()._configuration(B.boundary())
            BC = B.boundary()._configuration(C.boundary())
            AC = A.boundary()._configuration(C.boundary())

            required = False

            if AB == "convex":
                if BC == "concave":
                    assert AC in ["equal", "concave"]
                    required = True

                elif BC == "convex":
                    BC = B.boundary()._intersection(C.boundary())
                    required = AC == "negative" or (BC in A and BC not in A.boundary())

                elif BC == "negative":
                    required = True

                elif BC == "anti-parallel":
                    required = True

                else:
                    raise NotImplementedError(f"B and C are in unsupported configuration: {BC}")

            elif AB == "negative":
                required = True

            elif AB == "anti-parallel":
                required = True

            elif AB == "concave":
                required = True

            else:
                raise NotImplementedError(f"A and B are in unsupported configuration: {AB}")

            if required:
                required_half_spaces.append(B)
            elif len(required_half_spaces) > 1:
                half_spaces.append(required_half_spaces.pop())

        min = 0
        for i, half_space in enumerate(required_half_spaces):
            if HyperbolicHalfSpace._less_than(half_space, required_half_spaces[min]):
                min = i

        return self.parent().polygon(required_half_spaces[min:] + required_half_spaces[:min], check=False, assume_sorted=True, assume_minimal=True)

    def _normalize_drop_unit_disk_redundant(self):
        # TODO: Check docstring.
        r"""
        Return the intersection of the Euclidean ``half_spaces`` with the unit
        disk.

        The ``half_spaces`` must be minimal to describe their intersection in
        the Euclidean plane. If that intersection does not intersect the unit
        disk, then return the :meth:`empty_set`.

        Otherwise, return a minimal sublist of ``half_spaces`` that describes
        the intersection inside the unit disk.

        ALGORITHM:

        When passing to the Klein model, i.e., intersecting the polygon with the
        unit disk, some of the edges of the (possibly unbounded) polygon
        described by the ``half_spaces`` are unnecessary because they are not
        intersecting the unit disk.

        If none of the edges intersect the unit disk, then the polygon has
        empty intersection with the unit disk.

        Otherwise, we can drop the half spaces describing the edges that do not
        intersect the unit disk.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        An intersection which is a single point on the boundary of the unit
        disk::

            sage: polygon(*H.infinity()._half_spaces())._normalize_drop_unit_disk_redundant()
            ∞

        An intersection which is a segment outside of the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_drop_unit_disk_redundant()
            {}

        An intersection which is a polygon outside of the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(1, -2, 0, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_drop_unit_disk_redundant()
            {}

        An intersection which is an (unbounded) polygon touching the unit disk::

            sage: polygon(
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.vertical(1).right_half_space())._normalize_drop_unit_disk_redundant()
            ∞

        An intersection which is a segment touching the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.geodesic(-1, -2).right_half_space())._normalize_drop_unit_disk_redundant()
            ∞

        An intersection which is a polygon inside the unit disk::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space())._normalize_drop_unit_disk_redundant()
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) + x ≥ 0}

        A polygon which has no vertices inside the unit disk but intersects the unit disk::

            sage: polygon(
            ....:     H.geodesic(2, 3).left_half_space(),
            ....:     H.geodesic(-3, -2).left_half_space(),
            ....:     H.geodesic(-1/2, -1/3).left_half_space(),
            ....:     H.geodesic(1/3, 1/2).left_half_space())._normalize_drop_unit_disk_redundant()
            {6*(x^2 + y^2) - 5*x + 1 ≥ 0} ∩ {(x^2 + y^2) - 5*x + 6 ≥ 0} ∩ {(x^2 + y^2) + 5*x + 6 ≥ 0} ∩ {6*(x^2 + y^2) + 5*x + 1 ≥ 0}

        A single half plane::

            sage: polygon(H.vertical(0).left_half_space())._normalize_drop_unit_disk_redundant()
            {x ≤ 0}

        A pair of anti-parallel half planes::

            sage: polygon(
            ....:     H.geodesic(1/2, 2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).right_half_space())._normalize_drop_unit_disk_redundant()
            {2*(x^2 + y^2) - 5*x + 2 ≥ 0} ∩ {2*(x^2 + y^2) + 5*x + 2 ≥ 0}

        A segment in the unit disk with a superfluous half plane at infinity::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space())._normalize_drop_unit_disk_redundant()
            {x = 0}

        A polygon in the unit disk with several superfluous half planes::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.vertical(2).left_half_space(),
            ....:     H.geodesic(0, 1/2).left_half_space())._normalize_drop_unit_disk_redundant()
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) + x ≥ 0}

        A segment touching the inside of the unit disk::

            sage: polygon(
            ....:   H.vertical(1).left_half_space(),
            ....:   H.half_circle(0, 2).left_half_space(),
            ....:   H.vertical(0).left_half_space(),
            ....:   H.vertical(0).right_half_space(),
            ....: )._normalize_drop_unit_disk_redundant()
            {x = 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        An unbounded polygon touching the unit disk from the inside::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....: )._normalize_drop_unit_disk_redundant()
            {x - 1 ≤ 0} ∩ {x + 1 ≥ 0}

        A segment inside the unit disk::

            sage: polygon(
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(0).left_half_space(),
            ....:     H.geodesic(-2, 2).right_half_space(),
            ....:     H.geodesic(-1/2, 1/2).left_half_space(),
            ....: )._normalize_drop_unit_disk_redundant()
            {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {4*(x^2 + y^2) - 1 ≥ 0}

        """
        # TODO: Make all assumptions clear in the interface.

        required_half_spaces = []

        maybe_empty = True
        maybe_point = True
        maybe_segment = True

        for i in range(len(self._halfspaces)):
            A = self._halfspaces[i - 1]
            B = self._halfspaces[i]
            C = self._halfspaces[(i + 1) % len(self._halfspaces)]

            AB = None if A.boundary()._configuration(B.boundary()) == "concave" else A.boundary()._intersection(B.boundary())
            BC = None if B.boundary()._configuration(C.boundary()) == "concave" else B.boundary()._intersection(C.boundary())

            segment = self.parent().segment(B.boundary(), AB, BC, check=False)

            if isinstance(segment, HyperbolicEmptySet):
                pass
            elif isinstance(segment, HyperbolicPoint):
                maybe_empty = False

                if not maybe_point:
                    continue

                if maybe_point is True:
                    maybe_point = segment
                elif maybe_point != segment:
                    assert not maybe_point.is_finite()
                    assert not segment.is_finite()

                    maybe_point = False
            else:
                maybe_empty = False
                maybe_point = False

                if maybe_segment is True:
                    maybe_segment = segment
                elif maybe_segment == -segment:
                    # For the intersection to be only segment, we must see the
                    # segment twice, once from both sides.
                    return segment
                else:
                    maybe_segment = False

                required_half_spaces.append(B)

        if maybe_empty:
            return self.parent().empty_set()

        if maybe_point:
            return maybe_point

        if len(required_half_spaces) == 0:
            raise NotImplementedError("there is no convex set to represent the full space yet")

        if len(required_half_spaces) == 1:
            return required_half_spaces[0]

        return self.parent().polygon(required_half_spaces, check=False, assume_sorted=True, assume_minimal=True)

    def _euclidean_boundary(self):
        r"""
        TODO: This documentation is not entirely accurate anymore.

        Return a half space whose (Euclidean) boundary intersects the boundary
        of the intersection of ``half_spaces`` in more than a point.

        Consider the half spaces in the Klein model. Ignoring the unit disk,
        they also describe half spaces in the Euclidean plane.

        If their intersection contains a segment it must be on the boundary of
        one of the ``half_spaces`` which is returned by this method.

        If this is not the case, and the intersection is empty in the
        hyperbolic plane, return the :meth:`empty_set`. Otherwise, if the
        intersection is a point in the hyperbolic plane, return that point.

        The ``half_spaces`` must be sorted with respect to
        :meth:`HyperbolicHalfSpace._less_than`.

        ALGORITHM:

        We initially ignore the hyperbolic structure and just consider the half
        spaces of the Klein model as Euclidean half spaces.

        We use a relatively standard randomized optimization approach to find a
        point in the intersection: we randomly shuffle the half spaces and then
        optimize a segment on some boundary of the half spaces. The
        randomization makes this a linear time algorithm, see e.g.,
        https://www2.cs.arizona.edu/classes/cs437/fall07/Lecture4.prn.pdf

        If the only segment we can construct is a point, then the intersection
        is a single point in the Euclidean plane. The interesction in the
        hyperbolic plane might be a single point or empty.

        If not even a point exists, the intersection is empty in the Euclidean
        plane and therefore empty in the hyperbolic plane.

        Note that the segment returned might not be within the unit disk.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        Make the following randomized tests reproducible::

            sage: set_random_seed(0)

        An intersection that is already empty in the Euclidean plane::

            sage: polygon(
            ....:     H.geodesic(2, 1/2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).left_half_space()
            ....: )._euclidean_boundary()
            {}

        An intersection which in the Euclidean plane is a single point but
        outside the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(0, -1, 0, model="klein"),
            ....:     H.half_space(2, 2, -1, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....: )._euclidean_boundary()
            {}

        An intersection which is a single point inside the unit disk::

            sage: polygon(*H(I)._half_spaces())._euclidean_boundary()
            I

        An intersection which is a single point on the boundary of the unit
        disk::

            sage: polygon(*H.infinity()._half_spaces())._euclidean_boundary()
            {x - 1 ≥ 0}

        An intersection which is a segment outside of the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._euclidean_boundary()
            {x ≤ 0}

        An intersection which is a polygon outside of the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(1, -2, 0, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._euclidean_boundary()
            {9*(x^2 + y^2) + 32*x + 25 ≥ 0}

        An intersection which is an (unbounded) polygon touching the unit disk::

            sage: polygon(
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.vertical(1).right_half_space(),
            ....: )._euclidean_boundary()
            {x - 1 ≥ 0}

        An intersection which is a segment touching the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.geodesic(-1, -2).right_half_space(),
            ....: )._euclidean_boundary()
            {x ≥ 0}

        An intersection which is a polygon inside the unit disk::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....: )._euclidean_boundary()
            {(x^2 + y^2) - x ≥ 0}

        A polygon which has no vertices inside the unit disk but intersects the unit disk::

            sage: polygon(
            ....:     H.geodesic(2, 3).left_half_space(),
            ....:     H.geodesic(-3, -2).left_half_space(),
            ....:     H.geodesic(-1/2, -1/3).left_half_space(),
            ....:     H.geodesic(1/3, 1/2).left_half_space(),
            ....: )._euclidean_boundary()
            {6*(x^2 + y^2) - 5*x + 1 ≥ 0}

        A single half plane::

            sage: polygon(
            ....:     H.vertical(0).left_half_space()
            ....: )._euclidean_boundary()
            {x ≤ 0}

        A pair of anti-parallel half planes::

            sage: polygon(
            ....:     H.geodesic(1/2, 2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).right_half_space(),
            ....: )._euclidean_boundary()
            {2*(x^2 + y^2) - 5*x + 2 ≥ 0}

        TESTS:

        A case that caused problems at some point::

            sage: set_random_seed(1)

            sage: polygon(
            ....:    H.geodesic(300, 3389, -1166, model="half_plane").right_half_space(),
            ....:    H.geodesic(5, -24, -5, model="half_plane").left_half_space(),
            ....:    H.geodesic(182, -1135, 522, model="half_plane").left_half_space(),
            ....: )._euclidean_boundary()
            {5*(x^2 + y^2) - 24*x - 5 ≥ 0}

        """
        if len(self._halfspaces) == 0:
            raise ValueError("list of half spaces must not be empty")

        if len(self._halfspaces) == 1:
            return self._halfspaces[0]

        # Randomly shuffle the half spaces so the loop below runs in expected linear time.
        from sage.all import shuffle
        random_half_spaces = self._halfspaces[:]
        shuffle(random_half_spaces)

        # Move from the random starting point to a point that is contained in all half spaces.
        point = random_half_spaces[0].boundary().an_element()

        for half_space in random_half_spaces:
            if point in half_space:
                continue
            else:
                # The point is not in this half space. Find a point on the
                # boundary of half_space that is contained in all the half
                # spaces we have seen so far.
                boundary = half_space.boundary()

                # We parametrize the boundary points of half space, i.e., the
                # points that satisfy a + bx + cy = 0 by picking a base point B
                # and then writing points as (x, y) = B + λ(c, -b).

                # Each half space constrains the possible values of λ, starting
                # from (-∞,∞) to a smaller closed interval.
                from sage.all import RealSet, oo
                # Note that RealSet.real_line() would require SageMath 9.4
                interval = RealSet(-oo, oo)

                for constraining in random_half_spaces:
                    if constraining is half_space:
                        break

                    intersection = boundary._intersection(constraining.boundary())

                    if intersection is None:
                        # constraining is anti-parallel to half_space
                        if boundary.parametrize(0, model="euclidean", check=False) not in constraining:
                            return self.parent().empty_set()

                        # The intersection is non-empty, so this adds no further constraints.
                        continue

                    λ = boundary.parametrize(intersection, model="euclidean", check=False)

                    # Determine whether this half space constrains to (-∞, λ] or [λ, ∞).
                    if boundary.parametrize(λ + 1, model="euclidean", check=False) in constraining:
                        constraint = RealSet.unbounded_above_closed(λ)
                    else:
                        constraint = RealSet.unbounded_below_closed(λ)

                    interval = interval.intersection(constraint)

                    if interval.is_empty():
                        # The constraints leave no possibility for λ.
                        return self.parent().empty_set()

                # Construct a point from any of the λ in interval.
                λ = interval.an_element()

                point = boundary.parametrize(λ, model="euclidean", check=False)

        return self._extend_to_euclidean_boundary(point)

    def _extend_to_euclidean_boundary(self, point):
        half_spaces = [half_space for half_space in self._halfspaces if point in half_space.boundary()]

        if len(half_spaces) == 0:
            raise ValueError("point must be on the boundary of a defining half space")

        if len(half_spaces) < 3:
            return half_spaces[0]

        for (i, half_space) in enumerate(half_spaces):
            following = half_spaces[(i + 1) % len(half_spaces)]
            configuration = half_space.boundary()._configuration(following.boundary())

            if configuration == "convex":
                continue
            if configuration == "negative":
                return half_space

            if configuration == "concave":
                return half_space

            raise NotImplementedError(f"cannot extend point to segment when half spaces are in configuration {configuration}")

        if point.is_ultra_ideal():
            # There is no actual intersection in the hyperbolic plane.
            return self.parent().empty_set()

        return point

    def equations(self):
        r"""
        Return the equations describing the boundary of this polygon.

        The output is minimal and sorted by slope in the Klein model.
        """
        raise NotImplementedError

    def edges(self, as_segments=False):
        r"""
        Return the :class:`HyperbolicSegment`s and :class:`HyperbolicGeodesic`s defining this polygon.
        """
        edges = []

        boundaries = [half_space.boundary() for half_space in self._halfspaces]

        if len(boundaries) <= 1:
            return boundaries

        for i, B in enumerate(boundaries):
            A = boundaries[i - 1]
            C = boundaries[(i + 1) % len(boundaries)]

            AB = A._configuration(B)
            BC = B._configuration(C)

            start = None
            end = None

            if AB == "convex":
                start = A._intersection(B)
                if not start.is_finite():
                    start = None
            elif AB == "concave":
                pass
            else:
                raise NotImplementedError(f"cannot determine edges when boundaries are in configuration {AB}")

            if BC == "convex":
                end = B._intersection(C)
                if not end.is_finite():
                    end = None
            elif BC == "concave":
                pass
            else:
                raise NotImplementedError(f"cannot determine edges when boundaries are in configuration {BC}")

            edges.append(self.parent().segment(B, start=start, end=end, assume_normalized=as_segments, check=False))

        return edges

    def vertices(self):
        r"""
        Return the vertices of this polygon, i.e., the end points of the
        :meth:`edges`, in counterclockwise order.
        """
        raise NotImplementedError

    def _half_spaces(self):
        return self._halfspaces

    def _repr_(self):
        return " ∩ ".join([repr(half_space) for half_space in self._halfspaces])

    def plot(self, model="half_plane", **kwds):
        kwds.setdefault("color", "#efffff")
        kwds.setdefault("edgecolor", "#d1d1d1")

        if len(self._halfspaces) == 0:
            raise NotImplementedError("cannot plot full space")

        edges = self.edges(as_segments=True)

        pos = edges[0].start()

        commands = [BezierPath.Command("MOVETO", [pos])]

        for edge in edges:
            if edge.start() != pos:
                commands.append(BezierPath.Command("MOVETO", [edge.start()]))

            commands.append(BezierPath.Command("LINETO", [edge.end()]))
            pos = edge.end()

        if pos != edges[0].start():
            commands.append(BezierPath.Command("MOVETO", [edges[0].start()]))

        return self._enhance_plot(hyperbolic_path(commands, model=model, **kwds), model=model)

    def change_ring(self, ring):
        return HyperbolicPlane(ring).polygon([half_space.change_ring(ring) for half_space in self._halfspaces], check=False, assume_sorted=True, assume_minimal=True)

    def _richcmp_(self, other, op):
        # TODO: Pass to normalization
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicConvexPolygon):
                return False
            return self._halfspaces == other._halfspaces


class HyperbolicSegment(HyperbolicConvexSet):
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.
    """

    def __init__(self, parent, geodesic, start=None, end=None):
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicGeodesic):
            raise TypeError("geodesic must be a hyperbolic geodesic")

        if start is not None and not isinstance(start, HyperbolicPoint):
            raise TypeError("start must be a hyperbolic point")

        if end is not None and not isinstance(end, HyperbolicPoint):
            raise TypeError("start must be a hyperbolic point")

        self._geodesic = geodesic
        self._start = start
        self._end = end

    def _check(self, require_normalized=True):
        start = self._start
        end = self._end

        if start is not None:
            if start not in self._geodesic:
                raise ValueError("start point must be on the geodesic")

        if end is not None:
            if end not in self._geodesic:
                raise ValueError("end point must be on the geodesic")

        # TODO: Check end >= start and end > start if required_normalized.

        # TODO: Check start & end finite if require_normalized.

    def _normalize(self):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

        TESTS:

        We define a helper method for easier testing::

            sage: segment = lambda *args, **kwds: H.segment(*args, **kwds, check=False, assume_normalized=True)

        ::

            sage: segment(H.vertical(-1), start=H.infinity(), end=H.infinity())._normalize()
            ∞

        ::

            sage: segment(H.vertical(0), start=H.infinity(), end=None)._normalize()
            ∞

            sage: segment(H.vertical(0), start=None, end=H.infinity())._normalize()
            {-x = 0}

            sage: segment(-H.vertical(0), start=H.infinity(), end=None)._normalize()
            {x = 0}

            sage: segment(-H.vertical(0), start=None, end=H.infinity())._normalize()
            ∞

        ::

            sage: segment(H.vertical(0), start=I, end=H.infinity())._normalize()
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

            sage: segment(-H.vertical(0), start=H.infinity(), end=I)._normalize()
            {x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        """
        if self._geodesic.is_ultra_ideal():
            return self.parent().empty_set()

        start = self._start
        end = self._end

        if start is not None:
            λ_start = self._geodesic.parametrize(start, model="euclidean")

        if end is not None:
            λ_end = self._geodesic.parametrize(end, model="euclidean")

        if start is not None:
            if not start.is_finite():
                if λ_start > 0:
                    return self.parent().empty_set() if start.is_ultra_ideal() else start
                start = None

        if end is not None:
            if not end.is_finite():
                if λ_end < 0:
                    return self.parent().empty_set() if end.is_ultra_ideal() else end
                end = None

        if start is None and end is None:
            return self._geodesic

        assert (start is None or not start.is_ultra_ideal()) and (end is None or not end.is_ultra_ideal())

        if start == end:
            return start

        return self.parent().segment(self._geodesic, start=start, end=end, check=False, assume_normalized=True)

    @classmethod
    def _normalize_start(cls, geodesic, start, end):
        if start is None:
            return None, end

        start = geodesic.parent()(start)

        if not isinstance(start, HyperbolicPoint):
            raise TypeError("endpoint of segment must be a hyperbolic point")

        if start not in geodesic:
            raise ValueError("endpoint of segment must be on the geodesic")

        if not start.is_finite():
            if geodesic.parametrize(start, model="euclidean") < 0:
                return None, end
            if end is None:
                return start, start
            if end != start:
                raise ValueError("end point of segment cannot be before starting point")

        return start, end

    def _is_valid(self):
        if not self._geodesic._is_valid():
            return False

        if self._start is not None:
            if not self._start._is_valid():
                return False

            if self._start not in self._geodesic:
                return False

        if self._end is not None:
            if not self._end._is_valid():
                return False

            if self._end not in self._geodesic:
                return False

        if self._start is not None and self._end is not None:
            if self._start == self._end:
                return False

        # TODO: Check that the endpoints are ordered correctly.

        return True

    def _half_spaces(self):
        half_spaces = self._geodesic._half_spaces()

        if self._start is not None:
            half_spaces.append(self._start_half_space())

        if self._end is not None:
            half_spaces.append(self._end_half_space())

        return HyperbolicHalfSpace._merge_sorted(*[[half_space] for half_space in half_spaces])

    def _start_half_space(self):
        x, y = self._start.coordinates(model="klein")
        b, c = (self._geodesic._c, -self._geodesic._b)
        return self.parent().half_space(-b * x - c * y, b, c, model="klein", check=False)

    def _end_half_space(self):
        x, y = self._end.coordinates(model="klein")
        b, c = (-self._geodesic._c, self._geodesic._b)
        return self.parent().half_space(-b * x - c * y, b, c, model="klein", check=False)

    def _repr_(self):
        bounds = [repr(self._geodesic)]

        if self._start is not None:
            bounds.append(repr(self._start_half_space()))

        if self._end is not None:
            bounds.append(repr(self._end_half_space()))

        return " ∩ ".join(bounds)

    def _neg_(self):
        return self.parent().segment(-self._geodesic, self._end, self._start, check=False, assume_normalized=True)

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicSegment):
                return False
            return self._geodesic == other._geodesic and self._start == other._start and self._end == other._end

    def change_ring(self, ring):
        start = self._start.change_ring(ring) if self._start is not None else None
        end = self._end.change_ring(ring) if self._end is not None else None

        return HyperbolicPlane(ring).segment(self._geodesic.change_ring(ring), start=start, end=end, check=False, assume_normalized=True)

    def start(self, finite=False):
        if self._start is not None:
            return self._start

        return self._geodesic.start()

    def end(self, finite=False):
        if self._end is not None:
            return self._end

        return self._geodesic.end()

    def plot(self, model="half_plane", **kwds):
        from sage.all import RR

        self = self.change_ring(RR)
        plot = hyperbolic_path([BezierPath.Command("MOVETO", [self.start()]), BezierPath.Command("LINETO", [self.end()])], model=model, **kwds)

        return self._enhance_plot(plot, model=model)

    def configuration(self, other):
        if self._geodesic == other._geodesic:
            if self == other:
                return "equal"
            raise NotImplementedError("cannot determine configuration of segments on the same geodesic")

        if self._geodesic == -other._geodesic:
            if self == -other:
                return "negative"
            raise NotImplementedError("cannot determine configuration of segments on the same geodesic")

        intersection = self.intersection(other)

        if intersection is None:
            raise NotImplementedError("cannot determine configuration of segments that do not intersect")

        if intersection.is_finite():
            if intersection == self.end(finite=True) and intersection == other.start(finite=True):
                return "join"

            raise NotImplementedError("cannot determine configuration of segments that intersect in a finite point")

        if intersection == self.end():
            if intersection == other.start():
                return "join"
            if intersection == other.end():
                return "join-at-ends"

        if intersection == self.start():
            if intersection == other.start():
                return "join-at-starts"
            if intersection == other.end():
                return "join-reversed"

        raise NotImplementedError("cannot determine configuration of segments that intersect in an infinite point")


class HyperbolicEmptySet(HyperbolicConvexSet):
    r"""
    The empty subset of the hyperbolic plane.
    """

    def __init__(self, parent):
        super().__init__(parent)

    def _richcmp_(self, other, op):
        r"""
        Return how this set compares to ``other``.

        See :meth:`HyperbolicConvexSet._richcmp_`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: H.empty_set() == H.empty_set()
            True

        """
        from sage.structure.richcmp import rich_to_bool

        if isinstance(other, HyperbolicEmptySet):
            return rich_to_bool(op, 0)

        return rich_to_bool(op, -1)

    def _repr_(self):
        return "{}"

    def _apply_isometry_klein(self, isometry):
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: S = HyperbolicPlane(QQ).empty_set()
            sage: S.apply_isometry(matrix(2, [2, 1, 1, 1])) is S
            True
        """
        return self

    def plot(self, model="half_plane", **kwds):
        from sage.all import Graphics
        return self._enhance_plot(Graphics(), model=model)

    def change_ring(self, ring):
        return HyperbolicPlane(ring).empty_set()

    def dimension(self):
        from sage.all import ZZ
        return ZZ(-1)


def sl2_to_so12(m):
    r"""
    Return the lift of the 2x2 matrix ``m`` inside ``SO(1,2)``.
    """
    from sage.matrix.constructor import matrix

    if m.nrows() != 2 or m.ncols() != 2:
        raise ValueError('invalid matrix')
    a, b, c, d = m.list()
    return matrix(3, [a*d + b*c, a*c - b*d, a*c + b*d,
                      a*b - c*d, (a**2 - b**2 - c**2 + d**2) / 2, (a**2 + b**2 - c**2 - d**2) / 2,
                      a*b + c*d, (a**2 - b**2 + c**2 - d**2) / 2, (a**2 + b**2 + c**2 + d**2) / 2])


class BezierPath(GraphicPrimitive):
    # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
    def __init__(self, commands, options=None):
        options = options or {}

        valid_options = self._allowed_options()
        for option in options:
            if option not in valid_options:
                raise RuntimeError(f"option {option} not valid")

        self._commands = commands
        super().__init__(options)

    def _allowed_options(self):
        r"""
        Return the options that are supported by a path.

        We support all the options that are understood by a SageMath polygon.

        """
        from sage.plot.polygon import Polygon
        return Polygon([], [], {})._allowed_options()

    def _render_on_subplot(self, subplot):
        r"""
        Render this path on the subplot.

        Matplotlib was not really made to draw things that extend to infinity.
        The trick here is to register a callback that redraws whenever the
        viewbox of the plot changes, e.g., as more objects are added to the
        plot.
        """
        # Rewrite options to only contain matplotlib compatible entries
        matplotlib_options = {
            key: value for (key, value) in self.options().items()
            if key not in {'alpha', 'legend_color', 'legend_label', 'linestyle', 'rgbcolor', 'thickness'}
        }

        from matplotlib.path import Path
        fill_path = Path([(0, 0)])
        edge_path = Path([(0, 0)])

        from matplotlib.patches import PathPatch
        fill_patch = PathPatch(fill_path, **matplotlib_options)
        edge_patch = PathPatch(edge_path, **matplotlib_options)

        options = self.options()
        fill = options.pop('fill')
        if fill:
            subplot.axes.add_patch(fill_patch)
        subplot.axes.add_patch(edge_patch)

        # Translate SageMath options to matplotlib style.
        fill_patch.set_linewidth(float(options['thickness']))
        if 'linestyle' in options:
            fill_patch.set_linestyle(options['linestyle'])
        fill_patch.set_alpha(float(options['alpha']))

        from sage.plot.colors import to_mpl_color
        color = to_mpl_color(options.pop('rgbcolor'))

        fill_patch.set_fill(True)
        edge_patch.set_fill(False)

        edge_color = options.pop('edgecolor')
        if fill:
            if edge_color is None:
                fill_patch.set_color(color)
                edge_patch.set_color(color)
            else:
                fill_patch.set_facecolor(color)
                fill_patch.set_edgecolor(color)
                edge_patch.set_edgecolor(to_mpl_color(edge_color))
        else:
            edge_patch.set_edgecolor(edge_color or color)

        fill_patch.set_label(options['legend_label'])

        def redraw(_=None):
            r"""
            Redraw after the viewport has been rescaled to make sure that
            infinite rays reach the end of the viewport.
            """
            self._redraw_on_subplot(subplot, fill_patch, fill=True)
            self._redraw_on_subplot(subplot, edge_patch, fill=False)

        subplot.axes.callbacks.connect('ylim_changed', redraw)
        subplot.axes.callbacks.connect('xlim_changed', redraw)
        redraw()

    def _redraw_on_subplot(self, subplot, patch, fill):
        # Note that we throw RuntimeError (instead of NotImplementedError)
        # since some errors are silently consumed by matplotlib (or SageMath?)
        # it seems.

        from matplotlib.path import Path

        xlim = subplot.axes.get_xlim()
        ylim = subplot.axes.get_ylim()
        # TODO: Drop debug.
        # xlim = [-2, 10]
        # ylim = [-1, 3]

        commands = list(reversed(self._commands))

        command = commands.pop()

        def vertex(pos, direction):
            from sage.all import vector
            direction = vector(direction)
            pos = vector(pos)

            from sage.all import infinity
            if direction[0]:
                λx = max((xlim[0] - pos[0])/direction[0], (xlim[1] - pos[0])/direction[0])
            else:
                λx = infinity

            if direction[1]:
                λy = max((ylim[0] - pos[1])/direction[1], (ylim[1] - pos[1])/direction[1])
            else:
                λy = infinity

            λ = min(λx, λy)

            # Additionally, we now move out a full plot size so we are sure
            # that no artifacts coming from any sweeps (see below) show up in
            # the final plot.
            plot_size = (xlim[1] - xlim[0]) + (ylim[1] - ylim[0])
            λ += plot_size / direction.norm()

            return pos + λ * direction

        def extend(path):
            vertices.extend(path.vertices[1:])
            codes.extend(path.codes[1:])

        if command.code == "MOVETO":
            pos, = command.args
            direction = None
            vertices = [pos]
        elif command.code == "MOVETOINFINITY":
            pos, direction = command.args
            vertices = [vertex(pos, direction)]
        else:
            raise RuntimeError(f"path must not start with a {command.code} command")

        codes = [Path.MOVETO]

        while commands:
            command = commands.pop()

            if command.code == "LINETO":
                target, = command.args

                if direction is not None and pos != target:
                    raise ValueError(f"Cannot execute LINETO from infinite point at {pos} + λ {direction} when going to {target}")

                direction = None
                pos = target

                vertices.append(pos)
                codes.append(Path.LINETO)
            elif command.code == "LINETOINFINITY":
                assert direction is None

                base, direction = command.args
                # TODO: check up to epsilon
                # assert base == pos
                pos = base

                vertices.append(vertex(pos, direction))
                codes.append(Path.LINETO)
            elif command.code in "ARCTO":
                target, center = command.args

                assert direction is None

                extend(self._arc_path(center, pos, target))

                pos = target
            elif command.code == "RARCTO":
                target, center = command.args

                assert direction is None

                extend(self._arc_path(center, target, pos, reverse=True))

                pos = target
            elif command.code == "MOVETOINFINITY":
                assert direction is not None

                start = vertex(pos, direction)

                pos, direction = command.args
                end = vertex(pos, direction)

                # Sweep the bounding box counterclockwise from start to end
                from sage.all import vector
                # TODO: Is this the correct center?
                center = vector(((start[0] + end[0]) / 2, (start[1] + end[1]) / 2))

                extend(self._arc_path(center, start, end))

                vertices.append(end)
                codes.append(Path.LINETO if fill else Path.MOVETO)
            else:
                raise RuntimeError(f"cannot draw {command.code} yet")

        from matplotlib.path import Path
        patch.set_path(Path(vertices, codes))

    def get_minmax_data(self):
        try:
            from matplotlib.transforms import Bbox
            bbox = Bbox.null()

            pos = None
            for command in self._commands:
                if command.code in ["MOVETO", "LINETO"]:
                    pos, = command.args
                    bbox.update_from_data_xy([pos], ignore=False)
                elif command.code == "ARCTO":
                    target, center = command.args
                    bbox = bbox.union([bbox, self._arc_path(center, pos, target).get_extents()])
                    pos = target
                elif command.code == "RARCTO":
                    target, center = command.args
                    bbox = bbox.union([bbox, self._arc_path(center, target, pos, reverse=True).get_extents()])
                    pos = target
                elif command.code in ["LINETOINFINITY", "MOVETOINFINITY"]:
                    # TODO: Add a bit to the bounding box so these are always visible.
                    pass
                else:
                    raise NotImplementedError(f"cannot determine bounding box for {command.code} command")

            from sage.plot.plot import minmax_data

            return minmax_data(bbox.intervalx, bbox.intervaly, dict=True)
        except Exception as e:
            raise RuntimeError(e)

    @dataclass
    class Command:
        code: str
        args: tuple

    @classmethod
    def _arc_path(cls, center, start, end, reverse=False):
        from matplotlib.path import Path
        from math import atan2, pi
        from sage.all import vector
        # TODO: How many segments do we need?
        unit_arc = Path.arc(atan2(start[1] - center[1], start[0] - center[0]) / pi * 180, atan2(end[1] - center[1], end[0] - center[0]) / pi * 180, n=32)

        # Scale and translate the arc
        arc_vertices = unit_arc.vertices * vector((start[0] - center[0], start[1] - center[1])).norm() + center

        if reverse:
            arc_vertices = arc_vertices[::-1]

        return Path(arc_vertices, unit_arc.codes)

    @classmethod
    def hyperbolic_path(cls, commands, model, **kwds):
        if len(commands) < 2:
            raise ValueError("a path must contain at least two points")

        commands.reverse()

        command = commands.pop()
        if command.code == "MOVETO":
            pos, = command.args
            if model == "half_plane" and pos == pos.parent().infinity():
                next = commands[-1].args[0]
                bezier_commands = [BezierPath.Command("MOVETOINFINITY", [next.coordinates(model=model), (0, 1)])]
            else:
                from sage.all import RR
                bezier_commands = [BezierPath.Command("MOVETO", [pos.change_ring(RR).coordinates(model=model)])]
        else:
            raise ValueError("path must start with MOVETO or MOVETOINFINITY command")

        while commands:
            command = commands.pop()

            if command.code == "LINETO":
                target, = command.args
                bezier_commands.extend(cls._hyperbolic_segment(pos, target, model=model))
                pos = target
            elif command.code == "MOVETO":
                target, = command.args
                bezier_commands.extend(cls._hyperbolic_move(pos, target, model=model))
                pos = target
            else:
                raise NotImplementedError("unsuported hyperbolic plotting code")

        if bezier_commands[-1].code == "LINETOINFINITY" and bezier_commands[-1].args[-1] == (1, 0):
            assert bezier_commands[0].code == "MOVETOINFINITY" and bezier_commands[1].args[-1] == (0, 1)
            bezier_commands.append(bezier_commands[0])

        assert len(bezier_commands) >= 2

        return BezierPath(bezier_commands, kwds)

    @classmethod
    def _hyperbolic_segment(cls, start, end, model):
        if start == end:
            raise ValueError("cannot draw segment from point to itself")

        if model == "half_plane":
            if start == start.parent().infinity():
                return [
                        BezierPath.Command("MOVETOINFINITY", [end.coordinates(), (0, 1)]),
                        BezierPath.Command("LINETO", [end.coordinates()])]

            from sage.all import RR
            if end == end.parent().infinity():
                return [BezierPath.Command("LINETOINFINITY", [start.change_ring(RR).coordinates(model="half_plane"), (0, 1)])]

            p = start.change_ring(RR).coordinates(model="half_plane")
            q = end.change_ring(RR).coordinates(model="half_plane")

            if (p[0] - q[0]).abs() < (p[1] - q[1]).abs() * 1e-6:
                # This segment is (almost) vertical. We plot it as if it were
                # vertical to avoid numeric issus.
                return [BezierPath.Command("LINETO", [q])]

            geodesic = start.change_ring(RR).parent().geodesic(start, end)
            center = ((geodesic.start().coordinates()[0] + geodesic.end().coordinates()[0])/2, 0)

            return [BezierPath.Command("RARCTO" if p[0] < q[0] else "ARCTO", [q, center])]
        elif model == "klein":
            from sage.all import RR
            return [BezierPath.Command("LINETO", [end.change_ring(RR).coordinates(model="klein")])]
        else:
            raise NotImplementedError("cannot draw segment in this model")

    @classmethod
    def _hyperbolic_move(cls, start, end, model):
        if start == end:
            raise ValueError("cannot move from point to itself")

        if start.is_finite() or end.is_finite():
            raise ValueError("endpoints of move must be ideal")

        if model == "half_plane":
            if start == start.parent().infinity():
                from sage.all import RR
                return [
                    BezierPath.Command("MOVETOINFINITY", [end.change_ring(RR).coordinates(), (0, 1)]),
                    BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()]),
                ]

            if end == end.parent().infinity():
                from sage.all import RR
                return [BezierPath.Command("LINETOINFINITY", [start.change_ring(RR).coordinates(), (1, 0)])]

            from sage.all import RR
            if start.change_ring(RR).coordinates()[0] < end.change_ring(RR).coordinates()[0]:
                return [BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()])]
            else:
                return [
                    BezierPath.Command("LINETOINFINITY", [start.change_ring(RR).coordinates(), (1, 0)]),
                    BezierPath.Command("MOVETOINFINITY", [end.change_ring(RR).coordinates(), (-1, 0)]),
                    BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()]),
                ]

            raise NotImplementedError(f"cannot move from {start} to {end} in half plane model")
        elif model == "klein":
            # TODO: The actual arc should not be visible.
            from sage.all import RR
            return [BezierPath.Command("ARCTO", [end.change_ring(RR).coordinates(model="klein"), (0, 0)])]
        else:
            raise NotImplementedError("cannot move in this model")

    # def half_plane_segment(p, q):
    #     if p == p.parent().infinity():
    #         return reversed(half_plane_segment(q, p))

    #     if q == q.parent().infinity():
    #         from sage.all import infinity
    #         return [p, p, (infinity, 0, 1), (infinity, 0, 1)]

    #     p = p.coordinates(model="half_plane")
    #     q = q.coordinates(model="half_plane")

    #     if ((p[0] - q[0]) / (p[1] - q[1])).abs() < 1e-6:
    #         # This segment is (almost) vertical. We plot it as if it were
    #         # vertical to avoid numeric issus.
    #         return [p, q, p, q]

    #     raise NotImplementedError("cannot plot arcs yet")

    # def klein_segment(p, q):
    #     raise NotImplementedError

    # for p, q in zip(points, points[1:] + points[:1]):
    #     if p == q:
    #         raise ValueError("path must not contain repeated points")

    #     if model == "half_plane":
    #         control_points.extend(half_plane_segment(p, q))

    #     elif model == "klein":
    #         control_points.extend(klein_segment(p, q))

    #     else:
    #         raise NotImplementedError("cannot plot path in this model yet")

    # return bezier_path(control_points, **kwds)


@rename_keyword(color='rgbcolor')
@options(alpha=1, rgbcolor=(0, 0, 1), edgecolor=None, thickness=1, legend_label=None, legend_color=None, aspect_ratio=1.0, fill=True)
def hyperbolic_path(commands, model="half_plane", **options):
    if options["thickness"] is None:
        if options["fill"] and options["edgecolor"] is None:
            options["thickness"] = 0
        else:
            options["thickness"] = 1

    from sage.plot.all import Graphics

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))

    try:
        g.add_primitive(BezierPath.hyperbolic_path(commands[:], model=model, **options))
    except Exception as e:
        raise RuntimeError(f"Failed to render hyperbolic path {commands}", e)

    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    return g

# TODO: Ensure that every method has an INPUT section.
