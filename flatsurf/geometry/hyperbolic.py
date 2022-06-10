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
    sage: ab = H.geodesic(a, b)
    sage: ab
    {-x = 0}

Note that such a geodesic is oriented. The orientation is such that when we
replace the ``=`` in the above representation with a ``≥``, we obtain the half
space on its left::

    sage: H.geodesic(a, b).left_half_space()
    {x ≤ 0}

We can pass explicitly to the unoriented geodesic. Note that the oriented and
the unoriented version of a geodesic are not considered equal::

    sage: ab.unoriented()
    {x = 0}
    sage: ab == ab.unoriented()
    False
    sage: ab.is_subset(ab.unoriented())
    True
    sage: ab.unoriented().is_subset(ab)
    True

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

.. WARNING::

    Our implementation was not conceived with inexact rings in mind. Due to
    popular demand, we do allow inexact base rings but many operations have
    not been tuned for numerical stability, yet.

    To make our implementation work for a variety of (inexact) base rings, we
    delegate some numerically critical bits to a separate "geometry" class that
    can be specified when creating a hyperbolic plane. For exact base rings,
    this defaults to the :class:`HyperbolicExactGeometry` which uses exact text
    book algorithms.

    Over inexact rings, we implement a :class:`HyperbolicEpsilonGeometry` which
    considers two numbers two be equal if they differ only by a small
    (relative) error. This should work reasonably well for inexact base rings
    that have no denormalized numbers, i.e., it will work well for ``RR`` and
    general ``RealField``.

    sage: HyperbolicPlane(RR)
    Hyperbolic Plane over Real Field with 53 bits of precision

    There is currently no implementation that works well with ``RDF``. It
    should be easy to adapt :class:`HyperbolicEpsilonGeometry` for that purpose
    to take into account denormalized numbers::

        sage: HyperbolicPlane(RDF)
        Traceback (most recent call last):
        ...
        ValueError: geometry must be specified for HyperbolicPlane over inexact rings

    There is currently no implementation that works for ball arithmetic::

        sage: HyperbolicPlane(RBF)
        Traceback (most recent call last):
        ...
        ValueError: geometry must be specified for HyperbolicPlane over inexact rings

.. NOTE::

    This module implements different kinds of convex subsets as different
    classes. The alternative would have been to represent all subsets as
    collections of (in)equalities in some hyperbolic model. There is for
    example a :class:`HyperbolicUnorientedSegment` and a
    :class:`HyperbolicConvexPolygon` even though the former could in principle
    be expressed as the latter. The advantage of this approach is that we can
    provide a more natural user interface, e.g., a segment has a single
    underlying geodesic whereas the corresponding convex polygon would have
    four (or three.) Similarly, an oriented geodesic (which cannot really be
    expressed as a convex polygon due to the orientation) has a left and a
    right associated half spaces.

    Sometimes it can, however, be beneficial to treat each subset as a convex
    polygon. In such a case, one can explicitly create polygons from subsets by
    intersecting their :meth:`HyperbolicConvexSet.half_spaces`::

        sage: g = H.vertical(0)
        sage: P = H.polygon(g.half_spaces(), check=False, assume_minimal=True)
        sage: P
        {x ≤ 0} ∩ {x ≥ 0}

    Note that such an object might not be fully functional since some methods
    may assume that the object is an actual polygon::

        sage: P.dimension()
        2

    Similarly, a geodesic can be treated as a segment without endpoints::

        sage: H.segment(g, start=None, end=None, check=False, assume_normalized=True)
        {-x = 0}

.. NOTE::

    This implementation is an alternative to the one that comes with SageMath.
    The one in SageMath has a number of issues, see e.g.
    https://trac.sagemath.org/ticket/32400. The implementation here tries very
    hard to perform all operations over the same base ring, have the best
    complexities possible, keep all objects in the same (Klein) model, is not
    using any symbolic expressions, and tries to produce better plots.

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
from sage.misc.cachefunc import cached_method
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

    - ``geometry`` -- an implementation of the geometric primitives specified
      by :class:`HyperbolicExactGeometry`. If the ``base_ring`` is exact, this
      defaults to :class:`HyperbolicExactGeometry` over that base ring. If the
      base ring is ``RR`` or ``RealField``, this defaults to the
      :class:`HyperbolicEpsilonGeometry` over that ring. For other rings, a
      geometry must be explicitly provided.

    - ``category`` -- the category for this object; if not specified, defaults
      to sets. Note that we do not use metric spaces here since the elements
      are convex subsets of the hyperbolic plane and not just points and do
      therefore not satisfy the assumptions of a metric space.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

        sage: HyperbolicPlane()
        Hyperbolic Plane over Rational Field

    ::

        sage: HyperbolicPlane(AA)
        Hyperbolic Plane over Algebraic Real Field

    ::

        sage: HyperbolicPlane(RR)
        Hyperbolic Plane over Real Field with 53 bits of precision

    """

    @staticmethod
    def __classcall__(cls, base_ring=None, geometry=None, category=None):
        r"""
        Create the hyperbolic plane with normalized arguments to make it a
        unique SageMath parent.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicExactGeometry

            sage: HyperbolicPlane() is HyperbolicPlane(QQ)
            True

            sage: HyperbolicPlane() is HyperbolicPlane(geometry=HyperbolicExactGeometry(QQ))
            True

        """
        from sage.all import QQ

        base_ring = base_ring or QQ

        if geometry is None:
            from sage.all import RR
            if base_ring is RR:
                geometry = HyperbolicExactGeometry(QQ).change_ring(base_ring)
            elif base_ring.is_exact():
                geometry = HyperbolicExactGeometry(base_ring)
            else:
                raise ValueError("geometry must be specified for HyperbolicPlane over inexact rings")

        from sage.categories.all import Sets

        category = category or Sets()

        return super(HyperbolicPlane, cls).__classcall__(
            cls, base_ring=base_ring, geometry=geometry, category=category
        )

    def __init__(self, base_ring, geometry, category):
        r"""
        Create the hyperbolic plane over ``base_ring``.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: TestSuite(HyperbolicPlane(QQ)).run()
            sage: TestSuite(HyperbolicPlane(AA)).run()
            sage: TestSuite(HyperbolicPlane(RR)).run()

        """
        from sage.all import RR

        if geometry.base_ring() is not base_ring:
            raise ValueError(f"geometry base ring must be base ring of hyperbolic plane but {geometry.base_ring()} is not {base_ring}")

        if not RR.has_coerce_map_from(geometry.base_ring()):
            # We should check that the coercion is an embedding but this is not possible currently.
            raise ValueError("base ring must embed into the reals")

        super().__init__(category=category)
        self._base_ring = geometry.base_ring()
        self.geometry = geometry

    def change_ring(self, ring, geometry=None):
        r"""
        Return the hyperbolic plane over a different base ``ring``.

        INPUT:

        - ``ring`` -- a ring or ``None``; if ``None``, uses the current
          :meth:`base_ring`.

        - ``geometry`` -- a geometry or ``None``; if ``None``, tries to convert
          the existing geometry to ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane(QQ).change_ring(AA) is HyperbolicPlane(AA)
            True

        When changing to the ring ``RR`` and no geometry has been specified
        explicitly, the :class:`HyperbolicExactGeometry` changes to the
        :class:`HyperbolicEpsilonGeometry`, see
        :meth:`HyperbolicExactGeometry.change_ring`::

            sage: HyperbolicPlane(QQ).change_ring(RR) is HyperbolicPlane(RR)
            True

        In the opposite direction, the geometry cannot be determined automatically::

            sage: HyperbolicPlane(RR).change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: cannot change_ring() to an exact ring

        So the geometry has to be specified explicitly::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicExactGeometry
            sage: HyperbolicPlane(RR).change_ring(QQ, geometry=HyperbolicExactGeometry(QQ)) is HyperbolicPlane(QQ)
            True

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.change_ring` or more generally
            :meth:`HyperbolicConvexSet.change` to change the ring and geometry
            a set is defined over.

        """
        if ring is None and geometry is None:
            return self

        if ring is None:
            ring = self.base_ring()

        if geometry is None:
            geometry = self.geometry.change_ring(ring)

        return HyperbolicPlane(ring, geometry)

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

        elements = [
            self.empty_set(),
            # Points
            self.infinity(),
            self.real(0),
            self.real(1),
            self.real(-1),
            self.geodesic(0, 2).start(),
            # Oriented Geodesics
            self.vertical(1),
            self.half_circle(0, 1),
            self.half_circle(1, 3),
            # Unoriented Geodesics
            self.vertical(-1).unoriented(),
            self.half_circle(-1, 2).unoriented(),
            # Half spaces
            self.vertical(0).left_half_space(),
            self.half_circle(0, 2).left_half_space(),
        ]

        # The intersection algorithm is only implemented over exact rings.
        elements += [
            # An unbounded polygon
            self.vertical(1)
            .left_half_space()
            .intersection(self.vertical(-1).right_half_space()),
            # An unbounded polygon which is bounded in the Euclidean plane
            self.vertical(1)
            .left_half_space()
            .intersection(self.vertical(-1).right_half_space())
            .intersection(self.geodesic(0, 1).left_half_space())
            .intersection(self.geodesic(0, -1).right_half_space()),
            # A bounded polygon
            self.geodesic(-ZZ(1) / 3, 2)
            .left_half_space()
            .intersection(self.geodesic(ZZ(1) / 3, -2).right_half_space())
            .intersection(self.geodesic(-ZZ(2) / 3, 3).right_half_space())
            .intersection(self.geodesic(ZZ(2) / 3, -3).left_half_space()),
            # An unbounded oriented segment
            self.vertical(0).intersection(self.geodesic(-1, 1).left_half_space()),
            # A bounded oriented segment
            self.vertical(0)
            .intersection(self.geodesic(-2, 2).right_half_space())
            .intersection(self.geodesic(-ZZ(1) / 2, ZZ(1) / 2).left_half_space()),
            # An unbounded unoriented segment
            self.vertical(0)
            .intersection(self.geodesic(-1, 1).left_half_space())
            .unoriented(),
            # A bounded unoriented segment
            self.vertical(0)
            .intersection(self.geodesic(-2, 2).right_half_space())
            .intersection(self.geodesic(-ZZ(1) / 2, ZZ(1) / 2).left_half_space())
            .unoriented(),
        ]

        return elements

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

            TestSuite(x).run(
                verbose=tester._verbose,
                prefix=tester._prefix + "  ",
                raise_on_failure=is_sub_testsuite,
            )
            tester.info(tester._prefix + " ", newline=False)

    def random_element(self, kind=None):
        r"""
        Return a random convex subset of this hyperbolic plane.

        INPUT:

        - ``kind`` -- one of ``"empty_set"```, ``"point"```, ``"oriented geodesic"```,
          ``"unoriented geodesic"```, ``"half_space"```, ``"oriented segment"``,
          ``"unoriented segment"``, ``"polygon"``; the kind of set to produce.
          If not specified, the kind of set is chosen
          randomly.

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

            sage: H.random_element("oriented geodesic")
            {-12*(x^2 + y^2) + 1144*x + 1159 = 0}

            sage: H.random_element("unoriented geodesic")
            {648*(x^2 + y^2) + 1654*x + 85 = 0}

            sage: H.random_element("half_space")
            {3*(x^2 + y^2) - 5*x - 1 ≥ 0}

            sage: H.random_element("oriented segment")
            {-3*(x^2 + y^2) + x + 3 = 0} ∩ {9*(x^2 + y^2) - 114*x + 28 ≥ 0} ∩ {(x^2 + y^2) + 12*x - 1 ≥ 0}

            sage: H.random_element("unoriented segment")
            {16*(x^2 + y^2) - x - 16 = 0} ∩ {(x^2 + y^2) + 64*x - 1 ≥ 0} ∩ {496*(x^2 + y^2) - 1056*x + 529 ≥ 0}

            sage: H.random_element("polygon")
            {56766100*(x^2 + y^2) - 244977117*x + 57459343 ≥ 0} ∩ {822002048*(x^2 + y^2) - 3988505279*x + 2596487836 ≥ 0} ∩ {464*(x^2 + y^2) + 9760*x + 11359 ≥ 0} ∩ {4*(x^2 + y^2) + 45*x + 49 ≥ 0}

        .. SEEALSO::

            :meth:`some_elements` for a curated list of representative subsets.

        """
        kinds = ["empty_set", "point", "oriented geodesic", "unoriented geodesic", "half_space", "oriented segment", "unoriented segment", "polygon"]

        if kind is None:
            from sage.all import randint

            kind = kinds[randint(0, len(kinds) - 1)]

        if kind == "empty_set":
            set = self.empty_set()
        elif kind == "point":
            set = self.point(
                self.base_ring().random_element(),
                self.base_ring().random_element().abs(),
                model="half_plane",
                check=False,
            )
        elif "geodesic" in kind:
            a = self.random_element("point")
            b = self.random_element("point")
            while b == a:
                b = self.random_element("point")

            set = self.geodesic(a, b)
        elif kind == "half_space":
            set = self.random_element("geodesic").left_half_space()
        elif "segment" in kind:
            a = self.random_element("point")
            b = self.random_element("point")
            while a == b or (not a.is_finite() and not b.is_finite()):
                b = self.random_element("point")

            set = self.segment(self.geodesic(a, b), start=a, end=b)
        elif kind == "polygon":
            from sage.all import ZZ

            interior_points = [
                self.random_element("point")
                for i in range(ZZ.random_element().abs() + 3)
            ]

            half_spaces = []

            while len(half_spaces) < len(interior_points):
                half_space = self.random_element("half_space")

                for p in interior_points:
                    if p in half_space:
                        continue

                    a, b, c = half_space.equation(model="klein")

                    x, y = p.coordinates(model="klein")

                    a = -(b * x + c * y)

                    half_space = self.half_space(a, b, c, model="klein")

                    assert p in half_space

                half_spaces.append(half_space)

            set = self.polygon(half_spaces)
        else:
            raise ValueError(f"kind must be one of {kinds}")

        if "unoriented" in kind:
            set = set.unoriented()

        return set

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

        TESTS::

            sage: H(-I)
            Traceback (most recent call last):
            ...
            ValueError: point (0, -1) not in the upper half plane

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
                raise NotImplementedError(
                    "cannot create a hyperbolic point from an element in a number field that does not contain the imaginary unit"
                )

            return self.point(x.real(), x.imag(), model="half_plane")

        from sage.all import SR

        if parent is SR:
            return self.point(x.real(), x.imag(), model="half_plane")

        from sage.categories.all import Rings

        if parent in Rings():
            raise ValueError(
                f"cannot convert this element in {parent} to the hyperbolic plane over {self.base_ring()}"
            )

        raise NotImplementedError(
            "cannot create a subset of the hyperbolic plane from this element yet."
        )

    def base_ring(self):
        r"""
        Return the base ring over which objects in the plane are defined.

        More specifically, all geodesics must have an equation `a + bx + cy =
        0` in the Klein model with coefficients in this ring, and all points
        must have coordinates in this ring when written in the Klein model, or
        be the end point of a geodesic. All other objects are built from these
        primitives.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().base_ring()
            Rational Field

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.change_ring` to change the ring a set is defined over

        """
        return self._base_ring

    def infinity(self):
        r"""
        Return the point at infinity in the Poincaré half plane model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: p = H.infinity()
            sage: p
            ∞

            sage: p == H(oo)
            True
            sage: p.is_ideal()
            True

        .. SEEALSO::

            :meth:`point` to create points in general.

        """
        return self.projective(1, 0)

    def real(self, r):
        r"""
        Return the ideal point ``r`` on the real axis in the Poincaré half
        plane model.

        INPUT:

        - ``r`` -- an element of the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: p = H.real(-2)
            sage: p
            -2

            sage: p == H(-2)
            True
            sage: p.is_ideal()
            True

        .. SEEALSO::

            :meth:`point` to create points in general.

        """
        return self.projective(r, 1)

    def projective(self, p, q):
        r"""
        Return the ideal point with projective coordinates ``[p: q]`` in the
        upper half plane model.

        INPUT:

        - ``p`` -- an element of the :meth:`base_ring`.

        - ``q`` -- an element of the :meth:`base_ring`.

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

        .. SEEALSO::

            :meth:`point` to create points in general.

        """
        p = self.base_ring()(p)
        q = self.base_ring()(q)

        # TODO: Turn this into a proper predicate.
        if self.geometry.zero(p) and self.geometry.zero(q):
            raise ValueError("one of p and q must not be zero")

        # TODO: Turn this into a proper predicate.
        if self.geometry.zero(q):
            return self.point(0, 1, model="klein", check=False)

        return self.point(p / q, 0, model="half_plane", check=False)

    def start(self, geodesic):
        r"""
        Return the ideal starting point of ``geodesic``.

        INPUT:

        - ``geodesic`` -- an oriented geodesic

        .. NOTE::

            This method exists to keep all the methods that actually create
            hyperbolic sets on the lowest level in the
            :class:`HyperbolicPlane`. It is otherwise identical to
            :meth:`HyperbolicOrientedGeodesic.start`.

        EXAMPLES:

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H.start(H.vertical(0))
            sage: p
            0

            sage: H.vertical(0).start() == p
            True

        Points created this way might have coordinates that cannot be
        represented in the base ring::

            sage: p = H.half_circle(0, 2).start()
            sage: p.coordinates(model="klein")
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 ...

            sage: p.coordinates(model="half_plane")
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 ...

        """
        geodesic = self(geodesic)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("geodesic must be an oriented geodesic")

        return self.__make_element_class__(HyperbolicPointFromGeodesic)(self, geodesic)

    def point(self, x, y, model, check=True):
        r"""
        Return the point with coordinates (x, y) in the given model.

        When ``model`` is ``"half_plane"``, return the point `x + iy` in the upper half plane.

        When ``model`` is ``"klein"``, return the point (x, y) in the Klein model.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`

        - ``y`` -- an element of the :meth:`base_ring`

        - ``model`` -- one of ``"half_plane"`` or ``"klein"``

        - ``check`` -- whether to validate the inputs (default: ``True``); set
          this to ``False``, to create an ultra-ideal point, i.e., a point
          outside the unit circle in the Klein model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.point(0, 1, model="half_plane")
            I

            sage: H.point(1, 2, model="half_plane")
            1 + 2*I

            sage: H.point(0, 1, model="klein")
            ∞

        An ultra-ideal point::

            sage: H.point(2, 3, model="klein")
            Traceback (most recent call last):
            ...
            ValueError: point (2, 3) is not in the unit disk in the Klein model

            sage: H.point(2, 3, model="klein", check=False)
            (2, 3)

        .. SEEALSO::

            :meth:`HyperbolicOrientedGeodesic.start` and
            :meth:`HyperbolicOrientedGeodesic.end` to generate points that do
            not have coordinates over the base ring.
            :meth:`infinity`, :meth:`real`, and :meth:`projective` as shortcuts
            to generate ideal ponits.

        """
        x = self.base_ring()(x)
        y = self.base_ring()(y)

        if model == "klein":
            point = self.__make_element_class__(HyperbolicPointFromCoordinates)(self, x, y)
        elif model == "half_plane":
            # TODO: Turn this into a proper predicate.
            if self.geometry.sgn(y) < 0:
                raise ValueError(f"point {x, y} not in the upper half plane")

            denominator = 1 + x * x + y * y
            return self.point(
                x=2 * x / denominator,
                y=(-1 + x * x + y * y) / denominator,
                model="klein",
                check=check,
            )
        else:
            raise NotImplementedError("unsupported model")

        if check:
            point._check()

        return point

    def half_circle(self, center, radius_squared):
        # TODO: Check documentation.
        # TODO: Check for doctests
        r"""
        Return the geodesic centered around the real ``center`` and with
        ``radius_squared`` in the Poincaré half plane model. The geodesic is
        oriented such that the point at infinity is to its left.

        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.

        INPUT:

        - ``center`` -- an element of the :meth:`base_ring`, the center of the
          half circle on the real axis

        - ``radius_squared`` -- a positive element of the :meth:`base_ring`,
          the square of the radius of the half circle

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

        .. SEEALSO::

            :meth:`vertical` to get an oriented vertical in the
            half plane model and :meth:`geodesic` for the most general
            interface, producing a geodesic from an equation.

        """
        center = self.base_ring()(center)
        radius_squared = self.base_ring()(radius_squared)

        # TODO: Turn this into a proper predicate.
        if self.geometry.sgn(radius_squared) <= 0:
            raise ValueError("radius must be positive")

        # Represent this geodesic as a(x^2 + y^2) + b*x + c = 0
        a = 1
        b = -2 * center
        c = center * center - radius_squared

        return self.geodesic(a, b, c, model="half_plane")

    def vertical(self, real):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

    def geodesic(self, a, b, c=None, model=None, oriented=True, check=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
            ValueError: square root of 32 ...

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

            return self.geodesic(A, B, C, model="klein", oriented=oriented, check=check)

        if model is None:
            raise ValueError(
                "a model must be specified when specifying a geodesic with coefficients"
            )

        if model == "half_plane":
            # Convert to the Klein model.
            return self.geodesic(
                a + c, b, a - c, model="klein", oriented=oriented, check=check
            )

        if model == "klein":
            a = self.base_ring()(a)
            b = self.base_ring()(b)
            c = self.base_ring()(c)

            geodesic = self.__make_element_class__(
                HyperbolicOrientedGeodesic if oriented else HyperbolicUnorientedGeodesic
            )(self, a, b, c)

            if check:
                geodesic = geodesic._normalize()
                geodesic._check()

            return geodesic

        raise NotImplementedError(
            "cannot create geodesic from coefficients in this model"
        )

    def half_space(self, a, b, c, model, check=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

    # TODO: Look into orientation. When there is start & end then the result is
    # oriented. When there is only start or end, then geodesic must be oriented
    # an the result is oriented. If there is neither start nor end, we also
    # require geodesic to be oriented.
    def segment(
        self,
        geodesic,
        start=None,
        end=None,
        oriented=None,
        check=True,
        assume_normalized=False,
    ):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the segment on the ``geodesic`` bounded by ``start`` and ``end``.

        INPUT:

        - ``geodesic`` -- a :meth:`geodesic` in this space.

        - ``start`` -- ``None`` or a :meth:`point` on the ``geodesic``, e.g.,
          obtained from the :meth:`HyperbolicOrientedGeodesic.intersection` of
          ``geodesic`` with another geodesic. If ``None``, the segment starts
          at the infinite :meth:`HyperbolicOrientedGeodesic.start` point of the
          geodesic.

        - ``end`` -- ``None`` or a :meth:`point` on the ``geodesic``, same as
          ``start``; must be later on ``geodesic`` than ``start``.

        - ``oriented`` -- whether to produce an oriented segment or an
          unoriented segment. The default (``None``) is to produce an oriented
          segment iff ``geodesic`` is oriented or both ``start`` and ``end``
          are provided so the orientation can be deduced from their order.

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
            ValueError: square root of 32 ...

        The produced segment is oriented if the ``geodesic`` is oriented::

            sage: H.segment(H.vertical(0)).is_oriented()
            True

            sage: H.segment(H.vertical(0).unoriented()).is_oriented()
            False

        The segment is oriented if both ``start`` and ``end`` are provided::

            sage: H.segment(H.vertical(0).unoriented(), start=0, end=oo).is_oriented()
            True
            sage: H.segment(H.vertical(0).unoriented(), start=2*I, end=I).is_oriented()
            True
            sage: H.segment(H.vertical(0).unoriented(), start=I) != H.segment(H.vertical(0).unoriented(), end=I)
            True

        .. SEEALSO::

            :meth:`HyperbolicPoint.segment`

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

        if oriented is None:
            oriented = geodesic.is_oriented() or (start is not None and end is not None)

        if not geodesic.is_oriented():
            geodesic = geodesic.change(oriented=True)

        segment = self.__make_element_class__(
            HyperbolicOrientedSegment if oriented else HyperbolicUnorientedSegment
        )(self, geodesic, start, end)

        if check:
            segment._check(require_normalized=False)

        if not assume_normalized:
            segment = segment._normalize()

        if check:
            segment._check(require_normalized=True)

        return segment

    def polygon(
        self, half_spaces, check=True, assume_sorted=False, assume_minimal=False
    ):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
            sage: H.polygon(minimal.half_spaces(), check=False, assume_minimal=True)
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

            sage: H.polygon(H.infinity().half_spaces(), assume_sorted=True)
            ∞

        """
        half_spaces = [self.coerce(half_space) for half_space in half_spaces]

        half_spaces = HyperbolicHalfSpaces(half_spaces, assume_sorted=assume_sorted)

        polygon = self.__make_element_class__(HyperbolicConvexPolygon)(
            self, half_spaces
        )

        if check:
            polygon._check(require_normalized=False)

        if check or not assume_minimal:
            polygon = polygon._normalize()

        if check:
            polygon._check()

        return polygon

    def intersection(self, *subsets):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the intersection of convex ``subsets``.

        ALGORITHM:

        We compute the intersection of the
        :meth:`HyperbolicConvexSet.half_spaces` that make up the ``subsets``.
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

        if len(subsets) == 0:
            raise NotImplementedError(
                "the full hyperbolic space cannot be created as an intersection"
            )

        if len(subsets) == 1:
            return subsets[0].unoriented()

        half_spaces = sum(
            [subset.half_spaces() for subset in subsets], HyperbolicHalfSpaces([])
        )

        return self.polygon(
            half_spaces, assume_sorted=True, assume_minimal=False, check=False
        ).unoriented()

    def empty_set(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return an empty subset of this space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().empty_set()
            {}

        """
        return self.__make_element_class__(HyperbolicEmptySet)(self)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a printable representation of this hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: HyperbolicPlane(AA)
            Hyperbolic Plane over Algebraic Real Field

        """
        return f"Hyperbolic Plane over {repr(self.base_ring())}"


# TODO: Define more "hyperbolic" predicates instead.
class HyperbolicGeometry:
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    def __init__(self, ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        self._ring = ring

    def base_ring(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self._ring

    def zero(self, x):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.cmp(x, 0) == 0

    # TODO: Do not mix the vector case into this.
    def equal(self, x, y):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        try:
            return all(self._equal(a, b) for a, b in zip(x, y))
        except TypeError:
            return self._equal(x, y)

    def cmp(self, x, y):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self._equal(x, y):
            return 0
        if x < y:
            return -1

        assert x > y
        return 1

    def sgn(self, x):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.cmp(x, 0)

    def _equal(self, x, y):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        raise NotImplementedError

    def change_ring(ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        raise NotImplementedError


class HyperbolicExactGeometry(UniqueRepresentation, HyperbolicGeometry):

    def _equal(self, x, y):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return x == y

    def change_ring(self, ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import RR
        if ring is RR:
            return HyperbolicEpsilonGeometry(ring, 1e-6)

        if not ring.is_exact():
            raise ValueError("cannot change_ring() to an inexact ring")

        return HyperbolicExactGeometry(ring)


class HyperbolicEpsilonGeometry(UniqueRepresentation, HyperbolicGeometry):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    def __init__(self, ring, epsilon):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(ring)
        self._epsilon = ring(epsilon)

    def _equal(self, x, y):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # see https://floating-point-gui.de/errors/comparison/ but we do not
        # need to worry about denormalized numbers in RR; apparently the
        # exponent can get arbitrarily large there.
        # TODO: Test that this really does what we want it to do.
        if x == 0 or y == 0:
            return abs(x - y) < self._epsilon

        return abs(x - y) <= (abs(x) + abs(y)) * self._epsilon

    def change_ring(self, ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if ring.is_exact():
            raise ValueError("cannot change_ring() to an exact ring")

        return HyperbolicEpsilonGeometry(ring, self._epsilon)


# TODO: Change richcmp to match the description below.
class HyperbolicConvexSet(Element):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    Base class for convex subsets of :class:`HyperbolicPlane`.

    .. NOTE::

        Concrete subclasses should apply the following rules.

        There should only be a single type to describe a certain subset:
        normally, a certain subset, say a point, should only be described by a
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

        Comparison with ``==`` should mean "is essentially indistinguishable from":
        Implementing == to mean anything else would get us into trouble in the long run. In
        particular we cannot implement <= to mean "is subset of" since then an
        oriented and an unoriented geodesic would be `==`. So, objects of a
        different type are never equal. We do, however, treat objects as equal
        that only differ in their exact representation such as the geodesic x =
        1 and the geodesic 2x = 2.

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicConvexSet
        sage: H = HyperbolicPlane()

        sage: isinstance(H(0), HyperbolicConvexSet)
        True

    """

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a minimal set of half spaces whose intersection is this convex set.

        The half spaces are ordered by :meth:`HyperbolicHalfSpace._less_than`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).left_half_space().half_spaces()
            {{x ≤ 0},}

            sage: H.vertical(0).half_spaces()
            {{x ≤ 0}, {x ≥ 0}}

            sage: H(0).half_spaces()
            {{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}}

        """
        raise NotImplementedError(f"{type(self)} does not implement half_spaces()")

    def _test_half_spaces(self, **options):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Verify that this convex set implements :meth:`half_spaces` correctly.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_half_spaces()

        """
        tester = self._tester(**options)

        half_spaces = self.half_spaces()

        tester.assertEqual(self.parent().intersection(*half_spaces), self.unoriented())

        tester.assertTrue(isinstance(half_spaces, HyperbolicHalfSpaces))

        for a, b in zip(list(half_spaces), list(half_spaces)[1:]):
            tester.assertTrue(HyperbolicHalfSpaces._lt_(a, b))

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Validate this convex subset.

        If ``require_normalized``, we also check that the object has the
        correct implementation class, e.g., that a point is a
        :class:`HyperbolicPoint` and not say a
        :class:`HyperbolicOrientedSegment` of length zero.
        """
        pass

    def _normalize(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return this set possibly rewritten in a simpler form.

        This method is only relevant for sets created with ``check=False``.
        Such sets might have been created in a non-canonical way, e.g., when
        creating a :class:`HyperbolicOrientedSegment` whose start and end point are ideal,
        then this is actually a geodesic and it shuold be described as such.
        """
        return self

    def unoriented(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the non-oriented version of this set.

        Some sets such as geodesics and segments can have an explicit
        orientation. This method returns the underlying set without any
        explicit orientation.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.vertical(0).unoriented()
            {x = 0}

        """
        return self.change(oriented=False)

    def _test_unoriented(self, **options):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Verify that :meth:`unoriented` is implemented correctly.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_unoriented()

        """
        tester = self._tester(**options)

        tester.assertEqual(self.unoriented(), self.unoriented().unoriented())

    def intersection(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the intersection with the ``other`` convex set.
        """
        return self.parent().intersection(self, other)

    def __contains__(self, point):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return whether ``point`` is contained in this set.
        """
        for half_space in self.half_spaces():
            if point not in half_space:
                return False

        return True

    def is_finite(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return whether all points in this set are finite.
        """
        # TODO: This could be implemented generically.
        raise NotImplementedError(
            f"this {type(self)} does not support checking finiteness"
        )

    def change_ring(self, ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return this set as an element of the hyperbolic plane over ``ring``.
        """
        return self.change(ring=ring)

    def _test_change_ring(self, **options):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Verify that this set implements :meth:`change_ring`.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_change_ring()

        """
        tester = self._tester(**options)
        tester.assertEqual(self, self.change_ring(self.parent().base_ring()))

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a modified copy of this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: geodesic = H.geodesic(0, 1)

        We can change the base ring over which this set is defined::

            sage: geodesic.change(ring=AA)
            {(x^2 + y^2) - x = 0}

        We can drop the explicit orientation of a set::

            sage: geodesic.change(oriented=False)
            {(x^2 + y^2) - x = 0}

        """
        raise NotImplementedError(f"this {type(self)} does not implement change()")

    def _test_change(self, **options):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Verify that the full interface of :meth:`change` has been implemented.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.an_element()._test_change()

        """
        tester = self._tester(**options)

        # The ring parameter is supported
        tester.assertEqual(self, self.change(ring=self.parent().base_ring()))

        # The oriented parameter is supported
        tester.assertEqual(
            self.change(oriented=False),
            self.change(oriented=False).change(oriented=False),
        )
        if self != self.change(oriented=False):
            tester.assertEqual(self, self.change(oriented=True))

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a plot of this subset.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).plot()
            Graphics object consisting of 1 graphics primitive

        """
        # TODO: This could be implemented generically.
        raise NotImplementedError(f"this {type(self)} does not support plotting")

    def _test_plot(self, **options):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Turn this into a proper predicate.
        predicates = self.parent().geometry

        # TODO: check that isometry is actually a matrix over the right ring?
        if (
            isometry.nrows() != 3
            or isometry.ncols() != 3
            or not self.parent().base_ring().has_coerce_map_from(isometry.base_ring())
        ):
            raise ValueError("invalid isometry")

        from sage.matrix.special import diagonal_matrix

        D = isometry.transpose() * diagonal_matrix([1, 1, -1]) * isometry
        for i, row in enumerate(D):
            for j, entry in enumerate(row):
                if (i == j) == predicates.zero(entry):
                    raise ValueError("invalid isometry")

        if not predicates.equal(D[0, 0], D[1, 1]) or not predicates.equal(
            D[0, 0], -D[2, 2]
        ):
            raise ValueError("invalid isometry")

    def _apply_isometry_klein(self, isometry):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: This can be implemented generically.
        raise NotImplementedError

    def apply_isometry(self, isometry, model="half_plane"):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

        raise NotImplementedError(
            "applying isometry not supported in this hyperbolic model"
        )

    def _acted_upon_(self, x, self_on_left):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the result of acting upon this set with ``x``.

        EXAMPLES:

        The Möbius transformation that sends `z` to `(1 + 2z)/(3 + 4z)`::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: p = H(I)
            sage: matrix([[1, 2], [3, 4]]) * p
            11/25 + 2/25*I

        """
        return self.apply_isometry(x)

    def is_subset(self, other):
        r"""
        Return whether this set is a subset of ``other``.

        INPUT:

        - ``other`` -- another hyperbolic convex set

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H(I).is_subset(H.vertical(0))
            True
            sage: H.vertical(0).is_subset(H(I))
            False

        """
        if self.is_empty():
            return True

        other = self.parent()(other)

        if self.dimension() > other.dimension():
            return False

        vertices = self.vertices()

        for vertex in vertices:
            if vertex not in other:
                return False

        if self.dimension() > 1:
            raise NotImplementedError("cannot decide subset relation for this two-dimensional set")

        return True

    # TODO: Test that is_subset can compare all kinds of sets by inclusion.

    # TODO: Provide hashing.

    def an_element(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Test that everything implements an_element().
        return next(iter(self.vertices()))

    @classmethod
    def _enhance_plot(self, plot, model):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if model == "klein":
            from sage.all import circle

            plot = circle([0, 0], 1, fill=False, color="#d1d1d1", zorder=-1) + plot

        return plot

    def is_empty(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.dimension() < 0

    def __bool__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return not self.is_empty()

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Test that this is a ZZ integer.
        raise NotImplementedError(f"{type(self)} does not implement dimension() yet")

    def is_point(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.dimension() == 0

    def is_oriented(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Test that oriented sets implement _neg_ correctly.
        return isinstance(self, HyperbolicOrientedConvexSet)

    def edges(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Implement by iterating vertices.
        # TODO: Test that other implementation are consistent with vertices()
        raise NotImplementedError


class HyperbolicOrientedConvexSet(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        raise NotImplementedError


class HyperbolicHalfSpace(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A closed half space of the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)

        sage: H.half_circle(0, 1).left_half_space()
        {(x^2 + y^2) - 1 ≥ 0}

    """

    def __init__(self, parent, geodesic):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("geodesic must be an oriented geodesic")

        self._geodesic = geodesic

    def equation(self, model, gcd=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return an inequality for this half space as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is in the half space if it satisfies `a(x^2 + y^2) + bx + c \ge 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy \ge 0`.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.
        """
        return self._geodesic.equation(model=model, gcd=gcd)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
        # TODO: Turn this into a proper predicate.
        sgn = self.parent().geometry.sgn

        # Convert to the Poincaré half plane model as a(x^2 + y^2) + bx + c ≥ 0.
        a, b, c = self.equation(model="half_plane", gcd=None)

        # Remove any trailing - signs in the output.
        cmp = "≥"
        if sgn(a) < 0 or (sgn(a) == 0 and sgn(b) < 0):
            a *= -1
            b *= -1
            c *= -1
            cmp = "≤"

        from sage.all import PolynomialRing

        R = PolynomialRing(self.parent().base_ring(), names="x")
        if sgn(a) != 0:
            return (
                f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} {cmp} 0}}"
            )
        else:
            return f"{{{repr(R([c, b]))} {cmp} 0}}"

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: S = H.vertical(0).left_half_space()
            sage: [S] == list(S.half_spaces())
            True

        """
        return HyperbolicHalfSpaces([self], assume_sorted=True)

    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

    def boundary(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Turn this into a proper predicate.
        sgn = self.parent().geometry.sgn

        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        return sgn(a + b * x + c * y) >= 0

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicHalfSpace):
                return False
            return self._geodesic._richcmp_(other._geodesic, op)

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return (
            self.parent()
            .polygon([self], check=False, assume_minimal=True)
            .plot(model=model, **kwds)
        )

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if ring is not None or geometry is not None:
            self = self._geodesic.change(ring=ring, geometry=geometry).left_half_space()

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of half space")

        return self

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ(2)

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.boundary().vertices()


class HyperbolicGeodesic(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Change representation:
    # I would prefer geodesics to be encoded as dual to the quadratic form.
    # That is v = (a, b, c) should represent the half space {x in R^3 :  B(v,
    # x) >= 0} where B is the bilinear form associated to Q. This encoding
    # makes a much smoother transition between points, ideal points and
    # geodesics : they are all represented by elements of R^3 (points are Q >
    # 0, ideal points are Q = 0 and geodesics are Q < 0 (ultra-ideal points)).
    # See
    # https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/hyperbolic.20geometry/near/284722650
    # and implement what's written there.
    r"""
    A geodesic in the hyperbolic plane.

    This is the abstract base class of :class:`HyperbolicUnorientedGeodesic`
    and :class:`HyperbolicOrientedGeodesic`.

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicGeodesic
        sage: H = HyperbolicPlane()

        sage: geodesic = H.vertical(0)

        sage: isinstance(geodesic, HyperbolicGeodesic)
        True

        sage: isinstance(geodesic.unoriented(), HyperbolicGeodesic)
        True

    """

    def __init__(self, parent, a, b, c):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

    def __hash__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: This is wrong when there is a GCD. What about inexact rings?
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: hash(H.vertical(0)) == hash(H.vertical(0))
            True

        """
        return hash((self._a, self._b, self._c))

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # Convert to the Poincaré half plane model as a(x^2 + y^2) + bx + c = 0.
        a, b, c = self.equation(model="half_plane", gcd=None)

        from sage.all import PolynomialRing

        R = PolynomialRing(self.parent().base_ring(), names="x")
        # TODO: Turn this into a proper predicate.
        sgn = self.parent().geometry.sgn
        if sgn(a) != 0:
            return f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} = 0}}"
        else:
            return f"{{{repr(R([c, b]))} = 0}}"

    def equation(self, model, gcd=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return an equation for this geodesic as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is on the geodesic if it satisfies `a(x^2 + y^2) + bx + c = 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy = 0`.

        INPUT:

        - ``model`` -- the model in which this equation holds, either
          ``"half_plane"`` or ``"klein"``

        - ``gcd`` -- whether to divide the coefficients in the equation by
          their gcd. If ``None``, such a division is attempted but errors are
          silently ignored.

        If this geodesic :meth;`is_oriented`, then the sign of the coefficients
        is chosen to encode the orientation of this geodesic. The sign is such
        that the half plane obtained by replacing ``=`` with ``≥`` in above
        equationsis on the left of the geodesic.

        Note that the output does not uniquely describe the geodesic since the
        coefficients are only unique up to scaling.

        """
        a, b, c = self._a, self._b, self._c

        if model == "klein":
            a, b, c = a, b, c
        elif model == "half_plane":
            a, b, c = a + c, 2 * b, a - c
        else:
            raise NotImplementedError("cannot determine equation for this model yet")

        if gcd is not False:
            try:
                from sage.all import gcd as gcd_

                d = gcd_((a, b, c))
                assert d > 0
                a /= d
                b /= d
                c /= d
            except Exception:
                if gcd:
                    raise

        if not self.is_oriented():
            # TODO: Turn this into a proper predicate.
            sgn = self.parent().geometry.sgn
            if (
                sgn(a) < 0
                or (sgn(a) == 0 and b < 0)
                or (sgn(a) == 0 and sgn(b) == 0 and sgn(c) < 0)
            ):
                a *= -1
                b *= -1
                c *= -1

        return a, b, c

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H.vertical(0).half_spaces()
            {{x ≤ 0}, {x ≥ 0}}

        """
        self = self.change(oriented=True)
        return HyperbolicHalfSpaces([self.left_half_space(), self.right_half_space()])

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Create a plot of this geodesic in the hyperbolic ``model``.

        Additional arguments are passed on to the underlying SageMath plotting methods.
        """
        return (
            self.parent()
            .segment(
                self.change(oriented=True),
                start=None,
                end=None,
                check=False,
                assume_normalized=True,
            )
            .plot(model=model, **kwds)
        )

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            # TODO: Turn this into a proper predicate.
            equal = self.parent().geometry.equal
            sgn = self.parent().geometry.sgn
            if type(self) is not type(other):
                return False
            if sgn(self._b):
                return (
                    (not self.is_oriented() or sgn(self._b) == sgn(other._b))
                    and equal(self._a * other._b, other._a * self._b)
                    and equal(self._c * other._b, other._c * self._b)
                )
            else:
                return (
                    (not self.is_oriented() or sgn(self._c) == sgn(other._c))
                    and equal(self._a * other._c, other._a * self._c)
                    and equal(self._b * other._c, other._b * self._c)
                )

        super()._richcmp_(other, op)

    def __contains__(self, point):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        # TODO: Turn this into a proper predicate.
        zero = self.parent().geometry.zero
        return zero(a + b * x + c * y)

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ(1)

    def change(self, *, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a modified copy of this geodesic.

        EXAMPLES:

        The base ring over which this geodesic is defined can be changed::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(1).change_ring(QQ)
            {-x + 1 = 0}

            sage: H.vertical(AA(2).sqrt()).change(ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real ... to Rational

        We can forget the orientation of a geodesic::

            sage: v = H.vertical(0)
            sage: v.is_oriented()
            True
            sage: v = v.change(oriented=False)
            sage: v.is_oriented()
            False

        We can (somewhat randomly) pick the orientation of a geodesic::

            sage: v = v.change(oriented=True)
            sage: v.is_oriented()
            True

        """
        if ring is not None or geometry is not None:
            self = self.parent().change_ring(ring, geometry=geometry).geodesic(
                self._a,
                self._b,
                self._c,
                model="klein",
                check=False,
                oriented=self.is_oriented(),
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            self = self.parent().geodesic(
                self._a, self._b, self._c, model="klein", check=False, oriented=oriented
            )

        return self


class HyperbolicUnorientedGeodesic(HyperbolicGeodesic):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    An unoriented geodesic in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.vertical(0).unoriented()
        {x = 0}

    """

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.change(oriented=True).vertices()


class HyperbolicOrientedGeodesic(HyperbolicGeodesic, HyperbolicOrientedConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
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

    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the reversed geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: -H.vertical(0)
            {x = 0}

        """
        return self.parent().geodesic(
            -self._a, -self._b, -self._c, model="klein", check=False
        )

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        if self.is_ultra_ideal():
            raise ValueError(
                f"equation {self._a} + ({self._b})*x + ({self._c})*y = 0 does not define a chord in the Klein model"
            )

    def is_ultra_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        # TODO: Turn this into a proper predicate.
        cmp = self.parent().geometry.cmp
        return cmp(self._b * self._b + self._c * self._c, self._a * self._a) <= 0

    def start(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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
            ValueError: square root of 32 ...

        Passing to a bigger field, the coordinates can be represented::

            sage: K.<a> = QQ.extension(x^2 - 2, embedding=1.4)
            sage: H.half_circle(0, 2).change_ring(K).start()
            -a

        """
        return self.parent().start(self)

    def end(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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
            ValueError: square root of 32 ...

        Passing to a bigger field, the coordinates can be represented::

            sage: K.<a> = QQ.extension(x^2 - 2, embedding=1.4)
            sage: H.half_circle(0, 2).change_ring(K).end()
            a

        """
        return (-self).start()

    def left_half_space(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        Return the closed half space to the right of this (oriented) geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(0).right_half_space()
            {x ≥ 0}

        """
        return (-self).left_half_space()

    def _configuration(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        Return a classification of the angle between this
        geodesic and ``other`` in the Klein model.

        """
        # TODO: Can we make this public somehow?
        intersection = self._intersection(other)

        if intersection is None:
            sgn = self.parent().geometry.sgn
            orientation = sgn(self._b * other._b + self._c * other._c)

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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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
        if not isinstance(other, HyperbolicOrientedGeodesic):
            raise TypeError("can only intersect with another oriented geodesic")

        # TODO: Reference the trac ticket that says that solving for 2×2 matrices is very slow.
        det = self._b * other._c - self._c * other._b
        # TODO: Turn this into a proper predicate.
        zero = self.parent().geometry.zero

        if zero(det):
            return None

        x = (-other._c * self._a + self._c * other._a) / det
        y = (other._b * self._a - self._b * other._a) / det

        return self.parent().point(x, y, model="klein", check=False)

    def an_element(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return (
            self.parent()
            .geodesic(0, -self._c, self._b, model="klein", check=False)
            ._intersection(self)
        )

    def parametrize(self, point, model, check=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        if isinstance(point, HyperbolicPoint):
            if check and point not in self:
                raise ValueError("point must be on geodesic to be parametrized")

        if model == "euclidean":
            base = self.an_element().coordinates(model="klein")
            tangent = (self._c, -self._b)

            if isinstance(point, HyperbolicPoint):
                # TODO: Turn this into a proper predicate.
                zero = self.parent().geometry.zero
                coordinate = 0 if not zero(tangent[0]) else 1
                return (
                    point.coordinates(model="klein")[coordinate] - base[coordinate]
                ) / tangent[coordinate]

            λ = self.parent().base_ring()(point)

            return self.parent().point(
                x=base[0] + λ * tangent[0],
                y=base[1] + λ * tangent[1],
                model="klein",
                check=check,
            )

        raise NotImplementedError("cannot parametrize a geodesic over this model yet")

    def _apply_isometry_klein(self, isometry):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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

        b, c, a = (
            vector(self.parent().base_ring(), [self._b, self._c, self._a])
            * isometry.inverse()
        )
        return self.parent().geodesic(a, b, c, model="klein")

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return HyperbolicVertices([self.start(), self.end()])


class HyperbolicPoint(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A (possibly infinite) point in the :class:`HyperbolicPlane`.

    Internally, we typically represent a point as the Euclidean coordinates in
    the unit disk of the Klein model.

    Additionally, we allow points that are the (ideal) endpoints of geodesics
    even if these only have coordinates over a quadratic extension.
    """

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self.is_ultra_ideal():
            raise ValueError(f"point {self.coordinates(model='klein')} is not in the unit disk in the Klein model")

    def is_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        x, y = self.coordinates(model="klein")
        # TODO: Turn this into a proper predicate.
        return self.parent().geometry.cmp(x * x + y * y, 1) == 0

    def is_ultra_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        x, y = self.coordinates(model="klein")
        # TODO: Turn this into a proper predicate.
        return self.parent().geometry.cmp(x * x + y * y, 1) > 0

    def is_finite(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        x, y = self.coordinates(model="klein")
        # TODO: Turn this into a proper predicate.
        cmp = self.parent().geometry.cmp
        return cmp(x * x + y * y, 1) < 0

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H(I).half_spaces()
            {{(x^2 + y^2) + 2*x - 1 ≤ 0}, {x ≥ 0}, {(x^2 + y^2) - 1 ≥ 0}}

            sage: H(I + 1).half_spaces()
            {{x - 1 ≤ 0}, {(x^2 + y^2) - 3*x + 1 ≤ 0}, {(x^2 + y^2) - 2 ≥ 0}}

            sage: H.infinity().half_spaces()
            {{x ≤ 0}, {x - 1 ≥ 0}}

            sage: H(0).half_spaces()
            {{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}}

            sage: H(-1).half_spaces()
            {{x + 1 ≤ 0}, {(x^2 + y^2) - 1 ≤ 0}}

            sage: H(1).half_spaces()
            {{(x^2 + y^2) - x ≤ 0}, {(x^2 + y^2) - 1 ≥ 0}}

            sage: H(2).half_spaces()
            {{2*(x^2 + y^2) - 3*x - 2 ≥ 0}, {3*(x^2 + y^2) - 7*x + 2 ≤ 0}}

            sage: H(-2).half_spaces()
            {{(x^2 + y^2) - x - 6 ≥ 0}, {2*(x^2 + y^2) + 3*x - 2 ≤ 0}}

            sage: H(1/2).half_spaces()
            {{6*(x^2 + y^2) - x - 1 ≤ 0}, {2*(x^2 + y^2) + 3*x - 2 ≥ 0}}

            sage: H(-1/2).half_spaces()
            {{2*(x^2 + y^2) + 7*x + 3 ≤ 0}, {2*(x^2 + y^2) - 3*x - 2 ≤ 0}}

        For ideal endpoints of geodesics that do not have coordinates over the
        base ring, we cannot produce defining half spaces since these would
        require equations over a quadratic extension as well::

            sage: H.half_circle(0, 2).start().half_spaces()
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 ...

        ::

            sage: H = HyperbolicPlane(RR)

            sage: H(I).half_spaces()
            {{(x^2 + y^2) + 2.00000000000000*x - 1.00000000000000 ≤ 0},
             {x ≥ 0},
             {(x^2 + y^2) - 1.00000000000000 ≥ 0}}

            sage: H(0).half_spaces()
            {{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}}

        """
        x0, y0 = self.coordinates(model="klein")

        if self.is_finite():
            return HyperbolicHalfSpaces(
                # x ≥ x0
                [
                    self.parent().half_space(-x0, 1, 0, model="klein"),
                    # y ≥ y0
                    self.parent().half_space(-y0, 0, 1, model="klein"),
                    # x + y ≤ x0 + y0
                    self.parent().half_space(x0 + y0, -1, -1, model="klein"),
                ]
            )
        else:
            return HyperbolicHalfSpaces(
                # left of the line from (0, 0) to this point
                [
                    self.parent().half_space(0, -y0, x0, model="klein"),
                    # right of a line to this point with a starting point right of (0, 0)
                    self.parent().half_space(
                        -x0 * x0 - y0 * y0, y0 + x0, y0 - x0, model="klein"
                    ),
                ]
            )

    def coordinates(self, model="half_plane", ring=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return projective coordinates in the
        Poincaré half plane model.

        If ``model`` is ``"klein"``, return Euclidean coordinates in the Klein model.

        If no ``ring`` has been specified, an appropriate extension of the base
        ring of the :class:`HyperbolicPlane` is chosen where these coordinates
        live.
        """
        # TODO: Implement ring.
        # TODO: Fix documentation.

        if model == "half_plane":
            coordinates = self.coordinates(model="klein", ring=ring)
            if coordinates is None:
                return None

            x, y = coordinates

            if self == self.parent().infinity() or self.is_ultra_ideal():
                raise ValueError("point has no coordinates in the upper half plane")

            denominator = 1 - y

            if not self.is_finite():
                return (x / denominator, self.parent().base_ring().zero())

            square = 1 - x * x - y * y
            try:
                sqrt = square.sqrt()
                if sqrt not in self.parent().base_ring():
                    raise ValueError(f"square root of {square} not in {self.parent().base_ring()}")
            except ValueError:
                if ring == "try":
                    return None
                raise

            return (x / denominator, sqrt / denominator)

        raise NotImplementedError

    def real(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the real part of this point in the upper half plane model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: p = H(I + 2)
            sage: p.real()
            2

        """
        return self.coordinates(model="half_plane")[0]

    def imag(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the imaginary part of this point in the upper half plane model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: p = H(I + 2)
            sage: p.imag()
            1

        """
        return self.coordinates(model="half_plane")[1]

    def _apply_isometry_klein(self, isometry):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

    def segment(self, end, check=True, assume_normalized=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the oriented segment from this point to ``end``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: H(0).segment(I)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≤ 0}

        """
        return self.parent().segment(
            self.parent().geodesic(self, end),
            start=self,
            end=end,
            check=check,
            assume_normalized=assume_normalized,
        )

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ.zero()

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return HyperbolicVertices([self])


class HyperbolicPointFromCoordinates(HyperbolicPoint):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    def __init__(self, parent, x, y):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

        if x.parent() is not parent.base_ring():
            raise TypeError("x must be an element of the base ring")
        if y.parent() is not parent.base_ring():
            raise TypeError("y must be an element of the base ring")

        self._coordinates = (x, y)

    def coordinates(self, model="half_plane", ring=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return projective coordinates in the
        Poincaré half plane model.

        If ``model`` is ``"klein"``, return Euclidean coordinates in the Klein model.

        If no ``ring`` has been specified, an appropriate extension of the base
        ring of the :class:`HyperbolicPlane` is chosen where these coordinates
        live.
        """
        # TODO: Implement ring.
        # TODO: Fix documentation.
        if model == "klein":
            return self._coordinates

        return super().coordinates(model=model, ring=ring)

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            if isinstance(other, HyperbolicPointFromGeodesic):
                return other == self

            # TODO: Turn this into a proper predicate.
            equal = self.parent().geometry.equal
            return equal(
                self.coordinates(model="klein"), other.coordinates(model="klein")
            )

        super()._richcmp_(other, op)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self == self.parent().infinity():
            return "∞"

        if self.is_ultra_ideal():
            return repr(self.coordinates(model="klein"))

        x, y = self.coordinates(model="half_plane")

        from sage.all import PowerSeriesRing

        # We represent x + y*I in R[[I]] so we do not have to reimplement printing ourselves.
        return repr(PowerSeriesRing(self.parent().base_ring(), names="I")([x, y]))

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        def point(parent):
            return parent.point(*self._coordinates, model="klein", check=False)

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of a point")

        if ring is not None or geometry is not None:
            self = self.parent().change_ring(ring, geometry=geometry).point(*self._coordinates, model="klein", check=False)

        return self


class HyperbolicPointFromGeodesic(HyperbolicPoint):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    def __init__(self, parent, geodesic):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("x must be an oriented geodesic")

        self._geodesic = geodesic

    def is_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return True

    def is_ultra_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return False

    def is_finite(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return False

    @cached_method
    def coordinates(self, model="half_plane", ring=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Did we implement ring?
        if model == "klein":
            a, b, c = self._geodesic.equation(model="half_plane")

            # TODO: Turn this into a proper predicate.
            zero = self.parent().geometry.zero
            sgn = self.parent().geometry.sgn

            if zero(a):
                point = None
                if sgn(b) > 0:
                    point = self.parent().point(0, 1, model="klein", check=False)
                else:
                    point = self.parent().point(-c / b, 0, model="half_plane", check=False)

                return point.coordinates(model=model, ring=ring)
            else:
                discriminant = b * b - 4 * a * c
                try:
                    sqrt = discriminant.sqrt()
                    if sqrt not in self.parent().base_ring():
                        raise ValueError(f"square root of {discriminant} not in {self.parent().base_ring()}")
                except ValueError:
                    if ring == "try":
                        return None
                    raise

                endpoints = ((-b - sqrt) / (2 * a), (-b + sqrt) / (2 * a))

                return self.parent().point(
                        (min if sgn(a) > 0 else max)(endpoints),
                        0,
                        model="half_plane",
                        check=False,
                    ).coordinates(model=model, ring=ring)

        return super().coordinates(model=model, ring=ring)

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            if isinstance(other, HyperbolicPointFromGeodesic):
                if self._geodesic == other._geodesic:
                    return True
                if self._geodesic == -other._geodesic:
                    return False

                intersection = self._geodesic._intersection(other._geodesic)

                return self == intersection and other == intersection

            if not other.is_ideal():
                return False

            return other in self._geodesic and self._geodesic.parametrize(other, model="euclidean") < 0

        super()._richcmp_(other, op)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self == self.parent().infinity():
            return "∞"

        if not self.is_ultra_ideal():
            coordinates = self.coordinates(model="half_plane", ring="try")

            if coordinates is not None:
                return repr(self.parent().point(*coordinates, model="half_plane", check=False))

        coordinates = self.coordinates(model="klein", ring="try")

        if coordinates is not None:
            return repr(coordinates)

        from sage.all import RR
        return repr(self.change_ring(RR))

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of a point")

        if ring is not None or geometry is not None:
            self = self._geodesic.change(ring=ring, geometry=geometry).start()

        return self


class HyperbolicConvexPolygon(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`half spaces <HyperbolicHalfSpace>`.
    """

    def __init__(self, parent, half_spaces):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

        if not isinstance(half_spaces, HyperbolicHalfSpaces):
            raise TypeError("half_spaces must be HyperbolicHalfSpaces")

        self._half_spaces = half_spaces

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO
        pass

    # TODO: Add examples.
    def _normalize(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

        if not self._half_spaces:
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            sage: polygon(*(H.half_circle(0, 1).half_spaces() + H.half_circle(0, 2).half_spaces()))._normalize_drop_trivially_redundant()
            {(x^2 + y^2) - 1 ≤ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        """
        reduced = []

        for half_space in self._half_spaces:
            if reduced:
                a, b, c = half_space.equation(model="klein")
                A, B, C = reduced[-1].equation(model="klein")

                # TODO: Turn this into a proper predicate.
                equal = self.parent().geometry.equal
                sgn = self.parent().geometry.sgn
                if equal(c * B, C * b) and sgn(b) == sgn(B) and sgn(c) == sgn(C):
                    # The half spaces are parallel in the Euclidean plane. Since we
                    # assume spaces to be sorted by inclusion, we can drop this
                    # space.
                    continue

            reduced.append(half_space)

        return self.parent().polygon(
            reduced, check=False, assume_sorted=True, assume_minimal=True
        )

    def _normalize_drop_euclidean_redundant(self, boundary):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            sage: polygon(*H.infinity().half_spaces())._normalize_drop_euclidean_redundant(
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

        half_spaces = list(self._half_spaces)

        half_spaces = (
            half_spaces[half_spaces.index(boundary):]
            + half_spaces[: half_spaces.index(boundary)]
        )
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
                    raise NotImplementedError(
                        f"B and C are in unsupported configuration: {BC}"
                    )

            elif AB == "negative":
                required = True

            elif AB == "anti-parallel":
                required = True

            elif AB == "concave":
                required = True

            else:
                raise NotImplementedError(
                    f"A and B are in unsupported configuration: {AB}"
                )

            if required:
                required_half_spaces.append(B)
            elif len(required_half_spaces) > 1:
                half_spaces.append(required_half_spaces.pop())

        required_half_spaces = HyperbolicHalfSpaces(
            required_half_spaces, assume_sorted="rotated"
        )

        return self.parent().polygon(
            required_half_spaces, check=False, assume_sorted=True, assume_minimal=True
        )

    def _normalize_drop_unit_disk_redundant(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            sage: polygon(*H.infinity().half_spaces())._normalize_drop_unit_disk_redundant()
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

        for A, B, C in self._half_spaces.triples():
            AB = (
                None
                if A.boundary()._configuration(B.boundary()) == "concave"
                else A.boundary()._intersection(B.boundary())
            )
            BC = (
                None
                if B.boundary()._configuration(C.boundary()) == "concave"
                else B.boundary()._intersection(C.boundary())
            )

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
            raise NotImplementedError(
                "there is no convex set to represent the full space yet"
            )

        if len(required_half_spaces) == 1:
            return required_half_spaces[0]

        return self.parent().polygon(
            required_half_spaces, check=False, assume_sorted=True, assume_minimal=True
        )

    def _euclidean_boundary(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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

            sage: polygon(*H(I).half_spaces())._euclidean_boundary()
            I

        An intersection which is a single point on the boundary of the unit
        disk::

            sage: polygon(*H.infinity().half_spaces())._euclidean_boundary()
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
        if len(self._half_spaces) == 0:
            raise ValueError("list of half spaces must not be empty")

        if len(self._half_spaces) == 1:
            return next(iter(self._half_spaces))

        # Randomly shuffle the half spaces so the loop below runs in expected linear time.
        from sage.all import shuffle

        random_half_spaces = list(self._half_spaces)
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
                # TODO: This does not do any epsilons yet, i.e., it should probably use predicates somehow.
                interval = RealSet(-oo, oo)

                for constraining in random_half_spaces:
                    if constraining is half_space:
                        break

                    intersection = boundary._intersection(constraining.boundary())

                    if intersection is None:
                        # constraining is anti-parallel to half_space
                        if (
                            boundary.parametrize(0, model="euclidean", check=False)
                            not in constraining
                        ):
                            return self.parent().empty_set()

                        # The intersection is non-empty, so this adds no further constraints.
                        continue

                    λ = boundary.parametrize(
                        intersection, model="euclidean", check=False
                    )

                    # Determine whether this half space constrains to (-∞, λ] or [λ, ∞).
                    if (
                        boundary.parametrize(λ + 1, model="euclidean", check=False)
                        in constraining
                    ):
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        half_spaces = [
            half_space
            for half_space in self._half_spaces
            if point in half_space.boundary()
        ]

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

            raise NotImplementedError(
                f"cannot extend point to segment when half spaces are in configuration {configuration}"
            )

        if point.is_ultra_ideal():
            # There is no actual intersection in the hyperbolic plane.
            return self.parent().empty_set()

        return point

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ(2)

    def equations(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the equations describing the boundary of this polygon.

        The output is minimal and sorted by slope in the Klein model.
        """
        raise NotImplementedError

    def edges(self, as_segments=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Define in HyperbolicConvexSet
        r"""
        Return the :class:`segments <HyperbolicOrientedSegment>` and
        :class:`geodesics <HyperbolicOrientedGeodesic>` defining this polygon.
        """
        edges = []

        boundaries = [half_space.boundary() for half_space in self._half_spaces]

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
            elif AB == "negative":
                start = None
            else:
                raise NotImplementedError(
                    f"cannot determine edges when boundaries are in configuration {AB}"
                )

            if BC == "convex":
                end = B._intersection(C)
                if not end.is_finite():
                    end = None
            elif BC == "concave":
                pass
            elif BC == "negative":
                end = None
            else:
                raise NotImplementedError(
                    f"cannot determine edges when boundaries are in configuration {BC}"
                )

            edges.append(
                self.parent().segment(
                    B, start=start, end=end, assume_normalized=as_segments, check=False
                )
            )

        return edges

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Define in HyperbolicConvexSet
        r"""
        Return the vertices of this polygon, i.e., the (possibly ideal) end
        points of the :meth:`edges`, in counterclockwise order.
        """
        vertices = []

        vertex = None
        for edge in self.edges():
            start = edge.start()
            if vertex is not None and start != vertex:
                vertices.append(start)

            end = edge.end()
            vertices.append(end)
            vertex = end

        return HyperbolicVertices(vertices)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self._half_spaces

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return " ∩ ".join([repr(half_space) for half_space in self._half_spaces])

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        kwds.setdefault("color", "#efffff")
        kwds.setdefault("edgecolor", "#d1d1d1")

        if len(self._half_spaces) == 0:
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

        return self._enhance_plot(
            hyperbolic_path(commands, model=model, **kwds), model=model
        )

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if ring is not None or geometry is not None:
            self = self.parent().change_ring(ring, geometry=geometry).polygon(
                [half_space.change(ring=ring, geometry=geometry) for half_space in self._half_spaces],
                check=False,
                assume_sorted=True,
                assume_minimal=True,
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of a polygon")

        return self

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Pass to normalization
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicConvexPolygon):
                return False
            return self._half_spaces == other._half_spaces


class HyperbolicSegment(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A segment (possibly infinite) in the hyperbolic plane.

    This is an abstract base class of :class:`HyperbolicOrientedSegmented` and
    :class:`HyperbolicUnorientedSegment`.

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
        sage: H = HyperbolicPlane()

        sage: segment = H.segment(H.vertical(0), start=I)

        sage: isinstance(segment, HyperbolicSegment)
        True

        sage: isinstance(segment.unoriented(), HyperbolicSegment)
        True

    """

    def __init__(self, parent, geodesic, start=None, end=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("geodesic must be an oriented hyperbolic geodesic")

        if start is not None and not isinstance(start, HyperbolicPoint):
            raise TypeError("start must be a hyperbolic point")

        if end is not None and not isinstance(end, HyperbolicPoint):
            raise TypeError("end must be a hyperbolic point")

        self._geodesic = geodesic
        self._start = start
        self._end = end

    def _endpoint_half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        a, b, c = self._geodesic.equation(model="klein")

        if self._start is not None:
            x, y = self._start.coordinates(model="klein")
            yield self.parent().half_space(
                -c * x + b * y, c, -b, model="klein", check=False
            )

        if self._end is not None:
            x, y = self._end.coordinates(model="klein")
            yield self.parent().half_space(
                c * x - b * y, -c, b, model="klein", check=False
            )

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        bounds = [repr(self._geodesic)]
        bounds.extend(repr(half_space) for half_space in self._endpoint_half_spaces())

        return " ∩ ".join(bounds)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self._geodesic.half_spaces() + HyperbolicHalfSpaces(
            self._endpoint_half_spaces()
        )

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import RR

        kwds["fill"] = False

        self = self.change_ring(RR)
        plot = hyperbolic_path(
            [
                BezierPath.Command("MOVETO", [self.start()]),
                BezierPath.Command("LINETO", [self.end()]),
            ],
            model=model,
            **kwds,
        )

        return self._enhance_plot(plot, model=model)

    def start(self, finite=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self._start is not None:
            return self._start

        return self._geodesic.start()

    def end(self, finite=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if self._end is not None:
            return self._end

        return self._geodesic.end()

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Compares this segment to ``other`` with respect to ``op``.

        EXAMPLES:

        Oriented segments are equal if they have the same start and end points::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: H(I).segment(2*I) == H(2*I).segment(I)
            False

        For an unoriented segment the endpoints must be the same but order does not matter::

            sage: H(I).segment(2*I).unoriented() == H(2*I).segment(I).unoriented()
            True

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if type(self) is not type(other):
                return False
            return (
                self.geodesic() == other.geodesic()
                and self.vertices() == other.vertices()
            )

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if ring is not None or geometry is not None:
            start = self._start.change(ring=ring, geometry=geometry) if self._start is not None else None
            end = self._end.change(ring=ring, geometry=geometry) if self._end is not None else None

            self = self.parent().change_ring(ring=ring, geometry=geometry).segment(
                self._geodesic.change(ring=ring, geometry=geometry),
                start=start,
                end=end,
                check=False,
                assume_normalized=True,
                oriented=self.is_oriented(),
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            if not self.is_oriented():
                raise NotImplementedError("cannot orient unoriented segment")

            self = self.parent().segment(
                self._geodesic,
                start=self._start,
                end=self._end,
                check=False,
                assume_normalized=True,
                oriented=oriented,
            )

        return self

    def geodesic(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        geodesic = self._geodesic
        if not self.is_oriented():
            geodesic = geodesic.unoriented()
        return geodesic

    def vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return HyperbolicVertices([self.start(), self.end()])


class HyperbolicUnorientedSegment(HyperbolicSegment):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    An unoriented (possibly infinity) segment in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: segment = H.segment(H.vertical(0), start=I).unoriented()

    """


class HyperbolicOrientedSegment(HyperbolicSegment, HyperbolicOrientedConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.
    """

    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return self.parent().segment(
            -self._geodesic, self._end, self._start, check=False, assume_normalized=True
        )

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        start = self._start
        end = self._end

        if start is not None:
            if start not in self._geodesic:
                raise ValueError("start point must be on the geodesic")

        if end is not None:
            if end not in self._geodesic:
                raise ValueError("end point must be on the geodesic")

        # TODO: Check end >= start and end > start if require_normalized.

        # TODO: Check that end > start even if unoriented since otherwise the printing is wrong.

        # TODO: Check start & end finite if require_normalized.

    def _normalize(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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
            λ_start = self._geodesic.parametrize(start, model="euclidean", check=False)

        if end is not None:
            λ_end = self._geodesic.parametrize(end, model="euclidean", check=False)

        if start is not None:
            if not start.is_finite():
                # TODO: Turn this into a proper predicate.
                sgn = self.parent().geometry.sgn
                if sgn(λ_start) > 0:
                    return (
                        self.parent().empty_set() if start.is_ultra_ideal() else start
                    )
                start = None

        if end is not None:
            if not end.is_finite():
                # TODO: Turn this into a proper predicate.
                sgn = self.parent().geometry.sgn
                if sgn(λ_end) < 0:
                    return self.parent().empty_set() if end.is_ultra_ideal() else end
                end = None

        if start is None and end is None:
            return self._geodesic

        assert (start is None or not start.is_ultra_ideal()) and (
            end is None or not end.is_ultra_ideal()
        )

        if start == end:
            return start

        return self.parent().segment(
            self._geodesic, start=start, end=end, check=False, assume_normalized=True
        )

    def _is_valid(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
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

    def configuration(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        if self._geodesic == other._geodesic:
            if self == other:
                return "equal"
            raise NotImplementedError(
                "cannot determine configuration of segments on the same geodesic"
            )

        if self._geodesic == -other._geodesic:
            if self == -other:
                return "negative"
            raise NotImplementedError(
                "cannot determine configuration of segments on the same geodesic"
            )

        intersection = self.intersection(other)

        if intersection is None:
            raise NotImplementedError(
                "cannot determine configuration of segments that do not intersect"
            )

        if intersection.is_finite():
            if intersection == self.end(finite=True) and intersection == other.start(
                finite=True
            ):
                return "join"

            raise NotImplementedError(
                "cannot determine configuration of segments that intersect in a finite point"
            )

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

        raise NotImplementedError(
            "cannot determine configuration of segments that intersect in an infinite point"
        )


class HyperbolicEmptySet(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    The empty subset of the hyperbolic plane.
    """

    def __init__(self, parent):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        super().__init__(parent)

    def _richcmp_(self, other, op):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return "{}"

    def _apply_isometry_klein(self, isometry):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: S = HyperbolicPlane(QQ).empty_set()
            sage: S.apply_isometry(matrix(2, [2, 1, 1, 1])) is S
            True
        """
        return self

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import Graphics

        return self._enhance_plot(Graphics(), model=model)

    def change_ring(self, ring):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return HyperbolicPlane(ring).empty_set()

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ(-1)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return HyperbolicHalfSpaces(
            [
                self.parent().half_circle(-2, 1).right_half_space(),
                self.parent().half_circle(2, 1).right_half_space(),
            ]
        )

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if ring is not None:
            self = self.parent().change_ring(ring, geometry=geometry).empty_set()

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of empty set")

        return self


def sl2_to_so12(m):
    # TODO: Check documentation.
    # TODO: Check INPUT
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    Return the lift of the 2x2 matrix ``m`` inside ``SO(1,2)``.
    """
    from sage.matrix.constructor import matrix

    if m.nrows() != 2 or m.ncols() != 2:
        raise ValueError("invalid matrix")
    a, b, c, d = m.list()
    return matrix(
        3,
        [
            a * d + b * c,
            a * c - b * d,
            a * c + b * d,
            a * b - c * d,
            (a**2 - b**2 - c**2 + d**2) / 2,
            (a**2 + b**2 - c**2 - d**2) / 2,
            a * b + c * d,
            (a**2 - b**2 + c**2 - d**2) / 2,
            (a**2 + b**2 + c**2 + d**2) / 2,
        ],
    )


class SortedSet:
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A set of objects sorted by :meth:`SortedSet.cmp`.

    This is used to efficiently represent
    :meth:`HyperbolicConvexSet.half_spaces`,
    :meth:`HyperbolicConvexSet.vertices`, and
    :meth:`HyperbolicConvexSet.edges`. In particular, it allows us to create
    and merge such sets in linear time.
    """

    def __init__(self, entries, assume_sorted=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if assume_sorted is None:
            assume_sorted = isinstance(entries, SortedSet)

        if assume_sorted == "rotated":
            min = 0
            for i, entry in enumerate(entries):
                if self._lt_(entry, entries[min]):
                    min = i

            entries = entries[min:] + entries[:min]
            assume_sorted = True

        if not assume_sorted:
            entries = self._merge(*[[entry] for entry in entries])

        self._entries = tuple(entries)

    def _lt_(self, lhs, rhs):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        raise NotImplementedError

    def _merge(self, *sets):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the merge of sorted lists of ``sets``.

        Naturally, when there are a lot of small sets, such a merge sort takes
        quasi-linear time. However, when there are only a few sets, this runs
        in linear time.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: HyperbolicHalfSpaces([])._merge()
            []

            sage: HyperbolicHalfSpaces([])._merge(*[[half_space] for half_space in H.real(0).half_spaces()])
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpaces([])._merge(list(H.real(0).half_spaces()), list(H.real(0).half_spaces()))
            [{(x^2 + y^2) + x ≤ 0}, {(x^2 + y^2) + x ≤ 0}, {x ≥ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpaces([])._merge(*[[half_space] for half_space in list(H.real(0).half_spaces()) * 2])
            [{(x^2 + y^2) + x ≤ 0}, {(x^2 + y^2) + x ≤ 0}, {x ≥ 0}, {x ≥ 0}]

        """
        # A standard merge-sort implementation.
        count = len(sets)

        if count == 0:
            return []

        if count == 1:
            return sets[0]

        # The non-trivial base case.
        if count == 2:
            A = sets[0]
            B = sets[1]
            merged = []

            while A and B:
                if self._lt_(A[-1], B[-1]):
                    merged.append(B.pop())
                else:
                    merged.append(A.pop())

            merged.reverse()

            return A + B + merged

        # Divide & Conquer recursively.
        return self._merge(
            *(self._merge(*sets[: count // 2]), self._merge(*sets[count // 2:]))
        )

    def __eq__(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return whether this set is equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.vertical(0).vertices() == (-H.vertical(0)).vertices()
            True

        """
        if type(other) != type(self):
            return False

        return self._entries == other._entries

    def __ne__(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return whether this set is not equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.vertical(0).vertices() != H.vertical(1).vertices()
            True

        """
        return not (self == other)

    def __add__(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        entries = self._merge(list(self._entries), list(other._entries))
        return type(self)(entries, assume_sorted=True)

    def __repr__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return a printable representation of this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.half_circle(0, 1).vertices()
            {-1, 1}

        """
        return "{" + repr(self._entries)[1:-1] + "}"

    def __iter__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return iter(self._entries)

    def __len__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        return len(self._entries)

    def triples(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        for i in range(len(self._entries)):
            yield self._entries[i - 1], self._entries[i], self._entries[
                (i + 1) % len(self._entries)
            ]


class HyperbolicVertices(SortedSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    r"""
    A set of vertices on the boundary of a convex set in the hyperbolic plane,
    sorted in counterclockwise order.

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.vertices`

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)
        sage: V = H.vertical(0).vertices()
        sage: V
        {0, ∞}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicVertices
        sage: isinstance(V, HyperbolicVertices)
        True

    """

    def __init__(self, vertices, assume_sorted=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if len(vertices) == 0:
            raise ValueError("vertex set must not be empty")

        min = vertices[0]
        self._min = min.parent().point(0, 0, model="klein")

        for vertex in vertices[1:]:
            if self._lt_(vertex, min):
                min = vertex
            else:
                assert self._lt_(min, vertex), "vertices must not contain duplicates"

        self._min = min

        super().__init__(vertices, assume_sorted=assume_sorted)

    def _lt_(self, lhs, rhs):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: This is a bit hacky. And probably slow. (Do we need predicates here?)
        from sage.all import RR

        if lhs == rhs:
            return False

        if lhs.change_ring(RR).coordinates(model="klein") < rhs.change_ring(
            RR
        ).coordinates(model="klein"):
            return True
        elif lhs.change_ring(RR).coordinates(model="klein") == rhs.change_ring(
            RR
        ).coordinates(model="klein"):
            if lhs.coordinates(model="klein") < rhs.coordinates(model="klein"):
                return True

        return False


class HyperbolicHalfSpaces(SortedSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    @classmethod
    def _lt_(cls, lhs, rhs):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: This is essentially atan2.
        r"""
        Return whether the half space ``lhs`` is smaller than ``rhs`` in a cyclic
        ordering of normal vectors, i.e., in an ordering that half spaces
        whether their normal points to the left/right, the slope of the
        geodesic, and finally by containment.

        This ordering is such that :meth:`HyperbolicPlane.intersection` can be
        computed in linear time for two hyperbolic convex sets.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicHalfSpaces
            sage: H = HyperbolicPlane(QQ)

        A half space is equal to itself::

            sage: HyperbolicHalfSpaces._lt_(H.vertical(0).left_half_space(), H.vertical(0).left_half_space())
            False

        A half space whose normal in the Klein model points to the left is
        smaller than one whose normal points to the right::

            sage: HyperbolicHalfSpaces._lt_(H.vertical(1).left_half_space(), H.half_circle(0, 1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.vertical(0).left_half_space(), -H.vertical(0).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(-H.half_circle(-1, 1).left_half_space(), -H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(-H.half_circle(-1, 1).left_half_space(), -H.vertical(1/2).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.vertical(1).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.vertical(1/2).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True

        Half spaces are ordered by the slope of their normal in the Klein model::

            sage: HyperbolicHalfSpaces._lt_(H.vertical(-1).left_half_space(), H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(-H.half_circle(-1, 1).left_half_space(), H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.half_circle(-1, 1).left_half_space(), -H.vertical(1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.vertical(0).left_half_space(), H.vertical(1).left_half_space())
            True

        Parallel half spaces in the Klein model are ordered by inclusion::

            sage: HyperbolicHalfSpaces._lt_(-H.half_circle(-1, 1).left_half_space(), H.vertical(1/2).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(-H.vertical(1/2).left_half_space(), H.half_circle(-1, 1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.half_circle(0, 2).left_half_space(), H.half_circle(0, 1).left_half_space())
            True
            sage: HyperbolicHalfSpaces._lt_(H.half_circle(0, 1).right_half_space(), H.half_circle(0, 2).right_half_space())
            True

        Verify that comparisons are projective::

            sage: HyperbolicHalfSpaces._lt_(H.geodesic(5, -5, -1, model="half_plane").left_half_space(), H.geodesic(5/13, -5/13, -1/13, model="half_plane").left_half_space())
            False
            sage: HyperbolicHalfSpaces._lt_(H.geodesic(5/13, -5/13, -1/13, model="half_plane").left_half_space(), H.geodesic(5, -5, -1, model="half_plane").left_half_space())
            False

        """
        a, b, c = lhs.equation(model="klein")
        aa, bb, cc = rhs.equation(model="klein")

        # TODO: Should we use predicates here? If so, which?

        def normal_points_left(b, c):
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
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


class BezierPath(GraphicPrimitive):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
    def __init__(self, commands, options=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        options = options or {}

        valid_options = self._allowed_options()
        for option in options:
            if option not in valid_options:
                raise RuntimeError(f"option {option} not valid")

        self._commands = commands
        super().__init__(options)

    def _allowed_options(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Return the options that are supported by a path.

        We support all the options that are understood by a SageMath polygon.

        """
        from sage.plot.polygon import Polygon

        return Polygon([], [], {})._allowed_options()

    def _render_on_subplot(self, subplot):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        r"""
        Render this path on the subplot.

        Matplotlib was not really made to draw things that extend to infinity.
        The trick here is to register a callback that redraws whenever the
        viewbox of the plot changes, e.g., as more objects are added to the
        plot.
        """
        # Rewrite options to only contain matplotlib compatible entries
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

        fill_path = Path([(0, 0)])
        edge_path = Path([(0, 0)])

        from matplotlib.patches import PathPatch

        fill_patch = PathPatch(fill_path, **matplotlib_options)
        edge_patch = PathPatch(edge_path, **matplotlib_options)

        options = self.options()
        fill = options.pop("fill")
        if fill:
            subplot.axes.add_patch(fill_patch)
        subplot.axes.add_patch(edge_patch)

        # Translate SageMath options to matplotlib style.
        fill_patch.set_linewidth(float(options["thickness"]))
        if "linestyle" in options:
            fill_patch.set_linestyle(options["linestyle"])
        fill_patch.set_alpha(float(options["alpha"]))

        from sage.plot.colors import to_mpl_color

        color = to_mpl_color(options.pop("rgbcolor"))

        fill_patch.set_fill(True)
        edge_patch.set_fill(False)

        edge_color = options.pop("edgecolor")
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

        fill_patch.set_label(options["legend_label"])

        def redraw(_=None):
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
            r"""
            Redraw after the viewport has been rescaled to make sure that
            infinite rays reach the end of the viewport.
            """
            self._redraw_on_subplot(subplot, fill_patch, fill=True)
            self._redraw_on_subplot(subplot, edge_patch, fill=False)

        subplot.axes.callbacks.connect("ylim_changed", redraw)
        subplot.axes.callbacks.connect("xlim_changed", redraw)
        redraw()

    def _redraw_on_subplot(self, subplot, patch, fill):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
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
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
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

        def extend(path):
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
            vertices.extend(path.vertices[1:])
            codes.extend(path.codes[1:])

        if command.code == "MOVETO":
            (pos,) = command.args
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
                (target,) = command.args

                if direction is not None and pos != target:
                    raise ValueError(
                        f"Cannot execute LINETO from infinite point at {pos} + λ {direction} when going to {target}"
                    )

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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        try:
            from matplotlib.transforms import Bbox

            bbox = Bbox.null()

            pos = None
            for command in self._commands:
                if command.code in ["MOVETO", "LINETO"]:
                    (pos,) = command.args
                    bbox.update_from_data_xy([pos], ignore=False)
                elif command.code == "ARCTO":
                    target, center = command.args
                    bbox = bbox.union(
                        [bbox, self._arc_path(center, pos, target).get_extents()]
                    )
                    pos = target
                elif command.code == "RARCTO":
                    target, center = command.args
                    bbox = bbox.union(
                        [
                            bbox,
                            self._arc_path(
                                center, target, pos, reverse=True
                            ).get_extents(),
                        ]
                    )
                    pos = target
                elif command.code in ["LINETOINFINITY", "MOVETOINFINITY"]:
                    # TODO: Add a bit to the bounding box so these are always visible.
                    pass
                else:
                    raise NotImplementedError(
                        f"cannot determine bounding box for {command.code} command"
                    )

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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from matplotlib.path import Path
        from math import atan2, pi
        from sage.all import vector

        # TODO: How many segments do we need?
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

    @classmethod
    def hyperbolic_path(cls, commands, model, **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if len(commands) < 2:
            raise ValueError("a path must contain at least two points")

        commands.reverse()

        command = commands.pop()
        if command.code == "MOVETO":
            (pos,) = command.args
            if model == "half_plane" and pos == pos.parent().infinity():
                next = commands[-1].args[0]
                bezier_commands = [
                    BezierPath.Command(
                        "MOVETOINFINITY", [next.coordinates(model=model), (0, 1)]
                    )
                ]
            else:
                from sage.all import RR

                bezier_commands = [
                    BezierPath.Command(
                        "MOVETO", [pos.change_ring(RR).coordinates(model=model)]
                    )
                ]
        else:
            raise ValueError("path must start with MOVETO or MOVETOINFINITY command")

        while commands:
            command = commands.pop()

            if command.code == "LINETO":
                (target,) = command.args
                bezier_commands.extend(
                    cls._hyperbolic_segment(pos, target, model=model)
                )
                pos = target
            elif command.code == "MOVETO":
                (target,) = command.args
                bezier_commands.extend(cls._hyperbolic_move(pos, target, model=model))
                pos = target
            else:
                raise NotImplementedError("unsuported hyperbolic plotting code")

        if bezier_commands[-1].code == "LINETOINFINITY" and bezier_commands[-1].args[
            -1
        ] == (1, 0):
            assert bezier_commands[0].code == "MOVETOINFINITY" and bezier_commands[
                1
            ].args[-1] == (0, 1)
            bezier_commands.append(bezier_commands[0])

        assert len(bezier_commands) >= 2

        return BezierPath(bezier_commands, kwds)

    @classmethod
    def _hyperbolic_segment(cls, start, end, model):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if start == end:
            raise ValueError(f"cannot draw segment from point {start} to itself ({end})")

        if model == "half_plane":
            if start == start.parent().infinity():
                return [
                    BezierPath.Command("MOVETOINFINITY", [end.coordinates(), (0, 1)]),
                    BezierPath.Command("LINETO", [end.coordinates()]),
                ]

            from sage.all import RR

            if end == end.parent().infinity():
                return [
                    BezierPath.Command(
                        "LINETOINFINITY",
                        [start.change_ring(RR).coordinates(model="half_plane"), (0, 1)],
                    )
                ]

            p = start.change_ring(RR).coordinates(model="half_plane")
            q = end.change_ring(RR).coordinates(model="half_plane")

            if (p[0] - q[0]).abs() < (p[1] - q[1]).abs() * 1e-6:
                # This segment is (almost) vertical. We plot it as if it were
                # vertical to avoid numeric issus.
                return [BezierPath.Command("LINETO", [q])]

            geodesic = start.change_ring(RR).parent().geodesic(start, end, check=False)
            center = (
                (geodesic.start().coordinates()[0] + geodesic.end().coordinates()[0])
                / 2,
                0,
            )

            return [
                BezierPath.Command("RARCTO" if p[0] < q[0] else "ARCTO", [q, center])
            ]
        elif model == "klein":
            from sage.all import RR

            return [
                BezierPath.Command(
                    "LINETO", [end.change_ring(RR).coordinates(model="klein")]
                )
            ]
        else:
            raise NotImplementedError("cannot draw segment in this model")

    @classmethod
    def _hyperbolic_move(cls, start, end, model):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        if start == end:
            raise ValueError("cannot move from point to itself")

        if start.is_finite():
            raise ValueError("starting point of move must be ideal")

        if end.is_finite():
            raise ValueError("end of move must be ideal")

        if model == "half_plane":
            if start == start.parent().infinity():
                from sage.all import RR

                return [
                    BezierPath.Command(
                        "MOVETOINFINITY", [end.change_ring(RR).coordinates(), (0, 1)]
                    ),
                    BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()]),
                ]

            if end == end.parent().infinity():
                from sage.all import RR

                return [
                    BezierPath.Command(
                        "LINETOINFINITY", [start.change_ring(RR).coordinates(), (1, 0)]
                    )
                ]

            from sage.all import RR

            if (
                start.change_ring(RR).coordinates()[0]
                < end.change_ring(RR).coordinates()[0]
            ):
                return [
                    BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()])
                ]
            else:
                return [
                    BezierPath.Command(
                        "LINETOINFINITY", [start.change_ring(RR).coordinates(), (1, 0)]
                    ),
                    BezierPath.Command(
                        "MOVETOINFINITY", [end.change_ring(RR).coordinates(), (-1, 0)]
                    ),
                    BezierPath.Command("LINETO", [end.change_ring(RR).coordinates()]),
                ]

            raise NotImplementedError(
                f"cannot move from {start} to {end} in half plane model"
            )
        elif model == "klein":
            # TODO: The actual arc should not be visible.
            from sage.all import RR

            return [
                BezierPath.Command(
                    "ARCTO", [end.change_ring(RR).coordinates(model="klein"), (0, 0)]
                )
            ]
        else:
            raise NotImplementedError("cannot move in this model")


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
    # TODO: Check documentation.
    # TODO: Check INPUT
    # TODO: Check SEEALSO
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

    if options["legend_label"]:
        g.legend(True)
        g._legend_colors = [options["legend_color"]]
    return g
