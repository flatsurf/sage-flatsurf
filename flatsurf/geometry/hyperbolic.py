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

    sage: H_algebraic = HyperbolicPlane(AA)
    sage: H_algebraic(sqrt(2) + I)
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

    We usually display objects as if they were defined in the upper half plane
    model. However, internally, we store most objects in a representation in
    the Klein model. In that model it tends to be easier to perform
    computations without having to extend the base ring and we can also rely on
    standard algorithms for geometry in the Euclidean plane.

    For the Klein model, we use a unit disk centered at (0, 0). The map from
    the upper half plane sends the imaginary unit `i` to the center at the
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

    A geodesic in the upper half plane is given by an equation of the form

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

    in the upper half plane model.

    Note that the intersection of two geodesics defined by coefficients in a
    field `K` in the Klein model has coordinates in `K` in the Klein model.
    The corresponding statement is not true for the upper half plane model.

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

    def _coerce_map_from_(self, other):
        r"""
        Return a coercion map from ``other`` to this hyperbolic plane.

        EXAMPLES:

        Coercions between base rings induce coercion between hyperbolic planes,
        due to the :meth:`construction` functor::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: HyperbolicPlane(AA).has_coerce_map_from(H)
            True

        Base ring elements coerce as points on the real line::

            sage: H.has_coerce_map_from(QQ)
            True
            sage: H.has_coerce_map_from(ZZ)
            True

        Complex numbers do not coerce into the hyperbolic plane since that
        coercion would not be total::

            sage: H.has_coerce_map_from(CC)
            False
            sage: H.has_coerce_map_from(I.parent())
            False
            sage: HyperbolicPlane(RR).has_coerce_map_from(CC)
            False

        """
        if self.base_ring().has_coerce_map_from(other):
            return True

        if isinstance(other, HyperbolicPlane):
            return self.base_ring().has_coerce_map_from(other.base_ring())

        return False

    def __contains__(self, x):
        r"""
        Return whether the hyperboic plane contains ``x``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        We mostly rely on the standard implementation in SageMath, i.e., the
        interplay between the operator ``==`` and the coercion model::

            sage: 0 in H
            True

        However, we override that logic for complex numbers. There cannot be a
        coercion from the complex number to the hyperbolic plane since that map
        would not be total but we want the following to work::

            sage: I in H
            True

        We do not support such containment checks when conversion of the
        element could lead to a loss of precision::

            sage: CC(I) in H
            False

        .. NOTE::

            There is currently no way to check whether a point is in the
            interior of a set.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.__contains__` to check containment of a
            point in subsets of the hyperbolic plane.

        """
        from sage.categories.all import NumberFields

        import sage.structure.element
        from sage.structure.parent import Parent
        from sage.all import SR

        parent = sage.structure.element.parent(x)
        # Note that in old versions of SageMath (9.1 e.g.), I is not a number field element but a symbolic ring element.
        # The "parent is SR" part can probably removed at some point.
        if isinstance(parent, Parent) and parent in NumberFields() or parent is SR:
            if x.real() in self.base_ring() and x.imag() in self.base_ring() and x.imag() >= 0:
                return True

        return super().__contains__(x)

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

            sage: H = HyperbolicPlane()

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

            interior_points = []
            count = ZZ.random_element().abs() + 3

            while len(interior_points) < count:
                p = self.random_element("point")
                if p.is_ideal():
                    continue
                interior_points.append(p)

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

            sage: H = HyperbolicPlane()

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

    def is_exact(self):
        r"""
        Return whether hyperbolic subsets have exact coordinates.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: H.is_exact()
            True

            sage: H = HyperbolicPlane(RR)
            sage: H.is_exact()
            False

        """
        return self.base_ring().is_exact()

    def infinity(self):
        r"""
        Return the point at infinity in the upper half plane model.

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
        Return the ideal point ``r`` on the real axis in the upper half
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

        EXAMPLES::

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

        return self.geometry.projective(p, q, self.point)

    def start(self, geodesic, check=True):
        r"""
        Return the ideal starting point of ``geodesic``.

        INPUT:

        - ``geodesic`` -- an oriented geodesic

        - ``check`` -- whether to verify that ``geodesic`` is valid

        .. NOTE::

            This method exists to keep all the methods that actually create
            hyperbolic sets on the lowest level in the
            :class:`HyperbolicPlane`. It is otherwise identical to
            :meth:`HyperbolicOrientedGeodesic.start`.

        EXAMPLES::

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

        if check and geodesic.is_ultra_ideal():
            raise ValueError("geodesic does not intersect the Klein disk")

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
            to generate ideal points.

        """
        x = self.base_ring()(x)
        y = self.base_ring()(y)

        if model == "klein":
            point = self.__make_element_class__(HyperbolicPointFromCoordinates)(self, x, y)
        elif model == "half_plane":
            if self.geometry.classify_point(x, y, model="half_plane") < 0:
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
        r"""
        Return the geodesic centered around the real ``center`` and with
        ``radius_squared`` in the upper half plane model. The geodesic is
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
            half plane model and :meth:`geodesic` for the general
            interface, producing a geodesic from an equation.

        """
        center = self.base_ring()(center)
        radius_squared = self.base_ring()(radius_squared)

        return self.geometry.half_circle(center, radius_squared, self.geodesic)

    def vertical(self, real):
        r"""
        Return the vertical geodesic at the ``real`` ideal point in the
        upper half plane model. The geodesic is oriented such that it goes
        from ``real`` to the point at infinity.

        Use the ``-`` operator to pass to the geodesic with opposite
        orientation.

        Use :meth:`HyperbolicConvexSet.unoriented` to get the unoriented
        vertical.

        INPUT:

        - ``real`` -- an element of the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0)
            {-x = 0}

            sage: H.vertical(1)
            {-x + 1 = 0}

            sage: H.vertical(-1)
            {-x - 1 = 0}

        We can also create an unoriented geodesic::

            sage: v = H.vertical(0)
            sage: v.unoriented() == v
            False

        .. SEEALSO::

            :meth:`half_circle` to get an oriented geodesic that is not a
            vertical and :meth:`geodesic` for the general interface, producing
            a geodesic from an equation.

        """
        real = self.base_ring()(real)

        return self.geometry.vertical(real, self.geodesic)

    def geodesic(self, a, b, c=None, model=None, oriented=True, check=True):
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

        INPUT:

        - ``a`` -- a point in the hyperbolic plane or an element of the :meth:`base_ring`

        - ``b`` -- a point in the hyperbolic plane or an element of the :meth:`base_ring`

        - ``c`` -- ``None`` or an element of the :meth:`base_ring` (default: ``None``)

        - ``model`` -- ``None``, ``"half_plane"``, or ``"klein"`` (default:
          ``None``); when ``a``, ``b`` and ``c`` are elements of the
          :meth:`base_ring`, in which model they should be interpreted.

        - ``oriented`` -- whether the returned geodesic is oriented (default: ``True``)

        - ``check`` -- whether to verify that the arguments define a geodesic
          in the hyperbolic plane (default: ``True``)

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

        Geodesics cannot be defined from points whose coordinates are over a
        quadratic field extension::

            sage: H.geodesic(H.half_circle(0, 2).start(), H.half_circle(1, 2).end())
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not in Rational Field

        Except for some special cases::

            sage: H.geodesic(H.half_circle(0, 2).start(), H.half_circle(0, 2).end())
            {(x^2 + y^2) - 2 = 0}

        Disabling the ``check``, lets us define geodesics in the Klein model
        that lie outside the unit circle::

            sage: H.geodesic(2, 1, 0, model="klein")
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: geodesic = H.geodesic(2, 1, 0, model="klein", check=False)
            sage: geodesic
            {2 + x = 0}

            sage: geodesic.start()
            Traceback (most recent call last):
            ...
            ValueError: geodesic does not intersect the Klein disk

            sage: geodesic.end()
            Traceback (most recent call last):
            ...
            ValueError: geodesic does not intersect the Klein disk

        TESTS::

            sage: H.geodesic(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: points specifying a geodesic must be distinct

        ..SEEALSO::

            :meth:`half_circle` and :meth:`vertical`

        """
        if c is None:
            a = self(a)
            b = self(b)

            if a == b:
                raise ValueError("points specifying a geodesic must be distinct")

            if isinstance(a, HyperbolicPointFromGeodesic) and isinstance(b, HyperbolicPointFromGeodesic) and a._geodesic == -b._geodesic:
                return a._geodesic

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

        INPUT:

        - ``a`` -- an element of the :meth:`base_ring`

        - ``b`` -- an element of the :meth:`base_ring`

        - ``c`` -- an element of the :meth:`base_ring`

        - ``model`` -- one of ``"half_plane"``` or ``"klein"``

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_space(0, -1, 0, model="half_plane")
            {x ≤ 0}

        It is often easier to construct a half space as the space bounded by a geodesic::

            sage: H.vertical(0).left_half_space()
            {x ≤ 0}

        The half space y ≥ 0 given by its equation in the Klein model::

            sage: H.half_space(0, 0, 1, model="klein")
            {(x^2 + y^2) - 1 ≥ 0}

        ..SEEALSO::

            :meth:`HperbolicGeodesic.left_half_space`
            :meth:`HperbolicGeodesic.right_half_space`

        """
        geodesic = self.geodesic(a, b, c, model=model, check=check)

        return self.__make_element_class__(HyperbolicHalfSpace)(self, geodesic)

    def segment(
        self,
        geodesic,
        start=None,
        end=None,
        oriented=None,
        check=True,
        assume_normalized=False,
    ):
        r"""
        Return the segment on the ``geodesic`` bounded by ``start`` and ``end``.

        INPUT:

        - ``geodesic`` -- a :meth:`geodesic` in this space.

        - ``start`` -- ``None`` or a :meth:`point` on the ``geodesic``, e.g.,
          obtained from the :meth:`HyperbolicOrientedGeodesic.intersection` of
          ``geodesic`` with another geodesic. If ``None``, the segment starts
          at the infinite :meth:`HyperbolicOrientedGeodesic.start` point of the
          geodesic.

        - ``end`` -- ``None`` or a :meth:`point` on the ``geodesic``, as for
          ``start``; must be later on ``geodesic`` than ``start`` if the
          geodesic is oriented.

        - ``oriented`` -- whether to produce an oriented segment or an
          unoriented segment. The default (``None``) is to produce an oriented
          segment iff ``geodesic`` is oriented or both ``start`` and ``end``
          are provided so the orientation can be deduced from their order.

        - ``check`` -- boolean (default: ``True``), whether validation is
          performed on the arguments.

        - ``assume_normalized`` -- boolean (default: ``False``), if not set,
          the returned segment is normalized, i.e., if it is actually a
          geodesic, a :class:`HyperbolicGeodesic` is returned, if it is
          actually a point, a :class:`HyperbolicPoint` is returned.

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
            sage: H.segment(H.vertical(0), start=I) != H.segment(H.vertical(0), end=I)
            True

        TESTS:

        When only a ``start`` point is provided, we cannot deduce the orientation of the geodesic::

            sage: H.segment(H.vertical(0).unoriented(), start=0, oriented=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

            sage: H.segment(H.vertical(0).unoriented(), start=0, oriented=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

            sage: H.segment(H.vertical(0).unoriented(), start=I, oriented=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

        When only an ``end`` point is provided, we cannot deduce the orientation of the geodesic::

            sage: H.segment(H.vertical(0).unoriented(), end=0, oriented=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

            sage: H.segment(H.vertical(0).unoriented(), end=0, oriented=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

            sage: H.segment(H.vertical(0).unoriented(), end=I, oriented=True)
            Traceback (most recent call last):
            ...
            ValueError: cannot deduce segment from single endpoint on an unoriented geodesic

        When ``start`` and ``end`` are given, they must be ordered correctly::

            sage: H.segment(H.vertical(0), start=0, end=oo)
            {-x = 0}

            sage: H.segment(H.vertical(0).unoriented(), start=oo, end=0)
            {x = 0}

            sage: H.segment(H.vertical(0), start=oo, end=0)
            Traceback (most recent call last):
            ...
            ValueError: end point of segment must be after start point on the underlying geodesic

            sage: H.segment(H.vertical(0), start=I, end=2*I)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}

            sage: H.segment(H.vertical(0).unoriented(), start=2*I, end=I)
            {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

            sage: H.segment(H.vertical(0), start=2*I, end=I)
            Traceback (most recent call last):
            ...
            ValueError: end point of segment must be after start point on the underlying geodesic

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

            if start is None and end is None:
                # any orientation of the geodesic will do
                pass
            elif start is None or end is None or start == end:
                raise ValueError("cannot deduce segment from single endpoint on an unoriented geodesic")
            elif geodesic.parametrize(start, model="euclidean", check=False) > geodesic.parametrize(end, model="euclidean", check=False):
                geodesic = -geodesic

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
        self, half_spaces, check=True, assume_sorted=False, assume_minimal=False, marked_vertices=()
    ):
        r"""
        Return the convex polygon obtained by intersecting ``half_spaces``.

        INPUT:

        - ``half_spaces`` -- a non-empty iterable of
          :class:`HyperbolicHalfSpace`\ s of this hyperbolic plane.

        - ``check`` -- boolean (default: ``True``), whether the arguments are
          validated.

        - ``assume_sorted`` -- boolean (default: ``False``), whether to assume
          that the ``half_spaces`` are already sorted with respect to
          :meth:`HyperbolicHalfSpace._less_than`. When set, we omit sorting the
          half spaces explicitly, which is asymptotically the most exponsive
          part of the process of creating a polygon.

        - ``assume_minimal`` -- boolean (default: ``False``), whether to assume
          that the ``half_spaces`` provide a minimal representation of the
          polygon, i.e., removing any of them describes a different polygon.
          When set, we omit searching for a minimal subset of half spaces to
          describe the polygon.

        - ``marked_vertices`` -- an iterable of vertices (default: an empty
          tuple), the vertices are included in the
          :meth:`HyperbolicConvexPolygon.vertices` even if they are not in the set of
          minimal vertices describing this polygon.

        ALGORITHM:

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
        meet at ideal points::

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

        If we add a marked point to such a half space, the underlying type is a
        polygon again::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[I])
            sage: half_space
            {(x^2 + y^2) - 1 ≤ 0} ∪ {I}
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicConvexPolygon_with_category'>

        Marked points that coincide with vertices are ignored::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[-1])
            sage: half_space
            {(x^2 + y^2) - 1 ≤ 0}
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicHalfSpace_with_category'>

        Marked points must be on an edge of the polygon::

            sage: H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[-2])
            Traceback (most recent call last):
            ...
            ValueError: marked vertex must be on an edge of the polygon

            sage: H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[2*I])
            Traceback (most recent call last):
            ...
            ValueError: marked vertex must be on an edge of the polygon

        The intersection of the half spaces is computed in time quasi-linear in
        the number of half spaces. The limiting factor is sorting the half
        spaces by :meth:`HyperbolicHalfSpace._less_than`. If we know that the
        half spaces are already sorted like that, we can make the process run
        in linear time by setting ``assume_sorted``.

            sage: H.polygon(H.infinity().half_spaces(), assume_sorted=True)
            ∞

        .. SEEALSO::

            :meth:`intersection` to intersect arbitrary convex sets
            :meth:`convex_hull` to define a polygon by taking the convex hull
            of a union of convex sets

        """
        if not marked_vertices:
            marked_vertices = []

        half_spaces = [self(half_space) for half_space in half_spaces]
        marked_vertices = [self(vertex) for vertex in marked_vertices]

        half_spaces = HyperbolicHalfSpaces(half_spaces, assume_sorted=assume_sorted)

        polygon = self.__make_element_class__(HyperbolicConvexPolygon)(
            self, half_spaces, marked_vertices
        )

        if check:
            polygon._check(require_normalized=False)

        if check or not assume_minimal:
            polygon = polygon._normalize(marked_vertices=bool(marked_vertices))

        if check:
            polygon._check()

        return polygon

    def convex_hull(self, *subsets, marked_vertices=False):
        r"""
        Return the convex hull of the ``subsets``.

        INPUT:

        - ``subsets`` -- a sequence of subsets of this hyperbolic space.

        - ``marked_vertices`` -- a boolean (default: ``False``), whether to
          keep redundant vertices on the boundary.

        ALGORITHM:

        We use the standard Graham scan algorithm which runs in O(nlogn), see
        :meth:`HyperbolicHalfSpaces.convex_hull`.

        However, to get the unbounded bits of the convex hull right, we use a
        somewhat naive O(n²) algorithm which could probably be improved easily.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A polygon can also be created as the convex hull of its vertices::

            sage: H.convex_hull(I - 1, I + 1, 2*I - 1, 2*I + 1)
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 5 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        The vertices can also be infinite::

            sage: H.convex_hull(-1, 1, 2*I)
            {(x^2 + y^2) + 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        Redundant vertices are removed. However, they can be kept by setting
        ``marked_vertices``::

            sage: H.convex_hull(-1, 1, I, 2*I)
            {(x^2 + y^2) + 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

            sage: polygon = H.convex_hull(-1, 1, I, 2*I, marked_vertices=True)
            sage: polygon
            {(x^2 + y^2) + 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 3*x - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∪ {I}

            sage: polygon.vertices()
            {-1, I, 1, 2*I}

        The convex hull of a half space and a point::

            sage: H.convex_hull(H.half_circle(0, 1).right_half_space(), 0)
            {(x^2 + y^2) - 1 ≤ 0}

        To keep the additional vertices, again ``marked_vertices`` must be set::

            sage: H.convex_hull(H.half_circle(0, 1).left_half_space(), I, marked_vertices=True)
            {(x^2 + y^2) - 1 ≥ 0} ∪ {I}

            sage: H.convex_hull(H.vertical(0).left_half_space(), 0, I, oo, marked_vertices=True)
            {x ≤ 0} ∪ {I}

            sage: H.convex_hull(H.vertical(0).right_half_space(), I, marked_vertices=True)
            {x ≥ 0} ∪ {I}

        Note that this cannot be used to produce marked points on a geodesic::

            sage: H.convex_hull(-1, I, 1)
            {(x^2 + y^2) - 1 = 0}

            sage: H.convex_hull(-1, I, 1, marked_vertices=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot add marked vertices to low dimensional objects

        Note that this cannot be used to produce marked points on a segment::

            sage: H.convex_hull(I, 2*I, 3*I)
            {x = 0}

            sage: H.convex_hull(I, 2*I, 3*I, marked_vertices=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot add marked vertices to low dimensional objects

        The convex point of two polygons which contain infinitely many ideal points::

            sage: H.convex_hull(
            ....:     H.polygon([H.geodesic(-1, -1/4).left_half_space(), H.geodesic(0, -2).left_half_space()]),
            ....:     H.polygon([H.geodesic(4, 2).left_half_space(), H.geodesic(4, oo).left_half_space()]),
            ....:     H.polygon([H.geodesic(-1/2, 1/2).left_half_space(), H.geodesic(2, -2).left_half_space()])
            ....: )
            {2*(x^2 + y^2) - x ≥ 0} ∩ {(x^2 + y^2) - 2*x - 8 ≤ 0} ∩ {8*(x^2 + y^2) + 6*x + 1 ≥ 0}

        .. SEEALSO::

            :meth:`HyperbolicHalfSpaces.convex_hull` for the underlying implementation
            :meth:`intersection` to compute the intersection of convex sets

        """
        subsets = [self(subset) for subset in subsets]

        vertices = sum([list(subset.vertices()) for subset in subsets], [])

        polygon = self.polygon(HyperbolicHalfSpaces.convex_hull(vertices))

        half_spaces = []
        for subset in subsets:
            if subset.dimension() == 2:
                if isinstance(subset, HyperbolicHalfSpace):
                    half_spaces.append(subset)
                elif isinstance(subset, HyperbolicConvexPolygon):
                    # An infinity polygon is more than just the convex hull of
                    # its vertices, it may also contain entire half spaces,
                    # namely those that are between vertices that are not
                    # connected by an edge.
                    edges = subset.edges()
                    for (a, b) in subset.vertices().pairs():
                        if a.segment(b) not in edges:
                            half_spaces.append(self.geodesic(a, b).right_half_space())
                else:
                    raise NotImplementedError("cannot form convex hull of this kind of set yet")

        edges = []

        for edge in polygon.half_spaces():
            for half_space in half_spaces:
                if (-edge).is_subset(half_space):
                    break
            else:
                edges.append(edge)

        if marked_vertices:
            marked_vertices = [vertex for vertex in vertices if any(vertex in half_space.boundary() for half_space in edges)]

        polygon = self.polygon(edges, marked_vertices=marked_vertices)

        assert all(subset.is_subset(polygon) for subset in subsets), "convex hull does not contain all the sets is supposed to be the convex hull of"

        return polygon

    def intersection(self, *subsets):
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

        INPUT:

        - ``subsets`` -- a non-empty sequence of subsets of this hyperbolic space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.intersection(H.vertical(0).left_half_space())
            {x ≤ 0}

            sage: H.intersection(H.vertical(0).left_half_space(), H.vertical(0).right_half_space())
            {x = 0}

        We cannot form the intersection of no spaces yet::

            sage: H.intersection()
            Traceback (most recent call last):
            ...
            NotImplementedError: the full hyperbolic space cannot be created as an intersection

        .. SEEALSO::

            :meth:`HyperbolicPlane.polygon` for a specialized version for the intersection of half spaces
            :meth:`HyperbolicPLane.convex_hull` to compute the convex hull of subspaces

        """
        subsets = [self(subset) for subset in subsets]

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
        r"""
        Return an empty subset of this space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: HyperbolicPlane().empty_set()
            {}

        """
        return self.__make_element_class__(HyperbolicEmptySet)(self)

    def isometry(self, preimage, image, model="half_plane", on_right=False, normalized=False):
        r"""
        Return an isometry that maps ``preimage`` to ``image``.

        INPUT:

        - ``preimage`` -- a convex set in the hyperbolic plane or a list of
          such convex sets.

        - ``image`` -- a convex set in the hyperbolic plane or a list of such
          convex sets.

        - ``model`` -- one of ``"half_plane"`` and ``"klein"``, the model in
          which this isometry applies.

        - ``on_right`` -- a boolean (default: ``False``); whether the returned
          isometry maps ``preimage`` to ``image`` when multiplied from the
          right; otherwise from the left.

        - ``normalized`` -- a boolean (default: ``False``); whether the
          returned matrix has determinant ±1.

        OUTPUT:

        If ``model`` is ``"half_plane"``, returns a 2×2 matrix over the
        :meth:`base_ring`, if ``model`` is ``"klein"``, returns a 3×3 matrix
        over the base ring.
        See :meth:`apply_isometry` for meaning of this matrix.

        ALGORITHM:

        We compute an isometry with a very inefficient Gröbner basis approach.

        Essentially, we try to extract three points from the ``preimage`` with
        their prescribed images in ``image``, see :meth:`_isometry_conditions`
        and determine the unique isometry mapping the points by solving the
        corresponding polynomial system, see :meth:`_isometry_from_equations`
        for the hacky Gröbner basis bit.

        There are a lot of problems with this approach (apart from it being
        extremely slow.)

        Usually, we do not have access to meaningful points (the ideal end
        points of a geodesic do not typically live in the :meth:`base_ring`) so
        we have to map around geodesics instead. However, given two pairs of
        geodesics, there is in general not a unique isometry mapping one pair
        to the other, since there might be one isometry with positive and one
        with negative determinant with this property. This adds to the
        inefficiency because we have to try for both determinants and then
        check which isometry maps the actual objects in question correctly.

        Similarly, for more complex objects such as polygons, we do not know
        a-priori which edge of the preimage polygon can be mapped to which edge
        of the image polygon, so we have to try for all rotations of both
        polygons.

        Finally, the approach cannot work in general for certain
        undeterdetermined systems. We do not know how to determine an isometry
        that maps one geodesic to another geodesic. We can of course try to
        write down an isometry that maps a geodesic to the vertical at 0 but in
        general no such isometry is defined over the base ring. We can make the
        system determined by requiring that the midpoints of the geodesics map
        to each other but (apart from the midpoint not having coordinates in
        the base ring) that isometry might not be defined over the base ring
        even though there exists some isometry that maps the geodesics over the
        base ring.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        An isometry mapping one point to another; naturally this is not the
        unique isometry with this property::

            sage: m = H.isometry(0, 1)
            sage: H(0).apply_isometry(m)
            1
            sage: m
            [1 1]
            [0 1]

        ::

            sage: H.isometry(I, I+1)
            [1 1]
            [0 1]
            sage: H.isometry(0, oo)
            [ 0 -1]
            [-1  0]

        An isometry is uniquely determined by its image on three points::

            sage: H.isometry([0, 1, oo], [1, oo, 0])
            [ 0  1]
            [-1  1]

        It might be impossible to find an isometry with the prescribed mapping::

            sage: H.isometry(0, I)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

            sage: H.isometry([0, 1, oo, I], [0, 1, oo, I + 1])
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

        We can determine isometries by mapping more complex objects than
        points, e.g., geodesics::

            sage: H.isometry(H.geodesic(-1, 1), H.geodesic(1, -1))
            [ 0  1]
            [-1  0]

        We can also determine an isometry mapping polygons::

            sage: P = H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space(), H.geodesic(-1, 1).left_half_space()])
            sage: Q = H.polygon([H.geodesic(-1, 0).left_half_space(), H.geodesic(0, 1).left_half_space(), H.geodesic(1, -1).left_half_space()])
            sage: m = H.isometry(P, Q)
            sage: P.apply_isometry(m) == Q
            True

        When determining an isometry of polygons, marked vertices are mapped to
        marked vertices::

            sage: P = H.polygon(P.half_spaces(), marked_vertices=[1 + I])
            sage: Q = H.polygon(P.half_spaces(), marked_vertices=[])
            sage: H.isometry(P, Q)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

            sage: Q = H.polygon(P.half_spaces(), marked_vertices=[1 + 2*I])
            sage: H.isometry(P, Q)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

            sage: Q = H.polygon(P.half_spaces(), marked_vertices=[-1 + I])
            sage: H.isometry(P, Q)
            [ 1  0]
            [ 0 -1]

        We can explicitly ask for an isometry in the Klein model, given by a
        3×3 matrix::

            sage: H.isometry(P, Q, model="klein")
            [-1  0  0]
            [ 0  1  0]
            [ 0  0  1]

        The isometries are not returned as matrices of unit determinant since
        such an isometry might not exist without extending the base ring, we
        can, however, ask for an isometry of determinant ±1::

            sage: H.isometry(I, 2*I, normalized=True)
            Traceback (most recent call last):
            ...
            ValueError: not a perfect 2nd power

            sage: H.change_ring(AA).isometry(I, 2*I, normalized=True)
            [ 1.414213562373095?                   0]
            [                  0 0.7071067811865475?]

            sage: _.det()
            1.000000000000000?

        We can also explicitly ask for the isometry for the right action::

            sage: isometry = H.isometry(H.vertical(0), H.vertical(1), on_right=True)
            sage: isometry
            [ 1 -1]
            [ 0  1]
            sage: H.vertical(0).apply_isometry(isometry)
            {-x - 1 = 0}
            sage: H.vertical(0).apply_isometry(isometry, on_right=True)
            {-x + 1 = 0}

        TESTS:

        Points forming an ideal triangle with three ideal vertices::

            sage: preimage = H(-1), H(0), H(1)
            sage: image = H(0), H(1), H(2)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Points forming an ideal triangle with two ideal vertices::

            sage: preimage = H(-1), H(2*I), H(1)
            sage: image = H(0), H(1 + 2*I), H(2)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Points forming an ideal triangle with one ideal vertex::

            sage: preimage = H(I - 1), H(I + 1), H(oo)
            sage: image = H(I), H(I + 2), H(oo)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Points forming a finite triangle::

            sage: preimage = H(I), H(2*I), H(2*I + 1)
            sage: image = H(I + 1), H(2*I + 1), H(2*I + 2)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Points forming a degenerate ideal triangle::

            sage: preimage = H(-1), H(I), H(1)
            sage: image = H(0), H(I + 1), H(2)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Another degenerate ideal triangle::

            sage: preimage = H(0), H(I), H(oo)
            sage: image = H(1), H(I + 1), H(oo)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Points forming a finite degenerate triangle::

            sage: preimage = H(I), H(2*I), H(3*I)
            sage: image = H(I + 1), H(2*I + 1), H(3*I + 1)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        Impossible pairs of points (ideal and non-ideal) are detected::

            sage: preimage = H(-1), H(I), H(1)
            sage: image = H(0), H(1), H(2)
            sage: H.isometry(preimage, image)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

            sage: preimage = H(0), H(1), H(2)
            sage: image = H(-1), H(I), H(1)
            sage: H.isometry(preimage, image)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

        A case with d=0::

            sage: preimage = (H.geodesic(0, 1), H.geodesic(1, oo))
            sage: image = (H.geodesic(1, oo), H.geodesic(oo, 0))
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 1  0]

        A case with no solution, such an isometry would map four ideal points
        in an impossible way::

            sage: preimage = (H.geodesic(0, 1), H.geodesic(2, 3))
            sage: image = (H.geodesic(0, 1), H.geodesic(3, 4))
            sage: H.isometry(preimage, image, on_right=True)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

        An isometry that swaps end points but maps the corresponding oriented
        geodesics to themselves::

            sage: preimage = (-1, 1, -2, 2)
            sage: image = (1, -1, 2, -2)
            sage: H.isometry(preimage, image, on_right=True)
            [-1  0]
            [ 0  1]

        An example with negative determinant::

            sage: preimage = (I - 1, I + 1, I - 2, I + 2)
            sage: image = (I + 1, I - 1, I + 2, I - 2)
            sage: H.isometry(preimage, image, on_right=True)
            [-1  0]
            [ 0  1]

        A case that could initially not be solved over the rationals::

            sage: H.isometry((58*I, I + 1), (116*I - 1, 2*I + 1))
            [ 2 -1]
            [ 0  1]

        A case that caused problems at some point::

            sage: isometry = matrix([[1, 0], [0, -1/2]])
            sage: x = H(I/2 - 1)
            sage: y = H(I/3 - 1)
            sage: z = H(5/19 * I - 1/3)
            sage: H.isometry((x, y, z), (x.apply_isometry(isometry), y.apply_isometry(isometry), z.apply_isometry(isometry)))
            [-2  0]
            [ 0  1]

        ::

            sage: preimage = (I - 1, I + 1, I + 1, oo)
            sage: image = (I, I + 2, I + 2, oo)
            sage: H.isometry(preimage, image, on_right=True)
            [ 1 -1]
            [ 0  1]

        ::

            sage: preimage = (H.geodesic(I - 1, I + 1), H.geodesic(I - 2, I + 1))
            sage: image = (H.geodesic(I + 1, I - 1), H.geodesic(I + 1, I - 2))
            sage: H.isometry(preimage, image, on_right=True)
            [-1  2]
            [-1  1]

        ::

            sage: isometry = matrix([[2, 0], [1, 2]])
            sage: x = H(0)
            sage: y = H(I/2 - 1)
            sage: z = H(0)
            sage: H.isometry((x, y, z), (x.apply_isometry(isometry), y.apply_isometry(isometry), z.apply_isometry(isometry)))
            [  1   0]
            [1/2   1]

        ::

            sage: H.isometry(I, 2*I)
            [  1   0]
            [  0 1/2]

        An underdetermined case::

            sage: P = H.geodesic(126, 4447, 6387, model="half_plane").right_half_space()
            sage: isometry = matrix([[-1, 2], [2, -1/2]])
            sage: Q = P.apply_isometry(isometry)
            sage: H.isometry(P, Q)
            [  126/8579 -4321/8579]
            [         0          1]

        Here, there is also an isometry of negative determinant that maps some
        of the half spaces correctly::

            sage: P = H.polygon([
            ....:    H.geodesic(1, -5, 5, model="half_plane").left_half_space(),
            ....:    H.geodesic(247, -957, -5156, model="half_plane").right_half_space(),
            ....:    H.geodesic(1, -120, -137, model="half_plane").right_half_space()])
            sage: isometry = matrix([[0, -2], [1, 2]])
            sage: Q = P.apply_isometry(isometry)
            sage: H.isometry(P, Q)
            [  0  -1]
            [1/2   1]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.apply_isometry` to apply the returned
            isometry to a convex set.

        """
        if normalized:
            isometry = self.isometry(preimage=preimage, image=image, model=model, on_right=on_right, normalized=False)
            det = abs(isometry.det())
            λ = det.nth_root(isometry.nrows())
            return ~λ * isometry

        if model == "klein":
            isometry = self.isometry(preimage=preimage, image=image, model="half_plane", on_right=on_right, normalized=normalized)
            return gl2_to_sim12(isometry)
        elif model != "half_plane":
            raise NotImplementedError("unsupported model")

        if not on_right:
            isometry = self.isometry(preimage=preimage, image=image, model=model, on_right=True, normalized=normalized)
            if model == "half_plane":
                from sage.all import matrix
                # Pick a nice representative of the inverse matrix.
                isometry = matrix([[isometry[1][1], -isometry[0][1]], [-isometry[1][0], isometry[0][0]]])
            else:
                isometry = ~isometry

            return isometry

        # Normalize the arguments so that they are a list of convex sets.
        from collections.abc import Iterable
        if not isinstance(preimage, Iterable):
            preimage = [preimage]
        if not isinstance(image, Iterable):
            image = [image]

        preimage = [self(x) for x in preimage]
        image = [self(y) for y in image]

        if len(preimage) != len(image):
            raise ValueError("preimage and image must be the same size to determine an isometry between them")

        for (x, y) in zip(preimage, image):
            if x.dimension() != y.dimension():
                raise ValueError("preimage and image must be of the same dimensions to determine an isometry between them")
            if x.is_oriented() != y.is_oriented():
                raise ValueError("preimage and image must be oriented or unoriented consistently")

        # Drop empty sets from the preimage and image to make our lives easier
        preimage = [x for x in preimage if x.dimension() >= 0]
        image = [y for y in image if y.dimension() >= 0]

        # An isometry is uniquely determined by mapping three points.
        # In principle, we now just need to find three points in preimage, find
        # their corresponding images, construct the isometry, and then check
        # that it maps everything correctly.
        # However, there are objects that do not really define a mapping of
        # points. For example, when presented with an unoriented geodesic, we
        # can map its endpoints two possible ways (apart from that, an isometry
        # can swap the end points of a geodesic but map an oriented geodesic to
        # itself at the same time.) Similarly, when mapping a polygon, we can
        # permute the edges cyclically.

        # We need a mild form of backtracking to collect all possible triples
        # that define the isometry.
        def search_isometry(pairs):
            for conditions in self._isometry_conditions([], pairs):
                for isometry in self._isometry_from_primitives(conditions):
                    if isometry is None:
                        continue

                    if any([preimage.apply_isometry(isometry, on_right=True) != image for (preimage, image) in pairs]):
                        continue

                    return isometry

        isometry = search_isometry(list(zip(preimage, image)))

        if isometry is None:
            raise ValueError("no isometry can map these objects to each other")

        return isometry

    def _isometry_from_primitives(self, pairs):
        r"""
        Helper method for :meth:`isometry`.

        Return right isometries as 2x2 matrices that maps the elements of
        ``pairs`` to each other.

        INPUT:

        - ``pairs`` -- a sequence of pairs of geodesics or hyperbolic points

        OUTPUT:

        An iterator of matrices, see below.

        ALGORITHM:

        If ``pairs`` is a single pair of points, we construct an isometry with
        :meth:`_isometry_from_single_points`.

        If ``pairs`` is a single pair of geodesics, we make an attempt to
        construct an isometry with :meth:`_isometry_from_single_geodesics`.

        Otherwise, ``pairs`` up to three pairs of points and geodesics. These
        do not always uniquely define an isometry, e.g., often when presented
        with two geodesics. Namely, there could be one isometry of positive
        determinant and one isometry of negative determinant. We return both of
        them.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        An underdetermined system. Here we are lucky, the midpoint of the
        geodesics has coordinates in the base ring and the isometry mapping
        midpoints to each other is defined over the base ring::

            sage: preimage = H.geodesic(-126, -4447, -6387, model="half_plane")
            sage: image = H.geodesic(-8579, -13089, -4510, model="half_plane")
            sage: list(H._isometry_from_primitives([(preimage, image)]))
            [
            [        1 4321/8579]
            [        0  126/8579]
            ]

        Here, a non-trivial isometry maps these oriented geodesics to
        themselves, i.e., it stabilizes what's to the left of the geodesics::

           sage: g = H.geodesic(-1, 1)
           sage: h = H.geodesic(-2, 2)
           sage: g.apply_isometry(matrix([[-1, 0], [0, 1]])) == g
           True
           sage: h.apply_isometry(matrix([[-1, 0], [0, 1]])) == h
           True
           sage: list(H._isometry_from_primitives([(g, g), (h, h)]))
           [
           [1 0]  [-1  0]
           [0 1], [ 0  1]
           ]

        Note that that isometry swaps endpoints though::

            sage: H(-1).apply_isometry(matrix([[-1, 0], [0, 1]]))
            1

        """
        if len(pairs) == 0:
            from sage.all import matrix
            yield matrix(self.base_ring(), [[1, 0], [0, 1]])
            return

        if len(pairs) == 1 and pairs[0][0].dimension() == 0:
            yield self._isometry_from_single_points(pairs[0][0], pairs[0][1])
            return

        if len(pairs) == 1 and pairs[0][0].dimension() == 1:
            yield self._isometry_from_single_geodesics(pairs[0][0], pairs[0][1])
            return

        from sage.all import vector

        # Create polynomial equations that must be satisfied to map the pairs
        # to each other.
        def equations(isometry, λ):
            R = isometry.base_ring()

            equations = []

            for i, (preimage, image) in enumerate(pairs):
                equations.extend(preimage._isometry_equations(isometry, image, λ[i]))

            return equations

        # Create a predicate that can be used to check whether the isometry
        # correctly maps the pairs.
        def create_filter(det):
            def filter(isometry):
                if isometry.det().sign() != det:
                    return False
                for preimage, image in pairs:
                    if preimage.apply_isometry(isometry, on_right=True) != image:
                        return False
                return True

            return filter

        yield self._isometry_from_equations(equations, create_filter(1))
        yield self._isometry_from_equations(equations, create_filter(-1))

    def _isometry_untrivialize(self, preimage, image, defining):
        r"""
        Helper method for :meth:`isometry`.

        Return a pair of hyperbolic objects that describe the mapping of
        ``preimage`` to ``image`` or return ``None`` if that mapping is already
        captured by the pairs of objects in ``defining``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H._isometry_untrivialize(H(I), H(I), [])
            (I, I)

        There are many ways to map the geodesic between -1 and 1 to itself.
        Knowing that we need to fix `I` adds information::

            sage: H._isometry_untrivialize(H(I), H(I), [(H.geodesic(-1, 1), H.geodesic(-1, 1))])
            (I, I)

        However, we already know that ``-1`` must go to ``-1``. Well, actually
        that's not true, we could also be swapping the endpoints otherwise but
        we want to use this to build three conditions that are independent
        (this could probably be improved)::

            sage: H._isometry_untrivialize(H(-1), H(-1), [(H.geodesic(-1, 1), H.geodesic(-1, 1))]) is None
            True

        """
        existings = [x for (x, y) in defining]

        if preimage.dimension() == 0:
            if preimage in existings:
                return None
            if preimage.is_ideal():
                # Actually, an ideal point is not trivial if the existing
                # geodesic is unoriented. But it does not take a full
                # degree of freedom away, so we ignore it here.
                if any([preimage in existing for existing in existings]):
                    return None
            return (preimage, image)

        elif preimage.dimension() == 1:
            if preimage.unoriented() in [existing.unoriented() for existing in existings]:
                # Again, we ignore the distinction between oriented and
                # unoriented geodesics here.
                return None

            if preimage.start() in existings:
                return self._isometry_untrivialize((preimage.end(), image.end()))

            if preimage.end() in existings:
                return self._isometry_untrivialize((preimage.start(), image.start()))

            return (preimage, image)

        raise NotImplementedError

    def _isometry_conditions(self, defining, remaining):
        r"""
        Helper method for :meth:`isometry`.

        Return sequences (typically triples) of pairs of hyperbolic primitive
        objects (geodesics and points) that (almost) uniquely define a
        hyperbolic mapping.

        Build this sequence by extending ``defining`` with conditions extracted
        from ``remaining``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        Given four points, three points already uniquely determine the isometry::

            sage: conditions = H._isometry_conditions(defining=[], remaining=[(H(0), H(0)), (H(1), H(1)), (H(2), H(2)), (H(3), H(3))])
            sage: list(conditions)
            [[(0, 0), (1, 1), (2, 2)]]

        The data provided by ``remaining`` can contain redundancies::

            sage: conditions = H._isometry_conditions(defining=[], remaining=[(H(0), H(0)), (H(1), H(1)), (H(1), H(1)), (H(3), H(3))])
            sage: list(conditions)
            [[(0, 0), (1, 1), (3, 3)]]

        For more complex objects there might be lots of mappings possible (we
        could likely have a shorter list of possibilities here)::

            sage: P = H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space(), H.geodesic(-1, 1).left_half_space()], marked_vertices=[I + 1])
            sage: Q = H.polygon(P.half_spaces(), marked_vertices=[I - 1])
            sage: conditions = H._isometry_conditions([], [(P, Q)])
            sage: list(conditions)
            [[(-1, -1), (1, 1), (1 + I, ∞)],
             [(-1, 1), (1, ∞), (1 + I, -1 + I)],
             [(-1, ∞), (1, -1 + I), (1 + I, -1)],
             [(-1, -1 + I), (1, -1), (1 + I, 1)],
             [(-1, -1 + I), (1, ∞), (1 + I, 1)],
             [(-1, ∞), (1, 1), (1 + I, -1)],
             [(-1, 1), (1, -1), (1 + I, -1 + I)],
             [(-1, -1), (1, -1 + I), (1 + I, ∞)]]

        """
        def degree(preimage, image):
            if preimage.dimension() == 0:
                return 1
            if preimage.dimension() == 1:
                return 2
            assert False

        degree = sum([degree(preimage, image) for (preimage, image) in defining])

        # If we have three pairs of points, determine the unique isometry that maps them to each other.
        if degree >= 3:
            yield defining
            return

        # There are fewer than three points in "defining". Extend with more points.
        if remaining:
            # Extend by turning remaining[0] into a condition
            x, y = remaining[0]
            remaining = remaining[1:]

            # Extend with a pair of points in "remaining[0]"
            if x.dimension() == 0:
                assert y.dimension() == 0

                if pair := self._isometry_untrivialize(x, y, defining):
                    defining.append(pair)

                for conditions in self._isometry_conditions(defining[:], remaining):
                    yield conditions

            # Extend with a pair of geodesics in "remaining[0]"
            elif x.dimension() == 1:
                assert y.dimension() == 1

                f = x.geodesic()
                g = y.geodesic()

                if pair := self._isometry_untrivialize(f, g, defining):
                    defining.append(pair)

                for conditions in self._isometry_conditions(defining[:], remaining + [(x.start(), y.start()), (x.end(), y.end())]):
                    yield conditions

            # Extend with points coming from other hyperbolic objects in "remaining[0]"
            else:
                for pairs in x._isometry_conditions(y):
                    for conditions in self._isometry_conditions(defining[:], pairs + remaining):
                        yield conditions

        else:
            yield defining

    def _isometry_from_single_points(self, preimage, image):
        r"""
        Helper method for :meth:`isometry`.

        Return a right isometry that maps the point ``preimage`` to the point
        ``image`` or ``None`` when no such isometry exists.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H._isometry_from_single_points(H(I), H(I))
            [1 0]
            [0 1]

            sage: H._isometry_from_single_points(H(I), H(2*I))
            [1/2   0]
            [  0   1]

            sage: H._isometry_from_single_points(H(0), H(oo))
            [0 1]
            [1 0]

            sage: H._isometry_from_single_points(H(I), H(1)) is None
            True

        """
        if preimage.is_ideal() != image.is_ideal():
            return None

        from sage.all import MatrixSpace
        MS = MatrixSpace(self.base_ring(), 2, 2)

        isometry = MS(1)

        if preimage == image:
            return isometry

        if preimage.is_ideal():
            if preimage == self.infinity():
                isometry *= ~MS([[0, 1], [1, 0]])
            else:
                isometry *= ~MS([[1, -preimage.coordinates()[0]], [0, 1]])

            if image == self.infinity():
                isometry *= ~MS([[0, 1], [1, 0]])
            else:
                isometry *= ~MS([[1, image.coordinates()[0]], [0, 1]])

            return isometry
        else:
            isometry *= ~MS([[image.coordinates()[1] / preimage.coordinates()[1], 0], [0, 1]])
            isometry *= ~MS([[1, image.coordinates()[0] - preimage.apply_isometry(isometry, on_right=True).coordinates()[0]], [0, 1]])
            return isometry

    def _isometry_from_single_geodesics(self, preimage, image):
        r"""
        Helper method for :meth:`isometry`.

        Return a right isometry that maps the geodesic ``preimage`` to the
        geodesic ``image`` or ``None`` when no such isometry exists.

        ALGORITHM:

        We determine the isometry by forcing the midpoints of the geodesics to
        be mapped to each other. This might fail because the midpoints of the
        geodesics are not defined over the :meth:`base_ring`. Also, that
        isometry might not be defined over the base ring but some other
        isometry is.

        In general, this is not a good approach. There might be a much better
        way to determine such an isometry explicitly.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A case where we can actually map things::

            sage: H._isometry_from_single_geodesics(H.geodesic(-1, 1), H.geodesic(0, 2))
            [ 1 -1]
            [ 0  1]

        In many cases, we fail to find the isometry::

            sage: g = H.geodesic(1, 2, 3, model="klein")
            sage: h = g.apply_isometry(matrix([[1, 2], [3, 4]]), on_right=True)
            sage: H._isometry_from_single_geodesics(g, h)
            Traceback (most recent call last):
            ...
            ValueError: ...

        """
        for isometry in self._isometry_from_primitives([(preimage, image), (preimage.midpoint(), image.midpoint())]):
            if isometry is None:
                import warnings
                warnings.warn("Could not determine an isometry of geodesics over the base ring. There might still be one but the implementation failed to detect it.")
            return isometry

    def _isometry_from_equations(self, conditions, filter):
        r"""
        Helper method for :meth:`isometry`.

        Return an isometry that satisfies ``conditions`` and ``filter``.

        INPUT:

        - ``conditions`` -- a function that receives a (symbolic) isometry and
          some (symbolic) variables and creates polynomial equations that a
          concrete isometry must satisfy.

        - ``filter`` -- a function that receives a concrete isometry and
          returns whether it maps objects correctly.

        ALGORITHM:

        We guess determine the entries of a 2×2 matrix with entries a, b, c, d
        by guessing for each term that it is non-zero and then building the
        symbolic relations that the isometry must satisfy by invoking
        ``conditions``. These symbolic conditions contain free linear variables
        coming from the fact that points are encoded projectively and geodesics
        in the dual. We tune these variables so that all entries of the matrix
        are in the base ring. (Namely, so that we can take all the square roots
        that show up.)

        The whole process is very ad-hoc and very slow since it computes lots
        of Gröbner bases. It is very likely that the approach is not
        mathematically sound but it worked for many random inputs that we
        presented it with.

        """
        from sage.all import PolynomialRing, matrix

        # Over the reals, the equations are different depending on whether we
        # send a geodesic to -another geodesic or +another geodesic. We try
        # both cases to see which one yields a system that we can solve over
        # the base ring.
        for sgn in [self.base_ring().one(), -self.base_ring().one()]:
            # We try to determine the matrix describing the isometry assuming
            # that "variable" is non-zero.
            for variable in ["a", "d", "b", "c"]:
                variables = ["a", "b", "c", "d", "λ1", "λ2"]
                variables.remove(variable)
                variables.append(variable)
                # We use a term order that guarantees that we will see an
                # equation for "variable" in the Gröbner basis.
                R = PolynomialRing(self.base_ring(), names=variables, order="degrevlex(5), lex(1)")
                # We are going to run the same procedure twice. Once with λ0 =
                # ±1 and then with a λ0 tuned to a value so that we can
                # actually solve for "variable".
                λ0 = sgn
                λ1 = R("λ1")
                λ2 = R("λ2")
                a = R("a")
                b = R("b")
                c = R("c")
                d = R("d")
                variable = R(variable)

                # We keep track of whether we made any assumptions here that
                # mean that we might be ignoring solutions in this run.
                equivalence = True

                # The isometry as a symbolic 2×2 matrix.
                isometry = matrix([
                    [a, b],
                    [c, d]
                ])

                isometry = gl2_to_sim12(isometry)

                # Build equations for the symbolic variables and make sure that
                # the resulting variety is zero-dimensional.
                equations = conditions(isometry, (λ0, λ1, λ2))

                for λ in [λ1, λ2]:
                    if all([equation.degree(λ) <= 0 for equation in equations]):
                        # Force the unused variable λ to be =0
                        equations.append(λ)

                # The system of euations typically has no rational points.
                # We analyze the Gröbner basis to tune λ0 so that we get
                # rational points.
                J = list(R.ideal(equations).groebner_basis())

                if J == [1]:
                    # The equations are contradictory.
                    assert equivalence
                    return None

                # We extract an equation for "variable" from the Gröbner basis.
                equation = J[-1]

                assert equation.variables() == (variable,), f"expected Gröbner basis algorithm to yield an equation for {variable} but found {equation} instead"
                equation = equation.polynomial(variable)
                equation = equation.map_coefficients(self.base_ring(), new_base_ring=self.base_ring())
                variable = equation.parent().gen()

                # The variable must be zero. We continue with another
                # (non-zero) variable since we meant to deduce the value of the
                # scaling factor λ0.
                if equation == variable:
                    continue

                # The equation allows the case that the variable is zero.
                # We ignore this possibility here. If this is needed, then we
                # will find out in the loop for another variable.
                while not equation.constant_coefficient():
                    equivalence = False
                    equation >>= 1

                from sage.all import QQbar
                equation = equation.change_ring(QQbar)

                # We try to arrange things so that "variable" becomes an element of the base ring.
                # We look at the minpoly of variable and tune λ0 so that this
                # becomes a square root of something that is a square.
                for root in equation.roots(multiplicities=False):
                    minpoly = root.minpoly()
                    if minpoly.degree() == 1:
                        pass
                    elif minpoly.exponents() == [0, 2]:
                        λ0 = -sgn * ~self.base_ring()(minpoly.constant_coefficient())
                    else:
                        # We cannot change the equations for this root to show
                        # up over the base ring.
                        continue

                    assert λ0 in self.base_ring() and λ0 != 0, f"did not deduce a non-zero constant for {λ0=} from {equation}"

                    # We could now patch the existing Gröbner basis and solve directly, but
                    # we just solve again for the correct value of λ0.
                    equations = conditions(isometry, (λ0, λ1, λ2))

                    for λ in [λ1, λ2]:
                        if all([equation.degree(λ) <= 0 for equation in equations]):
                            # Force the unused variable λ to be =0
                            equations.append(λ)

                    solutions = R.ideal(equations).variety()

                    assert solutions, f"After tuning the constant of the equations describing the isometry, there should be a solution but we did not find any."

                    solutions = [matrix([[solution[a], solution[b]], [solution[c], solution[d]]]) for solution in solutions]

                    # Check which solutions do not only satisfy "conditions"
                    # but map the underlying objects correctly, i.e., they also
                    # pass "filter".
                    solutions = [solution for solution in solutions if filter(solution)]

                    if not solutions:
                        continue

                    # Prefer an isometry of determinant 1 and isometries with 1 entries.
                    return max(solutions, key=lambda isometry: (isometry.det(), isometry[1][1] == 1, isometry[1][0] == 1))

        # No luck with this approach. We hope that this means that no such
        # isometry exists over the base ring but it's not entirely clear
        # whether that is actually true.
        return None

    def _repr_(self):
        r"""
        Return a printable representation of this hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: HyperbolicPlane(AA)
            Hyperbolic Plane over Algebraic Real Field

        """
        return f"Hyperbolic Plane over {repr(self.base_ring())}"


class HyperbolicGeometry:
    r"""
    Predicates and primitive geometric constructions over a base ``ring``.

    This class and its subclasses implement the core underlying hyperbolic
    geometry that depends on the base ring. For example, when deciding whether
    two points in the hyperbolic plane are equal, we cannot just compare their
    coordinates if the base ring is inexact. Therefore, that predicate is
    implemented in this "geometry" class and is implemented differently by
    :class:`HyperbolicExactGeometry` for exact and
    :class:`HyperbolicEpsilonGeometry` for inexact rings.

    INPUT:

    - ``ring`` -- a ring, the ring in which coordinates in the hyperbolic plane
      will be represented

    .. NOTE::

        Abstract methods are not marked with `@abstractmethod` since we cannot
        use the ABCMeta metaclass to enforce their implementation; otherwise,
        our subclasses could not use the unique representation metaclasses.

    EXAMPLES:

    The specific hyperbolic geometry implementation is picked automatically,
    depending on whether the base ring is exact or not::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()
        sage: H.geometry
        Exact geometry over Rational Field
        sage: H(0) == H(1/1024)
        False

    However, we can explicitly use a different or custom geometry::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicEpsilonGeometry
        sage: H = HyperbolicPlane(QQ, HyperbolicEpsilonGeometry(QQ, 1/1024))
        sage: H.geometry
        Epsilon geometry with ϵ=1/1024 over Rational Field
        sage: H(0) == H(1/2048)
        True

    .. SEEALSO::

        :class:`HyperbolicExactGeometry`, :class:`HyperbolicEpsilonGeometry`

    """

    def __init__(self, ring):
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicGeometry
            sage: H = HyperbolicPlane()
            sage: isinstance(H.geometry, HyperbolicGeometry)
            True

        """
        self._ring = ring

    def base_ring(self):
        r"""
        Return the ring over which this geometry is implemented.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry.base_ring()
            Rational Field

        """
        return self._ring

    def _zero(self, x):
        r"""
        Return whether ``x`` should be considered zero in the
        :meth:`base_ring`.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry. Also, this predicate lacks
            the context of other elements; a proper predicate should also take
            other elements into account to decide this question relative to the
            other values.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)
            sage: H.geometry._zero(1)
            False
            sage: H.geometry._zero(1e-9)
            True

        """
        return self._cmp(x, 0) == 0

    def _cmp(self, x, y):
        r"""
        Return how ``x`` compares to ``y``.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`

        - ``y`` -- an element of the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry._cmp(0, 0)
            0
            sage: H.geometry._cmp(0, 1)
            -1
            sage: H.geometry._cmp(1, 0)
            1

        ::

            sage: H = HyperbolicPlane(RR)
            sage: H.geometry._cmp(0, 0)
            0
            sage: H.geometry._cmp(0, 1)
            -1
            sage: H.geometry._cmp(1, 0)
            1
            sage: H.geometry._cmp(1e-10, 0)
            0

        """
        if self._equal(x, y):
            return 0
        if x < y:
            return -1

        assert x > y, "Geometry over this ring must override _cmp since not (x == y) and not (x < y) does not imply x > y"
        return 1

    def _sgn(self, x):
        r"""
        Return the sign of ``x``.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry. Also, this predicate lacks
            the context of other elements; a proper predicate should also take
            other elements into account to decide this question relative to the
            other values.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)
            sage: H.geometry._sgn(1)
            1
            sage: H.geometry._sgn(-1)
            -1
            sage: H.geometry._sgn(1e-10)
            0

        """
        return self._cmp(x, 0)

    def _equal(self, x, y):
        r"""
        Return whether ``x`` and ``y`` should be considered equal in the :meth:`base_ring`.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`

        - ``y`` -- an element of the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)
            sage: H.geometry._equal(0, 1)
            False
            sage: H.geometry._equal(0, 1e-10)
            True

        """
        raise NotImplementedError("this geometry does not implement _equal()")

    def _determinant(self, a, b, c, d):
        r"""
        Return the determinant of the 2×2 matrix ``[[a, b], [c, d]]`` or
        ``None`` if the matrix is singular.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        INPUT:

        - ``a`` -- an element of the :meth:`base_ring`

        - ``b`` -- an element of the :meth:`base_ring`

        - ``c`` -- an element of the :meth:`base_ring`

        - ``d`` -- an element of the :meth:`base_ring`

        EXAMPLES:

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry._determinant(1, 2, 3, 4)
            -2
            sage: H.geometry._determinant(0, 10^-10, 1, 1)
            -1/10000000000

        """
        det = a * d - b * c
        if self._zero(det):
            return None
        return det

    def change_ring(ring):
        r"""
        Return this geometry with the :meth:`base_ring` changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry
            Exact geometry over Rational Field
            sage: H.geometry.change_ring(AA)
            Exact geometry over Algebraic Real Field

        """
        raise NotImplementedError("this geometry does not implement change_ring()")

    def projective(self, p, q, point):
        r"""
        Return the ideal point with projective coordinates ``[p: q]`` in the
        upper half plane model.

        INPUT:

        - ``p`` -- an element of the :meth:`base_ring`

        - ``q`` -- an element of the :meth:`base_ring`

        - ``point`` -- the :meth:`HyperbolicPlane.point` to create points

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.projective(1, 0, H.point)
            ∞
            sage: H.geometry.projective(0, 1, H.point)
            0

        """
        if self._zero(p) and self._zero(q):
            raise ValueError("one of p and q must not be zero")

        if self._zero(q):
            return point(0, 1, model="klein", check=False)

        return point(p / q, 0, model="half_plane", check=False)

    def half_circle(self, center, radius_squared, geodesic):
        r"""
        Return the geodesic around the real ``center`` and with
        ``radius_squared`` in the upper half plane.

        INPUT:

        - ``center`` -- an element of the :meth:`base_ring`, the center of the
          half circle on the real axis

        - ``radius_squared`` -- a positive element of the :meth:`base_ring`,
          the square of the radius of the half circle

        - ``geodesic`` -- the :meth:`HyperbolicPlane.geodesic` to create geodesics

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.half_circle(0, 1, H.geodesic)
            {(x^2 + y^2) - 1 = 0}

        Unfortunately, this does not work correctly over inexact fields yet::

            sage: H = HyperbolicPlane(RR)
            sage: H.geometry.half_circle(0, 1e-32, H.geodesic)
            Traceback (most recent call last):
            ...
            ValueError: radius must be positive

        """
        if self._sgn(radius_squared) <= 0:
            raise ValueError("radius must be positive")

        # Represent this geodesic as a(x^2 + y^2) + b*x + c = 0
        a = 1
        b = -2 * center
        c = center * center - radius_squared

        return geodesic(a, b, c, model="half_plane")

    def vertical(self, real, geodesic):
        r"""
        Return the vertical geodesic at the ``real`` ideal point in the upper
        half plane model.

        INPUT:

        - ``real`` -- an element of the :meth:`base_ring`

        - ``geodesic`` -- the :meth:`HyperbolicPlane.geodesic` to create geodesics

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.vertical(0, H.geodesic)
            {-x = 0}

        Unfortunately, this does not allow creation of verticals at large reals
        over inexact fields yet::

            sage: H = HyperbolicPlane(RR)
            sage: H.geometry.vertical(1e32, H.geodesic)
            Traceback (most recent call last):
            ...
            ValueError: equation ... does not define a chord in the Klein model

        """
        # Convert the equation -x + real = 0 to the Klein model.
        return geodesic(real, -1, -real, model="klein")

    def classify_point(self, x, y, model):
        r"""
        Return whether the point ``(x, y)`` is finite, ideal, or ultra-ideal.

        INPUT:

        - ``x`` -- an element of the :meth:`base_ring`

        - ``y`` -- an element of the :meth:`base_ring`

        - ``model`` -- a supported model, either ``"half_plane"`` or
          ``"klein"``

        OUTPUT: ``1`` if the point is finite, ``0`` if the point is ideal,
        ``-1`` if the point is neither of the two.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.classify_point(0, 1, model="half_plane")
            1
            sage: H.geometry.classify_point(0, 0, model="half_plane")
            0
            sage: H.geometry.classify_point(0, -1, model="half_plane")
            -1

        Unfortunately, over an inexact field, this detects points close to the
        real axis as being ultra-ideal::

            sage: H = HyperbolicPlane(RR)
            sage: H.geometry.classify_point(0, -1e32, model="half_plane")
            -1

        """
        if model == "half_plane":
            return self._sgn(y)

        if model == "klein":
            raise NotImplementedError("cannot classify points in the Klein model yet")

        raise NotImplementedError("unsupported model")

    def intersection(self, f, g):
        r"""
        Return the point of intersection between the Euclidean lines ``f`` and ``g``.

        INPUT:

        - ``f`` -- a triple of elements ``(a, b, c)`` of :meth:`base_ring`
          encoding the line `a + bx + cy = 0`

        - ``g`` -- a triple of elements ``(a, b, c)`` of :meth:`base_ring`
          encoding the line `a + bx + cy = 0`

        OUTPUT: A pair of elements of :meth:`base_ring`, the coordinates of the
        point of intersection, or ``None`` if the lines do not intersect.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.intersection((0, 1, 0), (0, 0, 1))
            (0, 0)

        """
        (fa, fb, fc) = f
        (ga, gb, gc) = g
        det = self._determinant(fb, fc, gb, gc)

        if det is None:
            return None

        x = (-gc * fa + fc * ga) / det
        y = (gb * fa - fb * ga) / det

        return (x, y)


class HyperbolicExactGeometry(UniqueRepresentation, HyperbolicGeometry):
    r"""
    Predicates and primitive geometric constructions over an exact base ring.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()
        sage: H.geometry
        Exact geometry over Rational Field

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicExactGeometry
        sage: isinstance(H.geometry, HyperbolicExactGeometry)
        True

    .. SEEALSO::

        :class:`HyperbolicEpsilonGeometry` for an implementation over inexact rings

    """

    def _equal(self, x, y):
        r"""
        Return whether the numbers ``x`` and ``y`` should be considered equal
        in exact geometry.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry._equal(0, 1)
            False
            sage: H.geometry._equal(0, 1/2**64)
            False
            sage: H.geometry._equal(0, 0)
            True

        """
        return x == y

    def change_ring(self, ring):
        r"""
        Return this geometry with the :meth:`base_ring` changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry.change_ring(QQ) == H.geometry
            True
            sage: H.geometry.change_ring(AA)
            Exact geometry over Algebraic Real Field

        When presented with the reals, we guess the epsilon for the
        :class:`HyperbolicEpsilonGeometry` to be consistent with the
        :class:`HyperbolicGeometry` constructor. (And also, because we use this
        frequently when plotting.)::

            sage: H.geometry.change_ring(RR)
            Epsilon geometry with ϵ=1.00000000000000e-6 over Real Field with 53 bits of precision

        """
        from sage.all import RR
        if ring is RR:
            return HyperbolicEpsilonGeometry(ring, 1e-6)

        if not ring.is_exact():
            raise ValueError("cannot change_ring() to an inexact ring")

        return HyperbolicExactGeometry(ring)

    def __repr__(self):
        r"""
        Return a printable representation of this geometry.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry
            Exact geometry over Rational Field

        """
        return f"Exact geometry over {self._ring}"


class HyperbolicEpsilonGeometry(UniqueRepresentation, HyperbolicGeometry):
    r"""
    Predicates and primitive geometric constructions over a base ``ring`` with
    "precision" ``epsilon``.

    This is an alternative to :class:`HyperbolicExactGeometry` over inexact
    rings. The exact meaning of the ``epsilon`` parameter is a bit fuzzy, but
    the basic idea is that two numbers are considered equal in this geometry if
    their relative difference is less than ``epsilon``, see :meth:`_equal` for
    details.

    INPUT:

    - ``ring`` -- a ring, the ring in which coordinates in the hyperbolic plane
      will be represented

    - ``epsilon`` -- an error bound

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicEpsilonGeometry
        sage: H = HyperbolicPlane(RR, HyperbolicEpsilonGeometry(RR, 1/1024))

    The ``epsilon`` affects the notion of equality in this geometry::

        sage: H(0) == H(1/2048)
        True

        sage: H(1/2048) == H(2/2048)
        False

    This geometry is meant for inexact rings, however, it can also be used in
    exact rings::

        sage: H = HyperbolicPlane(QQ, HyperbolicEpsilonGeometry(QQ, 1/1024))

    .. SEEALSO::

        :class:`HyperbolicExactGeometry`

    """

    def __init__(self, ring, epsilon):
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicEpsilonGeometry
            sage: H = HyperbolicPlane(RR)
            sage: isinstance(H.geometry, HyperbolicEpsilonGeometry)
            True

        """
        super().__init__(ring)
        self._epsilon = ring(epsilon)

    def _equal(self, x, y):
        r"""
        Return whether ``x`` and ``y`` should be considered equal numbers with
        respect to an ε error.

        .. NOTE::

            This method has not been tested much. Since this underlies much of
            the inexact geometry, we should probably do something better here,
            see e.g., https://floating-point-gui.de/errors/comparison/

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)

            sage: H.geometry._equal(1, 2)
            False
            sage: H.geometry._equal(1, 1 + 1e-32)
            True
            sage: H.geometry._equal(1e-32, 1e-32 + 1e-33)
            False
            sage: H.geometry._equal(1e-32, 1e-32 + 1e-64)
            True

        """
        if x == 0 or y == 0:
            return abs(x - y) < self._epsilon

        return abs(x - y) <= (abs(x) + abs(y)) * self._epsilon

    def _determinant(self, a, b, c, d):
        r"""
        Return the determinant of the 2×2 matrix ``[[a, b], [c, d]]`` or
        ``None`` if the matrix is singular.

        INPUT:

        - ``a`` -- an element of the :meth:`base_ring`

        - ``b`` -- an element of the :meth:`base_ring`

        - ``c`` -- an element of the :meth:`base_ring`

        - ``d`` -- an element of the :meth:`base_ring`

        EXAMPLES:

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)

            sage: H.geometry._determinant(1, 2, 3, 4)
            -2
            sage: H.geometry._determinant(1e-10, 0, 0, 1e-10)
            1.00000000000000e-20

        Unfortunately, we are not implementing any actual rank detecting
        algorithm (QR decomposition or such) here. So, we do not detect that
        this matrik is singular::

            sage: H.geometry._determinant(1e-127, 1e-128, 1, 1)
            9.00000000000000e-128

        """
        det = a * d - b * c
        if det == 0:
            # Note that we should instead numerically detect the rank here.
            return None
        return det

    def projective(self, p, q, point):
        r"""
        Return the ideal point with projective coordinates ``[p: q]`` in the
        upper half plane model.

        INPUT:

        - ``p`` -- an element of the :meth:`base_ring`

        - ``q`` -- an element of the :meth:`base_ring`

        - ``point`` -- the :meth:`HyperbolicPlane.point` to create points

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)

        The point ``[p: q]`` is the point at infinity if ``q`` is very small in
        comparison to ``p``::

            sage: H.geometry.projective(1, 0, H.point)
            ∞

            sage: H.geometry.projective(1e-8, 1e-16, H.point)
            ∞

            sage: H.geometry.projective(1e-8, -1e-16, H.point)
            ∞

        Even though ``q`` might be small, ``[p: q]`` is not the point at
        infinity if both coordinates are of similar size::

            sage: H.geometry.projective(1e-16, 1e-16, H.point)
            1.00000000000000

            sage: H.geometry.projective(-1e-16, 1e-16, H.point)
            -1.00000000000000

        """
        if self._zero(p) and self._zero(q):
            try:
                pq = p / q
            except ZeroDivisionError:
                return point(0, 1, model="klein", check=False)

            try:
                qp = q / p
            except ZeroDivisionError:
                return point(0, 0, model="half_plane", check=False)

            if self._zero(qp):
                return point(0, 1, model="klein", check=False)

            return point(pq, 0, model="half_plane", check=False)

        return super().projective(p, q, point)

    def change_ring(self, ring):
        r"""
        Return this geometry over ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)

            sage: H.geometry.change_ring(RR) is H.geometry
            True

            sage: H.geometry.change_ring(RDF)
            Epsilon geometry with ϵ=1e-06 over Real Double Field
 
        """
        if ring.is_exact():
            raise ValueError("cannot change_ring() to an exact ring")

        return HyperbolicEpsilonGeometry(ring, self._epsilon)

    def __repr__(self):
        r"""
        Return a printable representation of this geometry.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)
            sage: H.geometry
            Epsilon geometry with ϵ=1.00000000000000e-6 over Real Field with 53 bits of precision

        """
        return f"Epsilon geometry with ϵ={self._epsilon} over {self._ring}"


class HyperbolicConvexSet(Element):
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
        r"""
        Return a minimal set of half spaces whose intersection is this convex set.

        The half spaces are ordered by :meth:`HyperbolicHalfSpace._less_than`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).left_half_space().half_spaces()
            {{x ≤ 0},}

            sage: H.vertical(0).half_spaces()
            {{x ≤ 0}, {x ≥ 0}}

            sage: H(0).half_spaces()
            {{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}}

        """
        raise NotImplementedError(f"{type(self)} does not implement half_spaces()")

    def _test_half_spaces(self, **options):
        r"""
        Verify that this convex set implements :meth:`half_spaces` correctly.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_half_spaces()

        """
        tester = self._tester(**options)

        half_spaces = self.half_spaces()

        tester.assertEqual(self.parent().intersection(*half_spaces), self.unoriented())

        tester.assertTrue(isinstance(half_spaces, HyperbolicHalfSpaces))

        for a, b in zip(list(half_spaces), list(half_spaces)[1:]):
            tester.assertTrue(HyperbolicHalfSpaces._lt_(a, b))

    def _check(self, require_normalized=True):
        r"""
        Validate this convex subset.

        Subclasses run specific checks here that can be disabled when creating
        objects with ``check=False``.

        If ``require_normalized``, we also check that the object has the
        correct implementation class, e.g., that a point is a
        :class:`HyperbolicPoint` and not say a
        :class:`HyperbolicOrientedSegment` of length zero.

        EXAMPLES:

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: P = H.point(0, 0, model="klein")
            sage: P._check()
            sage: P = H.point(1, 1, model="klein", check=False)
            sage: P._check()
            Traceback (most recent call last):
            ...
            ValueError: point (1, 1) is not in the unit disk in the Klein model

        """
        pass

    def _normalize(self):
        r"""
        Return this set possibly rewritten in a simpler form.

        This method is only relevant for sets created with ``check=False``.
        Such sets might have been created in a non-canonical way, e.g., when
        creating a :class:`HyperbolicOrientedSegment` whose start and end point are ideal,
        then this is actually a geodesic and it shuold be described as such.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: segment = H.segment(H.vertical(-1), start=H.infinity(), end=H.infinity(), check=False, assume_normalized=True)
            sage: segment
            {-x - 1 = 0} ∩ {x - 1 ≥ 0} ∩ {x - 1 ≤ 0}
            sage: segment._normalize()
            ∞

        """
        return self

    def _test_normalize(self, **options):
        r"""
        Verify that normalization is idempotent.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: segment = H.segment(H.vertical(-1), start=H.infinity(), end=H.infinity(), check=False, assume_normalized=True)
            sage: segment._test_normalize()

        """
        tester = self._tester(**options)

        normalization = self._normalize()

        tester.assertEqual(normalization, normalization._normalize())

    def unoriented(self):
        r"""
        Return the non-oriented version of this set.

        Some sets such as geodesics and segments can have an explicit
        orientation. This method returns the underlying set without any
        explicit orientation.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).unoriented()
            {x = 0}

        """
        return self.change(oriented=False)

    def _test_unoriented(self, **options):
        r"""
        Verify that :meth:`unoriented` is implemented correctly.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_unoriented()

        """
        tester = self._tester(**options)

        tester.assertEqual(self.unoriented(), self.unoriented().unoriented())

    def intersection(self, other):
        r"""
        Return the intersection with the ``other`` convex set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).intersection(H.vertical(1))
            ∞

        ..SEEALSO::

            :meth:`HyperbolicPlane.intersection`

        """
        return self.parent().intersection(self, other)

    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(I) in H.empty_set()
            False

            sage: I in H.vertical(0)
            True

            sage: 2*I in H.half_circle(0, 1).left_half_space()
            True

            sage: I/2 in H.half_circle(0, 1).left_half_space()
            False

        .. NOTE::

            There is currently no way to check whether a point is in the
            interior of a set.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.is_subset` to check containment of
            arbitrary sets.

        """
        for half_space in self.half_spaces():
            if point not in half_space:
                return False

        return True

    def _test_contains(self, **options):
        r"""
        Verify that :meth:`__contains__` is implemented correctly.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0)._test_contains()

        """
        tester = self._tester(**options)

        try:
            vertices = self.vertices()
        except ValueError:
            # The vertex coordinates might not be representable in the base ring.
            return

        for vertex in vertices:
            tester.assertIn(vertex, self)

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices bounding this hyperbolic set.

        This returns both finite and ideal vertices.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``) whether to
          include marked vertices that are not actual cornerns of the convex
          set.

        OUTPUT:

        A set of points, namely :class:`HyperbolicVertices`. Iteration of this
        set happens incounterclockwise order (as seen from the inside of the
        convex set.)

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The empty set has no vertices::

            sage: H.empty_set().vertices()
            {}

        A point is its own vertex::

            sage: H(0).vertices()
            {0,}
            sage: H(I).vertices()
            {I,}
            sage: H(oo).vertices()
            {∞,}

        The vertices of a geodesic are its ideal end points::

            sage: H.vertical(0).vertices()
            {0, ∞}

        The vertices of a half space are the ideal end points of its boundary geodesic::

            sage: H.vertical(0).left_half_space().vertices()
            {0, ∞}

        The vertices a polygon can be finite and ideal::

            sage: P = H.polygon([H.vertical(0).left_half_space(), H.half_circle(0, 1).left_half_space()])
            sage: P.vertices()
            {-1, I, ∞}

        If a polygon has marked vertices they are included::

            sage: P = H.polygon([H.vertical(0).left_half_space(), H.half_circle(0, 1).left_half_space()], marked_vertices=[2*I])
            sage: P.vertices()
            {-1, I, 2*I, ∞}

            sage: P.vertices(marked_vertices=False)
            {-1, I, ∞}

        """
        return self.parent().polygon(self.half_spaces(), check=False, assume_sorted=True, assume_minimal=True).vertices()

    def is_finite(self):
        r"""
        Return whether all points in this set are finite.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().is_finite()
            True

            sage: H.vertical(0).is_finite()
            False

            sage: H.vertical(0).left_half_space().is_finite()
            False

            sage: H(I).segment(2*I).is_finite()
            True

            sage: H(0).segment(I).is_finite()
            False

            sage: P = H.polygon([
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 1).left_half_space(),
            ....:     H.half_circle(0, 4).right_half_space(),
            ....: ])
            sage: P.is_finite()
            False

            sage: P = H.polygon([
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space(),
            ....:     H.half_circle(0, 4).right_half_space(),
            ....: ])
            sage: P.is_finite()
            True

        """
        return all([vertex.is_finite() for vertex in self.vertices()])

    def change_ring(self, ring):
        r"""
        Return this set as an element of the hyperbolic plane over ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H(0)
            sage: q = p.change_ring(AA)
            sage: q.parent().base_ring()
            Algebraic Real Field

        Changing the base ring can provide coordinates for points::

            sage: p = H.half_circle(0, 2).start()
            sage: p.coordinates()
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: q = p.change_ring(AA)
            sage: q.coordinates()
            (-1.414213562373095?, 0)

        Note that changing the ring only works in relatively trivial cases::

            sage: q = HyperbolicPlane(AA).point(sqrt(2), 0, model="half_plane")

            sage: p = q.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: ...

        Most other sets also support changing the base ring::

            sage: g = H.half_circle(0, 2)
            sage: g.start().coordinates()
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: g.change_ring(AA).start().coordinates()
            (-1.414213562373095?, 0)

        .. SEEALSO::

            :meth:`change` for a more general interface to changing properties
            of hyperbolic sets.

            :meth:`HyperbolicPlane.change_ring` for the hyperbolic plane that
            the resulting objects lives in.

        """
        return self.change(ring=ring)

    def _test_change_ring(self, **options):
        r"""
        Verify that this set implements :meth:`change_ring`.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_change_ring()

        """
        tester = self._tester(**options)
        tester.assertEqual(self, self.change_ring(self.parent().base_ring()))

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return a modified copy of this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        r"""
        Verify that the full interface of :meth:`change` has been implemented.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        r"""
        Return a plot of this subset.

        Consult the implementation in the subclasses for a list supported
        keyword arguments, in particular :meth:`HyperbolicConvexPolygon.plot`.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"``

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).plot()
            ...Graphics object consisting of 1 graphics primitive

        """
        raise NotImplementedError(f"this {type(self)} does not support plotting")

    def _test_plot(self, **options):
        r"""
        Verify that this set implements :meth:`plot`.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?

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
                # TODO: Use a specialized predicate instead of the _method.
                if (i == j) == self.parent().geometry._zero(entry):
                    raise ValueError("invalid isometry")

        # TODO: Use a specialized predicate instead of the _method.
        if not self.parent().geometry._equal(D[0, 0], D[1, 1]):
            raise ValueError("invalid isometry")

        # TODO: Use a specialized predicate instead of the _method.
        if not self.parent().geometry._equal(
            D[0, 0], -D[2, 2]
        ):
            raise ValueError("invalid isometry")

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: This can be implemented generically.
        raise NotImplementedError

    def apply_isometry(self, isometry, model="half_plane", on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Explain what the isometry actually encodes, in particular if its determinant is negative.
        r"""
        Return the image of this set under the ``isometry``.

        INPUT:

        - ``isometry`` -- a 2×2 matrix in `GL(2,\mathbb{R})`, if ``model`` is
          ``"half_plane"``, or a 3×3 matrix giving a similitude that preserves
          a quadratic form of type `(1, 2)`, if ``model`` is ``"klein"``.

        - ``model`` -- a string (default: ``"half_plane"``); either
          ``"half_plane"`` or ``"klein"``

        - ``on_right`` -- a boolean (default: ``False``) whether to apply the
          right action.

        ALGORITHM:

        If ``model`` is ``"half_plane"``, the 2×2 matrix with entries `a, b, c,
        d` encodes a fractional linear transformation sending

        .. MATH::

            z \mapsto \frac{az + b}{cz + d}

        if the determinant is positive, and

        .. MATH::

            z \mapsto \frac{a\bar{z} + b}{a\bar{z} + d}

        if the determinant is negative. Note that these maps are invariant
        under scaling the matrix with a non-zero real.

        In any case, we convert the matrix to a corresponding 3×3 matrix, see
        :func:`gl2_to_sim12` and apply the isometry in the Klein model.

        TODO: Explain how the isometries in the Klein model work.

        REFERENCES:

        - Svetlana Katok, "Fuchsian Groups", Chicago University Press, Section
          1.3; for the isometries of the upper half plane.

        """
        if model == "half_plane":
            isometry = gl2_to_sim12(isometry)
            model = "klein"

        if model == "klein":
            self._check_isometry_klein(isometry)
            return self._apply_isometry_klein(isometry, on_right=on_right)

        raise NotImplementedError(
            "applying isometry not supported in this hyperbolic model"
        )

    def _acted_upon_(self, x, self_on_left):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return the result of acting upon this set with ``x``.

        EXAMPLES:

        The Möbius transformation that sends `z` to `(1 + 2z)/(3 + 4z)`::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: p = H(I)
            sage: m = matrix([[1, 2], [3, 4]])
            sage: m * p
            11/25 + 2/25*I
            sage: p * m
            -7/5 + 1/5*I

        TESTS::

            sage: m0 = matrix(2, [1, 2, 3, 4])
            sage: m1 = matrix(2, [1, 1, 0, 1])
            sage: p = HyperbolicPlane()(I + 1)
            sage: assert (m0 * m1) * p == m0 * (m1 * p)
            sage: assert p * (m0 * m1) == (p * m0) * m1
        """
        return self.apply_isometry(x, on_right=self_on_left)

    def is_subset(self, other):
        r"""
        Return whether this set is a subset of ``other``.

        INPUT:

        - ``other`` -- another hyperbolic convex set

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
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

        for vertex in self.vertices():
            if vertex not in other:
                return False

        if self.dimension() <= 1:
            return True

        # TODO: This is only correct for normalized sets?
        return self.intersection(other) == self

    # TODO: Test that is_subset can compare all kinds of sets by inclusion.

    def an_element(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Test that everything implements an_element().
        return next(iter(self.vertices()))

    @classmethod
    def _enhance_plot(self, plot, model):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if model == "klein":
            from sage.all import circle

            plot = circle([0, 0], 1, fill=False, color="#d1d1d1", zorder=-1) + plot

        return plot

    def is_empty(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self.dimension() < 0

    def __bool__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return not self.is_empty()

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Test that this is a ZZ integer.
        raise NotImplementedError(f"{type(self)} does not implement dimension() yet")

    def is_point(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self.dimension() == 0

    def is_oriented(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Test that oriented sets implement _neg_ correctly.
        return isinstance(self, HyperbolicOrientedConvexSet)

    def edges(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Implement by iterating vertices.
        # TODO: Test that other implementation are consistent with vertices()
        raise NotImplementedError

    def __hash__(self):
        r"""
        Return a hash value for this convex set.

        Specific sets should override this method.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: hash(H.empty_set())
            0

        Note that has values of sets over different base rings might not be
        consistent::

            sage: HyperbolicPlane(ZZ).half_circle(0, 2) == HyperbolicPlane(AA).half_circle(0, 2)
            True

            sage: hash(HyperbolicPlane(ZZ).half_circle(0, 2)) == hash(HyperbolicPlane(AA).half_circle(0, 2))
            False

        Sets over inexact base rings are not be hashable (since their hash
        would not be compatible with the notion of equality)::

            sage: hash(HyperbolicPlane(RR).vertical(0))
            Traceback (most recent call last):
            ...
            TypeError: cannot hash geodesic defined over inexact base ring

        """
        raise NotImplementedError(f"the set {self} is not hashable")

    def _test_hash(self, **options):
        r"""
        Verify that this set implements a good hash function.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_plot()

        Nothing is tested for unhashable sets::

            sage: H = HyperbolicPlane(RR)

            sage: H.an_element()._test_hash()

        """
        tester = self._tester(**options)

        # We refuse to hash inexact elements.
        if not self.parent().base_ring().is_exact():
            if not self.is_empty():
                with tester.assertRaises(TypeError):
                    hash(self)
            return

        tester.assertEqual(hash(self), hash(self))

        # We test that the hash function is good enough to distinguish this set
        # from some generic sets. While, there will be hash collisions
        # eventually, they should not show up with such non-random sets for a
        # good hash function.
        for subset in self.parent().some_elements():
            if subset != self:
                tester.assertNotEqual(hash(self), hash(subset))

    def _isometry_conditions(self, other):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        raise NotImplementedError

    def _isometry_equations(self, isometry, other, λ):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        raise NotImplementedError


class HyperbolicOrientedConvexSet(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        raise NotImplementedError


class HyperbolicHalfSpace(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    A closed half space of the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.half_circle(0, 1).left_half_space()
        {(x^2 + y^2) - 1 ≥ 0}

    """

    def __init__(self, parent, geodesic):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("geodesic must be an oriented geodesic")

        self._geodesic = geodesic

    def equation(self, model, normalization=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return an inequality for this half space as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is in the half space if it satisfies `a(x^2 + y^2) + bx + c \ge 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy \ge 0`.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.
        """
        return self._geodesic.equation(model=model, normalization=normalization)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return a printable representation of this half space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.half_circle(0, 1).right_half_space()

            sage: S
            {(x^2 + y^2) - 1 ≤ 0}

            sage: -S
            {(x^2 + y^2) - 1 ≥ 0}

        """
        # TODO: Use geodesic printing instead.
        # TODO: Use a specialized predicate instead of the _method.
        sgn = self.parent().geometry._sgn

        # Convert to the upper half plane model as a(x^2 + y^2) + bx + c ≥ 0.
        a, b, c = self.equation(model="half_plane", normalization=["gcd", None])

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
        # TODO: Benchmark?
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        r"""
        Return the closure of the complement of this half space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        r"""
        Return a geodesic on the boundary of this half space, oriented such that the half space is on its left.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.vertical(0).left_half_space()
            sage: S.boundary()
            {-x = 0}

        """
        return self._geodesic

    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this half space.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: h = H.vertical(0).left_half_space()
            sage: I in h
            True

            sage: I - 1 in h
            True

            sage: I + 1 in h
            False

            sage: oo in h
            True

        .. NOTE::

            The implementation is currently not very robust over inexact rings.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.is_subset` to check containment of
            arbitrary sets.

        """
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        # We should use a specialized predicate here to do something more
        # reasonable for points that are close to the boundary over inexact
        # rings.
        return self.parent().geometry._sgn(a + b * x + c * y) >= 0

    def _richcmp_(self, other, op):
        r"""
        Return how this half space compares to ``other`` with respect to the
        ``op`` operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether the two spaces are indistinguishable.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).left_half_space() == H.vertical(0).left_half_space()
            True

            sage: H.vertical(0).left_half_space() != H.vertical(0).right_half_space()
            True

            sage: H.vertical(0).left_half_space() != H.vertical(0)
            True

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicHalfSpace):
                return False
            return self._geodesic._richcmp_(other._geodesic, op)

        return super()._richcmp_(other, op)

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this half space in the hyperbolic ``model``.

        See :meth:`HyperbolicConvexPolygon.plot` for the supported keyword
        arguments.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: G = H.vertical(0).left_half_space().plot()

        In the half plane model, the half space is rendered as an infinite polygon::

            sage: G = H.vertical(0).left_half_space().plot()
            sage: G[0]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
                CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
                CartesianPathPlotCommand(code='RAYTO', args=(-1, 0)),
                CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])

        """
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
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        from sage.all import ZZ

        return ZZ(2)

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices bounding this half space.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          half space cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The vertices of a half space are the ideal end points of its boundary geodesic::

            sage: H.vertical(0).left_half_space().vertices()
            {0, ∞}

        Note that iteration in the set is not consistent with the orientation
        of the half space (it is chosen such that the subset relation on vertices
        can be checked quickly)::

            sage: h = H.vertical(0).left_half_space()
            sage: list(h.vertices())
            [0, ∞]

            sage: list((-h).vertices())
            [0, ∞]

        Use :meth:`HyperbolicGeodesic.start` and
        :meth:`HyperbolicGeodesic.end` on the :meth:`boundary` to get the end
        points in an order consistent with the orientation::

            sage: g = h.boundary()
            sage: g.start(), g.end()
            (0, ∞)

            sage: g = (-h).boundary()
            sage: g.start(), g.end()
            (∞, 0)

        Currently, vertices cannot be computed if some of them have coordinates
        which do not live over the :meth:`HyperbolicPlane.base_ring`; see
        :class:`HyperbolicVertices`::

            sage: H.half_circle(0, 2).left_half_space().vertices()
            Traceback (most recent call last):
            ...
            ValueError: ...

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return self.boundary().vertices()

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: This can be implemented generically.
        return self._geodesic.apply_isometry(isometry, model="klein", on_right=on_right).left_half_space()

    def __hash__(self):
        r"""
        Return a hash value for this half space

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since half spaces are hashable, they can be put in a hash table, such
        as a Python ``set``::

            sage: S = {H.vertical(0).left_half_space(), H.vertical(0).right_half_space()}
            sage: len(S)
            2

        """
        # Add the type to the hash value to distinguish the hash value from an
        # actual geodesic.
        return hash((type(self), self._geodesic))

    def _isometry_conditions(self, other):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        yield [(self.boundary(), other.boundary())]


class HyperbolicGeodesic(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
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
        # TODO: Benchmark?
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        # TODO: Make sure that all sets have is_finite, is_ideal, is_ultra_ideal.
        # TODO: Check that this also catches geodesics that touch the Klein disk.
        # TODO: Use a specialized predicate instead of the _method.
        return self.parent().geometry._cmp(self._b * self._b + self._c * self._c, self._a * self._a) <= 0

    def _repr_(self, model=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if model is None:
            model = "klein" if self.is_ultra_ideal() else "half_plane"

        if model == "half_plane":
            # Convert to the upper half plane model as a(x^2 + y^2) + bx + c = 0.
            a, b, c = self.equation(model="half_plane", normalization=["gcd", None])

            from sage.all import PolynomialRing

            R = PolynomialRing(self.parent().base_ring(), names="x")
            # TODO: Use a specialized method instead of the _method.
            if self.parent().geometry._sgn(a) != 0:
                return f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} = 0}}"
            else:
                return f"{{{repr(R([c, b]))} = 0}}"

        if model == "klein":
            a, b, c = self.equation(model="klein", normalization=["gcd", None])

            from sage.all import PolynomialRing

            R = PolynomialRing(self.parent().base_ring(), names=["x", "y"])
            polynomial_part = R({(1, 0): b, (0, 1): c})
            # TODO: Use a specialized predicate instead of the _method.
            if self.parent().geometry._sgn(a) != 0:
                return f"{{{repr(a)} + {repr(polynomial_part)} = 0}}"
            else:
                return f"{{{repr(polynomial_part)} = 0}}"

        raise NotImplementedError("printing not supported in this model")

    def equation(self, model, normalization=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Explain limitations when ultra ideal (and model=half_plane)
        r"""
        Return an equation for this geodesic as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is on the geodesic if it satisfies `a(x^2 + y^2) + bx + c = 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy = 0`.

        INPUT:

        - ``model`` -- the model in which this equation holds, either
          ``"half_plane"`` or ``"klein"``

        - ``normalization`` -- how to normalize the coefficients; the default
          ``None`` is not to normalize at all. Other options are ``gcd``, to
          divide the coefficients by their greatest common divisor, ``one``, to
          normalize the first non-zero coefficient to ±1. This can also be a
          list of such values which are then tried in order and exceptions are
          silently ignored unless they happen at the last option.

        If this geodesic :meth;`is_oriented`, then the sign of the coefficients
        is chosen to encode the orientation of this geodesic. The sign is such
        that the half plane obtained by replacing ``=`` with ``≥`` in above
        equationsis on the left of the geodesic.

        Note that the output might not uniquely describe the geodesic since the
        coefficients are only unique up to scaling.

        """
        normalization = normalization or [None]

        if isinstance(normalization, str):
            normalization = [normalization]

        from collections.abc import Sequence
        if not isinstance(normalization, Sequence):
            normalization = [normalization]

        normalization = list(normalization)
        normalization.reverse()

        a, b, c = self._a, self._b, self._c

        if model == "klein":
            a, b, c = a, b, c
        elif model == "half_plane":
            a, b, c = a + c, 2 * b, a - c
        else:
            raise NotImplementedError("cannot determine equation for this model yet")

        # TODO: Use a specialized predicate instead of the _method.
        sgn = self.parent().geometry._sgn
        sgn = -1 if (
                sgn(a) < 0
                or (sgn(a) == 0 and b < 0)
                or (sgn(a) == 0 and sgn(b) == 0 and sgn(c) < 0)
            ) else 1

        while normalization:
            strategy = normalization.pop()

            try:
                a, b, c = HyperbolicGeodesic._normalize_coefficients(a, b, c, sgn, strategy=strategy)
                break
            except Exception:
                if not normalization:
                    raise

        if not self.is_oriented():
            a *= sgn
            b *= sgn
            c *= sgn

        return a, b, c

    @classmethod
    def _normalize_coefficients(cls, a, b, c, sgn, strategy):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if strategy is None:
            return a, b, c

        if strategy == "gcd":
            def gcd(*coefficients):
                coefficients = [c for c in coefficients if c]
                if len(coefficients) == 1:
                    return sgn * coefficients[0]

                from sage.all import gcd
                return gcd(coefficients)

            d = gcd(a, b, c)
            assert d > 0
            return a / d, b / d, c / d

        if strategy == "one":
            if a:
                d = sgn * a
            elif b:
                d = sgn * b
            else:
                assert c
                d = sgn * c

            return a / d, b / d, c / d

        raise ValueError(f"unknown normalization {strategy}")

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).half_spaces()
            {{x ≤ 0}, {x ≥ 0}}

        """
        self = self.change(oriented=True)
        return HyperbolicHalfSpaces([self.left_half_space(), self.right_half_space()])

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this geodesic in the hyperbolic ``model``.

        See :meth:`HyperbolicSegment.plot` for the supported keyword arguments.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).plot()
            Graphics object consisting of 1 graphics primitive

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

    def pole(self):
        r"""
        Return the pole of this geodesic.

        ALGORITHM:

        The pole is the intersection of tangents of the Klein disk at
        the ideal endpoints of this geodesic, see `Wikipedia
        <https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model#Compass_and_straightedge_constructions>`.

        EXAMPLES:

        The pole of a geodesic is an ultra ideal point::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: p = H.vertical(2).pole(); p
            (1/2, 1)
            sage: p.is_ultra_ideal()
            True

        Computing the pole is only implemented if it is a finite point in the
        Euclidean plane::

            sage: H.half_circle(0, 1).pole()
            Traceback (most recent call last):
            ...
            NotImplementedError: can only compute pole if geodesic is a not a diameter in the Klein model

        The pole might not be defined without passing to a larger base ring::

            sage: H.half_circle(2, 2).pole()
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not in Rational Field

        """
        if self.is_diameter():
            raise NotImplementedError("can only compute pole if geodesic is a not a diameter in the Klein model")

        def tangent(endpoint):
            x, y = endpoint.coordinates(model="klein")
            return self.parent().geodesic(-(x*x + y*y), x, y, model="klein", check=False)

        A, B = self.vertices()

        pole = tangent(A)._intersection(tangent(B))

        assert pole is not None, "non-parallel lines must intersect"

        return pole

    def perpendicular(self, point_or_geodesic=None):
        r"""
        Return a geodesic that is perpendicular to this geodesic.

        If ``point_or_geodesic`` is a point, return a geodesic
        through that point.

        If ``point_or_geodesic`` is another geodesic, return a
        geodesic that is also perpendicular to that geodesic.

        ALGORITHM:

        We use the construction as explained on `Wikipedia
        <https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model#Compass_and_straightedge_constructions>`.

        INPUT:

        - ``point_or_geodesic`` -- a point or a geodesic in the
          hyperbolic plane or ``None`` (the default)

        EXAMPLES:

        Without parameters this method returns one of the many
        perpendicular geodesics::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: v = H.vertical(2)
            sage: v.perpendicular()
            {(x^2 + y^2) - 4*x - 1 = 0}

        We can request a perpendicular geodesic through a specific
        point::

            sage: v.perpendicular(2 + I)
            {(x^2 + y^2) - 4*x + 3 = 0}
            sage: v.perpendicular(I)
            {(x^2 + y^2) - 4*x - 1 = 0}

        In some cases, such a geodesic might not exist::

            sage: v.perpendicular(oo)
            Traceback (most recent call last):
            ...
            ValueError: ... does not define a chord in the Klein model

        We can request a geodesic that is also perpendicular to
        another geodesic::

            sage: v.perpendicular(H.half_circle(4, 1))
            {(x^2 + y^2) - 4*x + 1 = 0}

        In some cases, such a geodesic might not exist::

            sage: v.perpendicular(H.half_circle(2, 1))
            Traceback (most recent call last):
            ...
            ValueError: ... does not define a chord in the Klein model

        TESTS:

        Verify that this also works for geodesic that have no finite pole::

            sage: H.vertical(0).perpendicular()
            {(x^2 + y^2) - 1 = 0}
            sage: H.half_circle(0, 1).perpendicular()
            {x = 0}
            sage: H.half_circle(0, 1).perpendicular(0)
            {x = 0}
            sage: H.half_circle(0, 1).perpendicular(2)
            {2*(x^2 + y^2) - 5*x + 2 = 0}

            sage: H.vertical(0).perpendicular(H.half_circle(0, 1))
            Traceback (most recent call last):
            ...
            ValueError: no geodesic perpendicular to both {-x = 0} and {(x^2 + y^2) - 1 = 0}
            sage: H.half_circle(0, 1).perpendicular(H.vertical(0))
            Traceback (most recent call last):
            ...
            ValueError: no geodesic perpendicular to both {(x^2 + y^2) - 1 = 0} and {-x = 0}
            sage: H.vertical(0).perpendicular(H.vertical(0))
            {(x^2 + y^2) - 1 = 0}
            sage: H.vertical(0).perpendicular(H.half_circle(0, 1))
            Traceback (most recent call last):
            ...
            ValueError: no geodesic perpendicular to both {-x = 0} and {(x^2 + y^2) - 1 = 0}

        """
        # TODO: Orientation should be such that it is turning ccw from this geodesic.
        if point_or_geodesic is None:
            point_or_geodesic = self.an_element()

        point_or_geodesic = self.parent()(point_or_geodesic)

        if isinstance(point_or_geodesic, HyperbolicGeodesic):
            other = point_or_geodesic
            if self.unoriented() == other.unoriented():
                return self.perpendicular(self.an_element())

            if self.is_diameter() and other.is_diameter():
                raise ValueError(f"no geodesic perpendicular to both {self} and {other}")

            if other.is_diameter():
                return other.perpendicular(self)

            if self.is_diameter():
                # Construct the line a + bx + cy = 0 perpendicular to the
                # diameter through the pole of the other geodesic.
                b, c = (-self._c, self._b)
                x, y = other.pole().coordinates(model="klein")
                a = -(b*x + c*y)

                # The line might be not intersect the Klein disk. An error is
                # raised here in that case.
                return self.parent().geodesic(a, b, c, model="klein", oriented=False)

            # In the generic case, the perpendicular goes through both poles.
            # Throws an error if that line does not define a geodesic because
            # it's outside of the Klein disk.
            return self.parent().geodesic(self.pole(), other.pole(), oriented=False)
        else:
            point = point_or_geodesic
            if self.is_diameter():
                # Construct the line a + bx + cy = 0 perpendicular to the
                # diameter through the given point.
                b, c = (-self._c, self._b)
                x, y = point.coordinates(model="klein")
                a = -(b*x + c*y)

                perpendicular = self.parent().geodesic(a, b, c, model="klein", oriented=False)
            else:
                perpendicular = self.parent().geodesic(self.pole(), point, oriented=False)

            assert point in perpendicular

            return perpendicular

    def midpoint(self):
        r"""
        Return the fixed point of the (determinant one) Möbius transformation
        that interchanges the ideal endpoints of this geodesic.

        ALGORITHM:

        For the vertical connecting zero and infinity, the Möbius
        transformation sending z to `-1/z` has the imaginary unit as its fixed
        point. For a half circle centered at the origin its point on the
        imaginary axis must be the fixed point (due to symmetry or a direct
        computation.) All other geodesics, are just translated versions of
        these so we can just conjugate with a translation to determine the
        fixed point, i.e., the fixed point is a translate of one of the above.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).midpoint()
            I
            sage: H.vertical(1).midpoint()
            1 + I
            sage: H.half_circle(1, 1).midpoint()
            1 + I

        .. SEEALSO::

            :meth:`HyperbolicSegment.midpoint`

        """
        a, b = self.vertices()

        if self.is_vertical():
            if a == self.parent().infinity():
                a, b = b, a
            x, _ = a.coordinates(model="half_plane")

            return self.parent().point(x, 1, model="half_plane")

        ax, _ = a.coordinates(model="half_plane")
        bx, _ = b.coordinates(model="half_plane")

        m = (ax + bx) / 2
        r = abs(bx - ax) / 2

        return self.parent().point(m, r, model="half_plane")

    def is_diameter(self):
        r"""
        Return whether this geodesic is a diameter in the Klein model.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).is_diameter()
            True
            sage: H.vertical(1).is_diameter()
            False

        """
        return self.parent().point(0, 0, model="klein") in self

    def is_vertical(self):
        r"""
        Return whether this geodesic is a vertical in the upper half plane.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).is_vertical()
            True
            sage: H.half_circle(0, 1).is_vertical()
            False

        """
        # TODO: This check is probably done a few times in the code. We should use this predicate instead.
        return self.parent().infinity() in self

    def _richcmp_(self, other, op):
        r"""
        Return how this geodesic compares to ``other`` with respect to the
        ``op`` operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether the two geodesics are indistinguishable up to scaling
        of their defining equations.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0) == H.vertical(0)
            True

        We distinguish oriented and unoriented geodesics::

            sage: H.vertical(0).unoriented() == H.vertical(0)
            False

        We distinguish differently oriented geodesics::

            sage: H.vertical(0) == -H.vertical(0)
            False

        We do, however, identify geodesics whose defining equations differ by some scaling::

            sage: g = H.vertical(0)
            sage: g.equation(model="half_plane")
            (0, -2, 0)
            sage: h = H.geodesic(0, -4, 0, model="half_plane")
            sage: g.equation(model="half_plane") == h.equation(model="half_plane")
            False
            sage: g == h
            True

        .. NOTE::

            Over inexact rings, this method is not very reliable. To some
            extent this is inherent to the problem but also the implementation
            uses generic predicates instead of relying on a specialized
            implementation in the :class:`HyperbolicGeometry`.

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            # See note in the docstring. We should use specialized geometry here.
            equal = self.parent().geometry._equal
            sgn = self.parent().geometry._sgn

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

        return super()._richcmp_(other, op)

    def __contains__(self, point):
        r"""
        Return whether ``point`` lies on this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        ::

            sage: g = H.geodesic(20, -479, 858, model="half_plane")
            sage: g.start() in g
            True

        We can often decide containment for points coming from geodesics::

            sage: g = H.geodesic(-1, -1/2)
            sage: h = H.geodesic(1, 1/2)

            sage: g.start() in g
            True
            sage: g.start() in h
            False
            sage: g.end() in g
            True
            sage: g.end() in h
            False

        .. NOTE::

            The implementation is currently not very robust over inexact rings.

        """
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        if isinstance(point, HyperbolicPointFromGeodesic):
            # Short cut the most common case (that _intersection cannot handle.)
            if point._geodesic.unoriented() == self.unoriented():
                return True

            intersection = self._intersection(point._geodesic)

            if intersection is None:
                return False

            if intersection.is_ultra_ideal():
                return False

            return intersection == point

        x, y = point.coordinates(model="klein")
        a, b, c = self.equation(model="klein")

        # We should use a specialized predicate from the geometry class here to
        # handle points that are close to the geodesic in a more robust way.
        return self.parent().geometry._zero(a + b * x + c * y)

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import ZZ

        return ZZ(1)

    def change(self, *, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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

    def geodesic(self):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self

    def __hash__(self):
        r"""
        Return a hash value for this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since oriented geodesics are hashable, they can be put in a hash table, such as a
        Python ``set``::

            sage: S = {H.vertical(0), -H.vertical(0)}
            sage: len(S)
            2

        Same for unoriented geodesics::

            sage: {H.vertical(0).unoriented(), (-H.vertical(0)).unoriented()}
            {{x = 0}}

        Oriented and unoriented geodesics are distinct and so is their hash
        value (typically)::

            sage: hash(H.vertical(0)) != hash(H.vertical(0).unoriented())
            True

        We can also mix oriented and unoriented geodesics in hash tables::

            sage: S = {H.vertical(0), -H.vertical(0), H.vertical(0).unoriented(), (-H.vertical(0)).unoriented()}
            sage: len(S)
            3

        """
        if not self.parent().base_ring().is_exact():
            raise TypeError("cannot hash geodesic defined over inexact base ring")

        return hash((type(self), self.equation(model="klein", normalization=["one", "gcd"])))

    def _intersection(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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

        xy = self.parent().geometry.intersection((self._a, self._b, self._c), (other._a, other._b, other._c))

        if xy is None:
            return None

        return self.parent().point(*xy, model="klein", check=False)

    def _isometry_equations(self, isometry, image, λ):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        R = λ.parent()

        a, b, c = self.equation(model="klein")
        fa, fb, fc = image.equation(model="klein")
        from sage.all import vector
        condition = vector((b, c, a)) * isometry - λ * vector(R, (fb, fc, fa))
        return condition.list()


class HyperbolicUnorientedGeodesic(HyperbolicGeodesic):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    An unoriented geodesic in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.vertical(0).unoriented()
        {x = 0}

    """

    def vertices(self, marked_vertices=True):
        r"""
        Return the ideal end points of this unoriented geodesic.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          geodesic cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).unoriented().vertices()
            {0, ∞}

        Currently, vertices cannot be computed if some of them have coordinates
        which do not live over the :meth:`HyperbolicPlane.base_ring`; see
        :class:`HyperbolicVertices`::

            sage: H.half_circle(0, 2).unoriented().vertices()
            Traceback (most recent call last):
            ...
            ValueError: ...

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return self.change(oriented=True).vertices()

    def _isometry_conditions(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        self = self.half_spaces()[0].boundary()
        other = other.half_spaces()[0].boundary()
        yield [(self, other)]
        yield [(self, -other)]


class HyperbolicOrientedGeodesic(HyperbolicGeodesic, HyperbolicOrientedConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    An oriented geodesic in the hyperbolic plane.

    Internally, we represent geodesics as the chords satisfying the equation `a
    + bx + cy=0` in the unit disk of the Klein model.

    The geodesic is oriented such that the half space `a + bx + cy ≥ 0` is on
    its left.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        r"""
        Return the reversed geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: -H.vertical(0)
            {x = 0}

        """
        return self.parent().geodesic(
            -self._a, -self._b, -self._c, model="klein", check=False
        )

    def start(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        Return the ideal starting point of this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        Return the ideal end point of this geodesic.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
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
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        Return a classification of the angle between this
        geodesic and ``other`` in the Klein model.

        """
        # TODO: Can we make this public somehow?
        intersection = self._intersection(other)

        if intersection is None:
            # TODO: Use a specialized predicate instead of the _method.
            orientation = self.parent().geometry._sgn(self._b * other._b + self._c * other._c)

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

    def an_element(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        other = self.parent().geodesic(0, -self._c, self._b, model="klein", check=False)
        cross = other._intersection(self)
        assert cross
        return cross

    def parametrize(self, point, model, check=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        if isinstance(point, HyperbolicPoint):
            if check and point not in self:
                raise ValueError("point must be on geodesic to be parametrized")

        if model == "euclidean":
            base = self.an_element().coordinates(model="klein")
            tangent = (self._c, -self._b)

            if isinstance(point, HyperbolicPoint):
                # TODO: Use a specialized predicate instead of the _method.
                coordinate = 0 if not self.parent().geometry._zero(tangent[0]) else 1
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

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        INPUT:

        - ``isometry`` -- a 3 x 3 matrix

        - ``on_right`` -- an optional boolean which default to ``False``; set it to
          ``True`` if you want a right action instead.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: p0 = H(0)
            sage: p1 = H(1)
            sage: p2 = H(oo)
            sage: for (a, b, c, d) in [(2, 1, 1, 1), (1, 1, 0, 1), (1, 0, 1, 1), (2, 0, 0 , 1)]:
            ....:     for on_right in [True, False]:
            ....:         m = matrix(2, [a, b, c, d])
            ....:         q0 = p0.apply_isometry(m, on_right=on_right)
            ....:         q1 = p1.apply_isometry(m, on_right=on_right)
            ....:         q2 = p2.apply_isometry(m, on_right=on_right)
            ....:         assert H.geodesic(p0, p1).apply_isometry(m, on_right=on_right) == H.geodesic(q0, q1)
            ....:         assert H.geodesic(p1, p0).apply_isometry(m, on_right=on_right) == H.geodesic(q1, q0)
            ....:         assert H.geodesic(p1, p2).apply_isometry(m, on_right=on_right) == H.geodesic(q1, q2)
            ....:         assert H.geodesic(p2, p1).apply_isometry(m, on_right=on_right) == H.geodesic(q2, q1)
            ....:         assert H.geodesic(p2, p0).apply_isometry(m, on_right=on_right) == H.geodesic(q2, q0)
            ....:         assert H.geodesic(p0, p2).apply_isometry(m, on_right=on_right) == H.geodesic(q0, q2)

        """
        from sage.modules.free_module_element import vector

        if not on_right:
            isometry = isometry.inverse()

        b, c, a = (
            vector(self.parent().base_ring(), [self._b, self._c, self._a])
            * isometry
        )
        return self.parent().geodesic(a, b, c, model="klein")

    def vertices(self, marked_vertices=True):
        r"""
        Return the ideal end points of this oriented geodesic.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          geodesic cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).vertices()
            {0, ∞}

        Note that iteration in the set is not consistent with the orientation
        of the geodesic (it is chosen such that the subset relation on vertices
        can be checked quickly)::

            sage: v = H.vertical(0)
            sage: list(v.vertices())
            [0, ∞]

            sage: list((-v).vertices())
            [0, ∞]

        Use :meth:`HyperbolicGeodesic.start` and
        :meth:`HyperbolicGeodesic.end` to get the end points in an order that
        is consistent with orientation::

            sage: v.start(), v.end()
            (0, ∞)

            sage: (-v).start(), (-v).end()
            (∞, 0)

        Currently, vertices cannot be computed if some of them have coordinates
        which do not live over the :meth:`HyperbolicPlane.base_ring`; see
        :class:`HyperbolicVertices`::

            sage: H.half_circle(0, 2).vertices()
            Traceback (most recent call last):
            ...
            ValueError: ...

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return HyperbolicVertices([self.start(), self.end()])

    def _isometry_conditions(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        yield [(self.start(), other.start()), (self.end(), other.end())]


class HyperbolicPoint(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
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
        # TODO: Benchmark?
        if self.is_ultra_ideal():
            raise ValueError(f"point {self.coordinates(model='klein')} is not in the unit disk in the Klein model")

    def is_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        x, y = self.coordinates(model="klein")
        # TODO: Use a specialized predicate instead of the _method.
        return self.parent().geometry._cmp(x * x + y * y, 1) == 0

    def is_ultra_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        x, y = self.coordinates(model="klein")
        # TODO: Use a specialized predicate instead of the _method.
        return self.parent().geometry._cmp(x * x + y * y, 1) > 0

    def is_finite(self):
        r"""
        Return whether this point is finite.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).is_finite()
            False

            sage: H(I).is_finite()
            True

            sage: H.half_circle(0, 2).start().is_finite()
            False

        .. NOTE::

            Currently, the implementation is not robust over inexact rings.

        """
        x, y = self.coordinates(model="klein")
        # We should use specialized predicate from the geometry implementation
        # here to make this more robust over inexact rings.
        return self.parent().geometry._cmp(x * x + y * y, 1) < 0

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return projective coordinates in the
        upper half plane model.

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
                # TODO: We should throw a special exception instead and drop
                # this keyword. Otherwise, we need to write a lot of
                # boilerplate to check for None return values everywhere.
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
        # TODO: Benchmark?
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
        # TODO: Benchmark?
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

    def segment(self, end, check=True, assume_normalized=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return the oriented segment from this point to ``end``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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
        # TODO: Benchmark?
        if model == "half_plane":
            from sage.all import point

            # We need to wrap the coordinates into a list so they are not
            # interpreted as a list of complex numbers.
            plot = point([self.coordinates(model="half_plane")], **kwds)
        elif model == "klein":
            from sage.all import point

            # We need to wrap the coordinates into a list so they are not
            # interpreted as a list of complex numbers.
            plot = point([self.coordinates(model="klein")], **kwds)
        else:
            raise NotImplementedError

        return self._enhance_plot(plot, model=model)

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import ZZ

        return ZZ.zero()

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this point, i.e., this point.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          point cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(oo).vertices()
            {∞,}

        Currently, this does not work for points that have coordinates that do
        not live in the :meth:`HyperbolicPlane.base_ring`; see
        :class:`HyperbolicVertices`::

            sage: H.half_circle(0, 2).start().vertices()
            Traceback (most recent call last):
            ...
            ValueError: ...

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return HyperbolicVertices([self])

    def __hash__(self):
        r"""
        Return a hash value for this point.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since points are hashable, they can be put in a hash table, such as a
        Python ``set``::

            sage: S = {H(I), H(0)}
            sage: len(S)
            2

        ::

            sage: {H.half_circle(0, 1).start(), H.half_circle(-2, 1).end()}
            {-1}

        Also, endpoints of geodesics that have no coordinates in the base ring
        can be hashed, see :meth:`HyperbolicPointFromGeodesic.__hash__`::

            sage: S = {H.half_circle(0, 2).start()}

        """
        if not self.parent().base_ring().is_exact():
            raise TypeError("cannot hash a point defined over inexact base ring")

        return hash(self.coordinates(ring=self.parent().base_ring(), model="klein"))

    def _isometry_equations(self, isometry, image, λ):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        R = λ.parent()

        x, y, z = (*self.coordinates(model="klein"), R.one())
        fx, fy, fz = (*image.coordinates(model="klein"), R.one())

        from sage.all import vector
        equations = λ * vector((x, y, z)) - isometry * vector(R, (fx, fy, fz))
        return equations.list()


class HyperbolicPointFromCoordinates(HyperbolicPoint):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    def __init__(self, parent, x, y):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return projective coordinates in the
        upper half plane model.

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
        r"""
        Return how this point compares to ``other`` with respect to the ``op``
        operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two points are the same.

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

        .. NOTE::

            Over inexact rings, this method is not very reliable. To some
            extent this is inherent to the problem but also the implementation
            uses generic predicates instead of relying on a specialized
            implementation in the :class:`HyperbolicGeometry`.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.__contains__` to check for containment
            of a point in a set

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicPoint):
                return False

            if isinstance(other, HyperbolicPointFromGeodesic):
                return other == self

            # See note in the docstring. We should use specialized geometry
            # here to compare the coordinates simultaneously.
            return all(self.parent().geometry._equal(a, b) for (a, b) in zip(self.coordinates(model="klein"), other.coordinates(model="klein")))

        return super()._richcmp_(other, op)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: geodesic = H.geodesic(733833/5522174119, -242010/5522174119, -105111/788882017, model="klein")
            sage: geodesic.start()
            3.03625883227966

        """
        if self == self.parent().infinity():
            return "∞"

        if self.is_ultra_ideal():
            return repr(self.coordinates(model="klein"))

        try:
            x, y = self.coordinates(model="half_plane")
        except ValueError:
            # TODO: Use a more specific exception here.
            # TODO: Unify with the code above?
            return repr(self.coordinates(model="klein"))
        else:
            from sage.all import PowerSeriesRing

            # We represent x + y*I in R[[I]] so we do not have to reimplement printing ourselves.
            return repr(PowerSeriesRing(self.parent().base_ring(), names="I")([x, y]))

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        def point(parent):
            return parent.point(*self._coordinates, model="klein", check=False)

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of a point")

        if ring is not None or geometry is not None:
            self = self.parent().change_ring(ring, geometry=geometry).point(*self._coordinates, model="klein", check=False)

        return self

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        INPUT:

        - ``isometry`` -- a 3x3 matrix in `SO(1, 2)`

        - ``model`` (optional) -- either ``"half_plane"`` (default) or ``"klein"``

        - ``on_right`` (optional; default to ``False``) -- set it to ``True`` if you want
          the right action.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: for (a, b, c, d) in [(2, 1, 1, 1), (1, 1, 0, 1), (1, 0, 1, 1), (2, 0, 0 , 1)]:
            ....:     m = matrix(2, [a, b, c, d])
            ....:     assert H(0).apply_isometry(m) == H(b / d if d else oo)
            ....:     assert H(1).apply_isometry(m) == H((a + b) / (c + d) if c+d else oo)
            ....:     assert H(oo).apply_isometry(m) == H(a / c if c else oo)
        """
        from sage.modules.free_module_element import vector

        x, y = self.coordinates(model="klein")

        if on_right:
            isometry = isometry.inverse()

        x, y, z = isometry * vector(self.parent().base_ring(), [x, y, 1])
        return self.parent().point(x / z, y / z, model="klein")


class HyperbolicPointFromGeodesic(HyperbolicPoint):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    def __init__(self, parent, geodesic):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("x must be an oriented geodesic")

        self._geodesic = geodesic

    def is_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return True

    def is_ultra_ideal(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return False

    def is_finite(self):
        r"""
        Return whether this ideal point is finite, i.e., return ``False``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_circle(0, 2).start().is_finite()
            False

        """
        return False

    @cached_method
    def coordinates(self, model="half_plane", ring=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Did we implement ring?
        if model == "klein":
            a, b, c = self._geodesic.equation(model="half_plane")

            # TODO: Use specialized predicates instead of the _methods.
            if self.parent().geometry._zero(a):
                point = None
                # TODO: Use specialized predicates instead of the _methods.
                if self.parent().geometry._sgn(b) > 0:
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

                # TODO: Use specialized predicates instead of the _methods.
                return self.parent().point(
                        (min if self.parent().geometry._sgn(a) > 0 else max)(endpoints),
                        0,
                        model="half_plane",
                        check=False,
                    ).coordinates(model=model, ring=ring)

        return super().coordinates(model=model, ring=ring)

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other`` with respect to the ``op``
        operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two points are the same.

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

        .. NOTE::

            Over inexact rings, this method is not very reliable. To some
            extent this is inherent to the problem but also the implementation
            uses generic predicates instead of relying on a specialized
            implementation in the :class:`HyperbolicGeometry`.

            For points that are not defined by coordinates but merely as the
            starting points of a hyperbolic geodesic, this is probably not
            implemented to the full extent possible. Instead, this will often
            throw a ``ValueError`` even though we could say something about the
            equality of the points.

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

                # This is probably too complicated. If we can compute the
                # intersection, then it should have coordinates in the base
                # ring.
                intersection = self._geodesic._intersection(other._geodesic)

                return self == intersection and other == intersection

            if not other.is_ideal():
                return False

            # We should probably be a little bit more careful over inexact
            # rings here.
            return other in self._geodesic and self._geodesic.parametrize(other, model="euclidean") < 0

        return super()._richcmp_(other, op)

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of a point")

        if ring is not None or geometry is not None:
            self = self._geodesic.change(ring=ring, geometry=geometry).start()

        return self

    def __hash__(self):
        r"""
        Return a hash value for this point.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Points that are given as the endpoints of a geodesic may or may not
        have coordinates over the base ring::

            sage: H.half_circle(0, 1).start().coordinates(model="klein")
            (-1, 0)
            sage: H.half_circle(0, 2).start().coordinates(model="klein")
            Traceback (most recent call last):
            ...
            ValueError: square root of 32 not in Rational Field

        While they always have coordinates in a quadratic extension, the hash
        of the coordinates in the extension might not be consistent with has
        values in the base ring, so we cannot simply hash the coordinates over
        some field extension::

            sage: hash(QQ(1/2)) == hash(AA(1/2))
            False

        To obtain consistent hash values for sets over the same base ring, at
        least if that base ring is a field, we observe that a point whose
        coordinates are not in the base ring cannot be the starting point of
        two different geodesics with an equation in the base ring. Indeed, for
        it otherwise had coordinates in the base ring as it were the
        intersection of these two geodesics and whence a solution to a linear
        equation with coefficients in the base ring. So, for points that have
        no coordinates in the base ring, we can hash the equation of the
        oriented geodesic to obtain a hash value::

            sage: hash(H.half_circle(0, 2).start()) != hash(H.half_circle(0, 2).end())
            True

        """
        if self.coordinates(model="klein", ring="try") is not None:
            return super().__hash__()

        from sage.categories.all import Fields
        if self.parent().base_ring() not in Fields():
            raise NotImplementedError("cannot hash point defined by a geodesic over a non-field yet")
        return hash((type(self), self._geodesic))

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        image = self._geodesic.apply_isometry(isometry, model="klein", on_right=on_right)
        if isometry.det().sign() == -1:
            image = -image

        return image.start()


class HyperbolicConvexPolygon(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`half spaces <HyperbolicHalfSpace>`.
    """

    def __init__(self, parent, half_spaces, vertices):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        super().__init__(parent)

        if not isinstance(half_spaces, HyperbolicHalfSpaces):
            raise TypeError("half_spaces must be HyperbolicHalfSpaces")

        self._half_spaces = half_spaces
        self._marked_vertices = tuple(vertices)

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Check _marked_vertices not vertices of the polygon if require_normalized
        for vertex in self._marked_vertices:
            if not any([vertex in edge for edge in self.edges()]):
                raise ValueError("marked vertex must be on an edge of the polygon")

    # TODO: Add examples.
    def _normalize(self, marked_vertices=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        An instance that caused problems at some point::

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
        marked_vertices = self._marked_vertices if marked_vertices else []

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
        self = self._normalize_drop_unit_disk_redundant()

        marked_vertices = [vertex for vertex in marked_vertices if vertex not in self.vertices()]

        if marked_vertices:
            if self.dimension() < 2:
                raise NotImplementedError("cannot add marked vertices to low dimensional objects")

            self = self.parent().polygon(
                 self.half_spaces(), check=False, assume_sorted=True, assume_minimal=True, marked_vertices=marked_vertices
            )
            self = self._normalize_drop_marked_vertices()

        return self

    def _normalize_drop_trivially_redundant(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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

                # TODO: Use specialized predicates instead of the _methods.
                equal = self.parent().geometry._equal
                sgn = self.parent().geometry._sgn
                if equal(c * B, C * b) and sgn(b) == sgn(B) and sgn(c) == sgn(C):
                    # The half spaces are parallel in the Euclidean plane. Since we
                    # assume spaces to be sorted by inclusion, we can drop this
                    # space.
                    continue

            reduced.append(half_space)

        return self.parent().polygon(
             reduced, check=False, assume_sorted=True, assume_minimal=True, marked_vertices=False
        )

    def _normalize_drop_euclidean_redundant(self, boundary):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
            required_half_spaces, check=False, assume_sorted=True, assume_minimal=True, marked_vertices=False
        )

    def _normalize_drop_unit_disk_redundant(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
            required_half_spaces, check=False, assume_sorted=True, assume_minimal=True, marked_vertices=False,
        )

    def _euclidean_boundary(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        is a single point in the Euclidean plane. The intersection in the
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

                    # TODO: Fix RealSet in SageMath to work with number fields.
                    if λ.parent().is_exact():
                        from sage.all import AA
                        rλ = AA(λ)
                    else:
                        rλ = λ

                    # Determine whether this half space constrains to (-∞, λ] or [λ, ∞).
                    if (
                        boundary.parametrize(λ + 1, model="euclidean", check=False)
                        in constraining
                    ):
                        constraint = RealSet.unbounded_above_closed(rλ)
                    else:
                        constraint = RealSet.unbounded_below_closed(rλ)

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
        # TODO: Benchmark?
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

    def _normalize_drop_marked_vertices(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        vertices = [vertex for vertex in self._marked_vertices if vertex not in self.vertices(marked_vertices=False)]

        return self.parent().polygon(
            self._half_spaces, check=False, assume_sorted=True, assume_minimal=True, marked_vertices=vertices,
        )

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import ZZ

        return ZZ(2)

    def equations(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return the equations describing the boundary of this polygon.

        The output is minimal and sorted by slope in the Klein model.
        """
        raise NotImplementedError

    @cached_method
    def edges(self, as_segments=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Define in HyperbolicConvexSet
        # TODO: Check that this works for marked vertices.
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
            elif AB == "anti-parallel":
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
            elif BC == "anti-parallel":
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

        # TODO: Return a SortedSet
        return edges

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this polygon, i.e., the (possibly ideal) end
        points of the :meth:`edges`.

        INPUT:

        -- ``marked_vertices`` -- a boolean (default: ``True``) whether to
        inlude marked vertices in the output

        OUTPUT:

        Returns a set of points. Iteration over this set is in counterclockwise
        order.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A finite polygon with a marked vertex::

            sage: P = H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ], marked_vertices=[I])
            sage: P.vertices()
            {-1, I, 1, (2/5, 3/5), (-2/5, 3/5)}

            sage: P.vertices(marked_vertices=False)
            {-1, 1, (2/5, 3/5), (-2/5, 3/5)}

        Currently, vertices cannot be computed if some of them have coordinates
        which do not live over the :meth:`HyperbolicPlane.base_ring`; see
        :class:`HyperbolicVertices`::

            sage: P = H.polygon([
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 2).right_half_space(),
            ....: ])
            sage: P.vertices()
            Traceback (most recent call last):
            ...
            ValueError: ...

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        vertices = []

        edges = self.edges()

        end = edges[-1].end()
        for i, edge in enumerate(edges):
            start = edge.start()
            if start != end:
                vertices.append(start)

            end = edge.end()
            vertices.append(end)

        if marked_vertices:
            vertices.extend(self._marked_vertices)

        return HyperbolicVertices(vertices)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self._half_spaces

    def _repr_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        half_spaces = " ∩ ".join([repr(half_space) for half_space in self._half_spaces])
        vertices = ", ".join([repr(vertex) for vertex in self._marked_vertices])
        if vertices:
            return f"{half_spaces} ∪ {{{vertices}}}"
        return half_spaces

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this polygon in the hyperbolic ``model``.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"``

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A finite triangle::

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.geodesic(I, I + 1).left_half_space(),
            ....:     H.geodesic(I + 1, 2*I).left_half_space()
            ....: ])
            sage: P
            {(x^2 + y^2) - x - 1 ≥ 0} ∩ {(x^2 + y^2) + 2*x - 4 ≤ 0} ∩ {x ≥ 0}

        In the upper half plane model, this plots as a polygon bounded by a
        segment and two arcs::

            sage: P.plot("half_plane")[0]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 1.00000000000000)),
                               CartesianPathPlotCommand(code='RARCTO', args=((1.00000000000000, 1.00000000000000), (0.500000000000000, 0))),
                               CartesianPathPlotCommand(code='ARCTO', args=((0.000000000000000, 2.00000000000000), (-1.00000000000000, 0))),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000))])

        Technically, the plot consists of two parts, a filled layer with
        transparent stroke and a transparent layer with solid stroke. The
        latter shows the (finite) edges of the polygon::

            sage: P.plot("half_plane")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 1.00000000000000)),
                               CartesianPathPlotCommand(code='RARCTO', args=((1.00000000000000, 1.00000000000000), (0.500000000000000, 0))),
                               CartesianPathPlotCommand(code='ARCTO', args=((0.000000000000000, 2.00000000000000), (-1.00000000000000, 0))),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000))])

        In the Klein disk model, this plots as two identical Euclidean
        triangles (one for the background, one for the edges,) with an added
        circle representing the ideal points in that model::

            sage: P.plot("klein")[0]
            Circle defined by (0.0,0.0) with r=1.0
            sage: P.plot("klein")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.666666666666667, 0.333333333333333)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.600000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])
            sage: P.plot("klein")[2]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.666666666666667, 0.333333333333333)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.600000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])

        An ideal triangle plots the same way in the Klein model but now has two
        rays in the upper half plane model::

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....:     H.vertical(1).left_half_space()
            ....: ])
            sage: P
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x ≥ 0}

            sage: P.plot("half_plane")[0]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='RARCTO', args=((1.00000000000000, 0.000000000000000), (0.500000000000000, 0))),
                               CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])

            sage: P.plot("klein")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(0.000000000000000, -1.00000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, -1.00000000000000))])

        A polygon can contain infinitely many ideal points as is the case in
        this intersection of two half spaces::

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space()
            ....: ])
            sage: P
            {x - 1 ≤ 0} ∩ {x ≥ 0}

            sage: P.plot("half_plane")[0]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(1.00000000000000, 0.000000000000000))])

        The last part, the line connecting 0 and 1, is missing from the
        stroke plot since we only stroke finite edges::

            sage: P.plot("half_plane")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])

        Simalarly in the Klein model picture, the arc of infinite points is
        only part of the fill, not of the stroke::

            sage: P.plot("klein")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, -1.00000000000000)),
                               CartesianPathPlotCommand(code='ARCTO', args=((1.00000000000000, 0.000000000000000), (0, 0)))])

            sage: P.plot("klein")[2]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 1.00000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, -1.00000000000000))])

        If the polygon contains unbounded set of reals, we get a horizontal ray
        in the half plane picture::

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_circle(2, 1).left_half_space(),
            ....: ])
            sage: P
            {(x^2 + y^2) - 4*x + 3 ≥ 0} ∩ {x ≥ 0}

            sage: P.plot("half_plane")[0]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='RARCTO', args=((3.00000000000000, 0.000000000000000), (2.00000000000000, 0))),
                               CartesianPathPlotCommand(code='RAYTO', args=(1, 0)),
                               CartesianPathPlotCommand(code='RAYTO', args=(0, 1)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='LINETO', args=(1.00000000000000, 0.000000000000000))])

            sage: P.plot("half_plane")[1]
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(1.00000000000000, 0.000000000000000)),
                               CartesianPathPlotCommand(code='RARCTO', args=((3.00000000000000, 0.000000000000000), (2.00000000000000, 0))),
                               CartesianPathPlotCommand(code='MOVETOINFINITY', args=(0, 1)),
                               CartesianPathPlotCommand(code='LINETO', args=(0.000000000000000, 0.000000000000000))])

        """
        # TODO: Document keyword arguments.
        kwds.setdefault("color", "#efffff")
        kwds.setdefault("edgecolor", "blue")

        if len(self._half_spaces) == 0:
            raise NotImplementedError("cannot plot full space")

        edges = self.edges(as_segments=True)

        pos = edges[0].start()

        commands = [HyperbolicPathPlotCommand("MOVETO", pos)]

        for edge in edges:
            if edge.start() != pos:
                commands.append(HyperbolicPathPlotCommand("MOVETO", edge.start()))

            commands.append(HyperbolicPathPlotCommand("LINETO", edge.end()))
            pos = edge.end()

        if pos != edges[0].start():
            commands.append(HyperbolicPathPlotCommand("MOVETO", edges[0].start()))

        plot = hyperbolic_path(commands, model=model, **kwds)

        return self._enhance_plot(
            plot, model=model
        )

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        r"""
        Return how this polygon compares to ``other`` with respect to the
        ``op`` operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether polygons are essentially indistinguishable.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space()])
            sage: P == P
            True

        Marked vertices are taken into account to determine equality::

            sage: Q = H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space()], marked_vertices=[I + 1])
            sage: Q == Q
            True
            sage: P == Q
            False

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, HyperbolicConvexPolygon):
                return False
            return self._half_spaces == other._half_spaces and self._marked_vertices == other._marked_vertices

        return super()._richcmp_(other, op)

    def __hash__(self):
        r"""
        Return a hash value for this polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since polygons are hashable, they can be put in a hash table, such
        as a Python ``set``::

            sage: S = {H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space()])}

        """
        return hash((self._half_spaces, self._marked_vertices))

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        half_spaces = [h.apply_isometry(isometry, model="klein", on_right=on_right) for h in self._half_spaces]
        marked_vertices = [p.apply_isometry(isometry, model="klein", on_right=on_right) for p in self._marked_vertices]

        return self.parent().polygon(half_spaces=half_spaces, check=False, assume_minimal=True, marked_vertices=marked_vertices)

    def _isometry_conditions(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?

        # TODO: Return half spaces instead so we get orientations right. However, then we are missing the marked vertices which can also be cyclically permuted :(
        if self._marked_vertices or other._marked_vertices:
            self = list(self.vertices())
            other = list(other.vertices())

            if len(self) == len(other):
                for i in range(len(self)):
                    yield list(zip(self, other[i:] + other[:i]))

                other.reverse()

                for i in range(len(self)):
                    yield list(zip(self, other[i:] + other[:i]))
        else:
            self = list(self.half_spaces())
            other = list(other.half_spaces())

            if len(self) == len(other):
                for i in range(len(self)):
                    yield list(zip(self, other[i:] + other[:i]))

                other.reverse()

                for i in range(len(self)):
                    yield list(zip(self, other[i:] + other[:i]))


class HyperbolicSegment(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
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
        # TODO: Benchmark?
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

    def _check(self, require_normalized=True):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        # TODO: Should this be in the oriented class? Should there be an equivalent in the unoriented class?
        r"""
        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

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

        def λ(point):
            return self._geodesic.parametrize(point, model="euclidean", check=False)

        if start is not None and end is not None:
            if λ(start) > λ(end):
                raise ValueError("end point of segment must be after start point on the underlying geodesic")

        if start is not None:
            if not start.is_finite():
                # TODO: Use specialized predicate instead of _method.
                if self.parent().geometry._sgn(λ(start)) > 0:
                    return (
                        self.parent().empty_set() if start.is_ultra_ideal() else start
                    )
                start = None

        if end is not None:
            if not end.is_finite():
                # TODO: Use specialized predicate instead of _method.
                if self.parent().geometry._sgn(λ(end)) < 0:
                    return self.parent().empty_set() if end.is_ultra_ideal() else end
                end = None

        if start is None and end is None:
            segment = self._geodesic
            if not self.is_oriented():
                segment = segment.unoriented()
            return segment

        assert (start is None or not start.is_ultra_ideal()) and (
            end is None or not end.is_ultra_ideal()
        )

        if start == end:
            return start

        return self.parent().segment(
            self._geodesic, start=start, end=end, check=False, assume_normalized=True, oriented=self.is_oriented()
        )

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: segment = H.segment(H.vertical(0), start=I)
            sage: segment.apply_isometry(matrix(2, [1, 1, 0, 1]))
            {-x + 1 = 0} ∩ {2*(x^2 + y^2) - 3*x - 1 ≥ 0}

        """
        geodesic = self.geodesic()._apply_isometry_klein(isometry, on_right=on_right)
        start = self._start._apply_isometry_klein(isometry, on_right=on_right) if self._start is not None else None
        end = self._end._apply_isometry_klein(isometry, on_right=on_right) if self._end is not None else None
        return self.parent().segment(geodesic, start=start, end=end)

    def _endpoint_half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        bounds = [repr(self._geodesic)]
        bounds.extend(repr(half_space) for half_space in self._endpoint_half_spaces())

        return " ∩ ".join(bounds)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self._geodesic.half_spaces() + HyperbolicHalfSpaces(
            self._endpoint_half_spaces()
        )

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import RR

        kwds["fill"] = False

        self = self.change_ring(RR)
        plot = hyperbolic_path(
            [
                HyperbolicPathPlotCommand("MOVETO", self.start()),
                HyperbolicPathPlotCommand("LINETO", self.end()),
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
        # TODO: Benchmark?
        if self._start is not None:
            return self._start

        return self._geodesic.start()

    def end(self, finite=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if self._end is not None:
            return self._end

        return self._geodesic.end()

    def _richcmp_(self, other, op):
        r"""
        Return how this segment compares to ``other`` with respect to the
        ``op`` operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two segments are essentially indistinguishable.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()

        Oriented segments are equal if they have the same start and end points::

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

        return super()._richcmp_(other, op)

    def change(self, ring=None, geometry=None, oriented=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        geodesic = self._geodesic
        if not self.is_oriented():
            geodesic = geodesic.unoriented()
        return geodesic

    def vertices(self, marked_vertices=True):
        r"""
        Return the end poinst of this segment.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          segment cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s.vertices()
            {I, 2*I}

        Note that iteration in the set is not consistent with the orientation
        of the segment (it is chosen such that the subset relation on vertices
        can be checked quickly)::

            sage: (-s).vertices()
            {I, 2*I}

        Use :meth:`start` and :meth:`end` to get the vertices in an order that
        is consistent with the orientation::

            sage: s.start(), s.end()
            (I, 2*I)
            sage: (-s).start(), (-s).end()
            (2*I, I)

        Both finite and ideal end points of the segment are returned::

            sage: s = H(-1).segment(I)
            sage: s.vertices()
            {-1, I}

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return HyperbolicVertices([self.start(), self.end()])

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        from sage.all import ZZ

        return ZZ(1)

    def midpoint(self):
        r"""
        Return the midpoint of this segment.

        ALGORITHM:

        We use the construction as explained on `Wikipedia
        <https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model#Compass_and_straightedge_constructions>`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(4*I)
            sage: s.midpoint()
            2*I

        ::

            sage: K.<a> = NumberField(x^2 - 2, embedding=1.414)
            sage: H = HyperbolicPlane(K)
            sage: s = H(I - 1).segment(I + 1)
            sage: s.midpoint()
            a*I

        .. SEEALSO::

            :meth:`HyperbolicConvexPolygon.centroid` for a generalization of this
            :meth:`HyperbolicSegment.perpendicular` for the perpendicular bisector

        """
        start, end = self.vertices()

        if start == end:
            return start

        if not start.is_finite() and not end.is_finite():
            return self.geodesic().midpoint()

        if not start.is_finite() or not end.is_finite():
            raise NotImplementedError(f"cannot compute midpoint of unbounded segment {self}")

        for p in self.geodesic().perpendicular(start).vertices():
            for q in self.geodesic().perpendicular(end).vertices():
                line = self.parent().geodesic(p, q)
                intersection = self.intersection(line)
                if intersection:
                    return intersection

            # One of the two lines start at any p must intersect the segment
            # already. No need to check the other p.
            assert False, f"segment {self} must have a midpoint but the straightedge and compass construction did not yield any"

    def perpendicular(self, point=None):
        r"""
        Return the geodesic through ``point`` that is perpendicular to this
        segment.

        If no point is given, return the perpendicular bisector of this
        segment.

        ALGORITHM:

        See :meth:`HyperbolicGeodesic.perpendicular`.

        INPUT:

        - ``point`` -- a point on this segment or ``None`` (the default)

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicSegment
            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(4*I)
            sage: s.perpendicular()
            {(x^2 + y^2) - 4 = 0}
            sage: s.perpendicular(I)
            {(x^2 + y^2) - 1 = 0}

        ::

            sage: K.<a> = NumberField(x^2 - 2, embedding=1.414)
            sage: H = HyperbolicPlane(K)
            sage: s = H(I - 1).segment(I + 1)
            sage: s.perpendicular()
            {x = 0}
            sage: s.perpendicular(I - 1)
            {4/3*(x^2 + y^2) + 16/3*x + 8/3 = 0}

        """
        if point is None:
            point = self.midpoint()

        point = self.parent()(point)

        if point not in self:
            raise ValueError(f"point must be in the segment but {point} is not in {self}")

        return self.geodesic().perpendicular(point)


class HyperbolicUnorientedSegment(HyperbolicSegment):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    An unoriented (possibly infinity) segment in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: segment = H.segment(H.vertical(0), start=I).unoriented()

    """

    def __hash__(self):
        r"""
        Return a hash value for this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(2*I)

        Since an uriented segment is hashable, it can be put in a hash table, such as a
        Python ``set``::

            sage: {s.unoriented(), (-s).unoriented()}
            {{-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}}

        """
        return hash((frozenset([self._start, self._end]), self.geodesic()))

    def _isometry_conditions(self, other):
        # TODO: Check documentation
        # TODO: Check INPUTS
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        yield [(self.start(), other.start()), (self.end(), other.end())]
        yield [(self.start(), other.end()), (self.end(), other.start())]


class HyperbolicOrientedSegment(HyperbolicSegment, HyperbolicOrientedConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.
    """

    def _neg_(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return self.parent().segment(
            -self._geodesic, self._end, self._start, check=False, assume_normalized=True
        )

    def _is_valid(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
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

    def __hash__(self):
        r"""
        Return a hash value for this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(2*I)

        Since this set is hashable, it can be put in a hash table, such as a
        Python ``set``::

            sage: {s}
            {{-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}}

        """
        return hash((self._start, self._end, self.geodesic()))

    def _isometry_conditions(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        yield [(self.start(), other.start()), (self.end(), other.end())]


class HyperbolicEmptySet(HyperbolicConvexSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    The empty subset of the hyperbolic plane.
    """

    def __init__(self, parent):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        super().__init__(parent)

    def _richcmp_(self, other, op):
        r"""
        Return how this set compares to ``other`` with respect to ``op``.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether both sets are empty.

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
        # TODO: Benchmark?
        return "{}"

    def _apply_isometry_klein(self, isometry, on_right=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        INPUT:

        - ``isometry`` -- a 3 x 3 matrix

        - ``on_right`` -- an optional boolean which default to ``False``; set it to
          ``True`` if you want a right action instead.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: S = HyperbolicPlane().empty_set()
            sage: S.apply_isometry(matrix(2, [2, 1, 1, 1])) is S
            True

        """
        return self

    def plot(self, model="half_plane", **kwds):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import Graphics

        return self._enhance_plot(Graphics(), model=model)

    def dimension(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from sage.all import ZZ

        return ZZ(-1)

    def half_spaces(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        if ring is not None:
            self = self.parent().change_ring(ring, geometry=geometry).empty_set()

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot change orientation of empty set")

        return self

    def __hash__(self):
        r"""
        Return a hash value for this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since this set is hashable, it can be put in a hash table, such as a
        Python ``set``::

            sage: {H.empty_set()}
            {{}}

        """
        return 0

    def vertices(self, marked_vertices=True):
        r"""
        Retrun the vertices of this empty, i.e., an empty set of points.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().vertices()
            {}

        """
        return HyperbolicVertices([])


def gl2_to_sim12(m):
    # TODO: Check documentation.
    # TODO: Check INPUT
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    r"""
    Return the lift of the 2x2 matrix ``m`` inside ``SO(1,2)``.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import gl2_to_sim12
        sage: gl2_to_sim12(matrix(2, [1,1,0,1]))
        [   1   -1    1]
        [   1  1/2  1/2]
        [   1 -1/2  3/2]
        sage: gl2_to_sim12(matrix(2, [1,0,1,1]))
        [   1    1    1]
        [  -1  1/2 -1/2]
        [   1  1/2  3/2]
        sage: gl2_to_sim12(matrix(2, [2,0,0,1/2]))
        [   1    0    0]
        [   0 17/8 15/8]
        [   0 15/8 17/8]

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


def sim12_to_gl2(m):
    r"""
    Inverse of :func:`gl2_to_sim12`.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import sim12_to_gl2, gl2_to_sim12
        sage: sim12_to_gl2(matrix(3, [1, -1, 1, 1, 1/2, 1/2, 1, -1/2, 3/2]))
        [1 1]
        [0 1]
        sage: sim12_to_gl2(matrix(3, [1, 1, 1, -1, 1/2, -1/2, 1, 1/2, 3/2]))
        [1 0]
        [1 1]
        sage: sim12_to_gl2(matrix(3, [1, 0, 0, 0, 17/8, 15/8, 0, 15/8, 17/8]))
        [  2   0]
        [  0 1/2]
        sage: sim12_to_gl2(gl2_to_sim12(matrix([[-1, 0], [0, 1]])))
        [ 1  0]
        [ 0 -1]
    """
    from sage.matrix.constructor import matrix

    if m.nrows() != 3 or m.ncols() != 3:
        raise ValueError("invalid matrix")
    m00, m01, m02, m10, m11, m12, m20, m21, m22 = m.list()
    K = m.base_ring()
    two = K(2)

    so12_det = m.determinant()
    sl2_det = so12_det.nth_root(3)
    if so12_det.sign() != sl2_det.sign():
        sl2_det *= -1

    a2 = (m12 + m22 + m21 + m11) / two
    b2 = (m12 + m22 - m21 - m11) / two
    c2 = (m21 + m22 - m12 - m11) / two
    d2 = (m11 + m22 - m12 - m21) / two

    ab = (m10 + m20) / two
    ac = (m01 + m02) / two
    ad = (m00 + sl2_det) / two
    bc = (m00 - sl2_det) / two
    bd = (m02 - m01) / two
    cd = (m20 - m10) / two

    # we recover +/- a takings square roots
    pm_a = a2.sqrt()
    pm_b = b2.sqrt()
    pm_c = c2.sqrt()
    pm_d = d2.sqrt()

    # TODO: do something less naive as there is no need to iterate
    import itertools
    for sa, sb, sc, sd in itertools.product([1, -1], repeat=4):
        a = sa * pm_a
        b = sb * pm_b
        c = sc * pm_c
        d = sd * pm_d
        if (a * b == ab and a * c == ac and a * d == ad and b * c == bc and b * d == bd and c * d == cd):
            return matrix(K, 2, 2, [a, b, c, d])

    raise ValueError('no projection to GL(2, R) in the base ring')


# TODO: Rename to OrderedSet?
# TODO: Allow creating from generator while allowing to answer length immediately.
class SortedSet:
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Add a collections.abc class.
    # TODO: Benchmark?
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
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        raise NotImplementedError

    def _merge(self, *sets):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Benchmark?
        r"""
        Return whether this set is equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
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
        # TODO: Benchmark?
        r"""
        Return whether this set is not equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).vertices() != H.vertical(1).vertices()
            True

        """
        return not (self == other)

    def __hash__(self):
        r"""
        Return a hash value for this set that is consistent with :meth:`__eq__`
        and :meth:`__ne__`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: hash(H.vertical(0).vertices()) != hash(H.vertical(1).vertices())
            True

        """
        return hash(tuple(self._entries))

    def __add__(self, other):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if type(self) is not type(other):
            raise TypeError
        entries = self._merge(list(self._entries), list(other._entries))
        return type(self)(entries, assume_sorted=True)

    def __repr__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Return a printable representation of this set.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.half_circle(0, 1).vertices()
            {-1, 1}

        """
        return "{" + repr(self._entries)[1:-1] + "}"

    def __iter__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return iter(self._entries)

    def __len__(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        return len(self._entries)

    def pairs(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        for i in range(len(self._entries)):
            yield self._entries[i-1], self._entries[i]

    def triples(self):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        for i in range(len(self._entries)):
            yield self._entries[i - 1], self._entries[i], self._entries[
                (i + 1) % len(self._entries)
            ]


class HyperbolicVertices(SortedSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    # TODO: Do we need to disable _merge in most cases (probably fine if both sets of vertices define the same polygon?)
    r"""
    A set of vertices on the boundary of a convex set in the hyperbolic plane,
    sorted in counterclockwise order.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane
        sage: H = HyperbolicPlane()
        sage: V = H.vertical(0).vertices()
        sage: V
        {0, ∞}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicVertices
        sage: isinstance(V, HyperbolicVertices)
        True

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.vertices`

    """

    def __init__(self, vertices, assume_sorted=None):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?

        if len(vertices) != 0:
            # We sort vertices in counterclockwise order. We need to fix a starting
            # point consistently, namely we choose the leftmost point in the Klein
            # model and if there are ties the one with minimal y coordinate.
            # That way we can then order vertices by their slopes with the starting
            # point to get a counterclockwise walk.
            self._start = min(vertices, key=lambda vertex: vertex.coordinates(model="klein"))

        # _lt_ needs to know the global structure of the convex hull of the vertices.
        # The base class constructor will replace _entries with a sorted version of _entries.
        self._entries = tuple(vertices)

        super().__init__(vertices, assume_sorted=assume_sorted)

    def _slope(self, vertex):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        sx, sy = self._start.coordinates(model="klein")
        x, y = vertex.coordinates(model="klein")
        return (y - sy, x - sx)

    def _lt_(self, lhs, rhs):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        if lhs == self._start:
            return True
        if rhs == self._start:
            return False
        if lhs == rhs:
            return False

        dy_lhs, dx_lhs = self._slope(lhs)
        dy_rhs, dx_rhs = self._slope(rhs)

        assert dx_lhs >= 0 and dx_rhs >= 0, "all points must be to the right of the starting point due to chosen normalization"

        if dy_lhs * dx_rhs < dy_rhs * dx_lhs:
            return True

        if dy_lhs * dx_rhs > dy_rhs * dx_lhs:
            return False

        # The points (start, lhs, rhs) are collinear.
        # In general we cannot decide their order with only looking at start,
        # lhs, and rhs. We need to understand where the rest of the convex hull
        # lives.
        assert lhs in self._entries and rhs in self._entries, "cannot compare vertices that are not defining for the convex hull"

        # If there is any vertex with a bigger slope, then this line is at the
        # start of the walk in counterclockwise order.
        for vertex in self._entries:
            dy, dx = self._slope(vertex)
            if dy * dx_rhs > dy_rhs * dx:
                return dx_lhs < dx_rhs
            elif dy * dx_rhs < dy_rhs * dx:
                return dx_lhs > dx_rhs

        raise ValueError("cannot decide counterclockwise ordering of three collinear points")


class HyperbolicHalfSpaces(SortedSet):
    # TODO: Check documentation
    # TODO: Check INPUTS
    # TODO: Check SEEALSO
    # TODO: Check for doctests
    # TODO: Benchmark?
    @classmethod
    def _lt_(cls, lhs, rhs):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: This is essentially atan2.
        r"""
        Return whether the half space ``lhs`` is smaller than ``rhs`` in a cyclic
        ordering of normal vectors, i.e., order half spaces by whether their
        normal points to the left/right, the slope of the geodesic, and finally
        by containment.

        This ordering is such that :meth:`HyperbolicPlane.intersection` can be
        computed in linear time for two hyperbolic convex sets.

        TESTS::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

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
            # TODO: Benchmark?
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

    @classmethod
    def convex_hull(cls, vertices):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Can we use the sorting that HyperbolicVertices provides?
        # TODO: Make this work when the points are not on the convex hull
        vertices = [vertex for (i, vertex) in enumerate(vertices) if vertex not in vertices[i + 1:]]
        if not vertices:
            raise NotImplementedError("cannot convert convex hull of no points yet")

        if all(vertex == vertices[0] for vertex in vertices):
            return vertices[0].half_spaces()

        parent = vertices[0].parent()

        half_spaces = []

        # Determine the half spaces defining the convex hull of the
        # vertices by iterating through the vertices in counterclockwise
        # order (this is using the assumption that all the vertices are on
        # the boundary of the convex hull.)
        reference = min(vertices, key=lambda vertex: -vertex.coordinates(model="klein")[1])
        rays = []
        for vertex in vertices:
            if vertex == reference:
                continue
            ray = parent.geodesic(reference, vertex).left_half_space()
            rays.append(ray)
            ray.vertex = vertex

        rays = HyperbolicHalfSpaces(rays)

        start = reference
        for ray in rays:
            end = ray.vertex
            half_spaces.append(parent.geodesic(start, end).left_half_space())
            start = end
        half_spaces.append(parent.geodesic(start, reference).left_half_space())

        return HyperbolicHalfSpaces(half_spaces)


class CartesianPathPlot(GraphicPrimitive):
    r"""
    A plotted path in the hyperbolic plane, i.e., a sequence of commands and
    associated control points in the hyperbolic plane.

    The ``plot`` methods of most hyperbolic convex sets rely on such a path.
    Usually, such a path should not be produced directly.

    This can be considered a more generic version of
    :class:`sage.plot.line.Line` and :class:`sage.plot.polygon.Polygon` since
    this is not limited to finite line segments. At the same time this
    generalizes matplotlib's ``Path`` somewhat, again by allowing infinite rays
    and lines.

    INPUT:

    - ``commands`` -- a sequence of :class:`HyperbolicPath.Command` describing
      the path.

    - ``options`` -- a dict or ``None`` (the default), options to affect the
      plotting of the path; the options accepted are the same that
      :class:`sage.plot.polygon.Polygon` accepts.

    EXAMPLES:

    A geodesic plot as a single such path (wrapped in a SageMath graphics
    object)::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, CartesianPathPlot, CartesianPathPlotCommand
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

            sage: from flatsurf.geometry.hyperbolic import CartesianPathPlot
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

            sage: from flatsurf.geometry.hyperbolic import CartesianPathPlot, CartesianPathPlotCommand

            sage: P = CartesianPathPlot([
            ....:     CartesianPathPlotCommand("MOVETO", (-1, 0)),
            ....:     CartesianPathPlotCommand("LINETO", (0, 0)),
            ....: ])
            sage: P
            CartesianPathPlot([CartesianPathPlotCommand(code='MOVETO', args=(-1, 0)), CartesianPathPlotCommand(code='LINETO', args=(0, 0))])

        """
        return f"CartesianPathPlot({self._commands})"

    def _render_on_subplot(self, subplot):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        r"""
        Render this path on ``subplot``.

        Matplotlib was not really made to draw things that extend to infinity.
        The trick here is to register a callback that redraws whenever the
        viewbox of the plot changes, e.g., as more objects are added to the
        plot or as the plot is dragged around.
        """
        # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
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

        edge_color = None
        if "edgecolor" in options:
            edge_color = to_mpl_color(options.pop("edgecolor"))

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
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
            # TODO: Benchmark?
            r"""
            Redraw after the viewport has been rescaled to make sure that
            infinite rays reach the end of the viewport.
            """
            patch.set_path(self._create_path(subplot.axes.get_xlim(), subplot.axes.get_ylim(), fill=fill))

        subplot.axes.callbacks.connect("ylim_changed", redraw)
        subplot.axes.callbacks.connect("xlim_changed", redraw)
        redraw()

    def _create_path(self, xlim, ylim, fill):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?

        # Handle the first command. This controls how the path starts.
        command = self._commands[0]

        if command.code == "MOVETO":
            pos = command.args
            direction = None
            vertices = [pos]
        elif command.code == "MOVETOINFINITY":
            direction = command.args
            # TODO: What's pos?
            vertices = [self._infinity(pos, direction, xlim, ylim)]
        else:
            raise RuntimeError(f"path must not start with a {command.code} command")

        from matplotlib.path import Path

        codes = [Path.MOVETO]

        for command in self._commands[1:]:
            pos, direction = self._extend_path(vertices, codes, pos, direction, command, fill, xlim, ylim)

        return Path(vertices, codes)

    @staticmethod
    def _infinity(pos, direction, xlim, ylim):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        from matplotlib.path import Path

        def extend(path):
            # TODO: Check documentation.
            # TODO: Check INPUT
            # TODO: Check SEEALSO
            # TODO: Check for doctests
            # TODO: Benchmark?
            vertices.extend(path.vertices[1:])
            codes.extend(path.codes[1:])

        if command.code == "LINETO":
            target = command.args

            if direction is not None:
                vertices.append(CartesianPathPlot._infinity(target, direction, xlim, ylim))
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

                # TODO: Is this the correct center?
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
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
        try:
            from matplotlib.transforms import Bbox

            bbox = Bbox.null()

            pos = None
            for command in self._commands:
                if command.code in ["MOVETO", "LINETO"]:
                    pos = command.args
                    bbox.update_from_data_xy([pos], ignore=False)
                elif command.code in ["ARCTO", "RARCTO"]:
                    target, center = command.args
                    # We simplify the computation of the bounding box here to
                    # speed things up. Since these are hyperbolic arcs, their
                    # center must be at y=0 and the endpoints at y≥0.
                    # bbox = bbox.union(
                    #     [
                    #         bbox,
                    #         self._arc_path(
                    #             center, target, pos, reverse=command.code == "RARCTO"
                    #         ).get_extents(),
                    #     ]
                    # )
                    bbox = bbox.union([
                        bbox,
                        Bbox.from_bounds(*pos, 0, 0),
                        Bbox.from_bounds(*target, 0, 0),
                    ])
                    from sage.all import sgn
                    if sgn(pos[0] - center[0]) != sgn(target[0] - center[0]):
                        # The bounding box includes a point that is higher
                        # than the endpoints of the arc.
                        from math import sqrt
                        bbox = bbox.union([
                            bbox,
                            Bbox.from_bounds(center[0], sqrt((pos[0] - center[0]) ** 2 + pos[1]**2), 0, 0)
                        ])

                    pos = target
                elif command.code in ["RAYTO", "MOVETOINFINITY"]:
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

    @classmethod
    def _arc_path(cls, center, start, end, reverse=False):
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
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

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicPathPlotCommand
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

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicPathPlotCommand
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
                raise ValueError("when no previous cursor position is specified, command must be MOVETO")

            if model == "half_plane" and self.target == self.target.parent().infinity():
                return [CartesianPathPlotCommand("MOVETOINFINITY", (0, 1))]

            from sage.all import RR
            return [CartesianPathPlotCommand("MOVETO", self.target.change_ring(RR).coordinates(model=model))]

        if self.code == "LINETO":
            return HyperbolicPathPlotCommand.create_segment_cartesian(cursor, self.target, model=model)

        if self.code == "MOVETO":
            return HyperbolicPathPlotCommand.create_move_cartesian(cursor, self.target, model=model, fill=fill, stroke=stroke)

        raise NotImplementedError("cannot convert this command to a Cartesian plot command yet")

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

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicPathPlotCommand
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
            cartesian_commands.extend(command.cartesian(model=model, cursor=cursor, stroke=stroke, fill=fill))
            cursor = command.target

        while cartesian_commands and cartesian_commands[-1].code.startswith("MOVE"):
            cartesian_commands.pop()

        return cartesian_commands

    @staticmethod
    def create_segment_cartesian(start, end, model):
        r"""
        Return a sequence of :class:`CartesianPathPlotCommand` that represent
        the closed boundary of a :class:`HyperbolicConvexPolygon`, namely the segment
        to ``end`` (from the previous position ``start``.)

        This is a helper function for :meth:`cartesian`.

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPlane, HyperbolicPathPlotCommand
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
        :meth:`_hyperbolic_move` which contains the entire positive real axis::

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

        """
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
        if start == end:
            raise ValueError(f"cannot draw segment from point {start} to itself ({end})")

        if model == "half_plane":
            # TODO: Probably all use of RR should be replaced with ball arithmetic to make things work on a very small scale?
            from sage.all import RR
            if start != start.parent().infinity():
                start_x, start_y = start.change_ring(RR).coordinates(model="half_plane")

            if end != end.parent().infinity():
                end_x, end_y = end.change_ring(RR).coordinates(model="half_plane")

            if start == start.parent().infinity():
                return [
                    CartesianPathPlotCommand("LINETO", (end_x, end_y)),
                ]

            if end == end.parent().infinity():
                return [
                    CartesianPathPlotCommand(
                        "RAYTO",
                        (0, 1),
                    )
                ]

            # TODO: Should we be more careful here?
            if (start_x - end_x).abs() < (start_y - end_y).abs() * 1e-6:
                # This segment is (almost) vertical. We plot it as if it were
                # vertical to avoid numeric issus.
                return [CartesianPathPlotCommand("LINETO", (end_x, end_y))]

            real_hyperbolic_plane = HyperbolicPlane(RR)
            geodesic = real_hyperbolic_plane.geodesic(
                    real_hyperbolic_plane.point(start_x, start_y, model="half_plane"),
                    real_hyperbolic_plane.point(end_x, end_y, model="half_plane"))
            center = (
                (geodesic.start().coordinates()[0] + geodesic.end().coordinates()[0])
                / 2,
                0,
            )

            return [
                CartesianPathPlotCommand("RARCTO" if start_x < end_x else "ARCTO", ((end_x, end_y), center))
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
        Return a list of :class:`CartesianPathPlot.Command` that
        represent the open "segment" on the boundary of a polygon
        connecting ``start`` and ``end``.

        This is a helper function for :meth:`make_cartesian`.
        """
        # TODO: Check documentation.
        # TODO: Check INPUT
        # TODO: Check SEEALSO
        # TODO: Check for doctests
        # TODO: Benchmark?
        # TODO: Use sage's vector or matplotlib builtins more so we do not need to implement basic geometric primitives manually here.
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

                return [CartesianPathPlotCommand("MOVETO", end.change_ring(RR).coordinates())]

            if start == start.parent().infinity():
                return [
                    CartesianPathPlotCommand("RAYTO", (-1, 0)),
                    CartesianPathPlotCommand("LINETO", end.change_ring(RR).coordinates()),
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
                    CartesianPathPlotCommand("LINETO", end.change_ring(RR).coordinates())
                ]
            else:
                return [
                    CartesianPathPlotCommand(
                        "RAYTO", (1, 0)
                    ),
                    CartesianPathPlotCommand(
                        "RAYTO", (0, 1)
                    ),
                    CartesianPathPlotCommand(
                        "RAYTO", (-1, 0)
                    ),
                    CartesianPathPlotCommand("LINETO", end.change_ring(RR).coordinates()),
                ]

        elif model == "klein":
            from sage.all import RR

            if fill:
                return [
                    CartesianPathPlotCommand(
                        "ARCTO", (end.change_ring(RR).coordinates(model="klein"), (0, 0))
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


        sage: from flatsurf.geometry.hyperbolic import CartesianPathPlotCommand
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
        if options.get("fill", None):
            g.add_primitive(CartesianPathPlot(HyperbolicPathPlotCommand.make_cartesian(commands, model=model, fill=True, stroke=False), {**options, 'thickness': 0}))
        g.add_primitive(CartesianPathPlot(HyperbolicPathPlotCommand.make_cartesian(commands, model=model, fill=False, stroke=True), {**options, 'fill': False}))
    except Exception as e:
        raise RuntimeError(f"Failed to render hyperbolic path {commands}", e)

    if options["legend_label"]:
        g.legend(True)
        g._legend_colors = [options["legend_color"]]
    return g
