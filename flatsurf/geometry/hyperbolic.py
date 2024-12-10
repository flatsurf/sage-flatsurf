r"""
Two dimensional hyperbolic geometry.

.. jupyter-execute::
    :hide-code:

    # Allow jupyter-execute blocks in this module to contain doctests
    import jupyter_doctest_tweaks

EXAMPLES::

    sage: from flatsurf import HyperbolicPlane

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
    four (or three). Similarly, an oriented geodesic (which cannot really be
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

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2023 Julian Rüth
#                           2022 Sam Freedman
#                           2022 Vincent Delecroix
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

import collections.abc

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


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

        sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicExactGeometry

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
                raise ValueError(
                    "geometry must be specified for HyperbolicPlane over inexact rings"
                )

        from sage.categories.all import Sets

        category = category or Sets()

        return super().__classcall__(
            cls, base_ring=base_ring, geometry=geometry, category=category
        )

    def __init__(self, base_ring, geometry, category):
        r"""
        Create the hyperbolic plane over ``base_ring``.

        TESTS::

            sage: from flatsurf import HyperbolicPlane

            sage: TestSuite(HyperbolicPlane(QQ)).run()
            sage: TestSuite(HyperbolicPlane(AA)).run()  # long time (.5s)
            sage: TestSuite(HyperbolicPlane(RR)).run()

        """
        from sage.all import RR

        if geometry.base_ring() is not base_ring:
            raise ValueError(
                f"geometry base ring must be base ring of hyperbolic plane but {geometry.base_ring()} is not {base_ring}"
            )

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

        Coercions between base rings induce coercion between hyperbolic planes::

            sage: from flatsurf import HyperbolicPlane

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
        Return whether the hyperbolic plane contains ``x``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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

        # pylint does not see the Cython parent() so we disable the import check.
        # pylint: disable=c-extension-no-member
        parent = sage.structure.element.parent(x)
        # pylint: enable=c-extension-no-member

        # Note that in old versions of SageMath (e.g. Sage 9.1), I
        # is not a number field element but a symbolic ring element.
        # The "parent is SR" part can probably removed at some point.
        if isinstance(parent, Parent) and parent in NumberFields() or parent is SR:
            if (
                x.real() in self.base_ring()
                and x.imag() in self.base_ring()
                and x.imag() >= 0
            ):
                return True

        return super().__contains__(x)

    def change_ring(self, ring, geometry=None):
        r"""
        Return the hyperbolic plane over a different base ``ring``.

        INPUT:

        - ``ring`` -- a ring or ``None``; if ``None``, uses the current
          :meth:`~HyperbolicPlane.base_ring`.

        - ``geometry`` -- a geometry or ``None``; if ``None``, tries to convert
          the existing geometry to ``ring``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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
        Return an element of the hyperbolic plane (mostly for testing).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: HyperbolicPlane().an_element()
            0

        """
        return self.real(0)

    def some_subsets(self):
        r"""
        Return some subsets of the hyperbolic plane for testing.

        Some of the returned sets are elements of the hyperbolic plane (i.e.,
        points) some are parents themselves, e.g., polygons.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: HyperbolicPlane().some_elements()
            [∞, 0, 1, -1, ...]

        """
        from sage.all import ZZ

        elements = self.some_elements()

        elements += [
            self.empty_set(),
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

    def some_elements(self):
        r"""
        Return some representative elements, i.e., points of the hyperbolic
        plane for testing.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: HyperbolicPlane().some_elements()
            [∞, 0, 1, -1, ...]

        """
        return [
            self.infinity(),
            self.real(0),
            self.real(1),
            self.real(-1),
            self.geodesic(0, 2).start(),
            self.half_circle(0, 2).start(),
        ]

    def _test_some_subsets(self, tester=None, **options):
        r"""
        Run test suite on some representative convex subsets.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane
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
        kinds = {
            "empty_set": HyperbolicEmptySet,
            "point": HyperbolicPointFromCoordinates,
            "oriented geodesic": HyperbolicOrientedGeodesic,
            "unoriented geodesic": HyperbolicUnorientedGeodesic,
            "half_space": HyperbolicHalfSpace,
            "oriented segment": HyperbolicOrientedSegment,
            "unoriented segment": HyperbolicUnorientedSegment,
            "polygon": HyperbolicConvexPolygon,
        }

        if kind is None:
            from sage.all import randint

            kind = list(kinds.keys())[randint(0, len(kinds) - 1)]

        if kind not in kinds:
            raise ValueError(f"kind must be one of {kinds}")

        return kinds[kind].random_set(self)

    def __call__(self, x):
        r"""
        Return ``x`` as an element of the hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: H(1)
            1

        We need to override this method. The normal code path in SageMath
        requires the argument to be an Element but facade sets are not
        elements::

            sage: v = H.vertical(0)

            sage: Parent.__call__(H, v)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert HyperbolicOrientedGeodesic_with_category_with_category to sage.structure.element.Element

            sage: H(v)
            {-x = 0}

        """
        if isinstance(x, HyperbolicConvexFacade):
            return self._element_constructor_(x)

        return super().__call__(x)

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as an element of the hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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

        if isinstance(x, HyperbolicConvexSet):
            return x.change(ring=self.base_ring(), geometry=self.geometry)

        if x in self.base_ring():
            return self.real(x)

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

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane

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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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
            point = self.__make_element_class__(HyperbolicPointFromCoordinates)(
                self, x, y
            )
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            if (
                isinstance(a, HyperbolicPointFromGeodesic)
                and isinstance(b, HyperbolicPointFromGeodesic)
                and a._geodesic == -b._geodesic
            ):
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

            sage: from flatsurf import HyperbolicPlane
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
          obtained from the :meth:`HyperbolicGeodesic._intersection` of
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

            sage: from flatsurf import HyperbolicPlane
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
            ValueError: end point of segment must not be before start point on the underlying geodesic

            sage: H.segment(H.vertical(0), start=I, end=2*I)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}

            sage: H.segment(H.vertical(0).unoriented(), start=2*I, end=I)
            {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

            sage: H.segment(H.vertical(0), start=2*I, end=I)
            Traceback (most recent call last):
            ...
            ValueError: end point of segment must not be before start point on the underlying geodesic

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
                raise ValueError(
                    "cannot deduce segment from single endpoint on an unoriented geodesic"
                )
            elif geodesic.parametrize(
                start, model="euclidean", check=False
            ) > geodesic.parametrize(end, model="euclidean", check=False):
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
        self,
        half_spaces,
        check=True,
        assume_sorted=False,
        assume_minimal=False,
        marked_vertices=(),
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
          :meth:`HyperbolicHalfSpaces._lt_`. When set, we omit sorting the
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

            sage: from flatsurf import HyperbolicPlane
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
            <class 'flatsurf.geometry.hyperbolic.HyperbolicEmptySet_with_category_with_category'>

        ::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ])
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicHalfSpace_with_category_with_category'>

        If we add a marked point to such a half space, the underlying type is a
        polygon again::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[I])
            sage: half_space
            {(x^2 + y^2) - 1 ≤ 0} ∪ {I}
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicConvexPolygon_with_category_with_category'>

        Marked points that coincide with vertices are ignored::

            sage: half_space = H.polygon([
            ....:   H.half_circle(0, 1).right_half_space(),
            ....: ], marked_vertices=[-1])
            sage: half_space
            {(x^2 + y^2) - 1 ≤ 0}
            sage: type(half_space)
            <class 'flatsurf.geometry.hyperbolic.HyperbolicHalfSpace_with_category_with_category'>

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
        spaces by :meth:`HyperbolicHalfSpaces._lt_`. If we know that the
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

            sage: from flatsurf import HyperbolicPlane
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

        TESTS:

        A trivial case that did not work initially:

            sage: H.convex_hull(H(0), H(1), H(oo), H(I), H(I + 1))
            {(x^2 + y^2) - x ≥ 0} ∩ {x - 1 ≤ 0} ∩ {x ≥ 0}

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
                    for a, b in subset.vertices().pairs():
                        if a.segment(b) not in edges:
                            half_spaces.append(self.geodesic(a, b).right_half_space())
                else:
                    raise NotImplementedError(
                        "cannot form convex hull of this kind of set yet"
                    )

        edges = []

        for edge in polygon.half_spaces():
            for half_space in half_spaces:
                if (-edge).is_subset(half_space):
                    break
            else:
                edges.append(edge)

        if marked_vertices:
            marked_vertices = [
                vertex
                for vertex in vertices
                if any(vertex in half_space.boundary() for half_space in edges)
            ]

        polygon = self.polygon(edges, marked_vertices=marked_vertices)

        assert all(
            subset.is_subset(polygon) for subset in subsets
        ), "convex hull does not contain all the sets it is supposed to be the convex hull of"

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

            sage: from flatsurf import HyperbolicPlane
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
            :meth:`HyperbolicPlane.convex_hull` to compute the convex hull of subspaces

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

            sage: from flatsurf import HyperbolicPlane

            sage: HyperbolicPlane().empty_set()
            {}

        """
        return self.__make_element_class__(HyperbolicEmptySet)(self)

    def isometry(
        self, preimage, image, model="half_plane", on_right=False, normalized=False
    ):
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
        See :meth:`HyperbolicConvexSet.apply_isometry` for meaning of this
        matrix.

        ALGORITHM:

        We compute an isometry with a very inefficient Gröbner basis approach.

        Essentially, we try to extract three points from the ``preimage`` with
        their prescribed images in ``image``, see :meth:`_isometry_conditions`
        and determine the unique isometry mapping the points by solving the
        corresponding polynomial system, see :meth:`_isometry_from_equations`
        for the hacky Gröbner basis bit.

        There are a lot of problems with this approach (apart from it being
        extremely slow).

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

            sage: from flatsurf import HyperbolicPlane
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

            sage: H.isometry([0, 1, oo, I], [0, 1, oo, I + 1])  # long time (.4s)
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
            sage: H.isometry(P, Q)  # long time (1s)
            Traceback (most recent call last):
            ...
            ValueError: no isometry can map these objects to each other

            sage: Q = H.polygon(P.half_spaces(), marked_vertices=[-1 + I])
            sage: H.isometry(P, Q)  # long time (1s)
            [ 1  0]
            [ 0 -1]

        We can explicitly ask for an isometry in the Klein model, given by a
        3×3 matrix::

            sage: H.isometry(P, Q, model="klein")  # long time (1s)
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
            sage: H.isometry((x, y, z), (x.apply_isometry(isometry), y.apply_isometry(isometry), z.apply_isometry(isometry)))  # long time (.3s)
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
            sage: H.isometry(P, Q)  # long time (.4s)
            [  0  -1]
            [1/2   1]

        An isometry mapping unoriented segments, though not the most apparent
        one::

            sage: H.isometry(H(I).segment(2*I).unoriented(), H(2*I).segment(I).unoriented())
            [  0  -1]
            [1/2   0]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.apply_isometry` to apply the returned
            isometry to a convex set.

        """
        if normalized:
            isometry = self.isometry(
                preimage=preimage,
                image=image,
                model=model,
                on_right=on_right,
                normalized=False,
            )
            det = abs(isometry.det())
            λ = det.nth_root(isometry.nrows())
            return ~λ * isometry

        if model == "klein":
            isometry = self.isometry(
                preimage=preimage,
                image=image,
                model="half_plane",
                on_right=on_right,
                normalized=normalized,
            )
            return self._isometry_gl2_to_sim12(isometry)
        elif model != "half_plane":
            raise NotImplementedError("unsupported model")

        if not on_right:
            isometry = self.isometry(
                preimage=preimage,
                image=image,
                model=model,
                on_right=True,
                normalized=normalized,
            )
            if model == "half_plane":
                from sage.all import matrix

                # Pick a nice representative of the inverse matrix.
                isometry = matrix(
                    [
                        [isometry[1][1], -isometry[0][1]],
                        [-isometry[1][0], isometry[0][0]],
                    ]
                )
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
            raise ValueError(
                "preimage and image must be the same size to determine an isometry between them"
            )

        for x, y in zip(preimage, image):
            if x.dimension() != y.dimension():
                raise ValueError(
                    "preimage and image must be of the same dimensions to determine an isometry between them"
                )
            if x.is_oriented() != y.is_oriented():
                raise ValueError(
                    "preimage and image must be oriented or unoriented consistently"
                )

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
        # itself at the same time). Similarly, when mapping a polygon, we can
        # permute the edges cyclically.

        # We need a mild form of backtracking to collect all possible triples
        # that define the isometry.
        isometry = self._isometry_from_pairs(list(zip(preimage, image)))

        if isometry is None:
            raise ValueError("no isometry can map these objects to each other")

        return isometry

    def _isometry_gl2_to_sim12(self, isometry):
        r"""
        Return a lift of the ``isometry`` to a 3×3 matrix in the similitude
        group `\mathrm{Sim}(1, 2)` describing an isometry in hyperboloid
        model.

        This is a helper method for :meth:`isometry` and
        :meth:`HyperbolicConvexSet.apply_isometry` since isometries in the
        hyperboloid model can be directly applied to our objects which we
        represent in the Klein model.

        INPUT:

        - ``isometry`` -- a 2×2 matrix with non-zero determinant

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H._isometry_gl2_to_sim12(matrix(2, [1,1,0,1]))
            [   1   -1    1]
            [   1  1/2  1/2]
            [   1 -1/2  3/2]
            sage: H._isometry_gl2_to_sim12(matrix(2, [1,0,1,1]))
            [   1    1    1]
            [  -1  1/2 -1/2]
            [   1  1/2  3/2]
            sage: H._isometry_gl2_to_sim12(matrix(2, [2,0,0,1/2]))
            [   1    0    0]
            [   0 17/8 15/8]
            [   0 15/8 17/8]

        .. SEEALSO::

            :meth:`_isometry_sim12_to_gl2` for an inverse of this construction

        """
        from sage.matrix.constructor import matrix

        if isometry.dimensions() != (2, 2):
            raise ValueError(
                "matrix does not encode an isometry in the half plane model"
            )

        a, b, c, d = isometry.list()
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

    def _isometry_sim12_to_gl2(self, isometry):
        r"""
        Return an invertible 2×2 matrix that encodes the same isometry as
        ``isometry``.

        INPUT:

        - ``isometry`` -- a 3×3 matrix in `\mathrm{Sim}(1, 2)`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H._isometry_sim12_to_gl2(matrix(3, [1, -1, 1, 1, 1/2, 1/2, 1, -1/2, 3/2]))
            [1 1]
            [0 1]

            sage: H._isometry_sim12_to_gl2(matrix(3, [1, 1, 1, -1, 1/2, -1/2, 1, 1/2, 3/2]))
            [1 0]
            [1 1]

            sage: H._isometry_sim12_to_gl2(matrix(3, [1, 0, 0, 0, 17/8, 15/8, 0, 15/8, 17/8]))
            [  2   0]
            [  0 1/2]

            sage: H._isometry_sim12_to_gl2(H._isometry_gl2_to_sim12(matrix([[-1, 0], [0, 1]])))
            [ 1  0]
            [ 0 -1]

        .. SEEALSO::

            :meth:`_isometry_gl2_to_sim12` for an inverse of this construction

        """
        from sage.matrix.constructor import matrix

        if isometry.dimensions() != (3, 3):
            raise ValueError("matrix does not encode an isometry in the Klein model")

        m00, m01, m02, m10, m11, m12, m20, m21, m22 = isometry.list()
        K = isometry.base_ring()
        two = K(2)

        so12_det = isometry.determinant()
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

        # we recover +/- a taking square roots
        pm_a = a2.sqrt()
        pm_b = b2.sqrt()
        pm_c = c2.sqrt()
        pm_d = d2.sqrt()

        # We could do something less naive as there is no need to iterate
        import itertools

        for sa, sb, sc, sd in itertools.product([1, -1], repeat=4):
            a = sa * pm_a
            b = sb * pm_b
            c = sc * pm_c
            d = sd * pm_d
            if (
                a * b == ab
                and a * c == ac
                and a * d == ad
                and b * c == bc
                and b * d == bd
                and c * d == cd
            ):
                return matrix(K, 2, 2, [a, b, c, d])

        raise ValueError("no projection to GL(2, R) in the base ring")

    def _isometry_from_pairs(self, pairs):
        r"""
        Return a right isometry compatible with given (preimage, image) pairs.

        This is helper method for :meth:`isometry`.

        ALGORITHM:

        We extract pairs of primitive elements (points and geodesics) from
        ``pairs`` that necessarily need to map to each other for the mappings
        in ``pairs`` to be satisfied. (Based on the idea that there is only one
        isometry mapping three points in a prescribed way.) For these pairs of
        primitive elements, we determine the isometries mapping them (there
        might be more than one when we cannot extract three points) and check
        whether they map all ``pairs`` correctly.

        There is a bit of backtracking needed in :meth:`_isometry_conditions`
        since, e.g., there is no a unique isometry mapping two polygons to each
        other since we can, e.g., permute the edges of the polygon cyclically.

        INPUT:

        - ``pairs`` -- a sequence of pairs of hyperbolic sets

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H._isometry_from_pairs([(H(0), H(1)), (H(1), H(2)), (H(2), H(3))])
            [ 1 -1]
            [ 0  1]

        """
        for conditions in self._isometry_conditions([], pairs):
            for isometry in self._isometry_from_primitives(conditions):
                if isometry is None:
                    continue

                if any(
                    preimage.apply_isometry(isometry, on_right=True) != image
                    for (preimage, image) in pairs
                ):
                    continue

                return isometry

        return None

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

            sage: from flatsurf import HyperbolicPlane
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

        # Create polynomial equations that must be satisfied to map the pairs
        # to each other.
        def equations(isometry, λ):
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

            sage: from flatsurf import HyperbolicPlane
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
                if any(preimage in existing for existing in existings):
                    return None
            return (preimage, image)

        elif preimage.dimension() == 1:
            if preimage.unoriented() in [
                existing.unoriented() for existing in existings
            ]:
                # Again, we ignore the distinction between oriented and
                # unoriented geodesics here.
                return None

            if preimage.start() in existings:
                return self._isometry_untrivialize(
                    preimage.end(), image.end(), defining
                )

            if preimage.end() in existings:
                return self._isometry_untrivialize(
                    preimage.start(), image.start(), defining
                )

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

            sage: from flatsurf import HyperbolicPlane
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
            [[({-x + 1 = 0}, {x + 1 = 0}), ({x + 1 = 0}, {(x^2 + y^2) - 1 = 0})],
             [({-x + 1 = 0}, {(x^2 + y^2) - 1 = 0}), ({x + 1 = 0}, {-x + 1 = 0})],
             [({-x + 1 = 0}, {-x + 1 = 0}), ({x + 1 = 0}, {x + 1 = 0})],
             [({-x + 1 = 0}, {x + 1 = 0}), ({x + 1 = 0}, {-x + 1 = 0})],
             [({-x + 1 = 0}, {-x + 1 = 0}), ({x + 1 = 0}, {(x^2 + y^2) - 1 = 0})],
             [({-x + 1 = 0}, {(x^2 + y^2) - 1 = 0}), ({x + 1 = 0}, {x + 1 = 0})]]

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
            if isinstance(x, HyperbolicPoint):
                assert y.dimension() == 0

                pair = self._isometry_untrivialize(x, y, defining)
                if pair:
                    defining.append(pair)

                yield from self._isometry_conditions(defining[:], remaining)

            # Extend with a pair of geodesics in "remaining[0]"
            elif isinstance(x, HyperbolicOrientedGeodesic):
                assert y.dimension() == 1

                f = x.geodesic()
                g = y.geodesic()

                pair = self._isometry_untrivialize(f, g, defining)
                if pair:
                    defining.append(pair)

                yield from self._isometry_conditions(
                    defining[:],
                    remaining + [(x.start(), y.start()), (x.end(), y.end())],
                )

            # Extend with points coming from other hyperbolic objects in "remaining[0]"
            else:
                for pairs in x._isometry_conditions(y):
                    yield from self._isometry_conditions(defining[:], pairs + remaining)

        else:
            yield defining

    def _isometry_from_single_points(self, preimage, image):
        r"""
        Helper method for :meth:`isometry`.

        Return a right isometry that maps the point ``preimage`` to the point
        ``image`` or ``None`` when no such isometry exists.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
            isometry *= ~MS(
                [[image.coordinates()[1] / preimage.coordinates()[1], 0], [0, 1]]
            )
            isometry *= ~MS(
                [
                    [
                        1,
                        image.coordinates()[0]
                        - preimage.apply_isometry(
                            isometry, on_right=True
                        ).coordinates()[0],
                    ],
                    [0, 1],
                ]
            )
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

            sage: from flatsurf import HyperbolicPlane
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
        for isometry in self._isometry_from_primitives(
            [(preimage, image), (preimage.midpoint(), image.midpoint())]
        ):
            if isometry is None:
                import warnings

                warnings.warn(
                    "Could not determine an isometry of geodesics over the base ring. There might still be one but the implementation failed to detect it."
                )
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
                R = PolynomialRing(
                    self.base_ring(), names=variables, order="degrevlex(5), lex(1)"
                )
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
                isometry = matrix([[a, b], [c, d]])

                isometry = self._isometry_gl2_to_sim12(isometry)

                # Build equations for the symbolic variables and make sure that
                # the resulting variety is zero-dimensional.
                equations = conditions(isometry, (λ0, λ1, λ2))

                for λ in [λ1, λ2]:
                    if all(equation.degree(λ) <= 0 for equation in equations):
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

                assert equation.variables() == (
                    variable,
                ), f"expected Gröbner basis algorithm to yield an equation for {variable} but found {equation} instead"
                equation = equation.polynomial(variable)
                equation = equation.map_coefficients(
                    self.base_ring(), new_base_ring=self.base_ring()
                )
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

                    assert (
                        λ0 in self.base_ring() and λ0 != 0
                    ), f"did not deduce a non-zero constant for λ0={λ0} from {equation}"

                    # We could now patch the existing Gröbner basis and solve directly, but
                    # we just solve again for the correct value of λ0.
                    equations = conditions(isometry, (λ0, λ1, λ2))

                    for λ in [λ1, λ2]:
                        if all(equation.degree(λ) <= 0 for equation in equations):
                            # Force the unused variable λ to be =0
                            equations.append(λ)

                    solutions = R.ideal(equations).variety()

                    assert (
                        solutions
                    ), "After tuning the constant of the equations describing the isometry, there should be a solution but we did not find any."

                    solutions = [
                        matrix([[solution[a], solution[b]], [solution[c], solution[d]]])
                        for solution in solutions
                    ]

                    # Check which solutions do not only satisfy "conditions"
                    # but map the underlying objects correctly, i.e., they also
                    # pass "filter".
                    solutions = [solution for solution in solutions if filter(solution)]

                    if not solutions:
                        continue

                    # Prefer an isometry of determinant 1 and isometries with 1 entries.
                    return max(
                        solutions,
                        key=lambda isometry: (
                            isometry.det(),
                            isometry[1][1] == 1,
                            isometry[1][0] == 1,
                        ),
                    )

        # No luck with this approach. We hope that this means that no such
        # isometry exists over the base ring but it's not entirely clear
        # whether that is actually true.
        return None

    def _repr_(self):
        r"""
        Return a printable representation of this hyperbolic plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeometry
            sage: H = HyperbolicPlane()
            sage: isinstance(H.geometry, HyperbolicGeometry)
            True

        """
        self._ring = ring

    def base_ring(self):
        r"""
        Return the ring over which this geometry is implemented.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

        assert (
            x > y
        ), "Geometry over this ring must override _cmp since not (x == y) and not (x < y) does not imply x > y"
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

    def change_ring(self, ring):
        r"""
        Return this geometry with the :meth:`base_ring` changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

        sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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
        Return this geometry with the :meth:`~HyperbolicGeometry.base_ring`
        changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic import HyperbolicEpsilonGeometry
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

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicEpsilonGeometry
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

            sage: from flatsurf import HyperbolicPlane
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

        - ``a`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        - ``b`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        - ``c`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        - ``d`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        EXAMPLES:

            sage: from flatsurf import HyperbolicPlane
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

        - ``p`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        - ``q`` -- an element of the :meth:`~HyperbolicGeometry.base_ring`

        - ``point`` -- the :meth:`HyperbolicPlane.point` to create points

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(RR)
            sage: H.geometry
            Epsilon geometry with ϵ=1.00000000000000e-6 over Real Field with 53 bits of precision

        """
        return f"Epsilon geometry with ϵ={self._epsilon} over {self._ring}"


class HyperbolicConvexSet(SageObject):
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

        Comparison with ``==`` should mean "is essentially indistinguishable
        from": Implementing == to mean anything else would get us into trouble
        in the long run. In particular we cannot implement <= to mean "is
        subset of" since then an oriented and an unoriented geodesic would be
        `==`. So, objects of a different type should almost never be equal. A
        notable exception are objects that are indistinguishable to the end
        user but use different implementations: the starting point of the
        geodesic going from 0 to infinity, a
        :class:`HyperbolicPointFromGeodesic`, and the point with coordinates
        (0, 0) in the upper half plane model, a
        :class:`HyperbolicPointFromCoordinates`, are equal. Note that we also
        treat objects as equal that only differ in their exact representation
        such as the geodesic x = 1 and the geodesic 2x = 2.

    TESTS::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexSet
        sage: H = HyperbolicPlane()

        sage: isinstance(H(0), HyperbolicConvexSet)
        True

    """

    def half_spaces(self):
        r"""
        Return a minimal set of half spaces whose intersection is this convex set.

        Iteration of the half spaces is in counterclockwise order, i.e.,
        consistent with :meth:`HyperbolicHalfSpaces._lt_`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_half_spaces()

        """
        tester = self._tester(**options)

        try:
            half_spaces = self.half_spaces()
        except ValueError:
            # We cannot create half spaces defining a point which has no coordinates over the base ring.
            tester.assertIsInstance(self, HyperbolicPointFromGeodesic)
            return

        tester.assertEqual(self.parent().intersection(*half_spaces), self.unoriented())

        tester.assertTrue(isinstance(half_spaces, HyperbolicHalfSpaces))

        for a, b in zip(list(half_spaces), list(half_spaces)[1:]):
            tester.assertTrue(HyperbolicHalfSpaces._lt_(a, b))

    def _check(self, require_normalized=True):
        r"""
        Validate this convex subset.

        Subclasses run specific checks here that can be disabled when creating
        objects with ``check=False``.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); whether to
          include checks that assume that normalization has already happened

        EXAMPLES:

            sage: from flatsurf import HyperbolicPlane
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
        then this is actually a geodesic and it should be described as such.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).unoriented()
            {x = 0}

        """
        return self.change(oriented=False)

    def _test_unoriented(self, **options):
        r"""
        Verify that :meth:`unoriented` is implemented correctly.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_unoriented()

        """
        tester = self._tester(**options)

        tester.assertEqual(self.unoriented(), self.unoriented().unoriented())

    def intersection(self, other):
        r"""
        Return the intersection with the ``other`` convex set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0)._test_contains()

        """
        tester = self._tester(**options)

        for vertex in self.vertices():
            try:
                tester.assertIn(vertex, self)  # codespell:ignore assertin
            except ValueError:
                # Currently, containment can often not be decided when points
                # do not have coordinates over the base ring.
                tester.assertIsInstance(vertex, HyperbolicPointFromGeodesic)

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
        convex set).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        return (
            self.parent()
            .polygon(
                self.half_spaces(), check=False, assume_sorted=True, assume_minimal=True
            )
            .vertices()
        )

    def is_finite(self):
        r"""
        Return whether all points in this set are finite.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        .. SEEALSO::

            :meth:`is_ideal` and :meth:`is_ultra_ideal` for complementary notions

        """
        return all(
            vertex.is_finite() for vertex in self.vertices(marked_vertices=False)
        )

    def is_ideal(self):
        r"""
        Return whether all points in this set are ideal.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().is_ideal()
            True

            sage: H.vertical(0).is_ideal()
            False

            sage: H.vertical(0).start().is_ideal()
            True

            sage: H(I).is_ideal()
            False

        .. SEEALSO::

            :meth:`is_finite` and :meth:`is_ultra_ideal` for complementary notions

        """
        if self.dimension() >= 1:
            return False

        return all(vertex.is_ideal() for vertex in self.vertices(marked_vertices=False))

    def is_ultra_ideal(self):
        r"""
        Return whether all points in this set are ultra-ideal, i.e., the
        correspond to points outside the Klein disk.

        Note that it is normally not possible to create ultra ideal sets
        (except for the actual empty set). They only exist internally during
        geometric constructions in the Euclidean plane containing the Klein
        disk.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().is_ultra_ideal()
            True

            sage: H.vertical(0).is_ultra_ideal()
            False

        Normally, ultra-ideal objects are not permitted. They can often be created with the ``check=False`` keyword::

            sage: H.point(2, 0, check=False, model="klein").is_ultra_ideal()
            True

            sage: H.geodesic(2, 0, 1, check=False, model="klein").is_ultra_ideal()
            True

        .. SEEALSO::

            :meth:`is_finite` and :meth:`is_ultra_ideal` for complementary notions

        """
        if self.is_empty():
            return True

        if any(
            not vertex.is_ultra_ideal()
            for vertex in self.vertices(marked_vertices=False)
        ):
            return False

        raise NotImplementedError(
            f"{type(self)} does not implement is_ultra_ideal() yet"
        )

    def _test_is_finite(self, **options):
        r"""
        Verify that :meth;`is_finite`, :meth:`is_ideal`, and
        :meth:`is_ultra_ideal` are implemented correctly for this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set()._test_is_finite()

        """
        tester = self._tester(**options)

        finite = self.is_finite()
        ideal = self.is_ideal()
        ultra = self.is_ultra_ideal()

        if self.is_empty():
            tester.assertTrue(finite)
            tester.assertTrue(ideal)
            tester.assertTrue(ultra)
            return

        if finite:
            tester.assertFalse(ideal)
            tester.assertFalse(ultra)

        if ideal:
            tester.assertFalse(finite)
            tester.assertFalse(ultra)

        if ultra:
            tester.assertFalse(finite)
            tester.assertFalse(ideal)

    def change_ring(self, ring):
        r"""
        Return this set as an element of the hyperbolic plane over ``ring``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_change_ring()

        """
        tester = self._tester(**options)
        tester.assertEqual(self, self.change_ring(self.parent().base_ring()))

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this set.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the new set
          will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new set.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness) whether the new set will be explicitly oriented.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: geodesic = H.geodesic(0, 1)

        We can change the base ring over which this set is defined::

            sage: geodesic.change(ring=AA)
            {(x^2 + y^2) - x = 0}

        We can drop the explicit orientation of a set::

            sage: unoriented = geodesic.change(oriented=False)
            sage: unoriented.is_oriented()
            False

        We can also take an unoriented set and pick an orientation::

            sage: oriented = geodesic.change(oriented=True)
            sage: oriented.is_oriented()
            True

        .. SEEALSO::

            :meth:`is_oriented` for oriented an unoriented sets.

        """
        raise NotImplementedError(f"this {type(self)} does not implement change()")

    def _test_change(self, **options):
        r"""
        Verify that the full interface of :meth:`change` has been implemented.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_change()

        """
        tester = self._tester(**options)

        # The ring parameter is supported
        tester.assertEqual(self, self.change(ring=self.parent().base_ring()))

        # The geometry parameter is supported
        tester.assertEqual(self, self.change(geometry=self.parent().geometry))

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

        - ``model`` -- one of ``"half_plane"`` and ``"klein"`` (default: ``"half_plane"``)

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).plot()  # long time (.5s)
            ...Graphics object consisting of 1 graphics primitive

        """
        raise NotImplementedError(f"this {type(self)} does not support plotting")

    def _test_plot(self, **options):
        r"""
        Verify that this set implements :meth:`plot`.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_plot()

        """
        from sage.all import Graphics

        tester = self._tester(**options)
        if self != self.parent().infinity():
            tester.assertIsInstance(self.plot(), Graphics)
            tester.assertIsInstance(self.plot(model="half_plane"), Graphics)
        tester.assertIsInstance(self.plot(model="klein"), Graphics)

    def apply_isometry(self, isometry, model="half_plane", on_right=False):
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
        :meth:`HyperbolicPlane._isometry_gl2_to_sim12` and apply the isometry
        in the Klein model.

        To apply an isometry in the Klein model, we lift objects to the
        hyperboloid model, apply the isometry given by the 3×3 matrix there,
        and then project to the Klein model again.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The horizontal translation by 1 in the upper half plane model, a
        parabolic isometry::

            sage: isometry = matrix([[1, 1], [0, 1]])
            sage: H.vertical(0).apply_isometry(isometry)
            {-x + 1 = 0}

        The same isometry as an isometry of the hyperboloid model::

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: H.vertical(0).apply_isometry(isometry, model="klein")
            {-x + 1 = 0}

        An elliptic isometry::

            sage: isometry = matrix([[1, -1], [1, 1]])
            sage: H.vertical(0).apply_isometry(isometry)
            {(x^2 + y^2) - 1 = 0}

        A hyperbolic isometry::

            sage: isometry = matrix([[1, 0], [0, 1/2]])
            sage: H.vertical(0).apply_isometry(isometry)
            {-x = 0}
            sage: H(I).apply_isometry(isometry)
            2*I

        A reflection::

            sage: isometry = matrix([[-1, 0], [0, 1]])
            sage: H.vertical(0).apply_isometry(isometry)
            {x = 0}

        A glide reflection::

            sage: isometry = matrix([[-1, 0], [0, 1/2]])
            sage: H.vertical(0).apply_isometry(isometry)
            {x = 0}
            sage: H.vertical(1).apply_isometry(isometry)
            {x + 2 = 0}

        An isometry of the upper half plane must have non-zero determinant::

            sage: isometry = matrix([[1, 0], [1, 0]])
            sage: H.vertical(0).apply_isometry(isometry)
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry

        An isometry of the Klein model, must preserve a quadratic form of type
        `(1, 2)`::

            sage: isometry = matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            sage: H.vertical(0).apply_isometry(isometry, model="klein")
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry

        Isometries can be applied to half spaces::

            sage: isometry = matrix([[1, 1], [0, 1]])
            sage: H.vertical(0).left_half_space().apply_isometry(isometry)
            {x - 1 ≤ 0}

        Isometries can be applied to points::

            sage: H(I).apply_isometry(isometry)
            1 + I

        Isometries can be applied to polygons::

            sage: P = H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ], marked_vertices=[I])
            sage: P.apply_isometry(isometry)
            {(x^2 + y^2) - 2*x ≥ 0} ∩ {x - 2 ≤ 0} ∩ {(x^2 + y^2) - 2*x - 3 ≤ 0} ∩ {x ≥ 0} ∪ {1 + I}

        Isometries can be applied to segments::

            sage: segment = H(I).segment(2*I)
            sage: segment.apply_isometry(isometry)
            {-x + 1 = 0} ∩ {2*(x^2 + y^2) - 3*x - 1 ≥ 0} ∩ {(x^2 + y^2) - 3*x - 2 ≤ 0}

        REFERENCES:

        - Svetlana Katok, "Fuchsian Groups", Chicago University Press, Section
          1.3; for the isometries of the upper half plane.

        - James W. Cannon, William J. Floyd, Richard Kenyon, and Walter R.
          Parry, "Hyperbolic Geometry", Flavors of Geometry, MSRI Publications,
          Volume 31, 1997, Section 10; for the isometries as 3×3 matrices.

        """
        if model == "half_plane":
            isometry = self.parent()._isometry_gl2_to_sim12(isometry)
            model = "klein"

        if model == "klein":
            if isometry.dimensions() != (3, 3):
                raise ValueError(
                    "isometry in Klein model must be given as a 3×3 matrix"
                )

            isometry = isometry.change_ring(self.parent().base_ring())

            # Check that the matrix defines an isometry in the hyperboloid
            # model, see CFJK "Hyperbolic Geometry" Theorem 10.1
            from sage.all import diagonal_matrix

            D = (
                isometry.transpose()
                * diagonal_matrix(isometry.parent().base_ring(), [1, 1, -1])
                * isometry
            )

            # These checks should use a specialized predicate of the geometry
            # of this space so they are more robust over inexact rings.
            for i, row in enumerate(D):
                for j, entry in enumerate(row):
                    if (i == j) == self.parent().geometry._zero(entry):
                        raise ValueError("matrix does not define an isometry")

            if not self.parent().geometry._equal(D[0, 0], D[1, 1]):
                raise ValueError("matrix does not define an isometry")

            if not self.parent().geometry._equal(D[0, 0], -D[2, 2]):
                raise ValueError("matrix does not define an isometry")

            return self._apply_isometry_klein(isometry, on_right=on_right)

        raise NotImplementedError("cannot apply isometries in this model yet")

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the result of applying the ``isometry`` to this hyperbolic set.

        Helper methed for :meth:`apply_isometry`. Hyperbolic sets implement
        this method which is not implemented generically.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``); whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: isometry = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: H.empty_set()._apply_isometry_klein(isometry)
            {}

        """
        raise NotImplementedError(
            f"{type(self)} does not implement _apply_isometry_klein() yet"
        )

    def _acted_upon_(self, x, self_on_left):
        r"""
        Return the result of acting upon this set with ``x``.

        INPUT:

        - ``x`` -- an isometry encoded as a 2×2 matrix, see
          :meth:`apply_isometry`

        - ``self_on_left`` -- a boolean; whether this is the right action or
          the left action

        EXAMPLES:

        We apply the Möbius transformation that sends `z` to `(1 + 2z)/(3 +
        4z)`::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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

        if self.dimension() == 0:
            return self in other

        # Make sure that we do not get confused by marked vertices
        self = self.parent().intersection(*self.half_spaces())

        return self.intersection(other) == self

    def _test_is_subset(self, **options):
        r"""
        Verify that :meth:`is_subset` is implemented correctly.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.an_element()._test_is_subset()

        """
        tester = self._tester(**options)

        tester.assertTrue(self.is_subset(self))

        try:
            half_spaces = self.half_spaces()
        except ValueError:
            # We cannot determine half spaces for points whose coordinates are not over the base ring.
            tester.assertIsInstance(self, HyperbolicPointFromGeodesic)
        else:
            tester.assertTrue(self.is_subset(self.parent().intersection(*half_spaces)))

    def _an_element_(self):
        r"""
        Return a point of this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The returned element can be a finite point::

            sage: H.vertical(0).an_element()
            I

        But it can also be an infinite point::

            sage: H(0).segment(I).an_element()
            0

        An exception is raised when there are no elements in this set::

            sage: H.empty_set().an_element()
            Traceback (most recent call last):
            ...
            Exception: empty set has no points

        We get an element for geodesics without end points in the base ring,
        see :meth:`HyperbolicGeodesic._an_element_`::

            sage: H.half_circle(0, 2).an_element()
            (0, 1/3)

        """
        vertices = self.vertices()

        if not vertices:
            raise NotImplementedError(
                "cannot return points of this set yet because it has no vertices"
            )

        return next(iter(vertices))

    @classmethod
    def _enhance_plot(self, plot, model):
        r"""
        Modify the ``plot`` of this set to improve the resulting plot.

        Currently, this adds the unit circle to plots in the Klein disk model.

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H(0)
            sage: plot = plot([])

            sage: type(p)._enhance_plot(plot, model="half_plane")
            ...Graphics object consisting of 0 graphics primitives

        .. jupyter-execute::

            sage: type(p)._enhance_plot(plot, model="klein")
            ...Graphics object consisting of 1 graphics primitive

        """
        if model == "klein":
            from sage.all import circle

            plot = circle([0, 0], 1, fill=False, color="#d1d1d1", zorder=-1) + plot

        return plot

    def is_empty(self):
        r"""
        Return whether this set is empty.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).is_empty()
            False

            sage: H.empty_set().is_empty()
            True

        """
        return self.dimension() < 0

    def __bool__(self):
        r"""
        Return whether this set is non-empty.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: bool(H(0))
            True

            sage: bool(H.empty_set())
            False

        .. SEEALSO:

            :meth:`HyperbolicConvexSet.is_empty`

        """
        return not self.is_empty()

    def dimension(self):
        r"""
        Return the dimension of this set.

        OUTPUT:

        An integer, one of -1, 0, 1, 2.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We treat the empty set as -1-dimensional::

            sage: H.empty_set().dimension()
            -1

        Points are zero-dimensional::

            sage: H(0).dimension()
            0

        Geodesics and segments are one-dimensional::

            sage: H.vertical(0).dimension()
            1

            sage: H(I).segment(2*I).dimension()
            1

        Polygons and half spaces are two dimensional::

            sage: H.random_element("polygon").dimension()
            2

            sage: H.vertical(0).left_half_space().dimension()
            2

        """
        raise NotImplementedError(f"{type(self)} does not implement dimension() yet")

    def _test_dimension(self, **options):
        r"""
        Verify that :meth:`dimension` is implemented correctly.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set()._test_dimension()

        """
        tester = self._tester(**options)

        dimension = self.dimension()

        from sage.all import ZZ

        tester.assertEqual(dimension.parent(), ZZ)

        tester.assertIn(dimension, [-1, 0, 1, 2])  # codespell:ignore assertin

        tester.assertEqual(self.unoriented().dimension(), dimension)
        tester.assertEqual(self.is_point(), dimension == 0)
        tester.assertEqual(self.is_empty(), dimension == -1)

    def is_point(self):
        r"""
        Return whether this set is a single point.

        EXAMPLES:

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).is_point()
            True

            sage: H(I).is_point()
            True

            sage: H.empty_set().is_point()
            False

            sage: H.vertical(0).is_point()
            False

        .. SEEALSO::

            :meth:`is_ideal` to check whether this is a finite or an infinite
            point
        """
        return self.dimension() == 0

    def is_oriented(self):
        r"""
        Return whether this is a set with an explicit orientation.

        Some sets come in two flavors. There are oriented geodesics and
        unoriented geodesics. There are oriented segments and unoriented
        segments.

        This method answers whether a set is of the oriented kind if there is a
        choice.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        Normally, geodesics are oriented::

            sage: g = H.vertical(0)
            sage: g.is_oriented()
            True

            sage: g.start()
            0

            sage: g.end()
            ∞

        We can ask explicitly for an unoriented version::

            sage: h = g.unoriented()
            sage: h.is_oriented()
            False

            sage: h.start()
            Traceback (most recent call last):
            ...
            AttributeError: ... has no attribute 'start'...

        Segments are oriented::

            sage: s = H(I).segment(2*I)
            sage: s.is_oriented()
            True

            sage: s.start()
            I

            sage: s.end()
            2*I

        We can ask explicitly for an unoriented segment::

            sage: u = s.unoriented()
            sage: u.is_oriented()
            False

            sage: u.start()
            Traceback (most recent call last):
            ...
            AttributeError: ... has no attribute 'start'...

        Points are not oriented as there is no choice of orientation::

            sage: H(0).is_oriented()
            False

        Half spaces are not oriented:

            sage: H.vertical(0).left_half_space().is_oriented()
            False

        .. SEEALSO::

            :meth:`change` to pick an orientation on an unoriented set

            :meth:`HyperbolicHalfSpace.__neg__`,
            :meth:`HyperbolicOrientedGeodesic.__neg__`,
            :meth:`HyperbolicOrientedSegment.__neg__` i.e., the ``-`` operator,
            to invert the orientation of a set

        """
        return isinstance(self, HyperbolicOrientedConvexSet)

    def _test_is_oriented(self, **options):
        r"""
        Verify that this set implements :meth:`is_oriented` correctly.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set()._test_is_oriented()

        """
        tester = self._tester(**options)

        if self.is_oriented():
            tester.assertNotEqual(self, -self)

            # Verify that neg inverts the orientation of the set
            tester.assertEqual(self, -(-self))

    def edges(self, as_segments=False, marked_vertices=True):
        r"""
        Return the :class:`segments <HyperbolicOrientedSegment>` and
        :class:`geodesics <HyperbolicOrientedGeodesic>` that bound this set.

        INPUT:

        - ``as_segments`` -- a boolean (default: ``False``); whether to also
          return the geodesics as segments with ideal end points.

        - ``marked_vertices`` -- a boolean (default: ``True``); whether to
          report segments that start or end at redundant marked vertices or
          otherwise whether such marked vertices are completely ignored.

        OUTPUT:

        A set of segments and geodesics. Iteration through this set is in
        counterclockwise order with respect to the points of the set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The edges of a polygon::

            sage: P = H.intersection(
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.vertical(1).left_half_space(),
            ....:   H.half_circle(0, 2).left_half_space())

            sage: P.edges()
            {{-x + 1 = 0} ∩ {2*(x^2 + y^2) - 3*x - 1 ≥ 0}, {x + 1 = 0} ∩ {2*(x^2 + y^2) + 3*x - 1 ≥ 0}, {(x^2 + y^2) - 2 = 0} ∩ {(x^2 + y^2) + 3*x + 1 ≥ 0} ∩ {(x^2 + y^2) - 3*x + 1 ≥ 0}}

        The single edge of a half space:

            sage: H.vertical(0).left_half_space().edges()
            {{-x = 0},}

        A geodesic and a segment are bounded by two edges::

            sage: H.vertical(0).edges()
            {{-x = 0}, {x = 0}}

            sage: H(I).segment(2*I).edges()
            {{-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}, {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}}

        Lower dimensional objects have no edges:

            sage: H(I).edges()
            {}

            sage: H.empty_set().edges()
            {}

        """
        edges = []

        for v, w in self.vertices().pairs():
            edge = v.segment(w)

            if v.is_ideal() and w.is_ideal():
                if edge.right_half_space().is_subset(self):
                    continue

                if as_segments:
                    edge = self.parent().segment(edge, assume_normalized=True)
                    assert isinstance(edge, HyperbolicSegment)

            edges.append(edge)

        return HyperbolicEdges(edges)

    def _test_edges(self, **options):
        r"""
        Verify that this set implements :meth:`edges` correctly.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set()._test_edges()

        """
        tester = self._tester(**options)

        if not self.parent().base_ring().is_exact():
            # Inexact base rings typically do not allow hashing of coordinates
            # which is needed for the set below.
            return

        edges = self.edges()
        endpoints = set(
            [edge.start() for edge in edges] + [edge.end() for edge in edges]
        )

        if self.dimension() > 0:
            vertices = self.vertices()
            tester.assertEqual(endpoints, vertices)

    def area(self, numerical=True):
        r"""
        Return the area of this convex set divided by 2π.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: p = H.polygon([H.geodesic(0, 1).left_half_space(),
            ....:                 H.geodesic(1, Infinity).left_half_space(),
            ....:                 H.geodesic(Infinity, 0).left_half_space()])
            sage: p.area()
            0.5
            sage: p.area(numerical=False)
            1/2

        ::

            sage: a = H.point(0, 1, model='half_plane')
            sage: b = H.point(1, 1, model='half_plane')
            sage: c = H.point(1, 2, model='half_plane')
            sage: d = H.point(0, 2, model='half_plane')
            sage: p = H.polygon([H.geodesic(a, b).left_half_space(),
            ....:                H.geodesic(b, c).left_half_space(),
            ....:                H.geodesic(c, d).left_half_space(),
            ....:                H.geodesic(d, a).left_half_space()])
            sage: p.area()
            0.0696044872730639

        Zero and one-dimensional objects have area zero::

            sage: p = H(I)
            sage: p.area()
            0.0
            sage: p.area(numerical=False)
            0

        ::

            sage: g = H.geodesic(0, 1)
            sage: g.area()
            0.0
            sage: g.area(numerical=False)
            0

        A half space has infinite area::

            sage: h = g.left_half_space()
            sage: h.area()
            inf
            sage: h.area(numerical=False)
            +Infinity

        """
        from sage.all import QQ, Infinity

        if self.dimension() < 2:
            return 0.0 if numerical else QQ.zero()

        edges = [e.geodesic() for e in self.edges()]
        n = len(edges)
        if len(edges) <= 2 or any(
            edges[i].intersection(edges[(i + 1) % n]).is_empty() for i in range(n)
        ):
            return float("inf") if numerical else Infinity

        H = self.parent()
        return QQ((n - 2, 2)) - sum(
            (edges[(i + 1) % n]).angle(-edges[i], numerical=numerical) for i in range(n)
        )

    def __hash__(self):
        r"""
        Return a hash value for this convex set.

        Specific sets should override this method.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
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
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        When determining an isometry that maps sets to each other, we reduce to
        an isometry that maps points or geodesics to each other. Here, we
        produce such more primitive objects that map to each other.

        Sometimes, this mapping is not unique, e.g., when mapping polygons to
        each other, we may rotate the vertices of the polygon. Therefore, this
        returns an iterator that produces the possible mappings of primitive
        objects.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_circle(0, 4).left_half_space()],
            ....:     marked_vertices=[4*I])

            sage: conditions = P._isometry_conditions(P)
            sage: list(conditions)
            [[({x ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {x ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {x ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {x ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {x ≥ 0}),
              (4*I, 4*I)]]

        """
        raise NotImplementedError(
            f"this {type(self)} does not implement _isometry_conditions yet and cannot be used to compute an isometry"
        )

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random convex set.

        Concrete hyperbolic classes should override this method to provide random sets.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` the set should live in

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        Given a hyperbolic object, we can create another one of the same kind::

            sage: p = H(I)
            sage: type(p).random_set(H)  # random output
            7

        ::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexPolygon
            sage: P = HyperbolicConvexPolygon.random_set(H)
            sage: P.dimension()
            2

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        raise NotImplementedError(f"{cls} does not support producing random sets yet")


class HyperbolicOrientedConvexSet(HyperbolicConvexSet):
    r"""
    Base class for sets that have an explicit orientation.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic import HyperbolicOrientedConvexSet
        sage: H = HyperbolicPlane()

        sage: isinstance(H(0), HyperbolicOrientedConvexSet)
        False

        sage: isinstance(H.vertical(0), HyperbolicOrientedConvexSet)
        True

        sage: isinstance(H.vertical(0).unoriented(), HyperbolicOrientedConvexSet)
        False

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.is_oriented`

    """


class HyperbolicConvexFacade(HyperbolicConvexSet, Parent):
    r"""
    A convex subset of the hyperbolic plane that is itself a parent.

    This is the base class for all hyperbolic convex sets that are not points.
    This class solves the problem that we want convex sets to be "elements" of
    the hyperbolic plane but at the same time, we want these sets to live as
    parents in the category framework of SageMath; so they have be a Parent
    with hyperbolic points as their Element class.

    SageMath provides the (not very frequently used and somewhat flaky) facade
    mechanism for such parents. Such sets being a facade, their points can be
    both their elements and the elements of the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()
        sage: v = H.vertical(0)
        sage: p = H(0)
        sage: p in v
        True
        sage: p.parent() is H
        True
        sage: q = v.an_element()
        sage: q
        I
        sage: q.parent() is H
        True

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexFacade
        sage: isinstance(v, HyperbolicConvexFacade)
        True

    """

    def __init__(self, parent, category=None):
        Parent.__init__(self, facade=parent, category=category)

    def parent(self):
        r"""
        Return the hyperbolic plane this is a subset of.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: v = H.vertical(0)
            sage: v.parent()
            Hyperbolic Plane over Rational Field

        """
        return self.facade_for()[0]

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as a point of this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: v = H.vertical(0)
            sage: v(0)
            0
            sage: v(I)
            I
            sage: v(oo)
            ∞
            sage: v(2)
            Traceback (most recent call last):
            ...
            ValueError: point not contained in this set

            sage: 2 in v
            False

        """
        x = self.parent()(x)

        if isinstance(x, HyperbolicPoint):
            if not self.__contains__(x):
                raise ValueError("point not contained in this set")

        return x

    def base_ring(self):
        r"""
        Return the ring over which points of this set are defined.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: v = H.vertical(0)
            sage: v.base_ring()
            Rational Field

        """
        return self.parent().base_ring()


class HyperbolicHalfSpace(HyperbolicConvexFacade):
    r"""
    A closed half space of the hyperbolic plane.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` containing this half space.

    - ``geodesic`` -- the :class:`HyperbolicOrientedGeodesic` to whose left
      this half space lies.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.half_circle(0, 1).left_half_space()
        {(x^2 + y^2) - 1 ≥ 0}

    .. SEEALSO::

        :meth:`HyperbolicPlane.half_space`,
        :meth:`HyperbolicOrientedGeodesic.left_half_space`,
        :meth:`HyperbolicOrientedGeodesic.right_half_space` for the most common
        ways to create half spaces.

    """

    def __init__(self, parent, geodesic):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpace
            sage: H = HyperbolicPlane()

            sage: h = H.half_circle(0, 1).left_half_space()
            sage: isinstance(h, HyperbolicHalfSpace)
            True

            sage: TestSuite(h).run()

        """
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("geodesic must be an oriented geodesic")

        if not geodesic.parent() is parent:
            raise ValueError("geodesic must be in parent")

        self._geodesic = geodesic

    def equation(self, model, normalization=None):
        r"""
        Return an inequality for this half space as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is in the half space if it satisfies `a(x^2 + y^2) + bx + c \ge 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          `a + bx + cy \ge 0`.

        Note that the output is not unique since the coefficients can be scaled
        by a positive scalar.

        INPUT:

        - ``model`` -- either ``"half_plane"`` or ``"klein"``

        - ``normalization`` -- how to normalize the coefficients, see
          :meth:`HyperbolicGeodesic.equation` for details

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpace
            sage: H = HyperbolicPlane()

            sage: h = H.half_circle(0, 1).left_half_space()

            sage: h.equation(model="half_plane")
            (2, 0, -2)

            sage: H.half_space(2, 0, -2, model="half_plane") == h
            True

            sage: h.equation(model="klein")
            (0, 0, 2)

            sage: H.half_space(0, 0, 2, model="klein") == h
            True

        .. SEEALSO::

            :meth:`HyperbolicPlane.half_space` to build a half space from the
            coefficients returned by this method.

        """
        return self._geodesic.equation(model=model, normalization=normalization)

    def _repr_(self):
        r"""
        Return a printable representation of this half space.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.half_circle(0, 1).right_half_space()

            sage: S
            {(x^2 + y^2) - 1 ≤ 0}

            sage: -S
            {(x^2 + y^2) - 1 ≥ 0}

        """
        geodesic = repr(self.boundary())

        cmp = "≥"

        if geodesic.startswith("{-"):
            cmp = "≤"
            geodesic = repr(-self.boundary())

        return geodesic.replace("=", cmp)

    def half_spaces(self):
        r"""
        Return the half spaces defining this half space, i.e., this half space
        itself.

        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.vertical(0).left_half_space()
            sage: [S] == list(S.half_spaces())
            True

        """
        return HyperbolicHalfSpaces([self], assume_sorted=True)

    def __neg__(self):
        r"""
        Return the closure of the complement of this half space.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.half_circle(0, 1).left_half_space()
            sage: S
            {(x^2 + y^2) - 1 ≥ 0}

            sage: -S
            {(x^2 + y^2) - 1 ≤ 0}

        """
        return self._geodesic.right_half_space()

    def boundary(self):
        r"""
        Return a geodesic on the boundary of this half space, oriented such
        that the half space is on its left.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: S = H.vertical(0).left_half_space()
            sage: S.boundary()
            {-x = 0}

        .. SEEALSO::

            :meth:`HyperbolicOrientedGeodesic.left_half_space` to recover the
            half space from the oriented geodesic

        """
        return self._geodesic

    def __contains__(self, point):
        r"""
        Return whether ``point`` is contained in this half space.

        INPUT:

        - ``point`` -- a :class:`HyperbolicPoint`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        We can also check containment of ideal endpoints of geodesics::

            sage: g = H.half_circle(0, 2)
            sage: g.start() in h
            True

            sage: g.end() in h
            False

        ::

            sage: g = H.half_circle(-3, 2)
            sage: g.start() in h
            True

            sage: g.end() in h
            True

        ::

            sage: g = H.half_circle(3, 2)
            sage: g.start() in h
            False

            sage: g.end() in h
            False

        ::

            sage: h = H.half_circle(0, 2).left_half_space()
            sage: g = H.half_circle(0, 5)

            sage: g.start() in h
            True

            sage: g.start() in -h
            False

        .. NOTE::

            The implementation is currently not very robust over inexact rings.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.is_subset` to check containment of
            arbitrary sets.

        """
        point = self.parent()(point)

        if not isinstance(point, HyperbolicPoint):
            raise TypeError("point must be a point in the hyperbolic plane")

        try:
            x, y = point.coordinates(model="klein")
        except ValueError:
            # The point does not have coordinates in the base ring in the Klein model.
            # It is the starting point of a geodesic.
            assert point.is_ideal()

            if not isinstance(point, HyperbolicPointFromGeodesic):
                raise NotImplementedError(
                    "cannot decide whether this ideal point is contained in the half space yet"
                )

            boundary = self.boundary()
            intersection = boundary._intersection(point._geodesic)

            if intersection is None:
                # The boundary of the half space and the geodesic defining the
                # point do not intersect in a single point (not even in an
                # ultra-ideal point,) i.e., they are parallel in the Klein model.
                return point._geodesic.an_element() in self

            if not intersection.is_finite():
                # The intersection point between the geodesic defining the
                # point and the half space boundary is ideal or ultra-ideal.
                # Either the entire geodesic is in the half space or none of
                # it.
                return point._geodesic.an_element() in self

            # The intersection point between the geodesic defining the point
            # and the half space boundary is finite, i.e., inside the unit
            # circle in the Klein model. The segment between the starting point
            # of the geodesic and that intersection point is either completely
            # inside the half space or completely outside.
            ccw = boundary.ccw(point._geodesic)
            assert ccw != 0
            return ccw < 0

        a, b, c = self.equation(model="klein")

        # We should use a specialized predicate here to do something more
        # reasonable for points that are close to the boundary over inexact
        # rings.
        return self.parent().geometry._sgn(a + b * x + c * y) >= 0

    def __eq__(self, other):
        r"""
        Return whether this set is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: h = H.vertical(0).left_half_space()

            sage: h == H.vertical(0).left_half_space()
            True
            sage: h == H.vertical(0).right_half_space()
            False

        ::

            sage: h != H.vertical(0).left_half_space()
            False
            sage: h != H.vertical(0).right_half_space()
            True

        """
        if not isinstance(other, HyperbolicHalfSpace):
            return False
        return self._geodesic == other._geodesic

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this half space in the hyperbolic ``model``.

        See :meth:`HyperbolicConvexPolygon.plot` for the supported keyword
        arguments.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"`` (default: ``"half_plane"``)

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        r"""
        Return a modified copy of this half space.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the new half
          space will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new half space.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since half spaces cannot
          have an explicit orientation. See :meth:`~HyperbolicConvexSet.is_oriented`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: h = H.vertical(0).left_half_space()

        We change the base ring over which this space is defined::

            sage: h.change(ring=AA)
            {x ≤ 0}

        We cannot change the orientedness of a half space:

            sage: h.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: half spaces cannot have an explicit orientation

            sage: h.change(oriented=False)
            {x ≤ 0}

        """
        if ring is not None or geometry is not None:
            self = self._geodesic.change(ring=ring, geometry=geometry).left_half_space()

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("half spaces cannot have an explicit orientation")

        return self

    def dimension(self):
        r"""
        Return the dimension of this half space, i.e., 2.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).left_half_space().dimension()
            2

        """
        from sage.all import ZZ

        return ZZ(2)

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices bounding this half space.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          half space cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        Use :meth:`HyperbolicOrientedGeodesic.start` and
        :meth:`HyperbolicOrientedGeodesic.end` on the :meth:`boundary` to get
        the end points in an order consistent with the orientation::

            sage: g = h.boundary()
            sage: g.start(), g.end()
            (0, ∞)

            sage: g = (-h).boundary()
            sage: g.start(), g.end()
            (∞, 0)

        Vertices can be computed even though they do not have coordinates over
        the :meth:`HyperbolicPlane.base_ring`::

            sage: H.half_circle(0, 2).left_half_space().vertices()
            {-1.41421356237310, 1.41421356237310}

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return self.boundary().vertices()

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the image of this half space under ``isometry``.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: H.vertical(0).left_half_space()._apply_isometry_klein(isometry)
            {x - 1 ≤ 0}

        """
        return self._geodesic.apply_isometry(
            isometry, model="klein", on_right=on_right
        ).left_half_space()

    def __hash__(self):
        r"""
        Return a hash value for this half space

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: v = H.vertical(0).left_half_space()
            sage: w = H.vertical(1).right_half_space()

            sage: conditions = v._isometry_conditions(w)
            sage: list(conditions)
            [[({-x = 0}, {x - 1 = 0})]]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet._isometry_conditions` for a general description.

        r"""
        yield [(self.boundary(), other.boundary())]

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random half space.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the half plane

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpace
            sage: x = HyperbolicHalfSpace.random_set(H)

            sage: x.dimension()
            2

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        return HyperbolicOrientedGeodesic.random_set(parent).left_half_space()


class HyperbolicGeodesic(HyperbolicConvexFacade):
    r"""
    A geodesic in the hyperbolic plane.

    This is the abstract base class of :class:`HyperbolicUnorientedGeodesic`
    and :class:`HyperbolicOrientedGeodesic`.

    ALGORITHM:

    Internally, we represent geodesics as a triple `a, b, c` such that they satisfy the equation

    .. MATH::

        a + bx + cy = 0

    in the Klein disk model.

    Note that due to this representation we can always compute intersection
    points of geodesics but we cannot always get the coordinates of the ideal
    end points of a geodesic (since we would have to take a square root to
    solve for the points on the unit circle).

    It might be beneficial to store geodesics differently, see
    https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/hyperbolic.20geometry/near/284722650
    for a discussion.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` this geodesic lives in

    - ``a`` -- an element of :meth:`HyperbolicPlane.base_ring`

    - ``b`` -- an element of :meth:`HyperbolicPlane.base_ring`

    - ``c`` -- an element of :meth:`HyperbolicPlane.base_ring`

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.geodesic(1, 2, 3, model="klein")
        {2*(x^2 + y^2) + 2*x - 1 = 0}

    .. SEEALSO::

        :meth:`HyperbolicPlane.geodesic` for various ways of constructing geodesics

    """

    def __init__(self, parent, a, b, c):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: geodesic = H.vertical(0)

            sage: isinstance(geodesic, HyperbolicGeodesic)
            True

            sage: TestSuite(geodesic).run()

            sage: isinstance(geodesic.unoriented(), HyperbolicGeodesic)
            True

            sage: TestSuite(geodesic.unoriented()).run()

        """
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
        r"""
        Validate the equation defining this geodesic.

        Implements :meth:`HyperbolicConvexSet._check`.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); ignored

        EXAMPLES::


            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: geodesic = H.vertical(0)
            sage: geodesic._check()

            sage: geodesic = H.geodesic(2, 0, 1, model="klein", check=False)
            sage: geodesic._check()
            Traceback (most recent call last):
            ...
            ValueError: equation 2 + (0)*x + (1)*y = 0 does not define a chord in the Klein model

        .. SEEALSO::

            :meth:`is_ultra_ideal` to check whether a chord is completely outside the Klein disk

            :meth:`is_ideal` to check whether a chord touches the Klein disk

        """
        if self.is_ultra_ideal() or self.is_ideal():
            raise ValueError(
                f"equation {self._a} + ({self._b})*x + ({self._c})*y = 0 does not define a chord in the Klein model"
            )

    def is_ideal(self):
        r"""
        Return whether all hyperbolic points of this geodesic are ideal, i.e.,
        the defining equation of this geodesic in the Klein model only touches
        the Klein disk but does not intersect it.

        Note that it is normally not possible to create ideal geodesics. They
        only exist internally during constructions in the Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).is_ideal()
            False

            sage: geodesic = H.geodesic(1, 0, 1, model="klein", check=False)
            sage: geodesic.is_ideal()
            True

            sage: geodesic = H.geodesic(2, 0, 1, model="klein", check=False)
            sage: geodesic.is_ideal()
            False

        .. NOTE::

            The implementation of this predicate is not numerically robust over inexact rings.

        .. SEEALSO::

            :meth:`is_ultra_ideal` to detect whether a geodesic does not even
            touch the Klein disk

        """
        # We should probably use a specialized predicate of the geometry to
        # make this more robust over inexact rings.
        return self.parent().geometry._equal(
            self._b * self._b + self._c * self._c, self._a * self._a
        )

    def is_ultra_ideal(self):
        r"""
        Return whether the line given by the defining equation is completely
        outside the Klein disk, i.e., all "points" of this geodesic are ultra-ideal.

        Note that it is normally not possible to create ultra-ideal geodesics.
        They only exist internally during constructions in the Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).is_ultra_ideal()
            False

            sage: geodesic = H.geodesic(1, 0, 1, model="klein", check=False)
            sage: geodesic.is_ultra_ideal()
            False

            sage: geodesic = H.geodesic(2, 0, 1, model="klein", check=False)
            sage: geodesic.is_ultra_ideal()
            True

        .. NOTE::

            The implementation of this predicate is not numerically robust over inexact rings.

        .. SEEALSO::

            :meth:`is_ideal` to detect whether a geodesic touches the Klein disk

        """
        # We should probably use a specialized predicate of the geometry to
        # make this more robust over inexact rings.
        return (
            self.parent().geometry._cmp(
                self._b * self._b + self._c * self._c, self._a * self._a
            )
            < 0
        )

    def _repr_(self, model=None):
        r"""
        Return a printable representation of this geodesic.

        INPUT:

        - ``model`` -- ``"half_plane"`` or ``"klein"`` (default: ``None`` to
          use ``"half_plane"`` if possible); in which model of hyperbolic
          geometry the equation is realized

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        If the geodesic is oriented, we print a defining equation of the
        geodesic such that when replacing the ``=`` with a ``≥``, we get an
        equation for the half space to the left of the geodesic::

            sage: H.vertical(1)
            {-x + 1 = 0}
            sage: H.vertical(1).left_half_space()
            {x - 1 ≤ 0}

        Normally, geodesics are shown with their equation in the upper half
        plane model. We can also ask for an equation of the chord in the Klein
        disk::

            sage: H.vertical(1)._repr_(model="klein")
            '{1 + -x - y = 0}'

        This representation is also chosen automatically for objects that do
        not have a representation in upper half plane such as ultra-ideal
        geodesics::

            sage: H.geodesic(2, 0, 1, model="klein", check=False)
            {2 + y = 0}

        .. SEEALSO::

            :meth:`equation` to access the coefficients of the equation
            defining the geodesic

        """
        if model is None:
            model = (
                "klein" if self.is_ultra_ideal() or self.is_ideal() else "half_plane"
            )

        if model == "half_plane":
            # Convert to the upper half plane model as a(x^2 + y^2) + bx + c = 0.
            a, b, c = self.equation(model="half_plane", normalization=["gcd", None])

            from sage.all import PolynomialRing

            R = PolynomialRing(self.parent().base_ring(), names="x")
            if self.parent().geometry._sgn(a) != 0:
                return (
                    f"{{{repr(R([0, a]))[:-1]}(x^2 + y^2){repr(R([c, b, 1]))[3:]} = 0}}"
                )
            else:
                return f"{{{repr(R([c, b]))} = 0}}"

        if model == "klein":
            a, b, c = self.equation(model="klein", normalization=["gcd", None])

            from sage.all import PolynomialRing

            R = PolynomialRing(self.parent().base_ring(), names=["x", "y"])
            polynomial_part = R({(1, 0): b, (0, 1): c})
            if self.parent().geometry._sgn(a) != 0:
                return f"{{{repr(a)} + {repr(polynomial_part)} = 0}}"
            else:
                return f"{{{repr(polynomial_part)} = 0}}"

        raise NotImplementedError("printing not supported in this model")

    def equation(self, model, normalization=None):
        r"""
        Return an equation for this geodesic as a triple ``a``, ``b``, ``c`` such that:

        - if ``model`` is ``"half_plane"``, a point `x + iy` of the upper half
          plane is on the geodesic if it satisfies `a(x^2 + y^2) + bx + c = 0`.

        - if ``model`` is ``"klein"``, points `(x, y)` in the unit disk satisfy
          are on the geodesic if `a + bx + cy = 0`.

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

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: v = H.vertical(0)
            sage: v.equation(model="half_plane")
            (0, -2, 0)

            sage: v.equation(model="half_plane", normalization="gcd")
            (0, -1, 0)

            sage: v.equation(model="klein")
            (0, -1, 0)

        Sometimes, the desired normalization might not be possible (a more
        realistic example would be exact-real coefficients)::

            sage: H = HyperbolicPlane(ZZ)
            sage: g = H.geodesic(2, 3, -4, model="half_plane")

            sage: g.equation(model="half_plane", normalization="one")
            Traceback (most recent call last):
            ...
            TypeError: ...

        In this case, we can ask for the best of several normalization::

            sage: g.equation(model="half_plane", normalization=["one", "gcd", None])
            (2, 3, -4)

        For ultra-ideal geodesics, the equation in the half plane model is not
        very useful::

            sage: g = H.geodesic(2, 0, 1, model="klein", check=False)
            sage: g.equation(model="half_plane")  # i.e., 3*(x^2 + y^2) + 1 = 0
            (3, 0, 1)

        .. SEEALSO::

            :meth:`HyperbolicPlane.geodesic` to create a geodesic from an
            equation

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

        sgn = self.parent().geometry._sgn
        sgn = (
            -1
            if (
                sgn(a) < 0
                or (sgn(a) == 0 and b < 0)
                or (sgn(a) == 0 and sgn(b) == 0 and sgn(c) < 0)
            )
            else 1
        )

        while normalization:
            strategy = normalization.pop()

            try:
                a, b, c = HyperbolicGeodesic._normalize_coefficients(
                    a, b, c, strategy=strategy
                )
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
    def _normalize_coefficients(cls, a, b, c, strategy):
        r"""
        Return the normalized coefficients for the equation of a geodesic.

        Helper method for :meth:`equation`.

        INPUT:

        - ``a``, ``b``, ``c`` -- the coefficients of the geodesic equation

        - ``strategy`` -- one of ``"gcd"`` or ``"one"``

        EXAMPLES::

            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic

        Normalize the leading coefficient to be ±1::

            sage: HyperbolicGeodesic._normalize_coefficients(QQ(6), -QQ(10), -QQ(12), strategy="one")
            (1, -5/3, -2)
            sage: HyperbolicGeodesic._normalize_coefficients(-QQ(6), QQ(10), QQ(12), strategy="one")
            (-1, 5/3, 2)

        Divide the coefficients by their GCD::

            sage: HyperbolicGeodesic._normalize_coefficients(QQ(6), -QQ(10), -QQ(12), strategy="gcd")
            (3, -5, -6)
            sage: HyperbolicGeodesic._normalize_coefficients(-QQ(6), QQ(10), QQ(12), strategy="gcd")
            (-3, 5, 6)

        """
        R = a.parent()

        if strategy is None:
            return a, b, c

        if strategy == "gcd":

            def gcd(*coefficients):
                coefficients = [c for c in coefficients if c]
                if len(coefficients) == 1:
                    return coefficients[0]

                from sage.all import gcd

                return gcd(coefficients)

            d = gcd(a, b, c)
            if d < 0:
                d *= -1
            return R(a / d), R(b / d), R(c / d)

        if strategy == "one":
            if a:
                d = a
            elif b:
                d = b
            else:
                assert c
                d = c

            if d < 0:
                d *= -1

            return R(a / d), R(b / d), R(c / d)

        raise ValueError(f"unknown normalization {strategy}")

    def half_spaces(self):
        r"""
        Return the two half spaces whose intersection this geodesic is.

        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"`` (default: ``"half_plane"``)

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).plot()
            ...Graphics object consisting of 1 graphics primitive

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

            sage: from flatsurf import HyperbolicPlane
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
            raise NotImplementedError(
                "can only compute pole if geodesic is a not a diameter in the Klein model"
            )

        def tangent(endpoint):
            x, y = endpoint.coordinates(model="klein")
            return self.parent().geodesic(
                -(x * x + y * y), x, y, model="klein", check=False
            )

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

            sage: from flatsurf import HyperbolicPlane
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

        .. NOTE::

            Currently, the orientation of the returned geodesic is somewhat
            random. It should probably be counterclockwise to this geodesic.

        """
        if point_or_geodesic is None:
            point_or_geodesic = self.an_element()

        point_or_geodesic = self.parent()(point_or_geodesic)

        if isinstance(point_or_geodesic, HyperbolicGeodesic):
            other = point_or_geodesic
            if self.unoriented() == other.unoriented():
                return self.perpendicular(self.an_element())

            if self.is_diameter() and other.is_diameter():
                raise ValueError(
                    f"no geodesic perpendicular to both {self} and {other}"
                )

            if other.is_diameter():
                return other.perpendicular(self)

            if self.is_diameter():
                # Construct the line a + bx + cy = 0 perpendicular to the
                # diameter through the pole of the other geodesic.
                b, c = (-self._c, self._b)
                x, y = other.pole().coordinates(model="klein")
                a = -(b * x + c * y)

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
                a = -(b * x + c * y)

                perpendicular = self.parent().geodesic(
                    a, b, c, model="klein", oriented=False
                )
            else:
                perpendicular = self.parent().geodesic(
                    self.pole(), point, oriented=False
                )

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
        computation). All other geodesics, are just translated versions of
        these so we can just conjugate with a translation to determine the
        fixed point, i.e., the fixed point is a translate of one of the above.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
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

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0).is_vertical()
            True
            sage: H.half_circle(0, 1).is_vertical()
            False

        """
        return self.parent().infinity() in self

    def __eq__(self, other):
        r"""
        Return whether this geodesic is identical to other up to (orientation
        preserving) scaling of the defining equation.


        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.vertical(0) == H.vertical(0)
            True

        We distinguish oriented and unoriented geodesics::

            sage: H.vertical(0).unoriented() == H.vertical(0)
            False
            sage: H.vertical(0).unoriented() != H.vertical(0)
            True

        We distinguish differently oriented geodesics::

            sage: H.vertical(0) == -H.vertical(0)
            False
            sage: H.vertical(0) != -H.vertical(0)
            True

        We do, however, identify geodesics whose defining equations differ by some scaling::

            sage: g = H.vertical(0)
            sage: g.equation(model="half_plane")
            (0, -2, 0)
            sage: h = H.geodesic(0, -4, 0, model="half_plane")
            sage: g.equation(model="half_plane") == h.equation(model="half_plane")
            False
            sage: g == h
            True
            sage: g != h
            False

        .. NOTE::

            Over inexact rings, this method is not very reliable. To some
            extent this is inherent to the problem but also the implementation
            uses generic predicates instead of relying on a specialized
            implementation in the :class:`HyperbolicGeometry`.

        """
        if type(self) is not type(other):
            return False

        other = self.parent()(other)

        # See note in the docstring. We should use specialized geometry here.
        equal = self.parent().geometry._equal
        sgn = self.parent().geometry._sgn

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

    def __contains__(self, point):
        r"""
        Return whether ``point`` lies on this geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
            # Shortcut the most common case (that intersection cannot handle).
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
        r"""
        Return the dimension of this set, i.e., 1.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).dimension()
            1

        """
        from sage.all import ZZ

        return ZZ(1)

    def change(self, *, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this geodesic.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the new geodesic will be
          defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new geodesic.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); whether the new geodesic should be oriented.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

        The base ring over which this geodesic is defined can be changed::

            sage: H.vertical(1).change_ring(QQ)
            {-x + 1 = 0}

        But we cannot change the base ring if the geodesic cannot be expressed
        in the smaller ring::

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
            self = (
                self.parent()
                .change_ring(ring, geometry=geometry)
                .geodesic(
                    self._a,
                    self._b,
                    self._c,
                    model="klein",
                    check=False,
                    oriented=self.is_oriented(),
                )
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            self = self.parent().geodesic(
                self._a, self._b, self._c, model="klein", check=False, oriented=oriented
            )

        return self

    def geodesic(self):
        r"""
        Return the geodesic underlying this set, i.e., this geodesic itself.

        This method exists to unify the interface between segments and
        geodesics, see :meth:`HyperbolicSegment.geodesic`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).geodesic() == H.vertical(0)
            True

        """
        return self

    def __hash__(self):
        r"""
        Return a hash value for this geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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

        return hash(
            (type(self), self.equation(model="klein", normalization=["one", "gcd"]))
        )

    def _intersection(self, other):
        r"""
        Return the intersection of this geodesic and ``other``.

        Return ``None`` if they do not intersect in a point.

        INPUT:

        - ``other`` -- another :class:`HyperbolicGeodesic`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

        Geodesics can intersect in a finite point::

            sage: H.vertical(0)._intersection(H.half_circle(0, 1))
            I

        Geodesics can intersect in an ideal point::

            sage: H.vertical(1)._intersection(H.half_circle(0, 1))
            1

        Geodesics might intersect in an ultra ideal point::

            sage: H.half_circle(0, 1)._intersection(H.half_circle(1, 8))
            (-3, 0)

        Or they are parallel in the Klein model::

            sage: H.half_circle(0, 1)._intersection(H.half_circle(0, 4))

        Note that geodesics that overlap do not intersect in a point::

            sage: H.vertical(0)._intersection(H.vertical(0))

            sage: H.vertical(0)._intersection(-H.vertical(0))

        TESTS::

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

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.intersection` for intersection with more
            general sets.

            :meth:`HyperbolicPlane.intersection` for the generic
            intersection of convex sets.

        """
        if not isinstance(other, HyperbolicGeodesic):
            return super().intersection(other)

        xy = self.parent().geometry.intersection(
            (self._a, self._b, self._c), (other._a, other._b, other._c)
        )

        if xy is None:
            return None

        return self.parent().point(*xy, model="klein", check=False)

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the result of applying the ``isometry`` to this geodesic.

        Helper methed for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``); whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We apply a reflection to a geodesic::

            sage: isometry = matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: g = H.geodesic(-1, 1)

            sage: g._apply_isometry_klein(isometry)
            {(x^2 + y^2) - 1 = 0}

        Note how the reflection swaps the end points of the geodesic::

            sage: g.start().apply_isometry(isometry, model="klein")
            1
            sage: g.end().apply_isometry(isometry, model="klein")
            -1

        However, the isometry maps the oriented geodesic to itself since what's
        left of the geodesic is not changed::

            sage: g._apply_isometry_klein(isometry) == g
            True

        An isometry that changes the orientation of the geodesic::

            sage: isometry = matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
            sage: g._apply_isometry_klein(isometry)
            {-(x^2 + y^2) + 1 = 0}

        For an unoriented geodesic, the geodesic is unchanged though::

            sage: g.unoriented()._apply_isometry_klein(isometry) == g.unoriented()
            True

        TESTS::

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
            vector(self.parent().base_ring(), [self._b, self._c, self._a]) * isometry
        )
        return self.parent().geodesic(
            a, b, c, model="klein", oriented=self.is_oriented()
        )

    def _isometry_equations(self, isometry, image, λ):
        r"""
        Return equations that must be satisfied if this set converts to
        ``image`` under ``isometry`` using ``λ`` as a free variable for the
        scaling factor.

        Helper method for :meth:`HyperbolicPlane._isometry_from_primitives`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix describing a (right) isometry; typically
          not over the base ring but in symbolic variables

        - ``image`` -- a geodesic

        - ``λ`` -- a symbolic variable

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: R.<a, b, c, d, λ> = QQ[]
            sage: isometry = H._isometry_gl2_to_sim12(matrix([[a, b], [c, d]]))
            sage: isometry
            [                            b*c + a*d                             a*c - b*d                             a*c + b*d]
            [                            a*b - c*d 1/2*a^2 - 1/2*b^2 - 1/2*c^2 + 1/2*d^2 1/2*a^2 + 1/2*b^2 - 1/2*c^2 - 1/2*d^2]
            [                            a*b + c*d 1/2*a^2 - 1/2*b^2 + 1/2*c^2 - 1/2*d^2 1/2*a^2 + 1/2*b^2 + 1/2*c^2 + 1/2*d^2]

            sage: H.vertical(0)._isometry_equations(isometry, H.vertical(1), λ)
            [-b*c - a*d + λ, -a*c + b*d + λ, -a*c - b*d - λ]

        """
        R = λ.parent()

        a, b, c = self.equation(model="klein")
        fa, fb, fc = image.equation(model="klein")
        from sage.all import vector

        condition = vector((b, c, a)) * isometry - λ * vector(R, (fb, fc, fa))
        return condition.list()

    def _an_element_(self):
        r"""
        Return a finite point on this geodesic.

        ALGORITHM:

        We take the chord in the Klein model that intersects this geodesic
        perpendicularly and passes through the origin. The point of
        intersection must be a finite point.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).an_element()
            I

            sage: H.half_circle(0, 2).an_element()
            (0, 1/3)

        """
        other = self.parent().geodesic(0, -self._c, self._b, model="klein", check=False)
        cross = other._intersection(self)
        assert cross
        return cross

    def vertices(self, marked_vertices=True):
        r"""
        Return the ideal end points of this geodesic.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          geodesic cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        Use :meth:`HyperbolicOrientedGeodesic.start` and
        :meth:`HyperbolicOrientedGeodesic.end` to get the end points in an
        order that is consistent with orientation::

            sage: v.start(), v.end()
            (0, ∞)

            sage: (-v).start(), (-v).end()
            (∞, 0)

        The vertices can also be determined for an unoriented geodesic::

            sage: v.unoriented().vertices()
            {0, ∞}

        Vertices can be computed even if they do not have coordinates over the
        :meth:`HyperbolicPlane.base_ring`::

            sage: H.half_circle(0, 2).vertices()
            {-1.41421356237310, 1.41421356237310}

        .. NOTE::

            The implementation of this method is not robust over inexact rings.

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        if not self.is_oriented():
            self = self.change(oriented=True)

        return HyperbolicVertices([self.start(), self.end()])


class HyperbolicUnorientedGeodesic(HyperbolicGeodesic):
    r"""
    An unoriented geodesic in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic import HyperbolicUnorientedGeodesic
        sage: H = HyperbolicPlane()

        sage: H.vertical(0).unoriented()
        {x = 0}

    TESTS::

        sage: v = H.vertical(0).unoriented()
        sage: isinstance(v, HyperbolicUnorientedGeodesic)
        True

        sage: TestSuite(v).run()

    .. SEEALSO::

        :meth:`HyperbolicPlane.geodesic` to create geodesics from points or equations

    """

    def _isometry_conditions(self, other):
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: v = H.vertical(0).unoriented()
            sage: w = H.vertical(1).unoriented()

            sage: conditions = v._isometry_conditions(w)
            sage: list(conditions)
            [[({-x = 0}, {-x + 1 = 0})], [({-x = 0}, {x - 1 = 0})]]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet._isometry_conditions` for a general description.

        r"""
        self = self.change(oriented=True)
        other = other.change(oriented=True)
        yield [(self, other)]
        yield [(self, -other)]

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random unoriented geodesic.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the geodesic

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicUnorientedGeodesic
            sage: x = HyperbolicUnorientedGeodesic.random_set(H)

            sage: x.dimension()
            1

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        return HyperbolicOrientedGeodesic.random_set(parent).unoriented()


class HyperbolicOrientedGeodesic(HyperbolicGeodesic, HyperbolicOrientedConvexSet):
    r"""
    An oriented geodesic in the hyperbolic plane.

    Internally, we represent geodesics as the chords satisfying the equation

    .. MATH::

        a + bx + cy = 0

    in the unit disk of the Klein model.

    The geodesic is oriented such that the half space

    .. MATH::

        a + bx + cy ≥ 0

    is on its left.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` this geodesic lives in

    - ``a`` -- an element of :meth:`HyperbolicPlane.base_ring`

    - ``b`` -- an element of :meth:`HyperbolicPlane.base_ring`

    - ``c`` -- an element of :meth:`HyperbolicPlane.base_ring`

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.vertical(0)
        {-x = 0}

        sage: H.half_circle(0, 1)
        {(x^2 + y^2) - 1 = 0}

        sage: H.geodesic(H(I), 0)
        {x = 0}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicOrientedGeodesic

        sage: g = H.vertical(0)

        sage: isinstance(g, HyperbolicOrientedGeodesic)
        True

        sage: TestSuite(g).run()

    .. SEEALSO::

        :meth:`HyperbolicPlane.geodesic` for the most common ways to construct
        geodesics.

        :class:`HyperbolicUnorientedGeodesic` for geodesics without an explicit
        orientation and :class:`HyperbolicGeodesic` for shared functionality of
        all geodesics.

    """

    def __neg__(self):
        r"""
        Return this geodesic with its orientation reversed.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: -H.vertical(0)
            {x = 0}

        """
        return self.parent().geodesic(
            -self._a, -self._b, -self._c, model="klein", check=False
        )

    def start(self):
        r"""
        Return the ideal starting point of this geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        r"""
        Return the ideal end point of this geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        r"""
        Return the closed half space to the left of this (oriented) geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(0).left_half_space()
            {x ≤ 0}

        .. SEEALSO::

            :meth:`right_half_space` for the corresponding half space on the other side

            :meth:`HyperbolicPlane.half_space` for another method to create half spaces.

        """
        return self.parent().half_space(self._a, self._b, self._c, model="klein")

    def right_half_space(self):
        r"""
        Return the closed half space to the right of this (oriented) geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(AA)

            sage: H.vertical(0).right_half_space()
            {x ≥ 0}

        .. SEEALSO::

            :meth:`left_half_space` for the corresponding half space on the other side

            :meth:`HyperbolicPlane.half_space` for another method to create half spaces.

        """
        return (-self).left_half_space()

    def _configuration(self, other):
        r"""
        Return a classification of the (Euclidean) angle between this geodesic
        and ``other`` in the Klein model.

        This is a helper method for :meth:`HyperbolicConvexPolygon._normalize`.

        INPUT:

        - ``other`` -- another :class:`HyperbolicOrientedGeodesic`

        OUTPUT:

        A string explaining how the two geodesics are oriented.

        .. NOTE::

            This check is not robust over inexact rings and should be improved
            for that use case.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        Two geodesics can be equal::

            sage: H.vertical(0)._configuration(H.vertical(0))
            'equal'

        They can be equal with reversed orientation::

            sage: H.vertical(0)._configuration(-H.vertical(0))
            'negative'

        They can be parallel in the Klein model::

            sage: H.vertical(0)._configuration(H.geodesic(1/2, 2))
            'parallel'

        They can be parallel but with reversed orientation::

            sage: H.vertical(0)._configuration(H.geodesic(2, 1/2))
            'anti-parallel'

        Or they can intersect. We can distinguish the case that ``other``
        crosses over from left-to-right or from right-to-left::

            sage: H.vertical(0)._configuration(H.geodesic(1/3, 2))
            'concave'

            sage: H.vertical(0)._configuration(H.geodesic(1/2, 3))
            'convex'

        .. SEEALSO::

            :meth:`~HyperbolicGeodesic._intersection` to compute the (ultra-ideal) intersection of geodesics

        """
        if not isinstance(other, HyperbolicOrientedGeodesic):
            raise TypeError("other must be an oriented geodesic")

        intersection = self._intersection(other)

        if intersection is None:
            # We should use a specialized method of geometry here to make this
            # more robust over inexact rings.
            orientation = self.parent().geometry._sgn(
                self._b * other._b + self._c * other._c
            )

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

        # Probably convex and concave are not the best terms here.
        if orientation > 0:
            return "convex"

        # We probably would not need to consider the concave case if we always
        # packed all geodesics into a bounding box that contains the unit disk.
        return "concave"

    def parametrize(self, point, model, check=True):
        r"""
        Return the value of ``point`` in a linear parametrization of this
        geodesic.

        INPUT:

        - ``point`` -- a :class:`HyperbolicPoint` on this geodesic

        - ``model`` -- a string; currently only ``"euclidean"`` is supported

        - ``check`` -- a boolean (default: ``True``); whether to ensure that
          ``point`` is actually a point on the geodesic

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We can parametrize points on a geodesic such that the order of the
        points corresponds to the parameters we get::

            sage: g = H.vertical(0)

            sage: g.parametrize(0, model="euclidean")
            -1

            sage: g.parametrize(I, model="euclidean")
            0

            sage: g.parametrize(2*I, model="euclidean")
            3/5

            sage: g.parametrize(oo, model="euclidean")
            1

        .. NOTE::

            This method is not robust for points over inexact rings and should
            be improved.

        .. SEEALSO::

            :meth:`unparametrize` for the recovers the point from the parameter

        """
        if check:
            point = self.parent()(point)
            if point not in self:
                raise ValueError("point must be on geodesic to be parametrized")

        if model == "euclidean":
            base = self.an_element().coordinates(model="klein")
            tangent = (self._c, -self._b)

            # We should use a specialized predicate here to make this work
            # better over inexact rings.
            coordinate = 0 if not self.parent().geometry._zero(tangent[0]) else 1
            return (
                point.coordinates(model="klein")[coordinate] - base[coordinate]
            ) / tangent[coordinate]

        raise NotImplementedError("cannot parametrize a geodesic over this model yet")

    def unparametrize(self, λ, model, check=True):
        r"""
        Return the point parametrized by ``λ`` on this geodesic.

        INPUT:

        - ``λ`` -- an element of :meth:`HyperbolicPlane.base_ring`

        - ``model`` -- a string; currently only ``"euclidean"`` is supported

        - ``check`` -- a boolean (default: ``True``); whether to ensure that
          the returned point is not ultra ideal

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        This method is the inverse of :meth;`parametrize`::

            sage: g = H.vertical(0)

            sage: g.unparametrize(g.parametrize(0, model="euclidean"), model="euclidean")
            0

            sage: g.unparametrize(g.parametrize(I, model="euclidean"), model="euclidean")
            I

            sage: g.unparametrize(g.parametrize(2*I, model="euclidean"), model="euclidean")
            2*I

            sage: g.unparametrize(g.parametrize(oo, model="euclidean"), model="euclidean")
            ∞

        """
        if model == "euclidean":
            base = self.an_element().coordinates(model="klein")
            tangent = (self._c, -self._b)

            λ = self.parent().base_ring()(λ)

            return self.parent().point(
                x=base[0] + λ * tangent[0],
                y=base[1] + λ * tangent[1],
                model="klein",
                check=check,
            )

        raise NotImplementedError("cannot parametrize a geodesic over this model yet")

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random oriented geodesic.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the geodesic

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicOrientedGeodesic
            sage: x = HyperbolicOrientedGeodesic.random_set(H)

            sage: x.dimension()
            1

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        a = HyperbolicPointFromCoordinates.random_set(parent)
        b = a
        while b == a:
            b = HyperbolicPointFromCoordinates.random_set(parent)

        return parent.geodesic(a, b)

    def ccw(self, other, euclidean=False):
        r"""
        Return +1, -1, or zero if the turn from this geodesic to ``other`` is
        counterclockwise, clockwise, or if the geodesics are parallel at their
        point of intersection, respectively.

        If this geodesic and ``other`` do not intersect, a ``ValueError`` is
        raised.

        INPUT:

        - ``other`` -- another geodesic intersecting this geodesic

        - ``euclidean`` -- a boolean (default: ``False``); if ``True``, the
          geodesics are treated as lines in the Euclidean plane coming from
          their representations in the Klein model, i.e., geodesics
          intersecting at a non-finite point are reporting their Euclidean ccw.
          Otherwise, geodesics intersecting at an ideal point have a ``ccw`` of
          zero, and geodesics intersecting at an ultra-ideal point, are
          producing a ``ValueError``.

        EXAMPLES:

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: v = H.vertical(0)
            sage: g = H.half_circle(0, 2)

            sage: v.ccw(g)
            -1
            sage: v.ccw(-g)
            1

        Two verticals meet intersect with an angle 0 at infinity::

            sage: w = H.vertical(1)
            sage: v.ccw(w)
            0

        However, we can also get the way the chords are turning in the Klein model::

            sage: v.ccw(w, euclidean=True)
            1

        Geodesics that do not intersect produce an error::

            sage: g = H.half_circle(2, 1)
            sage: v.ccw(g)
            Traceback (most recent call last):
            ...
            ValueError: geodesics do not intersect

        However, these geodesics do intersect at an ultra-ideal point and we
        can get their orientation at that point::

            sage: v.ccw(g, euclidean=True)
            1

        """
        if not isinstance(other, HyperbolicOrientedGeodesic):
            raise NotImplementedError("can only compute ccw between oriented geodesics")

        other = self.parent()(other)

        intersection = self._intersection(other)

        if intersection is None:
            if self.unoriented() == other.unoriented():
                return 0
            raise ValueError("geodesics do not intersect")

        if intersection.is_ideal():
            if not euclidean:
                return 0

        if intersection.is_ultra_ideal():
            if not euclidean:
                raise ValueError("geodesics do not intersect")

        sgn = self.parent().geometry._sgn

        from flatsurf.geometry.euclidean import ccw

        return sgn(ccw((self._c, -self._b), (other._c, -other._b)))

    def _test_ccw(self, **options):
        r"""
        Verify that :meth:`ccw` has been implemented correctly

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0)._test_ccw()

        """
        tester = self._tester(**options)

        for euclidean in [False, True]:
            tester.assertEqual(self.ccw(self, euclidean=euclidean), 0)
            tester.assertEqual(self.ccw(-self, euclidean=euclidean), 0)

            for other in self.parent().some_subsets():
                if other.dimension() != 1:
                    continue

                other = other.geodesic().change(oriented=True)

                intersection = self._intersection(other)
                if intersection is None:
                    continue

                if intersection.is_ultra_ideal() and not euclidean:
                    continue

                tester.assertEqual(
                    self.ccw(other, euclidean=euclidean),
                    -other.ccw(self, euclidean=euclidean),
                )
                tester.assertEqual(
                    self.ccw(-other, euclidean=euclidean),
                    other.ccw(self, euclidean=euclidean),
                )

    def angle(self, other, numerical=True):
        r"""
        Compute the angle between this geodesic and ``other`` divided by 2π,
        i.e., as a number in `[0, 1)`.

        If this geodesic and ``other`` do not intersect, a ``ValueError`` is
        raised.

        INPUT:

        - ``other`` -- a geodesic intersecting this geodesic in a finite or
          ideal point.

        - ``numerical`` -- a boolean (default: ``True``); whether a numerical
          approximation of the angle is returned or whether we attempt to
          render it as a rational number.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)

            sage: g0 = H.geodesic(-1, 1)
            sage: g1 = H.geodesic(0, 2)
            sage: g2 = H.geodesic(1, 2)

            sage: g0.angle(g1, numerical=False)
            1/6

            sage: assert g0.angle(g1, numerical=False) == (-g0).angle(-g1, numerical=False)
            sage: assert g1.angle(g2, numerical=False) == (-g2).angle(-g1, numerical=False)
            sage: assert g0.angle(g1, numerical=False) + g1.angle(-g0, numerical=False) == 1/2
            sage: assert g1.angle(g2, numerical=False) + g1.angle(-g2, numerical=False) == 1/2

        ::

            sage: H.geodesic(0, 1).angle(H.geodesic(1, Infinity))
            0.5

            sage: H.geodesic(0, 1).angle(H.geodesic(Infinity, 1))
            0.0

        ::

            sage: m = matrix(2, [2, 1, 1, 1])

            sage: g0.apply_isometry(m).angle(g1.apply_isometry(m), numerical=False)
            1/6

        ::

            sage: a = H.point(0, 1, model='half_plane')
            sage: b = H.point(1, 1, model='half_plane')
            sage: c = H.point(1, 2, model='half_plane')

            sage: H.geodesic(a, b).angle(H.geodesic(b, c))
            0.32379180882521663

            sage: H.geodesic(b, a).angle(H.geodesic(c, b))
            0.32379180882521663

            sage: H.geodesic(a, b).angle(H.geodesic(b, c)) + H.geodesic(b, c).angle(H.geodesic(b, a))
            0.5

        ::

            sage: g3 = H.geodesic(2, 3)
            sage: g0.angle(g3)
            Traceback (most recent call last):
            ...
            ValueError: geodesics do not intersect

        """
        if not isinstance(other, HyperbolicOrientedGeodesic):
            raise NotImplementedError(
                "can only compute angle between oriented geodesics"
            )

        other = self.parent()(other)

        from sage.all import QQ

        ring = float if numerical else QQ

        ccw = self.ccw(other)

        if ccw == 0:
            if self.start() == other.end() or self.end() == other.start():
                return ring(0.5)

            assert (
                self.start() == other.start() or self.end() == other.end()
            ), "geodesics whose enclosed angle is zero, must have one ideal endpoint in common"

            return ring(0)

        a1 = self._a
        b1 = self._b
        c1 = self._c
        n1_squared = b1 * b1 + c1 * c1 - a1 * a1

        a2 = other._a
        b2 = other._b
        c2 = other._c
        n2_squared = b2 * b2 + c2 * c2 - a2 * a2

        import math

        n12 = math.sqrt(abs(n1_squared * n2_squared))

        cos_angle = (b1 * b2 + c1 * c2 - a1 * a2) / n12

        from flatsurf.geometry.euclidean import acos

        angle = acos(cos_angle, numerical=numerical)

        return angle if ccw > 0 else 1 - angle

    def _test_angle(self, **options):
        r"""
        Verify that :meth:`angle` has been implemented correctly.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0)._test_angle()

        """
        tester = self._tester(**options)

        tester.assertEqual(self.angle(self), 0)
        tester.assertEqual(self.angle(-self), 0.5)

        for other in self.parent().some_subsets():
            if other.dimension() != 1:
                continue

            other = other.geodesic().change(oriented=True)

            intersection = self._intersection(other)
            if intersection is None:
                continue

            tester.assertAlmostEqual(self.angle(other), (1 - other.angle(self)) % 1)
            tester.assertAlmostEqual(self.angle(-other), (0.5 - other.angle(self)) % 1)


class HyperbolicPoint(HyperbolicConvexSet, Element):
    r"""
    A (possibly infinite or even ultra-ideal) point in the
    :class:`HyperbolicPlane`.

    This is the abstract base class providing shared functionality for
    :class:`HyperbolicPointFromCoordinates` and
    :class:`HyperbolicPointFromGeodesic`.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

    A point with coordinates in the upper half plane::

        sage: p = H(0)

    The same point, created as an endpoint of a geodesic::

        sage: p = H.vertical(0).start()

    Another point on the same geodesic, a finite point::

        sage: p = H(I)

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicPoint
        sage: isinstance(p, HyperbolicPoint)
        True

        sage: TestSuite(p).run()

    .. SEEALSO::

        :meth:`HyperbolicPlane.point` for ways to create points

    """

    def _check(self, require_normalized=True):
        r"""
        Verify that this is a (possibly infinite) point in the hyperbolic plane.

        Implements :meth:`HyperbolicConvexSet._check`.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); ignored

        EXAMPLES::


            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: p = H(0)
            sage: p._check()

            sage: p = H.point(2, 0, model="klein", check=False)
            sage: p._check()
            Traceback (most recent call last):
            ...
            ValueError: point (2, 0) is not in the unit disk in the Klein model

        .. SEEALSO::

            :meth:`is_ultra_ideal` to see whether a point is outside of the Klein disk

        """
        if self.is_ultra_ideal():
            raise ValueError(
                f"point {self.coordinates(model='klein')} is not in the unit disk in the Klein model"
            )

    def is_ideal(self):
        r"""
        Return whether this point is at infinity.

        This implements :meth:`HyperbolicConvexSet.is_ideal`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: H(0).is_ideal()
            True

            sage: H(I).is_ideal()
            False

        .. NOTE::

            This check is not very robust over inexact rings and should be improved.

        .. SEEALSO::

            :meth:`is_ultra_ideal`, :meth:`is_finite`

        """
        x, y = self.coordinates(model="klein")
        # We should use a specialized predicate from the geometry here to be
        # more robust over inexact rings.
        return self.parent().geometry._cmp(x * x + y * y, 1) == 0

    def is_ultra_ideal(self):
        r"""
        Return whether this point is ultra-ideal, i.e., outside of the Klein
        disk.

        This implements :meth;`HyperbolicConvexSet.is_ultra_ideal`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicGeodesic
            sage: H = HyperbolicPlane()

            sage: H(0).is_ultra_ideal()
            False

            sage: H(I).is_ultra_ideal()
            False

            sage: H.point(2, 0, model="klein", check=False).is_ultra_ideal()
            True

        .. NOTE::

            This check is not very robust over inexact rings and should be improved.

        .. SEEALSO::

            :meth:`is_ideal`, :meth:`is_finite`

        """
        x, y = self.coordinates(model="klein")
        # We should use a specialized predicate from the geometry here to be
        # more robust over inexact rings.
        return self.parent().geometry._cmp(x * x + y * y, 1) > 0

    def is_finite(self):
        r"""
        Return whether this point is finite.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).is_finite()
            False

            sage: H(I).is_finite()
            True

            sage: H.half_circle(0, 2).start().is_finite()
            False

        .. NOTE::

            Currently, the implementation is not robust over inexact rings.

        .. SEEALSO::

            :meth:`is_ideal`, :meth:`is_ultra_ideal`

        """
        x, y = self.coordinates(model="klein")
        # We should use specialized predicate from the geometry implementation
        # here to make this more robust over inexact rings.
        return self.parent().geometry._cmp(x * x + y * y, 1) < 0

    def half_spaces(self):
        r"""
        Return a set of half spaces whose intersection is this point.

        Implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        .. SEEALSO::

            :meth:`HyperbolicPlane.intersection` and
            :meth:`HyperbolicPlane.polygon` to intersect half spaces

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
        r"""
        Return coordinates of this point in ``ring``.

        If ``model`` is ``"half_plane"``, return coordinates in the upper half
        plane model.

        If ``model`` is ``"klein"``, return Euclidean coordinates in the Klein
        model.

        INPUT:

        - ``model`` -- either ``"half_plane"`` or ``"klein"`` (default: ``"half_plane"``)

        - ``ring`` -- a ring, ``"maybe"``, or ``None`` (default: ``None``); in
          which ring the coordinates should be returned. If ``None``,
          coordinates are returned in the :meth:`HyperbolicPlane.base_ring`. If
          ``"maybe"``, same as ``None`` but instead of throwing an exception if
          the coordinates do not exist in the base ring, ``None`` is returned.

        .. NOTE::

            It would be good to add a ``"extension"`` mode here to
            automatically take a ring extension where the coordinates live.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(I).coordinates()
            (0, 1)

            sage: H(0).coordinates()
            (0, 0)

        The point at infinity has no coordinates in the upper half plane model::

            sage: H(oo).coordinates()
            Traceback (most recent call last):
            ...
            ValueError: point has no coordinates in the upper half plane

        Some points have coordinates in the Klein model but not in the upper
        half plane model::

            sage: p = H.point(1/2, 0, model="klein")

            sage: p.coordinates(model="klein")
            (1/2, 0)

            sage: p.coordinates()
            Traceback (most recent call last):
            ...
            ValueError: square root of 3/4 not in Rational Field

            sage: K.<a> = NumberField(x^2 - 3/4)
            sage: p.coordinates(ring=K)
            (1/2, a)

        Some points have no coordinates in either model unless we pass to a
        field extension::

            sage: p = H.half_circle(0, 2).start()

            sage: p.coordinates(model="half_plane")
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: p.coordinates(model="klein")
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: K.<a> = QuadraticField(2)
            sage: p.coordinates(ring=K)
            (-a, 0)

            sage: p.coordinates(ring=K, model="klein")
            (-2/3*a, 1/3)

        """
        coordinates = self._coordinates_klein(ring=ring)

        if coordinates is None:
            assert ring == "maybe"
            return coordinates

        if model == "klein":
            pass
        elif model == "half_plane":
            x, y = coordinates

            if self == self.parent().infinity() or self.is_ultra_ideal():
                raise ValueError("point has no coordinates in the upper half plane")

            denominator = 1 - y

            if not self.is_finite():
                return (x / denominator, self.parent().base_ring().zero())

            square = 1 - x * x - y * y
            try:
                sqrt = square.sqrt()
                if sqrt not in x.parent():
                    raise ValueError(f"square root of {square} not in {x.parent()}")
            except ValueError:
                if ring == "maybe":
                    return None
                raise

            coordinates = x / denominator, sqrt / denominator
        else:
            raise NotImplementedError("cannot determine coordinates in this model yet")

        assert coordinates is not None
        return coordinates

    def real(self):
        r"""
        Return the real part of this point in the upper half plane model.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H(I + 2)
            sage: p.real()
            2

        .. SEEALSO::

            :meth:`imag`

        """
        return self.coordinates(model="half_plane")[0]

    def imag(self):
        r"""
        Return the imaginary part of this point in the upper half plane model.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()

            sage: p = H(I + 2)
            sage: p.imag()
            1

        .. SEEALSO::

            :meth:`real`

        """
        return self.coordinates(model="half_plane")[1]

    def segment(self, end):
        r"""
        Return the oriented segment from this point to ``end``.

        INPUT:

        - ``end`` -- another :class:`HyperbolicPoint`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).segment(I)
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≤ 0}

        A geodesic is returned when both endpoints are ideal points::

            sage: H(0).segment(1) == H.geodesic(0, 1)
            True

        .. SEEALSO::

            :meth:`HyperbolicPlane.segment` for other ways to construct segments

        """
        end = self.parent()(end)

        geodesic = self.parent().geodesic(self, end)

        if self.is_ideal() and end.is_ideal():
            return geodesic

        return self.parent().segment(
            geodesic,
            start=self if self.is_finite() else None,
            end=end if end.is_finite() else None,
            assume_normalized=True,
        )

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this subset.

        See :meth:`HyperbolicConvexPolygon.plot` for the supported keyword
        arguments.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"`` (default: ``"half_plane"``)

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(I).plot()
            ...Graphics object consisting of 1 graphics primitive

        """
        coordinates = self.coordinates(model=model, ring="maybe")
        if not coordinates:
            from sage.all import RR

            coordinates = self.coordinates(model=model, ring=RR)

        from sage.all import point

        # We need to wrap the coordinates into a list so they are not
        # interpreted as a list of complex numbers.
        plot = point([coordinates], **kwds)

        return self._enhance_plot(plot, model=model)

    def dimension(self):
        r"""
        Return the dimension of this point, i.e., 0.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(0).dimension()
            0

        """
        from sage.all import ZZ

        return ZZ.zero()

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this point, i.e., this point.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored since a
          point cannot have marked vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(oo).vertices()
            {∞,}

        Vertices can be computed even if they do not have coordinates over the
        :meth:`HyperbolicPlane.base_ring`::

            sage: H.half_circle(0, 2).start().vertices()
            {-1.41421356237310,}

        .. SEEALSO::

            :meth:`HyperbolicConvexSet.vertices` for more details.

        """
        return HyperbolicVertices([self])

    def __hash__(self):
        r"""
        Return a hash value for this point.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

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
        r"""
        Return equations that must be satisfied if this set converts to
        ``image`` under ``isometry`` using ``λ`` as a free variable for the
        scaling factor.

        Helper method for :meth:`HyperbolicPlane._isometry_from_primitives`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix describing a (right) isometry; typically
          not over the base ring but in symbolic variables

        - ``image`` -- a point

        - ``λ`` -- a symbolic variable

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: R.<a, b, c, d, λ> = QQ[]
            sage: isometry = H._isometry_gl2_to_sim12(matrix([[a, b], [c, d]]))
            sage: isometry
            [                            b*c + a*d                             a*c - b*d                             a*c + b*d]
            [                            a*b - c*d 1/2*a^2 - 1/2*b^2 - 1/2*c^2 + 1/2*d^2 1/2*a^2 + 1/2*b^2 - 1/2*c^2 - 1/2*d^2]
            [                            a*b + c*d 1/2*a^2 - 1/2*b^2 + 1/2*c^2 - 1/2*d^2 1/2*a^2 + 1/2*b^2 + 1/2*c^2 + 1/2*d^2]

            sage: H(I)._isometry_equations(isometry, H(2*I), λ)
            [-8/5*a*c - 2/5*b*d,
             -4/5*a^2 - 1/5*b^2 + 4/5*c^2 + 1/5*d^2,
             -4/5*a^2 - 1/5*b^2 - 4/5*c^2 - 1/5*d^2 + λ]

        """
        R = λ.parent()

        x, y, z = (*self.coordinates(model="klein"), R.one())
        fx, fy, fz = (*image.coordinates(model="klein"), R.one())

        from sage.all import vector

        equations = λ * vector((x, y, z)) - isometry * vector(R, (fx, fy, fz))
        return equations.list()

    def __contains__(self, point):
        r"""
        Return whether the set comprised of this point contains ``point``,
        i.e., whether the points are equal.

        This implements :meth:`HyperbolicConvexSet.__contains__` without
        relying on :meth:`HyperbolicConvexSet.half_spaces` which can not be
        computed for points without coordinates in the
        :meth:`HyperbolicPlane.base_ring`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H(0)
            sage: q = H.half_circle(0, 2).start()

            sage: p in p
            True

            sage: p in q
            False

            sage: q in p
            False

            sage: q in q
            True

        """
        return self == point


class HyperbolicPointFromCoordinates(HyperbolicPoint):
    r"""
    A :class:`HyperbolicPoint` with explicit coordinates in the Klein model.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic import HyperbolicPointFromCoordinates
        sage: H = HyperbolicPlane()

        sage: H.point(0, 0, model="klein")
        I

    Points with coordinates in the half plane model are also stored as points
    in the Klein model::

        sage: p = H.point(0, 0, model="half_plane")

        sage: isinstance(p, HyperbolicPointFromCoordinates)
        True

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` containing this point

    - ``x`` -- an element of :meth:`HyperbolicPlane.base_ring`

    - ``y`` -- an element of :meth:`HyperbolicPlane.base_ring`

    .. SEEALSO::

        Use :meth;`HyperbolicPlane.point` to create points from coordinates

    """

    def __init__(self, parent, x, y):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H.point(0, 0, model="klein")
            sage: TestSuite(p).run()

        """
        super().__init__(parent)

        if x.parent() is not parent.base_ring():
            raise TypeError("x must be an element of the base ring")
        if y.parent() is not parent.base_ring():
            raise TypeError("y must be an element of the base ring")

        self._coordinates = (x, y)

    def _coordinates_klein(self, ring):
        r"""
        Return the coordinates of this point in the Klein model.

        This is a helper method for :meth:`HyperbolicPoint.coordinates`.

        INPUT:

        - ``ring`` -- see :meth:`HyperbolicPoint.coordinates` for this
          parameter

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H.point(0, 0, model="klein")
            sage: p._coordinates_klein(ring=None)
            (0, 0)

        """
        x, y = self._coordinates

        from sage.categories.all import Rings

        if ring is None or ring == "maybe":
            pass
        elif ring in Rings():
            x, y = ring(x), ring(y)
        else:
            raise NotImplementedError("cannot produce coordinates for this ring yet")

        return x, y

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other`` with respect to the ``op``
        operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two points are the same.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
            return all(
                self.parent().geometry._equal(a, b)
                for (a, b) in zip(
                    self.coordinates(model="klein"), other.coordinates(model="klein")
                )
            )

        return super()._richcmp_(other, op)

    def _repr_(self):
        r"""
        Return a printable representation of this point.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We represent points in the upper half plane model if possible::

            sage: H.point(0, 0, model="klein")
            I

        For some points this is not possible without extending the coordinate
        ring. Then we show their coordinates in the Klein model::

            sage: H.point(1/2, 0, model="klein")
            (1/2, 0)

        """
        if self == self.parent().infinity():
            return "∞"

        if self.is_ultra_ideal():
            return repr(self.coordinates(model="klein"))

        coordinates = self.coordinates(model="half_plane", ring="maybe")
        if coordinates is None:
            return repr(self.coordinates(model="klein"))

        from sage.all import PowerSeriesRing

        # We represent x + y*I in R[[I]] so we do not have to reimplement printing ourselves.
        return repr(
            PowerSeriesRing(self.parent().base_ring(), names="I")(list(coordinates))
        )

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this point.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the point
          will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          point.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since points cannot have
          an orientation.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We change the base ring over which the point is defined::

            sage: p = H(0)
            sage: p.change(ring=AA)
            0

        We cannot make a point oriented::

            sage: p.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot make a point oriented

            sage: p.change(oriented=False) == p
            True

        """

        def point(parent):
            return parent.point(*self._coordinates, model="klein", check=False)

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("cannot make a point oriented")

        if ring is not None or geometry is not None:
            self = (
                self.parent()
                .change_ring(ring, geometry=geometry)
                .point(*self._coordinates, model="klein", check=False)
            )

        return self

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the result of applying the ``isometry`` to this hyperbolic
        point.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
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

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random hyperbolic point.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the point

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicPointFromCoordinates
            sage: x = HyperbolicPointFromCoordinates.random_set(H)

            sage: x.dimension()
            0

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        return parent.point(
            parent.base_ring().random_element(),
            parent.base_ring().random_element().abs(),
            model="half_plane",
            check=False,
        )


class HyperbolicPointFromGeodesic(HyperbolicPoint):
    r"""
    An ideal :class:`HyperbolicPoint`, the end point of a geodesic.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.vertical(0).start()
        0

    This class is necessary because not all points have coordinates in the
    :meth:`HyperbolicPlane.base_ring`::

        sage: p = H.half_circle(0, 2).start()

        sage: p.coordinates(model="klein")
        Traceback (most recent call last):
        ...
        ValueError: ...

        sage: p.coordinates(model="half_plane")
        Traceback (most recent call last):
        ...
        ValueError: ...

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` containing this point

    - ``geodesic`` -- to :class:`HyperbolicOrientedGeodesic` whose :meth:`HyperbolicOrientedGeodesic.start` this point is

    .. SEEALSO::

        Use :meth:`HyperbolicOrientedGeodesic.start` and
        :meth:`HyperbolicOrientedGeodesic.end` to create endpoints of geodesics

    """

    def __init__(self, parent, geodesic):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicPointFromGeodesic
            sage: H = HyperbolicPlane()

            sage: p = H.half_circle(0, 2).start()
            sage: isinstance(p, HyperbolicPointFromGeodesic)
            True

            sage: TestSuite(p).run()

        """
        super().__init__(parent)

        if not isinstance(geodesic, HyperbolicOrientedGeodesic):
            raise TypeError("x must be an oriented geodesic")

        self._geodesic = geodesic

    def is_ideal(self):
        r"""
        Return whether this is an infinite point, i.e., return ``True`` since
        this is an endpoint of a geodesic.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).start().is_ideal()
            True

        """
        return True

    def is_ultra_ideal(self):
        r"""
        Return whether this is an ultra ideal point, i.e., return ``False``
        since this end point of a geodesic is not outside of the unit disk in
        the Klein model.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).start().is_ultra_ideal()
            False

        """
        return False

    def is_finite(self):
        r"""
        Return whether this is a finite point, i.e., return ``False`` since
        this end point of a geodesic is infinite.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_circle(0, 2).start().is_finite()
            False

        """
        return False

    @cached_method
    def _coordinates_klein(self, ring):
        r"""
        Return the coordinates of this point in the Klein model.

        This is a helper method for :meth:`HyperbolicPoint.coordinates`.

        INPUT:

        - ``ring`` -- see :meth:`HyperbolicPoint.coordinates` for this
          parameter


        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: p = H.point(0, 0, model="klein")
            sage: p._coordinates_klein(ring=None)
            (0, 0)

        Since the coordinates are given by intersection a line with the unit
        circle, they might only exist over a quadratic extension::

            sage: p = H.half_circle(0, 2).start()
            sage: p._coordinates_klein(ring=None)
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: K.<a> = QuadraticField(2)
            sage: p._coordinates_klein(ring=K)
            (-2/3*a, 1/3)

        .. NOTE::

            The implementation of this predicate is not numerically robust over
            inexact rings and should be improved.

        """
        a, b, c = self._geodesic.equation(model="half_plane")

        # We should probably use a specialized predicate of the geometry to make this
        # more robust over inexact rings.
        if self.parent().geometry._zero(a):
            point = None
            # We should probably use a specialized predicate of the geometry to make this
            # more robust over inexact rings.
            if self.parent().geometry._sgn(b) > 0:
                point = self.parent().point(0, 1, model="klein", check=False)
            else:
                point = self.parent().point(-c / b, 0, model="half_plane", check=False)

            return point.coordinates(model="klein", ring=ring)
        else:
            discriminant = b * b - 4 * a * c

            if ring is None or ring == "maybe":
                sqrt_ring = self.parent().base_ring()
            else:
                sqrt_ring = ring

            discriminant = sqrt_ring(discriminant)

            try:
                sqrt = discriminant.sqrt()
                if sqrt not in sqrt_ring:
                    raise ValueError(
                        f"square root of {discriminant} not in {sqrt_ring}"
                    )
            except ValueError:
                if ring == "maybe":
                    return None
                raise

            endpoints = ((-b - sqrt) / (2 * a), (-b + sqrt) / (2 * a))

            return (
                self.parent()
                .change_ring(sqrt_ring)
                .point(
                    # We should probably use a specialized predicate of the geometry to make this
                    # more robust over inexact rings.
                    (min if self.parent().geometry._sgn(a) > 0 else max)(endpoints),
                    0,
                    model="half_plane",
                    check=False,
                )
                .coordinates(model="klein", ring=ring)
            )

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other`` with respect to the ``op``
        operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two points are the same.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
            return (
                other in self._geodesic
                and self._geodesic.parametrize(other, model="euclidean") < 0
            )

        return super()._richcmp_(other, op)

    def _repr_(self):
        """
        Return a printable representation of this point.

        EXAMPLES:

        We try to represent this point in the upper half plane::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).start()
            0

        When this is not possible, we show approximate coordinates::

            sage: H.half_circle(0, 2).start()
            -1.41421356237310

        TESTS::

            sage: geodesic = H.geodesic(733833/5522174119, -242010/5522174119, -105111/788882017, model="klein")
            sage: geodesic.start()
            3.03625883227966

        """
        if self == self.parent().infinity():
            return "∞"

        if not self.is_ultra_ideal():
            coordinates = self.coordinates(model="half_plane", ring="maybe")

            if coordinates is not None:
                return repr(
                    self.parent().point(*coordinates, model="half_plane", check=False)
                )

        from sage.all import RR

        return repr(self.change_ring(RR))

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this point.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the point
          will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          point.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since points cannot have
          an orientation.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We change the base ring over which the point is defined::

            sage: p = H.half_circle(0, 2).start()
            sage: p.change(ring=AA)
            -1.414213562373095?

        We cannot make a point oriented::

            sage: p.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot change orientation of a point

            sage: p.change(oriented=False) == p
            True

        """
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

            sage: from flatsurf import HyperbolicPlane

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
        if self.coordinates(model="klein", ring="maybe") is not None:
            return super().__hash__()

        from sage.categories.all import Fields

        if self.parent().base_ring() not in Fields():
            raise NotImplementedError(
                "cannot hash point defined by a geodesic over a non-field yet"
            )
        return hash((type(self), self._geodesic))

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the result of applying the ``isometry`` to this hyperbolic
        point.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: point = H.half_circle(0, 2).start()
            sage: point
            -1.41421356237310

        We apply an isometry of positive determinant::

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: point._apply_isometry_klein(isometry)
            -0.414213562373095

        We apply an isometry of negative determinant::

            sage: isometry = matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: point._apply_isometry_klein(isometry)
            1.41421356237310

        """
        image = self._geodesic.apply_isometry(
            isometry, model="klein", on_right=on_right
        )

        if isometry.det().sign() == -1:
            # An isometry of negative determinant swaps the end points of the
            # geodesic
            image = -image

        return image.start()


class HyperbolicConvexPolygon(HyperbolicConvexFacade):
    r"""
    A (possibly unbounded) closed polygon in the :class:`HyperbolicPlane`,
    i.e., the intersection of a finite number of :class:`half spaces
    <HyperbolicHalfSpace>`.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` of which this is a subset

    - ``half_spaces`` -- the :class:`HyperbolicHalfSpace` of which this is an intersection

    - ``vertices`` -- marked vertices that should additionally be kept track of

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: P = H.polygon([
        ....:     H.vertical(0).right_half_space(),
        ....:     H.vertical(1).left_half_space(),
        ....:     H.half_circle(0, 2).left_half_space()])

    .. SEEALSO::

        Use :meth:`HyperbolicPlane.polygon` and
        :meth:`HyperbolicPlane.intersection` to create polygons in the
        hyperbolic plane.

    """

    def __init__(self, parent, half_spaces, vertices, category=None):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexPolygon
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space()])

            sage: isinstance(P, HyperbolicConvexPolygon)
            True

            sage: TestSuite(P).run()

        """
        if category is None:
            from flatsurf.geometry.categories import HyperbolicPolygons

            category = HyperbolicPolygons(parent.base_ring()).Convex().Simple()

        super().__init__(parent, category=category)

        if not isinstance(half_spaces, HyperbolicHalfSpaces):
            raise TypeError("half_spaces must be HyperbolicHalfSpaces")

        self._half_spaces = half_spaces
        self._marked_vertices = tuple(vertices)

    def _check(self, require_normalized=True):
        r"""
        Verify that the marked vertices of this polygon are actually on the edges of the polygon.

        This implements :meth:`HyperbolicConvexSet._check`.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); whether to
          assume that normalization has already happened, i.e., marked vertices
          that are actual vertices have already been removed

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexPolygon
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space()],
            ....:     marked_vertices=[4], check=False)

            sage: P._check()
            Traceback (most recent call last):
            ...
            ValueError: marked vertex must be on an edge of the polygon

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space()],
            ....:     marked_vertices=[oo], assume_minimal=True, check=False)

            sage: P._check()
            Traceback (most recent call last):
            ...
            ValueError: marked vertex must not be a non-marked vertex of the polygon

        """
        for vertex in self._marked_vertices:
            if not any(vertex in edge for edge in self.edges(marked_vertices=False)):
                raise ValueError("marked vertex must be on an edge of the polygon")

        if require_normalized:
            if any(
                vertex in self.vertices(marked_vertices=False)
                for vertex in self._marked_vertices
            ):
                raise ValueError(
                    "marked vertex must not be a non-marked vertex of the polygon"
                )

    def _normalize(self, marked_vertices=False):
        r"""
        Return a convex set describing the intersection of the half spaces
        underlying this polygon.

        This implements :meth:`HyperbolicConvexSet._normalize`.

        The half spaces are assumed to be already sorted respecting
        :meth:`HyperbolicHalfSpaces._lt_`.

        ALGORITHM:

        We compute the intersection of the half spaces in the Klein model in
        several steps:

        * Drop trivially redundant half spaces, e.g., repeated ones.
        * Handle the case that the intersection is empty or a single point, see
          :meth:`_normalize_euclidean_boundary`.
        * Compute the intersection of the corresponding half spaces in the
          Euclidean plane, see :meth:`_normalize_drop_euclidean_redundant`.
        * Remove redundant half spaces that make no contribution for the unit
          disk of the Klein model, see
          :meth:`_normalize_drop_unit_disk_redundant`.
        * Determine of which nature (point, segment, line, polygon) the
          intersection of half spaces is and return the resulting set.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``False``); whether to
          keep marked vertices when normalizing

        .. NOTE::

            Over inexact rings, this is probably mostly useless.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
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
        boundary = self._normalize_euclidean_boundary()

        if not isinstance(boundary, HyperbolicHalfSpace):
            # When there was no such segment, i.e., the intersection is empty
            # or just a point, we are done.
            return boundary

        # Compute a minimal subset of the half spaces that defines the intersection in the Euclidean plane.
        self = self._normalize_drop_euclidean_redundant(boundary)

        # Remove half spaces that make no contribution when restricting to the unit disk of the Klein model.
        self = self._normalize_drop_unit_disk_redundant()

        marked_vertices = [
            vertex for vertex in marked_vertices if vertex not in self.vertices()
        ]

        if marked_vertices:
            if self.dimension() < 2:
                raise NotImplementedError(
                    "cannot add marked vertices to low dimensional objects"
                )

            self = self.parent().polygon(
                self.half_spaces(),
                check=False,
                assume_sorted=True,
                assume_minimal=True,
                marked_vertices=marked_vertices,
            )
            self = self._normalize_drop_marked_vertices()

        return self

    def _normalize_drop_trivially_redundant(self):
        r"""
        Return a sublist of the ``half_spaces`` defining this polygon without
        changing their intersection by removing some trivially redundant half
        spaces.

        The ``half_spaces`` are assumed to be sorted consistent with
        :meth:`HyperbolicHalfSpaces._lt_`.

        This is a helper method for :meth:`_normalize`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

                equal = self.parent().geometry._equal
                sgn = self.parent().geometry._sgn
                if equal(c * B, C * b) and sgn(b) == sgn(B) and sgn(c) == sgn(C):
                    # The half spaces are parallel in the Euclidean plane. Since we
                    # assume spaces to be sorted by inclusion, we can drop this
                    # space.
                    continue

            reduced.append(half_space)

        return self.parent().polygon(
            reduced,
            check=False,
            assume_sorted=True,
            assume_minimal=True,
            marked_vertices=False,
        )

    def _normalize_drop_euclidean_redundant(self, boundary):
        r"""
        Return a minimal sublist of the ``half_spaces`` defining this polygon
        that describe their intersection as half spaces of the Euclidean plane.

        Consider the half spaces in the Klein model. Ignoring the unit disk,
        they also describe half spaces in the Euclidean plane.

        The half space ``boundary`` must be one of the ``half_spaces`` that
        defines a boundary edge of the intersection polygon in the Euclidean
        plane.

        This is a helper method for :meth:`_normalize`.

        ALGORITHM:

        We use an approach similar to gift-wrapping (but from the inside) to remove
        redundant half spaces from the input list. We start from the
        ``boundary`` which is one of the minimal half spaces and extend to the
        full intersection by walking the sorted half spaces.

        Since we visit each half space once, this reduction runs in linear time
        in the number of half spaces.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
            sage: P._normalize_drop_euclidean_redundant(boundary=P._normalize_euclidean_boundary())
            {(x^2 + y^2) - x ≥ 0} ∩ {2*x - 1 ≥ 0} ∩ {x - 1 ≥ 0}

        .. NOTE::

            There are some additional assumptions on the input than what is
            stated here. Please refer to the implementation.

        """
        half_spaces = list(self._half_spaces)

        half_spaces = (
            half_spaces[half_spaces.index(boundary) :]
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
            required_half_spaces,
            check=False,
            assume_sorted=True,
            assume_minimal=True,
            marked_vertices=False,
        )

    def _normalize_drop_unit_disk_redundant(self):
        r"""
        Return the intersection of the Euclidean ``half_spaces`` defining this
        polygon with the unit disk.

        The ``half_spaces`` must be minimal to describe their intersection in
        the Euclidean plane. If that intersection does not intersect the unit
        disk, then return the :meth:`HyperbolicPlane.empty_set`.

        Otherwise, return a minimal sublist of ``half_spaces`` that describes
        the intersection inside the unit disk.

        This is a helper method for :meth:`_normalize`.

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

            sage: from flatsurf import HyperbolicPlane
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

        .. NOTE::

            There are some additional assumptions on the input than what is
            stated here. Please refer to the implementation.

        """
        required_half_spaces = []

        maybe_empty = True
        maybe_point = True
        maybe_segment = True

        for A, B, C in self._half_spaces.triples(repeat=True):
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
                    # Unsurprisingly, pylint gets confused by maybe_point being
                    # both a boolean and a point at times. The code should
                    # probably be cleaned up. But here, it must be a point so
                    # the call is save.
                    # pylint: disable=no-member
                    assert not maybe_point.is_finite()
                    # pylint: enable=no-member
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
            required_half_spaces,
            check=False,
            assume_sorted=True,
            assume_minimal=True,
            marked_vertices=False,
        )

    def _normalize_euclidean_boundary(self):
        r"""
        Return a half space whose (Euclidean) boundary intersects the boundary
        of the intersection of the ``half_spaces`` defining this polygon in
        more than a point.

        Consider the half spaces in the Klein model. Ignoring the unit disk,
        they also describe half spaces in the Euclidean plane.

        If their intersection contains a segment it must be on the boundary of
        one of the ``half_spaces`` which is returned by this method.

        If this is not the case, and the intersection is empty in the
        hyperbolic plane, return the :meth:`HyperbolicPlane.empty_set`.
        Otherwise, if the intersection is a point in the hyperbolic plane,
        return that point.

        The ``half_spaces`` must already be sorted with respect to
        :meth:`HyperbolicHalfSpaces._lt_`.

        This is a helper method for :meth:`_normalize`.

        ALGORITHM:

        We initially ignore the hyperbolic structure and just consider the half
        spaces of the Klein model as Euclidean half spaces.

        We use a relatively standard randomized optimization approach to find a
        point in the intersection: we randomly shuffle the half spaces and then
        optimize a segment on some boundary of the half spaces. The
        randomization makes this a linear time algorithm, see e.g.,
        https://www2.cs.arizona.edu/classes/cs437/spring21/Lecture4.pdf.

        If the only segment we can construct is a point, then the intersection
        is a single point in the Euclidean plane. The intersection in the
        hyperbolic plane might be a single point or empty.

        If not even a point exists, the intersection is empty in the Euclidean
        plane and therefore empty in the hyperbolic plane.

        Note that the segment returned might not be within the unit disk.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        Make the following randomized tests reproducible::

            sage: set_random_seed(0)

        An intersection that is already empty in the Euclidean plane::

            sage: polygon(
            ....:     H.geodesic(2, 1/2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).left_half_space()
            ....: )._normalize_euclidean_boundary()
            {}

        An intersection which in the Euclidean plane is a single point but
        outside the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(0, -1, 0, model="klein"),
            ....:     H.half_space(2, 2, -1, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....: )._normalize_euclidean_boundary()
            {}

        An intersection which is a single point inside the unit disk::

            sage: polygon(*H(I).half_spaces())._normalize_euclidean_boundary()
            I

        An intersection which is a single point on the boundary of the unit
        disk::

            sage: polygon(*H.infinity().half_spaces())._normalize_euclidean_boundary()
            {x - 1 ≥ 0}

        An intersection which is a segment outside of the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_euclidean_boundary()
            {x ≤ 0}

        An intersection which is a polygon outside of the unit disk::

            sage: polygon(
            ....:     H.half_space(0, 1, 0, model="klein"),
            ....:     H.half_space(1, -2, 0, model="klein"),
            ....:     H.half_space(-2, -2, 1, model="klein"),
            ....:     H.half_space(17/8, 2, -1, model="klein"),
            ....: )._normalize_euclidean_boundary()
            {9*(x^2 + y^2) + 32*x + 25 ≥ 0}

        An intersection which is an (unbounded) polygon touching the unit disk::

            sage: polygon(
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.vertical(1).right_half_space(),
            ....: )._normalize_euclidean_boundary()
            {x - 1 ≥ 0}

        An intersection which is a segment touching the unit disk::

            sage: polygon(
            ....:     H.vertical(0).left_half_space(),
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(-1).left_half_space(),
            ....:     H.geodesic(-1, -2).right_half_space(),
            ....: )._normalize_euclidean_boundary()
            {x ≥ 0}

        An intersection which is a polygon inside the unit disk::

            sage: polygon(
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.geodesic(0, -1).right_half_space(),
            ....:     H.geodesic(0, 1).left_half_space(),
            ....: )._normalize_euclidean_boundary()
            {(x^2 + y^2) - x ≥ 0}

        A polygon which has no vertices inside the unit disk but intersects the unit disk::

            sage: polygon(
            ....:     H.geodesic(2, 3).left_half_space(),
            ....:     H.geodesic(-3, -2).left_half_space(),
            ....:     H.geodesic(-1/2, -1/3).left_half_space(),
            ....:     H.geodesic(1/3, 1/2).left_half_space(),
            ....: )._normalize_euclidean_boundary()
            {6*(x^2 + y^2) - 5*x + 1 ≥ 0}

        A single half plane::

            sage: polygon(
            ....:     H.vertical(0).left_half_space()
            ....: )._normalize_euclidean_boundary()
            {x ≤ 0}

        A pair of anti-parallel half planes::

            sage: polygon(
            ....:     H.geodesic(1/2, 2).left_half_space(),
            ....:     H.geodesic(-1/2, -2).right_half_space(),
            ....: )._normalize_euclidean_boundary()
            {2*(x^2 + y^2) - 5*x + 2 ≥ 0}

        TESTS:

        A case that caused problems at some point::

            sage: set_random_seed(1)

            sage: polygon(
            ....:    H.geodesic(300, 3389, -1166, model="half_plane").right_half_space(),
            ....:    H.geodesic(5, -24, -5, model="half_plane").left_half_space(),
            ....:    H.geodesic(182, -1135, 522, model="half_plane").left_half_space(),
            ....: )._normalize_euclidean_boundary()
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
                interval = RealSet(-oo, oo)

                for constraining in random_half_spaces:
                    if constraining is half_space:
                        break

                    intersection = boundary._intersection(constraining.boundary())

                    if intersection is None:
                        # constraining is anti-parallel to half_space
                        if (
                            boundary.unparametrize(0, model="euclidean", check=False)
                            not in constraining
                        ):
                            return self.parent().empty_set()

                        # The intersection is non-empty, so this adds no further constraints.
                        continue

                    λ = boundary.parametrize(
                        intersection, model="euclidean", check=False
                    )

                    # RealSet in SageMath does not like number fields. We move
                    # everything through AA (which might not always work) to
                    # work around this problem.
                    if λ.parent().is_exact():
                        from sage.all import AA

                        rλ = AA(λ)
                    else:
                        rλ = λ

                    # Determine whether this half space constrains to (-∞, λ] or [λ, ∞).
                    if (
                        boundary.unparametrize(λ + 1, model="euclidean", check=False)
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

                point = boundary.unparametrize(λ, model="euclidean", check=False)

        return self._normalize_extend_to_euclidean_boundary(point)

    def _normalize_extend_to_euclidean_boundary(self, point):
        r"""
        Extend ``point`` to a (Euclidean) half space which intersects the
        intersection of the ``half_spaces`` defining this polygon in more than
        one point.

        This is a helper method for :meth:`_normalize_euclidean_boundary`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        A helper to create non-normalized polygons for testing::

            sage: polygon = lambda *half_spaces: H.polygon(half_spaces, check=False, assume_sorted=False, assume_minimal=True)

        We extend from a single point to half space::

            sage: P = polygon(*H.infinity().half_spaces())

            sage: P._normalize_extend_to_euclidean_boundary(H.infinity())
            {x ≤ 0}

        """
        half_spaces = [
            half_space
            for half_space in self._half_spaces
            if point in half_space.boundary()
        ]

        if len(half_spaces) == 0:
            raise ValueError("point must be on the boundary of a defining half space")

        if len(half_spaces) < 3:
            return half_spaces[0]

        for i, half_space in enumerate(half_spaces):
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
        r"""
        Return a copy of this polygon with marked vertices removed that are
        already vertices of the polygon anyway.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space()],
            ....:     marked_vertices=[oo], assume_minimal=True, check=False)
            sage: P
            {x - 1 ≤ 0} ∩ {x ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0} ∪ {∞}

            sage: P._normalize_drop_marked_vertices()
            {x - 1 ≤ 0} ∩ {x ≥ 0} ∩ {(x^2 + y^2) - 2 ≥ 0}

        """
        vertices = [
            vertex
            for vertex in self._marked_vertices
            if vertex not in self.vertices(marked_vertices=False)
        ]

        return self.parent().polygon(
            self._half_spaces,
            check=False,
            assume_sorted=True,
            assume_minimal=True,
            marked_vertices=vertices,
        )

    def dimension(self):
        r"""
        Return the dimension of this polygon, i.e., 2.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        Note that this also returns 2 if the actual dimension of the polygon is
        smaller. This is, however, only possible for polygons created with
        :meth:`HyperbolicPlane.polygon` setting ``check=False``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 1).left_half_space(),
            ....:     H.half_circle(0, 4).right_half_space(),
            ....: ])
            sage: P.dimension()
            2

        """
        from sage.all import ZZ

        return ZZ(2)

    @cached_method
    def edges(self, as_segments=False, marked_vertices=True):
        r"""
        Return the :class:`segments <HyperbolicOrientedSegment>` and
        :class:`geodesics <HyperbolicOrientedGeodesic>` defining this polygon.

        This implements :meth:`HyperbolicConvexSet.edges` for polygons.

        INPUT:

        - ``as_segments`` -- a boolean (default: ``False``); whether to also
          return the geodesics as segments with ideal end points.

        - ``marked_vertices`` -- a boolean (default: ``True``); if set, edges
          with end points at a marked vertex are reported, otherwise, marked
          vertices are completely ignored.

        OUTPUT:

        A set of segments and geodesics. Iteration through this set is in
        counterclockwise order with respect to the points of the set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        The edges of a polygon::

            sage: P = H.intersection(
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.vertical(1).left_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space())

            sage: P.edges()
            {{-x + 1 = 0} ∩ {2*(x^2 + y^2) - 5*x - 3 ≤ 0}, {-(x^2 + y^2) + 4 = 0} ∩ {(x^2 + y^2) - 5*x + 1 ≥ 0} ∩ {(x^2 + y^2) + 5*x + 1 ≥ 0}, {x + 1 = 0} ∩ {2*(x^2 + y^2) + 5*x - 3 ≤ 0}, {(x^2 + y^2) - 1 = 0}}

            sage: [type(e) for e in P.edges()]
            [<class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedGeodesic_with_category_with_category'>]

            sage: [type(e) for e in P.edges(as_segments=True)]
            [<class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>,
             <class 'flatsurf.geometry.hyperbolic.HyperbolicOrientedSegment_with_category_with_category'>]

        The edges of a polygon with marked vertices::

            sage: P = H.convex_hull(-1, 1, I, 2*I, marked_vertices=True)
            sage: P.edges()
            {{-(x^2 + y^2) - 3*x + 4 = 0} ∩ {3*(x^2 + y^2) - 25*x - 12 ≤ 0}, {-(x^2 + y^2) + 3*x + 4 = 0} ∩ {3*(x^2 + y^2) + 25*x - 12 ≤ 0}, {(x^2 + y^2) - 1 = 0} ∩ {x ≤ 0}, {(x^2 + y^2) - 1 = 0} ∩ {x ≥ 0}}
            sage: P.edges(marked_vertices=False)
            {{-(x^2 + y^2) - 3*x + 4 = 0} ∩ {3*(x^2 + y^2) - 25*x - 12 ≤ 0}, {-(x^2 + y^2) + 3*x + 4 = 0} ∩ {3*(x^2 + y^2) + 25*x - 12 ≤ 0}, {(x^2 + y^2) - 1 = 0}}

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

        if marked_vertices and self._marked_vertices:
            edges_without_marked_vertices = edges
            edges = []
            for edge in edges_without_marked_vertices:
                vertices_on_edge = [
                    vertex for vertex in self._marked_vertices if vertex in edge
                ]

                def key(vertex, geodesic=edge.geodesic()):
                    return geodesic.parametrize(vertex, model="euclidean")

                vertices_on_edge.sort(key=key)
                vertices_on_edge.append(edge.end())

                start = edge.start()
                for vertex in vertices_on_edge:
                    edges.append(
                        self.parent().segment(
                            edge.geodesic(),
                            start=start,
                            end=vertex,
                            assume_normalized=as_segments,
                            check=False,
                        )
                    )
                    start = vertex

        return HyperbolicEdges(edges)

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this polygon, i.e., the (possibly ideal) end
        points of the :meth:`edges`.

        INPUT:

        -- ``marked_vertices`` -- a boolean (default: ``True``) whether to
        include marked vertices in the output

        OUTPUT:

        Returns a set of points. Iteration over this set is in counterclockwise
        order.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        edges = self.edges(marked_vertices=False)

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
        r"""
        Return a minimal set of half spaces whose intersection this polygon is.

        This implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        Marked vertices are not encoded in the half spaces in any way::

            sage: P = H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ], marked_vertices=[I + 1])
            sage: P
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∪ {1 + I}

            sage: H.polygon(P.half_spaces())
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        """
        return self._half_spaces

    def _repr_(self):
        r"""
        Return a printable representation of this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:   H.vertical(1).left_half_space(),
            ....:   H.vertical(-1).right_half_space(),
            ....:   H.half_circle(0, 1).left_half_space(),
            ....:   H.half_circle(0, 4).right_half_space(),
            ....: ], marked_vertices=[I + 1])
            sage: P
            {x - 1 ≤ 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∪ {1 + I}

        """
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

        - ``color`` -- a string (default: ``"#efffff"``); the fill color of
          polygons

        - ``edgecolor`` -- a string (default: ``"blue"``); the color of
          geodesics and segments

        See :func:`flatsurf.graphical.hyperbolic.hyperbolic_path` for additional keyword arguments to
        customize the plot.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

        Similarly in the Klein model picture, the arc of infinite points is
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
        kwds.setdefault("color", "#efffff")
        kwds.setdefault("edgecolor", "blue")

        if len(self._half_spaces) == 0:
            raise NotImplementedError("cannot plot full space")

        edges = self.edges(as_segments=True, marked_vertices=False)

        pos = edges[0].start()

        from flatsurf.graphical.hyperbolic import HyperbolicPathPlotCommand

        commands = [HyperbolicPathPlotCommand("MOVETO", pos)]

        for edge in edges:
            if edge.start() != pos:
                commands.append(HyperbolicPathPlotCommand("MOVETO", edge.start()))

            commands.append(HyperbolicPathPlotCommand("LINETO", edge.end()))
            pos = edge.end()

        if pos != edges[0].start():
            commands.append(HyperbolicPathPlotCommand("MOVETO", edges[0].start()))

        from flatsurf.graphical.hyperbolic import hyperbolic_path

        plot = hyperbolic_path(commands, model=model, **kwds)

        return self._enhance_plot(plot, model=model)

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this polygon.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the polygon
          will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          polygon.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since polygons cannot
          have an explicit orientation. See :meth:`~HyperbolicConvexSet.is_oriented`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We change the ring over which a polygon is defined::

            sage: P = H.polygon([
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space(),
            ....:     H.half_circle(0, 1).left_half_space()],
            ....:     marked_vertices=[I])

            sage: P.change(ring=AA)
            {x - 1 ≤ 0} ∩ {x + 1 ≥ 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∪ {I}

        We cannot give a polygon an explicit orientation::

            sage: P.change(oriented=False) == P
            True

            sage: P.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: polygons cannot have an explicit orientation

        """
        if ring is not None or geometry is not None:
            self = (
                self.parent()
                .change_ring(ring, geometry=geometry)
                .polygon(
                    [
                        half_space.change(ring=ring, geometry=geometry)
                        for half_space in self._half_spaces
                    ],
                    check=False,
                    assume_sorted=True,
                    assume_minimal=True,
                    marked_vertices=[
                        vertex.change(ring=ring, geometry=geometry)
                        for vertex in self._marked_vertices
                    ],
                )
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("polygons cannot have an explicit orientation")

        return self

    def __eq__(self, other):
        r"""
        Return whether this polygon is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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
        if not isinstance(other, HyperbolicConvexPolygon):
            return False

        return (
            self._half_spaces == other._half_spaces
            and self._marked_vertices == other._marked_vertices
        )

    def __hash__(self):
        r"""
        Return a hash value for this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since polygons are hashable, they can be put in a hash table, such
        as a Python ``set``::

            sage: S = {H.polygon([H.vertical(1).left_half_space(), H.vertical(-1).right_half_space()])}

        """
        return hash((self._half_spaces, self._marked_vertices))

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the image of this polygon under ``isometry``.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_circle(0, 4).left_half_space()],
            ....:     marked_vertices=[4*I])
            sage: P._apply_isometry_klein(isometry)
            {(x^2 + y^2) - 2*x - 3 ≥ 0} ∩ {x - 1 ≥ 0} ∪ {1 + 4*I}

        """
        half_spaces = [
            h.apply_isometry(isometry, model="klein", on_right=on_right)
            for h in self._half_spaces
        ]
        marked_vertices = [
            p.apply_isometry(isometry, model="klein", on_right=on_right)
            for p in self._marked_vertices
        ]

        return self.parent().polygon(
            half_spaces=half_spaces,
            check=False,
            assume_minimal=True,
            marked_vertices=marked_vertices,
        )

    def _isometry_conditions(self, other):
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: P = H.polygon([
            ....:     H.vertical(0).right_half_space(),
            ....:     H.half_circle(0, 4).left_half_space()],
            ....:     marked_vertices=[4*I])

            sage: conditions = P._isometry_conditions(P)
            sage: list(conditions)
            [[({x ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {x ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {x ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {x ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              (4*I, 4*I)],
             [({x ≥ 0}, {(x^2 + y^2) - 4 ≥ 0}),
              ({(x^2 + y^2) - 4 ≥ 0}, {x ≥ 0}),
              (4*I, 4*I)]]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet._isometry_conditions` for a general description.

        r"""
        # We are likely returning too many conditions here.
        # In particular there are many more ways to determine that no isometry
        # can possibly exist here.

        for reverse in [False, True]:
            preimage = list(self.half_spaces())
            image = list(other.half_spaces())

            if len(preimage) != len(image):
                continue

            if reverse:
                image.reverse()

            for i in range(len(image)):
                image = image[1:] + image[:1]

                if self._marked_vertices:
                    preimage_vertices = list(self._marked_vertices)
                    image_vertices = list(other._marked_vertices)

                    if len(preimage_vertices) != len(image_vertices):
                        continue

                    for j in range(len(image_vertices)):
                        image_vertices = image_vertices[1:] + image_vertices[:1]

                        yield list(
                            zip(preimage + preimage_vertices, image + image_vertices)
                        )
                else:
                    yield list(zip(preimage, image))

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random hyperbolic convex polygon.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the polygon

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexPolygon
            sage: x = HyperbolicConvexPolygon.random_set(H)

            sage: x.dimension()
            2

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        while True:
            from sage.all import ZZ

            interior_points = []
            count = ZZ.random_element().abs() + 3

            while len(interior_points) < count:
                p = HyperbolicPointFromCoordinates.random_set(parent)
                if p.is_ideal():
                    continue
                interior_points.append(p)

            half_spaces = []

            while len(half_spaces) < len(interior_points):
                half_space = HyperbolicHalfSpace.random_set(parent)

                for p in interior_points:
                    if p in half_space:
                        continue

                    a, b, c = half_space.equation(model="klein")

                    x, y = p.coordinates(model="klein")

                    a = -(b * x + c * y)

                    half_space = parent.half_space(a, b, c, model="klein")

                    assert p in half_space

                half_spaces.append(half_space)

            polygon = parent.polygon(half_spaces)

            if isinstance(polygon, HyperbolicConvexPolygon):
                return polygon

    def is_degenerate(self):
        r"""
        Return whether this is considered to be a degenerate polygon.

        EXAMPLES:

        We consider polygons of area zero as degenerate::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: P = H.polygon([
            ....:     H.vertical(0).left_half_space(),
            ....:     H.half_circle(0, 1).left_half_space(),
            ....:     H.half_circle(0, 2).right_half_space(),
            ....:     H.vertical(0).right_half_space()
            ....: ], check=False, assume_minimal=True)
            sage: P.is_degenerate()
            True

        We also consider polygons with marked points as degenerate::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: P = H.polygon([
            ....:     H.vertical(1).left_half_space(),
            ....:     H.half_circle(0, 2).left_half_space(),
            ....:     H.half_circle(0, 4).right_half_space(),
            ....:     H.vertical(-1).right_half_space()
            ....: ], marked_vertices=[2*I])
            sage: P.is_degenerate()
            True

            sage: H.polygon(P.half_spaces()).is_degenerate()
            False

        Finally, we consider polygons with ideal points as degenerate::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: P = H.polygon([
            ....:     H.vertical(1).left_half_space(),
            ....:     H.vertical(-1).right_half_space()
            ....: ])
            sage: P.is_degenerate()
            True

        .. NOTE::

            This is not a terribly meaningful notion. This exists mostly
            because degenerate polygons have a more obvious meaning in
            Euclidean geometry where this check is used when rendering a
            polygon as a string.

        """
        if self.parent().polygon(self.half_spaces()) != self:
            return True

        return not self.is_finite()


class HyperbolicSegment(HyperbolicConvexFacade):
    r"""
    A segment (possibly infinite) in the hyperbolic plane.

    This is an abstract base class of :class:`HyperbolicOrientedSegment` and
    :class:`HyperbolicUnorientedSegment`.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` containing this segment

    - ``geodesic`` -- the :class:`HyperbolicGeodesic` of which this segment is a subset

    - ``start`` -- a :class:`HyperbolicPoint` or ``None`` (default: ``None``);
      the finite endpoint of the segment. If ``None``, then the segment extends
      all the way to the ideal starting point of the geodesic.

    - ``end`` -- a :class:`HyperbolicPoint` or ``None`` (default: ``None``);
      the finite endpoint of the segment. If ``None``, then the segment extends
      all the way to the ideal end point of the geodesic.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.segment(H.vertical(0), start=I)
        {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        sage: H(I).segment(oo)
        {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

    .. SEEALSO::

        Use :meth:`HyperbolicPlane.segment` or :meth:`HyperbolicPoint.segment`
        to create segments.

    """

    def __init__(self, parent, geodesic, start=None, end=None):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: segment = H.segment(H.vertical(0), start=I)
            sage: isinstance(segment, HyperbolicSegment)
            True

            sage: TestSuite(segment).run()

            sage: segment = segment.unoriented()
            sage: isinstance(segment, HyperbolicSegment)
            True

            sage: TestSuite(segment).run()

        """
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
        r"""
        Verify that this is a valid segment.

        This implements :meth:`HyperbolicConvexSet._check`.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); whether to
          assume that normalization has already happened, i.e., segments with
          no finite end points have been rewritten into geodesics.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicConvexPolygon
            sage: H = HyperbolicPlane()

        The end points of the segment must be on the defining geodesic::

            sage: s = H.segment(H.vertical(0), start=1, check=False, assume_normalized=True)
            sage: s._check()
            Traceback (most recent call last):
            ...
            ValueError: start point must be on the geodesic

            sage: s = H.segment(H.vertical(0), end=1, check=False, assume_normalized=True)
            sage: s._check()
            Traceback (most recent call last):
            ...
            ValueError: end point must be on the geodesic

        The end points must be ordered correctly::

            sage: s = H.segment(H.vertical(0), start=2*I, end=I, check=False, assume_normalized=True)
            sage: s._check()
            Traceback (most recent call last):
            ...
            ValueError: end point of segment must not be before start point on the underlying geodesic

        The end points must be distinct::

            sage: s = H.segment(H.vertical(0), start=I, end=I, check=False, assume_normalized=True)
            sage: s._check(require_normalized=False)

            sage: s = H.segment(H.vertical(0), start=I, end=I, check=False, assume_normalized=True)
            sage: s._check(require_normalized=True)
            Traceback (most recent call last):
            ...
            ValueError: end point of segment must be after start point on the underlying geodesic

        """
        if self._start is not None:
            if self._start not in self._geodesic:
                raise ValueError("start point must be on the geodesic")

        if self._end is not None:
            if self._end not in self._geodesic:
                raise ValueError("end point must be on the geodesic")

        start = (
            None
            if self._start is None
            else self._geodesic.parametrize(self._start, model="euclidean", check=False)
        )
        end = (
            None
            if self._end is None
            else self._geodesic.parametrize(self._end, model="euclidean", check=False)
        )

        if start is not None and end is not None:
            if end < start:
                raise ValueError(
                    "end point of segment must not be before start point on the underlying geodesic"
                )
            if require_normalized and end <= start:
                raise ValueError(
                    "end point of segment must be after start point on the underlying geodesic"
                )

    def _normalize(self):
        r"""
        Return this set possibly rewritten in a simpler form.

        This implements :meth:`HyperbolicConvexSet._normalize`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We define a helper method for easier testing::

            sage: segment = lambda *args, **kwds: H.segment(*args, **kwds, check=False, assume_normalized=True)

        A segment that consists of an ideal point, is just that point::

            sage: segment(H.vertical(-1), start=H.infinity(), end=H.infinity())._normalize()
            ∞

            sage: segment(H.vertical(0), start=H.infinity(), end=None)._normalize()
            ∞

            sage: segment(-H.vertical(0), start=None, end=H.infinity())._normalize()
            ∞

        A segment that has two ideal end points is a geodesic::

            sage: segment(H.vertical(0), start=None, end=H.infinity())._normalize()
            {-x = 0}

            sage: segment(-H.vertical(0), start=H.infinity(), end=None)._normalize()
            {x = 0}

        Segments that remain segments in normalization::

            sage: segment(H.vertical(0), start=I, end=H.infinity())._normalize()
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

            sage: segment(-H.vertical(0), start=H.infinity(), end=I)._normalize()
            {x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        .. NOTE::

            This method is not numerically robust and should be improve over inexact rings.

        """
        if self._geodesic.is_ultra_ideal():
            return self.parent().empty_set()

        start = self._start
        end = self._end

        def λ(point):
            return self._geodesic.parametrize(point, model="euclidean", check=False)

        if start is not None and end is not None:
            if λ(start) > λ(end):
                raise ValueError(
                    "end point of segment must be after start point on the underlying geodesic"
                )

        if start is not None:
            if not start.is_finite():
                # We should use a specialized predicate of geometry to make
                # this more robust over inexact rings.
                if self.parent().geometry._sgn(λ(start)) > 0:
                    return (
                        self.parent().empty_set() if start.is_ultra_ideal() else start
                    )
                start = None

        if end is not None:
            if not end.is_finite():
                # We should use a specialized predicate of geometry to make
                # this more robust over inexact rings.
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
            self._geodesic,
            start=start,
            end=end,
            check=False,
            assume_normalized=True,
            oriented=self.is_oriented(),
        )

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the image of this segment under ``isometry``.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We apply an isometry to an oriented segment::

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: segment = H(I).segment(2*I)
            sage: segment._apply_isometry_klein(isometry)
            {-x + 1 = 0} ∩ {2*(x^2 + y^2) - 3*x - 1 ≥ 0} ∩ {(x^2 + y^2) - 3*x - 2 ≤ 0}

        We apply an isometry of negative determinant to an oriented segment::

            sage: isometry = matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: segment._apply_isometry_klein(isometry) == segment
            True

            sage: segment.start().apply_isometry(isometry, model="klein") == segment.start()
            True

            sage: segment.end().apply_isometry(isometry, model="klein") == segment.end()
            True

        Note that this behavior is different from how the start and end point
        of a geodesic behave under such an isometry::

            sage: segment.geodesic().apply_isometry(isometry, model="klein") == segment.geodesic()
            False

            sage: segment.geodesic().apply_isometry(isometry, model="klein").start() == segment.geodesic().end()
            True

            sage: segment.geodesic().apply_isometry(isometry, model="klein").end() == segment.geodesic().start()
            True

        We can also apply an isometry to an unoriented geodesic::

            sage: segment.unoriented()._apply_isometry_klein(isometry) == segment.unoriented()
            True

        """
        geodesic = self.geodesic()._apply_isometry_klein(isometry, on_right=on_right)

        if isometry.det().sign() == -1 and geodesic.is_oriented():
            geodesic = -geodesic

        start = (
            self._start.apply_isometry(isometry, model="klein", on_right=on_right)
            if self._start is not None
            else None
        )
        end = (
            self._end.apply_isometry(isometry, model="klein", on_right=on_right)
            if self._end is not None
            else None
        )

        return self.parent().segment(
            geodesic, start=start, end=end, oriented=self.is_oriented()
        )

    def _endpoint_half_spaces(self):
        r"""
        Return the half spaces that stop the segment at its endpoints.

        A helper method for printing and :meth:`half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: s = H.segment(H.half_circle(0, 1), end=I)
            sage: list(s._endpoint_half_spaces())
            [{x ≤ 0}]

        """
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
        r"""
        Return a printable representation of this segment.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: H.segment(H.half_circle(0, 1), end=I)
            {(x^2 + y^2) - 1 = 0} ∩ {x ≤ 0}

        """
        bounds = [repr(self._geodesic)]
        bounds.extend(repr(half_space) for half_space in self._endpoint_half_spaces())

        return " ∩ ".join(bounds)

    def half_spaces(self):
        r"""
        Return a minimal set of half spaces whose intersection is this segment.

        This implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: segment = H.segment(H.half_circle(0, 1), end=I)
            sage: segment.half_spaces()
            {{x ≤ 0}, {(x^2 + y^2) - 1 ≤ 0}, {(x^2 + y^2) - 1 ≥ 0}}

        """
        return self._geodesic.half_spaces() + HyperbolicHalfSpaces(
            self._endpoint_half_spaces()
        )

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of this segment.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` or ``"klein"`` (default:
          ``"half_plane"``); in which model to produce the plot

        See :func:`flatsurf.graphical.hyperbolic.hyperbolic_path` for additional keyword arguments to
        customize the plot.

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

            sage: segment = H.segment(H.half_circle(0, 1), end=I)
            sage: segment.plot()  # long time (.25s)
            ...Graphics object consisting of 1 graphics primitive

        """
        self = self.change(oriented=True)

        from sage.all import RR

        kwds["fill"] = False

        self = self.change_ring(RR)

        from flatsurf.graphical.hyperbolic import (
            hyperbolic_path,
            HyperbolicPathPlotCommand,
        )

        plot = hyperbolic_path(
            [
                HyperbolicPathPlotCommand("MOVETO", self.start()),
                HyperbolicPathPlotCommand("LINETO", self.end()),
            ],
            model=model,
            **kwds,
        )

        return self._enhance_plot(plot, model=model)

    def __eq__(self, other):
        r"""
        Return whether this segment is indistinguishable from ``other`` (except
        for scaling in the defining geodesic's equation).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
            sage: H = HyperbolicPlane()

        Oriented segments are equal if they have the same start and end points::

            sage: H(I).segment(2*I) == H(2*I).segment(I)
            False

        For an unoriented segment the endpoints must be the same but order does not matter::

            sage: H(I).segment(2*I).unoriented() == H(2*I).segment(I).unoriented()
            True

        """
        if type(self) is not type(other):
            return False
        return (
            self.geodesic() == other.geodesic() and self.vertices() == other.vertices()
        )

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this segment.

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~HyperbolicPlane.base_ring`); the ring over which the new half
          space will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new half space.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since half spaces cannot
          have an explicit orientation. See :meth:`~HyperbolicConvexSet.is_oriented`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        We change the ring over which the segment is defined::

            sage: s = H(I).segment(2*I)

            sage: s.change(ring=AA)
            {-x = 0} ∩ {3/5*(x^2 + y^2) - 3/5 ≥ 0} ∩ {6/25*(x^2 + y^2) - 24/25 ≤ 0}

        We make the segment unoriented::

            sage: s.change(oriented=False).is_oriented()
            False

        We pick a (somewhat) random orientation of an unoriented segment::

            sage: s.unoriented().change(oriented=True).is_oriented()
            True

        """
        if ring is not None or geometry is not None:
            start = (
                self._start.change(ring=ring, geometry=geometry)
                if self._start is not None
                else None
            )
            end = (
                self._end.change(ring=ring, geometry=geometry)
                if self._end is not None
                else None
            )

            self = (
                self.parent()
                .change_ring(ring=ring, geometry=geometry)
                .segment(
                    self._geodesic.change(ring=ring, geometry=geometry),
                    start=start,
                    end=end,
                    check=False,
                    assume_normalized=True,
                    oriented=self.is_oriented(),
                )
            )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
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
        r"""
        Return the geodesic on which this segment lies.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s.geodesic()
            {-x = 0}

        Since the segment is oriented, the geodesic is also oriented::

            sage: s.is_oriented()
            True

            sage: s.geodesic().is_oriented()
            True

            sage: s.unoriented().geodesic().is_oriented()
            False

        .. SEEALSO::

            geodesics also implement this method so that segments and geodesics
            can be treated uniformly, see :meth:`HyperbolicGeodesic.geodesic`

        """
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s.vertices()
            {I, 2*I}

        Note that iteration in the set is not consistent with the orientation
        of the segment (it is chosen such that the subset relation on vertices
        can be checked quickly)::

            sage: (-s).vertices()
            {I, 2*I}

        Use :meth:`~HyperbolicOrientedSegment.start` and
        :meth:`~HyperbolicOrientedSegment.end` to get the vertices in an order
        that is consistent with the orientation::

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
        self = self.change(oriented=True)
        return HyperbolicVertices([self.start(), self.end()])

    def dimension(self):
        r"""
        Return the dimension of this segment, i.e., 1.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        Note that this also returns 1 if the actual dimension of the segment is
        smaller. This is, however, only possible for segments created with
        :meth:`HyperbolicPlane.segment` setting ``check=False``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H(I).segment(2*I).dimension()
            1

        """
        from sage.all import ZZ

        return ZZ(1)

    def midpoint(self):
        r"""
        Return the midpoint of this segment.

        ALGORITHM:

        We use the construction as explained on `Wikipedia
        <https://en.wikipedia.org/wiki/Beltrami%E2%80%93Klein_model#Compass_and_straightedge_constructions>`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
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

            :meth:`HyperbolicSegment.perpendicular` for the perpendicular bisector

        """
        start, end = self.vertices()

        if start == end:
            return start

        if not start.is_finite() and not end.is_finite():
            return self.geodesic().midpoint()

        if not start.is_finite() or not end.is_finite():
            raise NotImplementedError(
                f"cannot compute midpoint of unbounded segment {self}"
            )

        for p in self.geodesic().perpendicular(start).vertices():
            for q in self.geodesic().perpendicular(end).vertices():
                line = self.parent().geodesic(p, q)
                intersection = self.intersection(line)
                if intersection:
                    return intersection

            # One of the two lines start at any p must intersect the segment
            # already. No need to check the other p.
            assert (
                False
            ), f"segment {self} must have a midpoint but the straightedge and compass construction did not yield any"

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

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicSegment
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
            raise ValueError(
                f"point must be in the segment but {point} is not in {self}"
            )

        return self.geodesic().perpendicular(point)


class HyperbolicUnorientedSegment(HyperbolicSegment):
    r"""
    An unoriented (possibly infinity) segment in the hyperbolic plane.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: segment = H.segment(H.vertical(0), start=I).unoriented()

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicUnorientedSegment
        sage: isinstance(segment, HyperbolicUnorientedSegment)
        True

    .. SEEALSO::

        Use :meth:`HyperbolicPlane.segment` or
        :meth:`~HyperbolicConvexSet.unoriented` to construct unoriented
        segments.

    """

    def __hash__(self):
        r"""
        Return a hash value for this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(2*I)

        Since an oriented segment is hashable, it can be put in a hash table,
        such as a Python ``set``::

            sage: {s.unoriented(), (-s).unoriented()}
            {{-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}}

        """
        return hash((frozenset([self._start, self._end]), self.geodesic()))

    def _isometry_conditions(self, other):
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I).unoriented()

            sage: conditions = s._isometry_conditions(s)
            sage: list(conditions)
            [[({-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0},
               {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0})],
             [({-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0},
               {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0})]]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet._isometry_conditions` for a general description.

        r"""
        self = self.change(oriented=True)
        other = other.change(oriented=True)

        yield [(self, other)]
        yield [(self, -other)]

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random unoriented segment.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the segment

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicUnorientedSegment
            sage: x = HyperbolicUnorientedSegment.random_set(H)

            sage: x.dimension()
            1

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        return HyperbolicOrientedSegment.random_set(parent).unoriented()


class HyperbolicOrientedSegment(HyperbolicSegment, HyperbolicOrientedConvexSet):
    r"""
    An oriented (possibly infinite) segment in the hyperbolic plane such as a
    boundary edge of a :class:`HyperbolicConvexPolygon`.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: s = H(I).segment(2*I)

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicOrientedSegment
        sage: isinstance(s, HyperbolicOrientedSegment)
        True

    .. SEEALSO::

        Use :meth:`HyperbolicPlane.segment` or :meth:`HyperbolicPoint.segment`
        to construct oriented segments.

    """

    def __neg__(self):
        r"""
        Return this segment with its orientation reversed.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s
            {-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}

            sage: -s
            {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}

        """
        return self.parent().segment(
            -self._geodesic, self._end, self._start, check=False, assume_normalized=True
        )

    def __hash__(self):
        r"""
        Return a hash value for this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()
            sage: s = H(I).segment(2*I)

        Since this set is hashable, it can be put in a hash table, such as a
        Python ``set``::

            sage: {s}
            {{-x = 0} ∩ {(x^2 + y^2) - 1 ≥ 0} ∩ {(x^2 + y^2) - 4 ≤ 0}}

        """
        return hash((self._start, self._end, self.geodesic()))

    def _isometry_conditions(self, other):
        r"""
        Return an iterable of primitive pairs that must map to each other in an
        isometry that maps this set to ``other``.

        Helper method for :meth:`HyperbolicPlane._isometry_conditions`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)

            sage: conditions = s._isometry_conditions(s)
            sage: list(conditions)
            [[(I, I), (2*I, 2*I)]]

        .. SEEALSO::

            :meth:`HyperbolicConvexSet._isometry_conditions` for a general description.

        r"""
        yield [(self.start(), other.start()), (self.end(), other.end())]

    def start(self):
        r"""
        Return the start point of this segment.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s.start()
            I

            sage: s.start().is_finite()
            True

        The start point can also be an ideal point::

            sage: s = H(0).segment(2*I)
            sage: s.start()
            0

            sage: s.start().is_ideal()
            True

        .. SEEALSO::

            :meth:`end`

        """
        if self._start is not None:
            return self._start

        return self._geodesic.start()

    def end(self):
        r"""
        Return the end point of this segment.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: s = H(I).segment(2*I)
            sage: s.end()
            2*I

            sage: s.end().is_finite()
            True

        The end point can also be an ideal point::

            sage: s = H(I).segment(oo)
            sage: s.end()
            ∞

            sage: s.end().is_ideal()
            True

        .. SEEALSO::

            :meth:`start`

        """
        if self._end is not None:
            return self._end

        return self._geodesic.end()

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random oriented segment.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` containing the segment

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicOrientedSegment
            sage: x = HyperbolicOrientedSegment.random_set(H)

            sage: x.dimension()
            1

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        a = HyperbolicPointFromCoordinates.random_set(parent)
        b = HyperbolicPointFromCoordinates.random_set(parent)
        while a == b or (not a.is_finite() and not b.is_finite()):
            b = HyperbolicPointFromCoordinates.random_set(parent)

        return parent.segment(parent.geodesic(a, b), start=a, end=b)


class HyperbolicEmptySet(HyperbolicConvexFacade):
    r"""
    The empty subset of the hyperbolic plane.

    INPUT:

    - ``parent`` -- the :class:`HyperbolicPlane` this is the empty set of

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: H.empty_set()
        {}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicEmptySet

        sage: ø = H.empty_set()
        sage: isinstance(ø, HyperbolicEmptySet)
        True

        sage: TestSuite(ø).run()

    .. SEEALSO::

        Use :meth:`HyperbolicPlane.empty_set` to construct the empty set.

    """

    def __eq__(self, other):
        r"""
        Return whether this empty set is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.empty_set() == H.empty_set()
            True
            sage: H.empty_set() == HyperbolicPlane(AA).empty_set()
            False

        """
        return isinstance(other, HyperbolicEmptySet) and self.parent() == other.parent()

    def some_elements(self):
        r"""
        Return some representative points of this set for testing.

        EXAMPLES:

        Since this set is empty, there are no points to return::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.empty_set().some_elements()
            []

        """
        return []

    def _test_an_element(self, **options):
        r"""
        Do not run tests on an element of this empty set (disabling the generic
        tests run by all parents otherwise).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H._test_an_element()

        """

    def _test_elements(self, **options):
        r"""
        Do not run any tests on the elements of this empty set (disabling the
        generic tests run by all parents otherwise).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H._test_elements()

        """

    def _repr_(self):
        r"""
        Return a printable representation of the empty set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set()
            {}

        """
        return "{}"

    def _apply_isometry_klein(self, isometry, on_right=False):
        r"""
        Return the result of applying the ``isometry`` to the empty set, i.e.,
        the empty set.

        Helper method for :meth:`HyperbolicConvexSet.apply_isometry`.

        INPUT:

        - ``isometry`` -- a 3×3 matrix over the
          :meth:`~HyperbolicPlane.base_ring` describing an isometry in the
          hyperboloid model.

        - ``on_right`` -- a boolean (default: ``False``) whether to return the
          result of the right action.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: isometry = matrix([[1, -1, 1], [1, 1/2, 1/2], [1, -1/2, 3/2]])
            sage: H.empty_set()._apply_isometry_klein(isometry) == H.empty_set()
            True

        """
        return self

    def plot(self, model="half_plane", **kwds):
        r"""
        Return a plot of the empty set.

        INPUT:

        - ``model`` -- one of ``"half_plane"`` and ``"klein"`` (default: ``"half_plane"``)

        Any keyword arguments are ignored.

        EXAMPLES:

        .. jupyter-execute::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().plot()
            ...Graphics object consisting of 0 graphics primitives

        """
        from sage.all import Graphics

        return self._enhance_plot(Graphics(), model=model)

    def dimension(self):
        r"""
        Return the dimension of this set; returns -1 for the empty set.

        This implements :meth:`HyperbolicConvexSet.dimension`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().dimension()
            -1

        """
        from sage.all import ZZ

        return ZZ(-1)

    def half_spaces(self):
        r"""
        Return a minimal set of half spaces whose intersection is empty.

        This implements :meth:`HyperbolicConvexSet.half_spaces`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().half_spaces()
            {{(x^2 + y^2) + 4*x + 3 ≤ 0}, {(x^2 + y^2) - 4*x + 3 ≤ 0}}

        """
        return HyperbolicHalfSpaces(
            [
                self.parent().half_circle(-2, 1).right_half_space(),
                self.parent().half_circle(2, 1).right_half_space(),
            ]
        )

    def change(self, ring=None, geometry=None, oriented=None):
        r"""
        Return a copy of the empty set.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`HyperbolicPlane.base_ring`); the ring over which the empty set
          will be defined.

        - ``geometry`` -- a :class:`HyperbolicGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          empty set.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since the empty set
          cannot have an explicit orientation.

        EXAMPLES:

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().change(ring=AA) == HyperbolicPlane(AA).empty_set()
            True

            sage: H.empty_set().change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot change orientation of empty set

        """
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

            sage: from flatsurf import HyperbolicPlane

            sage: H = HyperbolicPlane()

        Since this set is hashable, it can be put in a hash table, such as a
        Python ``set``::

            sage: {H.empty_set()}
            {{}}

        """
        return 0

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this empty, i.e., an empty set of points.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``), ignored

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().vertices()
            {}

        """
        return HyperbolicVertices([])

    def _an_element_(self):
        """
        Return a point in this set, i.e., raise an exception since there are no
        points.

        See :meth:`HyperbolicConvexSet._an_element_` for more interesting
        examples of this method.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.empty_set().an_element()
            Traceback (most recent call last):
            ...
            Exception: empty set has no points

        """
        raise Exception("empty set has no points")

    @classmethod
    def random_set(cls, parent):
        r"""
        Return a random empty set, i.e., the empty set.

        This implements :meth:`HyperbolicConvexSet.random_set`.

        INPUT:

        - ``parent`` -- the :class:`HyperbolicPlane` this is the empty set of.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: from flatsurf.geometry.hyperbolic import HyperbolicEmptySet
            sage: x = HyperbolicEmptySet.random_set(H)

            sage: x.dimension()
            -1

        .. SEEALSO::

            :meth:`HyperbolicPlane.random_element`

        """
        return parent.empty_set()


class OrderedSet(collections.abc.Set):
    r"""
    A set of objects sorted by :meth:`OrderedSet._lt_`.

    This is used to efficiently represent
    :meth:`HyperbolicConvexSet.half_spaces`,
    :meth:`HyperbolicConvexSet.vertices`, and
    :meth:`HyperbolicConvexSet.edges`. In particular, it allows us to create
    and merge such sets in linear time.

    This is an abstract base class for specialized sets such as
    :class:`HyperbolicHalfSpaces`, :class:`HyperbolicVertices`, and
    :class:`HyperbolicEdges`.

    INPUT:

    - ``entries`` -- an iterable, the elements of this set

    - ``assume_sorted`` -- a boolean or ``"rotated"`` (default: ``True``); whether to assume
      that the ``entries`` are already sorted with respect to :meth:`_lt_`. If
      ``"rotated"``, we assume that the entries are sorted modulo a cyclic permutation.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: segment = H(I).segment(2*I)
        sage: segment.vertices()
        {I, 2*I}

    """

    def __init__(self, entries, assume_sorted=None):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import OrderedSet
            sage: H = HyperbolicPlane()

            sage: vertices = H(I).segment(2*I).vertices()

            sage: isinstance(vertices, OrderedSet)
            True

        """
        if assume_sorted is None:
            assume_sorted = isinstance(entries, OrderedSet)

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
        r"""
        Return whether ``lhs`` should come before ``rhs`` in this set.

        Subclasses must implement this method.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import OrderedSet
            sage: H = HyperbolicPlane()

            sage: vertices = H(I).segment(2*I).vertices()

            sage: vertices._lt_(vertices[0], vertices[1])
            True

        """
        raise NotImplementedError

    def _merge(self, *sets):
        r"""
        Return the merge of sorted lists of ``sets`` using merge sort.

        Note that this set itself is not part of the merge.

        Naturally, when there are a lot of small sets, such a merge sort takes
        quasi-linear time. However, when there are only a few sets, this runs
        in linear time.

        INPUT:

        - ``sets`` -- iterables that are sorted with respect to :meth:`_lt_`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: HyperbolicHalfSpaces([])._merge()
            []

            sage: HyperbolicHalfSpaces([])._merge([H.vertical(0).left_half_space()], [H.vertical(0).left_half_space()])
            [{x ≤ 0}]

            sage: HyperbolicHalfSpaces([])._merge(*[[half_space] for half_space in H.real(0).half_spaces()])
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpaces([])._merge(list(H.real(0).half_spaces()), list(H.real(0).half_spaces()))
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

            sage: HyperbolicHalfSpaces([])._merge(*[[half_space] for half_space in list(H.real(0).half_spaces()) * 2])
            [{(x^2 + y^2) + x ≤ 0}, {x ≥ 0}]

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
                elif self._lt_(B[-1], A[-1]):
                    merged.append(A.pop())
                else:
                    # Drop duplicate from set
                    A.pop()

            merged.reverse()

            return A + B + merged

        # Divide & Conquer recursively.
        return self._merge(
            *(self._merge(*sets[: count // 2]), self._merge(*sets[count // 2 :]))
        )

    def __eq__(self, other):
        r"""
        Return whether this set is equal to ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).vertices() == (-H.vertical(0)).vertices()
            True

        """
        if type(other) is not type(self):
            other = type(self)(other)

        return self._entries == other._entries

    def __ne__(self, other):
        r"""
        Return whether this set is not equal to ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
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

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: hash(H.vertical(0).vertices()) != hash(H.vertical(1).vertices())
            True

        """
        return hash(tuple(self._entries))

    def __add__(self, other):
        r"""
        Return the :meth:`_merge` of this set and ``other``.

        INPUT:

        - ``other`` -- another set of the same kind

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.vertical(0).half_spaces() + H.vertical(1).half_spaces()
            {{x ≤ 0}, {x - 1 ≤ 0}, {x ≥ 0}, {x - 1 ≥ 0}}

        """
        if type(self) is not type(other):
            raise TypeError("both sets must be of the same type")
        entries = self._merge(list(self._entries), list(other._entries))
        return type(self)(entries, assume_sorted=True)

    def __repr__(self):
        r"""
        Return a printable representation of this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.half_circle(0, 1).vertices()
            {-1, 1}

        """
        return "{" + repr(self._entries)[1:-1] + "}"

    def __iter__(self):
        r"""
        Return an iterator of this set.

        Iteration happens in sorted order, consistent with :meth:`_lt_`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: list(vertices)
            [-1, 1]

        """
        return iter(self._entries)

    def __len__(self):
        r"""
        Return the cardinality of this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: len(vertices)
            2

        """
        return len(self._entries)

    def pairs(self, repeat=False):
        r"""
        Return an iterable that iterates over all consecutive pairs of elements
        in this set; including the pair formed by the last element and the
        first element.

        INPUT:

        - ``repeat`` -- a boolean (default: ``False``); whether to produce pair
          consisting of the first element twice if there is only one element

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: list(vertices.pairs())
            [(1, -1), (-1, 1)]

        .. SEEALSO::

            :meth:`triples`

        """
        if len(self._entries) <= 1 and not repeat:
            return

        for i in range(len(self._entries)):
            yield self._entries[i - 1], self._entries[i]

    def triples(self, repeat=False):
        r"""
        Return an iterable that iterates over all consecutive triples of
        elements in this set; including the triples formed by wrapping around
        the end of the set.

        INPUT:

        - ``repeat`` -- a boolean (default: ``False``); whether to produce
          triples by wrapping around even if there are fewer than three
          elements

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

        There must be at least three elements to form any triples::

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: list(vertices.triples())
            []

            sage: half_spaces = H(I).segment(2*I).half_spaces()
            sage: list(half_spaces.triples())
            [({(x^2 + y^2) - 1 ≥ 0}, {x ≤ 0}, {(x^2 + y^2) - 4 ≤ 0}),
             ({x ≤ 0}, {(x^2 + y^2) - 4 ≤ 0}, {x ≥ 0}),
             ({(x^2 + y^2) - 4 ≤ 0}, {x ≥ 0}, {(x^2 + y^2) - 1 ≥ 0}),
             ({x ≥ 0}, {(x^2 + y^2) - 1 ≥ 0}, {x ≤ 0})]

        However, we can force triples to be produced by wrapping around with
        ``repeat``::

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: list(vertices.triples(repeat=True))
            [(1, -1, 1), (-1, 1, -1)]

        .. SEEALSO::

            :meth:`pairs`

        """
        if len(self._entries) <= 2 and not repeat:
            return

        for i in range(len(self._entries)):
            yield self._entries[i - 1], self._entries[i], self._entries[
                (i + 1) % len(self._entries)
            ]

    def __getitem__(self, *args, **kwargs):
        r"""
        Return items from this set by index.

        INPUT:

        Any arguments that can be used to access a tuple are accepted.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: vertices = H.half_circle(0, 1).vertices()
            sage: vertices[0]
            -1

        """
        return self._entries.__getitem__(*args, **kwargs)

    def __contains__(self, x):
        r"""
        Return whether this set contains ``x``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: vertices = H.half_circle(0, 1).vertices()

            sage: H(0) in vertices
            False

            sage: H(1) in vertices
            True

        .. NOTE::

            Presently, this method is not used. It only exists to satisfy the
            conditions of the Python abc for sets. It could be implemented more
            efficiently.

        """
        return x in self._entries


class HyperbolicVertices(OrderedSet):
    r"""
    A set of vertices on the boundary of a convex set in the hyperbolic plane,
    sorted in counterclockwise order.

    INPUT:

    - ``vertices`` -- an iterable of :class:`HyperbolicPoint`, the vertices of this set

    - ``assume_sorted`` -- a boolean or ``"rotated"`` (default: ``True``);
      whether to assume that the ``vertices`` are already sorted with respect
      to :meth:`_lt_`. If ``"rotated"``, we assume that the vertices are sorted
      modulo a cyclic permutation.

    ALGORITHM:

    We keep vertices sorted in counterclockwise order relative to a fixed
    reference vertex (the leftmost and bottommost in the Klein model).

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()
        sage: V = H.vertical(0).vertices()
        sage: V
        {0, ∞}

    Note that in this example, ``0`` is chosen as the reference point::

        sage: V._start
        0

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicVertices
        sage: isinstance(V, HyperbolicVertices)
        True

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.vertices` to obtain such a set

    """

    def __init__(self, vertices, assume_sorted=None):
        r"""
        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicVertices
            sage: H = HyperbolicPlane()
            sage: V = H.vertical(0).vertices()

            sage: isinstance(V, HyperbolicVertices)
            True

        """
        vertices = list(vertices)

        if len(vertices) != 0:
            if assume_sorted or len(vertices) == 1:
                self._start = vertices[0]
            else:
                # We sort vertices in counterclockwise order. We need to fix a starting
                # point consistently, namely we choose the leftmost point in the Klein
                # model and if there are ties the one with minimal y coordinate.
                # That way we can then order vertices by their slopes with the starting
                # point to get a counterclockwise walk.

                # A very common special case is that we are presented with two
                # end points of a geodesic. This is a hack but getting the
                # consistent ordering to work with a mix of points without
                # using their coordinates is a lot of work.
                if (
                    len(vertices) == 2
                    and isinstance(vertices[0], HyperbolicPointFromGeodesic)
                    and isinstance(vertices[1], HyperbolicPointFromGeodesic)
                    and vertices[0]._geodesic == -vertices[1]._geodesic
                ):
                    geodesic = vertices[0]._geodesic
                    # Note that we should use more robust predicates from the "geometry" to
                    # make this work more reliably over inexact rings.
                    if geodesic.parent().geometry._zero(geodesic._c):
                        # This is a vertical in the Klein model. Both end points have the
                        # same x coordinate.
                        if geodesic._b > 0:
                            # This vertical is oriented downwards, so the end point has the
                            # minimal y coordinate.
                            vertices.reverse()
                    elif geodesic._c < 0:
                        # This geodesic points right-to-left in the Klein model. The
                        # end point has minimal x coordinate.
                        vertices.reverse()

                    self._start = vertices[0]
                    assume_sorted = True
                else:
                    self._start = min(
                        vertices, key=lambda vertex: vertex.coordinates(model="klein")
                    )

        # _lt_ needs to know the global structure of the convex hull of the vertices.
        # The base class constructor will replace _entries with a sorted version of _entries.
        self._entries = tuple(vertices)

        super().__init__(vertices, assume_sorted=assume_sorted)

    def _merge(self, *sets):
        r"""
        Return the merge of sorted lists of ``sets``.

        Note that this set itself is not part of the merge (but its reference
        point is used).

        INPUT:

        - ``sets`` -- iterables that are sorted with respect to :meth:`_lt_`.

        .. WARNING::

            For this to work correctly, the result of the merge must eventually
            have the reference point of this set as its reference point.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: V = H.vertical(0).vertices()

            sage: V._merge([H(1)], [H(0)], [H(oo)])
            [0, 1, ∞]

        """
        return super()._merge(*sets)

    def _slope(self, vertex):
        r"""
        Return the slope of ``vertex`` with respect to the chosen reference
        vertex of this set as a tuple (Δy, Δx).

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: V = H.vertical(0).vertices()

        We compute the Euclidean slope from 0 to 1 in the Klein model::

            sage: V._slope(H(1))
            (1, 1)

        """
        sx, sy = self._start.coordinates(model="klein")
        x, y = vertex.coordinates(model="klein")
        return (y - sy, x - sx)

    def _lt_(self, lhs, rhs):
        r"""
        Return whether ``lhs`` should come before ``rhs`` in this set.

        INPUT:

        - ``lhs`` -- a :class:`HyperbolicPoint`

        - ``rhs`` -- a :class:`HyperbolicPoint`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: V = H.vertical(0).vertices()

        We find that we go counterclockwise from 1 to ∞ when seen from 0 in the
        Klein model::

            sage: V._lt_(H(oo), H(1))
            False

        """
        if lhs == self._start:
            return True
        if rhs == self._start:
            return False
        if lhs == rhs:
            return False

        dy_lhs, dx_lhs = self._slope(lhs)
        dy_rhs, dx_rhs = self._slope(rhs)

        assert (
            dx_lhs >= 0 and dx_rhs >= 0
        ), "all points must be to the right of the starting point due to chosen normalization"

        if dy_lhs * dx_rhs < dy_rhs * dx_lhs:
            return True

        if dy_lhs * dx_rhs > dy_rhs * dx_lhs:
            return False

        # The points (start, lhs, rhs) are collinear.
        # In general we cannot decide their order with only looking at start,
        # lhs, and rhs. We need to understand where the rest of the convex hull
        # lives.
        assert (
            lhs in self._entries and rhs in self._entries
        ), "cannot compare vertices that are not defining for the convex hull"

        # If there is any vertex with a bigger slope, then this line is at the
        # start of the walk in counterclockwise order.
        for vertex in self._entries:
            dy, dx = self._slope(vertex)
            if dy * dx_rhs > dy_rhs * dx:
                return dx_lhs < dx_rhs
            elif dy * dx_rhs < dy_rhs * dx:
                return dx_lhs > dx_rhs

        raise ValueError(
            "cannot decide counterclockwise ordering of exactly three collinear points"
        )


class HyperbolicHalfSpaces(OrderedSet):
    r"""
    A set of half spaces in the hyperbolic plane ordered counterclockwise.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: half_spaces = H.vertical(0).half_spaces()
        sage: half_spaces
        {{x ≤ 0}, {x ≥ 0}}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces

        sage: isinstance(half_spaces, HyperbolicHalfSpaces)
        True

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.half_spaces` to obtain such a set
    """

    @classmethod
    def _lt_(cls, lhs, rhs):
        r"""
        Return whether the half space ``lhs`` is smaller than ``rhs`` in a cyclic
        ordering of normal vectors, i.e., order half spaces by whether their
        normal points to the left/right, the slope of the geodesic, and finally
        by containment.

        This ordering is such that :meth:`HyperbolicPlane.intersection` can be
        computed in linear time for two hyperbolic convex sets.

        INPUT:

        - ``lhs`` -- a :class:`HyperbolicHalfSpace`

        - ``rhs`` -- a :class:`HyperbolicHalfSpace`

        .. NOTE::

            The implementation is not very robust over inexact rings and should be improved for that use case.

        TESTS::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
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

    @staticmethod
    def convex_hull(vertices):
        r"""
        Return the convex hull of ``vertices`` as a ordered set of half spaces.

        INPUT:

        - ``vertices`` -- a sequence of :class:`HyperbolicPoint`

        ALGORITHM:

        We use the classical Euclidean Graham scan algorithm in the Klein
        model.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic import HyperbolicHalfSpaces
            sage: H = HyperbolicPlane()

            sage: HyperbolicHalfSpaces.convex_hull([H(0), H(1), H(oo)])
            {{(x^2 + y^2) - x ≥ 0}, {x - 1 ≤ 0}, {x ≥ 0}}

            sage: HyperbolicHalfSpaces.convex_hull([H(0), H(1), H(I), H(oo)])
            {{(x^2 + y^2) - x ≥ 0}, {x - 1 ≤ 0}, {x ≥ 0}}

            sage: HyperbolicHalfSpaces.convex_hull([H(0), H(1), H(I), H(I + 1), H(oo)])
            {{(x^2 + y^2) - x ≥ 0}, {x - 1 ≤ 0}, {x ≥ 0}}

            sage: HyperbolicHalfSpaces.convex_hull([H(1/2), H(-1/2), H(1), H(I), H(I + 1), H(oo)])
            {{2*(x^2 + y^2) - 3*x + 1 ≥ 0}, {x - 1 ≤ 0}, {2*x + 1 ≥ 0}, {4*(x^2 + y^2) - 1 ≥ 0}}

        .. SEEALSO::

            :meth:`HyperbolicPlane.convex_hull`

        """
        if not vertices:
            # We cannot return empty_set() because we do not know the HyperbolicPlane this lives in.
            raise NotImplementedError(
                "cannot compute convex hull of empty set of vertices"
            )

        H = vertices[0].parent()

        vertices = set(vertices)

        if len(vertices) == 1:
            return next(iter(vertices)).half_spaces()

        vertices = [vertex.coordinates(model="klein") for vertex in vertices]
        reference = min(vertices)

        class Slope:
            def __init__(self, xy):
                self.dx = xy[0] - reference[0]
                self.dy = xy[1] - reference[1]

            def __eq__(self, other):
                if self.dx == 0 and self.dy == 0:
                    return other.dx == 0 and other.dy == 0

                # Return whether the two points have the same slope relative to "reference"
                return self.dy * other.dx == other.dy * self.dx

            def __lt__(self, other):
                # Return whether the self has smaller slope relative to "reference" or
                if self.dy * other.dx < other.dy * self.dx:
                    return True
                if self.dy * other.dx > other.dy * self.dx:
                    return False
                # if slopes are the same, sort by distance
                if self.dx**2 + self.dy**2 < other.dx**2 + other.dy**2:
                    return True
                return False

        vertices.sort(key=Slope)

        assert vertices[0] == reference

        # Drop collinear points
        filtered = []
        for i, vertex in enumerate(vertices):
            if i + 1 == len(vertices):
                filtered.append(vertex)
                continue

            slope = Slope(vertex)
            next_slope = Slope(vertices[i + 1])

            if slope == next_slope:
                continue

            filtered.append(vertex)
            continue

        hull = []

        def ccw(A, B, C):
            r"""
            Return whether the vectors A->B, B->C describe a counter-clockwise turn.
            """
            return (B[0] - A[0]) * (C[1] - B[1]) > (C[0] - B[0]) * (B[1] - A[1])

        for vertex in filtered:
            while len(hull) >= 2 and not ccw(hull[-2], hull[-1], vertex):
                hull.pop()
            hull.append(vertex)

        assert hull[0] == reference

        hull = [H.point(*xy, model="klein") for xy in hull]

        half_spaces = []
        for i in range(len(hull)):
            half_spaces.append(H.geodesic(hull[i - 1], hull[i]).left_half_space())

        return HyperbolicHalfSpaces(half_spaces)


class HyperbolicEdges(OrderedSet):
    r"""
    A set of hyperbolic segments and geodesics ordered counterclockwise.

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane()

        sage: edges = H.vertical(0).edges()
        sage: edges
        {{-x = 0}, {x = 0}}

    TESTS::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicEdges
        sage: isinstance(edges, HyperbolicEdges)
        True

    .. SEEALSO::

        :meth:`HyperbolicConvexSet.edges` to obtain such a set

    """

    @classmethod
    def _lt_(cls, lhs, rhs):
        r"""
        Return whether ``lhs`` should come before ``rhs`` in the ordering of this set.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: edges = H.vertical(0).edges()
            sage: edges._lt_(edges[0], edges[1])
            True

        Segments on the same edge are ordered correctly::

            sage: segments = [
            ....:   H(0).segment(H(I)),
            ....:   H(I).segment(H(2*I)),
            ....:   H(2*I).segment(H(oo))
            ....: ]

            sage: edges._lt_(segments[0], segments[1])
            True
            sage: edges._lt_(segments[0], segments[2])
            True
            sage: edges._lt_(segments[1], segments[2])
            True

            sage: edges._lt_(segments[2], segments[1])
            False
            sage: edges._lt_(segments[2], segments[0])
            False
            sage: edges._lt_(segments[1], segments[0])
            False

        """
        lhs_geodesic = lhs
        if isinstance(lhs, HyperbolicOrientedSegment):
            lhs_geodesic = lhs.geodesic()

        rhs_geodesic = rhs
        if isinstance(rhs, HyperbolicOrientedSegment):
            rhs_geodesic = rhs.geodesic()

        if HyperbolicHalfSpaces._lt_(
            lhs_geodesic.left_half_space(), rhs_geodesic.left_half_space()
        ):
            return True

        if lhs_geodesic != rhs_geodesic:
            return False

        if lhs == rhs:
            return False

        # The geodesics containing the edges are the same but they are not the
        # same segments. We compare the finite points on the segment to decide
        # which edge comes first in counterclockwise order.
        geodesic = lhs_geodesic

        if lhs.start().is_ideal():
            if rhs.start().is_ideal():
                assert (
                    not lhs.end().is_ideal() and not rhs.end().is_ideal()
                ), "edges in a set of HyperbolicEdges must be sortable"
                assert (
                    lhs.end() != rhs.end()
                ), "edges were found to be different as segments but they are actually the same"

                return geodesic.parametrize(
                    lhs.end(), model="euclidean"
                ) < geodesic.parametrize(rhs.end(), model="euclidean")

            return True

        if rhs.start().is_ideal():
            return False

        return geodesic.parametrize(
            lhs.start(), model="euclidean"
        ) < geodesic.parametrize(rhs.start(), model="euclidean")
