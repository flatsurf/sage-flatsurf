# TODO: Benchmark how constructions here compare to constructions before we introduced the EuclideanPlane.


r"""
Two dimensional Euclidean geometry.

EXAMPLES::

    sage: from flatsurf import EuclideanPlane
    sage: E = EuclideanPlane(QQ)
    sage: E.circle((0, 0), radius=1)
    { x² + y² = 1 }

.. NOTE::

    Most functionality in this module is also implemented in
    SageMath in the context of linear programming/polyhedra/convex
    geometry. However, that implementation is much more general
    (higher dimensions) and therefore quite inefficient in two
    dimensions.

.. NOTE::

    Explicitly creating Python objects for everything in the Euclidean plane
    comes with a certain overhead. In practice, this overhead is negligible in
    comparison to the cost of performing computations with these objects.
    However, most classes here expose their core algorithms as static methods
    so they can be called with the bare coordinates instead of on explicit
    objects to gain a tiny bit of a speedup where this makes a difference.

.. SEEALSO

    :mod:`flatsurf.geometry.hyperbolic` for the geometry in the hyperbolic plane.

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2013-2020 Vincent Delecroix
#                      2013-2019 W. Patrick Hooper
#                      2020-2025 Julian Rüth
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
from typing import List

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.geometry import Geometry, ExactGeometry, EpsilonGeometry


class EuclideanPlane(Parent, UniqueRepresentation):
    r"""
    The Euclidean plane.

    All objects in the plane must be specified over the given base ring.

    The implemented objects of the plane are mostly convex (points, circles,
    segments, rays, convex polygons.) But some are also non-convex such as
    non-convex polygons.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane(QQ)

    TESTS::

        sage: isinstance(E, EuclideanPlane)
        True

        sage: TestSuite(E).run()

    """

    @staticmethod
    def __classcall__(cls, base_ring=None, geometry=None, category=None):
        r"""
        Create the Euclidean plane with normalized arguments to make it a
        unique SageMath parent.

        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: from flatsurf.geometry.euclidean import EuclideanExactGeometry

            sage: EuclideanPlane() is EuclideanPlane(QQ)
            True

            sage: EuclideanPlane() is EuclideanPlane(QQ, EuclideanExactGeometry(QQ))
            True

        """
        from sage.all import QQ

        base_ring = base_ring or QQ

        if geometry is None:
            if not base_ring.is_exact():
                raise ValueError("geometry must be specified over inexact rings")

            geometry = EuclideanExactGeometry(base_ring)

        from sage.categories.all import Sets

        category = category or Sets()

        return super().__classcall__(  # pyright: ignore
            cls, base_ring=base_ring, geometry=geometry, category=category
        )

    def __init__(self, base_ring, geometry, category=None):
        r"""
        Create the Euclidean plane over ``base_ring``.

        TESTS::

            sage: from flatsurf import EuclideanPlane

            sage: TestSuite(EuclideanPlane(QQ)).run()
            sage: TestSuite(EuclideanPlane(AA)).run()

        """
        from sage.all import RR

        if geometry.base_ring() is not base_ring:
            raise ValueError(
                f"geometry base ring must be base ring of Euclidean plane but {geometry.base_ring()} is not {base_ring}"
            )

        if type(geometry.base_ring()).__name__ == "ExactReals_with_category":
            # TODO: Fix this properly in exact-real?
            pass
        elif not RR.has_coerce_map_from(geometry.base_ring()):
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

        if isinstance(other, EuclideanPlane):
            # TODO: Should we not take geometry into account?
            return self.base_ring().has_coerce_map_from(other.base_ring())

        return False

    def similarity_group(self):
        from flatsurf.geometry.similarity import SimilarityGroup
        return SimilarityGroup(self.base_ring())

    def change_ring(self, ring, geometry=None):
        r"""
        Return the Euclidean plane over a different base ``ring``.

        INPUT:

        - ``ring`` -- a ring or ``None``; if ``None``, uses the current
          :meth:`~EuclideanPlane.base_ring`.

        - ``geometry`` -- a geometry or ``None``; if ``None``; trues to convert
          the existing geometry to ``ring``.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: EuclideanPlane(QQ).change_ring(AA) is EuclideanPlane(AA)
            True

        """
        if ring is None and geometry is None:
            return self

        if ring is None:
            ring = self.base_ring()

        if geometry is None:
            geometry = self.geometry.change_ring(ring)

        return EuclideanPlane(ring, geometry)

    def _an_element_(self):
        r"""
        Return a typical point of the Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()
            sage: E.an_element()
            (0, 0)

        """
        return self.point(0, 0)

    def some_subsets(self):
        # TODO
        raise NotImplementedError

    def some_elements(self):
        r"""
        Return some representative elements, i.e., points in the plane for
        testing.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: EuclideanPlane().some_elements()
            [(0, 0), (1, 0), (0, 1), ...]

        """
        from sage.all import QQ

        return [
            self((0, 0)),
            self((1, 0)),
            self((0, 1)),
            self((-QQ(1) / 2, QQ(1) / 2)),
        ]

    def _test_some_subsets(self, tester=None, **options):
        r"""
        Run test suite on some representative subsets of the Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: EuclideanPlane()._test_some_subsets()

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
        # TODO
        raise NotImplementedError

    def __call__(self, x):
        r"""
        Return ``x`` as an element of the Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()

            sage: E((1, 0))
            (1, 0)

        We need to override this method. The normal code path in SageMath
        requires the argument to be an Element but facade sets are not
        elements::

            sage: c = E.circle((0, 0), radius=1)
            sage: Parent.__call__(E, c)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert EuclideanCircle_with_category_with_category to sage.structure.element.Element

            sage: E(c)
            { x² + y² = 1 }

        """
        if isinstance(x, EuclideanFacade):
            return self._element_constructor_(x)

        return super().__call__(x)

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as an element of the plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()

            sage: E(E.an_element()) in E
            True

        Coordinates can be converted to points::

            sage: E((1, 2))
            (1, 2)

        Elements can be converted between planes with compatible base rings::

            sage: EuclideanPlane(AA)(E((0, 0)))
            (0, 0)

        TESTS::

            sage: E(0)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot convert this element in Integer Ring to Euclidean Plane over Rational Field

        """
        from sage.all import parent

        parent = parent(x)

        if parent is self:
            return x

        if isinstance(x, EuclideanSet):
            return x.change(ring=self.base_ring(), geometry=self.geometry)

        if self.vector_space().has_coerce_map_from(parent):
            x = tuple(self.vector_space()(x))
        else:
            try:
                if len(x) == 2:
                    x = tuple(x)
                    assert len(x) == 2
            except TypeError:
                # There's no easy reliable way to detect if something is a vector in SageMath afaik, so we just try to convert it to a list of two elements.
                pass

        if isinstance(x, tuple):
            if len(x) == 2:
                return self.point(*x)

            raise ValueError("coordinate tuple must have length 2")

        raise NotImplementedError(f"cannot convert this element in {parent} to {self}")

    def base_ring(self):
        r"""
        Return the base ring over which objects in the plane are defined.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: EuclideanPlane().base_ring()
            Rational Field

        """
        return self._base_ring

    def is_exact(self):
        r"""
        Return whether subsets have exact coordinates.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.is_exact()
            True

            sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
            sage: E = EuclideanPlane(RR, geometry=EuclideanEpsilonGeometry(RR, 1e-6))
            sage: E.is_exact()
            False

        """
        return self.base_ring().is_exact()

    @cached_method
    def vector_space(self):
        r"""
        Return the two dimensional standard vector space describing vectors in
        this Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.vector_space()
            Vector space of dimension 2 over Rational Field

        """
        # TODO: Remove this method. It's unclear if this goes to the fraction field or not.
        return self.base_ring() ** 2

    def point(self, x, y, z=None):
        r"""
        Return the point in the Euclidean plane with coordinates ``x`` and
        ``y``.

        If ``z`` is provided, return the projective point with homogeneous
        coordinates ``x``, ``y`` and ``z``.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.point(1, 2)
            (1, 2)
            sage: E.point(-1, 3, 2)
            (-1/2, 3/2)
            sage: E.point(2, 0, 2)
            (1, 0)
            sage: E.point(0, 3, 2)
            (0, 3/2)
            sage: E.point(2, 2, 0)
            [1:1:0]
            sage: E.point(0, 2, 0)
            [0:1:0]
            sage: E.point(2, 0, 0)
            [1:0:0]

        ::

            sage: E.point(sqrt(2), sqrt(3))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(2) to a rational

        ::

            sage: EuclideanPlane(ZZ).point(3, 0, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: no available coordinate normalization

        """
        x = self._base_ring(x)
        y = self._base_ring(y)

        if z is None:
            z = self._base_ring.one()
        else:
            z = self._base_ring(z)
            if not z:
                if not y:
                    x = self._base_ring.one()
                else:
                    try:
                        yi = y.inverse_of_unit()
                    except Exception:
                        raise NotImplementedError("no available coordinate normalization")
                    x *= yi
                    y = self._base_ring.one()

            else:
                try:
                    zi = z.inverse_of_unit()
                except Exception:
                    raise NotImplementedError("no available coordinate normalization")
                x *= zi
                y *= zi
                z = self._base_ring.one()


        if not x and not y and not z:
            raise ValueError("invalid coordinates to create a point")

        point = self.__make_element_class__(EuclideanPoint)(self, x, y, z)
        return point

    def circle(self, center, *, radius=None, radius_squared=None, check=True):
        r"""
        Return the circle around ``center`` with ``radius`` or
        ``radius_squared``.

        INPUT:

        - ``center`` -- a point in the Euclidean plane

        - ``radius`` or ``radius_squared`` -- exactly one of the parameters
          must be specified.

        - ``check`` -- whether to verify that the arguments define a circrle in the Euclidean plane (default: ``True``)

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.circle((0, 0), radius=2)
            { x² + y² = 4 }
            sage: E.circle((0, 0), radius_squared=4)
            { x² + y² = 4 }

        A circle with radius zero is a point::

            sage: E.circle((0, 0), radius=0)
            (0, 0)

        We can explicity create a circle with radius zero by setting ``check``
        to ``False``::

            sage: E.circle((0, 0), radius=0, check=False)
            { x² + y² = 0 }

        """
        # TODO: Allow radius to be an EuclideanDistance (and convert to one.)
        if (radius is None) == (radius_squared is None):
            raise ValueError(
                "exactly one of radius or radius_squared must be specified"
            )

        if radius is not None:
            if radius < 0:
                raise ValueError("radius must not be negative")
            radius_squared = radius**2

        center = self(center)
        radius_squared = self.base_ring()(radius_squared)

        circle = self.__make_element_class__(EuclideanCircle)(
            self, center, radius_squared
        )
        if check:
            circle = circle._normalize()
            circle._check()

        return circle

    # TODO: segment should accept the more straightforward syntax segment(finite_point_start, finite_point_end)
    def segment(
        self,
        line,
        start=None,
        end=None,
        oriented=None,
        check=True,
        assume_normalized=False,
    ):
        r"""
        Return the segment on the `line`` bounded by ``start`` and ``end``.

        INPUT:

        - ``line`` -- a line in the plane

        - ``start`` -- ``None`` or a :meth:`point` on the line, e.g., obtained
          as the :meth:`EuclideanLine.intersection` of ``line`` with another
          line. If ``None``, a ray is returned, unbounded on one side.

        - ``end`` -- ``None`` or a :meth:`point` on the line, e.g., obtained
          as the :meth:`EuclideanLine.intersection` of ``line`` with another
          line. If ``None``, a ray is returned, unbounded on one side. If
          ``start`` is also ``None``, the ``line`` is returned.

        - ``oriented`` -- whether to produce an oriented segment or an
          unoriented segment. The default (``None``) is to produce an oriented
          segment iff ``line`` is oriented or both ``start`` and ``end``
          are provided so the orientation can be deduced from their order.

        - ``check`` -- boolean (default: ``True``), whether validation is
          performed on the arguments.

        - ``assume_normalized`` -- boolean (default: ``False``), if not set,
          the returned segment is normalized, i.e., if it is actually a point,
          a :class:`EuclideanPoint` is returned.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: line = E.line((0, 0), (1, 1))
            sage: E.segment(line, (0, 0), (1, 1))
            (0, 0) → (1, 1)

            sage: E.segment(line, (0, 0), (1, 1), oriented=False)
            (0, 0) — (1, 1)

        A segment that consists only of a single point gets returned as a
        point::

            sage: E.segment(line, (0, 0), (0, 0))
            (0, 0)

        To explicitly obtain a non-normalized segment in such cases, we can set
        ``assume_normalized=True``::

            sage: E.segment(line, (0, 0), (0, 0), assume_normalized=True)
            (0, 0) → (0, 0)

        A segment with a single endpoint is a ray::

            sage: E.segment(line, start=(0, 0))
            Ray from (0, 0) in direction (1, 1)
            sage: E.segment(line, end=(0, 0))
            Ray to (0, 0) from direction (1, 1)

        A segment without endpoints is a line::

            sage: E.segment(line)
            {-x + y = 0}

        .. SEEALSO::

            :meth:`EuclideanPoint.segment` to create a segment from its two
            endpoints (without specifying a line.)

        """
        line = self(line)

        if not isinstance(line, EuclideanLine):
            raise TypeError("line must be a line")

        if start is not None:
            start = self(start)
            if not isinstance(start, EuclideanPoint):
                raise TypeError("start must be a point")

        if end is not None:
            end = self(end)
            if not isinstance(end, EuclideanPoint):
                raise TypeError("end must be a point")

        if oriented is None:
            oriented = line.is_oriented() or (start is not None and end is not None)

        if not line.is_oriented():
            line = line.change(oriented=True)

            if start is None and end is None:
                # any orientation of the line will do
                pass
            elif start is None or end is None or start == end:
                raise ValueError(
                    "cannot deduce segment from single endpoint on an unoriented line"
                )
            elif line.parametrize(start, check=False) > line.parametrize(
                end, check=False
            ):
                line = -line

        segment = self.__make_element_class__(
            EuclideanOrientedSegment if oriented else EuclideanUnorientedSegment
        )(self, line, start, end)

        if check:
            segment._check(require_normalized=False)

        if not assume_normalized:
            segment = segment._normalize()

        if check:
            segment._check(require_normalized=True)

        return segment

    def line(self, a, b, c=None, oriented=True, check=True):
        r"""
        Return a line in the Euclidean plane.

        If only ``a`` and ``b`` are given, return the line going through the
        points ``a`` and then ``b``.

        If ``c`` is specified, return the line given by the equation

        .. MATH::

            a + bx + cy = 0

        oriented such that the half plane

        .. MATH::

            a + bx + cy \ge 0

        is to its left.

        INPUT:

        - ``a`` -- a point or an element of the :meth:`base_ring`

        - ``b`` -- a point or an element of the :meth:`base_ring`

        - ``c`` -- ``None`` or an element of the :meth:`base_ring` (default: ``None``)

        - ``oriented`` -- whether the returned line is oriented (default: ``True``)

        - ``check`` -- whether to verify that the arguments actually define a
          line (default: ``True``)

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: E.line((0, 0), (1, 1))
            {-x + y = 0}

            sage: E.line(0, -1, 1)
            {-x + y = 0}

        """
        if c is None:
            a = self(a)
            b = self(b)

            if a == b:
                raise ValueError("points specifying a line must be distinct")

            ax, ay = a
            bx, by = b

            C = bx - ax
            B = ay - by
            A = -(B * ax + C * ay)

            return self.line(A, B, C, oriented=oriented, check=check)

        a = self.base_ring()(a)
        b = self.base_ring()(b)
        c = self.base_ring()(c)

        line = self.__make_element_class__(
            EuclideanOrientedLine if oriented else EuclideanUnorientedLine
        )(self, a, b, c)

        if check:
            line = line._normalize()
            line._check()

        return line

    geodesic = line

    def ray(self, base, direction, check=True):
        r"""
        Return a ray from ``base`` in ``direction``.

        INPUT:

        - ``base`` -- a point in the Euclidean plane

        - ``direction`` -- a vector in the two dimensional space over the
          :meth:`base_ring`, see :meth:`vector_space`

        - ``check`` -- a boolean (default: ``True``); whether to validate the
          parameters

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: ray = E.ray((0, 0), (1, 1))

        The base point is contained in the ray::

            sage: E.point(0, 0) in ray
            True

        The direction must be non-zero::

            sage: E.ray((0, 0), (0, 0))
            Traceback (most recent call last):
            ...
            ValueError: direction must be distinguishable from the zero vector

        """
        base = self(base)
        direction = self.vector_space()(direction)

        if self.geometry._zero(direction):
            raise ValueError("direction must be distinguishable from the zero vector")

        if not isinstance(base, EuclideanPoint):
            raise TypeError("base must be a point")

        line = self.line(base, base.translate(direction), check=check)
        return self.segment(line, base, None, check=check)

    def polygon(
        self,
        vertices=None,
        edges=None,
        angles=None,
        lengths=None,
        check=True,
        category=None,
    ):
        r"""
        Return a polygon from the given ``vertices``, ``edges``, or ``angles``.

        INPUT:

        - ``vertices`` -- a sequence of vertices or ``None`` (default:
          ``None``); the vertices of the polygon as points in the Euclidean
          plane

        # TODO: Allow edges to be segments and lines and not just vectors.

        - ``edges`` -- a sequence of vectors, segments, rays, or lines, or
          ``None`` (default: ``None``); the paths connecting the vertices of
          the polygon

        - ``angles`` -- a sequence of numbers that prescribe the inner angles
          of the polygon or ``None`` (default: ``None``); the angles are
          rescaled so that their sum matches the sum of the angles in an ngon.

        # TODO: Allow lengths to be norms?

        - ``lengths`` -- a sequence of numbers that prescribe the lengths of
          the edges of the polygon or ``None`` (default: ``None``)

        # TODO: Replace with the standard meaning of check in the EuclideanPlane.
        - ``check`` -- a boolean (default: ``True``); whether to check the
          consistency of the parameters or blindly trust them. Setting this to
          ``False`` allows creation of degenerate polygons in some cases. While
          they might be somewhat functional, no guarantees are made about such
          polygons.

        # TODO: Document category.

        EXAMPLES:

        A right triangle::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(QQ)
            sage: E.polygon(vertices=[(0, 0), (1, 0), (0, 1)])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        A right triangle that is not based at the origin::

            sage: E.polygon(vertices=[(1, 0), (2, 0), (1, 1)])
            Polygon(vertices=[(1, 0), (2, 0), (1, 1)])

        Vertices can also be actual points of the Euclidean plane::

            sage: E.polygon(vertices=[E.point(1, 0), (-1, 1), E.point(0, -1)])
            Polygon(vertices=[(1, 0), (-1, 1), (0, -1)])

        Edges can also be segments in the Euclidean plane::

            sage: E.polygon(edges=[
            ....:     E.point(0, 0).segment(E.point(1, 0)),
            ....:     E.point(1, 0).segment(E.point(0, 1)),
            ....:     E.point(0, 1).segment(E.point(0, 0))])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        Or a mix of segments and vectors::

            sage: E.polygon(edges=[
            ....:     E.point(0, 0).segment(E.point(1, 0)),
            ....:     (-1, 1),
            ....:     E.point(0, 1).segment(E.point(0, 0))])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        A right triangle at the origin, specified by giving the edge vectors::

            sage: E.polygon(edges=[(1, 0), (-1, 1), (0, -1)])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        When redundant information is given, it is checked for consistency::

            sage: E.polygon(vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
            sage: E.polygon(vertices=[(1, 0), (2, 0), (1, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
            Polygon(vertices=[(1, 0), (2, 0), (1, 1)])
            sage: E.polygon(vertices=[(0, 0), (2, 0), (1, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
            Traceback (most recent call last):
            ...
            ValueError: vertices and edges are not compatible

        Polygons given by edges must be closed (in particular we do not add an
        edge automatically to close things up since this is often not what the
        user wanted)::

            sage: E.polygon(edges=[(1, 0), (0, 1), (1, 1)])
            Traceback (most recent call last):
            ...
            ValueError: polygon is not closed

        A polygon with prescribed angles::

            sage: E.polygon(angles=[2, 1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        Again, if vertices and edges are also specified, they must be compatible
        with the angles::

            sage: E.polygon(angles=[2, 1, 1], vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

            sage: E.polygon(angles=[1, 2, 3], vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
            Traceback (most recent call last):
            ...
            ValueError: polygon does not have the prescribed angles

        When angles are specified, side lengths can also be prescribed::

            sage: E = EuclideanPlane(QuadraticField(3))
            sage: E.polygon(angles=[1, 1, 1], lengths=[1, 1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*a)])

        The function will deduce lengths if one or two are missing::

            sage: E.polygon(angles=[1, 1, 1, 1], lengths=[1, 1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            sage: E.polygon(angles=[1, 1, 1, 1], lengths=[1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            sage: E.polygon(angles=[1, 1, 1, 1], lengths=[1])
            Traceback (most recent call last):
            ...
            NotImplementedError: could only construct a digon from this data but not a quadrilateral

        Equally, we deduce vertices or edges::

            sage: E.polygon(angles=[1, 1, 1, 1], vertices=[(0, 0), (1, 0), (1, 1)])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            sage: E.polygon(angles=[1, 1, 1, 1], edges=[(1, 0), (0, 1)])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        When the angles are incompatible with the data, an error is reported (that
        might be somewhat cryptic at times)::

            sage: E.polygon(angles=[1, 1, 1, 1], edges=[(1, 0), (0, 1), (1, 2)])
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot recover a rational angle from these numerical results

        When lengths are given in addition to vertices or edges, they are checked for consistency::

            sage: E.polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)], lengths=[1, 1, 1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

            sage: E.polygon(vertices=[(0, 0), (1, 0), (0, 1)], lengths=[1, 1, 1])
            Traceback (most recent call last):
            ...
            ValueError: polygon does not have the prescribed lengths

        Currently, we cannot create a polygon from just lengths::

            sage: E.polygon(lengths=[1, 1, 1])
            Traceback (most recent call last):
            ...
            NotImplementedError: one of vertices, edges, or angles must be set

        Polygons do not have to be convex::

            sage: p = E.polygon(vertices=[(0, 0), (1, 1), (2, 0), (2, 4), (0, 4)])
            sage: p.describe_polygon()
            ('a', 'non-convex pentagon', 'non-convex pentagons')

        Polygons must be positively oriented::

            sage: E.polygon(vertices=[(0, 0), (0, 1), (1, 0)])
            Traceback (most recent call last):
            ...
            ValueError: polygon has negative area; probably the vertices are not in counter-clockwise order

        Polygons must have at least three sides::

            sage: E.polygon(vertices=[(0, 0), (1, 0)])  # TODO: This should be an unoriented segment
            Traceback (most recent call last):
            ...
            ValueError: polygon has zero area

            sage: E.polygon(vertices=[(0, 0), (1, 0), (2, 0)])  # TODO: Should this be supported? What do we do in the hyperbolic world?
            Traceback (most recent call last):
            ...
            ValueError: polygon has zero area

        Currently, polygons must not self-intersect::

            sage: p = E.polygon(vertices=[(0, 0), (2, 0), (0, 1), (1, -1), (2, 1)])  # not tested  # TODO: Enable this test again.
            Traceback (most recent call last):
            ...
            NotImplementedError: polygon self-intersects

        Currently, all angles must be less than 2π::

            sage: p = E.polygon(angles=[14, 1, 1, 1, 1])
            Traceback (most recent call last):
            ...
            NotImplementedError: each angle must be in (0, 2π)

        """
        # Determine the number of sides of this polygon.
        if angles:
            n = len(angles)
        elif edges:
            n = len(edges)
        elif vertices:
            n = len(vertices)
        else:
            raise NotImplementedError("one of vertices, edges, or angles must be set")

        choice = False

        # Determine the category of the polygon
        choice, category = self._polygon_category(
            n=n, angles=angles, choice=choice
        )

        # Write the vertices as points in the Euclidean plane
        choice, vertices = self._polygon_normalize_vertices(
            n=n, vertices=vertices, choice=choice
        )

        # Write the edges as segments, rays, and lines
        choice, edges = self._polygon_normalize_edges(
            n=n, vertices=vertices, edges=edges, choice=choice
        )

        choice, angles = self._polygon_normalize_angles(
            n=n, angles=angles, choice=choice
        )

        choice, lengths = self._polygon_normalize_lengths(
            n=n, lengths=lengths, choice=choice
        )

        # Deduce edges from vertices
        choice, edges, vertices = self._polygon_edges_from_vertices(
            n=n, vertices=vertices, edges=edges, choice=choice
        )

        # Deduce edges from angles (and possibly lengths)
        choice, edges, vertices, angles, lengths = self._polygon_edges_from_angles(
            n=n, edges=edges, vertices=vertices, angles=angles, lengths=lengths, choice=choice
        )

        # Add an edge if the polygon is not closed.
        choice, edges = self._polygon_edges_close(
            n=n, edges=edges, choice=choice
        )

        polygon = EuclideanPolygon(self, edges=tuple(edges), category=category)

        if check:
            polygon = polygon._normalize()
            self._polygon_check_n(polygon, n=n)
            polygon._check()
            self._polygon_check_vertices(polygon, vertices=vertices)
            self._polygon_check_edges(polygon, edges=edges)
            self._polygon_check_angles(polygon, angles=angles)
            self._polygon_check_lengths(polygon, lengths=lengths)

        return polygon

    def rectangle(self, width, height, **kwargs):
        from flatsurf.geometry.polygon import polygons
        return polygons.rectangle(width, height, parent=self, **kwargs)

    def square(self, side, **kwargs):
        from flatsurf.geometry.polygon import polygons
        return polygons.square(side, parent=self, **kwargs)

    # TODO: Add triangle
    # TODO: Add regular_ngon
    # TODO: Add right_triangle

    # TODO: All these helpers shoud probably go into the EuclideanPolygon class.

    def _polygon_check_n(self, polygon, n):
        if len(polygon.sides()) != n:
            from flatsurf.geometry.categories.polygons import Polygons
            ngon = Polygons._describe_polygon(n)
            mgon = Polygons._describe_polygon(len(polygon.sides()))
            raise NotImplementedError(f"could only construct {mgon[0]} {mgon[1]} from this data but not {ngon[0]} {ngon[1]}")

    def _polygon_check_vertices(self, polygon, vertices):
        if not vertices:
            return

        if polygon.vertices() != tuple(v.vector() for v in vertices):
            raise ValueError("vertices and edges are not compatible")

    def _polygon_check_edges(self, polygon, edges):
        if not edges:
            return

        if polygon.sides() != tuple(edges):
            raise ValueError

    def _polygon_check_angles(self, polygon, angles):
        if not angles:
            return

        # Check that the polygon has the prescribed angles
        from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
            EuclideanPolygonsWithAngles,
        )
        from flatsurf.geometry.categories.euclidean_polygons import (
            EuclideanPolygons,
        )

        # Use EuclideanPolygon's angle() so we do not use the precomputed angles set by the category.
        from flatsurf.geometry.categories.euclidean_polygons_with_angles import EuclideanPolygonsWithAngles
        if EuclideanPolygonsWithAngles._normalize_angles(angles) != tuple(
            EuclideanPolygons.ParentMethods.angle(polygon, i)
            for i in range(len(polygon.vertices()))
        ):
            raise ValueError("polygon does not have the prescribed angles")

    def _polygon_check_lengths(self, polygon, lengths):
        if not lengths:
            return

        for edge, length in zip(polygon.edges(), lengths):
            if edge.norm() != length:
                raise ValueError("polygon does not have the prescribed lengths")

    def _polygon_normalize_vertices(self, n: int, vertices, choice):
        if vertices is None:
            return choice, vertices

        vertices = [self(v) for v in vertices]
        return choice, vertices

    def _polygon_normalize_edges(self, n: int, vertices: List["EuclideanPoint"], edges, choice):
        if edges is None:
            return choice, edges

        edges = list(edges)

        pos = None
        if vertices:
            pos = vertices[0]
        for e in range(len(edges)):
            if not isinstance(edges[e], EuclideanFacade):
                # the edge is a vector
                if pos is None and e == 0:
                    choice = True
                    pos = self.point(0, 0)

                if pos is None:
                    raise ValueError("cannot use edge vector when there is no base point for an edge")

                edges[e] = pos.segment(pos.translate(edges[e]))

            if not edges[e].is_oriented():
                raise ValueError("edges must be oriented segments, rays, or lines")

            if isinstance(edges[e], EuclideanSegment):
                pos = edges[e].end()
            else:
                pos = None

        return choice, edges

    def _polygon_normalize_angles(self, n, angles, choice):
        if angles is not None:
            angles = list(angles)

        return choice, angles

    def _polygon_normalize_lengths(self, n, lengths, choice):
        if lengths is not None:
            lengths = list(lengths)

        return choice, lengths

    def _polygon_edges_from_vertices(self, n, vertices, edges, choice):
        if edges is not None:
            return choice, edges, vertices

        if vertices is None:
            return choice, edges, vertices

        # Intentionally, we omit the edge joining the last to the first vertex.
        # There might be two vertices missing if angles have been specified.
        # Otherwise _polygon_edges_close() is going to add that edge later.
        edges = [vertices[i].segment(vertices[i + 1]) for i in range(len(vertices) - 1)]
        return choice, edges, None

    def _polygon_edges_from_angles(self, n, edges, vertices, angles, lengths, choice):
        # TODO: Untangle this mess
        if edges is None:
            if not angles:
                return choice, edges, vertices, angles, lengths

            from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons
            category = EuclideanPolygons(self.base_ring()).Simple().WithAngles(angles)

            if lengths:
                edges = []
                for slope, length in zip(category.slopes(), lengths):
                    scale = self.base_ring()((length**2 / (slope[0] ** 2 + slope[1] ** 2)).sqrt())
                    edges.append(scale * slope)

                if len(edges) == n:
                    angles = None

                lengths = None
            else:
                choice = True

                # We pick the edges such that they form a closed polygon with the
                # prescribed angles. However, there might be self-intersection which
                # currently leads to an error.
                edges = [
                    length * slope
                    for (length, slope) in zip(
                        sum(r.vector() for r in category.lengths_polytope().rays()),
                        category.slopes(),
                    )
                ]

                angles = None

            choice, edges = self._polygon_normalize_edges(n, vertices, edges, choice)
        
        if angles and len(angles) == n and len(edges) == n - 2:
            # We do not use category.slopes() since the matrix formed by such
            # slopes might not be invertible (because exact-reals do not have a
            # fraction field implemented).
            from flatsurf.geometry.polygon import EuclideanPolygonsWithAngles
            slopes = EuclideanPolygonsWithAngles(angles).slopes()

            # We do not use solve_left() because the vertices might not live in
            # a ring that has a fraction field implemented (such as an
            # exact-real ring).
            from sage.all import matrix
            s, t = (edges[0].start().vector() - edges[-1].end().vector()) * matrix(
                [slopes[-1], slopes[n - 2]]
            ).inverse()
            assert edges[0].start().vector() - s * slopes[-1] == edges[-1].end().vector() + t * slopes[n - 2]

            if s <= 0 or t <= 0:
                raise (NotImplementedError if choice else ValueError)(
                    "cannot determine polygon with these angles from the given data"
                )

            edges.append(edges[-1].end().segment(edges[0].start().vector() - s * slopes[-1]))

        return choice, edges, vertices, angles, lengths

    def _polygon_edges_close(self, n, edges, choice):
        assert edges is not None

        if len(edges) < n and edges[-1].end() != edges[0].start():
            edges.append(edges[-1].end().segment(edges[0].start()))

        return choice, edges

    def _polygon_category(self, n, angles, choice):
        from flatsurf.geometry.categories import EuclideanPolygons

        # Currently, all polygons are assumed to be without self-intersection, i.e., simple.
        from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons
        category = EuclideanPolygons(self.base_ring()).Simple()
        if angles:
            category = category.WithAngles(angles)

        if n == 3:
            category = category.Convex()

        return choice, category

    def _polygon_normalize_arguments(self, category, n, vertices, edges, angles, lengths):
        r"""
        Return the normalized arguments defining a polygon. Additionally, a flag is
        returned that indicates whether we made a choice in normalizing these
        arguments.

        This is a helper function for :func:`Polygon`.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import _Polygon_normalize_arguments
            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: category = EuclideanPolygons(AA)
            sage: _Polygon_normalize_arguments(category=category, n=3, vertices=[(0, 0), (1, 0), (0, 1)], edges=None, angles=None, lengths=None)
            (False, [(0, 0), (1, 0), (0, 1)], None, None, None)
            sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=[(1, 0), (-1, 1), (0, -1)], angles=None, lengths=None)
            (False, [(0, 0), (1, 0), (0, 1)], None, None, None)

            sage: category = category.WithAngles([1, 1, 1])
            sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=None, angles=[1, 1, 1], lengths=None)
            (True, [(0, 0), (1, 0), (1/2, 0.866025403784439?)], None, None, None)
            sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=None, angles=[1, 1, 1], lengths=[AA(2).sqrt(), 1])
            (False,
             [(0, 0), (1.414213562373095?, 0), (0.9142135623730951?, 0.866025403784439?)],
             None,
             [1, 1, 1],
             None)

        """
        # TODO: Delete. Unused.
        base_ring = category.base_ring()

        # Track whether we made a choice that possibly is the reason that we fail
        # to find a polygon with the given data.
        choice = False

        # Rewrite angles and lengths as angles and edges.
        if angles and lengths and not edges:
            edges = []
            for slope, length in zip(category.slopes(), lengths):
                scale = base_ring((length**2 / (slope[0] ** 2 + slope[1] ** 2)).sqrt())
                edges.append(scale * slope)

            if len(edges) == n:
                angles = 0

            lengths = None

        # Deduce edges if only angles are given
        if angles and not edges and not vertices:
            assert not lengths

            choice = True

            # We pick the edges such that they form a closed polygon with the
            # prescribed angles. However, there might be self-intersection which
            # currently leads to an error.
            edges = [
                length * slope
                for (length, slope) in zip(
                    sum(r.vector() for r in category.lengths_polytope().rays()),
                    category.slopes(),
                )
            ]

            angles = None

        if vertices:
            vertices = [self.point(*v) for v in vertices]

        # Rewrite vertices as edges
        if vertices and not edges:
            edges = [v.segment(w) for v, w in zip(vertices, vertices[1:])]
            vertices = None

        # TODO: Make sure that everything are lists until the very end when we turn everything into tuples.
        edges = list(edges)

        # Rewrite edges as segments (and lines)
        pos = None
        if vertices:
            pos = self(vertices[0])
        for e in range(len(edges)):
            if not isinstance(edges[e], EuclideanFacade):
                # edge is a vector, turn it into a segment
                if pos is None:
                    choice = True
                    pos = self.point(0, 0)
                
                edges[e] = pos.segment(pos.translate(edges[e]))

            if not edges[e].is_oriented():
                raise ValueError("edges must be oriented segments, rays, or lines")

            pos = edges[e].end()

        return choice, vertices, edges, angles, lengths

    def _polygon_check(self, polygon, vertices, edges, angles, lengths):
        r"""
        Verify that ``p`` is a valid polygon and that it satisfies the constraints
        given.

        This is a helper function for :func:`Polygon`.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import _Polygon_check, Polygon
            sage: p = Polygon(angles=[1, 1, 1])
            sage: _Polygon_check(p, vertices=None, edges=None, angles=[1, 1, 1], lengths=None, convex=None)

        """
        # TODO: Unused. Make sure that everything is still in the new check functions.
        # Check that the polygon satisfies the assumptions of EuclideanPolygon
        area = polygon.area()

        if area < 0:
            raise ValueError(
                "polygon has negative area; probably the vertices are not in counter-clockwise order"
            )

        if area == 0:
            raise ValueError("polygon has zero area")

        if any(edge == 0 for edge in polygon.edges()):
            raise ValueError("polygon has zero edge")

        for i in range(len(polygon.edges())):
            from flatsurf.geometry.euclidean import is_anti_parallel

            if is_anti_parallel(polygon.edge(i), polygon.edge(i + 1)):
                raise ValueError("polygon has anti-parallel edges")

        from flatsurf.geometry.categories import EuclideanPolygons

        if not EuclideanPolygons.ParentMethods.is_simple(polygon):
            raise NotImplementedError("polygon self-intersects")

        # Check that any redundant data is compatible
        if vertices:
            if len(edges) != len(vertices):
                raise ValueError("vertices and edges must have the same length")

            for v, e in zip(vertices, polygon.sides()):
                if v != e.start():
                    raise ValueError("vertices and edges are not compatible")

        if angles:
            # Check that the polygon has the prescribed angles
            from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
                EuclideanPolygonsWithAngles,
            )
            from flatsurf.geometry.categories.euclidean_polygons import (
                EuclideanPolygons,
            )

            # Use EuclideanPolygon's angle() so we do not use the precomputed angles set by the category.
            from flatsurf.geometry.categories.euclidean_polygons_with_angles import EuclideanPolygonsWithAngles
            if EuclideanPolygonsWithAngles._normalize_angles(angles) != tuple(
                EuclideanPolygons.ParentMethods.angle(polygon, i)
                for i in range(len(polygon.vertices()))
            ):
                raise ValueError("polygon does not have the prescribed angles")

        if lengths:
            for edge, length in zip(polygon.edges(), lengths):
                if edge.norm() != length:
                    raise ValueError("polygon does not have the prescribed lengths")

    def _polygon_complete_edges(self, n, edges, angles, choice):
        r"""
        Return edges that define a polygon by completing the ``edges`` to an
        ``n``-gon with ``angles``.

        This is a helper function for :func:`Polygon`.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(AA)

            sage: E._polygon_complete_edges(3, [E((0, 0)).segment((1, 0))], [1, 1, 1], choice=False)
            [(0, 0) → (1, 0),
             (1, 0) → (1/2, 0.866025403784439?),
             (1/2, 0.866025403784439?) → (0, 0)]

        """
        if angles:
            if len(angles) == n and  len(edges) == n - 2:
                # We do not use category.slopes() since the matrix formed by such
                # slopes might not be invertible (because exact-reals do not have a
                # fraction field implemented).
                from flatsurf.geometry.polygon import EuclideanPolygonsWithAngles
                slopes = EuclideanPolygonsWithAngles(angles).slopes()

                # We do not use solve_left() because the vertices might not live in
                # a ring that has a fraction field implemented (such as an
                # exact-real ring).
                from sage.all import matrix
                s, t = (edges[0].start().vector() - edges[-1].end().vector()) * matrix(
                    [slopes[-1], slopes[n - 2]]
                ).inverse()
                assert edges[0].start().vector() - s * slopes[-1] == edges[-1].end().vector() + t * slopes[n - 2]

                if s <= 0 or t <= 0:
                    raise (NotImplementedError if choice else ValueError)(
                        "cannot determine polygon with these angles from the given data"
                    )

                edges.append(edges[-1].end().segment(edges[0].start().vector() - s * slopes[-1]))
            if len(edges) < n - 1:
                from flatsurf.geometry.categories import Polygons

                raise NotImplementedError(
                    f"cannot construct {' '.join(Polygons._describe_polygon(n)[:2])} from {n} angles and {len(edges)} edges"
                )

        if len(edges) == n - 1:
            edges.append(edges[-1].end().segment(edges[0].start()))

        return edges

    @cached_method
    def norm(self):
        r"""
        Return the Euclidean norm on this plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: norm = E.norm()
            sage: norm(2)
            2
            sage: norm.from_norm_squared(2)
            1.414213562373095?
            sage: norm.from_vector((1, 1))
            1.414213562373095?

        """
        return EuclideanDistances(self)

    def empty_set(self):
        # TODO
        raise NotImplementedError

    def _repr_(self):
        r"""
        Return a printable representation of this Euclidean plane.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: EuclideanPlane(AA)
            Euclidean Plane over Algebraic Real Field

        """
        return f"Euclidean Plane over {repr(self.base_ring())}"


class EuclideanGeometry(Geometry):
    r"""
    Predicates and primitive geometric constructions over a base ``ring``.

    This class and its subclasses implement the core underlying Euclidean
    geometry that depends on the base ring. For example, when deciding whether
    two points in the plane are equal, we cannot just compare their coordinates
    if the base ring is inexact. Therefore, that predicate is implemented in
    this "geometry" class and is implemented differently by
    :class:`EuclideanExactGeometry` for exact and
    :class:`EuclideanEpsilonGeometry` for inexact rings.

    INPUT:

    - ``ring`` -- a ring, the ring in which coordinates in the Euclidean plane
      will be represented

    .. NOTE::

        Abstract methods are not marked with `@abstractmethod` since we cannot
        use the ABCMeta metaclass to enforce their implementation; otherwise,
        our subclasses could not use the unique representation metaclasses.

    EXAMPLES:

    The specific Euclidean geometry implementation is picked automatically,
    depending on whether the base ring is exact or not::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: E.geometry
        Exact geometry over Rational Field
        sage: E((0, 0)) == E((1/1024, 0))
        False

    However, we can explicitly use a different or custom geometry::

        sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
        sage: E = EuclideanPlane(QQ, EuclideanEpsilonGeometry(QQ, 1/1024))
        sage: E.geometry
        Epsilon geometry with ϵ=1/1024 over Rational Field
        sage: E((0, 0)) == E((1/2048, 0))
        True

    .. SEEALSO::

        :class:`EuclideanExactGeometry`, :class:`EuclideanEpsilonGeometry`
    """

    def _equal_vector(self, v, w):
        r"""
        Return whether the vectors ``v`` and ``w`` should be considered equal.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        INPUT:

        - ``v`` -- a vector over the :meth:`base_ring`

        - ``w`` -- a vector over the :meth:`base_ring`

        EXAMPLES::

            sage: from flatsurf.geometry.euclidean import EuclideanExactGeometry
            sage: G = EuclideanExactGeometry(QQ)
            sage: G._equal_vector((0, 0), (0, 0))
            True
            sage: G._equal_vector((0, 0), (0, 1/1024))
            False

        """
        if len(v) != len(w):
            raise TypeError("v and w must be vectors in the same vector space")

        return all(self._equal(vv, ww) for (vv, ww) in zip(v, w))

    def _equal_point(self, p, q):
        r"""
        Return whether the points ``p`` and ``q`` should be considered equal.

        .. NOTE::

            This predicate should not be used directly in geometric
            constructions since it does not specify the context in which this
            question is asked. This makes it very difficult to override a
            specific aspect in a custom geometry.

        INPUT:

        - ``p`` -- a point in the Euclidean plane

        - ``q`` -- a point in the Euclidean plane

        EXAMPLES::

            sage: from flatsurf.geometry.euclidean import EuclideanPlane, EuclideanExactGeometry
            sage: E = EuclideanPlane(QQ)
            sage: G = EuclideanExactGeometry(QQ)
            sage: G._equal_point(E.point(0, 0), E.point(0, 0))
            True
            sage: G._equal_point(E.point(0, 0), E.point(0, 1/1024))
            False

        """
        return self._equal_vector((p._x, p._y, p._z), (q._x, q._y, q._z))


class EuclideanExactGeometry(UniqueRepresentation, EuclideanGeometry, ExactGeometry):
    r"""
    Predicates and primitive geometric constructions over an exact base ring.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: E.geometry
        Exact geometry over Rational Field

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanExactGeometry
        sage: isinstance(E.geometry, EuclideanExactGeometry)
        True

    .. SEEALSO::

        :class:`EuclideanEpsilonGeometry` for an implementation over inexact rings

    """

    def change_ring(self, ring):
        r"""
        Return this geometry with the :meth:`~EuclideanGeometry.base_ring`
        changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.geometry.change_ring(QQ) == E.geometry
            True
            sage: E.geometry.change_ring(AA)
            Exact geometry over Algebraic Real Field

        """
        if not ring.is_exact():
            raise ValueError("cannot change_ring() to an inexact ring")

        return EuclideanExactGeometry(ring)


class EuclideanEpsilonGeometry(
    UniqueRepresentation, EuclideanGeometry, EpsilonGeometry
):
    r"""
    Predicates and primitive geometric constructions over an inexact base ring.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
        sage: E = EuclideanPlane(RR, geometry=EuclideanEpsilonGeometry(RR, 1e-6))
        sage: E.geometry
        Epsilon geometry with ϵ=1.00000000000000e-6 over Real Field with 53 bits of precision

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
        sage: isinstance(E.geometry, EuclideanEpsilonGeometry)
        True

    .. SEEALSO::

        :class:`EuclideanExactGeometry` for an implementation over exact rings

    """

    def _equal_vector(self, v, w):
        r"""
        Return whether the vectors ``v`` and ``w`` should be considered equal.

        Implements :meth:`EuclideanGeometry._equal_vector` by comparing the
        Euclidean distance of the points to this geometry's epsilon.

        EXAMPLES::

            sage: from flatsurf.geometry.euclidean import EuclideanPlane, EuclideanEpsilonGeometry
            sage: G = EuclideanEpsilonGeometry(RR, 1e-3)
            sage: E = EuclideanPlane(RR, geometry=G)

            sage: G._equal_point(E.point(0, 0), E.point(0, 0))
            True
            sage: G._equal_point(E.point(0, 0), E.point(0, 1/1024))
            True

            sage: G._equal_vector((0, 0), (0, 0))
            True
            sage: G._equal_vector((0, 0), (0, 1/1024))
            True

        """
        if len(v) != len(w):
            raise TypeError("vectors must have same length")

        return sum((vv - ww) ** 2 for (vv, ww) in zip(v, w)) < self._epsilon**2

    def change_ring(self, ring):
        r"""
        Return this geometry with the :meth:`~EuclideanGeometry.base_ring`
        changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
            sage: G = EuclideanEpsilonGeometry(RR, 1e-3)
            sage: G.change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: cannot change_ring() to an exact ring
            sage: G.change_ring(RDF)
            Epsilon geometry with ϵ=0.001 over Real Double Field

        """
        if ring.is_exact():
            raise ValueError("cannot change_ring() to an exact ring")

        return EuclideanEpsilonGeometry(ring, self._epsilon)


class EuclideanSet(SageObject):
    r"""
    Base class for subsets of :class:`EuclideanPlane`.

    .. NOTE::

        Concrete subclasses should apply the following rules.

        There should only be a single type to describe a certain subset:
        normally, a certain subset, say a point, should only be described by a
        single class, namely :class:`Point`. Of course, one could
        describe a point as a polygon delimited by some edges that all
        intersect in that single point, such objects should be avoided. Namely,
        the methods that create a subset, say :meth:`EuclideanPlane.polygon`
        take care of this by calling a sets
        :meth:`EuclideanSet._normalize` to rewrite a set in its most natural
        representation. To get the denormalized representation, we can always
        set `check=False` when creating the object. For this to work, the
        `__init__` should not take care of any such normalization and accept
        any input that can possibly be made sense of.

        Comparison with ``==`` should mean "is essentially indistinguishable
        from": Implementing == to mean anything else would get us into trouble
        in the long run. In particular we cannot implement <= to mean "is
        subset of" since then an oriented and an unoriented geodesic would be
        `==`. So, objects of a different type should almost never be equal. A
        notable exception are objects that are indistinguishable to the end
        user but use different implementations.

    TESTS::

        sage: from flatsurf import EuclideanPlane
        sage: from flatsurf.geometry.euclidean import EuclideanSet
        sage: E = EuclideanPlane()

        sage: isinstance(E((0, 0)), EuclideanSet)
        True

    """

    def _check(self, require_normalized=True):
        r"""
        Validate this convex set.

        Subclasses run specific checks here that can be disabled when creating
        objects with ``check=False``.

        INPUT:

        - ``require_normalized`` -- a boolean (default: ``True``); whether to
          include checks that assume that normalization has already happened

        EXAMPLES:

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: P = E.point(0, 0)
            sage: P._check()

        """
        pass

    def _normalize(self):
        r"""
        Return this set possibly rewritten in a simpler form.

        This method is only relevant for sets created with ``check=False``.
        Such sets might have been created in a non-canonical way, e.g., when
        creating a :class:`OrientedSegment` whose start and end point is
        identical.

        """
        return self

    def _test_normalize(self, **options):
        r"""
        Verify that normalization is idempotent.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: segment = E((0, 0)).segment((1, 0))
            sage: segment._test_normalize()

        """
        tester = self._tester(**options)

        normalization = self._normalize()

        tester.assertEqual(normalization, normalization._normalize())

    def __contains__(self, point):
        r"""
        Return whether this set contains the point ``point``.

        INPUT:

        - ``point`` -- a point in the Euclidean plane

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: c = E.circle((0, 0), radius=1)
            sage: E((0, 0)) in c
            False

        """
        raise NotImplementedError(
            "this subset of the Euclidean plane cannot decide whether it contains a given point yet"
        )

    # TODO: Add a _test_contains test.

    # TODO: Add is_bounded() and a test method.

    def change(self, *, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this set.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~EuclideanPlane.base_ring`); the ring over which the new set
          will be defined.

        - ``geometry`` -- a :class:`EuclideanGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new set.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness) whether the new set will be explicitly oriented.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: segment = E((0, 0)).segment((1, 1))

        We can change the base ring over which this set is defined::

            sage: segment.change(ring=AA)
            (0, 0) → (1, 1)

        We can drop the explicit orientation of a set::

            sage: unoriented = segment.change(oriented=False)
            sage: unoriented.is_oriented()
            False

        We can also take an unoriented set and pick an orientation::

            sage: oriented = unoriented.change(oriented=True)
            sage: oriented.is_oriented()
            True

        .. SEEALSO::

            :meth:`is_oriented` to determine whether a set is oriented.

        """
        raise NotImplementedError(f"this {type(self)} does not implement change()")

    def change_ring(self, ring):
        r"""
        Return this set as an element of the Euclidean plane over ``ring``.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: p = E((0, 0))
            sage: p.change_ring(AA)
            (0, 0)

        """
        return self.change(ring=ring)

    # TODO: Add change_ring, change and test methods.

    # TODO: Add plot and test_plot()

    def _acted_upon_(self, g, self_on_left):
        r"""
        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: p = E.point(1, 1)

            sage: matrix(ZZ, 2, 2, [2, 3, 5, 7]) * p
            (5, 12)
            sage: matrix(ZZ, 2, 3, [0, 0, 3, 1, 1, 1]) * p
            (3, 3)
            sage: matrix(ZZ, 3, 3, [1, 2, 3, 1, 1, -1, 2, 0, 1]) * p
            (2, 1/3)
        """
        # NOTE: it might have been preferable to implement a proper action on subsets
        # of the hyperbolic plane via sage.categories.action.Action. However, it does
        # not work in our context since subsets are (facade) parents.
        from sage.structure.element import parent
        from sage.matrix.matrix_space import MatrixSpace
        E = self.parent()
        R = E.base_ring()
        P = parent(g)

        if P is R:
            return self._apply_scalar(g)

        if P is E.similarity_group():
            return self._apply_similarity(g)

        M2 = MatrixSpace(R, 2)
        if P is M2:
            return self._apply_2x2_matrix(g)

        M3 = MatrixSpace(R, 3)
        if P is M3:
            # NOTE: call the method without underscore that performs a sanity
            # check on orientation
            return self.apply_3x3_matrix(g)

        M23 = MatrixSpace(R, 2, 3)
        if P is M23:
            return self.apply_3x3_matrix(to_3x3_matrix(g))

        if R.has_coerce_map_from(P):
            return self._apply_scalar(R(g))

        if M2.has_coerce_map_from(P):
            return self._apply_2x2_matrix(M2(g))

        if M3.has_coerce_map_from(P):
            # NOTE: call the method without underscore that performs a sanity
            # check on orientation
            return self._apply_3x3_matrix(M3(g))

        if M23.has_coerce_map_from(P):
            return self._apply_3x3_matrix(to_3x3_matrix(M23(g)))

    def _apply_scalar(self, r):
        raise NotImplementedError

    # TODO: move TESTS as pytest
    # TODO: carefully implement similarity reversing orientations and tests
    def apply_similarity(self, g):
        r"""
        Apply the similarity ``g`` on this Euclidean set.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()
            sage: G = E.similarity_group()
            sage: g = G(1, -1, 3, 2)

        Action on a point::

            sage: point = E.point(1, -2)
            sage: point.apply_similarity(g)
            (2, -1)

        Note that this function is also registered as an action so that one can also use the multiplicative syntax::

            sage: A = get_coercion_model().get_action(G, E)
            sage: A
            Left action by Similarity group over Rational Field on Euclidean Plane over Rational Field
            sage: g * point
            (2, -1)
            sage: A(g, point)
            (2, -1)

        Action on a line::

            sage: line = E.line((1, 1), (3, -2))
            sage: line.apply_similarity(g)
            {-23 + 5*x - y = 0}
            sage: g * line
            {-23 + 5*x - y = 0}

        Action on segments::

            sage: segment = E.segment(line, start=(1, 1), end=(3, -2))
            sage: segment.apply_similarity(g)
            (5, 2) → (4, -3)
            sage: g * segment
            (5, 2) → (4, -3)

            sage: right_ray = E.segment(line, start=(3, -2))
            sage: right_ray.apply_similarity(g)
            Ray from (4, -3) in direction (-1, -5)
            sage: g * right_ray
            Ray from (4, -3) in direction (-1, -5)

            sage: left_ray = E.segment(line, end=(3, -2))
            sage: left_ray.apply_similarity(g)
            Ray to (4, -3) from direction (-1, -5)
            sage: g * left_ray
            Ray to (4, -3) from direction (-1, -5)

        Note that the action does not apply on more complicated subsets than
        points due to SageMath limitations::

            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(g, line)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(g, segment)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(g, right_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(g, left_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map

        TESTS::

            sage: point0 = E.point(0, 0)
            sage: point1 = E.point(1, 0)
            sage: point2 = E.point(0, 1)
            sage: point3 = E.point(1, 1)
            sage: point4 = E.point(-1, 0)
            sage: point5 = E.point(0, -1)
            sage: point6 = E.point(1, 0, 0)
            sage: point7 = E.point(1, 1, 0)
            sage: point8 = E.point(1, -1, 0)
            sage: line0 = E.line(point0, point1)
            sage: line1 = E.line(point1, point2)
            sage: line2 = E.line(point0, point2)
            sage: segment0 = E.segment(line0, start=point0, end=point1)
            sage: segment1 = E.segment(line0, start=point0)
            sage: segment2 = E.segment(line0, end=point0)
            sage: polygon0 = E.polygon(vertices=[point0, point1, point2])

            sage: points = [point0, point1, point2, point3, point4, point5, point6, point7, point8]
            sage: lines = [line0, line1, line2]
            sage: segments = [segment0, segment1, segment2]
            sage: polygons = [polygon0]

            sage: g0 = G(1, 1, 0, 1)
            sage: g1 = G(1, -1, 3, 1)
            sage: g2 = G(1, 0, 2, -3)
            sage: g3 = G(0, 1, 0, 0)
            sage: similarities = [g0, g1, g2, g3]

            sage: for obj in points + lines + segments + polygons:
            ....:     assert g0 * (g1 * obj) == (g0 * g1) * obj
            ....:     assert g0 * (g1 * (g2 * obj)) == (g0 * g1 * g2) * obj, obj
            ....:     assert g2 * (g3 * (g0 * obj)) == (g2 * g3 * g0) * obj, obj

            sage: from itertools import product
            sage: for p, obj, g in product(points, lines + segments + polygons, similarities):
            ....:     assert ((g * p) in (g * obj)) == (p in obj), (p, obj, g)
        """
        return self._apply_similarity(g)

    def _apply_similarity(self, g):
        r"""
        Concrete implementation of :meth:`apply_similarity` that must be
        overriden in subclasses.
        """
        raise NotImplementedError

    # TODO: move TESTS as pytest
    # TODO: carefully implement similarity reversing orientations and tests
    def apply_2x2_matrix(self, m):
        r"""
        Apply the 2x2 matrix ``m`` on this Euclidean set.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()
            sage: M = MatrixSpace(QQ, 2)
            sage: m = M([2, 1, 1, 1])

        Action on a point::

            sage: point = E.point(1, -2)
            sage: point.apply_2x2_matrix(m)
            (0, -1)

        Note that this function is also registered as an action so that one can also use the multiplicative syntax::

            sage: A = get_coercion_model().get_action(M, E)
            sage: A
            Left action by Full MatrixSpace of 2 by 2 dense matrices over Rational Field on Euclidean Plane over Rational Field
            sage: m * point
            (0, -1)
            sage: A(m, point)
            (0, -1)

        Action on an ideal point::

            sage: ideal_point = E.point(2, 1, 0)
            sage: ideal_point.apply_2x2_matrix(m)
            [5/3:1:0]
            sage: m * ideal_point
            [5/3:1:0]

        Action on a line::

            sage: line = E.line((1, 1), (3, -2))
            sage: line.apply_2x2_matrix(m)
            {-5 + x + y = 0}
            sage: m * line
            {-5 + x + y = 0}

        Action on segments::

            sage: segment = E.segment(line, start=(1, 1), end=(3, -2))
            sage: segment.apply_2x2_matrix(m)
            (3, 2) → (4, 1)
            sage: m * segment
            (3, 2) → (4, 1)

            sage: right_ray = E.segment(line, start=(3, -2))
            sage: right_ray.apply_2x2_matrix(m)
            Ray from (4, 1) in direction (1, -1)
            sage: m * right_ray
            Ray from (4, 1) in direction (1, -1)

            sage: left_ray = E.segment(line, end=(3, -2))
            sage: left_ray.apply_2x2_matrix(m)
            Ray to (4, 1) from direction (1, -1)
            sage: m * left_ray
            Ray to (4, 1) from direction (1, -1)

        Note that the action does not apply on more complicated subsets than
        points due to SageMath limitations::

            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, line)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, segment)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, right_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, left_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map

        TESTS::

            sage: point0 = E.point(0, 0)
            sage: point1 = E.point(1, 0)
            sage: point2 = E.point(0, 1)
            sage: point3 = E.point(1, 1)
            sage: point4 = E.point(-1, 0)
            sage: point5 = E.point(0, -1)
            sage: point6 = E.point(1, 0, 0)
            sage: point7 = E.point(1, 1, 0)
            sage: point8 = E.point(1, -1, 0)
            sage: line0 = E.line(point0, point1)
            sage: line1 = E.line(point1, point2)
            sage: line2 = E.line(point0, point2)
            sage: segment0 = E.segment(line0, start=point0, end=point1)
            sage: segment1 = E.segment(line0, start=point0)
            sage: segment2 = E.segment(line0, end=point0)
            sage: polygon0 = E.polygon(vertices=[point0, point1, point2])

            sage: points = [point0, point1, point2, point3, point4, point5, point6, point7, point8]
            sage: lines = [line0, line1, line2]
            sage: segments = [segment0, segment1, segment2]
            sage: polygons = [polygon0]

            sage: m0 = M([2, 1, 1, 1])
            sage: m1 = M([2, 0, 0, 2])
            sage: m2 = M([1, -1, 0, 1])
            sage: m3 = M([0, 1, -1, 0])
            sage: matrices = [m0, m1, m2, m3]

            sage: for obj in points + lines + segments + polygons:
            ....:     assert m0 * (m1 * obj) == (m0 * m1) * obj, obj
            ....:     assert m0 * (m1 * (m2 * obj)) == (m0 * m1 * m2) * obj, obj

            sage: from itertools import product
            sage: for p, obj, m in product(points, lines + segments + polygons, matrices):
            ....:     assert ((m * p) in (m * obj)) == (p in obj), (p, obj, m)

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: p = E.polygon(vertices = [(1,0),(0,1),(-1,-1)])
            sage: p
            Polygon(vertices=[(1, 0), (0, 1), (-1, -1)])
            sage: matrix(ZZ, [[0, 1], [1, 0]]) * p
            Polygon(vertices=[(0, 1), (-1, -1), (1, 0)])
            sage: matrix(ZZ,[[2, 0], [0, 1]]) * p
            Polygon(vertices=[(2, 0), (0, 1), (-2, -1)])
        """
        return self._apply_2x2_matrix(m)

    def _apply_2x2_matrix(self, m):
        raise NotImplementedError

    def apply_3x3_matrix(self, m):
        r"""
        Apply the 3x3 matrix ``m`` on this Euclidean set.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane

            sage: E = EuclideanPlane()
            sage: M = MatrixSpace(QQ, 3)
            sage: m = M([2, 1, 1, 1, 1, 1, 0, 0, 1])

        Action on a point::

            sage: point = E.point(1, -2)
            sage: point.apply_3x3_matrix(m)
            (1, 0)

        Note that this function is also registered as an action so that one can also use the multiplicative syntax::

            sage: A = get_coercion_model().get_action(M, E)
            sage: A
            Left action by Full MatrixSpace of 3 by 3 dense matrices over Rational Field on Euclidean Plane over Rational Field
            sage: m * point
            (1, 0)
            sage: A(m, point)
            (1, 0)

        Action on an ideal point::

            sage: ideal_point = E.point(1, 1, 0)
            sage: ideal_point.apply_3x3_matrix(m)
            [3/2:1:0]
            sage: m * ideal_point
            [3/2:1:0]

        Action on a line::

            sage: line = E.line((1, 1), (3, -2))
            sage: line.apply_3x3_matrix(m)
            {-7 + x + y = 0}
            sage: m * line
            {-7 + x + y = 0}

        Action on segments::

            sage: segment = E.segment(line, start=(1, 1), end=(3, -2))
            sage: segment.apply_3x3_matrix(m)
            (4, 3) → (5, 2)
            sage: m * segment
            (4, 3) → (5, 2)

            sage: right_ray = E.segment(line, start=(3, -2))
            sage: right_ray.apply_3x3_matrix(m)
            Ray from (5, 2) in direction (1, -1)
            sage: m * right_ray
            Ray from (5, 2) in direction (1, -1)

            sage: left_ray = E.segment(line, end=(3, -2))
            sage: left_ray.apply_3x3_matrix(m)
            Ray to (5, 2) from direction (1, -1)
            sage: m * left_ray
            Ray to (5, 2) from direction (1, -1)

        Note that the action does not apply on more complicated subsets than
        points due to SageMath limitations::

            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, line)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, segment)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, right_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map
            sage: A(m, left_ray)
            Traceback (most recent call last):
            ...
            TypeError: embedding must be a parent or map

        TESTS::

            sage: point0 = E.point(0, 0)
            sage: point1 = E.point(1, 0)
            sage: point2 = E.point(0, 1)
            sage: point3 = E.point(1, 1)
            sage: point4 = E.point(-1, 0)
            sage: point5 = E.point(0, -1)
            sage: point6 = E.point(1, 0, 0)
            sage: point7 = E.point(1, 1, 0)
            sage: point8 = E.point(1, -1, 0)
            sage: line0 = E.line(point0, point1)
            sage: line1 = E.line(point1, point2)
            sage: line2 = E.line(point0, point2)
            sage: segment0 = E.segment(line0, start=point0, end=point1)
            sage: segment1 = E.segment(line0, start=point0)
            sage: segment2 = E.segment(line0, end=point0)
            sage: polygon0 = E.polygon(vertices=[point0, point1, point2])

            sage: points = [point0, point1, point2, point3, point4, point5, point6, point7, point8]
            sage: lines = [line0, line1, line2]
            sage: segments = [segment0, segment1, segment2]
            sage: polygons = [polygon0]

            sage: m0 = M([2, 1, 1, 1, 1, -2, 0, 0, 1])
            sage: m1 = M([2, 0, 0, 0, 2, 0, 0, 0, 1])
            sage: m2 = M([1, -1, 1, 0, 1, 0, 0, 0, 1])
            sage: m3 = M([0, 1, 0, -1, 0, 1, 0, 0, 1])
            sage: matrices = [m0, m1, m2, m3]

            sage: for obj in points + lines:
            ....:     assert m0 * (m1 * obj) == (m0 * m1) * obj, obj
            ....:     assert m0 * (m1 * (m2 * obj)) == (m0 * m1 * m2) * obj, obj

            sage: from itertools import product
            sage: for p, obj, m in product(points, lines, matrices):
            ....:     assert ((m * p) in (m * obj)) == (p in obj), (p, obj, m)
        """
        if self.is_oriented():
            if m[2, 0] or m[2, 1] or m[2, 2] < 0:
                raise ValueError("ambiguous action of 3x3 matrix on oriented Euclidean set")
        return self._apply_3x3_matrix(m)

    def _apply_3x3_matrix(self, m):
        raise NotImplementedError

    # TODO: Add is_subset() and a test method.

    # TODO: Add _an_element_ and some_elements()

    # TODO: Add is_empty() and __bool__

    # TODO: Add is_point()

    def is_oriented(self):
        r"""
        Return whether this is a set with an explicit orientation.

        Some sets come in two flavors. There are oriented segments and
        unoriented segments.

        This method answers whether a set is in the oriented kind if there is a
        choice.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

        Segments are normally oriented::

            sage: s = E((0, 0)).segment((1, 0))
            sage: s.is_oriented()
            True

        We can explicitly ask for an unoriented segment::

            sage: u = s.unoriented()
            sage: u.is_oriented()
            False

        Points are not oriented, there is no choice of orientation::

            sage: p = E((0, 0))
            sage: p.is_oriented()
            False

        """
        return isinstance(self, EuclideanOrientedSet)

    # TODO: Add __hash__ and test method

    # TODO: Add random_set for testing.


class EuclideanOrientedSet(EuclideanSet):
    r"""
    Base class for sets that have an explicit orientation.

    .. SEEALSO::

        :meth:`EuclideanSet.is_oriented`

    """
    def __contains__(self, point):
        return point in self.unoriented()

    def unoriented(self):
        return self.change(oriented=False)


class EuclideanFacade(EuclideanSet, Parent):
    r"""
    A subset of the Euclidean plane that is itself a parent.

    This is the base class for all Euclidean sets that are not points.
    This class solves the problem that we want sets to be "elements" of the
    Euclidean plane but at the same time, we want these sets to live as parents
    in the category framework of SageMath; so they have a Parent with Euclidean
    points as their Element class.

    SageMath provides the (not very frequently used and somewhat flaky) facade
    mechanism for such parents. Such sets being a facade, their points can be
    both their elements and the elements of the Euclidean plane.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: c = E.circle((0, 0), radius=1)
        sage: p = c.center()
        sage: p in c
        False
        sage: p.parent() is E
        True
        sage: q = c.an_element()
        sage: q
        (1, 0)
        sage: q in c
        True
        sage: q.parent() is E
        True

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanFacade
        sage: isinstance(c, EuclideanFacade)
        True

    """

    def __init__(self, parent, category=None):
        Parent.__init__(self, facade=parent, category=category)

    def parent(self):
        r"""
        Return the Euclidean plane this is a subset of.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: c = E.circle((0, 0), radius=1)
            sage: c.parent()
            Euclidean Plane over Rational Field

        """
        return self.facade_for()[0]

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as a point of this set.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: c = E.circle((0, 0), radius=1)
            sage: c((1, 0))
            (1, 0)
            sage: c((0, 0))
            Traceback (most recent call last):
            ...
            ValueError: point not contained in this set

        """
        x = self.parent()(x)

        if isinstance(x, EuclideanPoint):
            if not self.__contains__(x):
                raise ValueError("point not contained in this set")

        return x

    def base_ring(self):
        r"""
        Return the ring over which points of this set are defined.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: c = E.circle((0, 0), radius=1)
            sage: c.base_ring()
            Rational Field

        """
        return self.parent().base_ring()


class EuclideanCircle(EuclideanFacade):
    r"""
    A circle in the Euclidean plane.

    INPUT:

    - ``parent`` -- the :class:`EuclideanPlane` containing this circle

    - ``center`` -- the :class:`EuclideanPoint`` at the center of this circle

    - ``radius_squared`` -- the square of the radius of this circle

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: c = EuclideanPlane().circle((0, 0), radius=1)

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanCircle
        sage: isinstance(c, EuclideanCircle)
        True
        sage: TestSuite(c).run()

    .. SEEALSO::

        :meth:`EuclideanPlane.circle` for a method to create circles

    """

    def __init__(self, parent, center, radius_squared):
        super().__init__(parent)

        self._center = center
        self._radius_squared = radius_squared

    def __eq__(self, other):
        if not isinstance(other, EuclideanCircle):
            return False

        if self.parent() is not other.parent():
            return False

        return self._center == other._center and self._radius_squared == other._radius_squared

    def __hash__(self):
        return hash((self._center, self._radius_squared))

    def _an_element_(self):
        # TODO: Maybe try to be smarter to construct an element without taking a field extension here.
        return self.parent()(self._center.vector() + self.parent().vector_space()((self._radius_squared.sqrt(), 0)))

    def _repr_(self):
        r"""
        Return a printable representation of this circle.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: c = EuclideanPlane().circle((0, 0), radius=1)
            sage: c
            { x² + y² = 1 }

        """
        x, y = self._center

        x = f"(x - {x})" if x else "x"
        y = f"(y - {y})" if y else "y"

        return f"{{ {x}² + {y}² = {self._radius_squared} }}"

    def change(self, *, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this circle.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~EuclideanPlane.base_ring`); the ring over which the new
          circle will be defined.

        - ``geometry`` -- a :class:`EuclideanGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new circle.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since circles cannot
          have an explicit orientation. See :meth:`~EuclideanSet.is_oriented`.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: c = E.circle((0, 0), radius=1)

        We change the base ring over which this circle is defined::

            sage: c.change(ring=AA)
            { x² + y² = 1 }

        We cannot change the orientation of a circle::

            sage: c.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: circles cannot have an explicit orientation

            sage: c.change(oriented=False)
            { x² + y² = 1 }

        """
        if ring is not None or geometry is not None:
            P = self.parent()
            Q = P.change_ring(ring, geometry=geometry)
            if P is not Q:
                self = Q.circle(self._center, radius_squared=self._radius_squared, check=False)

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("circles cannot have an explicit orientation")

        return self

    def _normalize(self):
        r"""
        Return this set possibly rewritten in a simpler form.

        This implements :meth:`EuclideanSet._normalize`.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: circle = E.circle((0, 0), radius=0, check=False)
            sage: circle
            { x² + y² = 0 }
            sage: circle._normalize()
            (0, 0)

        """
        if self.parent().geometry._zero(self._radius_squared):
            return self._center

        return self

    def center(self):
        r"""
        Return the point at the center of the circle.
        """
        return self._center

    def radius_squared(self):
        r"""
        Return the square of the radius of the circle.
        """
        return self._radius_squared

    def point_position(self, point):
        r"""
        Return 1 if point lies in the circle, 0 if the point lies on the circle,
        and -1 if the point lies outide the circle.
        """
        # TODO: Deprecate?
        value = (
            (point[0] - self._center[0]) ** 2
            + (point[1] - self._center[1]) ** 2
            - self._radius_squared
        )

        if value > 0:
            return -1
        if value < 0:
            return 1
        return 0

    def closest_point_on_line(self, point, direction_vector):
        r"""
        Consider the line through the provided point in the given direction.
        Return the closest point on this line to the center of the circle.
        """
        # TODO: Rewrite or deprecate and generalize
        V3 = self.parent().base_ring() ** 3
        V2 = self.parent().base_ring() ** 2

        cc = V3((self._center[0], self._center[1], 1))
        # point at infinite orthogonal to direction_vector:
        dd = V3((direction_vector[1], -direction_vector[0], 0))
        l1 = cc.cross_product(dd)

        pp = V3((point[0], point[1], 1))
        # direction_vector pushed to infinity
        ee = V3((direction_vector[0], direction_vector[1], 0))
        l2 = pp.cross_product(ee)

        # This is the point we want to return
        rr = l1.cross_product(l2)
        try:
            return V2((rr[0] / rr[2], rr[1] / rr[2]))
        except ZeroDivisionError:
            raise ValueError(
                "Division by zero error. Perhaps direction is zero. "
                + "point="
                + str(point)
                + " direction="
                + str(direction_vector)
                + " circle="
                + str(self)
            )

    # TODO: Not used anywhere.
    ## def line_position(self, point, direction_vector):
    ##     r"""
    ##     Consider the line through the provided point in the given direction.
    ##     We return 1 if the line passes through the circle, 0 if it is tangent
    ##     to the circle and -1 if the line does not intersect the circle.
    ##     """
    ##     return self.point_position(self.closest_point_on_line(point, direction_vector))

    # TODO: Create Segment class.
    def line_segment_position(self, p, q):
        r"""
        Consider the open line segment pq.We return 1 if the line segment
        enters the interior of the circle, zero if it touches the circle
        tangentially (at a point in the interior of the segment) and
        and -1 if it does not touch the circle or its interior.
        """
        if self.point_position(p) == 1:
            return 1
        if self.point_position(q) == 1:
            return 1
        r = self.closest_point_on_line(p, q - p)
        pos = self.point_position(r)
        if pos == -1:
            return -1
        # This checks if r lies in the interior of pq
        if p[0] == q[0]:
            if (p[1] < r[1] and r[1] < q[1]) or (p[1] > r[1] and r[1] > q[1]):
                return pos
        elif (p[0] < r[0] and r[0] < q[0]) or (p[0] > r[0] and r[0] > q[0]):
            return pos
        # It does not lie in the interior.
        return -1

    # TODO: Not used anywhere.
    ## def tangent_vector(self, point):
    ##     r"""
    ##     Return a vector based at the provided point (which must lie on the circle)
    ##     which is tangent to the circle and points in the counter-clockwise
    ##     direction.

    ##     EXAMPLES::

    ##         sage: from flatsurf.geometry.circle import Circle
    ##         sage: c=Circle(vector((0,0)), 2, base_ring=QQ)
    ##         sage: c.tangent_vector(vector((1,1)))
    ##         (-1, 1)
    ##     """
    ##     if not self.point_position(point) == 0:
    ##         raise ValueError("point not on circle.")
    ##     return vector((self._center[1] - point[1], point[0] - self._center[0]))

    # TODO: Not used anywhere.
    ## def other_intersection(self, p, v):
    ##     r"""
    ##     Consider a point p on the circle and a vector v. Let L be the line
    ##     through p in direction v. Then L intersects the circle at another
    ##     point q. This method returns q.

    ##     Note that if p and v are both in the field of the circle,
    ##     then so is q.

    ##     EXAMPLES::

    ##         sage: from flatsurf.geometry.circle import Circle
    ##         sage: c=Circle(vector((0,0)), 25, base_ring=QQ)
    ##         sage: c.other_intersection(vector((3,4)),vector((1,2)))
    ##         (-7/5, -24/5)
    ##     """
    ##     pp = self._V3((p[0], p[1], self._base_ring.one()))
    ##     vv = self._V3((v[0], v[1], self._base_ring.zero()))
    ##     L = pp.cross_product(vv)
    ##     cc = self._V3((self._center[0], self._center[1], self._base_ring.one()))
    ##     vvperp = self._V3((-v[1], v[0], self._base_ring.zero()))
    ##     # line perpendicular to L through center:
    ##     Lperp = cc.cross_product(vvperp)
    ##     # intersection of L and Lperp:
    ##     rr = L.cross_product(Lperp)
    ##     r = self._V2((rr[0] / rr[2], rr[1] / rr[2]))
    ##     return self._V2((2 * r[0] - p[0], 2 * r[1] - p[1]))

    def __rmul__(self, similarity):
        r"""
        Apply a similarity to the circle.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.square_torus()
            sage: c = s.polygon(0).circumscribed_circle()
            sage: c
            Circle((1/2, 1/2), 1/2)
            sage: s.edge_transformation(0,2)
            (x, y) |-> (x, y - 1)
            sage: s.edge_transformation(0,2) * c
            Circle((1/2, -1/2), 1/2)
        """
        # TODO: Implement this properly for this parent
        from .similarity import SimilarityGroup

        SG = SimilarityGroup(self.parent().base_ring())
        s = SG(similarity)
        return self.parent().circle(
            s(self._center), radius_squared=s.det() * self._radius_squared
        )

    ## def __str__(self):
    ##     return (
    ##         "circle with center "
    ##         + str(self._center)
    ##         + " and radius squared "
    ##         + str(self._radius_squared)
    ##     )

    def __contains__(self, point):
        if point.is_ideal():
            return False
        v = self._center.vector() - point.vector()
        return v.dot_product(v) == self._radius_squared


class EuclideanPoint(EuclideanSet, Element):
    r"""
    A point in the :class:`EuclideanPlane`.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()

        sage: p = E.point(0, 0)

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanPoint
        sage: isinstance(p, EuclideanPoint)
        True

        sage: TestSuite(p).run()

    .. SEEALSO::

        :meth:`EuclideanPlane.point` for ways to create points

    """

    def __init__(self, parent, x, y, z):
        super().__init__(parent)

        self._x = x
        self._y = y
        self._z = z

        # TODO: handle non-normalized coordinates
        assert self._z.is_one() or (self._z.is_zero() and self._y.is_one() or (self._y.is_zero() and self._x.is_one()))

    def is_ideal(self):
        return not self._z

    def __iter__(self):
        r"""
        Return an iterator over the coordinates of this point.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: p = E.point(1, 2)
            sage: list(p)
            [1, 2]

        """
        if not self._z:
            raise ValueError("ideal point")
        yield self._x
        yield self._y

    def vector(self):
        if not self._z:
            raise ValueError("ideal point")
        v = self.parent().vector_space()((self._x, self._y))
        v.set_immutable()
        return v

    def translate(self, v):
        if self.is_ideal():
            return self
        v = self.parent().vector_space()(v)
        return self.parent().point(*(self.vector() + v))

    def _apply_scalar(self, r):
        if not self._z:
            return self
        return self.parent().point(r * self._x, r * self._y, self._z)

    def _apply_similarity(self, g):
        return self._apply_3x3_matrix(g.matrix())

    def _apply_2x2_matrix(self, m):
        V = self.parent().base_ring() ** 2
        x, y = m * V((self._x, self._y))
        return self.parent().point(x, y, self._z)

    def _apply_3x3_matrix(self, m):
        V = self.parent().base_ring() ** 3
        return self.parent().point(*(m * V((self._x, self._y, self._z))))

    def _richcmp_(self, other, op):
        r"""
        Return how this point compares to ``other`` with respect to the ``op``
        operator.

        This is only implemented for the operators ``==`` and ``!=``. It
        returns whether two points are the same.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: E((0, 0)) == E((0, 0))
            True

        .. SEEALSO::

            :meth:`EuclideanSet.__contains__` to check for containment of a
            point in a set

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        if op == op_EQ:
            if not isinstance(other, EuclideanPoint):
                return False

            return self.parent().geometry._equal_point(self, other)

        return super()._richcmp_(other, op)

    def _repr_(self):
        r"""
        Return a printable representation of this point.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: p = E.point(0, 0)
            sage: p
            (0, 0)
            sage: E.point(1, 0, 0)
            [1:0:0]

        """
        if not self._z:
            return f"[{self._x}:{self._y}:0]"
        return repr(tuple(self))

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coordinate of this point.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: p = E.point(1, 2)
            sage: p[0]
            1
            sage: p[1]
            2
            sage: p[2]
            Traceback (most recent call last):
            ...
            KeyError

        """
        if i == 0:
            return self._x
        if i == 1:
            return self._y

        raise KeyError

    def change(self, *, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this point.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~EuclideanPlane.base_ring`); the ring over which the new
          point will be defined.

        - ``geometry`` -- a :class:`EuclideanGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new point.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); must be ``None`` or ``False`` since points cannot
          have an explicit orientation. See :meth:`~EuclideanSet.is_oriented`.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()

            sage: p = E((0, 0))

        We change the base ring over which this point is defined::

            sage: p.change_ring(ring=AA)
            (0, 0)

        We cannot change the orientation of a point:

            sage: p.change(oriented=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: points cannot have an explicit orientation

            sage: p.change(oriented=False)
            (0, 0)

        Ideal points::

            sage: E.point(1, 0, 0).change_ring(AA)
            [1:0:0]

        TESTS:

        Test that trivial changes do not trigger a copy::

            sage: point = E.point(1, 0)
            sage: point.change_ring(QQ) is point
            True
            sage: point.change(geometry=point.parent().geometry) is point
            True

        """
        if ring is not None or geometry is not None:
            P = self.parent()
            Q = P.change_ring(ring, geometry=geometry)
            if P is not Q:
                self = Q.point(self._x, self._y, self._z)

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            raise NotImplementedError("points cannot have an explicit orientation")

        return self

    def segment(self, end):
        end = self.parent()(end)

        line = self.parent().line(self, end)

        return self.parent().segment(line, start=self, end=end)

    def __hash__(self):
        r"""
        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: len(set(map(hash, [E.point(0, 0), E.point(1, 0), E.point(0, 1), E.point(1, 0, 0), E.point(0, 1, 0)]))) == 5
            True
        """
        return hash((self._x, self._y, self._z))


class EuclideanLine(EuclideanFacade):
    r"""
    A line in the Euclidean plane.

    This is a common base class for oriented and unoriented lines, see
    :class:`EuclideanOrientedLine` and :class:`EuclideanUnorientedLine`.

    Internally, we represent a line by its equation, i.e., the ``a``,
    ``b``, ``c`` such that points on the line satisfy

    .. MATH::

        a + bx + cy = 0

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: line = E.line((0, 0), (1, 1))

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanLine
        sage: isinstance(line, EuclideanLine)
        True
        sage: TestSuite(line).run()

    .. SEEALSO::

        :meth:`EuclideanPlane.line`

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

    def is_compact(self):
        return False

    def is_ideal(self):
        return not self._b and not self._c

    def dimension(self):
        return 1

    def change(self, *, ring=None, geometry=None, oriented=None):
        r"""
        Return a modified copy of this line.

        INPUT:

        - ``ring`` -- a ring (default: ``None`` to keep the current
          :meth:`~EuclideanPlane.base_ring`); the ring over which the new line
          will be defined.

        - ``geometry`` -- a :class:`EuclideanGeometry` (default: ``None`` to
          keep the current geometry); the geometry that will be used for the
          new line.

        - ``oriented`` -- a boolean (default: ``None`` to keep the current
          orientedness); whether the new line should be oriented.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(AA)

        The base ring over which this line is defined can be changed::

            sage: E.line((0, 0), (1, 1)).change_ring(QQ)
            {-x + y = 0}

        But we cannot change the base ring if the line's equation cannot be
        expressed in the smaller ring::

            sage: E.line((0, 0), (1, AA(2).sqrt())).change_ring(QQ)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real ... to Rational

        We can forget the orientation of a line::

            sage: line = E.line((0, 0), (1, 1))
            sage: line.is_oriented()
            True
            sage: line = line.change(oriented=False)
            sage: line.is_oriented()
            False

        We can (somewhat randomly) pick the orientation of a line::

            sage: line = line.change(oriented=True)
            sage: line.is_oriented()
            True


        TESTS:

        Test that trivial changes do not trigger a copy::

            sage: line = E.line((0, 0), (1, 1))
            sage: line.change_ring(AA) is line
            True
            sage: line.change(geometry=line.parent().geometry) is line
            True
        """
        if ring is not None or geometry is not None:
            P = self.parent()
            Q = P.change_ring(ring, geometry=geometry)
            if P is not Q:
                self = Q.line(
                    self._a,
                    self._b,
                    self._c,
                    check=False,
                    oriented=self.is_oriented(),
                )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            self = self.parent().line(
                self._a, self._b, self._c, check=False, oriented=oriented
            )

        return self

    def _repr_(self):
        r"""
        Return a printable representation of this line.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(AA)
            sage: E.line((0, 0), (1, 1))
            {-x + y = 0}
            sage: E.line((0, 1), (1, 1))
            {-1 + y = 0}
            sage: E.line((1, 0), (1, 1))  # TODO: Fix printing
            {1 + -x = 0}
            sage: E.line(1, 0, 0)  # TODO: should we keep this for the line at infinity?
            {1 + 0 = 0}

        """
        a, b, c = self.equation(normalization=["gcd", None])

        from sage.all import PolynomialRing

        R = PolynomialRing(self.parent().base_ring(), names=["x", "y"])
        polynomial_part = R({(1, 0): b, (0, 1): c})
        if self.parent().geometry._sgn(a) != 0:
            return f"{{{repr(a)} + {repr(polynomial_part)} = 0}}"
        else:
            return f"{{{repr(polynomial_part)} = 0}}"

    def equation(self, normalization=None):
        r"""
        Return an equation for this line as a triple ``a``, ``b``, ``c`` such
        that the line is given by the points satisfying

        .. MATH::

            a + bx + cy = 0

        INPUT:

        - ``normalization`` -- how to normalize the coefficients; the default
          ``None`` is not to normalize at all. Other options are ``gcd``, to
          divide the coefficients by their greatest common divisor, ``one``, to
          normalize the first non-zero coefficient to ±1. This can also be a
          list of such values which are then tried in order and exceptions are
          silently ignored unless they happen at the last option.

        If this line :meth:`is_oriented`, then the sign of the coefficients
        is chosen to encode the orientation of this line. The sign is such
        that the half plane obtained by replacing ``=`` with ``≥`` in the
        equation is on the left of the line.

        Note that the output might not uniquely describe the line. The
        coefficients are only unique up to scaling.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(AA)
            sage: E.line((0, 0), (1, 1)).equation()
            (0, -1, 1)
            sage: E.line((0, 1), (1, 1)).equation()
            (-1, 0, 1)
            sage: E.line((1, 0), (1, 1)).equation()
            (1, -1, 0)

        Some normalizations might not be possible over some base rings::

            sage: E = EuclideanPlane(ZZ)
            sage: line = E.line((1, 3), (6, 8))
            sage: line
            {-2 + -x + y = 0}
            sage: line.equation()
            (-10, -5, 5)
            sage: line.equation(normalization="one")
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: line.equation(normalization="gcd")
            (-2, -1, 1)

        In such cases, we can also use a list of normalizations to select the
        best one possible::

            sage: line.equation(normalization=["one", "gcd", None])
            (-2, -1, 1)

        .. SEEALSO::

            :meth:`HyperbolicGeodesic.equation` which does essentially the same
            for geodesics in the hyperbolic plane

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

            # TODO: Do not call into the hyperbolic machinery (but maybe the other way round.)
            from flatsurf.geometry.hyperbolic import HyperbolicGeodesic

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

    def _an_element_(self):
        if self._b:
            return self.parent().point(-self._a / self._b, 0)

        assert self._c
        return self.parent().point(-self._a / self._c, 0)

    def _apply_scalar(self, r):
        r"""
        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane(QQ)
            sage: E.line(1, 1, 1)
            {1 + x + y = 0}
            sage: 2 * E.line(1, 1, 1)
            {2 + x + y = 0}
        """
        if not r:
            return self.parent().point(0, 0)
        return self.parent().line(self._a, self._b / r, self._c / r, oriented=self.is_oriented(), check=False)

    def _apply_similarity(self, g):
        A = self._a
        B = self._b
        C = self._c
        g = ~g
        a = g._a
        b = g._b
        s = g._s
        t = g._t
        sign = g._sign
        if sign.is_one():
            AA = A + B*s + C*t
            BB = B*a + C*b
            CC = -B*b + C*a
        else:
            raise NotImplementedError
        return self.parent().line(AA, BB, CC, oriented=self.is_oriented(), check=False)

    def _apply_2x2_matrix(self, m):
        A = self._a
        B = self._b
        C = self._c
        m = m.inverse_of_unit()
        a, b = m[0]
        c, d = m[1]
        # TODO: check that it does actually make sense for matrices reversing orientation
        AA = A
        BB = B * m[0, 0] + C * m[1, 0]
        CC = B * m[0, 1] + C * m[1, 1]
        return self.parent().line(AA, BB, CC, oriented=self.is_oriented(), check=False)

    def _apply_3x3_matrix(self, m):
        V = self.parent().base_ring() ** 3
        BB, CC, AA = V((self._b, self._c, self._a)) * m.inverse_of_unit()
        return self.parent().line(AA, BB, CC, oriented=self.is_oriented(), check=False)

    def projection(self, point):
        # Move the line to the origin, i.e., instead of a + bx + cy = 0,
        # consider bx + cy = 0.
        # Let v be a vector parallel to this line.
        shift = self.an_element()

        point = point.translate(-shift.vector())

        (x, y) = point
        v = (self._c, -self._b)
        vv = ~(v[0] ** 2 + v[1] ** 2)

        p = (v[0] ** 2 * x + v[0] * v[1] * y, v[0] * v[1] * x + v[1] * v[1] * y)
        p = (p[0] * vv, p[1] * vv)

        assert p[0] * self._b + p[1] * self._c == 0

        return self.parent().point(*p).translate(shift.vector())

    def __contains__(self, point):
        return (self._a * point._z + self._b * point._x + self._c * point._y).is_zero()

    def __eq__(self, other):
        # TODO: This is literally identical to __eq__ of hyperbolic geodesics. One should call the other instead.
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

    def __hash__(self):
        if not self.parent().base_ring().is_exact():
            raise TypeError("cannot hash geodesic defined over inexact base ring")

        return hash(
            (type(self), self.equation(normalization=["one", "gcd"]))
        )
        
    def translate(self, v):
        return self.parent().line(self._a - self._b * v[0] - self._c * v[1], self._b, self._c, oriented=self.is_oriented(), check=False)


class EuclideanOrientedLine(EuclideanLine, EuclideanOrientedSet):
    r"""
    A line in the Euclidean plane with an explicit orientation.

    Internally, we represent a line by its equation, i.e., the ``a``,
    ``b``, ``c`` such that points on the line satisfy

    .. MATH::

        a + bx + cy = 0

    The orientation of that line is such that

    .. MATH::

        a + bx + cy \ge 0

    is to its left.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: line = E.line((0, 0), (1, 1))

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanOrientedLine
        sage: isinstance(line, EuclideanOrientedLine)
        True
        sage: TestSuite(line).run()

    """

    def direction(self, normalization=None):
        r"""
        Return a vector pointing in the direction of this oriented line.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: line = E.line((0, 0), (1, 1))
            sage: line.direction()
            (1, 1)
            sage: line = E.line((0, 0), (2, 2))
            sage: line.direction()
            (2, 2)
            sage: line.direction(normalization="gcd")
            (1, 1)

        """
        a, b, c = self.equation(normalization=normalization)
        return self.parent().vector_space()((c, -b))

    def ccw(self, other):
        # TODO: Check for intersection and make sure other is an oriented line. This is identical to geodesic ccw.
        sgn = self.parent().geometry._sgn
        return sgn(ccw((self._c, -self._b), (other._c, -other._b)))

    def __neg__(self):
        r"""
        Return this line with reversed orientation.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: line = E.line((0, 0), (1, 1))
            sage: -line
            {x - y = 0}

        """
        return self.parent().line(-self._a, -self._b, -self._c, check=False)

    def start(self):
        # TODO: Return an ideal point.
        return None

    def end(self):
        # TODO: Return an ideal point.
        return None


class EuclideanUnorientedLine(EuclideanLine):
    pass


class EuclideanSegment(EuclideanFacade):
    r"""
    A line segment in the Euclidean plane.

    This is a common base class for oriented and unoriented segments, see
    :class:`EuclideanOrientedSegment` and :class:`EuclideanUnorientedSegment`.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: start = E((0, 0))
        sage: end = E((1, 1))
        sage: segment = start.segment(end)

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanSegment
        sage: isinstance(segment, EuclideanSegment)
        True
        sage: TestSuite(segment).run()

    .. SEEALSO::

        :meth:`EuclideanPlane.segment`
        :meth:`EuclideanPoint.segment`

    """

    def __init__(self, parent, line, start, end):
        super().__init__(parent)

        if not isinstance(line, EuclideanLine):
            raise TypeError("line must be a Euclidean line")
        if start is not None and not (isinstance(start, EuclideanPoint) and start in line):
            raise TypeError("start must be a Euclidean point on line")
        if end is not None and not (isinstance(end, EuclideanPoint) and end in line):
            raise TypeError("end must be a Euclidean point on line")

        self._line = line
        self._start = start
        self._end = end

    def _normalize(self):
        r"""
        Return a normalized version of this segment.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: line = E.line((0, 0), (1, 1))

        A segment with two identical endpoints, is a point::

            sage: segment = E.segment(line, start=(0, 0), end=(0, 0), assume_normalized=True, check=False)
            sage: segment
            (0, 0) → (0, 0)
            sage: segment._normalize()
            (0, 0)

        A segment without endpoints is a line::

            sage: segment = E.segment(line, start=None, end=None, assume_normalized=True, check=False)
            sage: segment  # TODO: Should we fix this printing?
            Ray to None from direction (1, 1)
            sage: segment._normalize()
            {-x + y = 0}

        """
        line = self._line
        start = self._start
        end = self._end

        if start is None and end is None:
            return self._line.change(oriented=self.is_oriented())

        if start == end:
            return start

        if start is None:
            if not self.is_oriented():
                start, end = end, start
                line = -line

        return self

    def _apply_scalar(self, r):
        if not r:
            return self.parent().point(0, 0)
        line = self._line._apply_scalar(r)
        start = None if self._start is None else self._start._apply_scalar(r)
        end = None if self._end is None else self._end._apply_scalar(r)
        if r < 0:
            start, end = end, start
        return self.parent().segment(line, start, end, oriented=self.is_oriented(), check=False)

    def _apply_similarity(self, g):
        line = self._line._apply_similarity(g)
        start = None if self._start is None else self._start._apply_similarity(g)
        end = None if self._end is None else self._end._apply_similarity(g)
        if not g.sign().is_one():
            start, end = end, start
        return self.parent().segment(line, start, end, oriented=self.is_oriented(), check=False)

    def _apply_2x2_matrix(self, m):
        line = self._line._apply_2x2_matrix(m)
        start = None if self._start is None else self._start._apply_2x2_matrix(m)
        end = None if self._end is None else self._end._apply_2x2_matrix(m)
        if m.det() < 0:
            start, end = end, start
        return self.parent().segment(line, start, end, oriented=self.is_oriented(), check=False)

    def _apply_3x3_matrix(self, m):
        line = self._line._apply_3x3_matrix(m)
        start = None if self._start is None else self._start._apply_3x3_matrix(m)
        end = None if self._end is None else self._end._apply_3x3_matrix(m)
        if m.det() < 0:
            start, end = end, start
        return self.parent().segment(line, start, end, oriented=self.is_oriented(), check=False)

    def distance(self, point):
        point = self.parent()(point)

        # To compute the distance from the point to the segment, we compute the
        # distance from the point to the line containing the segment.
        # If the closest point on the line is on the segment, that's the
        # distance to the segment. If not, the minimum distance is at one of
        # the endpoints of the segment.
        norm = self.parent().norm()
        p = self._line.projection(point)
        if p in self:
            return norm.from_vector(point.vector() - p.vector())
        return min(
            norm.from_vector(point.vector() - self._start.vector()),
            norm.from_vector(point.vector() - self._end.vector()),
        )

    def change(self, ring=None, geometry=None, oriented=None):
        if ring is not None or geometry is not None:
            P = self.parent()
            Q = P.change_ring(ring, geometry=geometry)
            if P is not Q:
                self = Q.segment(
                        self._line,
                        self._start,
                        self._end,
                        check=False,
                        oriented=self.is_oriented(),
                    )

        if oriented is None:
            oriented = self.is_oriented()

        if oriented != self.is_oriented():
            self = self.parent().segment(
                self._line, self._start, self._end, check=False, oriented=oriented
            )

        return self

    def translate(self, v):
        return self.parent().segment(self._line.translate(v), None if self._start is None else self._start.translate(v), None if self._end is None else self._end.translate(v), oriented=self.is_oriented(), check=False)

    def __neg__(self):
        return self.parent().segment(-self._line, self._end, self._start, oriented=self.is_oriented(), check=False)


class EuclideanOrientedSegment(EuclideanSegment, EuclideanOrientedSet):
    r"""
    An oriented line segment in the Euclidean plane going from a ``start``
    point to an ``end`` point.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: start = E((0, 0))
        sage: end = E((1, 0))
        sage: segment = start.segment(end)

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanOrientedSegment
        sage: isinstance(segment, EuclideanOrientedSegment)
        True
        sage: TestSuite(segment).run()

    """

    def _repr_(self):
        r"""
        Return a printable representation of this segment.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: start = E((0, 0))
            sage: end = E((1, 0))
            sage: segment = start.segment(end)
            sage: segment
            (0, 0) → (1, 0)

        """
        # TODO: we reproduce the call in EuclideanLine._repr_, maybe this could be abstracted
        d = self._line.direction(normalization=["gcd", None])
        if self._start is None:
            return f"Ray to {self._end!r} from direction {d!r}"
        if self._end is None:
            return f"Ray from {self._start!r} in direction {d!r}"

        return f"{self._start!r} → {self._end!r}"

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        return (
            self._line == other._line and self._start == other._start and self._end == other._end
        )

    def __hash__(self):
        return hash((self._line, self._start, self._end))

    def start(self):
        # TODO: The EuclideanPlane should Model RP², so this should not be None but an ideal point.
        return self._start

    def end(self):
        # TODO: The EuclideanPlane should Model RP², so this should not be None but an ideal point.
        return self._end

    def __neg__(self):
        return self.parent().segment(-self._line, self._end, self._start, check=False)

    def direction(self):
        return self._line.direction()

    def is_compact(self):
        return self._start is not None and self._end is not None

    def dimension(self):
        return 1


class EuclideanUnorientedSegment(EuclideanSegment):
    r"""
    An unoriented line segment in the Euclidean plane connecting ``start`` and
    ``end``.

    EXAMPLES::

        sage: from flatsurf import EuclideanPlane
        sage: E = EuclideanPlane()
        sage: start = E((0, 0))
        sage: end = E((1, 0))
        sage: segment = start.segment(end)
        sage: segment = segment.unoriented()

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanUnorientedSegment
        sage: isinstance(segment, EuclideanUnorientedSegment)
        True
        sage: TestSuite(segment).run()

    """

    def _repr_(self):
        r"""
        Return a printable representation of this segment.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: start = E((0, 0))
            sage: end = E((1, 0))
            sage: segment = start.segment(end).unoriented()
            sage: segment
            (0, 0) — (1, 0)

        """
        return f"{self._start!r} — {self._end!r}"

    def __contains__(self, point):
        if point not in self._line:
            return False

        if point.is_ideal():
            return self._start is None or self._end is None

        if self._start is not None:
            t, _ = time_on_ray(self._start.vector(), self._line.direction(), point)
            if t < 0:
                return False
        if self._end is not None:
            t, _ = time_on_ray(self._end.vector(), -self._line.direction(), point)
            if t < 0:
                return False

        return True

    def __eq__(self, other):
        if type(self) is not type(other):
            return False
        return (
            self._line.unoriented() == other._line.unoriented() and {self._start, self._end} == {other._start, other._end}
        )

    def __hash__(self):
        return hash((self._line.unoriented(), {self._start, self._end}))


# TODO: should we allow the "complement" of a polygon to be a proper polygon?
# (ie the complement of a unit square)
class EuclideanPolygon(EuclideanFacade):
    r"""
    A (possibly non-convex) simple polygon in the plane `\mathbb{R}^2`.

    EXAMPLES::

        sage: from flatsurf import polygons, Polygon
        sage: s = polygons.square()
        sage: s
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanPolygon
        sage: isinstance(s, EuclideanPolygon)
        True

        sage: TestSuite(s).run()

    """

    def __init__(self, parent, edges, category=None):
        from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons
        if category is None:
            category = EuclideanPolygons(parent.base_ring())

        category &= EuclideanPolygons(parent.base_ring())

        super().__init__(parent, category=category)

        self._edges = edges

        if "Convex" not in category.axioms() and self.is_convex():
            self._refine_category_(category.Convex())

        # The category is not refined automatically to the WithAngles()
        # subcategory since computation of angles can be very costly.
        # The category gets further refined when angles() is invoked.

    def is_compact(self):
        return all(s.is_compact() for s in self.sides())

    def _an_element_(self):
        r"""
        Return a point of this polygon.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()
            sage: s.an_element()
            (0, 0)

        """
        return self.parent()(self.vertices()[0])

    @cached_method
    def __hash__(self):
        return hash(self._edges)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf import polygons, Polygon
            sage: p1 = polygons.square()
            sage: p2 = Polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)], base_ring=AA)
            sage: p1 == p2
            True

            sage: p3 = Polygon(edges=[(2, 0), (-1, 1), (-1, -1)])
            sage: p1 == p3
            False

        TESTS::

            sage: from flatsurf import Polygon, polygons
            sage: p1 = polygons.square()
            sage: p2 = Polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)], base_ring=AA)
            sage: p1 != p2
            False

            sage: p3 = Polygon(edges=[(2, 0), (-1, 1), (-1, -1)])
            sage: p1 != p3
            True

        """
        if not isinstance(other, EuclideanPolygon):
            return False

        return self._edges == other._edges

    # TODO: Do we need this?
    def cmp(self, other):
        r"""
        Implement a total order on polygons
        """
        if not isinstance(other, EuclideanPolygon):
            raise TypeError("__cmp__ only implemented for ConvexPolygons")
        if not self.base_ring() == other.base_ring():
            raise ValueError(
                "__cmp__ only implemented for ConvexPolygons defined over the same base_ring"
            )
        sign = len(self.vertices()) - len(other.vertices())
        if sign > 0:
            return 1
        if sign < 0:
            return -1
        sign = self.area() - other.area()
        if sign > self.base_ring().zero():
            return 1
        if sign < self.base_ring().zero():
            return -1
        for v in range(1, len(self.vertices())):
            p = self.vertex(v)
            q = other.vertex(v)
            sign = p[0] - q[0]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
            sign = p[1] - q[1]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
        return 0

    def translate(self, u):
        r"""
        Return a copy of this polygon that has been translated by ``u``.

        TESTS::

            sage: from flatsurf import Polygon
            sage: Polygon(vertices=[(0,0), (2,0), (1,1)]).translate((3,-2))
            Polygon(vertices=[(3, -2), (5, -2), (4, -1)])

        """
        # TODO: Require in base class. Implement partially in base class. Probably more generally as apply_matrix or something like that.
        return self.parent().polygon(
            edges=[e.translate(u) for e in self._edges],
            check=False
        )

    def _apply_scalar(self, r):
        if not r:
            return self.parent().point()(0, 0)
        elif r > 0:
            return self.parent().polygon(
                edges=[e._apply_scalar(r) for e in self._edges],
                check=False)
        else:
            return self.parent().polygon(
                edges=[e._apply_scalar(r) for e in reversed(self._edges)],
                check=False)

    def _apply_similarity(self, g):
        if g.sign().is_one():
            return self.parent().polygon(
                edges=[e._apply_similarity(g) for e in self._edges],
                check=False)
        else:
            return self.parent().polygon(
                edges=[e._apply_similarity(g) for e in reversed(self._edges)],
                check=False)

    def _apply_2x2_matrix(self, m):
        if m.det() > 0:
            return self.parent().polygon(
                edges=[e._apply_2x2_matrix(m) for e in self._edges],
                check=False)
        else:
            # TODO: set check=False below
            return self.parent().polygon(
                edges=[e._apply_2x2_matrix(m) for e in reversed(self._edges)],
                check=True)

    def _apply_3x3_matrix(self, m):
        if m.det() > 0:
            return self.parent().polygon(
                edges=[e._apply_3x3_matrix(m) for e in self._edges],
                check=False)
        else:
            # TODO: set check=False below
            return self.parent().polygon(
                edges=[e._apply_3x3_matrix(m) for e in reversed(self._edges)],
                check=True)

    def _repr_(self):
        r"""
        Return a printable representation of this polygon.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: S
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        try:
            return f"Polygon(vertices={repr(list(self.vertices()))})"
        except (AttributeError, NotImplementedError):
            return f"Polygon(edges={repr(list(self.sides()))})"

    def vertices(self, marked_vertices=True):
        # TODO: Expose as points as well.
        r"""
        Return the vertices of this polygon in counterclockwise order as
        vectors in the real plane.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``); whether to
          include vertices with a π angle in the output.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()
            sage: s.vertices()
            ((0, 0), (1, 0), (1, 1), (0, 1))

        """
        if not marked_vertices:
            return tuple(
                vertex
                # TODO: Test this code path.
                for (vertex, slope) in zip(self.vertices(), self.slopes(relative=True))
                if slope[1] != 0
            )

        vertices = [e.start() for e in self.sides()]
        if any(v is None for v in vertices):
            raise NotImplementedError("some sides of the polygon do not end in a finite point")

        return tuple(v.vector() for v in vertices)

    def side(self, e):
        return self._edges[e % len(self._edges)]

    def sides(self):
        return self._edges

    def corners(self):
        # TODO: Make this faster by caching?
        corners = []

        n = len(self._edges)
        for i in range(n):
            vertex = None
            if isinstance(self._edges[i], EuclideanSegment):
                # TODO: There must be a public way to implement this without reaching into the implementation details.
                vertex = self._edges[i]._start
            corners.append((-self._edges[i - 1].direction(), vertex, self._edges[i].direction()))

        return corners

    def _check(self, require_normalized=True):
        r"""
        Verify that this is a valid polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import _Polygon_check, Polygon
            sage: p = Polygon(angles=[1, 1, 1])
            sage: _Polygon_check(p, vertices=None, edges=None, angles=[1, 1, 1], lengths=None, convex=None)

        """
        # Check that the polygon satisfies the assumptions of EuclideanPolygon
        area = self.area()

        if area < 0:
            raise ValueError(
                "polygon has negative area; probably the vertices are not in counter-clockwise order"
            )

        if require_normalized:
            if area == 0:
                raise ValueError("polygon has zero area")

            if any(edge.dimension() == 0 for edge in self.sides()):
                raise ValueError("polygon has zero edge")

            for v, corner, w in self.corners():
                from flatsurf.geometry.euclidean import is_parallel

                if corner is None:
                    continue

                if is_parallel(v, w):
                    raise ValueError("polygon has anti-parallel edges")

        sides = self.sides()
        for i in range(len(sides)):
            if sides[i - 1].end() != sides[i].start():
                raise ValueError("polygon is not closed")

        from flatsurf.geometry.categories import EuclideanPolygons

        if not EuclideanPolygons.ParentMethods.is_simple(self):
            raise NotImplementedError("polygon self-intersects")

    def change(self, ring=None, geometry=None, oriented=None):
        if ring is not None or geometry is not None:
            P = self.parent()
            Q = P.change_ring(ring, geometry=geometry)
            if P is not Q:
                # TODO: this is not good enough to rebuild a polygon
                self = Q.polygon(vertices=self.vertices())

        if oriented:
            raise ValueError("polygons cannot be oriented")

        return self

    def __contains__(self, point):
        # TODO: make a proper implementation
        if point.is_ideal():
            # TODO: include point at infinity in the polygon
            return False
        return self.get_point_position(point).is_inside()


class EuclideanDistance_base(Element):
    def _acted_upon_(self, other, self_on_left):
        return self.scale(other)

    def scale(self, scalar):
        scalar = self.parent().base_ring()(scalar)

        if scalar < 0:
            raise ValueError("cannot scale a distance by a negative scalar")

        return self._scale(scalar)

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE, op_LE, op_LT, op_GE, op_GT

        if op == op_LT:
            return self <= other and not (self >= other)
        if op == op_LE:
            return self._le_(other)
        if op == op_EQ:
            return self._eq_(other)
        if op == op_NE:
            return not self == other
        if op == op_GT:
            return self >= other and not (self <= other)
        if op == op_GE:
            return self._ge_(other)

        raise NotImplementedError("Operator not implemented for this distance")

    def _le_(self, other):
        if self.is_finite() and not other.is_finite():
            return True

        return self.norm_squared() <= other.norm_squared()

    def _ge_(self, other):
        if self.is_finite() and not other.is_finite():
            return False

        return self.norm_squared() >= other.norm_squared()

    def _eq_(self, other):
        return self.norm_squared() == other.norm_squared()

    def _div_(self, other):
        return self.parent().from_quotient(self, other)

    def __float__(self):
        from math import sqrt

        f = float(self.norm_squared())
        if f < 0:
            print("Bug https://github.com/sagemath/sage/issues/37983.")
            return float(0)
        return sqrt(f)


class EuclideanDistance_squared(EuclideanDistance_base):
    def __init__(self, parent, norm_squared):
        super().__init__(parent)

        self._norm_squared = parent.base_ring()(norm_squared)

    def norm_squared(self):
        return self._norm_squared

    def norm(self):
        from sage.all import AA

        return AA(self._norm_squared).sqrt()

    def _scale(self, scalar):
        return self.parent().from_norm_squared(scalar**2 * self._norm_squared)

    def is_finite(self):
        return True

    def _add_(self, other):
        return self.parent().from_sum(self, other)

    def _repr_(self):
        try:
            return repr(self.norm())
        except (NotImplementedError, ValueError):
            pass

        return f"√{self.norm_squared()}"


class EuclideanDistance_infinite(EuclideanDistance_base):
    def _scale(self, scalar):
        if scalar == 0:
            raise NotImplementedError

        return self

    def norm_squared(self):
        from sage.all import oo

        return oo

    def is_finite(self):
        return False

    def _repr_(self):
        return "∞"


class EuclideanDistances(Parent):
    def __init__(self, euclidean_plane, category=None):
        # TODO: Pick a better category.
        from sage.categories.all import Sets

        super().__init__(euclidean_plane.base_ring(), category=category or Sets())
        self._euclidean_plane = euclidean_plane

    def _repr_(self):
        return f"Euclidean Norm on {self._euclidean_plane}"

    @cached_method
    def infinite(self):
        r"""
        Return an infinite distance.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.norm().infinite()
            ∞

        """
        return self.__make_element_class__(EuclideanDistance_infinite)(self)

    def zero(self):
        return self.from_norm_squared(0)

    def from_norm_squared(self, x):
        return self.__make_element_class__(EuclideanDistance_squared)(self, x)

    def from_vector(self, v):
        v = self._euclidean_plane.vector_space()(v)
        return self.from_norm_squared(v.dot_product(v))

    def _element_constructor_(self, x):
        from sage.all import parent

        if self.base_ring().has_coerce_map_from(parent(x)):
            return self.from_norm_squared(x**2)

        return self.from_vector(x)


### TODO: PRE-EUCLIDEAN-PLANE CODE HERE


def is_cosine_sine_of_rational(cos, sin, scaled=False):
    r"""
    Check whether the given pair is a cosine and sine of a same rational angle.

    INPUT:

    - ``cos`` -- a number

    - ``sin`` -- a number

    - ``scaled`` -- a boolean (default: ``False``); whether to allow ``cos``
      and ``sin`` to be scaled by the same positive algebraic number

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_cosine_sine_of_rational

        sage: c = s = AA(sqrt(2))/2
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: c = AA(sqrt(3))/2
        sage: s = AA(1/2)
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: c = AA(sqrt(5)/2)
        sage: s = (1 - c**2).sqrt()
        sage: c**2 + s**2
        1.000000000000000?
        sage: is_cosine_sine_of_rational(c, s)
        False

        sage: c = (AA(sqrt(5)) + 1)/4
        sage: s = (1 - c**2).sqrt()
        sage: is_cosine_sine_of_rational(c, s)
        True

        sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
        sage: is_cosine_sine_of_rational(K.zero(), -K.one())
        True

    TESTS::

        sage: from pyexactreal import ExactReals # optional: pyexactreal  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: R = ExactReals() # optional: pyexactreal
        sage: is_cosine_sine_of_rational(R.one(), R.zero()) # optional: pyexactreal
        True

    """
    from sage.all import AA

    # We cannot check in AA due to https://github.com/flatsurf/exact-real/issues/172
    # We just trust that non-algebraic elements won't allow conversion to AA.
    # if cos not in AA:
    #     return False
    # if sin not in AA:
    #     return False

    if not scaled:
        if cos**2 + sin**2 != 1:
            return False

    try:
        cos = AA(cos)
    except ValueError:
        # This is a replacement for the "in AA" checked disabled above.
        return False

    cos = cos.as_number_field_element(embedded=True)
    # We need an explicit conversion to the number field due to https://github.com/sagemath/sage/issues/35613
    cos = cos[0](cos[1])

    try:
        sin = AA(sin)
    except ValueError:
        # This is a replacement for the "in AA" checked disabled above.
        return False
    sin = sin.as_number_field_element(embedded=True)
    # We need an explicit conversion to the number field due to https://github.com/sagemath/sage/issues/35613
    sin = sin[0](sin[1])

    from sage.all import ComplexBallField

    CBF = ComplexBallField(53)

    x = CBF(cos) + CBF.gen(0) * CBF(sin)
    xN = x

    # Suppose that (cos, sin) are indeed sine and cosine of a rational angle.
    # Then x = cos + I*sin generates a cyclotomic field C and for some N we
    # have x^N = ±1. Since C is contained in the compositum of K=Q(cos) and
    # L=Q(i*sin) and Q(cos) and Q(sin) are both contained in C, the degree of C
    # is bounded from above by twice (accounting for the imaginary unit) the
    # degrees of K and L. The degree of C is the totient of N which is bounded
    # from below by n / (e^γ loglog n + 3 / loglog n) [cf. wikipedia].
    degree_bound = 2 * cos.minpoly().degree() * sin.minpoly().degree()

    from itertools import count

    for n in count(2):
        xN *= x

        c = xN.real()
        s = xN.imag()

        if xN.real().contains_zero() or xN.imag().contains_zero():
            c, s = cos, sin
            for _ in range(n - 1):
                c, s = c * cos - s * sin, s * cos + c * sin

            if c == 0 or s == 0:
                return True

            CBF = ComplexBallField(CBF.precision() * 2)
            x = CBF(cos) + CBF.gen(0) * CBF(sin)
            xN = x**n

        from math import log

        if n / (2.0 * log(log(n)) + 3 / log(log(n))) > 2 * degree_bound:
            return False


def acos(cos_angle, numerical=False):
    r"""
    Return the arccosine of ``cos_angle`` as a multiple of 2π, i.e., as a value
    between 0 and 1/2.

    INPUT:

    - ``cos_angle`` -- a floating point number, the cosine of an angle

    - ``numerical`` -- a boolean (default: ``False``); whether to return a
      numerical approximation of the arccosine or try to reconstruct an exact
      rational value for the arccosine (in radians.)

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import acos

        sage: acos(1)
        0
        sage: acos(.5)
        1/6
        sage: acos(0)
        1/4
        sage: acos(-.5)
        1/3
        sage: acos(-1)
        1/2

        sage: acos(.25)
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot recover a rational angle from these numerical results
        sage: acos(.25, numerical=True)
        0.2097846883724169

    """
    import math

    angle = math.acos(cos_angle) / (2 * math.pi)

    assert 0 <= angle <= 0.5

    if numerical:
        return angle

    # fast and dirty way using floating point approximation
    from sage.all import RR

    angle_rat = RR(angle).nearby_rational(0.00000001)
    if angle_rat.denominator() > 256:
        raise NotImplementedError(
            "cannot recover a rational angle from these numerical results"
        )
    return angle_rat


def angle(u, v, numerical=False):
    r"""
    Return the angle between the vectors ``u`` and ``v`` divided by `2 \pi`.

    INPUT:

    - ``u``, ``v`` - vectors in the plane

    - ``numerical`` - boolean (default: ``False``), whether to return floating
      point numbers

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import angle

    As the implementation is dirty, we at least check that it works for all
    denominator up to 20::

        sage: u = vector((AA(1),AA(0)))
        sage: for n in xsrange(1,20):       # long time  (1.5s)
        ....:     for k in xsrange(1,n):
        ....:         v = vector((AA(cos(2*k*pi/n)), AA(sin(2*k*pi/n))))
        ....:         assert angle(u,v) == k/n

    The numerical version (working over floating point numbers)::

        sage: import math
        sage: u = (1, 0)
        sage: for n in xsrange(1,20):
        ....:     for k in xsrange(1,n):
        ....:         a = 2 * k * math.pi / n
        ....:         v = (math.cos(a), math.sin(a))
        ....:         assert abs(angle(u,v,numerical=True) * 2 * math.pi - a) < 1.e-10

    If the angle is not rational, then the method returns an element in the real
    lazy field::

        sage: v = vector((AA(sqrt(2)), AA(sqrt(3))))
        sage: a = angle(u, v)
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot recover a rational angle from these numerical results
        sage: a = angle(u, v, numerical=True)
        sage: a    # abs tol 1e-14
        0.14102355421224375
        sage: exp(2*pi.n()*CC(0,1)*a)
        0.632455532033676 + 0.774596669241483*I
        sage: v / v.norm()
        (0.6324555320336758?, 0.774596669241484?)

    """
    import math

    u0 = float(u[0])
    u1 = float(u[1])
    v0 = float(v[0])
    v1 = float(v[1])

    cos_uv = (u0 * v0 + u1 * v1) / math.sqrt((u0 * u0 + u1 * u1) * (v0 * v0 + v1 * v1))
    if cos_uv < -1.0:
        assert cos_uv > -1.0000001
        cos_uv = -1.0
    elif cos_uv > 1.0:
        assert cos_uv < 1.0000001
        cos_uv = 1.0

    angle = acos(cos_uv, numerical=numerical)
    return 1 - angle if u0 * v1 - u1 * v0 < 0 else angle


def ccw(v, w):
    r"""
    Return a positive number if the turn from ``v`` to ``w`` is
    counterclockwise, a negative number if it is clockwise, and zero if the two
    vectors are collinear.

    .. NOTE::

        This function is sometimes also referred to as the wedge product or
        simply the determinant. We chose the more customary name ``ccw`` from
        computational geometry here.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import ccw
        sage: ccw((1, 0), (0, 1))
        1
        sage: ccw((1, 0), (-1, 0))
        0
        sage: ccw((1, 0), (0, -1))
        -1
        sage: ccw((1, 0), (1, 0))
        0

    """
    return v[0] * w[1] - v[1] * w[0]


def is_parallel(v, w):
    r"""
    Return whether the vectors ``v`` and ``w`` are parallel (but not
    anti-parallel.)

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_parallel
        sage: is_parallel((0, 1), (0, 1))
        True
        sage: is_parallel((0, 1), (0, 2))
        True
        sage: is_parallel((0, 1), (0, -2))
        False
        sage: is_parallel((0, 1), (0, 0))
        False
        sage: is_parallel((0, 1), (1, 0))
        False

    TESTS::

        sage: V = QQ**2

        sage: is_parallel(V((0,1)), V((0,2)))
        True
        sage: is_parallel(V((1,-1)), V((2,-2)))
        True
        sage: is_parallel(V((4,-2)), V((2,-1)))
        True
        sage: is_parallel(V((1,2)), V((2,4)))
        True
        sage: is_parallel(V((0,2)), V((0,1)))
        True

        sage: is_parallel(V((1,1)), V((1,2)))
        False
        sage: is_parallel(V((1,2)), V((2,1)))
        False
        sage: is_parallel(V((1,2)), V((1,-2)))
        False
        sage: is_parallel(V((1,2)), V((-1,-2)))
        False
        sage: is_parallel(V((2,-1)), V((-2,1)))
        False

    """
    if ccw(v, w) != 0:
        # vectors are not collinear
        return False

    return v[0] * w[0] + v[1] * w[1] > 0


def is_anti_parallel(v, w):
    r"""
    Return whether the vectors ``v`` and ``w`` are anti-parallel, i.e., whether
    ``v`` and ``-w`` are parallel.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_anti_parallel
        sage: V = QQ**2

        sage: is_anti_parallel(V((0,1)), V((0,-2)))
        True
        sage: is_anti_parallel(V((1,-1)), V((-2,2)))
        True
        sage: is_anti_parallel(V((4,-2)), V((-2,1)))
        True
        sage: is_anti_parallel(V((-1,-2)), V((2,4)))
        True

        sage: is_anti_parallel(V((1,1)), V((1,2)))
        False
        sage: is_anti_parallel(V((1,2)), V((2,1)))
        False
        sage: is_anti_parallel(V((0,2)), V((0,1)))
        False
        sage: is_anti_parallel(V((1,2)), V((1,-2)))
        False
        sage: is_anti_parallel(V((1,2)), V((-1,2)))
        False
        sage: is_anti_parallel(V((2,-1)), V((-2,-1)))
        False

    """
    return is_parallel(v, -w)


def line_intersection(l, m):
    r"""
    Return the point of intersection between the lines ``l`` and ``m``. If the
    lines do not have a single point of intersection, returns None.

    INPUT:

    - ``l`` -- a line in the plane given by two points (as vectors) in the plane

    - ``m`` -- a line in the plane given by two points (as vectors) in the plane

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import line_intersection
        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((0, -1)), vector((0, 1))))
        (0, 0)

    Parallel lines have no single point of intersection::

        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((-1, 1)), vector((1, 1))))

    Identical lines have no single point of intersection::

        sage: line_intersection((vector((-1, 0)), vector((1, 0))), (vector((-2, 0)), vector((2, 0))))

    """
    Δl = l[1] - l[0]
    Δm = m[1] - m[0]

    # We solve the linear system determining the time when l[0] + t * Δl hits
    # m, i.e., l[0] + t Δl == m[0] + s Δm.
    a, b, c, d = Δl[0], -Δm[0], Δl[1], -Δm[1]
    rhs = m[0] - l[0]

    det = a * d - b * c
    if det == 0:
        # The lines are parallel
        return None

    # Solve the linear system. We only need t (and not s.)
    t = (d * rhs[0] - b * rhs[1]) / det

    return l[0] + t * Δl


def time_on_ray(p, direction, q):
    if direction[0]:
        dim = 0
    else:
        dim = 1

    delta = q[dim] - p[dim]
    length = direction[dim]
    if length < 0:
        delta *= -1
        length *= -1

    return delta, length


def time_on_segment(segment, p):
    if p == segment[0]:
        return 0
    if not is_parallel(p - segment[0], segment[1] - segment[0]):
        return None

    delta, length = time_on_ray(segment[0], segment[1] - segment[0], p)
    if delta > length:
        return None
    if delta < 0:
        return None

    return delta / length


def ray_segment_intersection(p, direction, segment):
    r"""
    Return the intersection of the ray from ``p`` in ``direction`` with
    ``segment``.

    If the segment and the ray intersect in a point, return that point as a
    vector.

    If the segment and the ray overlap in a segment, return the end points of
    that segment (in order.)

    If the segment and the ray do not intersect, return ``None``.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import ray_segment_intersection
        sage: V = QQ**2

        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((1, -1)), V((1, 1))))
        (1, 0)
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((0, 0)), V((1, 0))))
        ((0, 0), (1, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((1, 0)), V((2, 0))))
        ((1, 0), (2, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((-1, 0)), V((1, 0))))
        ((0, 0), (1, 0))
        sage: ray_segment_intersection(V((0, 0)), V((1, 0)), (V((-1, -1)), V((-1, 1))))

    TESTS::

        sage: ray_segment_intersection(V((0, 0)), V((5, 1)), (V((3, 2)), V((3, 3))))

    """
    intersection = line_intersection((p, p + direction), segment)

    if intersection is None:
        # ray and segment are parallel.
        if ccw(direction, segment[0] - p) != 0:
            # ray and segment are not on the same line
            return None

        t0, length = time_on_ray(p, direction, segment[0])
        t1, _ = time_on_ray(p, direction, segment[1])

        if t1 < t0:
            t0, t1 = t1, t0

        if t1 < 0:
            return None

        if t1 == 0:
            return p

        t0 /= length
        t1 /= length

        if t0 < 0:
            return (p, p + t1 * direction)

        return (p + t0 * direction, p + t1 * direction)

    if time_on_ray(p, direction, intersection)[0] < 0:
        return None

    if ccw(segment[0] - p, direction) * ccw(segment[1] - p, direction) > 0:
        return None

    return intersection


def is_box_intersecting(b, c):
    r"""
    Return whether the (bounding) boxes ``b`` and ``c`` intersect.

    INPUT:

    - ``b`` -- a pair of corners
    - ``c`` -- a pair of corners

    OUTPUT:

    - ``0`` -- do not intersect
    - ``1`` -- intersect in a point
    - ``2`` -- intersect in a segment
    - ``3`` -- intersection has interior points

    """

    def normalize(b):
        x_inverted = b[0][0] > b[1][0]
        y_inverted = b[0][1] > b[1][1]

        return (
            (b[1][0] if x_inverted else b[0][0], b[1][1] if y_inverted else b[0][1]),
            (b[0][0] if x_inverted else b[1][0], b[0][1] if y_inverted else b[1][1]),
        )

    b = normalize(b)
    c = normalize(c)

    if b[0][0] > c[1][0]:
        return 0
    if c[0][0] > b[1][0]:
        return 0
    if b[0][1] > c[1][1]:
        return 0
    if c[0][1] > b[1][1]:
        return 0

    # TODO: This is not the full algorithm yet.
    return 3


def is_segment_intersecting(s, t):
    r"""
    Return whether the segments ``s`` and ``t`` intersect.

    INPUT:

    - ``s`` -- a segment given as a pair of endpoints (given as vectors in the plane.)

    - ``t`` -- a segment given as a pair of endpoints (given as vectors in the plane.)

    OUTPUT:

    - ``0`` - do not intersect
    - ``1`` - exactly one endpoint in common
    - ``2`` - non-trivial intersection

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_segment_intersecting
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, 1)), vector((0, 3))))
        0
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, 0)), vector((0, 3))))
        1
        sage: is_segment_intersecting((vector((0, 0)), vector((1, 0))), (vector((0, -1)), vector((0, 3))))
        2
        sage: is_segment_intersecting((vector((-1, -1)), vector((1, 1))), (vector((0, 0)), vector((2, 2))))
        2
        sage: is_segment_intersecting((vector((-1, -1)), vector((1, 1))), (vector((1, 1)), vector((2, 2))))
        1

    """
    if not is_box_intersecting(s, t):
        return 0

    Δs = s[1] - s[0]
    turn_from_s = ccw(Δs, t[0] - s[0]) * ccw(Δs, t[1] - s[0])
    if turn_from_s > 0:
        # Both endpoints of t are on the same side of s
        return 0

    Δt = t[1] - t[0]
    turn_from_t = ccw(Δt, s[0] - t[0]) * ccw(Δt, s[1] - t[0])
    if turn_from_t > 0:
        # Both endpoints of s are on the same side of t
        return 0

    if ccw(Δs, Δt) == 0:
        # Segments are parallel
        if is_anti_parallel(Δs, Δt):
            # Ensure segments are oriented the same way.
            t = (t[1], t[0])
            Δt = -Δt

        time_to_t0, length = time_on_ray(s[0], Δs, t[0])
        time_to_t1, _ = time_on_ray(s[0], Δs, t[1])

        if time_to_t0 < 0 and time_to_t1 < 0:
            # Both endpoints of t are earlier than s
            return 0

        if time_to_t0 > length and time_to_t1 > length:
            # Both endpoints of t are later than s
            return 0

        if time_to_t0 == length:
            return 1

        if time_to_t1 == 0:
            return 1

        return 2

    if turn_from_t == 0 and turn_from_s == 0:
        return 1

    return 2


def is_between(begin, end, v):
    r"""
    Check whether the vector ``v`` is strictly in the sector formed by the vectors
    ``begin`` and ``end`` (in counter-clockwise order).

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import is_between
        sage: is_between((1, 0), (1, 1), (2, 1))
        True

        sage: from itertools import product
        sage: vecs = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
        sage: for (i, vi), (j, vj), (k, vk) in product(enumerate(vecs), repeat=3):
        ....:     assert is_between(vi, vj, vk) == ((i == j and i != k) or i < k < j or k < j < i or j < i < k), ((i, vi), (j, vj), (k, vk))
    """
    if begin[0] * end[1] > end[0] * begin[1]:
        # positive determinant
        # [ begin[0] end[0] ]^-1 = [ end[1] -end[0] ]
        # [ begin[1] end[1] ]      [-begin[1]  begin[0] ]
        # v[0] * end[1] - end[0] * v[1] > 0
        # - v[0] * begin[1] + begin[0] * v[1] > 0
        return end[1] * v[0] > end[0] * v[1] and begin[0] * v[1] > begin[1] * v[0]
    elif begin[0] * end[1] == end[0] * begin[1]:
        # aligned vector
        if begin[0] * end[0] >= 0 and begin[1] * end[1] >= 0:
            return (
                v[0] * begin[1] != v[1] * begin[0]
                or v[0] * begin[0] < 0
                or v[1] * begin[1] < 0
            )
        else:
            return begin[0] * v[1] > begin[1] * v[0]
    else:
        # negative determinant
        # [ end[0] begin[0] ]^-1 = [ begin[1] -begin[0] ]
        # [ end[1] begin[1] ]      [-end[1]  end[0] ]
        # v[0] * begin[1] - begin[0] * v[1] > 0
        # - v[0] * end[1] + end[0] * v[1] > 0
        return begin[1] * v[0] < begin[0] * v[1] or end[0] * v[1] < end[1] * v[0]


def solve(x, u, y, v):
    r"""
    Return (a,b) so that: x + au = y + bv

    INPUT:

    - ``x``, ``u``, ``y``, ``v`` -- two dimensional vectors

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import solve
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: V = VectorSpace(K,2)
        sage: x = V((1,-sqrt2))
        sage: y = V((1,1))
        sage: a = V((0,1))
        sage: b = V((-sqrt2, sqrt2+1))
        sage: u = V((0,1))
        sage: v = V((-sqrt2, sqrt2+1))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True

        sage: u = V((1,1))
        sage: v = V((1,sqrt2))
        sage: a, b = solve(x,u,y,v)
        sage: x + a*u == y + b*v
        True

    """
    d = -u[0] * v[1] + u[1] * v[0]
    if d.is_zero():
        raise ValueError("parallel vectors")
    a = v[1] * (x[0] - y[0]) + v[0] * (y[1] - x[1])
    b = u[1] * (x[0] - y[0]) + u[0] * (y[1] - x[1])
    return (a / d, b / d)


def projectivization(x, y, signed=True, denominator=None):
    r"""
    Return a simplified version of the projective coordinate [x: y].

    If ``signed`` (the default), the second coordinate is made non-negative;
    otherwise the coordinates keep their signs.

    If ``denominator`` is ``False``, returns [x/y: 1] up to sign. Otherwise,
    the returned coordinates have no denominator and no non-unit gcd.

    TESTS::

        sage: from flatsurf.geometry.euclidean import projectivization

        sage: projectivization(2/3, -3/5, signed=True, denominator=True)
        (10, -9)
        sage: projectivization(2/3, -3/5, signed=False, denominator=True)
        (-10, 9)
        sage: projectivization(2/3, -3/5, signed=True, denominator=False)
        (10/9, -1)
        sage: projectivization(2/3, -3/5, signed=False, denominator=False)
        (-10/9, 1)

        sage: projectivization(-1/2, 0, signed=True, denominator=True)
        (-1, 0)
        sage: projectivization(-1/2, 0, signed=False, denominator=True)
        (1, 0)
        sage: projectivization(-1/2, 0, signed=True, denominator=False)
        (-1, 0)
        sage: projectivization(-1/2, 0, signed=False, denominator=False)
        (1, 0)

    """
    from sage.all import Sequence

    parent = Sequence([x, y]).universe()
    if y:
        z = x / y
        if denominator is True or (denominator is None and hasattr(z, "denominator")):
            d = parent(z.denominator())
        else:
            d = parent(1)
        if signed and y < 0:
            d *= -1
        return (z * d, d)
    elif signed and x < 0:
        return (parent(-1), parent(0))
    else:
        return (parent(1), parent(0))


def slope(a, rotate=1):
    r"""
    Return either ``1`` (positive slope) or ``-1`` (negative slope).

    If ``rotate`` is set to 1 then consider the edge as if it was rotated counterclockwise
    infinitesimally.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import slope
        sage: slope((1, 1))
        1
        sage: slope((-1, 1))
        -1
        sage: slope((-1, -1))
        1
        sage: slope((1, -1))
        -1

        sage: slope((1, 0))
        1
        sage: slope((0, 1))
        -1
        sage: slope((-1, 0))
        1
        sage: slope((0, -1))
        -1

        sage: slope((1, 0), rotate=-1)
        -1
        sage: slope((0, 1), rotate=-1)
        1
        sage: slope((-1, 0), rotate=-1)
        -1
        sage: slope((0, -1), rotate=-1)
        1

        sage: slope((1, 0), rotate=0)
        0
        sage: slope((0, 1), rotate=0)
        0
        sage: slope((-1, 0), rotate=0)
        0
        sage: slope((0, -1), rotate=0)
        0

        sage: slope((0, 0))
        Traceback (most recent call last):
        ...
        ValueError: zero vector
    """
    x, y = a
    if not x and not y:
        raise ValueError("zero vector")
    if (x > 0 and y > 0) or (x < 0 and y < 0):
        return 1
    elif (x > 0 and y < 0) or (x < 0 and y > 0):
        return -1
    if rotate == 0:
        return 0
    if rotate == 1:
        return 1 if x else -1
    if rotate == -1:
        return 1 if y else -1
    raise ValueError("invalid argument rotate={}".format(rotate))


def to_3x3_matrix(g):
    from sage.matrix.constructor import matrix
    from sage.structure.element import Matrix
    from flatsurf.geometry.similarity import Similarity
    if isinstance(g, Similarity):
        return g.matrix()
    elif isinstance(g, Matrix):
        if g.nrows() == g.ncols() == 2:
            return matrix(g.base_ring(), 3, [g[0, 0], g[0, 1], 0, g[1, 0], g[1, 1], 0, 0, 0, 1])
        elif g.nrows() == 2 and g.ncols() == 3:
            return matrix(g.base_ring(), 3, [g[0, 0], g[0, 1], g[0, 2], g[1, 0], g[1, 1], g[1, 2], 0, 0, 1])
        elif g.nrows() != 3 or g.ncols() != 3:
            raise ValueError(f"invalid element g={g}")
    else:
        raise TypeError(f"invalid element g={g} of type={type(g)}")
