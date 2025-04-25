r"""
Generic geometric notions that are independent of a concrete model of geometry
such as the hyperbolic or Euclidean geometry.

EXAMPLES:

This module provides generic geometries such as an "exact" geometry and a
numerical geometry. Concrete geometries then inherit some functionality from
these generic geometries::

    sage: from flatsurf.geometry.hyperbolic import HyperbolicExactGeometry
    sage: from flatsurf.geometry.geometry import ExactGeometry

    sage: G = HyperbolicExactGeometry(QQ)
    sage: isinstance(G, ExactGeometry)
    True

    sage: G._zero(1/2**64)
    False

::

    sage: from flatsurf.geometry.hyperbolic import HyperbolicEpsilonGeometry
    sage: from flatsurf.geometry.geometry import EpsilonGeometry

    sage: G = HyperbolicEpsilonGeometry(QQ, 1e-9)
    sage: isinstance(G, EpsilonGeometry)
    True

    sage: G._zero(1/2**64)
    True

"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023-2025 Julian Rüth
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


class Geometry:
    r"""
    Predicates and primitive geometric constructions over a base ``ring``.

    This is an abstract base class to collect shared functionality for concrete
    geometries such as the :class:`EuclideanGeometry` and the
    :class:`HyperbolicGeometry`.

    INPUT:

    - ``ring`` -- a ring, the ring in which object in this geometry will be
      represented

    TESTS::

        sage: from flatsurf.geometry.euclidean import EuclideanExactGeometry, Geometry
        sage: geometry = EuclideanExactGeometry(QQ)
        sage: isinstance(geometry, Geometry)
        True

    """

    def __init__(self, ring):
        r"""
        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: from flatsurf.geometry.euclidean import EuclideanGeometry
            sage: E = EuclideanPlane()
            sage: isinstance(E.geometry, EuclideanGeometry)
            True

        """
        self._ring = ring

    def base_ring(self):
        r"""
        Return the ring over which this geometry is implemented.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.geometry.base_ring()
            Rational Field

        """
        return self._ring

    def change_ring(self, ring):
        r"""
        Return this geometry with the :meth:`base_ring` changed to ``ring``.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.geometry
            Exact geometry over Rational Field
            sage: E.geometry.change_ring(AA)
            Exact geometry over Algebraic Real Field

        ::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry
            Exact geometry over Rational Field
            sage: H.geometry.change_ring(AA)
            Exact geometry over Algebraic Real Field

        """
        raise NotImplementedError("this geometry does not implement change_ring()")

    def interval_intersection(self, i, j):
        r"""
        Return the intersection of the two intervals of the real line ``i`` and
        ``j``.

        INPUT:

        - ``i`` -- a closed interval of the real line encoded as ``None`` for
          the empty interval, or a pair ``(left, right)`` where an entry
          ``None`` encodes that the interval is unbounded.

        - ``j`` -- a closed interval encoded like ``i``

        """
        if i is None or j is None:
            return None

        if i[0] is None:
            left = j[0]
        elif j[0] is None:
            left = i[0]
        else:
            left = max(i[0], j[0])

        if i[1] is None:
            right = j[1]
        elif j[1] is None:
            right = i[1]
        else:
            right = min(i[1], j[1])

        if left is not None and right is not None:
            if left > right:
                return None

        return (left, right)

    def interval_element(self, i):
        r"""
        Return an element in the interval ``i``.

        INPUT:

        - ``i`` -- a closed interval of the real line, see
          :meth:`interval_intersection` for details.

        """
        if i is None:
            raise ValueError("empty interval has no elements")

        if i[0] is None:
            return i[1]

        if i[1] is None:
            return i[0]

        return (i[0] + i[1]) / 2

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

    def _segment_intersection(self, f, g):
        Δf = (f[1][0] - f[0][0], f[1][1] - f[0][1])
        # TODO: Bring the speedups from is_segment_intersecting here.
        ccwfg0 = self._ccw(Δf, (g[0][0] - f[0][0], g[0][1] - f[0][1]))
        ccwfg1 = self._ccw(Δf, (g[1][0] - f[0][0], g[1][1] - f[0][1]))

        if ccwfg0 * ccwfg1 > 0:
            return None

        Δg = (g[1][0] - g[0][0], g[1][1] - g[0][1])
        ccwgf0 = self._ccw(Δg, (f[0][0] - g[0][0], f[0][1] - g[0][1]))
        ccwgf1 = self._ccw(Δg, (f[1][0] - g[0][0], f[1][1] - g[0][1]))

        if ccwgf0 * ccwgf1 > 0:
            return None

        if ccwfg0 == ccwfg1:
            # The segments are parallel
            assert ccwfg0 == ccwfg1 == ccwgf0 == ccwgf1 == 0

            t0 = self._time_on_ray(f[0], Δf, g[0])
            t1 = self._time_on_ray(f[0], Δf, g[1])

            if t0[0] * t1[1] > t1[0] * t0[1]:
                g = g[1], g[0]
                Δg = -Δg[0], -Δg[1]
                t0, t1 = t1, t0

            if t0[0] <= 0:
                if t1[0] < 0:
                    return None

                a = f[0]
                if t1[0] <= t1[1]:
                    b = g[1]
                else:
                    b = f[1]

            else:
                if t0[0] > t0[1]:
                    return None

                a = g[0]
                if t1[0] < t1[1]:
                    b = g[1]
                else:
                    b = f[1]

            if a == b:
                return ((a[0], 1), (a[1], 1))

            return (a[0], a[1], b[0], b[1])

        fb, fc = (-(f[1][1] - f[0][1]), f[1][0] - f[0][0])
        gb, gc = (-(g[1][1] - g[0][1]), g[1][0] - g[0][0])

        fa = -(fb * f[0][0] + fc * f[0][1])
        ga = -(gb * g[0][0] + gc * g[0][1])

        return self._line_intersection((fa, fb, fc), (ga, gb, gc))

    def _time_on_ray(self, p, direction, q):
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

    def _ccw(self, v, w):
        r"""
        Return a positive number if the turn from ``v`` to ``w`` is
        counterclockwise, a negative number if it is clockwise, and zero if the two
        vectors are collinear.

        .. NOTE::

            This function is sometimes also referred to as the wedge product or
            simply the determinant. We chose the more customary name ``ccw`` from
            computational geometry here.

        # TODO
        # EXAMPLES::

        #     sage: from flatsurf.geometry.euclidean import ccw
        #     sage: ccw((1, 0), (0, 1))
        #     1
        #     sage: ccw((1, 0), (-1, 0))
        #     0
        #     sage: ccw((1, 0), (0, -1))
        #     -1
        #     sage: ccw((1, 0), (1, 0))
        #     0

        """
        # TODO: Use sgn to figure out whether this is too close to zero.
        return v[0] * w[1] - v[1] * w[0]

    def segment_intersection(self, f, g):
        r"""
        Return the intersection between the Euclidean segments ``f`` and ``g``.

        INPUT:

        - ``f`` -- a pair of points ``((x0, y0), (x1, y1))``, the end points of
          the segment.

        - ``g`` -- a pair of points ``((x0, y0), (x1, y1))``, the end points of
          the segment.

        OUTPUT: ``None`` if the segments do not intersect, a pair of
        coordinates if the segments intersect in a single point, or the end
        points of a segment ``(x0, y0, x1, y1)`` if the segments intersect in a
        segment.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.segment_intersection(((-1, 0), (1, 0)), ((0, -1), (0, 1)))
            (0, 0)

            sage: H.geometry.segment_intersection(((-1, 0), (1, 0)), ((1, -1), (1, 1)))
            (1, 0)

            sage: H.geometry.segment_intersection(((-1, 0), (1, 0)), ((2, -1), (2, 1))) is None
            True

            sage: H.geometry.segment_intersection(((-1, 0), (1, 0)), ((-2, 0), (0, 0)))
            (-1, 0, 0, 0)

        """
        xy = self._segment_intersection(f, g)

        if xy is None:
            return None

        if len(xy) == 2:
            x, y = xy
            return x[0] / x[1], y[0] / y[1]

        if len(xy) == 4:
            return xy

        assert "unexpected output from _segment_intersection"

    def _line_intersection(self, f, g):
        (fa, fb, fc) = f
        (ga, gb, gc) = g
        det = self._determinant(fb, fc, gb, gc)

        if det is None:
            if self._determinant(fa, fb, ga, gb) is None and self._determinant(fa, fc, ga, gc) is None:
                # The (unoriented) lines are identical
                return f

            return None

        x = (-gc * fa + fc * ga), det
        y = (gb * fa - fb * ga), det

        return (x, y)

    def line_intersection(self, f, g):
        r"""
        Return the intersection between the Euclidean lines ``f`` and ``g``.

        INPUT:

        - ``f`` -- a triple of elements ``(a, b, c)`` of :meth:`base_ring`
          encoding the line `a + bx + cy = 0`

        - ``g`` -- a triple of elements ``(a, b, c)`` of :meth:`base_ring`
          encoding the line `a + bx + cy = 0`

        OUTPUT: ``None`` if they lines do not intersect, ``f`` if the lines
        overlap, or a pair of elements of :meth:`base_ring`, the coordinates of
        the point of intersection.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()

            sage: H.geometry.line_intersection((0, 1, 0), (0, 0, 1))
            (0, 0)

            sage: H.geometry.line_intersection((0, 0, 1), (0, 0, 1))
            (0, 0, 1)

            sage: H.geometry.line_intersection((0, 0, 1), (1, 0, 1)) is None
            True

        """
        xy = self._line_intersection(f, g)
        if xy is None:
            return None

        if len(xy) == 2:
            x, y = xy
            return x[0] / x[1], y[0] / y[1]

        return xy


class ExactGeometry(Geometry):
    r"""
    Shared base class for predicates and geometric constructions over exact rings.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import EuclideanExactGeometry
        sage: geometry = EuclideanExactGeometry(QQ)

    TESTS::

        sage: from flatsurf.geometry.geometry import ExactGeometry
        sage: isinstance(geometry, ExactGeometry)
        True

    ::

        sage: from flatsurf.geometry.hyperbolic import HyperbolicExactGeometry
        sage: isinstance(geometry, ExactGeometry)
        True

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

    def __repr__(self):
        r"""
        Return a printable representation of this geometry.

        EXAMPLES::

            sage: from flatsurf import EuclideanPlane
            sage: E = EuclideanPlane()
            sage: E.geometry
            Exact geometry over Rational Field

        ::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane()
            sage: H.geometry
            Exact geometry over Rational Field

        """
        return f"Exact geometry over {self._ring}"


class EpsilonGeometry(Geometry):
    r"""
    Shared base class for predicates and primitive geometric constructions over
    an inexact base ring.

    EXAMPLES::

        sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
        sage: geometry = EuclideanEpsilonGeometry(RR, 1e-6)

    TESTS::

        sage: from flatsurf.geometry.geometry import EpsilonGeometry
        sage: isinstance(geometry, EpsilonGeometry)
        True

    """

    def __init__(self, ring, epsilon):
        r"""
        TESTS::

            sage: from flatsurf import EuclideanPlane
            sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
            sage: E = EuclideanPlane(RR, EuclideanEpsilonGeometry(RR, 1e-6))
            sage: isinstance(E.geometry, EuclideanEpsilonGeometry)
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

    def __repr__(self):
        r"""
        Return a printable representation of this geometry.

        EXAMPLES::

            sage: from flatsurf.geometry.euclidean import EuclideanEpsilonGeometry
            sage: EuclideanEpsilonGeometry(RR, 1e-6)
            Epsilon geometry with ϵ=1.00000000000000e-6 over Real Field with 53 bits of precision

        """
        return f"Epsilon geometry with ϵ={self._epsilon} over {self._ring}"


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


