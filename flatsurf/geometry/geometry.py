# TODO: Document module.


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
