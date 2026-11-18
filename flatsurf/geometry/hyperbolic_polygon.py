r"""
Compact polygons in the hyperbolic plane, given by their vertices.

A :class:`HyperbolicPolygon` is a possibly non-convex polygon of the
:class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane` whose vertices are
finite points. Its vertices are stored in the Klein model, where geodesics are
straight lines: the polygon is a Euclidean polygon inside the disk of the
quadratic form of the plane, so that its combinatorics (orientation,
convexity, simplicity, position of a point) are those of that Euclidean
polygon. Its metric quantities (side lengths and angles) are decided exactly
in the base ring through that quadratic form, the Minkowski bilinear form of
the hyperboloid model for the unit disk.

EXAMPLES::

    sage: from flatsurf import HyperbolicPlane
    sage: H = HyperbolicPlane(QQ)
    sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
    sage: P
    Hyperbolic polygon with vertices [I, (1/2, 0), (1/2, 1/2), (0, 1/2)]
    sage: P.angles(numerical=True)  # abs tol 1e-12
    (0.25, 0.25, 0.195913276015304, 0.25)
    sage: P.area(numerical=True)  # abs tol 1e-12
    0.0540867239846964

.. SEEALSO::

    :class:`~flatsurf.geometry.hyperbolic.HyperbolicConvexPolygon` for the
    (possibly unbounded) convex polygons obtained as intersections of half
    spaces

"""

######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Vincent Delecroix
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

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent


def _minkowski(u, v, form=None):
    r"""
    Return the Minkowski product `u_0 v_0 + u_1 v_1 - u_2 v_2` of two vectors
    of dimension three, or `u^T Q v` for the diagonal matrix ``form`` `Q`.

    The form is not validated: use
    :func:`~flatsurf.geometry.hyperbolic_isometry._minkowski_form` to build it.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_polygon import _minkowski
        sage: _minkowski((1, 2, 3), (1, 1, 1))
        0
        sage: _minkowski((1, 1, 1), (1, 1, 1), diagonal_matrix([1, 2, -3]))
        0

    """
    if form is None:
        return u[0] * v[0] + u[1] * v[1] - u[2] * v[2]

    return (
        form[0, 0] * u[0] * v[0] + form[1, 1] * u[1] * v[1] + form[2, 2] * u[2] * v[2]
    )


def _form_repr(form):
    r"""
    Return the printed form `x^2 + b y^2 - c z^2` of the diagonal matrix
    ``form`` `\mathrm{diag}(1, b, -c)`.

    A coefficient is put in parentheses when it prints with a sum, a
    difference or a leading minus sign, so that the signs of the form stay
    readable.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_polygon import _form_repr
        sage: _form_repr(diagonal_matrix(QQ, [1, 2, -3]))
        'x^2 + 2*y^2 - 3*z^2'
        sage: _form_repr(diagonal_matrix(QQ, [1, 1, -1]))
        'x^2 + y^2 - z^2'

        sage: K.<a> = NumberField(x^2 - 2, embedding=-1.414)
        sage: _form_repr(diagonal_matrix(K, [1, -a, 2*a]))
        'x^2 + (-a)*y^2 - (-2*a)*z^2'

    """

    def monomial(coefficient, variable):
        if coefficient == 1:
            return variable
        coefficient = repr(coefficient)
        if "+" in coefficient or " - " in coefficient or coefficient.startswith("-"):
            coefficient = f"({coefficient})"
        return f"{coefficient}*{variable}"

    return f"x^2 + {monomial(form[1, 1], 'y^2')} - {monomial(-form[2, 2], 'z^2')}"


def _dual_minkowski(u, v, form=None):
    r"""
    Return the product `u^T Q^{-1} v` of the dual form of the diagonal
    matrix ``form`` `Q` (default: `\mathrm{diag}(1, 1, -1)`) on two triples of
    coefficients.

    The line `a + bx + cy = 0` of the Klein model is the triple `(b, c, a)`,
    so that it vanishes on the lift `(x, y, 1)` of its points. The dual
    product of a geodesic with itself is positive, zero or negative when
    the line crosses, touches or misses the disk of the form; two geodesics
    are perpendicular when their dual product vanishes.

    The form is not validated: use
    :func:`~flatsurf.geometry.hyperbolic_isometry._minkowski_form` to build it.
    Without a form, no division takes place.

    EXAMPLES:

    For `Q = \mathrm{diag}(1, 2, -3)`, the line `x = 1/3` crosses the disk,
    the line `-3 + x + 2y = 0` touches it, and the line `x = 2` misses it::

        sage: from flatsurf.geometry.hyperbolic_polygon import _dual_minkowski
        sage: Q = diagonal_matrix(QQ, [1, 2, -3])
        sage: _dual_minkowski((1, 0, -1/3), (1, 0, -1/3), Q)
        26/27
        sage: _dual_minkowski((1, 2, -3), (1, 2, -3), Q)
        0
        sage: _dual_minkowski((1, 0, -2), (1, 0, -2), Q)
        -1/3
        sage: _dual_minkowski((0, 1, -1), (0, 1, -1), Q)
        1/6
        sage: _dual_minkowski((1, -1, 0), (1, -1, 0), Q)
        3/2
        sage: _dual_minkowski((1, 2, -3/2), (1, 2, -3/2), Q)
        9/4

    The lines `x = 1/3` and `y = 0` are perpendicular, and so are the lines
    `x - y = 0` and `x + 2y = 3/2`::

        sage: _dual_minkowski((1, 0, -1/3), (0, 1, 0), Q)
        0
        sage: _dual_minkowski((1, -1, 0), (1, 2, -3/2), Q)
        0

    Without a form, this is the product of the standard plane::

        sage: _dual_minkowski((1, 2, 3), (1, 1, 1))
        0

    """
    if form is None:
        return u[0] * v[0] + u[1] * v[1] - u[2] * v[2]

    return (
        u[0] * v[0] / form[0, 0] + u[1] * v[1] / form[1, 1] + u[2] * v[2] / form[2, 2]
    )


def _lift(v):
    r"""
    Return the lift `(x, y, 1)` of the point ``v`` of the Klein model to the
    cone of timelike vectors of Minkowski space.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_polygon import _lift
        sage: _lift((1/2, 0))
        (1/2, 0, 1)

    """
    return (v[0], v[1], v[0].parent().one())


def _cosh2(p, q, form=None):
    r"""
    Return the square of the hyperbolic cosine of the distance between the
    points ``p`` and ``q`` of the Klein model, for the diagonal matrix
    ``form`` (default: `\mathrm{diag}(1, 1, -1)`).

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_polygon import _cosh2
        sage: V = QQ**2
        sage: _cosh2(V((0, 0)), V((1/2, 0)))
        4/3
        sage: Q = diagonal_matrix([1, 2, -3])
        sage: _cosh2(V((0, 0)), V((1, 0)), Q)
        3/2
        sage: _cosh2(V((1, 0)), V((0, 1)), Q)
        9/2

    """
    P = _lift(p)
    R = _lift(q)
    return _minkowski(P, R, form) ** 2 / (
        _minkowski(P, P, form) * _minkowski(R, R, form)
    )


def _rational_arccos(sign, cos2):
    r"""
    Return the angle in `(0, 1/2)`, as a multiple of `2\pi`, whose cosine has
    the sign ``sign`` and the square ``cos2``, if it is rational; return
    ``None`` otherwise.

    INPUT:

    - ``sign`` -- one of `-1`, `0`, `1`

    - ``cos2`` -- an element of a subfield of the real algebraic numbers in
      `[0, 1)`

    ALGORITHM:

    If the angle is `k/m` with `k` and `m` coprime and `m \ge 3`, then its
    cosine has degree `\varphi(m)/2` over `\QQ`. That degree is at most
    twice the degree `d` of ``cos2``, so `\varphi(m) \le 4d`, and
    `\varphi(m) \ge \sqrt{m/2}` gives `m \le 32 d^2 = B`. Two distinct
    rationals with denominators at most `B` are at least `1/B^2` apart, so
    an approximation of the angle within `1/(2B^2)` of `k/m` has `k/m` as
    its best rational approximation with denominator at most `B`. That
    rational is the only candidate; it is then checked exactly against the
    cosine of `2\pi k/m` in the real algebraic numbers, so that the
    precision only matters for not missing a rational angle.

    The cosine is computed with `p = 64 + 4 b` bits, `b` the number of bits
    of `B`, so that its error is `\varepsilon \le 4 \cdot 2^{-p}`. The
    arccosine is ill-conditioned near `\pm 1`, but `|\arccos u - \arccos v|
    \le \pi \sqrt{|u - v|/2}` (the worst case is at an endpoint, where
    `\arccos(1 - \delta) = 2 \arcsin\sqrt{\delta/2}`). The error on the angle
    as a multiple of `2\pi` is thus at most `\sqrt{\varepsilon/8} \le
    2^{-(p + 1)/2} \le 2^{-32} B^{-2}`, well within `1/(2B^2)` even with
    the rounding of the arccosine and of the division.

    EXAMPLES::

        sage: from flatsurf.geometry.hyperbolic_polygon import _rational_arccos
        sage: _rational_arccos(1, 1/4)
        1/6
        sage: _rational_arccos(-1, 1/2)
        3/8
        sage: _rational_arccos(0, 0)
        1/4
        sage: _rational_arccos(1, 1/9) is None
        True

    """
    from sage.all import AA, QQ, QQbar, ZZ, RealField

    if sign == 0:
        return QQ((1, 4))

    cos2 = AA(cos2)
    degree = cos2.minpoly().degree()
    bound = ZZ(32 * degree**2)

    R = RealField(64 + 4 * bound.nbits())
    cos = R(sign) * R(cos2).sqrt()
    approximation = cos.arccos() / (2 * R.pi())
    candidate = QQ(approximation.nearby_rational(max_denominator=bound))

    if not 0 < candidate < QQ((1, 2)):
        return None

    exact = AA((QQbar.zeta(candidate.denominator()) ** candidate.numerator()).real())
    if exact.sign() != sign or exact**2 != cos2:
        return None

    return candidate


def hyperbolic_polygon(plane, vertices, model="klein", check=True):
    r"""
    Return the polygon of the hyperbolic ``plane`` with the given
    ``vertices``.

    This is a helper function for
    :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.polygon`.

    INPUT:

    - ``plane`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane`

    - ``vertices`` -- a sequence of points of ``plane`` or of pairs of
      coordinates in ``model``, in counterclockwise order; each must be a
      finite point

    - ``model`` -- ``"klein"`` or ``"half_plane"`` (default: ``"klein"``);
      the model in which pairs of coordinates are given

    - ``check`` -- a boolean (default: ``True``); whether to verify that the
      vertices are in counterclockwise order and that the polygon does not
      self-intersect

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: from flatsurf.geometry.hyperbolic_polygon import hyperbolic_polygon
        sage: H = HyperbolicPlane(QQ)
        sage: hyperbolic_polygon(H, [I, 1 + I, 2*I], model="half_plane")
        Hyperbolic polygon with vertices [I, 1 + I, 2*I]

    TESTS::

        sage: hyperbolic_polygon(H, [(0, 0), (1, 0), (0, 1/2)])
        Traceback (most recent call last):
        ...
        ValueError: vertex (1, 0) is not a finite point of the hyperbolic plane
        sage: hyperbolic_polygon(H, [(0, 0), (1/2, 0), (0, 1/2)], model="poincare")
        Traceback (most recent call last):
        ...
        NotImplementedError: unsupported model 'poincare'

    The vertices must lie in the disk of the quadratic form of the plane,
    here `x^2 + 2y^2 < 3`::

        sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
        sage: hyperbolic_polygon(H3, [(0, 0), (1, 0), (0, 1/2)])
        Hyperbolic polygon with vertices [(0, 0), (1, 0), (0, 1/2)]
        sage: hyperbolic_polygon(H3, [(0, 0), (2, 0), (0, 1/2)])
        Traceback (most recent call last):
        ...
        ValueError: vertex (2, 0) is not a finite point of the hyperbolic plane

    """
    from flatsurf.geometry.hyperbolic import HyperbolicPoint

    if model not in ("klein", "half_plane"):
        raise NotImplementedError(f"unsupported model {model!r}")

    K = plane.base_ring()

    klein = []
    for vertex in vertices:
        if isinstance(vertex, HyperbolicPoint):
            vertex = plane(vertex)
            if not vertex.is_finite():
                raise ValueError(
                    f"vertex {vertex} is not a finite point of the hyperbolic plane"
                )
            x, y = vertex.coordinates(model="klein")
        else:
            x, y = vertex
            if model == "half_plane":
                point = plane.point(x, y, model="half_plane")
                if not point.is_finite():
                    raise ValueError(
                        f"vertex {point} is not a finite point of the hyperbolic plane"
                    )
                x, y = point.coordinates(model="klein")
        x = K(x)
        y = K(y)
        if not _minkowski(_lift((x, y)), _lift((x, y)), plane._form) < 0:
            raise ValueError(
                f"vertex {(x, y)} is not a finite point of the hyperbolic plane"
            )
        klein.append((x, y))

    polygon = HyperbolicPolygon(plane, klein)

    if check:
        polygon._check()

    return polygon


class HyperbolicPolygon(Parent):
    r"""
    A compact polygon of the hyperbolic plane with finite vertices, possibly
    non-convex.

    The polygon is a facade of its hyperbolic plane: its elements are points
    of the plane. Its vertices are stored in the Klein model, in the order
    given, which is counterclockwise for a valid polygon.

    INPUT:

    - ``parent`` -- a :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane`

    - ``vertices`` -- a sequence of pairs of elements of the base ring, the
      Klein coordinates of the vertices, each strictly inside the disk of the
      quadratic form of ``parent``, the unit disk for `x^2 + y^2 - z^2`

    - ``category`` -- a category or ``None`` (default: ``None``); refined to
      convex and simple polygons when the vertices describe such a polygon

    EXAMPLES::

        sage: from flatsurf import HyperbolicPlane
        sage: H = HyperbolicPlane(QQ)
        sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
        sage: P
        Hyperbolic polygon with vertices [(0, -1/2), (1/2, 1/2), I, (-1/2, 1/2)]
        sage: P.category()
        Category of facade simple hyperbolic polygons over Rational Field

    TESTS::

        sage: from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon
        sage: isinstance(P, HyperbolicPolygon)
        True
        sage: TestSuite(P).run()

        sage: Q = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
        sage: Q.category()
        Category of facade convex simple hyperbolic polygons over Rational Field
        sage: TestSuite(Q).run()

    """

    def __init__(self, parent, vertices, category=None):
        from flatsurf.geometry.categories import HyperbolicPolygons
        from flatsurf.geometry.euclidean import is_convex_polygon, is_simple_polygon

        K = parent.base_ring()
        V = K**2
        self._v = tuple(V(vertex) for vertex in vertices)
        for v in self._v:
            v.set_immutable()

        if category is None:
            category = HyperbolicPolygons(K)

        category &= HyperbolicPolygons(K)

        if "Simple" not in category.axioms() and is_simple_polygon(self._v):
            category = category.Simple()

        if "Convex" not in category.axioms() and is_convex_polygon(self._v):
            category = category.Convex()

        Parent.__init__(self, base=K, facade=parent, category=category)

    def _check(self):
        r"""
        Verify that this polygon is valid: at least three vertices, in
        counterclockwise order, without zero or backtracking edges and
        without self-intersection.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (0, 1/2), (1/2, 0)])
            Traceback (most recent call last):
            ...
            ValueError: polygon has negative area; probably the vertices are not in counter-clockwise order
            sage: H.polygon(vertices=[(0, 0), (1/2, 0)])
            Traceback (most recent call last):
            ...
            ValueError: polygon must have at least three sides
            sage: H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0)])
            Traceback (most recent call last):
            ...
            ValueError: polygon has zero area
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 0), (0, 1/2)])
            Traceback (most recent call last):
            ...
            ValueError: polygon has zero edge
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/4, 0), (0, 1/2)])
            Traceback (most recent call last):
            ...
            ValueError: polygon has anti-parallel edges
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/4), (1/4, -1/4), (1/2, 1/4)])
            Traceback (most recent call last):
            ...
            ValueError: polygon self-intersects

        """
        from flatsurf.geometry.euclidean import is_anti_parallel, is_simple_polygon

        n = len(self._v)
        if n < 3:
            raise ValueError("polygon must have at least three sides")

        area = self._klein_area()
        if area < 0:
            raise ValueError(
                "polygon has negative area; probably the vertices are not in counter-clockwise order"
            )
        if area == 0:
            raise ValueError("polygon has zero area")

        edges = [self._v[(i + 1) % n] - self._v[i] for i in range(n)]
        if any(edge == 0 for edge in edges):
            raise ValueError("polygon has zero edge")

        if any(is_anti_parallel(edges[i - 1], edges[i]) for i in range(n)):
            raise ValueError("polygon has anti-parallel edges")

        if not is_simple_polygon(self._v):
            raise ValueError("polygon self-intersects")

    def _klein_area(self):
        r"""
        Return the signed Euclidean area of this polygon in the Klein model.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P._klein_area()
            1/4

        """
        n = len(self._v)
        total = sum(
            self._v[i][0] * self._v[(i + 1) % n][1]
            - self._v[(i + 1) % n][0] * self._v[i][1]
            for i in range(n)
        )
        return total / 2

    def parent(self):
        r"""
        Return the hyperbolic plane this polygon is a subset of.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P.parent()
            Hyperbolic Plane over Rational Field

        """
        return self.facade_for()[0]

    def _element_constructor_(self, x):
        r"""
        Return ``x`` as a point of this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P(I)
            I
            sage: P(2*I)
            Traceback (most recent call last):
            ...
            ValueError: point not contained in this polygon

        """
        x = self.parent()(x)

        if not self.get_point_position(x).is_inside():
            raise ValueError("point not contained in this polygon")

        return x

    def _an_element_(self):
        r"""
        Return a point of this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P.an_element()
            I

        """
        return self.vertex(0)

    def some_elements(self):
        r"""
        Return some points of this polygon, its vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P.some_elements()
            (I, (1/2, 0), (0, 1/2))

        """
        return self.vertices()

    def __eq__(self, other):
        r"""
        Return whether this polygon lies in the same plane and has the same
        vertices, in the same order, as ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P == H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            True
            sage: P == H.polygon(vertices=[(1/2, 0), (0, 1/2), (0, 0)])
            False
            sage: P != H.polygon(vertices=[(1/2, 0), (0, 1/2), (0, 0)])
            True
            sage: P == P.change_ring(AA)
            False

        Polygons with the same Klein coordinates in planes of different forms
        differ::

            sage: P == HyperbolicPlane(QQ, form=[1, 2, -3]).polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            False

        """
        if not isinstance(other, HyperbolicPolygon):
            return False

        return self.parent() == other.parent() and self._v == other._v

    def __ne__(self, other):
        r"""
        Return whether this polygon differs from ``other``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P != H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            False

        """
        return not self == other

    def __hash__(self):
        r"""
        Return a hash value for this polygon, compatible with ``==``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: hash(P) == hash(H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            True

        """
        return hash((self.parent(), self._v))

    def _repr_(self):
        r"""
        Return a printable representation of this polygon.

        The vertices print as points of the hyperbolic plane: in the upper
        half plane model when their coordinates there are in the base ring,
        in the Klein model otherwise.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[I, 1 + I, 2*I], model="half_plane")
            Hyperbolic polygon with vertices [I, 1 + I, 2*I]

        """
        return f"Hyperbolic polygon with vertices {list(self.vertices())!r}"

    def vertices(self, marked_vertices=True):
        r"""
        Return the vertices of this polygon, in counterclockwise order
        starting from the start vertex.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``); whether to
          include the vertices with an angle π

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0), (0, 1/2)])
            sage: P.vertices()
            (I, (1/4, 0), (1/2, 0), (0, 1/2))
            sage: P.vertices(marked_vertices=False)
            (I, (1/2, 0), (0, 1/2))

        """
        vertices = self._v
        if not marked_vertices:
            from flatsurf.geometry.euclidean import polygon_corners

            vertices = polygon_corners(vertices)

        H = self.parent()
        return tuple(H.point(x, y, model="klein", check=False) for (x, y) in vertices)

    def vertex(self, i):
        r"""
        Return the ``i``-th vertex of this polygon, with ``i`` taken modulo
        the number of vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P.vertex(1)
            (1/2, 0)
            sage: P.vertex(-1)
            (0, 1/2)
            sage: P.vertex(3)
            I

        """
        x, y = self._v[i % len(self._v)]
        return self.parent().point(x, y, model="klein", check=False)

    def edge_segment(self, i):
        r"""
        Return the oriented geodesic segment from the vertex ``i`` to the
        vertex ``i + 1``, with ``i`` taken modulo the number of vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[I, 1 + I, 2*I], model="half_plane")
            sage: s = P.edge_segment(2)
            sage: s
            {x = 0} ∩ {(x^2 + y^2) - 4 ≤ 0} ∩ {(x^2 + y^2) - 1 ≥ 0}
            sage: s.start(), s.end()
            (2*I, I)
            sage: P.edge_segment(-1) == s
            True
            sage: P.edge_segment(3).start()
            I

        """
        return self.vertex(i).segment(self.vertex(i + 1))

    def edge_segments(self):
        r"""
        Return the oriented geodesic segments joining consecutive vertices of
        this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[I, 1 + I, 2*I], model="half_plane")
            sage: [(s.start(), s.end()) for s in P.edge_segments()]
            [(I, 1 + I), (1 + I, 2*I), (2*I, I)]

        """
        return [self.edge_segment(i) for i in range(len(self._v))]

    def is_convex(self, strict=False):
        r"""
        Return whether this polygon is convex.

        INPUT:

        - ``strict`` -- a boolean (default: ``False``); whether to check for
          strict convexity, i.e., a polygon with a π angle is not considered
          convex

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0), (0, 1/2)])
            sage: P.is_convex()
            True
            sage: P.is_convex(strict=True)
            False
            sage: H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)]).is_convex()
            False

        """
        from flatsurf.geometry.euclidean import is_convex_polygon

        return is_convex_polygon(self._v, strict=strict)

    def is_simple(self):
        r"""
        Return whether this polygon does not self-intersect.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)]).is_simple()
            True
            sage: H.polygon(vertices=[(0, 0), (1/2, 1/2), (1/2, 0), (0, 1/2)], check=False).is_simple()
            False

        """
        from flatsurf.geometry.euclidean import is_simple_polygon

        return is_simple_polygon(self._v)

    def is_degenerate(self):
        r"""
        Return whether this polygon is degenerate, i.e., has zero area or
        marked vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]).is_degenerate()
            False
            sage: H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0), (0, 1/2)]).is_degenerate()
            True
            sage: H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0)], check=False).is_degenerate()
            True

        """
        from flatsurf.geometry.euclidean import polygon_corners

        if self._klein_area() == 0:
            return True

        return len(polygon_corners(self._v)) != len(self._v)

    def get_point_position(self, point):
        r"""
        Return the combinatorial position of ``point`` with respect to this
        polygon.

        INPUT:

        - ``point`` -- a point of the hyperbolic plane, or a pair of
          coordinates in the Klein model

        OUTPUT:

        A :class:`~flatsurf.geometry.polygon.PolygonPosition`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P.get_point_position((1/4, 1/4))
            point positioned in interior of polygon
            sage: P.get_point_position(I)
            point positioned on vertex 0 of polygon
            sage: P.get_point_position(H.point(1/2, 1/4, model="klein"))
            point positioned on interior of edge 1 of polygon
            sage: P.get_point_position((3/4, 0))
            point positioned outside polygon
            sage: P.get_point_position(0)
            point positioned outside polygon

        """
        from flatsurf.geometry.euclidean import polygon_point_position

        return polygon_point_position(self._v, self._klein_vector(point))

    def _klein_vector(self, point):
        r"""
        Return the Klein coordinates of ``point`` as a vector.

        INPUT:

        - ``point`` -- a point of the hyperbolic plane, a pair of coordinates
          in the Klein model, or anything the hyperbolic plane converts to a
          point

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P._klein_vector((1/4, 0))
            (1/4, 0)
            sage: P._klein_vector(I)
            (0, 0)
            sage: P._klein_vector(H.vertical(0))
            Traceback (most recent call last):
            ...
            TypeError: not a point of the hyperbolic plane

        """
        from flatsurf.geometry.hyperbolic import HyperbolicPoint

        V = self.base_ring() ** 2

        if not isinstance(point, HyperbolicPoint):
            try:
                pair = len(point) == 2
            except TypeError:
                pair = False

            if pair:
                return V(point)

            point = self.parent()(point)

            if not isinstance(point, HyperbolicPoint):
                raise TypeError("not a point of the hyperbolic plane")

        return V(self.parent()(point).coordinates(model="klein"))

    def __contains__(self, point):
        r"""
        Return whether ``point`` is in this (closed) polygon.

        INPUT:

        - ``point`` -- a point of the hyperbolic plane, or a pair of
          coordinates in the Klein model

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: I in P
            True
            sage: (1/8, 1/8) in P
            True
            sage: (1/2, 1/2) in P
            False
            sage: H.vertical(0) in P
            False

        """
        try:
            return self.get_point_position(point).is_inside()
        except TypeError:
            return False

    def change_ring(self, ring):
        r"""
        Return this polygon as a polygon of the hyperbolic plane over
        ``ring``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: Q = P.change_ring(AA)
            sage: Q.base_ring()
            Algebraic Real Field
            sage: Q.parent()
            Hyperbolic Plane over Algebraic Real Field

        The quadratic form of the plane is kept::

            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: W = H3.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: W.change_ring(AA).parent()
            Hyperbolic Plane over Algebraic Real Field with form x^2 + 2*y^2 - 3*z^2

        """
        if ring is self.base_ring():
            return self

        return HyperbolicPolygon(self.parent().change_ring(ring), self._v)

    def erase_marked_vertices(self):
        r"""
        Return a copy of this polygon without its vertices of angle π.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/4, 0), (1/2, 0), (0, 1/2)])
            sage: P.erase_marked_vertices()
            Hyperbolic polygon with vertices [I, (1/2, 0), (0, 1/2)]

        """
        from flatsurf.geometry.euclidean import polygon_corners

        return HyperbolicPolygon(self.parent(), polygon_corners(self._v))

    def cyclic_quotient(self, m):
        r"""
        Return a fundamental domain of the rotation of order ``m`` that
        leaves this polygon invariant.

        Write `v_0, \ldots, v_{n-1}` for the vertices of this polygon, marked
        vertices included, and `k = n/m`. The rotation `g` is the orientation
        preserving isometry that sends the edge `v_0 v_1` onto the edge `v_k
        v_{k+1}`; it must send every `v_i` to `v_{i+k}` (indices modulo `n`)
        and is then of order `m`. Its centre `c` is its fixed point.

        For `m \ge 3`, the result is the sector `(c, v_0, v_1, \ldots, v_k)`:
        gluing its edge `0` (from `c` to `v_0`) to its edge `k + 1` (from `v_k`
        to `c`) gives the quotient of this polygon by `g`, with a cone point
        of angle `1/m` (in multiples of `2\pi`) at `c`. For `m = 2`, the result
        is the half `(v_0, v_1, \ldots, v_k)`: gluing its last edge (from
        `v_k` to `v_0`, through `c`) to itself gives the quotient, the
        half-turn `g` exchanging the two halves of that edge about its
        midpoint `c`.

        The rotation permutes the vertices, so its matrix in
        `\mathrm{SO}(2, 1)^+`, or in the group of the quadratic form of the
        plane, and its centre `c` are defined over the base ring: so is the
        result.

        INPUT:

        - ``m`` -- an integer at least 2 that divides the number of vertices

        OUTPUT:

        A polygon over the base ring of this polygon, of the same hyperbolic
        plane. A :class:`ValueError` is raised when ``m`` is not an integer at
        least 2, when `m \nmid n`, or when no rotation sends every `v_i` to
        `v_{i+k}`. For a non-convex polygon, the sector may be no valid
        polygon, in which case the polygon constructor raises its error.

        EXAMPLES:

        A square cell centred at the Klein point `(3/5, 0)`, i.e., `3/5 +
        4/5 i` in the upper half plane::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(11/13, 0), (3/5, 2/5), (1/7, 0), (3/5, -2/5)])
            sage: Q = P.cyclic_quotient(4)
            sage: Q
            Hyperbolic polygon with vertices [3/5 + 4/5*I, (11/13, 0), (3/5, 2/5)]

        Its edges 0 and 2 glued, the centre is a cone point of angle `1/4`::

            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: S = MutableOrientedHyperbolicSurface(QQ)
            sage: S.add_polygon(Q)
            0
            sage: S.glue((0, 0), (0, 2))
            sage: S.set_immutable()
            sage: S(0, 0).angle()
            1/4
            sage: S(0, 0).angle(numerical=True)
            0.25

        The quotient by the half-turn is half of the square, its last edge
        glued to itself::

            sage: P.cyclic_quotient(2)
            Hyperbolic polygon with vertices [(11/13, 0), (3/5, 2/5), (1/7, 0)]

        TESTS::

            sage: P.cyclic_quotient(3)
            Traceback (most recent call last):
            ...
            ValueError: polygon is not invariant under a rotation of order 3
            sage: P.cyclic_quotient(1)
            Traceback (most recent call last):
            ...
            ValueError: m must be an integer at least 2
            sage: P.cyclic_quotient(1/2)
            Traceback (most recent call last):
            ...
            ValueError: m must be an integer at least 2

        A rectangle has a half-turn but no quarter turn::

            sage: R = H.polygon(vertices=[(-1/2, -1/4), (1/2, -1/4), (1/2, 1/4), (-1/2, 1/4)])
            sage: R.cyclic_quotient(2)
            Hyperbolic polygon with vertices [(-1/2, -1/4), (1/2, -1/4), (1/2, 1/4)]
            sage: R.cyclic_quotient(4)
            Traceback (most recent call last):
            ...
            ValueError: polygon is not invariant under a rotation of order 4

        A Lambert quadrilateral has no rotation::

            sage: L = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: L.cyclic_quotient(2)
            Traceback (most recent call last):
            ...
            ValueError: polygon is not invariant under a rotation of order 2

        The rotation is an isometry of the quadratic form of the plane: this
        triangle is invariant under a rotation of order 3 for the form `x^2 +
        3y^2 - z^2`, not for `x^2 + y^2 - z^2`::

            sage: H31 = HyperbolicPlane(QQ, form=[1, 3, -1])
            sage: T = H31.polygon(vertices=[(1/2, 0), (-1/4, 1/4), (-1/4, -1/4)])
            sage: T._cosh2_sides()
            (9/4, 9/4, 9/4)
            sage: T.cyclic_quotient(3)
            Hyperbolic polygon with vertices [(0, 0), (1/2, 0), (-1/4, 1/4)]
            sage: H.polygon(vertices=[(1/2, 0), (-1/4, 1/4), (-1/4, -1/4)]).cyclic_quotient(3)
            Traceback (most recent call last):
            ...
            ValueError: polygon is not invariant under a rotation of order 3

        """
        from sage.all import ZZ

        from flatsurf.geometry.hyperbolic_isometry import (
            HyperbolicIsometryGroup,
            segment_isometry,
        )

        try:
            m = ZZ(m)
        except (TypeError, ValueError) as error:
            raise ValueError("m must be an integer at least 2") from error
        if m < 2:
            raise ValueError("m must be an integer at least 2")

        refused = f"polygon is not invariant under a rotation of order {m}"

        vertices = self._v
        n = len(vertices)
        if n % m:
            raise ValueError(refused)
        k = n // m

        K = self.base_ring()
        form = self.parent()._form
        try:
            M = segment_isometry(
                (vertices[0], vertices[1]),
                (vertices[k], vertices[(k + 1) % n]),
                ring=K,
                form=form,
            )
        except ValueError as error:
            raise ValueError(refused) from error

        for i, (x, y) in enumerate(vertices):
            X = M[0, 0] * x + M[0, 1] * y + M[0, 2]
            Y = M[1, 0] * x + M[1, 1] * y + M[1, 2]
            Z = M[2, 0] * x + M[2, 1] * y + M[2, 2]
            if (X / Z, Y / Z) != tuple(vertices[(i + k) % n]):
                raise ValueError(refused)

        if m == 2:
            return self.parent().polygon(vertices=vertices[: k + 1])

        centre = HyperbolicIsometryGroup(K, form=form)(M).fixed_point()
        return self.parent().polygon(
            vertices=[centre.coordinates(model="klein")] + list(vertices[: k + 1])
        )

    def _transformation_group(self):
        r"""
        Return the group of the orientation preserving isometries of the
        hyperbolic plane, over the base ring of this polygon, with the
        quadratic form of its plane.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])._transformation_group()
            Group of orientation preserving isometries of the hyperbolic plane over Rational Field
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])._transformation_group()
            Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2

        """
        from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup

        return HyperbolicIsometryGroup(
            self.parent().base_ring(), form=self.parent()._form
        )

    def _mul_(self, g, switch_sides=None):
        r"""
        Return the image of this polygon under the hyperbolic isometry ``g``.

        A matrix is refused since it is ambiguous between the models; use
        :meth:`apply_isometry` for it.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: M = matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            sage: h = HyperbolicIsometryGroup(QQ)(M)
            sage: h * P
            Hyperbolic polygon with vertices [I, (0, 1/2), (-1/2, 0)]
            sage: P * h
            Hyperbolic polygon with vertices [I, (0, 1/2), (-1/2, 0)]

            sage: K.<a> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
            sage: hK = HyperbolicIsometryGroup(K)(matrix(K, [[a/2, -a/2, 0], [a/2, a/2, 0], [0, 0, 1]]))
            sage: hK * P
            Hyperbolic polygon with vertices [I, (1/4*a, 1/4*a), (-1/4*a, 1/4*a)]

            sage: M * P
            Traceback (most recent call last):
            ...
            TypeError: cannot multiply a hyperbolic polygon with an object of type <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>; use apply_isometry() to apply a matrix

        """
        from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometry

        if isinstance(g, HyperbolicIsometry):
            return g(self)

        raise TypeError(
            "cannot multiply a hyperbolic polygon with an object of type {}; use apply_isometry() to apply a matrix".format(
                type(g)
            )
        )

    def apply_isometry(self, isometry, model="klein"):
        r"""
        Return the image of this polygon under the orientation preserving
        ``isometry``.

        INPUT:

        - ``isometry`` -- a
          :class:`~flatsurf.geometry.hyperbolic_isometry.HyperbolicIsometry`,
          or a matrix as in
          :meth:`~flatsurf.geometry.hyperbolic.HyperbolicConvexSet.apply_isometry`:
          a 2×2 matrix of positive determinant if ``model`` is
          ``"half_plane"``, a 3×3 matrix preserving the quadratic form of the
          plane up to a positive factor if ``model`` is ``"klein"``

        - ``model`` -- ``"klein"`` or ``"half_plane"`` (default:
          ``"klein"``); the model in which a matrix ``isometry`` acts, ignored
          for a :class:`~flatsurf.geometry.hyperbolic_isometry.HyperbolicIsometry`

        OUTPUT:

        The polygon of this hyperbolic plane whose vertices are the images of
        the vertices of this polygon, in the same order.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P.apply_isometry(matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
            Hyperbolic polygon with vertices [I, (0, 1/2), (-1/2, 0)]
            sage: P.apply_isometry(matrix([[2, 0], [0, 1]]), model="half_plane")
            Hyperbolic polygon with vertices [2*I, (2/5, 3/5), (0, 11/13)]

            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: g = HyperbolicIsometryGroup(QQ)(matrix([[2, 0], [0, 1]]))
            sage: P.apply_isometry(g) == g(P) == P.apply_isometry(matrix([[2, 0], [0, 1]]), model="half_plane")
            True

        The isometry must preserve the orientation::

            sage: P.apply_isometry(matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: isometry must preserve the orientation
            sage: P.apply_isometry(matrix([[1, 0], [0, -1]]), model="half_plane")
            Traceback (most recent call last):
            ...
            ValueError: isometry must preserve the orientation
            sage: P.apply_isometry(matrix([[1, 2, 0], [0, 1, 0], [0, 0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry
            sage: P.apply_isometry(matrix([[1, 0], [0, 0]]), model="half_plane")
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry
            sage: P.apply_isometry(matrix([[1, 0], [0, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: isometry in Klein model must be given as a 3×3 matrix

        On a plane with another form, the isometry preserves that form::

            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: W = H3.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: B = matrix([[2, 0, 3], [0, 1, 0], [1, 0, 2]])
            sage: W.apply_isometry(B)
            Hyperbolic polygon with vertices [(3/2, 0), (8/5, 0), (8/5, 1/5), (3/2, 1/4)]
            sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
            sage: W.apply_isometry(B) == G3(B)(W) == G3(B) * W
            True
            sage: P.apply_isometry(G3(B))
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry
            sage: W.apply_isometry(g)
            Traceback (most recent call last):
            ...
            ValueError: matrix does not define an isometry

        .. SEEALSO::

            :meth:`~flatsurf.geometry.hyperbolic.HyperbolicPlane.isomorphism`
            to map a polygon to a hyperbolic plane of another form

        """
        from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometry

        H = self.parent()
        K = H.base_ring()

        if isinstance(isometry, HyperbolicIsometry):
            M = isometry.matrix()
        elif model == "half_plane":
            if isometry.dimensions() != (2, 2):
                raise ValueError(
                    "isometry in half plane model must be given as a 2×2 matrix"
                )
            if isometry.determinant() == 0:
                raise ValueError("matrix does not define an isometry")
            M = H._isometry_gl2_to_sim12(isometry)
        elif model == "klein":
            if isometry.dimensions() != (3, 3):
                raise ValueError(
                    "isometry in Klein model must be given as a 3×3 matrix"
                )
            M = isometry
        else:
            raise NotImplementedError(f"unsupported model {model!r}")

        M = M.change_ring(K)
        J = H.quadratic_form()
        D = M.transpose() * J * M
        if not D[0, 0] > 0 or D != D[0, 0] * J:
            raise ValueError("matrix does not define an isometry")

        images = []
        for x, y in self._v:
            X = M[0, 0] * x + M[0, 1] * y + M[0, 2]
            Y = M[1, 0] * x + M[1, 1] * y + M[1, 2]
            Z = M[2, 0] * x + M[2, 1] * y + M[2, 2]
            images.append((X / Z, Y / Z))

        polygon = HyperbolicPolygon(H, images)

        # An isometry is a diffeomorphism of the Klein disk, so it reverses
        # the orientation of every polygon or of none.
        if polygon._klein_area() < 0:
            raise ValueError("isometry must preserve the orientation")

        return polygon

    def edge_cosh_squared(self, i):
        r"""
        Return the square of the hyperbolic cosine of the length of the edge
        ``i`` of this polygon.

        The square lies in the base ring of this polygon, which the
        hyperbolic cosine itself does not in general.

        INPUT:

        - ``i`` -- an integer, the index of the edge from the vertex ``i`` to
          the vertex ``i + 1``, with `0 \le i <` the number of vertices

        OUTPUT:

        An element of the base ring.

        EXAMPLES:

        The edges of the quadrilateral whose vertices are those of a
        Euclidean square in the Klein model::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: [P.edge_cosh_squared(i) for i in range(4)]
            [4/3, 3/2, 3/2, 4/3]

        The edge of a regular octagon, from the vertex 0 to the vertex 1::

            sage: from flatsurf import hyperbolic_polygons
            sage: P = hyperbolic_polygons.regular_ngon(8, 1/8)
            sage: P.edge_cosh_squared(0)
            40*a^2 + 17
            sage: P.edge_cosh_squared(0).parent() is P.base_ring()
            True

        The edge 1 of a triangle with a right angle at the vertex 0::

            sage: T = hyperbolic_polygons.triangle(QQ(1)/4, QQ(1)/12, QQ(1)/12)
            sage: T.base_ring()
            Number Field in a with defining polynomial y^2 - 6 with a = 2.449489742783178?
            sage: T.edge_cosh_squared(1)
            9

        TESTS::

            sage: T.edge_cosh_squared(3)
            Traceback (most recent call last):
            ...
            ValueError: the edge 3 is not an edge of a polygon with 3 vertices
            sage: T.edge_cosh_squared(-1)
            Traceback (most recent call last):
            ...
            ValueError: the edge -1 is not an edge of a polygon with 3 vertices
            sage: T.edge_cosh_squared(1/2)
            Traceback (most recent call last):
            ...
            ValueError: the edge 1/2 is not an edge of a polygon with 3 vertices

        """
        from sage.all import ZZ

        n = len(self._v)
        if i not in ZZ or not 0 <= i < n:
            raise ValueError(
                f"the edge {i} is not an edge of a polygon with {n} vertices"
            )
        i = int(i)
        return _cosh2(self._v[i], self._v[(i + 1) % n], self.parent()._form)

    def _cosh2_sides(self):
        r"""
        Return the squares of the hyperbolic cosines of the side lengths of
        this polygon.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P._cosh2_sides()
            (4/3, 3/2, 3/2, 4/3)

        """
        n = len(self._v)
        form = self.parent()._form
        return tuple(_cosh2(self._v[i], self._v[(i + 1) % n], form) for i in range(n))

    def _angle_key(self, i):
        r"""
        Return a triple of exact data that determines the interior angle at
        the vertex ``i``.

        The triple is the sign of the Euclidean turn at the vertex in the
        Klein model (`1` at a convex corner, `-1` at a reflex one, `0` at a
        marked vertex), the sign of the cosine of the angle `\theta_0 \in (0,
        \pi)` between the two sides, and the square of that cosine.

        ALGORITHM:

        With the vertex lifted to `Q = (x, y, 1)`, the tangent vector of the
        hyperboloid at `Q` pointing to another point `R` is the projection
        `R - \langle R, Q\rangle / \langle Q, Q\rangle Q` of `R` to the
        Minkowski orthogonal of `Q`; the cosine of `\theta_0` is the
        normalized Minkowski product of the two tangent vectors. The Klein
        chart preserves the orientation, so the interior angle is
        `\theta_0` at a convex corner and `2\pi - \theta_0` at a reflex one.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
            sage: P._angle_key(2)
            (-1, 0, 0)
            sage: P._angle_key(0)
            (1, 1, 169/361)

        """
        from flatsurf.geometry.euclidean import ccw

        n = len(self._v)
        p = self._v[(i - 1) % n]
        q = self._v[i % n]
        r = self._v[(i + 1) % n]

        turn = ccw(q - p, r - q)
        turn = (turn > 0) - (turn < 0)

        P = _lift(p)
        Q = _lift(q)
        R = _lift(r)

        form = self.parent()._form

        QQ_ = _minkowski(Q, Q, form)
        RQ = _minkowski(R, Q, form)
        PQ = _minkowski(P, Q, form)

        product = _minkowski(R, P, form) - RQ * PQ / QQ_
        norm_r = _minkowski(R, R, form) - RQ * RQ / QQ_
        norm_p = _minkowski(P, P, form) - PQ * PQ / QQ_

        sign = (product > 0) - (product < 0)

        return (turn, sign, product**2 / (norm_r * norm_p))

    def angle(self, i, numerical=None):
        r"""
        Return the interior angle at the vertex ``i`` of this polygon, as a
        multiple of `2\pi`.

        INPUT:

        - ``i`` -- an integer, taken modulo the number of vertices

        - ``numerical`` -- a boolean or ``None`` (default: ``None``); whether
          to return a floating point approximation instead of an exact
          rational number. When ``None``, the angle is exact if the polygon
          is :meth:`rational <is_rational>` and numerical with a warning
          otherwise.

        OUTPUT:

        A number in `(0, 1)`.

        EXAMPLES:

        A Lambert quadrilateral, with three right angles::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P.angle(0, numerical=False)
            1/4
            sage: P.angle(2, numerical=True)  # abs tol 1e-12
            0.195913276015304
            sage: a = P.angle(2)
            doctest:warning
            ...
            UserWarning: for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the angle
            sage: a  # abs tol 1e-12
            0.195913276015304
            sage: P.angle(2, numerical=False)
            Traceback (most recent call last):
            ...
            ValueError: angle is not a rational multiple of 2π; set numerical=True to get a numerical approximation

        A reflex angle::

            sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
            sage: P.angle(2, numerical=False)
            3/4

        """
        if numerical is None:
            numerical = not self.is_rational()

            if numerical:
                import warnings

                warnings.warn(
                    "for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the angle"
                )

        turn, sign, cos2 = self._angle_key(i)

        if numerical:
            import math

            if turn == 0:
                return 0.5

            cos = sign * math.sqrt(float(cos2))
            theta = math.acos(min(1.0, max(-1.0, cos))) / (2 * math.pi)
            return theta if turn > 0 else 1 - theta

        from sage.all import QQ

        if turn == 0:
            return QQ((1, 2))

        theta = _rational_arccos(sign, cos2)
        if theta is None:
            raise ValueError(
                "angle is not a rational multiple of 2π; set numerical=True to get a numerical approximation"
            )

        return theta if turn > 0 else 1 - theta

    def angles(self, numerical=None):
        r"""
        Return the interior angles of this polygon, as multiples of `2\pi`.

        INPUT:

        - ``numerical`` -- a boolean or ``None`` (default: ``None``); as in
          :meth:`angle`

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
            sage: P.angles(numerical=True)  # abs tol 1e-12
            (0.130073469146503, 0.0368407666897357, 0.75, 0.0368407666897357)

        """
        if numerical is None:
            numerical = not self.is_rational()

            if numerical:
                import warnings

                warnings.warn(
                    "for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the angles"
                )

        return tuple(self.angle(i, numerical=numerical) for i in range(len(self._v)))

    def area(self, numerical=None):
        r"""
        Return the hyperbolic area of this polygon divided by `2\pi`.

        INPUT:

        - ``numerical`` -- a boolean or ``None`` (default: ``None``); as in
          :meth:`angle`

        ALGORITHM:

        By the Gauss-Bonnet theorem, the area of a hyperbolic `n`-gon is
        `(n - 2)\pi` minus the sum of its interior angles.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)])
            sage: P.area(numerical=True)  # abs tol 1e-12
            0.0462449974740256

        The area agrees with the area of the convex hull of the vertices
        for a convex polygon::

            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P.area(numerical=True)  # abs tol 1e-12
            0.0540867239846964
            sage: H.convex_hull(*P.vertices()).area()  # abs tol 1e-12
            0.0540867239846964

        A Lambert quadrilateral of the plane of the form `x^2 + 2y^2 - 3z^2`,
        whose vertex `(1, 0)` lies outside of the unit disk::

            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: L = H3.polygon(vertices=[(0, 0), (1, 0), (1, 1/2), (0, 1/2)])
            sage: [L.edge_cosh_squared(i) for i in range(4)]
            [3/2, 4/3, 5/3, 6/5]
            sage: L.angles(numerical=True)  # abs tol 1e-12
            (0.25, 0.25, 0.198791808825217, 0.25)
            sage: L.area(numerical=True)  # abs tol 1e-12
            0.0512081911747834
            sage: H3.convex_hull(*L.vertices()).area()  # abs tol 1e-12
            0.0512081911747834

        """
        from sage.all import QQ

        if numerical is None:
            numerical = not self.is_rational()

            if numerical:
                import warnings

                warnings.warn(
                    "for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the area"
                )

        n = len(self._v)
        total = sum(self.angles(numerical=numerical))

        if numerical:
            return (n - 2) / 2 - float(total)

        return QQ((n - 2, 2)) - total

    @cached_method
    def is_rational(self):
        r"""
        Return whether all the interior angles of this polygon are rational
        multiples of `\pi`.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]).is_rational()
            False

        A triangle with angles `\pi/2`, `\pi/6` and `\pi/6`::

            sage: K.<a> = NumberField(x^2 - 6, embedding=2.4)
            sage: H = HyperbolicPlane(K)
            sage: P = H.polygon(vertices=[(0, 0), (a/3, 0), (0, a/3)])
            sage: P.is_rational()
            True
            sage: P.angles()
            (1/4, 1/12, 1/12)
            sage: P.area()
            1/12

        """
        n = len(self._v)
        for i in range(n):
            turn, sign, cos2 = self._angle_key(i)
            if turn != 0 and _rational_arccos(sign, cos2) is None:
                return False

        return True

    def is_equilateral(self):
        r"""
        Return whether all the sides of this polygon have the same hyperbolic
        length.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]).is_equilateral()
            False
            sage: H.polygon(vertices=[(1/2, 0), (0, 1/2), (-1/2, 0), (0, -1/2)]).is_equilateral()
            True

        """
        return len(set(self._cosh2_sides())) == 1

    def is_equiangular(self):
        r"""
        Return whether all the interior angles of this polygon are equal.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]).is_equiangular()
            False
            sage: H.polygon(vertices=[(1/2, 0), (0, 1/2), (-1/2, 0), (0, -1/2)]).is_equiangular()
            True

        """
        return len({self._angle_key(i) for i in range(len(self._v))}) == 1

    def is_right_triangle(self):
        r"""
        Return whether this polygon is a triangle with a right angle.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]).is_right_triangle()
            True
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2)]).is_right_triangle()
            True
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/4, 1/2)]).is_right_triangle()
            False
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]).is_right_triangle()
            False

        """
        if len(self._v) != 3:
            return False

        return any(self._angle_key(i)[1] == 0 for i in range(3))

    def is_isosceles_triangle(self):
        r"""
        Return whether this polygon is a triangle with two sides of the same
        hyperbolic length.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]).is_isosceles_triangle()
            True
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/3)]).is_isosceles_triangle()
            False
            sage: H.polygon(vertices=[(1/2, 0), (0, 1/2), (-1/2, 0), (0, -1/2)]).is_isosceles_triangle()
            False

        """
        if len(self._v) != 3:
            return False

        return len(set(self._cosh2_sides())) < 3

    def _is_right_angled(self):
        r"""
        Return whether all the interior angles of this polygon are right
        angles.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])._is_right_angled()
            False

        """
        return all(self._angle_key(i)[:2] == (1, 0) for i in range(len(self._v)))

    def describe_polygon(self):
        r"""
        Return a textual description of this polygon for generating
        human-readable messages.

        The description is returned as a triple (indeterminate article,
        singular, plural). The Euclidean names of quadrilaterals (square,
        rectangle, rhombus) are not used.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: H.polygon(vertices=[(1/2, 0), (0, 1/2), (-1/2, 0), (0, -1/2)]).describe_polygon()
            ('a', 'regular quadrilateral', 'regular quadrilaterals')
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)]).describe_polygon()
            ('a', 'quadrilateral', 'quadrilaterals')
            sage: H.polygon(vertices=[(0, -1/2), (1/2, 1/2), (0, 0), (-1/2, 1/2)]).describe_polygon()
            ('a', 'non-convex quadrilateral', 'non-convex quadrilaterals')
            sage: H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]).describe_polygon()
            ('an', 'isosceles triangle', 'isosceles triangles')
            sage: H.polygon(vertices=[(1/2, 0), (0, 1/2), (-1/2, 0), (-1/4, -1/4), (0, -1/2)]).describe_polygon()
            ('a', 'regular quadrilateral with a marked vertex', 'regular quadrilaterals with a marked vertex')

        """
        from flatsurf.geometry.categories.polygons import Polygons
        from flatsurf.geometry.euclidean import polygon_corners

        polygon = self
        marked_vertices = len(self._v) - len(polygon_corners(self._v))

        if marked_vertices and self._klein_area() != 0:
            polygon = self.erase_marked_vertices()

        n = len(polygon._v)

        properties = {
            "degenerate": polygon.is_degenerate(),
            "equilateral": polygon.is_equilateral(),
            "equiangular": polygon.is_equiangular(),
            "convex": polygon.is_convex(),
            "marked_vertices": marked_vertices,
            "euclidean": False,
            "right_angled": polygon._is_right_angled(),
        }

        if n == 3:
            properties["right"] = polygon.is_right_triangle()
            properties["isosceles"] = polygon.is_isosceles_triangle()

        return Polygons._describe_polygon(n, **properties)

    def _test_gauss_bonnet(self, **options):
        r"""
        Verify that the angles lie in `(0, 1)` and that the area obtained from
        them is positive, and, for a convex polygon, equal to the area of the
        convex hull of the vertices.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (1/2, 1/2), (0, 1/2)])
            sage: P._test_gauss_bonnet()

        """
        tester = self._tester(**options)

        angles = self.angles(numerical=True)
        tester.assertTrue(all(0 < angle < 1 for angle in angles))

        area = self.area(numerical=True)
        tester.assertGreater(area, 0)

        if self.is_convex():
            hull = self.parent().convex_hull(*self.vertices())
            tester.assertAlmostEqual(area, float(hull.area(numerical=True)))

    def _test_vertex(self, **options):
        r"""
        Verify that :meth:`vertex` and :meth:`vertices` are compatible.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P._test_vertex()

        """
        tester = self._tester(**options)

        vertices = self.vertices()
        n = len(vertices)
        for i in range(-1, n + 1):
            tester.assertEqual(self.vertex(i), vertices[i % n])

    def _test_edge_segments(self, **options):
        r"""
        Verify that the edge segments go from a vertex to the next one.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H = HyperbolicPlane(QQ)
            sage: P = H.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: P._test_edge_segments()

        """
        tester = self._tester(**options)

        segments = self.edge_segments()
        tester.assertEqual(len(segments), len(self.vertices()))
        for i, segment in enumerate(segments):
            tester.assertEqual(segment.start(), self.vertex(i))
            tester.assertEqual(segment.end(), self.vertex(i + 1))
