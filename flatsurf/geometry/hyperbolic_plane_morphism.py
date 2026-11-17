r"""
Isomorphisms between hyperbolic planes of different quadratic forms.

A :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane` with the quadratic
form `Q = \mathrm{diag}(1, b, -c)` is the Klein model of the cone of `Q`: the
disk `x^2 + b y^2 < c`. A 3×3 matrix `A` with `A^T Q_2 A = \lambda Q_1`,
`\lambda > 0`, sends the cone of `Q_1` onto the cone of `Q_2`, and the Klein
point `(x, y, 1)` to `A (x, y, 1)`; it is an isometry from the plane of `Q_1`
to the plane of `Q_2`. It keeps the orientation when `\det A > 0` and the last
coordinate of `A (0, 0, 1)` is positive.

EXAMPLES:

The planes of the forms `x^2 + 8y^2 - 12z^2` and `x^2 + 2y^2 - 3z^2` differ by
a rescaling of the axes::

    sage: from flatsurf import HyperbolicPlane
    sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
    sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
    sage: f = H8.isomorphism(H3)
    sage: f
    Hyperbolic plane morphism:
      From: Hyperbolic Plane over Rational Field with form x^2 + 8*y^2 - 12*z^2
      To:   Hyperbolic Plane over Rational Field with form x^2 + 2*y^2 - 3*z^2
      Defn: [1 0 0]
            [0 2 0]
            [0 0 2]
    sage: f(H8.point(1/4, 1/8, model="klein"))
    (1/8, 1/8)

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

from sage.categories.morphism import Morphism


class HyperbolicPlaneIsomorphism(Morphism):
    r"""
    An orientation preserving isometry from a hyperbolic plane to a
    hyperbolic plane over the same base ring, of another quadratic form.

    INPUT:

    - ``parent`` -- the homset ``Hom(H_1, H_2)`` of two
      :class:`~flatsurf.geometry.hyperbolic.HyperbolicPlane` over the same
      base ring `F`, of the forms `Q_1` and `Q_2`

    - ``matrix`` -- a 3×3 matrix `A` over `F` with `A^T Q_2 A = \lambda Q_1`
      and `\lambda > 0`, of positive determinant, the last coordinate of `A
      (0, 0, 1)` positive; it is defined up to a positive scalar

    - ``check`` -- a boolean (default: ``True``); whether to check the
      conditions on ``matrix``

    EXAMPLES:

    An isomorphism from the plane of `x^2 + 2y^2 - 3z^2` to the plane of `x^2
    + y^2 - z^2` over the rationals, which sends the isotropic vector `(1, 1,
    1)` to the isotropic vector `(1, 0, 1)`::

        sage: from flatsurf import HyperbolicPlane
        sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
        sage: H = HyperbolicPlane(QQ)
        sage: A3 = matrix([[25/4, 23/2, -69/4], [0, 6, -6], [-23/4, -25/2, 75/4]])
        sage: f = H3.isomorphism(H, matrix=A3)
        sage: A3.transpose() * H.quadratic_form() * A3 == 6 * H3.quadratic_form()
        True

    It sends a kite to a polygon with the same sides, angles and area::

        sage: P = H3.polygon(vertices=[(0, 0), (1/2, -1/2), (3/4, 0), (1/2, 1/2)])
        sage: Q = f(P)
        sage: [v.coordinates(model="klein") for v in Q.vertices()]
        [(-23/25, -8/25), (-53/59, -24/59), (-67/77, -32/77), (-67/77, -24/77)]
        sage: Q.parent() is H
        True
        sage: Q._cosh2_sides() == P._cosh2_sides()
        True
        sage: Q._cosh2_sides()
        (4/3, 49/39, 49/39, 4/3)
        sage: Q.angles(numerical=True)  # abs tol 1e-12
        (0.3040867239846964, 0.14102355421224375, 0.38101825537203543, 0.14102355421224375)
        sage: Q.area(numerical=True)  # abs tol 1e-12
        0.032847912218780584

    It conjugates the side pairing of the kite, from its edge 0 to its edge
    3, into the group of `x^2 + y^2 - z^2`::

        sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
        sage: G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
        sage: g = G3([[-1/3, -4/3, 0], [2/3, -1/3, 0], [0, 0, 1]])
        sage: f(g)
        [ -179/8   -71/6 -607/24]
        [   -7/2      -3    -9/2]
        [  181/8    73/6  617/24]
        sage: f(g).trace()
        1/3
        sage: f(g).parent()
        Group of orientation preserving isometries of the hyperbolic plane over Rational Field

    It sends the surface glued from the kite to a surface of the plane of `x^2
    + y^2 - z^2`::

        sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
        sage: S = MutableOrientedHyperbolicSurface(H3)
        sage: S.add_polygon(P)
        0
        sage: S.glue((0, 0), (0, 3))
        sage: S.glue((0, 1), (0, 2))
        sage: S.set_immutable()
        sage: T = f(S)
        sage: T.hyperbolic_plane() is H
        True
        sage: T.is_mutable()
        False
        sage: T.edge_transformation(0, 0) == f(g)
        True
        sage: sorted(T.angles(numerical=True))  # abs tol 1e-12
        [0.2820471084244875, 0.3040867239846964, 0.38101825537203543]

    TESTS::

        sage: from flatsurf.geometry.hyperbolic_plane_morphism import HyperbolicPlaneIsomorphism
        sage: isinstance(f, HyperbolicPlaneIsomorphism)
        True
        sage: TestSuite(T).run()

    """

    def __init__(self, parent, matrix, check=True):
        Morphism.__init__(self, parent)

        from sage.matrix.constructor import matrix as Matrix

        H1 = parent.domain()
        H2 = parent.codomain()
        F = H1.base_ring()

        if H2.base_ring() is not F:
            raise ValueError(
                f"the planes must have the same base ring but {F} is not {H2.base_ring()}"
            )

        A = Matrix(F, matrix)

        if check:
            self._check(H1, H2, A)

        A.set_immutable()
        self._matrix = A

    @staticmethod
    def _check(H1, H2, A):
        r"""
        Raise a ``ValueError`` if ``A`` is not the matrix of an orientation
        preserving isometry from ``H1`` to ``H2``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_plane_morphism import HyperbolicPlaneIsomorphism
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: HyperbolicPlaneIsomorphism._check(H8, H3, diagonal_matrix([1, 2, 2]))

        TESTS:

        Each condition on the matrix is checked::

            sage: HyperbolicPlaneIsomorphism._check(H8, H3, identity_matrix(3))
            Traceback (most recent call last):
            ...
            ValueError: matrix must send the form of the codomain to a positive multiple of the form of the domain
            sage: HyperbolicPlaneIsomorphism._check(H8, H3, -diagonal_matrix([1, 2, 2]))
            Traceback (most recent call last):
            ...
            ValueError: matrix must send the point (0, 0, 1) to a positive last coordinate
            sage: HyperbolicPlaneIsomorphism._check(H8, H3, diagonal_matrix([1, -2, 2]))
            Traceback (most recent call last):
            ...
            ValueError: matrix must have a positive determinant
            sage: HyperbolicPlaneIsomorphism._check(H8, H3, identity_matrix(2))
            Traceback (most recent call last):
            ...
            ValueError: matrix must be 3×3

        """
        if A.dimensions() != (3, 3):
            raise ValueError("matrix must be 3×3")

        Q1 = H1.quadratic_form()
        D = A.transpose() * H2.quadratic_form() * A
        # D has the signature of Q1 when it is a nonzero multiple of Q1, so the
        # multiple D[0, 0] cannot be negative.
        if not D[0, 0] > 0 or D != D[0, 0] * Q1:
            raise ValueError(
                "matrix must send the form of the codomain to a positive multiple of the form of the domain"
            )

        if not A[2, 2] > 0:
            raise ValueError(
                "matrix must send the point (0, 0, 1) to a positive last coordinate"
            )

        if not A.det() > 0:
            raise ValueError("matrix must have a positive determinant")

    def matrix(self):
        r"""
        Return the matrix of this isomorphism, defined up to a positive scalar.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H8.isomorphism(H3).matrix()
            [1 0 0]
            [0 2 0]
            [0 0 2]

        """
        return self._matrix

    def _repr_type(self):
        r"""
        Return the type of this morphism as printed.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.isomorphism(H3)._repr_type()
            'Hyperbolic plane'

        """
        return "Hyperbolic plane"

    def _repr_defn(self):
        r"""
        Return the definition of this morphism as printed, its matrix.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.isomorphism(H3)._repr_defn()
            '[1 0 0]\n[0 1 0]\n[0 0 1]'

        """
        return repr(self._matrix)

    def is_identity(self):
        r"""
        Return whether this isomorphism is the identity of its domain.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.isomorphism(H3).is_identity()
            True
            sage: H3.isomorphism(H3, matrix=diagonal_matrix([2, 2, 2])).is_identity()
            True
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H8.isomorphism(H3).is_identity()
            False

        """
        A = self._matrix
        return self.domain() is self.codomain() and A.is_scalar(A[0, 0])

    def _normalized_matrix(self):
        r"""
        Return the matrix of this isomorphism divided by the absolute value of
        its first nonzero entry.

        Two matrices that differ by a positive scalar have the same
        normalization.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.isomorphism(H3, matrix=diagonal_matrix([3, 3, 3]))._normalized_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        A = self._matrix
        pivot = next(entry for entry in A.list() if entry)
        if pivot < 0:
            pivot = -pivot
        A = A / pivot
        A.set_immutable()
        return A

    def _richcmp_(self, other, op):
        r"""
        Compare this isomorphism to ``other``, an isomorphism with the same
        domain and codomain; they are equal when their matrices differ by a
        positive scalar.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H = HyperbolicPlane(QQ)
            sage: A3 = matrix([[25/4, 23/2, -69/4], [0, 6, -6], [-23/4, -25/2, 75/4]])
            sage: H3.isomorphism(H, matrix=A3) == H3.isomorphism(H, matrix=2*A3)
            True
            sage: H3.isomorphism(H3) != H3.isomorphism(H3, matrix=diagonal_matrix([2, 2, 2]))
            False
            sage: H3.isomorphism(H, matrix=A3) == H.isomorphism(H)
            False

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op not in (op_EQ, op_NE):
            return NotImplemented

        equal = (
            isinstance(other, HyperbolicPlaneIsomorphism)
            and self.domain() is other.domain()
            and self.codomain() is other.codomain()
            and self._normalized_matrix() == other._normalized_matrix()
        )
        return equal == (op == op_EQ)

    def __hash__(self):
        r"""
        Return a hash value for this isomorphism, compatible with ``==``.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: hash(H3.isomorphism(H3)) == hash(H3.isomorphism(H3, matrix=diagonal_matrix([2, 2, 2])))
            True

        """
        return hash((self.domain(), self.codomain(), self._normalized_matrix()))

    def __invert__(self):
        r"""
        Return the inverse of this isomorphism, of matrix `A^{-1}`.

        EXAMPLES:

        Over `\QQ(\sqrt{3})`, an isomorphism to the plane of `x^2 + y^2 - z^2`
        and its inverse::

            sage: from flatsurf import HyperbolicPlane
            sage: K.<a> = NumberField(x^2 - 3, embedding=AA(3).sqrt())
            sage: H3K = HyperbolicPlane(K, form=[1, 2, -3])
            sage: HK = HyperbolicPlane(K)
            sage: f = H3K.isomorphism(HK, matrix=matrix(K, [[7/12*a, 0, -5/4], [0, 1, 0], [-5/12*a, 0, 7/4]]))
            sage: P = H3K.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])
            sage: Q = f(P)
            sage: [v.coordinates(model="klein") for v in Q.vertices()]
            [(-5/7, 0), (48/563*a - 385/563, 0), (-5/7, 2/7)]
            sage: Q.angles(numerical=True)  # abs tol 1e-12
            (0.25, 0.14875643774533862, 0.09122814004003268)
            sage: P.angles(numerical=True)  # abs tol 1e-12
            (0.25, 0.14875643774533862, 0.09122814004003268)
            sage: (~f)(Q) == P
            True
            sage: ~f
            Hyperbolic plane morphism:
              From: Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
              To:   Hyperbolic Plane over Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878? with form x^2 + 2*y^2 - 3*z^2
              Defn: [7/6*a     0 5/6*a]
                    [    0     1     0]
                    [  5/6     0   7/6]

        A polygon of the plane of `x^2 + 2y^2 - 3z^2` over `\QQ` is moved to
        the plane over `\QQ(\sqrt{3})` first::

            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f(H3.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)])) == Q
            True

        The composition with the inverse is the identity::

            sage: (~f) * f == H3K.isomorphism(H3K)
            True
            sage: f * (~f) == HK.isomorphism(HK)
            True

        """
        from sage.categories.homset import Hom

        A = self._matrix.inverse().change_ring(self.domain().base_ring())
        return HyperbolicPlaneIsomorphism(
            Hom(self.codomain(), self.domain()), A, check=False
        )

    def section(self):
        r"""
        Return the inverse of this isomorphism.

        This is the same as ``~self``; surface morphisms spell their inverse
        ``section()``, so that code written for both works.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)
            sage: f.section() == ~f
            True
            sage: f.section().matrix()
            [  1   0   0]
            [  0 1/2   0]
            [  0   0 1/2]

        """
        return ~self

    def _composition_(self, right, homset):
        r"""
        Return the composition of this isomorphism with ``right``, ``right``
        applied first.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H = HyperbolicPlane(QQ)
            sage: A3 = matrix([[25/4, 23/2, -69/4], [0, 6, -6], [-23/4, -25/2, 75/4]])
            sage: g = H3.isomorphism(H, matrix=A3) * H8.isomorphism(H3)
            sage: g.matrix()
            [ 25/4    23 -69/2]
            [    0    12   -12]
            [-23/4   -25  75/2]
            sage: g.domain() is H8 and g.codomain() is H
            True

        """
        if isinstance(right, HyperbolicPlaneIsomorphism):
            return HyperbolicPlaneIsomorphism(
                homset, self._matrix * right._matrix, check=False
            )

        return super()._composition_(right, homset)

    def _call_(self, x):
        r"""
        Return the image of the point or convex set ``x`` of the domain.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)
            sage: f._call_(H8.point(1/4, 1/8, model="klein"))
            (1/8, 1/8)

        """
        return x._apply_isometry_klein(self._matrix, codomain=self.codomain())

    def __call__(self, x):
        r"""
        Return the image of ``x`` under this isomorphism.

        INPUT:

        - ``x`` -- a point or a convex set of the domain, a
          :class:`~flatsurf.geometry.hyperbolic_polygon.HyperbolicPolygon`, an
          isometry of the group of the domain, or a hyperbolic surface whose
          polygons live in the domain; an object of a plane with the form of
          the domain over a subring of the base ring is moved to the domain
          first

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)

        Points and convex sets::

            sage: p = H8.point(1/4, 1/8, model="klein")
            sage: f(p)
            (1/8, 1/8)
            sage: f(p).parent() is H3
            True
            sage: q = H8.point(0, 0, model="klein")
            sage: f(p.segment(q)) == f(p).segment(f(q))
            True
            sage: f(H8.geodesic(p, q)) == H3.geodesic(f(p), f(q))
            True
            sage: f(H8.geodesic(p, q).left_half_space()) == H3.geodesic(f(p), f(q)).left_half_space()
            True
            sage: f(H8.empty_set()) == H3.empty_set()
            True

        Polygons::

            sage: f(H8.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            Hyperbolic polygon with vertices [(0, 0), (1/4, 0), (0, 1/2)]

        TESTS::

            sage: f(42)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot map Integer through a hyperbolic plane isomorphism
            sage: f(H3.point(0, 0, model="klein"))
            Traceback (most recent call last):
            ...
            ValueError: (0, 0) does not lie in the domain of this isomorphism

        """
        from flatsurf.geometry.categories import HyperbolicConeSurfaces
        from flatsurf.geometry.hyperbolic import HyperbolicConvexSet
        from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometry
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon
        from flatsurf.geometry.surface import Surface_base

        H1 = self.domain()
        F = H1.base_ring()

        if isinstance(x, (HyperbolicConvexSet, HyperbolicPolygon)):
            self._require_domain(x, x.parent())
            if isinstance(x, HyperbolicPolygon):
                return self._image_polygon(x.change_ring(F))
            return self._call_(H1(x))

        if isinstance(x, HyperbolicIsometry):
            self._require_domain(x, x.parent().hyperbolic_plane())
            return self._image_isometry(x)

        if isinstance(x, Surface_base) and x in HyperbolicConeSurfaces().Oriented():
            self._require_domain(x, x.hyperbolic_plane())
            if x.base_ring() is not F:
                x = x.change_ring(F)
            return self._image_surface(x)

        raise NotImplementedError(
            f"cannot map {type(x).__name__} through a hyperbolic plane isomorphism"
        )

    def _require_domain(self, x, plane):
        r"""
        Raise a ``ValueError`` if the object ``x`` of the hyperbolic ``plane``
        does not coerce into the domain of this isomorphism.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H3.isomorphism(H3)
            sage: f._require_domain("x", HyperbolicPlane(QQ))
            Traceback (most recent call last):
            ...
            ValueError: x does not lie in the domain of this isomorphism

        """
        if not self.domain().has_coerce_map_from(plane):
            raise ValueError(f"{x} does not lie in the domain of this isomorphism")

    def _image_polygon(self, P):
        r"""
        Return the image of the polygon ``P`` of the domain: the polygon of
        the codomain whose vertices are the images of its vertices, in the
        same order.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)
            sage: f._image_polygon(H8.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            Hyperbolic polygon with vertices [(0, 0), (1/4, 0), (0, 1/2)]

        """
        from flatsurf.geometry.hyperbolic_polygon import HyperbolicPolygon

        A = self._matrix
        images = []
        for x, y in P._v:
            X = A[0, 0] * x + A[0, 1] * y + A[0, 2]
            Y = A[1, 0] * x + A[1, 1] * y + A[1, 2]
            Z = A[2, 0] * x + A[2, 1] * y + A[2, 2]
            images.append((X / Z, Y / Z))

        # The map is projective with Z > 0 on the disk and det A > 0, so it
        # keeps the orientation, the convexity and the simplicity of P.
        return HyperbolicPolygon(self.codomain(), images)

    def _image_isometry(self, g):
        r"""
        Return the image `A g A^{-1}` of the isometry ``g`` of the plane of
        the domain, in the group of the codomain.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)
            sage: g = HyperbolicIsometryGroup(QQ, form=[1, 8, -12]).some_elements()[2]
            sage: f._image_isometry(g)
            [-7/9  8/9    0]
            [-4/9 -7/9    0]
            [   0    0    1]

        """
        G1, G2 = self._isometry_groups()
        A = self._matrix
        return G2(A * G1(g).matrix() * A.inverse())

    def _isometry_groups(self):
        r"""
        Return the groups of orientation preserving isometries of the domain
        and of the codomain.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: H3.isomorphism(H3)._isometry_groups()
            (Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2,
             Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2)

        """
        from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup

        F = self.domain().base_ring()
        return (
            HyperbolicIsometryGroup(F, form=self.domain()._form),
            HyperbolicIsometryGroup(F, form=self.codomain()._form),
        )

    def isometry_group_morphism(self):
        r"""
        Return the morphism of groups `g \mapsto A g A^{-1}` from the
        isometries of the domain to the isometries of the codomain.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: f = H8.isomorphism(H3)
            sage: phi = f.isometry_group_morphism()
            sage: phi.domain()
            Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 8*y^2 - 12*z^2
            sage: phi.codomain()
            Group of orientation preserving isometries of the hyperbolic plane over Rational Field with form x^2 + 2*y^2 - 3*z^2
            sage: g = phi.domain().an_element()
            sage: phi(g) == f(g)
            True
            sage: phi(g).trace() == g.trace()
            True

        """
        from sage.categories.homset import Hom
        from sage.categories.morphism import SetMorphism

        G1, G2 = self._isometry_groups()
        return SetMorphism(Hom(G1, G2), self._image_isometry)

    def _image_surface(self, S):
        r"""
        Return the image of the hyperbolic surface ``S`` whose polygons live
        in the domain: the surface of the codomain with the same labels, the
        images of the polygons and the same gluings, immutable when ``S`` is.

        EXAMPLES::

            sage: from flatsurf import HyperbolicPlane
            sage: from flatsurf.geometry.hyperbolic_surface import MutableOrientedHyperbolicSurface
            sage: H8 = HyperbolicPlane(QQ, form=[1, 8, -12])
            sage: H3 = HyperbolicPlane(QQ, form=[1, 2, -3])
            sage: S = MutableOrientedHyperbolicSurface(H8)
            sage: S.add_polygon(H8.polygon(vertices=[(0, 0), (1/2, 0), (0, 1/2)]))
            0
            sage: T = H8.isomorphism(H3)._image_surface(S)
            sage: T.polygon(0)
            Hyperbolic polygon with vertices [(0, 0), (1/4, 0), (0, 1/2)]
            sage: T.is_mutable()
            True

        """
        from flatsurf.geometry.hyperbolic_surface import (
            MutableOrientedHyperbolicSurface,
        )

        T = MutableOrientedHyperbolicSurface(self.codomain())
        for label in S.labels():
            T.add_polygon(self._image_polygon(S.polygon(label)), label=label)

        glued = set()
        for x, y in S.gluings():
            if x in glued:
                continue
            T.glue(x, y)
            glued.add(x)
            glued.add(y)

        if not S.is_mutable():
            T.set_immutable()

        return T
