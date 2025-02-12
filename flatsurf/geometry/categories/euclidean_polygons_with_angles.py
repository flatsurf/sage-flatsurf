r"""
The category of polygons in the real plane with fixed rational
angles.

This module provides a common structure for all polygons with certain fixed
angles.

See :mod:`flatsurf.geometry.categories` for a general description of the
category framework in sage-flatsurf.

Normally, you won't create this (or any other) category directly. The correct
category is automatically determined for polygons.

EXAMPLES:

The category of rectangles::

    sage: from flatsurf.geometry.categories import EuclideanPolygons
    sage: C = EuclideanPolygons(QQ).WithAngles([1, 1, 1, 1])

It is often tedious to create this category manually, since you need to
determine a base ring that can describe the coordinates of polygons with such
angles::

    sage: C = EuclideanPolygons(QQ).WithAngles([1, 1, 1])
    sage: C.slopes()
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce c to a rational

    sage: C = EuclideanPolygons(AA).WithAngles([1, 1, 1])
    sage: C.slopes()
    [(1, 0), (-0.5773502691896258?, 1), (-0.5773502691896258?, -1)]

Instead, we can use :func:`~.polygon.EuclideanPolygonsWithAngles` to create this category
over a minimal number field::

    sage: from flatsurf import EuclideanPolygonsWithAngles
    sage: C = EuclideanPolygonsWithAngles([1, 1, 1])
    sage: C.slopes()
    [(1, 0), (-c, 3), (-c, -3)]

The category of polygons is automatically determined when using
:func:`~.polygon.Polygon`::

    sage: from flatsurf import Polygon
    sage: p = Polygon(angles=(1, 1, 1))
    sage: p.category()
    Category of convex simple euclidean equilateral triangles over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

However, it can be very costly to determine that a polygon is rational and what
its actual angles are (the "equilateral" in the previous example). Therefore,
the category might get refined once these aspects have been determined::

    sage: p = Polygon(edges=[(1, 0), (0, 1), (-1, 0), (0, -1)])
    sage: p.category()
    Category of convex simple euclidean polygons over Rational Field
    sage: p.is_rational()
    True
    sage: p.category()
    Category of rational convex simple euclidean polygons over Rational Field
    sage: p.angles()
    (1/4, 1/4, 1/4, 1/4)
    sage: p.category()
    Category of convex simple euclidean rectangles over Rational Field

Note that SageMath applies the same strategy when determining whether the
integers modulo N are a field::

    sage: K = Zmod(1361)
    sage: K.category()
    Join of Category of finite commutative rings and Category of subquotients of monoids and Category of quotients of semigroups and Category of finite enumerated sets
    sage: K.is_field()
    True
    sage: K.category()
    Join of Category of finite enumerated fields and Category of subquotients of monoids and Category of quotients of semigroups

"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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
from sage.misc.cachefunc import cached_method, cached_function
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from flatsurf.geometry.categories.euclidean_polygons import EuclideanPolygons


class EuclideanPolygonsWithAngles(Category_over_base_ring):
    r"""
    The category of euclidean polygons with fixed rational angles.

    EXAMPLES::

        sage: from flatsurf.geometry.categories import EuclideanPolygons
        sage: C = EuclideanPolygons(QQ).WithAngles([1, 1, 1, 1])

    TESTS::

        sage: TestSuite(C).run()

    """

    def __init__(self, base_ring, angles):
        self._angles = angles

        super().__init__(base_ring)

    def super_categories(self):
        r"""
        Return the other categories such polygons are automatically members of,
        namely, the category of rational euclidean polygons polygons.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: C = EuclideanPolygons(QQ).WithAngles([1, 1, 1, 1])
            sage: C.super_categories()
            [Category of rational euclidean polygons over Rational Field]

        """
        return [EuclideanPolygons(self.base_ring()).Rational()]

    @staticmethod
    def _normalize_angles(angles):
        r"""
        Return ``angles`` normalized such that they sum to n/2-1 where ``n`` is
        the number of angles, i.e., they scale to the angle in an n-gon divided
        by π.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygonsWithAngles
            sage: EuclideanPolygonsWithAngles._normalize_angles([1, 1, 1, 1])
            (1/4, 1/4, 1/4, 1/4)
            sage: EuclideanPolygonsWithAngles._normalize_angles([1, 2, 3, 4])
            (1/10, 1/5, 3/10, 2/5)
            sage: EuclideanPolygonsWithAngles._normalize_angles([1, 2, 3])
            (1/12, 1/6, 1/4)

        """
        n = len(angles)
        if n < 3:
            raise ValueError("there must be at least three angles")

        from sage.all import QQ, ZZ

        angles = [QQ.coerce(a) for a in angles]
        if any(angle <= 0 for angle in angles):
            raise ValueError("angles must be positive rationals")

        # Store each angle as a multiple of 2π, i.e., normalize them such their sum is (n - 2)/2.
        angles = [a / sum(angles) for a in angles]
        angles = [a * ZZ(n - 2) / 2 for a in angles]
        if any(angle <= 0 or angle >= 1 for angle in angles):
            raise NotImplementedError("each angle must be in (0, 2π)")

        angles = tuple(angles)

        return angles

    @cached_method
    def _slopes(self):
        slopes = _slopes(self._angles)

        # We bring the slopes first into the minimal number field in which they
        # are defined since otherwise conversion from the cosines ring to the
        # base_ring might fail. E.g., when the base ring is the exact-reals
        # over the minimal base ring.
        minimal_base_ring = _base_ring(self._angles)
        slopes = [
            (minimal_base_ring(slope[0]), minimal_base_ring(slope[1]))
            for slope in slopes
        ]
        return [
            (self.base_ring()(slope[0]), self.base_ring()(slope[1])) for slope in slopes
        ]

    @cached_method
    def _cosines_ring(self):
        slopes = _slopes(self._angles)
        return slopes[0][0].parent()

    def _repr_object_names(self):
        r"""
        Helper method to create the name of this category.

        EXAMPLES::

            sage: from flatsurf import EuclideanPolygonsWithAngles
            sage: EuclideanPolygonsWithAngles([1/6, 1/6, 1/6])
            Category of simple euclidean equilateral triangles over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?
            sage: EuclideanPolygonsWithAngles([1/4, 1/4, 1/4, 1/4])
            Category of simple euclidean rectangles over Rational Field
            sage: EuclideanPolygonsWithAngles([1/10, 2/10, 3/10, 4/10])
            Category of simple euclidean quadrilaterals with angles (1/10, 1/5, 3/10, 2/5) over Number Field in c with defining polynomial x^4 - 5*x^2 + 5 with c = 1.902113032590308?

        """
        names = super()._repr_object_names()

        equiangular = len(set(self._angles)) == 1

        from flatsurf.geometry.categories.polygons import Polygons

        _, _, polygons = Polygons._describe_polygon(
            len(self._angles), equiangular=equiangular
        )

        with_angles = "" if equiangular else f" with angles {self.angles(False)}"

        return names.replace(" with angles", with_angles).replace("polygons", polygons)

    class ParentMethods:
        r"""
        Provides methods available to all polygons with known angles.

        If you want to add functionality to all such polygons, you probably
        want to put it here.
        """

        def is_convex(self, strict=False):
            r"""
            Return whether this is a convex polygon.

            INPUT:

            - ``strict`` -- whether to check for strict convexity, i.e., a
              polygon with a π angle is not considered convex.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.is_convex()
                True
                sage: S.is_convex(strict=True)
                True

            """
            return self.category().is_convex()

        def angle(self, e, numerical=None, assume_rational=None):
            r"""
            Return the angle at the beginning of the start point of the edge ``e``.

            EXAMPLES::

                sage: from flatsurf.geometry.polygon import polygons
                sage: polygons.square().angle(0)
                1/4
                sage: polygons.regular_ngon(8).angle(0)
                3/8

                sage: from flatsurf import Polygon
                sage: T = Polygon(vertices=[(0,0), (3,1), (1,5)])
                sage: [T.angle(i, numerical=True) for i in range(3)]  # abs tol 1e-13
                [0.16737532973071603, 0.22741638234956674, 0.10520828791971722]
                sage: sum(T.angle(i, numerical=True) for i in range(3))   # abs tol 1e-13
                0.5

            """
            if assume_rational is not None:
                import warnings

                warnings.warn(
                    "assume_rational has been deprecated as a keyword to angle() and will be removed from a future version of sage-flatsurf"
                )

            if numerical is None:
                numerical = not self.is_rational()

                if numerical:
                    import warnings

                    warnings.warn(
                        "the behavior of angle() has been changed in recent versions of sage-flatsurf; for non-rational polygons, numerical=True must be set explicitly to get a numerical approximation of the angle"
                    )

            angle = self.category().angles()[e]
            if numerical:
                from sage.all import RR

                angle = RR(angle)

            return angle

    class SubcategoryMethods:
        def is_convex(self, strict=False):
            r"""
            Return whether the polygons in this category are convex.

            INPUT:

            - ``strict`` -- a boolean (default: ``False``); whether only to
              consider polygons convex if all angles are <π.

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(1, 2, 5).is_convex()
                True
                sage: EuclideanPolygonsWithAngles(2, 2, 3, 13).is_convex()
                False

            ::

                sage: E = EuclideanPolygonsWithAngles([1, 1, 1, 1, 2])
                sage: E.angles()
                (1/4, 1/4, 1/4, 1/4, 1/2)
                sage: E.is_convex(strict=False)
                True
                sage: E.is_convex(strict=True)
                False

            """
            if strict:
                return all(2 * a < 1 for a in self.angles())

            return all(2 * a <= 1 for a in self.angles())

        def convexity(self):
            r"""
            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(1, 2, 5).convexity()
                doctest:warning
                ...
                UserWarning: convexity() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex() instead
                True
                sage: EuclideanPolygonsWithAngles(2, 2, 3, 13).convexity()
                False

            """
            import warnings

            warnings.warn(
                "convexity() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex() instead"
            )

            return self.is_convex()

        def strict_convexity(self):
            r"""
            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: E = EuclideanPolygonsWithAngles([1, 1, 1, 1, 2])
                sage: E.angles()
                (1/4, 1/4, 1/4, 1/4, 1/2)
                sage: E.convexity()
                True
                sage: E.strict_convexity()
                doctest:warning
                ...
                UserWarning: strict_convexity() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex(strict=True) instead
                False

            """
            import warnings

            warnings.warn(
                "strict_convexity() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex(strict=True) instead"
            )

            return self.is_convex(strict=True)

        def angles(self, integral=False):
            r"""
            Return the interior angles of this polygon as multiples of 2π.

            INPUT:

            - ``integral`` -- a boolean (default: ``False``); whether to return
              the angles not as multiples of 2π but rescaled so that they have
              no denominators.

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: E = EuclideanPolygonsWithAngles(1, 1, 1, 2, 6)
                sage: E.angles()
                (3/22, 3/22, 3/22, 3/11, 9/11)

            When ``integral`` is set, the output is scaled to eliminate
            denominators::

                sage: E.angles(integral=True)
                (1, 1, 1, 2, 6)

            """
            angles = self.__angles()
            if integral:
                from sage.all import lcm, ZZ, gcd

                C = lcm([a.denominator() for a in self.angles()]) / gcd(
                    [a.numerator() for a in self.angles()]
                )
                angles = tuple(ZZ(C * a) for a in angles)
            return angles

        @cached_method
        def __angles(self):
            r"""
            Helper method for :meth:`angles` to lookup the stored angles if
            this is a subcategory of :class:`EuclideanPolygonsWithAngles`.
            """
            if isinstance(self, EuclideanPolygonsWithAngles):
                return self._angles

            for category in self.all_super_categories():
                if isinstance(category, EuclideanPolygonsWithAngles):
                    return category._angles

            assert (
                False
            ), "EuclideanPolygonsWithAngles should be a supercategory of this category"

        def slopes(self, e0=(1, 0)):
            r"""
            Return the slopes of the edges as a list of vectors.

            INPUT:

            - ``e0`` -- the first slope returned (default: ``(1, 0)``)

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(1, 2, 1, 2).slopes()
                [(1, 0), (c, 3), (-1, 0), (-c, -3)]

            """
            V = self.base_ring() ** 2
            slopes = self.__slopes()
            n = len(slopes)
            cosines = [x[0] for x in slopes]
            sines = [x[1] for x in slopes]
            e = V(e0)
            edges = [e]
            for i in range(n - 1):
                e = (
                    -cosines[i + 1] * e[0] - sines[i + 1] * e[1],
                    sines[i + 1] * e[0] - cosines[i + 1] * e[1],
                )
                from flatsurf.geometry.euclidean import projectivization

                e = projectivization(*e)
                edges.append(V(e))
            return edges

        @cached_method
        def __slopes(self):
            r"""
            Helper method for :meth:`slopes` to lookup the stored slopes if
            this is a subcategory of :class:`EuclideanPolygonsWithAngles`.
            """
            if isinstance(self, EuclideanPolygonsWithAngles):
                return self._slopes()

            for category in self.all_super_categories():
                if isinstance(category, EuclideanPolygonsWithAngles):
                    return category._slopes()

            assert (
                False
            ), "EuclideanPolygonsWithAngles should be a supercategory of this category"

        # TODO: rather than lengths, it would be more convenient to have access
        # to the tangent space (that is the space of possible holonomies). However,
        # since it is not defined over the real numbers, there are several possible ways
        # to handle the data.
        # TODO: here we ignored the direction SO(2) which provides additional symmetry
        # in the tangent space
        @cached_method
        def lengths_polytope(self):
            r"""
            Return the polytope parametrizing the admissible vectors data.

            This polytope parametrizes the tangent space to the set of these
            equiangular polygons. Be careful that even though the lengths are
            admissible, they may not define a polygon without intersection.

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(1, 2, 1, 2).lengths_polytope()
                A 2-dimensional polyhedron in (Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?)^4 defined as the convex hull of 1 vertex and 2 rays
            """
            n = len(self.angles())
            slopes = self.slopes()
            eqns = [[0] + [s[0] for s in slopes], [0] + [s[1] for s in slopes]]
            ieqs = []
            for i in range(n):
                ieq = [0] * (n + 1)
                ieq[i + 1] = 1
                ieqs.append(ieq)

            from sage.geometry.polyhedron.constructor import Polyhedron

            return Polyhedron(eqns=eqns, ieqs=ieqs, base_ring=self.base_ring())

        def an_element(self):
            r"""
            Return a polygon in this category.

            Since currently polygons must not be self-intersecting, the
            construction used might fail.

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(4, 3, 4, 4, 3, 4).an_element()
                Polygon(vertices=[(0, 0),
                                  (1/22*c + 1, 0),
                                  (9*c^9 + 1/2*c^8 - 88*c^7 - 9/2*c^6 + 297*c^5 + 27/2*c^4 - 396*c^3 - 15*c^2 + 3631/22*c + 11/2, 1/2*c + 11),
                                  (16*c^9 + c^8 - 154*c^7 - 9*c^6 + 506*c^5 + 27*c^4 - 638*c^3 - 30*c^2 + 4841/22*c + 9, c + 22),
                                  (16*c^9 + c^8 - 154*c^7 - 9*c^6 + 506*c^5 + 27*c^4 - 638*c^3 - 30*c^2 + 220*c + 8, c + 22),
                                  (7*c^9 + 1/2*c^8 - 66*c^7 - 9/2*c^6 + 209*c^5 + 27/2*c^4 - 242*c^3 - 15*c^2 + 55*c + 7/2, 1/2*c + 11)])
            """
            from flatsurf import Polygon

            p = Polygon(angles=self.angles())

            if p not in self:  # pylint: disable=unsupported-membership-test
                raise NotImplementedError(
                    "cannot create an element in this category yet"
                )

            return p

        def random_element(self, ring=None, **kwds):
            r"""
            Return a random polygon in this category.

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: EuclideanPolygonsWithAngles(1, 1, 1, 2, 5).random_element()
                Polygon(vertices=[(0, 0), ...])
                sage: EuclideanPolygonsWithAngles(1,1,1,15,15,15).random_element()
                Polygon(vertices=[(0, 0), ...])
                sage: EuclideanPolygonsWithAngles(1,15,1,15,1,15).random_element()
                Polygon(vertices=[(0, 0), ...])

            """
            if ring is None:
                from sage.all import QQ

                ring = QQ

            rays = [r.vector() for r in self.lengths_polytope().rays()]

            def random_element():
                while True:
                    coeffs = []
                    while len(coeffs) < len(rays):
                        x = ring.random_element(**kwds)
                        while x < 0:
                            x = ring.random_element(**kwds)
                        coeffs.append(x)

                    sol = sum(c * r for c, r in zip(coeffs, rays))
                    if all(x > 0 for x in sol):
                        return coeffs, sol

            while True:
                coeffs, lengths = random_element()
                edges = [
                    length * slope for (length, slope) in zip(lengths, self.slopes())
                ]

                from flatsurf import Polygon

                p = Polygon(edges=edges, check=False)

                from flatsurf.geometry.categories import EuclideanPolygons

                if not EuclideanPolygons.ParentMethods.is_simple(p):
                    continue

                p = Polygon(edges=edges, angles=self.angles(), check=False)
                break

            if p not in self:  # pylint: disable=unsupported-membership-test
                raise NotImplementedError(
                    "cannot create a random element in this category yet"
                )

            return p

        def billiard_unfolding_angles(self, cover_type="translation"):
            r"""
            Return the angles of the unfolding rational, half-translation or translation surface.

            INPUT:

            - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
              ``"half-translation"`` or ``"translation"``

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles

                sage: E = EuclideanPolygonsWithAngles(1, 2, 5)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/8: 1, 1/4: 1, 5/8: 1}
                sage: (1/8 - 1) + (1/4 - 1) + (5/8 - 1)  # Euler characteristic (of the sphere)
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 3, 5/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 3, 5: 1}

                sage: E = EuclideanPolygonsWithAngles(1, 3, 1, 7)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/6: 2, 1/2: 1, 7/6: 1}
                sage: 2 * (1/6 - 1) + (1/2 - 1) + (7/6 - 1) # Euler characteristic
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 5, 7/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 5, 7: 1}

                sage: E = EuclideanPolygonsWithAngles(1, 3, 5, 7)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/8: 1, 3/8: 1, 5/8: 1, 7/8: 1}
                sage: (1/8 - 1) + (3/8 - 1) + (5/8 - 1) + (7/8 - 1) # Euler characteristic
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 1, 3/2: 1, 5/2: 1, 7/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 1, 3: 1, 5: 1, 7: 1}

                sage: E = EuclideanPolygonsWithAngles(1, 2, 8)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/11: 1, 2/11: 1, 8/11: 1}
                sage: (1/11 - 1) + (2/11 - 1) + (8/11 - 1) # Euler characteristic
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1: 1, 2: 1, 8: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 1, 2: 1, 8: 1}
            """
            rat_angles = {}
            for a in self.angles():
                if 2 * a in rat_angles:
                    rat_angles[2 * a] += 1
                else:
                    rat_angles[2 * a] = 1
            if cover_type == "rational":
                return rat_angles

            from sage.all import lcm

            N = lcm([x.denominator() for x in rat_angles])
            if N % 2:
                N *= 2

            cov_angles = {}
            for x, mult in rat_angles.items():
                y = x.numerator()
                d = x.denominator()
                if d % 2:
                    d *= 2
                else:
                    y = y / 2
                assert N % d == 0
                if y in cov_angles:
                    cov_angles[y] += mult * N // d
                else:
                    cov_angles[y] = mult * N // d

            if cover_type == "translation" and any(
                y.denominator() == 2 for y in cov_angles
            ):
                covcov_angles = {}
                for y, mult in cov_angles.items():
                    yy = y.numerator()
                    if yy not in covcov_angles:
                        covcov_angles[yy] = 0
                    covcov_angles[yy] += 2 // y.denominator() * mult
                return covcov_angles
            elif cover_type == "half-translation" or cover_type == "translation":
                return cov_angles
            else:
                raise ValueError("unknown 'cover_type' {!r}".format(cover_type))

        def billiard_unfolding_stratum(
            self, cover_type="translation", marked_points=False
        ):
            r"""
            Return the stratum of quadratic or Abelian differential obtained by
            unfolding a billiard in a polygon of this equiangular family.

            INPUT:

            - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
              ``"half-translation"`` or ``"translation"``

            - ``marked_poins`` (optional, default ``False``) - whether the stratum should
              have regular marked points

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles, similarity_surfaces

                sage: E = EuclideanPolygonsWithAngles(1, 2, 5)
                sage: E.billiard_unfolding_stratum("half-translation")
                Q_1(3, -1^3)
                sage: E.billiard_unfolding_stratum("translation")
                H_3(4)
                sage: E.billiard_unfolding_stratum("half-translation", True)
                Q_1(3, -1^3)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_3(4, 0^3)

                sage: E = EuclideanPolygonsWithAngles(1, 3, 1, 7)
                sage: E.billiard_unfolding_stratum("half-translation")
                Q_1(5, -1^5)
                sage: E.billiard_unfolding_stratum("translation")
                H_4(6)
                sage: E.billiard_unfolding_stratum("half-translation", True)
                Q_1(5, -1^5)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_4(6, 0^5)

                sage: P = E.an_element()
                sage: S = similarity_surfaces.billiard(P)
                sage: S.minimal_cover("half-translation").stratum()
                Q_1(5, -1^5)
                sage: S.minimal_cover("translation").stratum()
                H_4(6, 0^5)

                sage: E = EuclideanPolygonsWithAngles(1, 3, 5, 7)
                sage: E.billiard_unfolding_stratum("half-translation")
                Q_3(5, 3, 1, -1)
                sage: E.billiard_unfolding_stratum("translation")
                H_7(6, 4, 2)

                sage: P = E.an_element()
                sage: S = similarity_surfaces.billiard(P)
                sage: S.minimal_cover("half-translation").stratum()
                Q_3(5, 3, 1, -1)
                sage: S.minimal_cover("translation").stratum()
                H_7(6, 4, 2, 0)

                sage: E = EuclideanPolygonsWithAngles(1, 2, 8)
                sage: E.billiard_unfolding_stratum("half-translation")
                H_5(7, 1)
                sage: E.billiard_unfolding_stratum("translation")
                H_5(7, 1)

                sage: E.billiard_unfolding_stratum("half-translation", True)
                H_5(7, 1, 0)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_5(7, 1, 0)

                sage: E = EuclideanPolygonsWithAngles(9, 6, 3, 2)
                sage: p = E.an_element()
                sage: B = similarity_surfaces.billiard(p)
                sage: B.minimal_cover("half-translation").stratum()
                Q_4(7, 4, 1, 0)
                sage: E.billiard_unfolding_stratum("half-translation", True)
                Q_4(7, 4, 1, 0)
                sage: B.minimal_cover("translation").stratum()
                H_8(8, 2^3, 0^2)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_8(8, 2^3, 0^2)
            """
            from sage.all import ZZ

            angles = self.billiard_unfolding_angles(cover_type)

            from surface_dynamics import Stratum

            if all(a.is_integer() for a in angles):
                if not marked_points and len(angles) == 1 and 1 in angles:
                    return Stratum([0], 1)
                else:
                    return Stratum(
                        sum(
                            (
                                [ZZ(a - 1)] * mult
                                for a, mult in angles.items()
                                if marked_points or a != 1
                            ),
                            [],
                        ),
                        1,
                    )
            else:
                return Stratum(
                    sum(
                        (
                            [ZZ(2 * (a - 1))] * mult
                            for a, mult in angles.items()
                            if marked_points or a != 1
                        ),
                        [],
                    ),
                    2,
                )

        def billiard_unfolding_stratum_dimension(
            self, cover_type="translation", marked_points=False
        ):
            r"""
            Return the dimension of the stratum of quadratic or Abelian differential
            obtained by unfolding a billiard in a polygon of this equiangular family.

            INPUT:

            - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
              ``"half-translation"`` or ``"translation"``

            - ``marked_poins`` (optional, default ``False``) - whether the stratum should
              have marked regular points

            EXAMPLES::

                sage: from flatsurf import EuclideanPolygonsWithAngles

                sage: E = EuclideanPolygonsWithAngles(1, 1, 1)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                2
                sage: E.billiard_unfolding_stratum_dimension("translation")
                2
                sage: E.billiard_unfolding_stratum_dimension("half-translation", True)
                4
                sage: E.billiard_unfolding_stratum_dimension("translation", True)
                4

                sage: E = EuclideanPolygonsWithAngles(1, 2, 5)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                4
                sage: E.billiard_unfolding_stratum("half-translation").dimension()
                4
                sage: E.billiard_unfolding_stratum_dimension(cover_type="translation")
                6
                sage: E.billiard_unfolding_stratum("translation").dimension()
                6
                sage: E.billiard_unfolding_stratum_dimension("translation", True)
                9
                sage: E.billiard_unfolding_stratum("translation", True).dimension()
                9

                sage: E = EuclideanPolygonsWithAngles(1, 3, 5)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                6
                sage: E.billiard_unfolding_stratum("half-translation").dimension()
                6
                sage: E.billiard_unfolding_stratum_dimension("translation")
                6
                sage: E.billiard_unfolding_stratum("translation").dimension()
                6

                sage: E = EuclideanPolygonsWithAngles(1, 3, 1, 7)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                6

                sage: E = EuclideanPolygonsWithAngles(1, 3, 5, 7)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                8

                sage: E = EuclideanPolygonsWithAngles(1, 2, 8)
                sage: E.billiard_unfolding_stratum_dimension()
                11
                sage: E.billiard_unfolding_stratum().dimension()
                11
                sage: E.billiard_unfolding_stratum_dimension(marked_points=True)
                12
                sage: E.billiard_unfolding_stratum(marked_points=True).dimension()
                12
            """
            from sage.all import ZZ

            if cover_type == "rational":
                raise NotImplementedError
            if cover_type != "translation" and cover_type != "half-translation":
                raise ValueError

            angles = self.billiard_unfolding_angles(cover_type)
            if not marked_points:
                if 1 in angles:
                    del angles[1]
                if not angles:
                    angles[ZZ.one()] = ZZ.one()

            abelian = all(a.is_integer() for a in angles)
            s = sum(angles.values())
            chi = sum(mult * (a - 1) for a, mult in angles.items())
            assert chi.denominator() == 1
            chi = ZZ(chi)
            assert chi % 2 == 0
            g = chi // 2 + 1
            return 2 * g + s - 1 if abelian else 2 * g + s - 2

    class Simple(CategoryWithAxiom_over_base_ring):
        r"""
        The subcategory of simple Euclidean polygons with prescribed angles.

        EXAMPLES::

            sage: from flatsurf.geometry.categories import EuclideanPolygons
            sage: EuclideanPolygons(QQ).WithAngles([1, 1, 1, 1]).Simple()
            Category of simple euclidean rectangles over Rational Field

        """

        # Due to some limitations in SageMath this does not work. Apparently,
        # extra_super_categories cannot depend on parameters of the category.
        # With this code, some polygons are randomly declared as convex.
        # Unfortunately, this means that the category prints as "convex simple
        # rectangles" and not just "rectangles".
        # def extra_super_categories(self):
        #     r"""
        #     Return the categories that simple polygons with prescribed angles
        #     are additionally contained in; namely, in some cases the category
        #     of convex polygons.

        #     EXAMPLES::

        #         sage: from flatsurf.geometry.categories import EuclideanPolygons
        #         sage: C = EuclideanPolygons(QQ).Simple().WithAngles([1, 1, 1, 1])
        #         sage: "Convex" in C.axioms()
        #         True

        #         sage: C = EuclideanPolygons(QQ).Simple().WithAngles([2, 2, 1, 6, 1])
        #         sage: "Convex" in C.axioms()
        #         False

        #     """
        #     # We cannot call is_convex() yet because the SubcategoryMethods
        #     # have not been established yet.
        #     if self._base_category.is_convex():
        #         return (self._base_category.Convex(),)

        #     return ()

        def __call__(self, *lengths, normalized=False, base_ring=None):
            r"""
            Return a polygon with these angles from ``lengths``.

            TESTS::

                sage: from flatsurf import EuclideanPolygonsWithAngles
                sage: P = EuclideanPolygonsWithAngles(1, 2, 1, 2)
                sage: L = P.lengths_polytope()
                sage: r0, r1 = [r.vector() for r in L.rays()]
                sage: lengths = r0 + r1
                sage: P(*lengths[:-2])
                doctest:warning
                ...
                UserWarning: calling EuclideanPolygonsWithAngles() has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon(angles=[...], lengths=[...]) instead.
                To make the resulting polygon non-normalized, i.e., the lengths are not actual edge lengths but the multiple of slope vectors,
                use Polygon(edges=[length * slope for (length, slope) in zip(lengths, EuclideanPolygonsWithAngles(angles).slopes())]).
                Polygon(vertices=[(0, 0), (1, 0), (c + 1, 3), (c, 3)])

                sage: from flatsurf import Polygon, EuclideanPolygonsWithAngles
                sage: P = EuclideanPolygonsWithAngles([1, 2, 1, 2])
                sage: Polygon(angles=[1, 2, 1, 2], lengths=lengths[:-2])
                Polygon(vertices=[(0, 0), (1, 0), (3/2, 1/2*c), (1/2, 1/2*c)])
                sage: Polygon(angles=[1, 2, 1, 2], edges=[length * slope for (length, slope) in zip(lengths[:-2], P.slopes())])
                Polygon(vertices=[(0, 0), (1, 0), (c + 1, 3), (c, 3)])

                sage: P = EuclideanPolygonsWithAngles(2, 2, 3, 13)
                sage: r0, r1 = [r.vector() for r in P.lengths_polytope().rays()]
                sage: P(r0 + r1)
                Polygon(vertices=[(0, 0), (20, 0), (5, -15*c^3 + 60*c), (5, -5*c^3 + 20*c)])

                sage: P = EuclideanPolygonsWithAngles([2, 2, 3, 13])
                sage: Polygon(angles=[2, 2, 3, 13], lengths=r0 + r1)
                Traceback (most recent call last):
                ...
                ValueError: polygon not closed
                sage: Polygon(angles=[2, 2, 3, 13], edges=[length * slope for (length, slope) in zip(r0 + r1, P.slopes())])
                Polygon(vertices=[(0, 0), (20, 0), (5, -15*c^3 + 60*c), (5, -5*c^3 + 20*c)])

            """
            # __call__() cannot be properly inherited in subcategories since it
            # cannot be in SubcategoryMethods; that's why we want to get rid of it.
            import warnings

            warning = "calling EuclideanPolygonsWithAngles() has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon(angles=[...], lengths=[...]) instead."

            if not normalized:
                warning += (
                    " To make the resulting polygon non-normalized, i.e., the lengths are not actual edge lengths but the multiple of slope vectors, use "
                    "Polygon(edges=[length * slope for (length, slope) in zip(lengths, EuclideanPolygonsWithAngles(angles).slopes())])."
                )

            warnings.warn(warning)

            from sage.structure.element import Vector

            if len(lengths) == 1 and isinstance(lengths[0], (tuple, list, Vector)):
                lengths = lengths[0]

            n = len(self.angles())
            if len(lengths) != n - 2 and len(lengths) != n:
                raise ValueError(
                    f"must provide {n - 2} or {n} lengths but provided {len(lengths)}"
                )

            V = self.base_ring() ** 2
            slopes = self.slopes()
            if normalized:
                cosines_ring = (
                    self._without_axiom("Simple")
                    ._without_axiom("Convex")
                    ._cosines_ring()
                )
                V = V.change_ring(cosines_ring)
                for i, s in enumerate(slopes):
                    x, y = map(cosines_ring, s)
                    norm2 = (x**2 + y**2).sqrt()
                    slopes[i] = V((x / norm2, y / norm2))

            if base_ring is None:
                from sage.all import Sequence

                base_ring = Sequence(lengths).universe()

                from sage.categories.pushout import pushout

                if normalized:
                    base_ring = pushout(base_ring, cosines_ring)
                else:
                    base_ring = pushout(base_ring, self.base_ring())

            v = V((0, 0))
            vertices = [v]

            from sage.all import vector, matrix

            if len(lengths) == n - 2:
                for i in range(n - 2):
                    v += lengths[i] * slopes[i]
                    vertices.append(v)
                s, t = (
                    vector(vertices[0] - vertices[n - 2])
                    * matrix([slopes[-1], slopes[n - 2]]).inverse()
                )
                assert (
                    vertices[0] - s * slopes[-1] == vertices[n - 2] + t * slopes[n - 2]
                )
                if s <= 0 or t <= 0:
                    raise ValueError(
                        "the provided lengths do not give rise to a polygon"
                    )
                vertices.append(vertices[0] - s * slopes[-1])

            elif len(lengths) == n:
                for i in range(n):
                    v += lengths[i] * slopes[i]
                    vertices.append(v)
                if not vertices[-1].is_zero():
                    raise ValueError(
                        "the provided lengths do not give rise to a polygon"
                    )
                vertices.pop(-1)

            category = self.change_ring(base_ring)
            if self.is_convex():
                category = category.Convex()

            from flatsurf.geometry.polygon import Polygon

            return Polygon(base_ring=base_ring, vertices=vertices, category=category)


@cached_function
def _slopes(angles):
    r"""
    Return the slopes of the sides of a polygon with ``angles`` in a
    (possibly non-minimal) number field.

    .. NOTE::

        This function gets called a lot from the other functions here. We
        could refactor things and make sure that the function is only
        invoked once. However, it seems easier to just cache its output.
        Also, this speeds up the relative common case when multiple
        polygons with the same angles are created. (If lots of polygons
        with different angles get created then we pay this cache with RAM
        of course but in our experiments it has not been an issue.)

    .. NOTE::

        SageMath 9.1 does not support cached functions that are staticmethods.
        Therefore we define this function at the module level.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.euclidean_polygons_with_angles import _slopes
        sage: _slopes((1/6, 1/6, 1/6))
        [(c, 3), (c, 3), (c, 3)]
        sage: _slopes((1/4, 1/4, 1/4, 1/4))
        [(0, 1), (0, 1), (0, 1), (0, 1)]
        sage: _slopes((1/10, 2/10, 3/10, 4/10))
        [(c^3, 5), (3*c^3 - 10*c, 5), (-3*c^3 + 10*c, 5), (-c^3, 5)]

    """
    from sage.all import QQ, RIF, lcm, AA, NumberField
    from flatsurf.geometry.subfield import chebyshev_T, cos_minpoly

    # We determine the number field that contains the slopes of the sides,
    # i.e., the cosines and sines of the inner angles of the polygon.
    # Let us first write all angles as multiples of 2π/N with the smallest
    # possible common N.
    N = lcm(a.denominator() for a in angles)
    # The field containing the cosine and sine of 2π/N might be too small
    # to write down all the slopes when N is not divisible by 4.
    if N == 1:
        raise ValueError("there cannot be a polygon with all angles multiples of 2π")
    if N == 2:
        pass
    elif N % 4:
        while N % 4:
            N *= 2

    angles = [QQ(a * N) for a in angles]

    if N == 2:
        base_ring = QQ
        c = QQ.zero()
    else:
        # Construct the minimal polynomial f(x) of c = 2 cos(2π / N)
        f = cos_minpoly(N // 2)
        emb = AA.polynomial_root(f, 2 * (2 * RIF.pi() / N).cos())
        base_ring = NumberField(f, "c", embedding=emb)
        c = base_ring.gen()

    # Construct the cosine and sine of each angle as an element of our number field.
    def cosine(a):
        return chebyshev_T(abs(a), c) / 2

    def sine(a):
        # Use sin(x) = cos(π/2 - x)
        return cosine(N // 4 - a)

    slopes = [(cosine(a), sine(a)) for a in angles]

    assert all((x**2 + y**2).is_one() for x, y in slopes)

    from flatsurf.geometry.euclidean import projectivization

    return [projectivization(x, y) for x, y in slopes]


@cached_function
def _base_ring(angles):
    r"""
    Return a minimal number field containing all the :meth:`_slopes` of a
    polygon with ``angles``.

    .. NOTE::

        Internally, this uses
        :func:`~flatsurf.geometry.subfield.subfield_from_element` which is
        very slow. We therefore cache the result currently.

    .. NOTE::

        SageMath 9.1 does not support cached functions that are staticmethods.
        Therefore we define this function at the module level.

    EXAMPLES::

        sage: from flatsurf.geometry.categories.euclidean_polygons_with_angles import _base_ring
        sage: _base_ring((1/6, 1/6, 1/6))
        Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?
        sage: _base_ring((1/4, 1/4, 1/4, 1/4))
        Rational Field
        sage: _base_ring((1/10, 2/10, 3/10, 4/10))
        Number Field in c with defining polynomial x^4 - 5*x^2 + 5 with c = 1.902113032590308?

    """
    slopes = _slopes(angles)
    base_ring = slopes[0][0].parent()

    # It might be the case that the slopes generate a smaller field. For
    # now we use an ugly workaround via subfield_from_elements.
    old_slopes = []
    for v in slopes:
        old_slopes.extend(v)
    from flatsurf.geometry.subfield import subfield_from_elements

    L, _, _ = subfield_from_elements(base_ring, old_slopes)
    return L
