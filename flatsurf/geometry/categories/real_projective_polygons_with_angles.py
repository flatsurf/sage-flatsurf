from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from flatsurf.geometry.categories.real_projective_polygons import RealProjectivePolygons


class RealProjectivePolygonsWithAngles(Category_over_base_ring):
    # TODO: This is the docstring of EquiangularPolygons
    r"""
    Polygons with fixed (rational) angles.

    EXAMPLES::

        sage: from flatsurf import EquiangularPolygons

    The polygons with inner angles `\pi/4`, `\pi/2`, `5\pi/4`::

        sage: P = EquiangularPolygons(1, 2, 5)
        sage: P
        Category of real projective polygons with angles (1, 2, 5) over Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    Internally, polygons are given by their vertices' coordinates over some
    number field, in this case a quadratic field::

        sage: P.base_ring()
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    Polygons can also be defined over other number field implementations::

        sage: from pyeantic import RealEmbeddedNumberField # optional: eantic  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: K = RealEmbeddedNumberField(P.base_ring()) # optional: eantic
        sage: P(K(1)) # optional: eantic
        polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])
        sage: _.base_ring() # optional: eantic
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    However, specific instances of such polygons might be defined over another ring::

        sage: P(1)
        polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])
        sage: _.base_ring()
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

        sage: P(AA(1))
        polygon(vertices=[(0, 0), (1, 0), (0.7071067811865475?, 0.2928932188134525?)])
        sage: _.base_ring()
        Algebraic Real Field

    Polygons can also be defined over a module containing transcendent parameters::

        sage: from pyexactreal import ExactReals # optional: exactreal  # random output due to deprecation warnings with some versions of pkg_resources
        sage: R = ExactReals(P.base_ring()) # optional: exactreal
        sage: P(R(1)) # optional: exactreal
        polygon(vertices=[(0, 0), (1, 0), ((1/2*c0 ~ 0.70710678), (-1/2*c0+1 ~ 0.29289322))])
        sage: P(R(R.random_element([0.2, 0.3]))) # random output, optional: exactreal
        polygon(vertices=[(0, 0),])
                 (ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54), 0),
                 (((12*c0+17 ~ 33.970563)*ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54))/((17*c0+24 ~ 48.041631)),
                 ((5*c0+7 ~ 14.071068)*ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54))/((17*c0+24 ~ 48.041631)))
        sage: _.base_ring() # optional: exactreal
        Real Numbers as (Real Embedded Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?)-Module

    ::

        sage: L = P.lengths_polytope()    # polytope of admissible lengths for edges
        sage: L
        A 1-dimensional polyhedron in (Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?)^3 defined as the convex hull of 1 vertex and 1 ray
        sage: lengths = L.rays()[0].vector()
        sage: lengths
        (1, -1/2*c0 + 1, -1/2*c0 + 1)
        sage: p = P(*lengths)    # build one polygon with the given lengths
        sage: p
        polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])
        sage: p.angles()
        (1/16, 1/8, 5/16)
        sage: P.angles(integral=False)
        (1/16, 1/8, 5/16)

        sage: P = EquiangularPolygons(1, 2, 1, 2, 2, 1)
        sage: L = P.lengths_polytope()
        sage: L
        A 4-dimensional polyhedron in (Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?)^6 defined as the convex hull of 1 vertex and 6 rays
        sage: rays = [r.vector() for r in L.rays()]
        sage: rays
        [(1, 0, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (0, 1, 0, 0, c^2 - 3, c^2 - 2),
         (1/3*c^4 - 2*c^2 + 3, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (-c^4 + 4*c^2, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (0, 1/3*c^4 - 2*c^2 + 3, c^2 - 3, 0, 0, 1/3*c^4 - c^2),
         (0, -c^4 + 4*c^2, 0, c^2 - 3, 0, -c^4 + 5*c^2 - 3)]
        sage: lengths = 3*rays[0] + rays[2] + 2*rays[3] + rays[4]
        sage: p = P(*lengths)
        sage: p
        polygon(vertices=[(0, 0),
                          (-5/3*c^4 + 6*c^2 + 6, 0),
                          (3*c^5 - 5/3*c^4 - 16*c^3 + 6*c^2 + 18*c + 6, c^4 - 6*c^2 + 9),
                          (2*c^5 - 2*c^4 - 10*c^3 + 15/2*c^2 + 9*c + 5, -1/2*c^5 + c^4 + 5/2*c^3 - 3*c^2 - 2*c),
                          (2*c^5 - 10*c^3 - 3/2*c^2 + 9*c + 9, -3/2*c^5 + c^4 + 15/2*c^3 - 3*c^2 - 6*c),
                          (2*c^5 - 10*c^3 - 3*c^2 + 9*c + 12, -3*c^5 + c^4 + 15*c^3 - 3*c^2 - 12*c)])

        sage: p.angles()
        (2/9, 4/9, 2/9, 4/9, 4/9, 2/9)

        sage: EquiangularPolygons(1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2, 1)
        Category of real projective polygons with angles (1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2, 1) over Number Field in c with defining polynomial x^22 - 23*x^20 + 230*x^18 - 1311*x^16 + 4692*x^14 - 10948*x^12 + 16744*x^10 - 16445*x^8 + 9867*x^6 - 3289*x^4 + 506*x^2 - 23 with c = 1.995337538381079?

    A regular pentagon::

        sage: E = EquiangularPolygons(1, 1, 1, 1, 1)
        sage: E(1, 1, 1, 1, 1, normalized=True)
        polygon(vertices=[(0, 0), (1, 0), (1/2*c^2 - 1/2, 1/2*c), (1/2, 1/2*c^3 - c), (-1/2*c^2 + 3/2, 1/2*c)])
    """

    def super_categories(self):
        return [RealProjectivePolygons(self.base_ring()).Rational()]

    @staticmethod
    def _normalize_angles(angles):
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
            raise NotImplementedError("each angle must be > 0 and < 2 pi")

        angles = tuple(angles)

        return angles

    @staticmethod
    def _slopes(angles):
        from sage.all import ZZ, QQ, RIF, lcm, AA, NumberField
        from flatsurf.geometry.subfield import chebyshev_T, cos_minpoly

        # We determine the number field that contains the slopes of the sides,
        # i.e., the cosines and sines of the inner angles of the polygon.
        # Let us first write all angles as multiples of 2π/N with the smallest
        # possible common N.
        N = lcm(a.denominator() for a in angles)
        # The field containing the cosine and sine of 2π/N might be too small
        # to write down all the slopes when N is not divisible by 4.
        assert N != 1, "there cannot be a polygon with all angles multiples of 2π"
        if N == 2:
            pass
        elif N % 4:
            while N % 4:
                N *= 2

        angles = [ZZ(a * N) for a in angles]

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

        from flatsurf.geometry.polygon import projectivization
        return [projectivization(x, y) for x, y in slopes]

    @staticmethod
    def _base_ring(angles):
        slopes = RealProjectivePolygonsWithAngles._slopes(angles)
        base_ring = slopes[0][0].parent()

        # TODO: It might be the case that the slopes generate a smaller
        # field. For now we use an ugly workaround via subfield_from_elements.
        old_slopes = []
        for v in slopes:
            old_slopes.extend(v)
        from flatsurf.geometry.subfield import subfield_from_elements
        L, _, _ = subfield_from_elements(base_ring, old_slopes)
        return L

    def __init__(self, base_ring, angles):
        self._angles = angles
        self._slopes = RealProjectivePolygonsWithAngles._slopes(angles)
        self._cosines_ring = self._slopes[0][0].parent()
        self._slopes = [(base_ring(slope[0]), base_ring(slope[1])) for slope in self._slopes]

        super().__init__(base_ring)

    def _repr_object_names(self):
        names = super()._repr_object_names()
        return names.replace("with angles", f"with angles {self.angles(True)}")

    def __call__(self, *lengths, normalized=False, base_ring=None):
        r"""
        TESTS::

            sage: from flatsurf import EquiangularPolygons
            sage: P = EquiangularPolygons(1, 2, 1, 2)
            sage: L = P.lengths_polytope()
            sage: r0, r1 = [r.vector() for r in L.rays()]
            sage: lengths = r0 + r1
            sage: P(*lengths[:-2])
            polygon(vertices=[(0, 0), (1, 0), (c + 1, 3), (c, 3)])

            sage: P = EquiangularPolygons(2, 2, 3, 13)
            sage: r0, r1 = [r.vector() for r in P.lengths_polytope().rays()]
            sage: P(r0 + r1)
            polygon(vertices=[(0, 0), (20, 0), (5, -15*c^3 + 60*c), (5, -5*c^3 + 20*c)])
        """
        from sage.structure.element import Vector
        if len(lengths) == 1 and isinstance(lengths[0], (tuple, list, Vector)):
            lengths = lengths[0]

        n = len(self._angles)
        if len(lengths) != n - 2 and len(lengths) != n:
            raise ValueError(
                "must provide %d or %d lengths but provided %d"
                % (n - 2, n, len(lengths))
            )

        V = self.module()
        slopes = self.slopes()
        if normalized:
            V = V.change_ring(self._cosines_ring)
            for i, s in enumerate(slopes):
                x, y = map(self._cosines_ring, s)
                norm2 = (x**2 + y**2).sqrt()
                slopes[i] = V((x / norm2, y / norm2))

        if base_ring is None:
            from sage.all import Sequence

            base_ring = Sequence(lengths).universe()

            from sage.categories.pushout import pushout

            if normalized:
                base_ring = pushout(base_ring, self._cosines_ring)
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
            assert vertices[0] - s * slopes[-1] == vertices[n - 2] + t * slopes[n - 2]
            if s <= 0 or t <= 0:
                raise ValueError("the provided lengths do not give rise to a polygon")
            vertices.append(vertices[0] - s * slopes[-1])

        elif len(lengths) == n:
            for i in range(n):
                v += lengths[i] * slopes[i]
                vertices.append(v)
            if not vertices[-1].is_zero():
                raise ValueError("the provided lengths do not give rise to a polygon")
            vertices.pop(-1)

        category = RealProjectivePolygons(base_ring)
        if self.convexity():
            category = category.Convex()

        from flatsurf.geometry.polygon import EuclideanPolygon
        return EuclideanPolygon(base_ring, vertices=vertices, category=category)

    class SubcategoryMethods:
        def convexity(self):
            r"""
            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: EquiangularPolygons(1, 2, 5).convexity()
                True
                sage: EquiangularPolygons(2, 2, 3, 13).convexity()
                False
            """
            return all(2 * a <= 1 for a in self._angles)

        # def base_ring(self):
        #     r"""
        #     Return the number field over which the coordinates of the vertices of
        #     this family of polygons are represented internally.

        #     EXAMPLES::

        #         sage: from flatsurf import EquiangularPolygons
        #         sage: EquiangularPolygons(1, 2, 5).base_ring()
        #         Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

        #     """
        #     return self._base_ring

        def angles(self, integral=False):
            r"""
            Return the interior angles of this polygon as multiples 2π.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: E = EquiangularPolygons(1, 1, 1, 2, 6)
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
                C = lcm([a.denominator() for a in self._angles]) / gcd(
                    [a.numerator() for a in self._angles]
                )
                angles = tuple(ZZ(C * a) for a in angles)
            return angles

        @cached_method
        def __angles(self):
            if isinstance(self, RealProjectivePolygonsWithAngles):
                return self._angles

            for category in self.super_categories():
                if isinstance(category, RealProjectivePolygonsWithAngles):
                    return category._angles

            raise NotImplementedError

        def strict_convexity(self):
            r"""
            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: E = EquiangularPolygons([1, 1, 1, 1, 2])
                sage: E.angles()
                (1/4, 1/4, 1/4, 1/4, 1/2)
                sage: E.convexity()
                True
                sage: E.strict_convexity()
                False

            """
            return all(2 * a < 1 for a in self._angles)

        @cached_method
        def module(self):
            r"""
            Return the free module of rank 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: C = EquiangularPolygons(1, 2, 3)
                sage: C.module()
                Vector space of dimension 2 over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

            """
            from sage.all import FreeModule
            return FreeModule(self.base_ring(), 2)

        @cached_method
        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: C = EquiangularPolygons(1, 2, 3)
                sage: C.vector_space()
                Vector space of dimension 2 over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

            """
            from sage.all import VectorSpace
            return VectorSpace(self.base_ring().fraction_field(), 2)

        def slopes(self, e0=(1, 0)):
            r"""
            List of slopes of the edges as a list of vectors.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: EquiangularPolygons(1, 2, 1, 2).slopes()
                [(1, 0), (c, 3), (-1, 0), (-c, -3)]
            """
            V = self.module()
            slopes = self._slopes
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
                from flatsurf.geometry.polygon import projectivization
                e = projectivization(*e)
                edges.append(V(e))
            return edges

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

                sage: from flatsurf import EquiangularPolygons
                sage: EquiangularPolygons(1, 2, 1, 2).lengths_polytope()
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
            Return a polygon in this family.

            Note that this might fail due to intersection.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: EquiangularPolygons(4, 3, 4, 4, 3, 4).an_element()
                polygon(vertices=[(0, 0),
                                  (1/22*c + 1, 0),
                                  (9*c^9 + 1/2*c^8 - 88*c^7 - 9/2*c^6 + 297*c^5 + 27/2*c^4 - 396*c^3 - 15*c^2 + 3631/22*c + 11/2, 1/2*c + 11),
                                  (16*c^9 + c^8 - 154*c^7 - 9*c^6 + 506*c^5 + 27*c^4 - 638*c^3 - 30*c^2 + 4841/22*c + 9, c + 22),
                                  (16*c^9 + c^8 - 154*c^7 - 9*c^6 + 506*c^5 + 27*c^4 - 638*c^3 - 30*c^2 + 220*c + 8, c + 22),
                                  (7*c^9 + 1/2*c^8 - 66*c^7 - 9/2*c^6 + 209*c^5 + 27/2*c^4 - 242*c^3 - 15*c^2 + 55*c + 7/2, 1/2*c + 11)])
            """
            return self(sum(r.vector() for r in self.lengths_polytope().rays()))

        def random_element(self, ring=None, **kwds):
            r"""
            Return a random polygon.

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons
                sage: EquiangularPolygons(1, 1, 1, 2, 5).random_element()
                polygon(vertices=[(0, 0), ...])
                sage: EquiangularPolygons(1,1,1,15,15,15).random_element()
                polygon(vertices=[(0, 0), ...])
                sage: EquiangularPolygons(1,15,1,15,1,15).random_element()
                polygon(vertices=[(0, 0), ...])
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
                coeffs, r = random_element()
                try:
                    return self(*r)
                except ValueError:
                    pass

        def billiard_unfolding_angles(self, cover_type="translation"):
            r"""
            Return the angles of the unfolding rational, half-translation or translation surface.

            INPUT:

            - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
              ``"half-translation"`` or ``"translation"``

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons

                sage: E = EquiangularPolygons(1, 2, 5)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/8: 1, 1/4: 1, 5/8: 1}
                sage: (1/8 - 1) + (1/4 - 1) + (5/8 - 1)  # Euler characteristic (of the sphere)
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 3, 5/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 3, 5: 1}

                sage: E = EquiangularPolygons(1, 3, 1, 7)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/6: 2, 1/2: 1, 7/6: 1}
                sage: 2 * (1/6 - 1) + (1/2 - 1) + (7/6 - 1) # Euler characteristic
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 5, 7/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 5, 7: 1}

                sage: E = EquiangularPolygons(1, 3, 5, 7)
                sage: E.billiard_unfolding_angles(cover_type="rational")
                {1/8: 1, 3/8: 1, 5/8: 1, 7/8: 1}
                sage: (1/8 - 1) + (3/8 - 1) + (5/8 - 1) + (7/8 - 1) # Euler characteristic
                -2
                sage: E.billiard_unfolding_angles(cover_type="half-translation")
                {1/2: 1, 3/2: 1, 5/2: 1, 7/2: 1}
                sage: E.billiard_unfolding_angles(cover_type="translation")
                {1: 1, 3: 1, 5: 1, 7: 1}

                sage: E = EquiangularPolygons(1, 2, 8)
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

        def billiard_unfolding_stratum(self, cover_type="translation", marked_points=False):
            r"""
            Return the stratum of quadratic or Abelian differential obtained by
            unfolding a billiard in a polygon of this equiangular family.

            INPUT:

            - ``cover_type`` (optional, default ``"translation"``) - either ``"rational"``,
              ``"half-translation"`` or ``"translation"``

            - ``marked_poins`` (optional, default ``False``) - whether the stratum should
              have regular marked points

            EXAMPLES::

                sage: from flatsurf import EquiangularPolygons, similarity_surfaces

                sage: E = EquiangularPolygons(1, 2, 5)
                sage: E.billiard_unfolding_stratum("half-translation")
                Q_1(3, -1^3)
                sage: E.billiard_unfolding_stratum("translation")
                H_3(4)
                sage: E.billiard_unfolding_stratum("half-translation", True)
                Q_1(3, -1^3)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_3(4, 0^3)

                sage: E = EquiangularPolygons(1, 3, 1, 7)
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

                sage: E = EquiangularPolygons(1, 3, 5, 7)
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

                sage: E = EquiangularPolygons(1, 2, 8)
                sage: E.billiard_unfolding_stratum("half-translation")
                H_5(7, 1)
                sage: E.billiard_unfolding_stratum("translation")
                H_5(7, 1)

                sage: E.billiard_unfolding_stratum("half-translation", True)
                H_5(7, 1, 0)
                sage: E.billiard_unfolding_stratum("translation", True)
                H_5(7, 1, 0)

                sage: E = EquiangularPolygons(9, 6, 3, 2)
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
            if all(a.is_integer() for a in angles):
                from surface_dynamics import AbelianStratum

                if not marked_points and len(angles) == 1 and 1 in angles:
                    return AbelianStratum([0])
                else:
                    return AbelianStratum(
                        {
                            ZZ(a - 1): mult
                            for a, mult in angles.items()
                            if marked_points or a != 1
                        }
                    )
            else:
                from surface_dynamics import QuadraticStratum

                return QuadraticStratum(
                    {
                        ZZ(2 * (a - 1)): mult
                        for a, mult in angles.items()
                        if marked_points or a != 1
                    }
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

                sage: from flatsurf import EquiangularPolygons

                sage: E = EquiangularPolygons(1, 1, 1)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                2
                sage: E.billiard_unfolding_stratum_dimension("translation")
                2
                sage: E.billiard_unfolding_stratum_dimension("half-translation", True)
                4
                sage: E.billiard_unfolding_stratum_dimension("translation", True)
                4

                sage: E = EquiangularPolygons(1, 2, 5)
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

                sage: E = EquiangularPolygons(1, 3, 5)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                6
                sage: E.billiard_unfolding_stratum("half-translation").dimension()
                6
                sage: E.billiard_unfolding_stratum_dimension("translation")
                6
                sage: E.billiard_unfolding_stratum("translation").dimension()
                6

                sage: E = EquiangularPolygons(1, 3, 1, 7)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                6

                sage: E = EquiangularPolygons(1, 3, 5, 7)
                sage: E.billiard_unfolding_stratum_dimension("half-translation")
                8

                sage: E = EquiangularPolygons(1, 2, 8)
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
