r"""
TODO: Document this module.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Julian Rüth
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
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.misc.cachefunc import cached_method
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.all import ZZ
from dataclasses import dataclass


class HarmonicDifferential(Element):
    def __init__(self, parent, series):
        super().__init__(parent)
        self._series = series

    def _add_(self, other):
        r"""
        Return the sum of this harmonic differential and ``other`` by summing
        their underlying power series.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f); η  # tol 1e-6
            (0 - 1.*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
             0 - 1.*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

            sage: η + η  # tol 1e-6
            (0 - 2*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
             0 - 2*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

        """
        return self.parent()({
            triangle: self._series[triangle] + other._series[triangle]
            for triangle in self._series
        })

    def evaluate(self, triangle, Δ, derivative=0):
        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        return self._evaluate(C.evaluate(triangle, Δ, derivative=derivative))

    def _evaluate(self, expression):
        r"""
        Evaluate an expression by plugging in the coefficients of the power
        series defining this differential.

        This might not correspond to evaluating the actual power series somewhere.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

        Compute the sum of the constant coefficients::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 5)
            sage: R = C.symbolic_ring()
            sage: η._evaluate(R(C.gen(0, 0)) + R(C.gen(1, 0)))  # tol 1e-6
            0 - 2*I

        """
        coefficients = {}

        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), self.parent()._geometry)

        for gen in expression.variables():
            kind, triangle, k = C._describe_generator(gen)
            coefficient = self._series[triangle][k]

            if kind == "real":
                coefficients[gen] = coefficient.real()
            else:
                assert kind == "imag"
                coefficients[gen] = coefficient.imag()

        value = expression.parent()(C._subs(expression, coefficients))
        assert value.degree() <= 0
        return value.constant_coefficient()

    @cached_method
    def precision(self):
        precisions = set(series.precision_absolute() for series in self._series.values())
        assert len(precisions) == 1
        return next(iter(precisions))

    def integrate(self, cycle):
        r"""
        Return the integral of this differential along the homology class
        ``cycle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)

        Construct a differential form such that integrating it along `a` yields real part
        1, and along `b` real part 0::

            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: η.integrate(a).real()  # tol 1e-6
            1

            sage: η.integrate(b).real()  # tol 1e-6
            0

        """
        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), self.parent()._geometry)
        return self._evaluate(C.integrate(cycle))

    def _repr_(self):
        return repr(tuple(self._series.values()))


class HarmonicDifferentials(UniqueRepresentation, Parent):
    r"""
    The space of harmonic differentials on this surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
        sage: T.set_immutable()

        sage: Ω = HarmonicDifferentials(T); Ω
        Ω(TranslationSurface built from 2 polygons)

    ::

        sage: H = SimplicialCohomology(T)
        sage: Ω(H())
        (O(z0^10), O(z1^10))

    ::

        sage: a, b = H.homology().gens()
        sage: f = H({b: 1})
        sage: η = Ω(f)
        sage: η.integrate(a).real()  # tol 1e-6
        0
        sage: η.integrate(b).real()  # tol 1e-6
        1

    """
    Element = HarmonicDifferential

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, category=None):
        r"""
        Normalize parameters when creating the space of harmonic differentials.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: HarmonicDifferentials(T) is HarmonicDifferentials(T)
            True

        """
        from sage.all import RR
        return super().__classcall__(cls, surface, coefficients or RR, category or SetsWithPartialMaps())

    def __init__(self, surface, coefficients, category):
        if surface != surface.delaunay_triangulation():
            raise NotImplementedError("Surface must be Delaunay triangulated")

        Parent.__init__(self, category=category, base=coefficients)

        self._surface = surface
        # TODO: Coefficients must be reals of some sort?
        self._coefficients = coefficients

        self._geometry = GeometricPrimitives(surface)

    def surface(self):
        return self._surface

    @cached_method
    def power_series_ring(self, triangle):
        r"""
        Return the power series ring to write down the series describing a
        harmonic differential in a Voronoi cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: Ω = HarmonicDifferentials(T)
            sage: Ω.power_series_ring(1)
            Power Series Ring in z1 over Complex Field with 53 bits of precision
            sage: Ω.power_series_ring(2)
            Power Series Ring in z2 over Complex Field with 53 bits of precision

        """
        from sage.all import PowerSeriesRing
        # TODO: Should we use self._coefficients in some way?
        from sage.all import CC
        return PowerSeriesRing(CC, f"z{triangle}")

    def _repr_(self):
        return f"Ω({self._surface})"

    def _element_constructor_(self, x, *args, **kwargs):
        from flatsurf.geometry.cohomology import SimplicialCohomology
        cohomology = SimplicialCohomology(self._surface, self._coefficients)

        if not x:
            return self._element_from_cohomology(cohomology(), *args, **kwargs)

        if isinstance(x, dict):
            return self.element_class(self, x, *args, **kwargs)

        if x.parent() is cohomology:
            return self._element_from_cohomology(x, *args, **kwargs)

        raise NotImplementedError()

    def _element_from_cohomology(self, cocycle, /, prec=10, algorithm=["midpoint_derivatives", "area_upper_bound"], check=True):
        # TODO: In practice we could speed things up a lot with some smarter
        # caching. A lot of the quantities used in the computations only depend
        # on the surface & precision. When computing things for many cocycles
        # we do not need to recompute them. (But currently, we probably do
        # because they live in the constraints instance.)

        # We develop a consistent system of Laurent series at each vertex of the Voronoi diagram
        # to describe a differential.

        # Let η be the differential we are looking for. To describe η we will use ω, the differential
        # corresponding to the flat structure given by our triangulation on this Riemann surface.
        # Then f:=η/ω is a meromorphic function which we can develop locally into a Laurent series.
        # Away from the vertices of the triangulation, ω has no zeros, so f has no poles there and is
        # thus given by a power series.

        # At each vertex of the Voronoi diagram, write f=Σ a_k z^k + O(z^prec). Our task is now to determine
        # the a_k.

        constraints = PowerSeriesConstraints(self.surface(), prec=prec, geometry=self._geometry)

        # We use a variety of constraints. Which ones to use exactly is
        # determined by the "algorithm" parameter. If algorithm is a dict, it
        # can be used to configure aspects of the constraints.
        def get_parameter(alg, default):
            assert alg in algorithm
            if isinstance(algorithm, dict):
                return algorithm[alg]
            return default

        # (0) If two power series are developed around essentially the same
        # point (because the Delaunay triangulation is ambiguous) we force them
        # to coincide.
        constraints.require_equality()

        # (1) The radius of convergence of the power series is the distance from the vertex of the Voronoi
        # cell to the closest vertex of the triangulation (since we use a Delaunay triangulation, all vertices
        # are at the same distance in fact.) So the radii of convergence of two neigbhouring cells overlap
        # and the power series must coincide there. Note that this constraint is unrelated to the cohomology
        # class Φ.
        if "midpoint_derivatives" in algorithm:
            derivatives = get_parameter("midpoint_derivatives", prec//3)
            constraints.require_midpoint_derivatives(derivatives)

        # (1') TODO: Describe L2 optimization.
        if "L2" in algorithm:
            weight = get_parameter("L2", 1)
            constraints.optimize(weight * constraints._L2_consistency())

        # (2) We have that for any cycle γ, Re(∫fω) = Re(∫η) = Φ(γ). We can turn this into constraints
        # on the coefficients as we integrate numerically following the path γ as it intersects the radii of
        # convergence.
        constraints.require_cohomology(cocycle)

        # (3) Since the area ∫ η \wedge \overline{η} must be finite [TODO:
        # REFERENCE?] we optimize for a proxy of this quantity to be minimal.
        if "area_upper_bound" in algorithm:
            weight = get_parameter("area_upper_bound", 1)
            constraints.optimize(weight * constraints._area_upper_bound())

        # (3') We can also optimize for the exact quantity to be minimal but
        # this is much slower.
        if "area" in algorithm:
            weight = get_parameter("area", 1)
            constraints.optimize(weight * constraints._area())

        η = self.element_class(self, constraints.solve())

        # TODO: Factor this out so we can use it in reporting.
        if check:
            # Check whether this is actually a global differential:
            # (1) Check that the series are actually consistent where the Voronoi cells overlap.
            def check(actual, expected, message, abs_error_bound = 1e-9, rel_error_bound = 1e-6):
                abs_error = abs(expected - actual)
                if abs_error > abs_error_bound:
                    if expected == 0 or abs_error / abs(expected) > rel_error_bound:
                        print(f"{message}; expected: {expected}, got: {actual}")

            for (triangle, edge) in self._surface.edge_iterator():
                triangle_, edge_ = self._surface.opposite_edge(triangle, edge)
                for derivative in range(prec//3):
                    expected = η.evaluate(triangle, self._geometry.midpoint(triangle, edge), derivative)
                    other = η.evaluate(triangle_, self._geometry.midpoint(triangle_, edge_), derivative)
                    check(other, expected, f"power series defining harmonic differential are not consistent: {derivative}th derivate does not match between {(triangle, edge)} and {(triangle_, edge_)}")

            # (2) Check that differential actually integrates like the cohomology class.
            for gen in cocycle.parent().homology().gens():
                expected = cocycle(gen)
                actual = η.integrate(gen).real
                if callable(actual):
                    actual = actual()
                check(actual, expected, f"harmonic differential does not have prescribed integral at {gen}")

            # (3) Check that the area is finite.

        return η


class GeometricPrimitives:
    def __init__(self, surface):
        # TODO: Require immutable.
        self._surface = surface

    @cached_method
    def midpoint(self, triangle, edge):
        r"""
        Return the vector to go from the center of the circumcircle of
        ``triangle`` to the midpoint of ``edge``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import GeometricPrimitives
            sage: G = GeometricPrimitives(T)

            sage: G.midpoint(0, 0)
            0j
            sage: G.midpoint(0, 1)
            0.5j
            sage: G.midpoint(0, 2)
            (-0.5+0j)
            sage: G.midpoint(1, 0)
            0j
            sage: G.midpoint(1, 1)
            -0.5j
            sage: G.midpoint(1, 2)
            (0.5+0j)

        """
        polygon = self._surface.polygon(triangle)
        return -self.center(triangle) + complex(*polygon.vertex(edge)) + complex(*polygon.edge(edge)) / 2

    @cached_method
    def center(self, triangle):
        return complex(*self._surface.polygon(triangle).circumscribing_circle().center())


class PowerSeriesConstraints:
    r"""
    A collection of (linear) constraints on the coefficients of power series
    developed at the vertices of the Voronoi cells of a Delaunay triangulation.

    This is used to create harmonic differentials from cohomology classes.
    """
    @dataclass
    class Constraint:
        real: dict
        imag: dict
        lagrange: list
        value: complex

    def __init__(self, surface, prec, geometry=None):
        self._surface = surface
        self._prec = prec
        self._geometry = geometry or GeometricPrimitives(surface)
        self._constraints = []
        self._cost = self.symbolic_ring().zero()

    def __repr__(self):
        return repr(self._constraints)

    @cached_method
    def symbolic_ring(self, *triangles):
        r"""
        Return the polynomial ring in the coefficients of the power series at
        ``triangle``.

        If ``triangle`` is not set, return the polynomial ring with the
        coefficients for all the triangles.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.symbolic_ring(1)
            Multivariate Polynomial Ring in Re_a1_0, Re_a1_1, Re_a1_2, Im_a1_0, Im_a1_1, Im_a1_2 over Complex Field with 53 bits of precision

            sage: C.symbolic_ring()
            Multivariate Polynomial Ring in Re_a0_0, Re_a0_1, Re_a0_2, Im_a0_0, Im_a0_1, Im_a0_2, Re_a1_0, Re_a1_1, Re_a1_2, Im_a1_0, Im_a1_1, Im_a1_2 over Complex Field with 53 bits of precision

        """
        gens = []

        if not triangles:
            triangles = list(self._surface.label_iterator())

        for t in sorted(set(triangles)):
            gens += [f"Re_a{t}_{n}" for n in range(self._prec)]
            gens += [f"Im_a{t}_{n}" for n in range(self._prec)]

        # TODO: Should we use a better/configured base ring here?

        from sage.all import PolynomialRing, CC
        return PolynomialRing(CC, gens)

    @cached_method
    def gen(self, triangle, k, ring=None, conjugate=False):
        real = self.real(triangle, k, ring=ring)
        imag = self.imag(triangle, k, ring=ring)

        I = imag.parent().base_ring().gen()

        if conjugate:
            I = -I

        assert I*I == -1

        return real + I*imag


    @cached_method
    def real(self, triangle, k, ring=None):
        r"""
        Return the real part of the kth generator of the :meth:`symbolic_ring`
        for ``triangle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.real(0, 0)
            Re_a0_0
            sage: C.real(0, 1)
            Re_a0_1
            sage: C.real(1, 2)
            Re_a1_2

        """
        if k >= self._prec:
            raise ValueError("symbolic ring has no k-th generator")

        if ring is None:
            return self.symbolic_ring(triangle).gen(k)

        return ring.gen(ring.variable_names().index(f"Re_a{triangle}_{k}"))

    @cached_method
    def imag(self, triangle, k, ring=None):
        r"""
        Return the imaginary part of the kth generator of the :meth:`symbolic_ring`
        for ``triangle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.imag(0, 0)
            Im_a0_0
            sage: C.imag(0, 1)
            Im_a0_1
            sage: C.imag(1, 2)
            Im_a1_2

        """
        if k >= self._prec:
            raise ValueError("symbolic ring has no k-th generator")

        if ring is None:
            return self.symbolic_ring(triangle).gen(self._prec + k)

        return ring.gen(ring.variable_names().index(f"Im_a{triangle}_{k}"))

    @cached_method
    def _describe_generator(self, gen):
        r"""
        Return which kind of symbolic generator ``gen`` is.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C._describe_generator(C.imag(1, 2))
            ('imag', 1, 2)
            sage: C._describe_generator(C.real(2, 1))
            ('real', 2, 1)

        """
        gen = str(gen)
        if gen.startswith("Re_a"):
            kind = "real"
            gen = gen[4:]
        elif gen.startswith("Im_a"):
            kind = "imag"
            gen = gen[4:]

        triangle, k = gen.split('_')

        return kind, int(triangle), int(k)

    def project(self, x, part):
        r"""
        Return the ``"real"`` or ``"imag"```inary ``part`` of ``x``.
        """
        if part not in ["real", "imag"]:
            raise ValueError("part must be one of real or imag")

        # Return the real part of a complex number.
        if hasattr(x, part):
            x = getattr(x, part)
            if callable(x):
                x = x()
            return x

        return x.map_coefficients(lambda c: self.project(c, part))

    @staticmethod
    def _subs(polynomial, substitutions):
        r"""
        A faster version of multivariate polynomial's ``subs``.

        Unfortunately, ``subs`` is extremely slow for polynomials with lots of
        variables. Part of this are trivialities, namely, ``subs`` stringifies
        all of the generators of the polynomial ring. But also due to the
        evaluation algorithm that is not very fast when most variables are
        unchanged.
        """
        R = polynomial.parent()
        gens = R.gens()

        result = R.zero()

        substituted_generator_indexes = [i for i, gen in enumerate(gens) if gen in substitutions]

        for coefficient, monomial, exponents in zip(polynomial.coefficients(), polynomial.monomials(), polynomial.exponents()):
            for index in substituted_generator_indexes:
                if exponents[index]:
                    break
            else:
                # monomial is unaffected by this substitution
                result += coefficient * monomial
                continue

            for i, (gen, exponent) in enumerate(zip(gens, exponents)):
                if not exponent:
                    continue
                if gen in substitutions:
                    coefficient *= substitutions[gen] ** exponent

            result += coefficient

        return result

    def real_part(self, x):
        r"""
        Return the real part of ``x``.

        EXAMPLES::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: C = PowerSeriesConstraints(T, prec=3)

            sage: C.real_part(1 + I)  # tol 1e-9
            1
            sage: C.real_part(1.)  # tol 1e-9
            1

        ::

            sage: C.real_part(C.gen(0, 0))
            Re_a0_0
            sage: C.real_part(C.real(0, 0))
            Re_a0_0
            sage: C.real_part(C.imag(0, 0))
            Im_a0_0
            sage: C.real_part(2*C.gen(0, 0))  # tol 1e-9
            2*Re_a0_0
            sage: C.real_part(2*I*C.gen(0, 0))  # tol 1e-9
            (-2)*Im_a0_0

        """
        return self.project(x, "real")

    def imaginary_part(self, x):
        r"""
        Return the imaginary part of ``x``.

        EXAMPLES::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()
            sage: C = PowerSeriesConstraints(T, prec=3)

            sage: C.imaginary_part(1 + I)  # tol 1e-9
            1
            sage: C.imaginary_part(1.)  # tol 1e-9
            0

        ::

            sage: C.imaginary_part(C.gen(0, 0))
            Im_a0_0
            sage: C.imaginary_part(C.real(0, 0))
            0
            sage: C.imaginary_part(C.imag(0, 0))
            0
            sage: C.imaginary_part(2*C.gen(0, 0))  # tol 1e-9
            2*Im_a0_0
            sage: C.imaginary_part(2*I*C.gen(0, 0))  # tol 1e-9
            2*Re_a0_0

        """
        return self.project(x, "imag")

    def add_constraint(self, expression, value=ZZ(0), lagrange=[]):
        expression -= value

        if expression == 0:
            return

        if expression.total_degree() == 0:
            raise ValueError("cannot solve for constraints c == 0")

        if expression.total_degree() > 1:
            raise NotImplementedError("can only encode linear constraints")

        value = -expression.constant_coefficient()

        # We encode a constraint Σ c_i a_i = v as its real and imaginary part.
        # (Our solver can handle complex systems but we also want to add
        # constraints that only concern the real part of the a_i.)

        for part in [self.real_part, self.imaginary_part]:
            e = part(expression)

            real = {}
            imag = {}

            for gen in e.variables():
                kind, triangle, k = self._describe_generator(gen)

                assert kind in ["real", "imag"]

                bucket = real if kind == "real" else imag

                coefficients = bucket.setdefault(triangle, [])

                if len(coefficients) <= k:
                    coefficients.extend([0]*(k + 1 - len(coefficients)))

                coefficients[k] = e[gen]

            self._add_constraint(
                real=real,
                imag=imag,
                lagrange=[part(l) for l in lagrange],
                value=part(value))

    def _add_constraint(self, real, imag, value, lagrange=[]):
        # Simplify constraint by dropping zero coefficients.
        for triangle in real:
            while real[triangle] and not real[triangle][-1]:
                real[triangle].pop()

        for triangle in imag:
            while imag[triangle] and not imag[triangle][-1]:
                imag[triangle].pop()

        # Simplify constraints by dropping trivial constraints
        real = {triangle: real[triangle] for triangle in real if real[triangle]}
        imag = {triangle: imag[triangle] for triangle in imag if imag[triangle]}

        # Simplify lagrange constraints
        while lagrange and not lagrange[-1]:
            lagrange.pop()

        # Ignore trivial constraints
        if not real and not imag and not value and not lagrange:
            return

        constraint = PowerSeriesConstraints.Constraint(real=real, imag=imag, value=value, lagrange=lagrange)

        if constraint not in self._constraints:
            self._constraints.append(constraint)

    @cached_method
    def _formal_power_series(self, triangle, base_ring=None):
        if base_ring is None:
            base_ring = self.symbolic_ring(triangle)

        from sage.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, 'z')

        return R([self.gen(triangle, n, ring=base_ring) for n in range(self._prec)])

    def develop(self, triangle, Δ=0, base_ring=None):
        r"""
        Return the power series obtained by developing at z + Δ.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.develop(0)
            Re_a0_0 + 1.00000000000000*I*Im_a0_0 + (Re_a0_1 + 1.00000000000000*I*Im_a0_1)*z + (Re_a0_2 + 1.00000000000000*I*Im_a0_2)*z^2
            sage: C.develop(1, 1)
            Re_a1_0 + Re_a1_1 + Re_a1_2 + 1.00000000000000*I*Im_a1_0 + 1.00000000000000*I*Im_a1_1 + 1.00000000000000*I*Im_a1_2 + (Re_a1_1 + 2.00000000000000*Re_a1_2 + 1.00000000000000*I*Im_a1_1 + 2.00000000000000*I*Im_a1_2)*z + (Re_a1_2 + 1.00000000000000*I*Im_a1_2)*z^2

        """
        # TODO: Check that Δ is within the radius of convergence.
        f = self._formal_power_series(triangle, base_ring=base_ring)
        return f(f.parent().gen() + Δ)

    def integrate(self, cycle):
        r"""
        Return the linear combination of the power series coefficients that
        describe the integral of a differential along the homology class
        ``cycle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=5)

            sage: C.integrate(H())
            0

            sage: a, b = H.gens()
            sage: C.integrate(a)  # tol 2e-3
            (0.500 + 0.500*I)*Re_a0_0 + (-0.250)*Re_a0_1 + (0.0416 - 0.0416*I)*Re_a0_2 + (0.00625 + 0.00625*I)*Re_a0_4 + (-0.500 + 0.500*I)*Im_a0_0 + (-0.250*I)*Im_a0_1 + (0.0416 + 0.0416*I)*Im_a0_2 + (-0.00625 + 0.00625*I)*Im_a0_4 + (0.500 + 0.500*I)*Re_a1_0 + 0.250*Re_a1_1 + (0.0416 - 0.0416*I)*Re_a1_2 + (0.00625 + 0.00625*I)*Re_a1_4 + (-0.500 + 0.500*I)*Im_a1_0 + 0.250*I*Im_a1_1 + (0.0416 + 0.0416*I)*Im_a1_2 + (-0.00625 + 0.00625*I)*Im_a1_4
            sage: C.integrate(b)  # tol 2e-3
            (-0.500)*Re_a0_0 + 0.125*Re_a0_1 + (-0.0416)*Re_a0_2 + 0.0156*Re_a0_3 + (-0.00625)*Re_a0_4 + (-0.500*I)*Im_a0_0 + 0.125*I*Im_a0_1 + (-0.0416*I)*Im_a0_2 + 0.0156*I*Im_a0_3 + (-0.00625*I)*Im_a0_4 + (-0.500)*Re_a1_0 + (-0.125)*Re_a1_1 + (-0.0416)*Re_a1_2 + (-0.0156)*Re_a1_3 + (-0.00625)*Re_a1_4 + (-0.500*I)*Im_a1_0 + (-0.125*I)*Im_a1_1 + (-0.0416*I)*Im_a1_2 + (-0.0156*I)*Im_a1_3 + (-0.00625*I)*Im_a1_4

        """
        surface = cycle.surface()

        R = self.symbolic_ring(*[triangle for path in cycle.voronoi_path().monomial_coefficients().keys() for triangle, _ in path])

        expression = R.zero()

        for path, multiplicity in cycle.voronoi_path().monomial_coefficients().items():

            for S, T in zip((path[-1],) + path, path):
                # Integrate from the midpoint of the edge of S to the midpoint of the edge of T
                S = surface.opposite_edge(*S)
                assert S[0] == T[0], f"consecutive elements of a path must be attached to the same face in {path} but {S} and {T} do not have that property"

                # Namely we integrate the power series defined around the Voronoi vertex of S by symbolically integrating each monomial term.

                # The midpoints of the edges
                P = self._geometry.midpoint(*S)
                Q = self._geometry.midpoint(*T)

                P_power = P
                Q_power = Q

                for k in range(self._prec):
                    expression += multiplicity * self.gen(S[0], k, R) / (k + 1) * (Q_power - P_power)

                    P_power *= P
                    Q_power *= Q

        return expression

    def evaluate(self, triangle, Δ, derivative=0):
        r"""
        Return the value of the power series evaluated at Δ in terms of
        symbolic variables.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.evaluate(0, 0)
            Re_a0_0 + 1.00000000000000*I*Im_a0_0
            sage: C.evaluate(1, 0)
            Re_a1_0 + 1.00000000000000*I*Im_a1_0
            sage: C.evaluate(1, 2)
            Re_a1_0 + 2.00000000000000*Re_a1_1 + 4.00000000000000*Re_a1_2 + 1.00000000000000*I*Im_a1_0 + 2.00000000000000*I*Im_a1_1 + 4.00000000000000*I*Im_a1_2

        """
        # TODO: Check that Δ is within the radius of convergence.

        if derivative >= self._prec:
            raise ValueError

        from sage.all import factorial
        return self.develop(triangle=triangle, Δ=Δ)[derivative] * factorial(derivative)

    def require_equality(self):
        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            parent = self.symbolic_ring(triangle0, triangle1)

            Δ0 = self._geometry.midpoint(triangle0, edge0)
            Δ1 = self._geometry.midpoint(triangle1, edge1)

            if abs(Δ0) < 1e-6 and abs(Δ1) < 1e-6:
                # Force power series to be identical if the Delaunay triangulation is ambiguous at this edge.
                for k in range(self._prec):
                    self.add_constraint(self.gen(triangle0, k, parent) - self.gen(triangle1, k, parent))

    def require_midpoint_derivatives(self, derivatives):
        r"""
        The radius of convergence of the power series is the distance from the
        vertex of the Voronoi cell to the closest singularity of the
        triangulation (since we use a Delaunay triangulation, all vertices are
        at the same distance in fact.) So the radii of convergence of two
        neigbhouring cells overlap and the power series must coincide there.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

        This example is a bit artificial. Both power series are developed
        around the same point since a common edge is ambiguous in the Delaunay
        triangulation. Therefore, we require all coefficients to be identical.
        However, we also require the power series to be compatible with itself
        away from the midpoint which cannot be seen in this example because the
        precision is too low::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_midpoint_derivatives(1)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1], 1: [-1]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1], 1: [-1]}, lagrange=[], value=0)]

        If we add more coefficients, we get three pairs of contraints for the
        three edges surrounding a face; for the edge on which the centers of
        the Voronoi cells fall, we get a more restrictive constraint forcing
        all the coefficients to be equal (these are the first four constraints
        recorded below.)::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(1)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [0.0, -0.50], 1: [0.0, -0.50]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 0.50], 1: [0.0, 0.50]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [1.0, -0.50], 1: [-1.0, -0.50]}, imag={}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.50], 1: [-1.0, -0.50]}, lagrange=[], value=-0.0)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(2)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [0, -0.50], 1: [0, -0.50]}, lagrange=[], value=-0.0),
            PowerSeriesConstraints.Constraint(real={0: [0, 0.50], 1: [0, 0.50]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=-0.0),
            PowerSeriesConstraints.Constraint(real={0: [0, 1.0], 1: [0, -1.0]}, imag={}, lagrange=[], value=-0.0),
            PowerSeriesConstraints.Constraint(real={}, imag={0: [0, 1.0], 1: [0, -1.0]}, lagrange=[], value=-0.0),
            PowerSeriesConstraints.Constraint(real={0: [1.0, -0.50], 1: [-1.0, -0.50]}, imag={}, lagrange=[], value=-0.0),
            PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.50], 1: [-1.0, -0.50]}, lagrange=[], value=-0.0)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            parent = self.symbolic_ring(triangle0, triangle1)

            Δ0 = self._geometry.midpoint(triangle0, edge0)
            Δ1 = self._geometry.midpoint(triangle1, edge1)

            # TODO: Are these good constants?
            if abs(Δ0) < 1e-6 and abs(Δ1) < 1e-6:
                continue

            # Require that the 0th, ..., derivatives-1th derivatives are the same at the midpoint of the edge.
            for derivative in range(derivatives):
                self.add_constraint(
                    parent(self.evaluate(triangle0, Δ0, derivative)) - parent(self.evaluate(triangle1, Δ1, derivative)))

    def _L2_consistency(self):
        r"""
        For each pair of adjacent triangles meeting at and edge `e`, let `v` be
        the midpoint of `e`. We develop the power series coming from both
        triangles around that midpoint and check them for agreement. Namely, we
        integrate the square of their difference on the circle of maximal
        radius around `v` as a line integral. (If the power series agree, that
        integral should be zero.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

        This example is a bit artificial. Both power series are developed
        around the same point since a common edge is ambiguous in the Delaunay
        triangulation. Therefore, we require all coefficients to be identical.
        However, we also require the power series to be compatible with itself
        away from the midpoint::

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f, prec=6)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: consistency = PowerSeriesConstraints(T, η.precision())._L2_consistency()
            sage: η._evaluate(consistency)  # tol 1e-9
            0

        """
        R = self.symbolic_ring()

        cost = R.zero()

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            # The midpoint of the edge where the triangles meet with respect to
            # the center of the triangle.
            Δ0 = self._geometry.midpoint(triangle0, edge0)
            Δ1 = self._geometry.midpoint(triangle1, edge1)

            # TODO: Should we skip such steps here?
            # if abs(Δ0) < 1e-6 and abs(Δ1) < 1e-6:

            # Develop both power series around that midpoint, i.e., Taylor expand them.
            T0 = self.develop(triangle0, Δ0, base_ring=R)
            T1 = self.develop(triangle1, Δ1, base_ring=R)

            # Write b_n for the difference of the n-th coefficient of both power series.
            # We want to minimize the sum of |b_n|^2 r^2n where r is half the
            # length of the edge we are on.
            # TODO: What is the correct exponent here actually?
            b = (T0 - T1).list()
            edge = self._surface.polygon(triangle0).edges()[edge0]
            r2 = (edge[0]**2 + edge[1]**2) / 4

            for n, b_n in enumerate(b):
                cost += (self.real_part(b_n)**2 + self.imaginary_part(b_n)**2) * r2**(n + 1)

        return cost

    @cached_method
    def _elementary_line_integrals(self, triangle, n, m):
        r"""
        Return the integrals f(z)dx and f(z)dy where f(z) = z^n\overline{z}^m
        along the boundary of the ``triangle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 3)

            sage: C._elementary_line_integrals(0, 0, 0)
            (0.0, 0.0)
            sage: C._elementary_line_integrals(0, 1, 0)  # tol 1e-9
            (0 - 0.5*I, 0.5 + 0.0*I)
            sage: C._elementary_line_integrals(0, 0, 1)  # tol 1e-9
            (0.0 + 0.5*I, 0.5 - 0.0*I)
            sage: C._elementary_line_integrals(0, 1, 1)  # tol 1e-9
            (-0.1666666667, -0.1666666667)

        """
        from sage.all import I

        ix = 0
        iy = 0

        triangle = self._surface.polygon(triangle)
        center = triangle.circumscribing_circle().center()

        for v, e in zip(triangle.vertices(), triangle.edges()):
            Δx, Δy = e
            x0, y0 = -center + v

            def f(x, y):
                from sage.all import I
                return complex((x + I*y)**n * (x - I*y)**m)

            def fx(t):
                if abs(Δx) < 1e-6:
                    return complex(0)
                return f(x0 + t, y0 + t * Δy/Δx)

            def fy(t):
                if abs(Δy) < 1e-6:
                    return complex(0)
                return f(x0 + t * Δx/Δy, y0 + t)

            def integrate(value, t0, t1):
                from scipy.integrate import quad
                return quad(value, t0, t1)[0]

            ix += integrate(lambda t: fx(t).real, 0, Δx) + I * integrate(lambda t: fx(t).imag, 0, Δx)
            iy += integrate(lambda t: fy(t).real, 0, Δy) + I * integrate(lambda t: fy(t).imag, 0, Δy)

        return ix, iy

    @cached_method
    def _elementary_area_integral(self, triangle, n, m):
        r"""
        Return the integral of z^n\overline{z}^m on the ``triangle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 3)
            sage: C._elementary_area_integral(0, 0, 0)  # tol 1e-9
            0.5 + 0.0*I

            sage: C._elementary_area_integral(0, 1, 0)  # tol 1e-6
            -0.083333333 + 0.083333333*I

        """
        # Write f(n, m) for z^n\overline{z}^m.
        # Then 1/(2m + 1) [d/dx f(n, m+1) - d/dy -i f(n, m+1)] = f(n, m).

        # So we can use Green's theorem to compute this integral by integrating
        # on the boundary of the triangle:
        # -i/(2m + 1) f(n, m + 1) dx + 1/(2m + 1) f(n, m + 1) dy

        ix, iy = self._elementary_line_integrals(triangle, n, m+1)

        from sage.all import I
        return -I/(2*(m + 1)) * ix + 1/(2*(m + 1)) * iy

    def _area(self):
        r"""
        Return a formula for the area ∫ η \wedge \overline{η}.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f, prec=6)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: area = PowerSeriesConstraints(T, η.precision())._area()

            sage: η._evaluate(area)  # tol 1e-9
            1.0 + 0.0*I

        """
        R = self.symbolic_ring()

        area = R.zero()

        # We eveluate the integral of |f|^2 on each triangle where f is the
        # power series on the Voronoy cell containing that triangle.
        for triangle in self._surface.label_iterator():
            # We expand the integrand Σa_nz^n · \overline{Σa_mz^m} naively as
            # the sum of a_n \overline{a_m} z^n \overline{z^m}.
            for n in range(self._prec):
                for m in range(self._prec):
                    coefficient = self.gen(triangle, n, R) * self.gen(triangle, m, R, conjugate=True)
                    # Now we have to integrate z^n \overline{z^m} on the triangle.
                    area += coefficient * self._elementary_area_integral(triangle, n, m)

        return area

    def _area_upper_bound(self):
        r"""
        Return an upper bound for the area 1 /(2iπ) ∫ η \wedge \overline{η}.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: area = PowerSeriesConstraints(T, η.precision())._area_upper_bound()
            sage: area  # tol 2e-2
            0.250*Re_a0_0^2 + 0.125*Re_a0_1^2 + 0.0625*Re_a0_2^2 + 0.0312*Re_a0_3^2 + 0.0156*Re_a0_4^2 + 0.00781*Re_a0_5^2 + 0.00390*Re_a0_6^2 + 0.00195*Re_a0_7^2 + 0.000976*Re_a0_8^2 + 0.000488*Re_a0_9^2 + 0.250*Im_a0_0^2 + 0.125*Im_a0_1^2 + 0.0625*Im_a0_2^2 + 0.0312*Im_a0_3^2 + 0.0156*Im_a0_4^2 + 0.00781*Im_a0_5^2 + 0.00390*Im_a0_6^2 + 0.00195*Im_a0_7^2 + 0.000976*Im_a0_8^2 + 0.000488*Im_a0_9^2
            + 0.250*Re_a1_0^2 + 0.125*Re_a1_1^2 + 0.0625*Re_a1_2^2 + 0.0312*Re_a1_3^2 + 0.0156*Re_a1_4^2 + 0.00781*Re_a1_5^2 + 0.00390*Re_a1_6^2 + 0.00195*Re_a1_7^2 + 0.000976*Re_a1_8^2 + 0.000488*Re_a1_9^2 + 0.250*Im_a1_0^2 + 0.125*Im_a1_1^2 + 0.0625*Im_a1_2^2 + 0.0312*Im_a1_3^2 + 0.0156*Im_a1_4^2 + 0.00781*Im_a1_5^2 + 0.00390*Im_a1_6^2 + 0.00195*Im_a1_7^2 + 0.000976*Im_a1_8^2 + 0.000488*Im_a1_9^2

        The correct area would be 1/2π here. However, we are overcounting
        because we sum the single Voronoi cell twice. And also, we approximate
        the square with a circle for another factor π/2::

            sage: η._evaluate(area)  # tol 1e-2
            .5

        """
        # TODO: Verify that this is actually correct. I think the formula below
        # is offf. It's certainly not in sync with the documentation.

        # To make our lives easier, we do not optimize this value but instead
        # half the sum of the |a_k|^2·radius^(k+2) = (Re(a_k)^2 +
        # Im(a_k)^2)·radius^(k+2) which is a very rough upper bound for the
        # area.
        area = self.symbolic_ring().zero()

        for triangle in self._surface.label_iterator():
            R = float(self._surface.polygon(triangle).circumscribing_circle().radius_squared().sqrt())

            for k in range(self._prec):
                area += (self.real(triangle, k)**2 + self.imag(triangle, k)**2) * R**(2*k + 2)

        return area/2

    def optimize(self, f):
        r"""
        Add constraints that optimize the symbolic expression ``f``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

        ::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_midpoint_derivatives(1)
            sage: R = C.symbolic_ring()
            sage: f = 3*C.real(0, 0, R)^2 + 5*C.imag(0, 0, R)^2 + 7*C.real(1, 0, R)^2 + 11*C.imag(1, 0, R)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            ...
            PowerSeriesConstraints.Constraint(real={0: [6.00000000000000]}, imag={}, lagrange=[-1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={}, imag={0: [10.0000000000000]}, lagrange=[0, -1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={1: [14.0000000000000]}, imag={}, lagrange=[1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={}, imag={1: [22.0000000000000]}, lagrange=[0, 1.00000000000000], value=-0.000000000000000)]

        """
        self._cost += self.symbolic_ring()(f)

    def _optimize_cost(self):
        # We use Lagrange multipliers to rewrite this expression.
        # If we let
        #   L(Re(a), Im(a), λ) = f(Re(a), Im(a)) - Σ λ_i g_i(Re(a), Im(a))
        # and denote by g_i=0 all the affine linear conditions collected so
        # far, then we get two constraints for each a_k, one real, one
        # imaginary, namely that the partial derivative wrt Re(a_k) and Im(a_k)
        # vanishes. Note that we have one Lagrange multiplier for each affine
        # linear constraint collected so far.
        lagranges = len(self._constraints)

        g = self._constraints

        # We form the partial derivative with respect to the variables Re(a_k)
        # and Im(a_k).
        for triangle in range(self._surface.num_polygons()):
            for k in range(self._prec):
                for part in ["real", "imag"]:
                    gen = getattr(self, part)(triangle, k)

                    if self._cost.degree(gen) <= 0:
                        continue

                    gen = self._cost.parent()(gen)

                    from more_itertools import nth

                    lagrange = [nth(getattr(g[i], part).get(triangle, []), k, 0) for i in range(lagranges)]

                    self.add_constraint(self._cost.derivative(gen), lagrange=lagrange)

        # We form the partial derivatives with respect to the λ_i. This yields
        # the condition -g_i=0 which is already recorded in the linear system.

        # Prevent us from calling this method again.
        self._cost = None

    def require_cohomology(self, cocycle):
        r""""
        Create a constraint by integrating numerically following the paths that
        form a basis of homology.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()

        ::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_cohomology(H({a: 1}))
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [0.5, -0.25], 1: [0.5, 0.25]}, imag={0: [-0.5], 1: [-0.5]}, lagrange=[], value=1),
             PowerSeriesConstraints.Constraint(real={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, imag={}, lagrange=[], value=0)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [0.5, -0.25], 1: [0.5, 0.25]}, imag={0: [-0.5], 1: [-0.5]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, imag={}, lagrange=[], value=1)]

        """
        for cycle in cocycle.parent().homology().gens():
            self.add_constraint(self.real_part(self.integrate(cycle)), self.real_part(cocycle(cycle)))

    def matrix(self):
        lagranges = max(len(constraint.lagrange) for constraint in self._constraints)
        triangles = list(self._surface.label_iterator())

        prec = int(self._prec)

        import numpy

        A = numpy.zeros((len(self._constraints), 2*len(triangles)*prec + lagranges))
        b = numpy.zeros((len(self._constraints),))

        for row, constraint in enumerate(self._constraints):
            b[row] = constraint.value
            for block, triangle in enumerate(triangles):
                for column, value in enumerate(constraint.real.get(triangle, [])):
                    A[row, block * (2*prec) + column] = value
                for column, value in enumerate(constraint.imag.get(triangle, [])):
                    A[row, block * (2*prec) + prec + column] = value

            for column, value in enumerate(constraint.lagrange):
                A[row, 2*len(triangles)*prec + column] = value

        return A, b

    def solve(self):
        r"""
        Return a solution for the system of constraints with minimal error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C._add_constraint(real={0: [-1], 1: [1]}, imag={}, value=0)
            sage: C._add_constraint(imag={0: [-1], 1: [1]}, real={}, value=0)
            sage: C._add_constraint(real={0: [1]}, imag={}, value=1)
            sage: C.solve()
            {0: 1.00000000000000 + O(z0), 1: 1.00000000000000 + O(z1)}

        """
        self._optimize_cost()

        A, b = self.matrix()

        import scipy.linalg
        solution, residues, _, _ = scipy.linalg.lstsq(A, b, check_finite=False, overwrite_a=True, overwrite_b=True)

        lagranges = max(len(constraint.lagrange) for constraint in self._constraints)

        if lagranges:
            solution = solution[:-lagranges]

        solution = [solution[2*k*self._prec:2*(k+1)*self._prec] for k in range(self._surface.num_polygons())]

        from sage.all import CC
        return {
            triangle: HarmonicDifferentials(self._surface).power_series_ring(triangle)([CC(solution[k][p], solution[k][p + self._prec]) for p in range(self._prec)]).add_bigoh(self._prec)
            for (k, triangle) in enumerate(self._surface.label_iterator())
        }
