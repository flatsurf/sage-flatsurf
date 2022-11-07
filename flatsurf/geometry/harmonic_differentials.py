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
from sage.all import ZZ, CC
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement
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
            sage: η = Ω(f); η  # TODO: WRONG OUTPUT
            (-1.02907963995315e-15 + 0.992449699487969*I + (-1.05471187339390e-15 + 3.46944695195361e-16*I)*z0 + (1.22124532708767e-15 - 0.102080540518758*I)*z0^2 + (2.55351295663786e-15 + 4.08006961549745e-15*I)*z0^3 + (-1.55431223447522e-15 + 0.0286707690963231*I)*z0^4 + (-8.32667268468867e-15 - 5.82867087928207e-16*I)*z0^5 + (-8.16013923099490e-15 + 0.544429549433375*I)*z0^6 + (-5.32907051820075e-15 - 9.08995101411847e-15*I)*z0^7 + (7.04991620636974e-15 - 0.229366152770579*I)*z0^8 + (3.65263375101677e-14 + 4.88498130835069e-15*I)*z0^9 + O(z0^10), -4.02455846426619e-16 + 0.992449699487968*I + (2.49800180540660e-16*I)*z1 + (1.60982338570648e-15 - 0.102080540518759*I)*z1^2 + (1.91513471747840e-15 + 3.96904731303493e-15*I)*z1^3 + (-1.27675647831893e-15 + 0.0286707690963231*I)*z1^4 + (-7.88258347483861e-15 - 6.38378239159465e-16*I)*z1^5 + (-8.33361157859258e-15 + 0.544429549433375*I)*z1^6 + (-4.90579799006241e-15 - 8.95117313604032e-15*I)*z1^7 + (6.90593415786367e-15 - 0.229366152770578*I)*z1^8 + (3.65540930857833e-14 + 4.54324078358326e-15*I)*z1^9 + O(z1^10))
            sage: # (0 - 1.*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
            ....: #  0 - 1.*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

            sage: η + η  # TODO: WRONG OUTPUT
            (-2.05815927990630e-15 + 1.98489939897594*I + (-2.10942374678780e-15 + 6.93889390390723e-16*I)*z0 + (2.44249065417534e-15 - 0.204161081037516*I)*z0^2 + (5.10702591327572e-15 + 8.16013923099490e-15*I)*z0^3 + (-3.10862446895044e-15 + 0.0573415381926462*I)*z0^4 + (-1.66533453693773e-14 - 1.16573417585641e-15*I)*z0^5 + (-1.63202784619898e-14 + 1.08885909886675*I)*z0^6 + (-1.06581410364015e-14 - 1.81799020282369e-14*I)*z0^7 + (1.40998324127395e-14 - 0.458732305541158*I)*z0^8 + (7.30526750203353e-14 + 9.76996261670138e-15*I)*z0^9 + O(z0^10), -8.04911692853238e-16 + 1.98489939897594*I + (4.99600361081320e-16*I)*z1 + (3.21964677141295e-15 - 0.204161081037517*I)*z1^2 + (3.83026943495679e-15 + 7.93809462606987e-15*I)*z1^3 + (-2.55351295663786e-15 + 0.0573415381926463*I)*z1^4 + (-1.57651669496772e-14 - 1.27675647831893e-15*I)*z1^5 + (-1.66672231571852e-14 + 1.08885909886675*I)*z1^6 + (-9.81159598012482e-15 - 1.79023462720806e-14*I)*z1^7 + (1.38118683157273e-14 - 0.458732305541157*I)*z1^8 + (7.31081861715666e-14 + 9.08648156716652e-15*I)*z1^9 + O(z1^10))
            sage: # (0 - 2*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
            ....: #  0 - 2*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

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
            sage: η._evaluate(R(C.gen(0, 0)) + R(C.gen(1, 0)))  # TODO: WRONG OUTPUT, too much imprecision
            -1.43153548637977e-15 + 1.98489939897594*I
            sage: # 0 - 2*I

        """
        coefficients = {}

        for gen in expression.variables():
            triangle, k, kind = gen.gen()
            coefficient = self._series[triangle][k]

            if kind == "real":
                coefficients[gen] = coefficient.real()
            else:
                assert kind == "imag"
                coefficients[gen] = coefficient.imag()

        value = expression(coefficients)
        if isinstance(value, SymbolicCoefficientExpression):
            assert value.total_degree() <= 0
            value = value.constant_coefficient()

        return value

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
            def check(actual, expected, message, abs_error_bound=1e-9, rel_error_bound=1e-6):
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
    # TODO: Run test suite
    # TODO: Make sure that we never have zero coefficients as these would break degree computations.

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


class SymbolicCoefficientExpression(CommutativeRingElement):
    def __init__(self, parent, coefficients, constant):
        super().__init__(parent)

        self._coefficients = coefficients
        self._constant = constant

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not (self == other)

        if op == op_EQ:
            return self._constant == other._constant and self._coefficients == other._coefficients

        raise NotImplementedError

    def _repr_(self):
        terms = self.items()

        if self.is_constant():
            return repr(self._constant)

        def variable_name(triangle, k, kind):
            if kind == "real":
                prefix = "Re"
            elif kind == "imag":
                prefix = "Im"
            else:
                assert False

            return f"{prefix}_a{triangle}_{k}"

        variable_names = [variable_name(*gen) for gen in sorted(set(gen for (monomial, coefficient) in terms for gen in monomial))]

        from sage.all import PolynomialRing
        R = PolynomialRing(self.base_ring(), tuple(variable_names))

        def polynomial_monomial(monomial):
            from sage.all import prod
            return prod([R(variable_name(*gen))**exponent for (gen, exponent) in monomial.items()])

        f = sum(coefficient * polynomial_monomial(monomial) for (monomial, coefficient) in terms)

        return repr(f)

    def degree(self, gen):
        if not isinstance(gen, SymbolicCoefficientExpression):
            raise NotImplementedError

        if not gen.is_monomial():
            raise ValueError

        if self.is_zero():
            return -1

        key = next(iter(gen._coefficients))

        degree = 0

        while isinstance(self, SymbolicCoefficientExpression):
            self = self._coefficients.get(key, None)

            if self is None:
                break

            degree += 1

        return degree

    def is_monomial(self):
        return len(self._coefficients) == 1 and not self._constant and next(iter(self._coefficients.values())) == 1

    def is_constant(self):
        return not self._coefficients

    def _neg_(self):
        parent = self.parent()
        return parent.element_class(parent, {key: -coefficient for (key, coefficient) in self._coefficients.items()}, -self._constant)

    def _add_(self, other):
        parent = self.parent()

        if len(self._coefficients) < len(other._coefficients):
            self, other = other, self

        coefficients = self._coefficients | other._coefficients

        for key in other._coefficients:
            c = self._coefficients.get(key)
            if c is not None:
                coefficients[key] += c

        return parent.element_class(parent, coefficients, self._constant + other._constant)

    def _sub_(self, other):
        return self._add_(-other)

    def _mul_(self, other):
        parent = self.parent()

        if other.is_zero() or self.is_zero():
            return parent.zero()

        if other.is_one():
            return self

        if self.is_one():
            return other

        if other.is_constant():
            constant = other._constant
            return parent({key: constant * value for (key, value) in self._coefficients.items()}, constant * self._constant)

        if self.is_constant():
            return other * self

        value = parent.zero()

        for (monomial, coefficient) in self.items():
            for (monomial_, coefficient_) in other.items():
                if not monomial and not monomial_:
                    value += coefficient * coefficient_
                    continue

                from copy import copy
                monomial__ = copy(monomial)

                for (gen, exponent) in monomial_.items():
                    monomial__.setdefault(gen, 0)
                    monomial__[gen] += exponent

                coefficient__ = coefficient * coefficient_

                if not monomial__:
                    value += coefficient__
                    continue

                def unfold(monomial, coefficient):
                    if not monomial:
                        return coefficient

                    unfolded = {}
                    for gen, exponent in monomial.items():
                        unfolded[gen] = unfold({
                            g: e if g != gen else e - 1
                            for (g, e) in monomial.items()
                            if g != gen or e != 1
                        }, coefficient)

                    return parent(unfolded)

                value += unfold(monomial__, coefficient__)

        return value

    def _rmul_(self, right):
        return self._lmul_(right)

    def _lmul_(self, left):
        return self.parent()({key: left * value for (key, value) in self._coefficients.items()}, self._constant * left)

    def constant_coefficient(self):
        return self._constant

    def variables(self):
        return [self.parent()({variable: 1}) for variable in self._coefficients]

    def __getitem__(self, gen):
        if not gen.is_monomial():
            raise ValueError

        return self._coefficients.get(next(iter(gen._coefficients.keys())), 0)

    def __hash__(self):
        return hash((tuple(sorted(self._coefficients.items())), self._constant))

    def total_degree(self):
        if not self._coefficients:
            if not self._constant:
                return -1
            return 0

        degree = 1

        for key, coefficient in self._coefficients.items():
            if not isinstance(coefficient, SymbolicCoefficientExpression):
                continue
            degree = max(degree, 1 + coefficient.total_degree())

        return degree

    def derivative(self, gen):
        key = gen.gen()

        # Compute derivative with product rule
        value = self._coefficients.get(key, self.parent().zero())
        if value and isinstance(value, SymbolicCoefficientExpression):
            value += gen * value.derivative(gen)

        return value

    def gen(self):
        if not self.is_monomial():
            raise ValueError

        gen, coefficient = next(iter(self._coefficients.items()))

        if coefficient != 1:
            raise ValueError

        return gen

    def map_coefficients(self, f):
        def g(coefficient):
            if isinstance(coefficient, SymbolicCoefficientExpression):
                return coefficient.map_coefficients(f)
            return f(coefficient)

        return self.parent()({key: g(value) for (key, value) in self._coefficients.items()}, g(self._constant))

    def items(self):
        items = []

        def collect(element, prefix=()):
            if not isinstance(element, SymbolicCoefficientExpression):
                if element:
                    items.append((prefix, element))
                return

            if element._constant:
                items.append((prefix, element._constant))

            for key, value in element._coefficients.items():
                if prefix and key < prefix[-1]:
                    # Don't add the same monomial twice.
                    continue

                collect(value, prefix + (key,))

        collect(self)

        def monomial(gens):
            monomial = {}
            for gen in gens:
                monomial.setdefault(gen, 0)
                monomial[gen] += 1

            return monomial

        return [(monomial(gens), coefficient) for (gens, coefficient) in items]

    def __call__(self, values):
        parent = self.parent()

        values = {gen.gen(): value for (gen, value) in values.items()}

        from sage.all import prod
        return sum([
            coefficient * prod([
                (parent({gen: 1}) if gen not in values else values[gen])**e for (gen, e) in monomial.items()
            ]) for (monomial, coefficient) in self.items()])


class SymbolicCoefficientRing(UniqueRepresentation, CommutativeRing):
    def __init__(self, surface, base_ring=CC, category=None):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T)
            sage: R.has_coerce_map_from(CC)
            True

            sage: TestSuite(R).run()

        """
        self._surface = surface

        from sage.categories.all import CommutativeRings
        CommutativeRing.__init__(self, base_ring, category=category or CommutativeRings())
        self.register_coercion(base_ring)

    Element = SymbolicCoefficientExpression

    def _repr_(self):
        return r"Ring of Power Series Coefficients"

    def is_exact(self):
        return self.base_ring().is_exact()

    def base_ring(self):
        return self.base()

    def _element_constructor_(self, x, constant=None):
        if isinstance(x, tuple) and len(x) == 3:
            assert constant is None

            # x describes a monomial
            return self.element_class(self, {x: self.base_ring().one()}, self.base_ring().zero())

        if isinstance(x, dict):
            constant = constant or 0

            return self.element_class(self, x, constant)

        if x in self.base_ring():
            return self.element_class(self, {}, self.base_ring()(x))

        raise NotImplementedError(f"symbolic expression from {x}")

    @cached_method
    def imaginary_unit(self):
        from sage.all import I
        return self(self.base_ring()(I))

    def ngens(self):
        raise NotImplementedError


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
    def symbolic_ring(self):
        r"""
        Return the polynomial ring in the coefficients of the power series of
        the triangles.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.symbolic_ring()
            Ring of Power Series Coefficients

        """
        return SymbolicCoefficientRing(self._surface)

    @cached_method
    def gen(self, triangle, k, conjugate=False):
        real = self.real(triangle, k)
        imag = self.imag(triangle, k)

        I = self.symbolic_ring().imaginary_unit()

        if conjugate:
            I = -I

        return real + I*imag

    @cached_method
    def real(self, triangle, k):
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

        return self.symbolic_ring()((triangle, k, "real"))

    @cached_method
    def imag(self, triangle, k):
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

        return self.symbolic_ring()((triangle, k, "imag"))

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
            -2.0000000000000*Im_a0_0

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

        total_degree = expression.total_degree()

        if total_degree == 0:
            raise ValueError(f"cannot solve for constraint {expression} == 0")

        if total_degree > 1:
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
                triangle, k, kind = gen.gen()

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

        # We could deduplicate constraints here. But it turned out to be more
        # expensive to deduplicate then the price we pay for adding constraints
        # twice. (However, duplicate constraints also increase the weight of
        # that constraint for the lstsq solver which is probably beneficial.)
        self._constraints.append(constraint)

    @cached_method
    def _formal_power_series(self, triangle, base_ring=None):
        if base_ring is None:
            base_ring = self.symbolic_ring()

        from sage.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, 'z')

        return R([self.gen(triangle, n) for n in range(self._prec)])

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
            1.00000000000000*I*Im_a0_0 + Re_a0_0 + (1.00000000000000*I*Im_a0_1 + Re_a0_1)*z + (1.00000000000000*I*Im_a0_2 + Re_a0_2)*z^2
            sage: C.develop(1, 1)
            1.00000000000000*I*Im_a1_0 + Re_a1_0 + 1.00000000000000*I*Im_a1_1 + Re_a1_1 + 1.00000000000000*I*Im_a1_2 + Re_a1_2 + (1.00000000000000*I*Im_a1_1 + Re_a1_1 + 2.00000000000000*I*Im_a1_2 + 2.00000000000000*Re_a1_2)*z + (1.00000000000000*I*Im_a1_2 + Re_a1_2)*z^2

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
            0.000000000000000

            sage: a, b = H.gens()
            sage: C.integrate(a)  # tol 1e-6
            (0.500000000000000 - 0.500000000000000*I)*Im_a0_0 + (0.500000000000000 + 0.500000000000000*I)*Re_a0_0 + 0.250000000000000*I*Im_a0_1 + (-0.250000000000000)*Re_a0_1 + (-0.0416666666666667 - 0.0416666666666667*I)*Im_a0_2 + (0.0416666666666667 - 0.0416666666666667*I)*Re_a0_2 + (0.00625000000000000 - 0.00625000000000000*I)*Im_a0_4 + (0.00625000000000000 + 0.00625000000000000*I)*Re_a0_4 + (0.500000000000000 - 0.500000000000000*I)*Im_a1_0 + (0.500000000000000 + 0.500000000000000*I)*Re_a1_0 + (-0.250000000000000*I)*Im_a1_1 + 0.250000000000000*Re_a1_1 + (-0.0416666666666667 - 0.0416666666666667*I)*Im_a1_2 + (0.0416666666666667 - 0.0416666666666667*I)*Re_a1_2 + (0.00625000000000000 - 0.00625000000000000*I)*Im_a1_4 + (0.00625000000000000 + 0.00625000000000000*I)*Re_a1_4
            sage: C.integrate(b)  # tol 1e-6
            0.500000000000000*I*Im_a0_0 + (-0.500000000000000)*Re_a0_0 + (-0.125000000000000*I)*Im_a0_1 + 0.125000000000000*Re_a0_1 + 0.0416666666666667*I*Im_a0_2 + (-0.0416666666666667)*Re_a0_2 + (-0.0156250000000000*I)*Im_a0_3 + 0.0156250000000000*Re_a0_3 + 0.00625000000000000*I*Im_a0_4 + (-0.00625000000000000)*Re_a0_4 + 0.500000000000000*I*Im_a1_0 + (-0.500000000000000)*Re_a1_0 + 0.125000000000000*I*Im_a1_1 + (-0.125000000000000)*Re_a1_1 + 0.0416666666666667*I*Im_a1_2 + (-0.0416666666666667)*Re_a1_2 + 0.0156250000000000*I*Im_a1_3 + (-0.0156250000000000)*Re_a1_3 + 0.00625000000000000*I*Im_a1_4 + (-0.00625000000000000)*Re_a1_4

        """
        surface = cycle.surface()

        R = self.symbolic_ring()

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
            1.00000000000000*I*Im_a0_0 + Re_a0_0
            sage: C.evaluate(1, 0)
            1.00000000000000*I*Im_a1_0 + Re_a1_0
            sage: C.evaluate(1, 2)
            1.00000000000000*I*Im_a1_0 + Re_a1_0 + 2.00000000000000*I*Im_a1_1 + 2.00000000000000*Re_a1_1 + 4.00000000000000*I*Im_a1_2 + 4.00000000000000*Re_a1_2

        """
        # TODO: Check that Δ is within the radius of convergence.

        if derivative >= self._prec:
            raise ValueError

        parent = self.symbolic_ring()

        value = parent.zero()

        z = 1

        from sage.all import factorial
        factor = factorial(derivative)

        for k in range(derivative, self._prec):
            value += factor * self.gen(triangle, k) * z

            factor *= k + 1
            factor /= k - derivative + 1

            z *= Δ

        return value

    def require_equality(self):
        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            parent = self.symbolic_ring()

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
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000)]

        If we add more coefficients, we get three pairs of contraints for the
        three edges surrounding a face; for the edge on which the centers of
        the Voronoi cells fall, we get a more restrictive constraint forcing
        all the coefficients to be equal (these are the first four constraints
        recorded below.)::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={0: [0.000000000000000, -0.500000000000000], 1: [-0.000000000000000, -0.500000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [0.000000000000000, 0.500000000000000], 1: [-0.000000000000000, 0.500000000000000]}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [1.00000000000000, -0.500000000000000], 1: [-1.00000000000000, -0.500000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000, -0.500000000000000], 1: [-1.00000000000000, -0.500000000000000]}, lagrange=[], value=-0.000000000000000)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(2)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={0: [0.000000000000000, -0.500000000000000], 1: [-0.000000000000000, -0.500000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [0.000000000000000, 0.500000000000000], 1: [-0.000000000000000, 0.500000000000000]}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [0, 1.00000000000000], 1: [0, -1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [0, 1.00000000000000], 1: [0, -1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [1.00000000000000, -0.500000000000000], 1: [-1.00000000000000, -0.500000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000, -0.500000000000000], 1: [-1.00000000000000, -0.500000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [0, 1.00000000000000], 1: [0, -1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [0, 1.00000000000000], 1: [0, -1.00000000000000]}, lagrange=[], value=-0.000000000000000)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            parent = self.symbolic_ring()

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
                    coefficient = self.gen(triangle, n) * self.gen(triangle, m, conjugate=True)
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
            sage: area  # tol 1e-6
            0.250000000000000*Im_a0_0^2 + 0.250000000000000*Re_a0_0^2 + 0.125000000000000*Im_a0_1^2 + 0.125000000000000*Re_a0_1^2 + 0.0625000000000000*Im_a0_2^2 + 0.0625000000000000*Re_a0_2^2 + 0.0312500000000000*Im_a0_3^2 + 0.0312500000000000*Re_a0_3^2 + 0.0156250000000000*Im_a0_4^2 + 0.0156250000000000*Re_a0_4^2 + 0.00781250000000001*Im_a0_5^2 + 0.00781250000000001*Re_a0_5^2 + 0.00390625000000000*Im_a0_6^2 + 0.00390625000000000*Re_a0_6^2 + 0.00195312500000000*Im_a0_7^2 + 0.00195312500000000*Re_a0_7^2 + 0.000976562500000001*Im_a0_8^2 + 0.000976562500000001*Re_a0_8^2 + 0.000488281250000001*Im_a0_9^2 + 0.000488281250000001*Re_a0_9^2 + 0.250000000000000*Im_a1_0^2 + 0.250000000000000*Re_a1_0^2 + 0.125000000000000*Im_a1_1^2 + 0.125000000000000*Re_a1_1^2 + 0.0625000000000000*Im_a1_2^2 + 0.0625000000000000*Re_a1_2^2 + 0.0312500000000000*Im_a1_3^2 + 0.0312500000000000*Re_a1_3^2 + 0.0156250000000000*Im_a1_4^2 + 0.0156250000000000*Re_a1_4^2 + 0.00781250000000001*Im_a1_5^2 + 0.00781250000000001*Re_a1_5^2 + 0.00390625000000000*Im_a1_6^2 + 0.00390625000000000*Re_a1_6^2 + 0.00195312500000000*Im_a1_7^2 + 0.00195312500000000*Re_a1_7^2 + 0.000976562500000001*Im_a1_8^2 + 0.000976562500000001*Re_a1_8^2 + 0.000488281250000001*Im_a1_9^2 + 0.000488281250000001*Re_a1_9^2

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
            sage: f = 3*C.real(0, 0)^2 + 5*C.imag(0, 0)^2 + 7*C.real(1, 0)^2 + 11*C.imag(1, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [1.00000000000000], 1: [-1.00000000000000]}, imag={}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.00000000000000], 1: [-1.00000000000000]}, lagrange=[], value=-0.000000000000000),
             PowerSeriesConstraints.Constraint(real={0: [6.00000000000000]}, imag={}, lagrange=[1.00000000000000, 0, 1.00000000000000], value=0.000000000000000),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [10.0000000000000]}, lagrange=[0, 1.00000000000000, 0, 1.00000000000000], value=0.000000000000000),
             PowerSeriesConstraints.Constraint(real={1: [14.0000000000000]}, imag={}, lagrange=[-1.00000000000000, 0, -1.00000000000000], value=0.000000000000000), 
             PowerSeriesConstraints.Constraint(real={}, imag={1: [22.0000000000000]}, lagrange=[0, -1.00000000000000, 0, -1.00000000000000], value=0.000000000000000)]

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

                    def nth(L, n, default):
                        return (L[n:n+1] or [default])[0]

                    lagrange = [nth(getattr(g[i], part).get(triangle, []),k, 0) for i in range(lagranges)]

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
