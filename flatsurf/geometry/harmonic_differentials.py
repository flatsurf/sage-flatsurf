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
            (0 - 1*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
             0 - 1*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

            sage: η + η  # tol 1e-6
            (0 - 2*I + (0 + 0*I)*z0 + (0 + 0*I)*z0^2 + (0 + 0*I)*z0^3 + (0 + 0*I)*z0^4 + (0 + 0*I)*z0^5 + (0 + 0*I)*z0^6 + (0 + 0*I)*z0^7 + (0 + 0*I)*z0^8 + (0 + 0*I)*z0^9 + O(z0^10),
             0 - 2*I + (0 + 0*I)*z1 + (0 + 0*I)*z1^2 + (0 + 0*I)*z1^3 + (0 + 0*I)*z1^4 + (0 + 0*I)*z1^5 + (0 + 0*I)*z1^6 + (0 + 0*I)*z1^7 + (0 + 0*I)*z1^8 + (0 + 0*I)*z1^9 + O(z1^10))

        """
        # TODO: Some of the imaginary parts of the above output are not correct.

        return self.parent()({
            triangle: self._series[triangle] + other._series[triangle]
            for triangle in self._series
        })

    @staticmethod
    def _midpoint(surface, triangle, edge):
        r"""
        Return the (complex) vector from the center of the circumcircle to
        the center of the edge.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import HarmonicDifferential
            sage: HarmonicDifferential._midpoint(T, 0, 0)
            0j
            sage: HarmonicDifferential._midpoint(T, 0, 1)
            0.5j
            sage: HarmonicDifferential._midpoint(T, 0, 2)
            (-0.5+0j)
            sage: HarmonicDifferential._midpoint(T, 1, 0)
            0j
            sage: HarmonicDifferential._midpoint(T, 1, 1)
            -0.5j
            sage: HarmonicDifferential._midpoint(T, 1, 2)
            (0.5+0j)

        """
        P = surface.polygon(triangle)
        Δ = -P.circumscribing_circle().center() + P.vertex(edge) + P.edge(edge) / 2
        return complex(*Δ)

    def evaluate(self, triangle, Δ, derivative=0):
        C = PowerSeriesConstraints(self.parent().surface(), self.precision())
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
            sage: η._evaluate(C.gen(0, 0) + C.gen(1, 0))  # tol 1e-6
            0 - 2*I

        """
        coefficients = {}

        C = PowerSeriesConstraints(self.parent().surface(), self.precision())

        for gen in expression.variables():
            kind, triangle, k = C._describe_generator(gen)
            coefficient = self._series[triangle][k]

            if kind == "gen":
                coefficients[gen] = coefficient
            elif kind == "real":
                coefficients[gen] = coefficient.real()
            else:
                assert kind == "imag"
                coefficients[gen] = coefficient.imag()

        value = expression.parent()(expression.substitute(coefficients))
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
        C = PowerSeriesConstraints(self.parent().surface(), self.precision())
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

    def _element_from_cohomology(self, cocycle, /, prec=10, consistency=3):
        # We develop a consistent system of Laurent series at each vertex of the Voronoi diagram
        # to describe a differential.

        # Let η be the differential we are looking for. To describe η we will use ω, the differential
        # corresponding to the flat structure given by our triangulation on this Riemann surface.
        # Then f:=η/ω is a meromorphic function which we can develop locally into a Laurent series.
        # Away from the vertices of the triangulation, ω has no zeros, so f has no poles there and is
        # thus given by a power series.

        # At each vertex of the Voronoi diagram, write f=Σ a_k x^k + O(x^prec). Our task is now to determine
        # the a_k.

        constraints = PowerSeriesConstraints(self.surface(), prec=prec)

        # (1) The radius of convergence of the power series is the distance from the vertex of the Voronoi
        # cell to the closest vertex of the triangulation (since we use a Delaunay triangulation, all vertices
        # are at the same distance in fact.) So the radii of convergence of two neigbhouring cells overlap
        # and the power series must coincide there. Note that this constraint is unrelated to the cohomology
        # class Φ.
        constraints.require_consistency(consistency)

        # TODO: What should the rank be after this step?
        # from numpy.linalg import matrix_rank
        # A = constraints.matrix()[0]
        # print(len(A), len(A[0]), matrix_rank(A))

        # (2) We have that for any cycle γ, Re(∫fω) = Re(∫η) = Φ(γ). We can turn this into constraints
        # on the coefficients as we integrate numerically following the path γ as it intersects the radii of
        # convergence.
        constraints.require_cohomology(cocycle)

        # (3) Since the area 1 /(2iπ) ∫ η \wedge \overline{η} must be finite [TODO: REFERENCE?] we optimize for
        # this quantity to be minimal.
        constraints.require_finite_area()

        η = self.element_class(self, constraints.solve())

        # Check whether this is actually a global differential:
        # (1) Check that the series are actually consistent where the Voronoi cells overlap.

        def check(actual, expected, message, abs_error_bound = 1e-9, rel_error_bound = 1e-6):
            abs_error = abs(expected - actual)
            if abs_error > abs_error_bound:
                if expected == 0 or abs_error / abs(expected) > rel_error_bound:
                    print(f"{message}; expected: {expected}, got: {actual}")

        for (triangle, edge) in self._surface.edge_iterator():
            triangle_, edge_ = self._surface.opposite_edge(triangle, edge)
            for derivative in range(consistency):
                expected = η.evaluate(triangle, HarmonicDifferential._midpoint(self._surface, triangle, edge), derivative)
                other = η.evaluate(triangle_, HarmonicDifferential._midpoint(self._surface, triangle_, edge_), derivative)
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

        def get(self, gen):
            r"""
            Return the coefficients that are multiplied with the coefficient
            ``gen`` of the power series.
            """
            gen = str(gen)

            # TODO: This use of magic strings is not great.
            if gen.startswith("Re_a"):
                coefficients = self.real
            elif gen.startswith("Im_a"):
                coefficients = self.imag
            else:
                raise NotImplementedError

            gen = gen[4:]
            triangle, k = gen.split('_')
            triangle = int(triangle)
            k = int(k)

            coefficients = coefficients.get(triangle, [])[k:k+1]
            if not coefficients:
                from sage.all import ZZ
                return ZZ(0)
            return coefficients[0]

    def __init__(self, surface, prec):
        self._surface = surface
        self._prec = prec
        self._constraints = []

    def __repr__(self):
        return repr(self._constraints)

    @cached_method
    def symbolic_ring(self, triangle=None):
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
            sage: C.symbolic_ring(triangle=1)
            Multivariate Polynomial Ring in a1_0, a1_1, a1_2, Re_a1_0, Re_a1_1, Re_a1_2, Im_a1_0, Im_a1_1, Im_a1_2 over Complex Field with 53 bits of precision

            sage: C.symbolic_ring()
            Multivariate Polynomial Ring in a0_0, a0_1, a0_2, Re_a0_0, Re_a0_1, Re_a0_2, Im_a0_0, Im_a0_1, Im_a0_2, a1_0, a1_1, a1_2, Re_a1_0, Re_a1_1, Re_a1_2, Im_a1_0, Im_a1_1, Im_a1_2 over Complex Field with 53 bits of precision

        """
        gens = []

        for t in [triangle] if triangle else self._surface.label_iterator():
            gens += [f"a{t}_{n}" for n in range(self._prec)]
            gens += [f"Re_a{t}_{n}" for n in range(self._prec)]
            gens += [f"Im_a{t}_{n}" for n in range(self._prec)]

        # TODO: Should we use a better/configured base ring here?

        from sage.all import PolynomialRing, CC
        return PolynomialRing(CC, gens)

    @cached_method
    def gen(self, triangle, k):
        r"""
        Return the kth generator of the :meth:`symbolic_ring` for ``triangle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.gen(0, 0)
            a0_0
            sage: C.gen(0, 1)
            a0_1
            sage: C.gen(1, 2)
            a1_2

        """
        if k >= self._prec:
            raise ValueError("symbolic ring has no k-th generator")
        return self.symbolic_ring(triangle).gen(k)

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
        return self.symbolic_ring(triangle).gen(self._prec + k)

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
        return self.symbolic_ring(triangle).gen(2*self._prec + k)

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
            sage: C._describe_generator(C.gen(0, 0))
            ('gen', 0, 0)
            sage: C._describe_generator(C.imag(1, 2))
            ('imag', 1, 2)
            sage: C._describe_generator(C.real(2, 1))
            ('real', 2, 1)

        """
        gen = str(gen)
        if gen.startswith("a"):
            kind = "gen"
            gen = gen[1:]
        elif gen.startswith("Re_a"):
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

        # TODO: Is there a more generic ring than RR?
        from sage.all import RR
        if x in RR:
            if part == "real":
                # If this is just a constant, return it.
                return RR(x)
            elif part == "image":
                return RR.zero()

            assert False  # unreachable

        # Eliminate the generators a_k by rewriting them as Re(a_k) + I*Im(a_k)
        substitutions = {}
        for gen in x.variables():
            kind, triangle, k = self._describe_generator(gen)

            if kind != "gen":
                continue

            real = self.real(triangle, k)
            imag = self.imag(triangle, k)
            imag *= imag.parent().base_ring().gen()
            substitutions[gen] = real + imag

        if substitutions:
            x = x.substitute(substitutions)

        terms = []

        # We use Re(c*Re(a_k)) = Re(c) * Re(a_k) and Re(c*Im(a_k)) = Re(c) * Im(a_k)
        # and Im(c*Re(a_k)) = Im(c) * Re(a_k) and Im(c*Im(a_k)) = Im(c) * Im(a_k), respectively.
        for gen in x.variables():
            if x.degree(gen) > 1:
                raise NotImplementedError

            kind, triangle, k = self._describe_generator(gen)

            assert kind != "gen"

            terms.append(self.project(x[gen], part) * gen)

        terms.append(self.project(x.constant_coefficient(), part))

        return sum(terms)

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
    def _formal_power_series(self, triangle):
        from sage.all import PowerSeriesRing
        R = PowerSeriesRing(self.symbolic_ring(triangle=triangle), 'z')

        return R([self.gen(triangle, n) for n in range(self._prec)])

    def develop(self, triangle, Δ=0):
        r"""
        Return the power series obtained by developing at z + Δ.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, prec=3)
            sage: C.develop(0)
            a0_0 + a0_1*z + a0_2*z^2
            sage: C.develop(1, 1)  # tol 1e-9
            a1_0 + a1_1 + a1_2 + (a1_1 + 2*a1_2)*z + a1_2*z^2

        """
        # TODO: Check that Δ is within the radius of convergence.
        f = self._formal_power_series(triangle)
        return f(f.parent().gen() + Δ)

    def integrate(self, cycle):
        r"""
        Return the linear combination of the power series coefficients that
        decsribe the integral of a differential along the homology class
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
            sage: C.integrate(a)  # tol 1e-6
            (0.500000000000000 + 0.500000000000000*I)*a0_0 + (-0.250000000000000)*a0_1 + (0.0416666666666667 - 0.0416666666666667*I)*a0_2 + (0.00625000000000000 + 0.00625000000000000*I)*a0_4 + (0.500000000000000 + 0.500000000000000*I)*a1_0 + 0.250000000000000*a1_1 + (0.0416666666666667 - 0.0416666666666667*I)*a1_2 + (0.00625000000000000 + 0.00625000000000000*I)*a1_4
            sage: C.integrate(b)  # tol 1e-6
            (-0.500000000000000)*a0_0 + 0.125000000000000*a0_1 + (-0.0416666666666667)*a0_2 + 0.0156250000000000*a0_3 + (-0.00625000000000000)*a0_4 + (-0.500000000000000)*a1_0 + (-0.125000000000000)*a1_1 + (-0.0416666666666667)*a1_2 + (-0.0156250000000000)*a1_3 + (-0.00625000000000000)*a1_4

        """
        surface = cycle.surface()

        expression = 0

        for path, multiplicity in cycle.voronoi_path().monomial_coefficients().items():

            for S, T in zip((path[-1],) + path, path):
                # Integrate from the midpoint of the edge of S to the midpoint of the edge of T
                S = surface.opposite_edge(*S)
                assert S[0] == T[0], f"consecutive elements of a path must be attached to the same face in {path} but {S} and {T} do not have that property"

                # Namely we integrate the power series defined around the Voronoi vertex of S by symbolically integrating each monomial term.

                # The midpoints of the edges
                P = HarmonicDifferential._midpoint(surface, *S)
                Q = HarmonicDifferential._midpoint(surface, *T)

                for k in range(self._prec):
                    expression -= self.gen(S[0], k) * multiplicity * P**(k + 1) / (k + 1)
                    expression += self.gen(T[0], k) * multiplicity * Q**(k + 1) / (k + 1)

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
            a0_0
            sage: C.evaluate(1, 0)
            a1_0
            sage: C.evaluate(1, 2)  # tol 1e-9
            a1_0 + 2*a1_1 + 4*a1_2

        """
        # TODO: Check that Δ is within the radius of convergence.

        if derivative >= self._prec:
            raise ValueError

        from sage.all import factorial
        return self.develop(triangle=triangle, Δ=Δ)[derivative] * factorial(derivative)

    def require_consistency(self, derivatives):
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
        triangulation. Therefore, it would be better to just require all
        coefficients to be identical in the first place::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_consistency(1)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1], 1: [-1]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1], 1: [-1]}, lagrange=[], value=0)]

        If we add more coefficients, we get three pairs of contraints for the
        three edges surrounding a face; for the edge on which the centers of
        the Voronoi cells fall, we get a more restrictive constraint forcing
        all the coefficients to be equal (these are the first four constraints
        recorded below.)::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_consistency(1)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 1.0], 1: [0.0, -1.0]}, imag={}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [0.0, 1.0], 1: [0.0, -1.0]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [0.0, -0.50], 1: [0.0, -0.50]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 0.50], 1: [0.0, 0.50]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={0: [1.0, -0.50], 1: [-1.0, -0.50]}, imag={}, lagrange=[], value=-0.0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.50], 1: [-1.0, -0.50]}, lagrange=[], value=-0.0)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_consistency(2)
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0, 1.0], 1: [0, -1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [0, 1.0], 1: [0, -1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [-0.0, -0.5], 1: [0.0, -0.5]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 0.5], 1: [-0.0, 0.5]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0, -0.5], 1: [-1.0, -0.5]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.5], 1: [-1.0, -0.5]}, lagrange=[], value=0)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            Δ0 = HarmonicDifferential._midpoint(self._surface, triangle0, edge0)
            Δ1 = HarmonicDifferential._midpoint(self._surface, triangle1, edge1)

            if abs(Δ0) < 1e-6 and abs(Δ1) < 1e-6:
                # Force power series to be identical if the Delaunay triangulation is ambiguous at this edge.
                for k in range(self._prec):
                    self.add_constraint(self.gen(triangle0, k) - self.gen(triangle1, k))

                continue

            # Require that the 0th, ..., derivatives-1th derivatives are the same at the midpoint of the edge.
            for derivative in range(derivatives):
                self.add_constraint(
                    self.evaluate(triangle0, Δ0, derivative) - self.evaluate(triangle1, Δ1, derivative))

    def require_L2_consistency(self):
        r"""
        For each pair of adjacent triangles meeting at and edge `e`, let `v` be
        the midpoint of `e`. We develop the power series coming from both
        triangles around that midpoint and check them for agreement. Namely, we
        integrate the square of their difference on the circle of maximal
        radius around `v`. (If the power series agree, that integral should be
        zero.)

        TODO

        """
        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            # The midpoint of the edge where the triangles meet with respect to
            # the center of the triangle.
            Δ0 = HarmonicDifferential._midpoint(self._surface, triangle0, edge0)
            Δ1 = HarmonicDifferential._midpoint(self._surface, triangle1, edge1)

            # Develop both power series around that midpoint, i.e., Taylor expand them.
            T0 = HarmonicDifferential._taylor_symbolic(Δ0, self._prec)
            T1 = HarmonicDifferential._taylor_symbolic(Δ0, self._prec)

            # Write b_n for the difference of the n-th coefficient of both power series.
            # b = [...]
            # We want to minimize the sum of b_n^2 R^n where R is half the
            # length of the edge we are on.
            # Taking a square root of R^n, we merge R^n with b_n^2.
            raise NotImplementedError
            # To optimize this error, write it as Lagrange multipliers.
            raise NotImplementedError # TODO: We need some generic infrastracture to merge quadratic conditions such as this and require_finite_area by weighing both.

    def _area(self):
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
            sage: area = PowerSeriesConstraints(T, η.precision())._area()
            sage: area  # tol 2e-2
            0.250*Re_a0_0^2 + 0.176*Re_a0_1^2 + 0.125*Re_a0_2^2 + 0.0883*Re_a0_3^2 + 0.0625*Re_a0_4^2 + 0.0441*Re_a0_5^2 + 0.0312*Re_a0_6^2 + 0.0220*Re_a0_7^2 + 0.0156*Re_a0_8^2 + 0.0110*Re_a0_9^2
            + 0.250*Im_a0_0^2 + 0.176*Im_a0_1^2 + 0.125*Im_a0_2^2 + 0.0883*Im_a0_3^2 + 0.0625*Im_a0_4^2 + 0.0441*Im_a0_5^2 + 0.0312*Im_a0_6^2 + 0.0220*Im_a0_7^2 + 0.0156*Im_a0_8^2 + 0.0110*Im_a0_9^2
            + 0.250*Re_a1_0^2 + 0.176*Re_a1_1^2 + 0.125*Re_a1_2^2 + 0.0883*Re_a1_3^2 + 0.0625*Re_a1_4^2 + 0.0441*Re_a1_5^2 + 0.0312*Re_a1_6^2 + 0.0220*Re_a1_7^2 + 0.0156*Re_a1_8^2 + 0.0110*Re_a1_9^2
            + 0.250*Im_a1_0^2 + 0.176*Im_a1_1^2 + 0.125*Im_a1_2^2 + 0.0883*Im_a1_3^2 + 0.0625*Im_a1_4^2 + 0.0441*Im_a1_5^2 + 0.0312*Im_a1_6^2 + 0.0220*Im_a1_7^2 + 0.0156*Im_a1_8^2 + 0.0110*Im_a1_9^2


        The correct area would be 1/2π here. However, we are overcounting
        because we sum the single Voronoi cell twice. And also, we approximate
        the square with a circle for another factor π/2::

            sage: η._evaluate(area)  # tol 1e-2
            .5

        """
        # To make our lives easier, we do not optimize this value but instead
        # half the sum of the |a_k|^2·radius^(k+2) = (Re(a_k)^2 +
        # Im(a_k)^2)·radius^(k+2) which is a very rough upper bound for the
        # area.
        area = 0

        for triangle in range(self._surface.num_polygons()):
            R = float(self._surface.polygon(triangle).circumscribing_circle().radius_squared().sqrt())

            for k in range(self._prec):
                area += (self.real(triangle, k)**2 + self.imag(triangle, k)**2) * R**(k+2)

        return area/2

    def require_finite_area(self):
        r"""
        Since the area 1 /(2iπ) ∫ η \wedge \overline{η} must be finite [TODO:
        REFERENCE?] we can optimize for this quantity to be minimal.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

        ::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_consistency(1)

            sage: C.require_finite_area()
            sage: C  # tol 1e-9
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0]}, imag={}, lagrange=[-1.0], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0]}, lagrange=[0, -1.0], value=0),
             PowerSeriesConstraints.Constraint(real={1: [1.0]}, imag={}, lagrange=[1.0], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={1: [1.0]}, lagrange=[0, 1.0], value=0)]

        """
        self.optimize(self._area())

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
            sage: C.require_consistency(1)
            sage: f = 3*C.real(0, 0)^2 + 5*C.imag(0, 0)^2 + 7*C.real(1, 0)^2 + 11*C.imag(1, 0)^2
            sage: C.optimize(f)
            sage: C
            ...
            PowerSeriesConstraints.Constraint(real={0: [6.00000000000000]}, imag={}, lagrange=[-1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={}, imag={0: [10.0000000000000]}, lagrange=[0, -1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={1: [14.0000000000000]}, imag={}, lagrange=[1.00000000000000], value=-0.000000000000000),
            PowerSeriesConstraints.Constraint(real={}, imag={1: [22.0000000000000]}, lagrange=[0, 1.00000000000000], value=-0.000000000000000)]

        """
        # We cannot optimize if there is an unbound z in the expression.
        f = self.symbolic_ring()(f)

        # We rewrite a_k as Re(a_k) + i Im(a_k).
        for triangle in range(self._surface.num_polygons()):
            for k in range(self._prec):
                a_k = self.gen(triangle, k)
                if f.degree(a_k):
                    raise NotImplementedError("cannot rewrite a_k as Re(a_k) + i Im(a_k) yet")

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
                for gen in [self.real(triangle, k), self.imag(triangle, k)]:
                    if f.degree(gen) <= 0:
                        continue

                    gen = f.parent()(gen)

                    self.add_constraint(f.derivative(gen), lagrange=[-g[i].get(gen) for i in range(lagranges)], value=ZZ(0))

        # We form the partial derivatives with respect to the λ_i. This yields
        # the condition -g_i=0 which is already recorded in the linear system.

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
        A = []
        b = []

        lagranges = max(len(constraint.lagrange) for constraint in self._constraints)

        for constraint in self._constraints:
            A.append([])
            b.append(constraint.value)
            for triangle in self._surface.label_iterator():
                real = constraint.real.get(triangle, [])
                real.extend([ZZ(0)] * (self._prec - len(real)))
                A[-1].extend(real)

                imag = constraint.imag.get(triangle, [])
                imag.extend([ZZ(0)] * (self._prec - len(imag)))
                A[-1].extend(imag)

            lagrange = constraint.lagrange
            lagrange.extend([ZZ(0)] * (lagranges - len(lagrange)))

            A[-1].extend(lagrange)

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
        A, b = self.matrix()

        import scipy.linalg
        import numpy
        solution, residues, _, _ = scipy.linalg.lstsq(numpy.matrix(A), numpy.array(b))

        lagranges = max(len(constraint.lagrange) for constraint in self._constraints)

        if lagranges:
            solution = solution[:-lagranges]

        solution = [solution[2*k*self._prec:2*(k+1)*self._prec] for k in range(self._surface.num_polygons())]

        from sage.all import CC
        return {
            triangle: HarmonicDifferentials(self._surface).power_series_ring(triangle)([CC(solution[k][p], solution[k][p + self._prec]) for p in range(self._prec)]).add_bigoh(self._prec)
            for (k, triangle) in enumerate(self._surface.label_iterator())
        }
