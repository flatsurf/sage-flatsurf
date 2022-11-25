r"""
TODO: Document this module.

EXAMPLES:

We compute harmonic differentials on the square torus::

    sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
    sage: T.set_immutable()

    sage: H = SimplicialHomology(T)
    sage: a, b = H.gens()

First, the harmonic differentials that sends the diagonal `a` to 1 and the negative
horizontal `b` to zero. Note that integrating the differential given by the
power series Σa_n z^n along `a` yields `Re(a_0) - Im(a_0)` and along `b` we get
`-Re(a_0)`. Therefore, this must have Re(a_0) = 0 and Im(a_0) = -1::

    sage: H = SimplicialCohomology(T)
    sage: f = H({a: 1})
    sage: Ω = HarmonicDifferentials(T)
    sage: Ω(f)
    (-1.000000000000000000000*I + O(z0^10), -1.000000000000000000000*I + O(z1^10))

The harmonic differential that integrates as 0 along `a` but 1 along `b` must
similarly have Re(a_0) = -1 but Im(a_0) = -1::

    sage: g = H({b: 1})
    sage: Ω(g)
    (-1.000000000000000000000 - 1.000000000000000000000*I + O(z0^10), -1.000000000000000000000 - 1.000000000000000000000*I + O(z1^10))

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
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement


class HarmonicDifferential(Element):
    def __init__(self, parent, series, residue=None, cocycle=None):
        super().__init__(parent)
        self._series = series
        self._residue = residue
        self._cocycle = cocycle

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

            sage: Ω(f) + Ω(f) + Ω(H({a: -2}))
            (O(z0^10), O(z1^10))

        """
        return self.parent()({
            triangle: self._series[triangle] + other._series[triangle]
            for triangle in self._series
        })

    def error(self, kind=None, verbose=False, abs_tol=1e-6, rel_tol=1e-6):
        r"""
        Return whether this differential is likely inaccurate.

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
            sage: η.error()
            False

        """
        error = False

        def errors(expected, actual):
            abs_error = abs(expected - actual)
            rel_error = 0
            if abs(expected) > 1e-12:
                rel_error = abs_error / abs(expected)

            return abs_error, rel_error

        if kind is None or "residue" in kind:
            if self._residue is not None:
                report = f"Harmonic differential created by solving Ax=b with |Ax-b| = {self._residue}."
                if verbose:
                    print(report)
                if self._residue > abs_tol:
                    error = report
                    if not verbose:
                        return error

        if kind is None or "cohomology" in kind:
            if self._cocycle is not None:
                for gen in self._cocycle.parent().homology().gens():
                    expected = self._cocycle(gen)
                    actual = self.integrate(gen).real
                    if callable(actual):
                        actual = actual()

                    abs_error, rel_error = errors(expected, actual)

                    report = f"Integrating along cycle gives {actual} whereas the cocycle gave {expected}, i.e., an absolute error of {abs_error} and a relative error of {rel_error}."
                    if verbose:
                        print(report)

                    if abs_error > abs_tol or rel_error > rel_tol:
                        error = report
                        if not verbose:
                            return error

        if kind is None or "midpoint_derivatives" in kind:
            C = PowerSeriesConstraints(self.parent().surface(), self.precision())
            for (triangle, edge) in self.parent().surface().edge_iterator():
                triangle_, edge_ = self.parent().surface().opposite_edge(triangle, edge)
                for derivative in range(self.precision()//3):
                    expected = self.evaluate(triangle, C.complex_field()(*self.parent()._geometry.midpoint(triangle, edge)), derivative)
                    other = self.evaluate(triangle_, C.complex_field()(*self.parent()._geometry.midpoint(triangle_, edge_)), derivative)

                    abs_error, rel_error = errors(expected, other)

                    if abs_error > abs_tol or rel_error > rel_tol:
                        report = f"Power series defining harmonic differential are not consistent where triangles meet. {derivative}th derivative does not match between {(triangle, edge)} where it is {expected} and {(triangle_, edge_)} where it is {other}, i.e., there is an absolute error of {abs_error} and a relative error of {rel_error}."
                        if verbose:
                            print(report)

                        error = report
                        if not verbose:
                            return error

        if kind is None or "area" in kind:
            if verbose:
                C = PowerSeriesConstraints(self.parent().surface(), self.precision())
                area = self._evaluate(C._area_upper_bound())

                report = f"Area (upper bound) is {area}."
                print(report)

        if kind is None or "L2" in kind:
            C = PowerSeriesConstraints(self.parent().surface(), self.precision())
            abs_error = self._evaluate(C._L2_consistency())

            report = f"L2 norm of differential is {abs_error}."
            if verbose:
                print(report)

            if abs_error > abs_tol:
                error = report
                if not verbose:
                    return error

        return error

    def series(self, triangle):
        r"""
        Return the power series describing this harmonic differential at the
        center of the circumcircle of the given Delaunay ``triangle``.

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

            sage: η.series(0)
            -1.000000000000000000000*I + O(z0^10)

        """
        return self._series[triangle]

    def evaluate(self, triangle, Δ, derivative=0):
        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        return self._evaluate(C.evaluate(triangle, Δ, derivative=derivative))

    def roots(self):
        r"""
        Return the roots of this harmonic differential.

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
            sage: η.roots()
            []

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.regular_octagon().delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a = H.gens()[0]
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f, prec=16, check=False)  # TODO: There is a small discontinuity (rel error 1e-6.)
            sage: η.roots()  # TODO: Where should the roots be?
            []

        """
        roots = []

        surface = self.parent().surface()

        for triangle in surface.label_iterator():
            series = self._series[triangle]
            series = series.truncate(series.prec())
            for (root, multiplicity) in series.roots():
                if multiplicity != 1:
                    raise NotImplementedError

                root += root.parent()(*self.parent()._geometry.center(triangle))

                from sage.all import vector
                root = vector(root)

                # TODO: Keep roots that are within the circumcircle for sanity checking.
                # TODO: Make sure that roots on the edges are picked up by exactly one triangle.

                if not surface.polygon(triangle).contains_point(root):
                    continue

                roots.append((triangle, vector(root)))

        # TODO: Deduplicate roots.
        # TODO: Compute roots at the vertices.

        return roots

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
            sage: η._evaluate(R(C.gen(0, 0)) + R(C.gen(1, 0)))  # tol 1e-9
            0 - 2.0000000000000000*I

        """
        coefficients = {}

        for variable in expression.variables():
            kind, triangle, k = variable.describe()

            coefficient = self._series[triangle][k]

            if kind == "real":
                coefficients[variable] = coefficient.real()
            else:
                assert kind == "imag"
                coefficients[variable] = coefficient.imag()

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

    def cauchy_residue(self, vertex, n, angle=None):
        r"""
        Return the n-th coefficient of the power series expansion around the
        ``vertex`` in local coordinates of that vertex.

        Let ``d`` be the degree of the vertex and pick a (d+1)-th root z' of z such
        that this differential can be written as a power series Σb_n z'^n. This
        method returns the `b_n` by evaluating the Cauchy residue formula,
        i.e., 1/2πi ∫ f/z'^{n+1} integrating along a loop around the vertex.

        EXAMPLES:

        For the square torus, the vertices are no singularities, so harmonic
        differentials have no pole at the vertex::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: vertex = [(0, 1), (1, 0), (0, 2), (1, 1), (0, 0), (1, 2)]
            sage: angle = 1

            sage: η.cauchy_residue(vertex, 0, 1)  # tol 1e-9
            0 - 1.0*I
            sage: abs(η.cauchy_residue(vertex, -1, 1)) < 1e-9
            True
            sage: abs(η.cauchy_residue(vertex, -2, 1)) < 1e-9
            True
            sage: abs(η.cauchy_residue(vertex, 1, 1)) < 1e-9
            True
            sage: abs(η.cauchy_residue(vertex, 2, 1)) < 1e-9
            True

        """
        surface = self.parent().surface()

        if angle is None:
            vertex = set(vertex)
            for angle, v in surface.angles(return_adjacent_edges=True):
                if vertex & set(v):
                    assert vertex == set(v)
                    vertex = v
                    break
            else:
                raise ValueError("not a vertex in this surface")

        parts = []

        complex_field = self._series[0].parent().base_ring()

        # Integrate real & complex part independently.
        for part in ["real", "imag"]:
            integral = 0

            # TODO print(f"Building {part} part")

            # We integrate along a path homotopic to a (counterclockwise) unit
            # circle centered at the vertex.

            # We keep track of our last choice of a (d+1)-st root and pick the next
            # root that is closest while walking counterclockwise on that circle.
            arg = 0

            def root(z):
                if angle == 1:
                    return complex_field(z)

                roots = complex_field(z).nth_root(angle, all=True)

                positives = [root for root in roots if (root.arg() - arg) % (2*3.14159265358979) > -1e-6]

                return min(positives)

            for triangle, edge in vertex:
                # We integrate on the line segment from the midpoint of edge to
                # the midpoint of the previous edge in triangle, i.e., the next
                # edge in counterclockwise order walking around the vertex.
                edge_ = (edge - 1) % 3

                # TODO print(f"integrating across {triangle} from the midpoint of {edge} to {edge_}")

                P = complex_field(*self.parent()._geometry.midpoint(triangle, edge))
                Q = complex_field(*self.parent()._geometry.midpoint(triangle, edge_))

                # The vector from z=0, the center of the Voronoi cell, to the vertex.
                δ = P - complex_field(*surface.polygon(triangle).edge(edge)) / 2

                def at(t):
                    z = (1 - t) * P + t * Q
                    root_at_z = root(z - δ)
                    return z, root_at_z

                root_at_P = at(0)[1]
                arg = root_at_P.arg()

                # TODO print(f"in terms of the Voronoi cell midpoint, this is integrating from {P} to {Q} which lift to {at(0)[1]} and {at(1)[1]} relative to the vertex which is at {δ} from the center of the Voronoi cell")

                def integrand(t):
                    t = complex_field(t)
                    z, root_at_z = at(t)
                    denominator = root_at_z ** (n + 1)
                    numerator = self.evaluate(triangle, z)
                    integrand = numerator / denominator * (Q - P)
                    # TODO print(f"evaluating at {z} resp its root {root_at_z} produced ({numerator}) / ({denominator}) * ({Q - P}) = {integrand}")
                    integrand = getattr(integrand, part)()
                    return integrand

                from sage.all import numerical_integral
                integral_on_segment, error = numerical_integral(integrand, 0, 1)
                # TODO: What should be do about the error?
                integral += integral_on_segment

                root_at_Q = at(1)[1]
                arg = root_at_Q.arg()

            parts.append(integral)

        return complex_field(*parts) / complex_field(0, 2*3.14159265358979)

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
        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
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
    def __classcall__(cls, surface, category=None):
        r"""
        Normalize parameters when creating the space of harmonic differentials.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: HarmonicDifferentials(T) is HarmonicDifferentials(T)
            True

        """
        return super().__classcall__(cls, surface, category or SetsWithPartialMaps())

    def __init__(self, surface, category):
        if surface != surface.delaunay_triangulation():
            raise NotImplementedError("Surface must be Delaunay triangulated")

        Parent.__init__(self, category=category)

        self._surface = surface

        self._geometry = GeometricPrimitives(surface)

    def surface(self):
        return self._surface

    def _repr_(self):
        return f"Ω({self._surface})"

    def _element_constructor_(self, x, *args, **kwargs):
        if not x:
            η = self.element_class(self, None, *args, **kwargs)

        if isinstance(x, dict):
            return self.element_class(self, x, *args, **kwargs)

        return self._element_from_cohomology(x, *args, **kwargs)

    def _element_from_cohomology(self, cocycle, /, prec=10, algorithm=["L2"], check=True):
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

        solution, residue = constraints.solve()
        η = self.element_class(self, solution, residue=residue, cocycle=cocycle)

        if check:
            if report := η.error():
                raise ValueError(report)

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
            (0, 0)
            sage: G.midpoint(0, 1)
            (0, 1/2)
            sage: G.midpoint(0, 2)
            (-1/2, 0)
            sage: G.midpoint(1, 0)
            (0, 0)
            sage: G.midpoint(1, 1)
            (0, -1/2)
            sage: G.midpoint(1, 2)
            (1/2, 0)

        """
        polygon = self._surface.polygon(triangle)
        return -self.center(triangle) + polygon.vertex(edge) + polygon.edge(edge) / 2

    @cached_method
    def center(self, triangle):
        return self._surface.polygon(triangle).circumscribing_circle().center()


class SymbolicCoefficientExpression(CommutativeRingElement):
    # TODO: Make sure that we never have zero coefficients as these would break degree computations.

    def __init__(self, parent, coefficients):
        super().__init__(parent)

        self._coefficients = coefficients

    def _richcmp_(self, other, op):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a == a
            True
            sage: a == b
            False

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not (self == other)

        if op == op_EQ:
            return self._coefficients == other._coefficients

        raise NotImplementedError

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a
            Im(a0,0)
            sage: b
            Re(a0,0)
            sage: a + b
            Re(a0,0) + Im(a0,0)
            sage: a + b + 1
            Re(a0,0) + Im(a0,0) + 1.00000000000000

        """
        if self.is_constant():
            return repr(self.constant_coefficient())

        def decode(gen):
            if gen < 0:
                return -gen-1,

            kind = "Im" if gen % 2 else "Re"
            gen //= 2
            polygon = gen % self.parent()._surface.num_polygons()
            gen //= self.parent()._surface.num_polygons()
            k = gen

            return kind, polygon, k

        def variable_name(gen):
            gen = decode(gen)

            if len(gen) == 1:
                return f"λ{gen[0]}"

            kind, polygon, k = gen
            return f"{kind}__open__a{polygon}__comma__{k}__close__"

        def key(gen):
            gen = decode(gen)

            if len(gen) == 1:
                n = gen[0]
                return 1e9, n

            kind, polygon, k = gen
            return polygon, k, 0 if kind == "Re" else 1

        gens = list({gen for monomial in self._coefficients.keys() for gen in monomial})
        gens.sort(key=key)

        variable_names = tuple(variable_name(gen) for gen in gens)

        from sage.all import PolynomialRing
        R = PolynomialRing(self.base_ring(), variable_names)

        def monomial(variables):
            monomial = R.one()
            for variable in variables:
                monomial *= R.gen(gens.index(variable))
            return monomial

        f = sum(coefficient * monomial(gens) for (gens, coefficient) in self._coefficients.items())

        return repr(f).replace('__open__', '(').replace('__close__', ')').replace('__comma__', ',')

    def degree(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.degree(a)
            1
            sage: (a + b).degree(a)
            1
            sage: (a * b + a).degree(a)
            1
            sage: R.one().degree(a)
            0
            sage: R.zero().degree(a)
            -1

        """
        if not gen.is_variable():
            raise ValueError(f"gen must be a monomial not {gen}")

        variable = next(iter(gen._coefficients))[0]

        return max([monomial.count(variable) for monomial in self._coefficients], default=-1)

    def is_monomial(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.is_monomial()
            True
            sage: (a + b).is_monomial()
            False
            sage: R.one().is_monomial()
            False
            sage: R.zero().is_monomial()
            False
            sage: (a * a).is_monomial()
            True

        """
        if len(self._coefficients) != 1:
            return False

        ((key, value),) = self._coefficients.items()

        return bool(key) and value.is_one()

    def is_constant(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.is_constant()
            False
            sage: (a + b).is_constant()
            False
            sage: R.one().is_constant()
            True
            sage: R.zero().is_constant()
            True
            sage: (a * a).is_constant()
            False

        """
        coefficients = len(self._coefficients)

        if coefficients == 0:
            return True

        if coefficients > 1:
            return False

        monomial = next(iter(self._coefficients.keys()))

        return not monomial

    def norm(self, p=2):
        r"""
        Return the p-norm of the coefficient vector.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: x = R.gen(('imag', 0, 0)) + 1; x
            Im(a0,0) + 1.00000000000000
            sage: x.norm(1)
            2.00000000000000
            sage: x.norm(oo)
            1.00000000000000

        """
        from sage.all import vector
        return vector(self._coefficients.values()).norm(p)

    def _neg_(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: -a
            -Im(a0,0)
            sage: -(a + b)
            -Re(a0,0) - Im(a0,0)
            sage: -(a * a)
            -Im(a0,0)^2
            sage: -R.one()
            -1.00000000000000
            sage: -R.zero()
            0.000000000000000

        """
        parent = self.parent()
        return type(self)(parent, {key: -coefficient for (key, coefficient) in self._coefficients.items()})

    def _add_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a + 1
            Im(a0,0) + 1.00000000000000
            sage: a + (-a)
            0.000000000000000
            sage: a + b
            Re(a0,0) + Im(a0,0)
            sage: a * a + b * b
            Re(a0,0)^2 + Im(a0,0)^2

        """
        parent = self.parent()

        if len(self._coefficients) < len(other._coefficients):
            self, other = other, self

        coefficients = dict(self._coefficients)

        for monomial, coefficient in other._coefficients.items():
            assert coefficient
            if monomial not in coefficients:
                coefficients[monomial] = coefficient
            else:
                coefficients[monomial] += coefficient

                if not coefficients[monomial]:
                    del coefficients[monomial]

        return type(self)(parent, coefficients)

    def _sub_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a - 1
            Im(a0,0) - 1.00000000000000
            sage: a - a
            0.000000000000000
            sage: a * a - b * b
            -Re(a0,0)^2 + Im(a0,0)^2

        """
        return self._add_(-other)

    def _mul_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a * a
            Im(a0,0)^2
            sage: a * b
            Re(a0,0)*Im(a0,0)
            sage: a * R.one()
            Im(a0,0)
            sage: a * R.zero()
            0.000000000000000
            sage: (a + b) * (a - b)
            -Re(a0,0)^2 + Im(a0,0)^2

        """
        parent = self.parent()

        if other.is_zero() or self.is_zero():
            return parent.zero()

        if other.is_one():
            return self

        if self.is_one():
            return other

        coefficients = {}

        for self_monomial, self_coefficient in self._coefficients.items():
            assert self_coefficient
            for other_monomial, other_coefficient in other._coefficients.items():
                assert other_coefficient

                monomial = tuple(sorted(self_monomial + other_monomial))
                coefficient = self_coefficient * other_coefficient

                if monomial not in coefficients:
                    coefficients[monomial] = coefficient
                else:
                    coefficients[monomial] += coefficient
                    if not coefficients[monomial]:
                        del coefficients[monomial]

        return type(self)(self.parent(), coefficients)

    def _rmul_(self, right):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: a * 0
            0.000000000000000
            sage: a * 1
            Im(a0,0)
            sage: a * 2
            2.00000000000000*Im(a0,0)

        """
        return self._lmul_(right)

    def _lmul_(self, left):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: 0 * a
            0.000000000000000
            sage: 1 * a
            Im(a0,0)
            sage: 2 * a
            2.00000000000000*Im(a0,0)

        """
        return type(self)(self.parent(), {key: coefficient for (key, value) in self._coefficients.items() if (coefficient := left * value)})

    def constant_coefficient(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.constant_coefficient()
            0.000000000000000
            sage: (a + b).constant_coefficient()
            0.000000000000000
            sage: R.one().constant_coefficient()
            1.00000000000000
            sage: R.zero().constant_coefficient()
            0.000000000000000

        """
        return self._coefficients.get((), self.parent().base_ring().zero())

    def variables(self, kind=None):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.variables()
            {Im(a0,0)}
            sage: (a + b).variables()
            {Im(a0,0), Re(a0,0)}
            sage: (a * a).variables()
            {Im(a0,0)}
            sage: (a * b).variables()
            {Im(a0,0), Re(a0,0)}
            sage: R.one().variables()
            set()
            sage: R.zero().variables()
            set()

        """
        if kind == "lagrange":
            def filter(gen):
                return gen < 0
        elif kind is None:
            def filter(gen):
                return True
        else:
            raise ValueError("unsupported kind")

        return set(self.parent()((gen,)) for monomial in self._coefficients.keys() for gen in monomial if filter(gen))

    def polygon(self):
        r"""
        Return the label of the polygon affected by this variable.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.polygon()
            0
            sage: b.polygon()
            0
            sage: (a*b).polygon()
            Traceback (most recent call last):
            ...
            ValueError: element must be a variable

        """
        if not self.is_variable():
            raise ValueError("element must be a variable")

        variable = next(iter(self._coefficients.keys()))[0]

        if variable < 0:
            raise ValueError("Lagrange multipliers are not associated to a polygon")

        return variable % (2*self.parent()._surface.num_polygons()) // 2

    def is_variable(self):
        r"""
        Return the label of the polygon affected by this variable.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.is_variable()
            True
            sage: R.zero().is_variable()
            False
            sage: R.one().is_variable()
            False
            sage: (a + 1).is_variable()
            False
            sage: (a*b).is_variable()
            False

        """
        if not self.is_monomial():
            return False

        monomial = next(iter(self._coefficients.keys()))

        if len(monomial) != 1:
            return False

        return True

    def describe(self):
        r"""
        Return a tuple describing the nature of this variable.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a.describe()
            ('imag', 0, 0)
            sage: b.describe()
            ('real', 0, 0)
            sage: (a + b).describe()
            Traceback (most recent call last):
            ...
            ValueError: element must be a variable

        """
        if not self.is_variable():
            raise ValueError("element must be a variable")

        variable = next(iter(self._coefficients.keys()))[0]

        if variable < 0:
            return ("lagrange", -variable-1)

        triangle = self.polygon()
        k = variable // (2*self.parent()._surface.num_polygons())
        if variable % 2:
            return ("imag", triangle, k)
        else:
            return ("real", triangle, k)

    def real(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: c = (a + b)**2
            sage: c.real()
            Re(a0,0)^2 + 2.00000000000000*Re(a0,0)*Im(a0,0) + Im(a0,0)^2

        """
        return self.map_coefficients(lambda c: c.real(), self.parent().change_ring(self.parent().real_field()))

    def imag(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: c = (a + b)**2
            sage: c.imag()
            0.000000000000000

            sage: c = (I*a + b)**2
            sage: c.imag()
            2.00000000000000*Re(a0,0)*Im(a0,0)

        """
        return self.map_coefficients(lambda c: c.imag(), self.parent().change_ring(self.parent().real_field()))

    def __getitem__(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a[a]
            1.00000000000000
            sage: a[b]
            0.000000000000000
            sage: (a + b)[a]
            1.00000000000000
            sage: (a * b)[a]
            0.000000000000000
            sage: (a * b)[a * b]
            1.00000000000000

        """
        if not gen.is_monomial():
            raise ValueError

        return self._coefficients.get(next(iter(gen._coefficients.keys())), self.parent().base_ring().zero())

    def __hash__(self):
        return hash(tuple(sorted(self._coefficients.items())))

    def total_degree(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: R.zero().total_degree()
            -1
            sage: R.one().total_degree()
            0
            sage: a.total_degree()
            1
            sage: (a * a + b).total_degree()
            2

        """
        degrees = [len(monomial) for monomial in self._coefficients]
        return max(degrees, default=-1)

    def derivative(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: R.zero().derivative(a)
            0.000000000000000
            sage: R.one().derivative(a)
            0.000000000000000
            sage: a.derivative(a)
            1.00000000000000
            sage: a.derivative(b)
            0.000000000000000
            sage: c = (a + b) * (a - b)
            sage: c.derivative(a)
            2.00000000000000*Im(a0,0)
            sage: c.derivative(b)
            -2.00000000000000*Re(a0,0)

        """
        if not gen.is_variable():
            raise ValueError

        gen = next(iter(gen._coefficients.keys()))[0]

        derivative = self.parent().zero()

        for monomial, coefficient in self._coefficients.items():
            assert coefficient

            exponent = monomial.count(gen)

            if not exponent:
                continue

            monomial = list(monomial)
            monomial.remove(gen)
            monomial = tuple(monomial)

            derivative += self.parent()({monomial: exponent * coefficient})

        return derivative

    def map_coefficients(self, f, ring=None):
        if ring is None:
            ring = self.parent()

        return ring({key: image for key, value in self._coefficients.items() if (image := f(value))})

    def __call__(self, values):
        r"""
        Return the value of this symbolic expression.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: a = R.gen(('imag', 0, 0))
            sage: b = R.gen(('real', 0, 0))
            sage: a({a: 1})
            1.00000000000000
            sage: a({a: 2, b: 1})
            2.00000000000000
            sage: (2 * a * b)({a: 3, b: 5})
            30.0000000000000

        """
        def evaluate(monomial):
            product = self.parent().base_ring().one()

            for variable in monomial:
                product *= values[self.parent().gen(variable)]

            return product

        return sum([
            coefficient * evaluate(monomial) for (monomial, coefficient) in self._coefficients.items()
            ])


class SymbolicCoefficientRing(UniqueRepresentation, CommutativeRing):
    @staticmethod
    def __classcall__(cls, surface, base_ring, category=None):
        from sage.categories.all import CommutativeRings
        return super().__classcall__(cls, surface, base_ring, category or CommutativeRings())

    def __init__(self, surface, base_ring, category):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing
            sage: R = SymbolicCoefficientRing(T, CC)
            sage: R.has_coerce_map_from(CC)
            True

            sage: TestSuite(R).run()

        """
        self._surface = surface
        self._base_ring = base_ring

        CommutativeRing.__init__(self, base_ring, category=category, normalize=False)
        self.register_coercion(base_ring)

    Element = SymbolicCoefficientExpression

    def _repr_(self):
        return f"Ring of Power Series Coefficients over {self.base_ring()}"

    def change_ring(self, ring):
        return SymbolicCoefficientRing(self._surface, ring, category=self.category())

    def real_field(self):
        from sage.all import RealField
        return RealField(self._base_ring.prec())

    def is_exact(self):
        return self.base_ring().is_exact()

    def base_ring(self):
        return self._base_ring

    def _coerce_map_from_(self, other):
        if isinstance(other, SymbolicCoefficientRing):
            return self.base_ring().has_coerce_map_from(other.base_ring())

    def _element_constructor_(self, x):
        if isinstance(x, SymbolicCoefficientExpression):
            return x.map_coefficients(self.base_ring(), ring=self)

        if isinstance(x, tuple):
            # x describes a variable
            assert x == tuple(sorted(x))
            return self.element_class(self, {x: self._base_ring.one()})

        if isinstance(x, dict):
            return self.element_class(self, x)

        from sage.all import parent
        if parent(x) is self._base_ring:
            if not x:
                return self.element_class(self, {})
            return self.element_class(self, {(): x})

        raise NotImplementedError(f"symbolic expression from {x}")

    @cached_method
    def imaginary_unit(self):
        from sage.all import I
        return self(self._base_ring(I))

    def gen(self, n):
        if isinstance(n, tuple):
            if len(n) == 3:
                kind, polygon, k = n

                if kind == "real":
                    kind = 0
                elif kind == "imag":
                    kind = 1
                else:
                    raise NotImplementedError

                n = k * 2 * self._surface.num_polygons() + 2 * polygon + kind
            elif len(n) == 2:
                kind, k = n

                if kind != "lagrange":
                    raise ValueError
                if k < 0:
                    raise ValueError

                n = -k-1
            else:
                raise ValueError

        from sage.all import parent, ZZ
        if parent(n) is ZZ:
            n = int(n)

        if not isinstance(n, int):
            raise ValueError

        return self((n,))

    def ngens(self):
        raise NotImplementedError


class PowerSeriesConstraints:
    r"""
    A collection of (linear) constraints on the coefficients of power series
    developed at the vertices of the Voronoi cells of a Delaunay triangulation.

    This is used to create harmonic differentials from cohomology classes.
    """

    def __init__(self, surface, prec, bitprec=None, geometry=None):
        from sage.all import log, ceil, factorial

        self._surface = surface
        self._prec = prec
        self._bitprec = bitprec or ceil(log(factorial(prec), 2) + 53)
        self._geometry = geometry or GeometricPrimitives(surface)
        self._constraints = []
        self._cost = self.symbolic_ring().zero()

    def __repr__(self):
        return repr(self._constraints)

    @cached_method
    def symbolic_ring(self, base_ring=None):
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
            Ring of Power Series Coefficients over Complex Field with 56 bits of precision

        """
        return SymbolicCoefficientRing(self._surface, base_ring=base_ring or self.complex_field())

    @cached_method
    def complex_field(self):
        from sage.all import ComplexField
        return ComplexField(self._bitprec)

    @cached_method
    def real_field(self):
        from sage.all import RealField
        return RealField(self._bitprec)

    @cached_method
    def gen(self, triangle, k, /, conjugate=False):
        assert conjugate is True or conjugate is False
        real = self.real(triangle, k)
        imag = self.imag(triangle, k)

        i = self.symbolic_ring().imaginary_unit()

        if conjugate:
            i = -i

        return real + i*imag

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
            Re(a0,0)
            sage: C.real(0, 1)
            Re(a0,1)
            sage: C.real(1, 2)
            Re(a1,2)

        """
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator for this triangle")

        return self.symbolic_ring().gen(("real", triangle, k))

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
            Im(a0,0)
            sage: C.imag(0, 1)
            Im(a0,1)
            sage: C.imag(1, 2)
            Im(a1,2)

        """
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator for this triangle")

        return self.symbolic_ring().gen(("imag", triangle, k))

    @cached_method
    def lagrange(self, k):
        return self.symbolic_ring().gen(("lagrange", k))

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
            Re(a0,0)
            sage: C.real_part(C.real(0, 0))
            Re(a0,0)
            sage: C.real_part(C.imag(0, 0))
            Im(a0,0)
            sage: C.real_part(2*C.gen(0, 0))  # tol 1e-9
            2*Re(a0,0)
            sage: C.real_part(2*I*C.gen(0, 0))  # tol 1e-9
            -2.0000000000000*Im(a0,0)

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
            Im(a0,0)
            sage: C.imaginary_part(C.real(0, 0))
            0.0000000000000000
            sage: C.imaginary_part(C.imag(0, 0))
            0.0000000000000000
            sage: C.imaginary_part(2*C.gen(0, 0))  # tol 1e-9
            2*Im(a0,0)
            sage: C.imaginary_part(2*I*C.gen(0, 0))  # tol 1e-9
            2*Re(a0,0)

        """
        return self.project(x, "imag")

    def add_constraint(self, expression):
        total_degree = expression.total_degree()

        if total_degree == -1:
            return

        if total_degree == 0:
            raise ValueError(f"cannot solve for constraint {expression} == 0")

        if total_degree > 1:
            raise NotImplementedError("can only encode linear constraints")

        if expression.parent().base_ring() is self.real_field():
            # TODO: Should we scale?
            from sage.all import oo
            # self._constraints.append(expression / expression.norm(1))
            self._constraints.append(expression)
        elif expression.parent().base_ring() is self.complex_field():
            self.add_constraint(expression.real())
            self.add_constraint(expression.imag())
        else:
            raise NotImplementedError("cannot handle expressions over this base ring")

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
            Re(a0,0) + 1.000000000000000*I*Im(a0,0) + (Re(a0,1) + 1.000000000000000*I*Im(a0,1))*z + (Re(a0,2) + 1.000000000000000*I*Im(a0,2))*z^2
            sage: C.develop(1, 1)
            Re(a1,0) + 1.000000000000000*I*Im(a1,0) + Re(a1,1) + 1.000000000000000*I*Im(a1,1) + Re(a1,2) + 1.000000000000000*I*Im(a1,2) + (Re(a1,1) + 1.000000000000000*I*Im(a1,1) + 2.000000000000000*Re(a1,2) + 2.000000000000000*I*Im(a1,2))*z + (Re(a1,2) + 1.000000000000000*I*Im(a1,2))*z^2

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
            sage: C = PowerSeriesConstraints(T, prec=1)

            sage: C.integrate(H())
            0.000000000000000

            sage: a, b = H.gens()
            sage: C.integrate(a)
            (0.500000000000000 + 0.500000000000000*I)*Re(a0,0) + (0.500000000000000 - 0.500000000000000*I)*Im(a0,0) + (0.500000000000000 + 0.500000000000000*I)*Re(a1,0) + (0.500000000000000 - 0.500000000000000*I)*Im(a1,0)
            sage: C.integrate(b)
            (-0.500000000000000)*Re(a0,0) + 0.500000000000000*I*Im(a0,0) + (-0.500000000000000)*Re(a1,0) + 0.500000000000000*I*Im(a1,0)

            sage: C = PowerSeriesConstraints(T, prec=5)
            sage: C.integrate(a) + C.integrate(-a)
            0.00000000000000000
            sage: C.integrate(b) + C.integrate(-b)
            0.00000000000000000

        """
        surface = cycle.surface()

        R = self.symbolic_ring()

        expression = R.zero()

        for path, multiplicity in cycle.voronoi_path().monomial_coefficients().items():

            for S, T in zip((path[-1],) + path, path):
                # Integrate from the midpoint of the edge of S to the midpoint
                # of the edge of T by crossing over the face of T.
                S = surface.opposite_edge(*S)
                assert S[0] == T[0], f"consecutive elements of a path must be attached to the same face in {path} but {S} and {T} do not have that property"

                # Namely we integrate the power series defined around the Voronoi vertex of S by symbolically integrating each monomial term.

                # The midpoints of the edges
                P = self.complex_field()(*self._geometry.midpoint(*S))
                Q = self.complex_field()(*self._geometry.midpoint(*T))

                P_power = P
                Q_power = Q

                for k in range(self._prec):
                    expression += multiplicity * self.gen(S[0], k) / (k + 1) * (Q_power - P_power)

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
            Re(a0,0) + 1.000000000000000*I*Im(a0,0)
            sage: C.evaluate(1, 0)
            Re(a1,0) + 1.000000000000000*I*Im(a1,0)
            sage: C.evaluate(1, 2)
            Re(a1,0) + 1.000000000000000*I*Im(a1,0) + 2.000000000000000*Re(a1,1) + 2.000000000000000*I*Im(a1,1) + 4.000000000000000*Re(a1,2) + 4.000000000000000*I*Im(a1,2)

        """
        # TODO: Check that Δ is within the radius of convergence.

        if derivative >= self._prec:
            raise ValueError

        parent = self.symbolic_ring()

        value = parent.zero()

        z = self.complex_field().one()

        from sage.all import factorial
        factor = self.complex_field()(factorial(derivative))

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

            if abs(Δ0 - Δ1) < 1e-6:
                # Force power series to be identical if they have the same center of Voronoi cell.
                for k in range(self._prec):
                    self.add_constraint(self.gen(triangle0, k) - self.gen(triangle1, k))

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
            [Re(a0,0) - Re(a1,0), Im(a0,0) - Im(a1,0), Re(a0,0) - Re(a1,0), Im(a0,0) - Im(a1,0)]

        If we add more coefficients, we get three pairs of contraints for the
        three edges surrounding a face; for the edge on which the centers of
        the Voronoi cells fall, we get a more restrictive constraint forcing
        all the coefficients to be equal (these are the first four constraints
        recorded below.)::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a0,0) - 0.500000000000000*Im(a0,1) - Re(a1,0) - 0.500000000000000*Im(a1,1),
             Im(a0,0) + 0.500000000000000*Re(a0,1) - Im(a1,0) + 0.500000000000000*Re(a1,1),
             Re(a0,0) - 0.500000000000000*Re(a0,1) - Re(a1,0) - 0.500000000000000*Re(a1,1),
             Im(a0,0) - 0.500000000000000*Im(a0,1) - Im(a1,0) - 0.500000000000000*Im(a1,1)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_midpoint_derivatives(2)
            sage: C  # tol 1e-9
            [Re(a0,0) - 0.500000000000000*Im(a0,1) - Re(a1,0) - 0.500000000000000*Im(a1,1),
             Im(a0,0) + 0.500000000000000*Re(a0,1) - Im(a1,0) + 0.500000000000000*Re(a1,1),
             Re(a0,1) - Re(a1,1),
             Im(a0,1) - Im(a1,1),
             Re(a0,0) - 0.500000000000000*Re(a0,1) - Re(a1,0) - 0.500000000000000*Re(a1,1),
             Im(a0,0) - 0.500000000000000*Im(a0,1) - Im(a1,0) - 0.500000000000000*Im(a1,1),
             Re(a0,1) - Re(a1,1),
             Im(a0,1) - Im(a1,1)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            parent = self.symbolic_ring()

            Δ0 = self.complex_field()(*self._geometry.midpoint(triangle0, edge0))
            Δ1 = self.complex_field()(*self._geometry.midpoint(triangle1, edge1))

            # TODO: Are these good constants?
            if abs(Δ0 - Δ1) < 1e-6:
                continue

            # Require that the 0th, ..., derivatives-1th derivatives are the same at the midpoint of the edge.
            for derivative in range(derivatives):
                self.add_constraint(
                    parent(self.evaluate(triangle0, Δ0, derivative)) - parent(self.evaluate(triangle1, Δ1, derivative)))

    def _L2_consistency_edge(self, triangle0, edge0):
        cost = self.symbolic_ring(self.real_field()).zero()

        triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

        # The midpoint of the edge where the triangles meet with respect to
        # the center of the triangle.
        Δ0 = self.complex_field()(*self._geometry.midpoint(triangle0, edge0))
        Δ1 = self.complex_field()(*self._geometry.midpoint(triangle1, edge1))

        # Develop both power series around that midpoint, i.e., Taylor expand them.
        T0 = self.develop(triangle0, Δ0)
        T1 = self.develop(triangle1, Δ1)

        # Write b_n for the difference of the n-th coefficient of both power series.
        # We want to minimize the sum of |b_n|^2 r^2n where r is half the
        # length of the edge we are on.
        b = (T0 - T1).list()
        edge = self._surface.polygon(triangle0).edges()[edge0]
        r2 = (edge[0]**2 + edge[1]**2) / 4

        r2n = r2
        for n, b_n in enumerate(b):
            # TODO: In the article it says that it should be R^n as a
            # factor but R^{2n+2} is actually more reasonable. See
            # https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/Harmonic.20Differentials/near/308863913
            real = b_n.real()
            imag = b_n.imag()
            cost += (real * real + imag * imag) * r2n

            r2n *= r2

        return cost

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
        R = self.symbolic_ring(self.real_field())

        cost = R.zero()

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            cost += self._L2_consistency_edge(triangle0, edge0)

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
            (0.0000000000000000, 0.0000000000000000)
            sage: C._elementary_line_integrals(0, 1, 0)  # tol 1e-9
            (0 - 0.5*I, 0.5 + 0.0*I)
            sage: C._elementary_line_integrals(0, 0, 1)  # tol 1e-9
            (0.0 + 0.5*I, 0.5 - 0.0*I)
            sage: C._elementary_line_integrals(0, 1, 1)  # tol 1e-9
            (-0.1666666667, -0.1666666667)

        """
        ix = self.complex_field().zero()
        iy = self.complex_field().zero()

        triangle = self._surface.polygon(triangle)
        center = triangle.circumscribing_circle().center()

        for v, e in zip(triangle.vertices(), triangle.edges()):
            Δx, Δy = e
            x0, y0 = -center + v

            def f(x, y):
                from sage.all import I
                return self.complex_field()((x + I*y)**n * (x - I*y)**m)

            def fx(t):
                if abs(Δx) < 1e-6:
                    return self.complex_field().zero()
                return f(x0 + t, y0 + t * Δy/Δx)

            def fy(t):
                if abs(Δy) < 1e-6:
                    return self.complex_field().zero()
                return f(x0 + t * Δx/Δy, y0 + t)

            def integrate(value, t0, t1):
                from sage.all import numerical_integral
                # TODO: Should we do something about the error that is stored in [1]?
                return numerical_integral(value, t0, t1)[0]

            C = self.complex_field()
            ix += C(integrate(lambda t: fx(t).real(), 0, Δx), integrate(lambda t: fx(t).imag(), 0, Δx))
            iy += C(integrate(lambda t: fy(t).real(), 0, Δy), integrate(lambda t: fy(t).imag(), 0, Δy))

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
        C = self.complex_field()
        # Write f(n, m) for z^n\overline{z}^m.
        # Then 1/(2m + 1) [d/dx f(n, m+1) - d/dy -i f(n, m+1)] = f(n, m).

        # So we can use Green's theorem to compute this integral by integrating
        # on the boundary of the triangle:
        # -i/(2m + 1) f(n, m + 1) dx + 1/(2m + 1) f(n, m + 1) dy

        ix, iy = self._elementary_line_integrals(triangle, n, m+1)

        from sage.all import I
        return -I/(C(2)*(m + 1)) * ix + C(1)/(2*(m + 1)) * iy

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
        Return an upper bound for the area 1/π ∫ η \wedge \overline{η}.

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

        The correct area would be 1/π here. However, we are overcounting
        because we sum the single Voronoi cell twice. And also, we approximate
        the square with a circle of radius 1/sqrt(2) for another factor π/2::

            sage: η._evaluate(area)  # tol 1e-9
            1.0

        """
        # Our upper bound is integrating the quantity over the circumcircles of
        # the Delaunay triangles instead, i.e.,
        # we integrate the sum of the a_n \overline{a_m} z^n \overline{z}^m.
        # Integration in polar coordinates, namely
        # ∫ z^n\overline{z}^m = ∫_0^R ∫_0^{2π} ir e^{it(n-m)}
        # shows that only when n = m the value does not vanish and is
        # π/(n+1) |a_n|^2 R^(2n+2).
        area = self.symbolic_ring().zero()

        for triangle in self._surface.label_iterator():
            R2 = self.real_field()(self._surface.polygon(triangle).circumscribing_circle().radius_squared())

            for k in range(self._prec):
                area += (self.real(triangle, k)**2 + self.imag(triangle, k)**2) * R2**(k + 1)

        return area

    def optimize(self, f):
        r"""
        Add constraints that optimize the symbolic expression ``f``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_midpoint_derivatives(1)
            sage: R = C.symbolic_ring()
            sage: f = 3*C.real(0, 0)^2 + 5*C.imag(0, 0)^2 + 7*C.real(1, 0)^2 + 11*C.imag(1, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [Re(a0,0) - Re(a1,0),
             Im(a0,0) - Im(a1,0),
             Re(a0,0) - Re(a1,0),
             Im(a0,0) - Im(a1,0),
             6.00000000000000*Re(a0,0) + λ0 + λ2,
             10.0000000000000*Im(a0,0) + λ1 + λ3,
             14.0000000000000*Re(a1,0) - λ0 - λ2,
             22.0000000000000*Im(a1,0) - λ1 - λ3]

        """
        self._cost += f

    def _optimize_cost(self):
        # We use Lagrange multipliers to rewrite this expression.
        # If we let
        #   L(Re(a), Im(a), λ) = f(Re(a), Im(a)) + Σ λ_i g_i(Re(a), Im(a))
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

                    L = self._cost.derivative(gen)

                    for i in range(lagranges):
                        L += g[i][gen] * self.lagrange(i)

                    self.add_constraint(L)

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

        Integrating along the (negative horizontal) cycle `b`, produces
        something with `-Re(a_0)` in the real part.
        Integration along the diagonal `a`, produces essentially `Re(a_0) -
        Im(a_0)`. Note that the two variables ``a0`` and ``a1`` are the same
        because the centers of the Voronoi cells for the two triangles are
        identical::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.500000000000000*Re(a0,0) + 0.500000000000000*Im(a0,0) + 0.500000000000000*Re(a1,0) + 0.500000000000000*Im(a1,0),
            -0.500000000000000*Re(a0,0) - 0.500000000000000*Re(a1,0) - 1.00000000000000]

        If we increase precision, we see additional higher imaginary parts.
        These depend on the choice of base point of the integration and will be
        found to be zero by other constraints::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.500000000000000*Re(a0,0) + 0.500000000000000*Im(a0,0) - 0.250000000000000*Re(a0,1) + 0.500000000000000*Re(a1,0) + 0.500000000000000*Im(a1,0) + 0.250000000000000*Re(a1,1),
            -0.500000000000000*Re(a0,0) + 0.125000000000000*Re(a0,1) - 0.500000000000000*Re(a1,0) - 0.125000000000000*Re(a1,1) - 1.00000000000000]

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.500000000000000*Re(a0,0) + 0.500000000000000*Im(a0,0) - 0.250000000000000*Re(a0,1) + 0.500000000000000*Re(a1,0) + 0.500000000000000*Im(a1,0) + 0.250000000000000*Re(a1,1),
            -0.500000000000000*Re(a0,0) + 0.125000000000000*Re(a0,1) - 0.500000000000000*Re(a1,0) - 0.125000000000000*Re(a1,1) - 1.00000000000000]

        """
        for cycle in cocycle.parent().homology().gens():
            self.add_constraint(self.real_part(self.integrate(cycle)) - self.real_part(cocycle(cycle)))

    def lagrange_variables(self):
        return set(variable for constraint in self._constraints for variable in constraint.variables("lagrange"))

    def matrix(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_midpoint_derivatives(1)
            sage: C.matrix()
            (
            [ 1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000 -1.00000000000000]
            [ 1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000 -1.00000000000000], (0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000)
            )

        ::

            sage: R = C.symbolic_ring()
            sage: f = 3*C.real(0, 0)^2 + 5*C.imag(0, 0)^2 + 7*C.real(1, 0)^2 + 11*C.imag(1, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C.matrix()
            (
            [ 1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000]
            [ 1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000]
            [ 6.00000000000000 0.000000000000000 0.000000000000000 0.000000000000000  1.00000000000000 0.000000000000000  1.00000000000000 0.000000000000000]
            [0.000000000000000  10.0000000000000 0.000000000000000 0.000000000000000 0.000000000000000  1.00000000000000 0.000000000000000  1.00000000000000]
            [0.000000000000000 0.000000000000000  14.0000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 -1.00000000000000 0.000000000000000]
            [0.000000000000000 0.000000000000000 0.000000000000000  22.0000000000000 0.000000000000000 -1.00000000000000 0.000000000000000 -1.00000000000000], (0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000)
            )

        """
        triangles = list(self._surface.label_iterator())

        prec = int(self._prec)

        from sage.all import matrix, vector
        A = matrix(self.real_field(), len(self._constraints), 2*len(triangles)*prec + len(self.lagrange_variables()))
        b = vector(self.real_field(), len(self._constraints))

        for row, constraint in enumerate(self._constraints):
            for monomial, coefficient in constraint._coefficients.items():
                if not monomial:
                    b[row] = -coefficient
                    continue

                assert len(monomial) == 1

                gen = monomial[0]
                if gen < 0:
                    column = 2*len(triangles)*prec + (-gen-1)
                else:
                    column = gen

                A[row, column] = coefficient

        return A, b

    @cached_method
    def power_series_ring(self, triangle):
        r"""
        Return the power series ring to write down the series describing a
        harmonic differential in a Voronoi cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: Ω = PowerSeriesConstraints(T, 8)
            sage: Ω.power_series_ring(1)
            Power Series Ring in z1 over Complex Field with 69 bits of precision
            sage: Ω.power_series_ring(2)
            Power Series Ring in z2 over Complex Field with 69 bits of precision

        """
        from sage.all import PowerSeriesRing
        return PowerSeriesRing(self.complex_field(), f"z{triangle}")

    def solve(self):
        r"""
        Return a solution for the system of constraints with minimal error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.add_constraint(C.real(0, 0) - C.real(1, 0))
            sage: C.add_constraint(C.real(0, 0) - 1)
            sage: C.solve()
            ({0: 1.00000000000000 + O(z0), 1: 1.00000000000000 + O(z1)}, 0.000000000000000)

        """
        self._optimize_cost()

        A, b = self.matrix()

        from sage.all import ComplexBallField
        C = ComplexBallField(self.complex_field().prec())
        CA = A.change_ring(C)
        Cb = b.change_ring(C)

        solution = CA.solve_right(Cb)

        solution = solution.change_ring(self.real_field())

        residue = (A*solution -b).norm()

        lagranges = len(self.lagrange_variables())

        if lagranges:
            solution = solution[:-lagranges]

        P = self._surface.num_polygons()
        real_solution = [solution[2*p::2*P] for p in range(P)]
        imag_solution = [solution[2*p + 1::2*P] for p in range(P)]

        return {
            triangle: self.power_series_ring(triangle)([self.complex_field()(real_solution[p][k], imag_solution[p][k]) for k in range(self._prec)], self._prec)
            for (p, triangle) in enumerate(self._surface.label_iterator())
        }, residue
