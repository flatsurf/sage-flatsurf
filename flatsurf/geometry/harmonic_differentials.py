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
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.all import ZZ
from dataclasses import dataclass
from collections import namedtuple


class HarmonicDifferential(Element):
    def __init__(self, parent, series):
        super().__init__(parent)
        self._series = series

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

    @staticmethod
    def _integrate_symbolic(cycle, prec):
        r"""
        Return the linear combination of the power series coefficients that
        decsribe the integral of a differential along the homology class
        ``cycle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: from flatsurf.geometry.harmonic_differentials import HarmonicDifferential
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: HarmonicDifferential._integrate_symbolic(H(), prec=5)
            {}

            sage: a, b = H.gens()
            sage: HarmonicDifferential._integrate_symbolic(a, prec=5)
            {0: [(0.5+0.5j),
              (-0.25+0j),
              (0.041666666666666664-0.041666666666666664j),
              0j,
              (0.00625+0.00625j)],
             1: [(0.5+0.5j),
              (0.25+0j),
              (0.041666666666666664-0.041666666666666664j),
              0j,
              (0.00625+0.00625j)]}
            sage: HarmonicDifferential._integrate_symbolic(b, prec=5)
            {0: [(-0.5+0j),
              (0.125+0j),
              (-0.041666666666666664+0j),
              (0.015625+0j),
              (-0.00625+0j)],
             1: [(-0.5+0j),
              (-0.125+0j),
              (-0.041666666666666664+0j),
              (-0.015625+0j),
              (-0.00625+0j)]}

        """
        surface = cycle.surface()

        expression = {}

        for path, multiplicity in cycle.voronoi_path().monomial_coefficients().items():

            for S, T in zip((path[-1],) + path, path):
                # Integrate from the midpoint of the edge of S to the midpoint of the edge of T
                S = surface.opposite_edge(*S)
                assert S[0] == T[0], f"consecutive elements of a path must be attached to the same face in {path} but {S} and {T} do not have that property"

                # Namely we integrate the power series defined around the Voronoi vertex of S by symbolically integrating each monomial term.
                cell = S[0]
                coefficients = expression.get(cell, [0] * prec)

                # The midpoints of the edges
                P = HarmonicDifferential._midpoint(surface, *S)
                Q = HarmonicDifferential._midpoint(surface, *T)

                for k in range(prec):
                    coefficients[k] -= multiplicity * P**(k + 1) / (k + 1)
                    coefficients[k] += multiplicity * Q**(k + 1) / (k + 1)

                expression[cell] = coefficients

        return expression

    @staticmethod
    def _evaluate_symbolic(Δ, derivative, prec):
        from sage.all import ZZ, factorial
        return [ZZ(0) if k < derivative else factorial(k) / factorial(k - derivative) * Δ**(k - derivative) for k in range(prec)]

    def _evaluate(self, expression):
        r"""
        Evaluate an expression by plugging in the coefficients of the power
        series defining this differential.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f); η
            (-0.200000000000082 - 0.500000000000258*I + ... at 0, -0.199999999999957 - 0.500000000001160*I + ... at 1)

        Compute the sum of the constant coefficients::

            sage: η._evaluate({0: [1], 1: [1]})
            -0.400000000000039 - 1.00000000000142*I

        """
        return sum(c*a for face in expression for (c, a) in zip(expression.get(face, []), self._series[face]._coefficients))

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
        # TODO: The above outputs are wrong :(
        return self._evaluate(HarmonicDifferential._integrate_symbolic(cycle, max(series.precision() for series in self._series.values())))

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
        (O(z^5) at 0, O(z^5) at 1)

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

    def _repr_(self):
        return f"Ω({self._surface})"

    def _element_constructor_(self, x, *args, **kwargs):
        from flatsurf.geometry.cohomology import SimplicialCohomology
        cohomology = SimplicialCohomology(self._surface, self._coefficients)

        if not x:
            return self._element_from_cohomology(cohomology(), *args, **kwargs)

        if x.parent() is cohomology:
            return self._element_from_cohomology(x, *args, **kwargs)

        raise NotImplementedError()

    def _element_from_cohomology(self, cocycle, /, prec=5, consistency=2):
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
        for (triangle, edge) in self._surface.edge_iterator():
            triangle_, edge_ = self._surface.opposite_edge(triangle, edge)
            for derivative in range(consistency):
                expected = η._evaluate({triangle: HarmonicDifferential._evaluate_symbolic(HarmonicDifferential._midpoint(self._surface, triangle, edge), derivative, prec)})
                other = η._evaluate({triangle_: HarmonicDifferential._evaluate_symbolic(HarmonicDifferential._midpoint(self._surface, triangle_, edge_), derivative, prec)})
                abs_error = abs(expected - other)
                rel_error = abs(abs_error / expected)
                if abs_error > 1e-9 and rel_error > 1e-6:
                    print(f"power series defining harmonic differential are not consistent: {derivative}th derivate does not match between {(triangle, edge)} and {(triangle_, edge_)}; relative error is {reL_error:.6f}")

        # (2) Check that differential actually integrates like the cohomology class.
        for gen in cocycle.parent().homology().gens():
            expected = cocycle(gen)
            actual = η.integrate(gen)
            error = abs((expected - actual) / expected)
            if error > 1e-6:
                print(f"harmonic differential does not have prescribed integral at {gen}; relative error is {error:.6f}")

        # (3) Check that the area is finite.

        return η


class PowerSeries:
    r"""
    A power series developed around the center of a Voronoi cell.

    This class is used in the implementation of :class:`HormonicDifferential`.
    A harmonic differential is a power series developed at each Voronoi cell.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
        sage: T.set_immutable()

        sage: from flatsurf.geometry.harmonic_differentials import PowerSeries
        sage: PowerSeries(T, 0, [1, 2, 3, 0, 0])
        1 + 2*z + 3*z^2 + O(z^5) at 0

    """

    def __init__(self, surface, polygon, coefficients):
        self._surface = surface
        self._polygon = polygon
        self._coefficients = coefficients

    def precision(self):
        return len(self._coefficients)

    def __repr__(self):
        from sage.all import Sequence, PowerSeriesRing, O
        R = Sequence(self._coefficients, immutable=True).universe()
        R = PowerSeriesRing(R, 'z')
        f = R(self._coefficients) + O(R.gen()**len(self._coefficients))
        return f"{f} at {self._polygon}"


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

        Coefficient = namedtuple("Coefficient", ["real", "imag"])

        def get(self, triangle, k):
            r"""
            Return the coefficients that are multiplied with the coefficient
            a_k of the power series for the ``triangle`` in this linear
            constraint.
            """
            from sage.all import ZZ

            real = self.real.get(triangle, [])[k:k+1]
            imag = self.imag.get(triangle, [])[k:k+1]

            return PowerSeriesConstraints.Constraint.Coefficient(
                real=real[0] if real else ZZ(0),
                imag=imag[0] if imag else ZZ(0))

    def __init__(self, surface, prec):
        self._surface = surface
        self._prec = prec
        self._constraints = []

    def __repr__(self):
        return repr(self._constraints)

    def add_constraint(self, coefficients, value, real=True, imag=True):
        def _imag(x):
            x = x.imag
            if callable(x):
                x = x()
            return x

        def _real(x):
            if hasattr(x, "real"):
                x = x.real
            if callable(x):
                x = x()
            return x

        # We encode a constraint Σ c_i a_i = v as its real and imaginary part.
        # (Our solver can handle complex systems but we also want to add
        # constraints that only concern the real part of the a_i.)
        # We have Re(Σ c_i a_i) = Σ Re(c_i) Re(a_i) - Im(c_i) Im(a_i)
        #     and Im(Σ c_i a_i) = Σ Im(c_i) Re(a_i) + Re(c_i) Im(a_i).
        if real:
            self._add_constraint(
                real={triangle: [_real(c) for c in coefficients[triangle]] for triangle in coefficients.keys()},
                imag={triangle: [-_imag(c) for c in coefficients[triangle]] for triangle in coefficients.keys()},
                value=_real(value))
        if complex:
            self._add_constraint(
                real={triangle: [_imag(c) for c in coefficients[triangle]] for triangle in coefficients.keys()},
                imag={triangle: [_real(c) for c in coefficients[triangle]] for triangle in coefficients.keys()},
                value=_imag(value))

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

    def require_consistency(self, prec):
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
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0)]

        If we add more coefficients, we get three pairs of contraints for the
        three edges surrounding a face::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_consistency(1)
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [-0.0, -0.5], 1: [0.0, -0.5]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 0.5], 1: [-0.0, 0.5]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0, -0.5], 1: [-1.0, -0.5]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.5], 1: [-1.0, -0.5]}, lagrange=[], value=0)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_consistency(2)
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0, 1.0], 1: [0, -1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [0, 1.0], 1: [0, -1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={0: [-0.0, -0.5], 1: [0.0, -0.5]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0.0, 0.5], 1: [-0.0, 0.5]}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [1.0, -0.5], 1: [-1.0, -0.5]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0, -0.5], 1: [-1.0, -0.5]}, lagrange=[], value=0)]

        """
        if prec > self._prec:
            raise ValueError("prec must not exceed global precision")

        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            Δ0 = HarmonicDifferential._midpoint(self._surface, triangle0, edge0)
            Δ1 = HarmonicDifferential._midpoint(self._surface, triangle1, edge1)

            # Require that the 0th, ..., prec-1th derivatives are the same at the midpoint of the edge.
            # The series f(z) = Σ a_k z^k has derivative Σ k!/(k-d)! a_k z^{k-d}
            for d in range(prec):
                self.add_constraint({
                    triangle0: HarmonicDifferential._evaluate_symbolic(Δ0, d, self._prec),
                    triangle1: [-c for c in HarmonicDifferential._evaluate_symbolic(Δ1, d, self._prec)],
                }, ZZ(0))

    def require_finite_area(self):
        r"""
        Since the area 1 /(2iπ) ∫ η \wedge \overline{η} must be finite [TODO:
        REFERENCE?] we can optimize for this quantity to be minimal.

        EXPMALES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
            sage: T.set_immutable()

        EXAMPLES::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 1)
            sage: C.require_consistency(1)

            sage: C.require_finite_area()
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [1.0], 1: [-1.0]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [1.0], 1: [-1.0]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [2.0]}, imag={}, lagrange=[-1.0], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [2.0]}, lagrange=[0, -1.0], value=0),
             PowerSeriesConstraints.Constraint(real={1: [2.0]}, imag={}, lagrange=[1.0], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={1: [2.0]}, lagrange=[0, 1.0], value=0)]

        """
        # To make our lives easier, we do not optimize this value but instead
        # the sum of the |a_k|^2·radius^k = (Re(a_k)^2 + Im(a_k)^2)·radius^k
        # which are easier to compute. (TODO: Explain why this is reasonable to
        # do instead.) We rewrite this condition using Langrange multipliers
        # into a system of linear equations.

        # If we let
        #   L(Re(a), Im(a), λ) = f(Re(a), Im(a)) - Σ λ_i g_i(Re(a), Im(a))
        # and denote by f(Re(a), Im(a)) the above sum, and by g_i=0 all the
        # affine linear conditions collected so far, then we get two
        # constraints for each a_k, one real, one imaginary, namely that the
        # partial derivative wrt Re(a_k) and Im(a_k) vanishes. Note that we
        # have one Lagrange multiplier for each affine linear constraint
        # collected so far.
        lagranges = len(self._constraints)

        # if any(constraint.value for constraint in self._constraints):
        #     raise Exception("cannot add Lagrange Multiplier constraint when non-linear constraints have been added")

        g = self._constraints

        # We form the partial derivative with respect to the variables Re(a_k)
        # and Im(a_k):
        for triangle in range(self._surface.num_polygons()):
            R = float(self._surface.polygon(triangle).circumscribing_circle().radius_squared().sqrt())

            for k in range(self._prec):
                # We get a constraint by forming dL/dRe(a_k) = 0, namely
                # 2 Re(a_k) R^k - Σ λ_i dg_i/dRe(a_k) = 0
                self._add_constraint(
                    real={triangle: [ZZ(0) if j != k else 2*R**k for j in range(k+1)]},
                    imag={},
                    lagrange=[-g[i].get(triangle, k).real for i in range(lagranges)],
                    value=ZZ(0))
                # We get another constraint by forming dL/dIm(a_k) = 0, namely
                # 2 Im(a_k) R^k + λ^t dg/dIm(a_k) = 0
                self._add_constraint(
                    real={},
                    imag={triangle: [ZZ(0) if j != k else 2*R**k for j in range(k+1)]},
                    lagrange=[-g[i].get(triangle, k).imag for i in range(lagranges)],
                    value=ZZ(0))

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
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [0.5, -0.25], 1: [0.5, 0.25]}, imag={0: [-0.5], 1: [-0.5]}, lagrange=[], value=1),
             PowerSeriesConstraints.Constraint(real={0: [0.5], 1: [0.5]}, imag={0: [0.5, -0.25], 1: [0.5, 0.25]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, imag={}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, lagrange=[], value=0)]

        ::

            sage: C = PowerSeriesConstraints(T, 2)
            sage: C.require_cohomology(H({b: 1}))
            sage: C
            [PowerSeriesConstraints.Constraint(real={0: [0.5, -0.25], 1: [0.5, 0.25]}, imag={0: [-0.5], 1: [-0.5]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [0.5], 1: [0.5]}, imag={0: [0.5, -0.25], 1: [0.5, 0.25]}, lagrange=[], value=0),
             PowerSeriesConstraints.Constraint(real={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, imag={}, lagrange=[], value=1),
             PowerSeriesConstraints.Constraint(real={}, imag={0: [-0.5, 0.125], 1: [-0.5, -0.125]}, lagrange=[], value=0)]

        """
        for cycle in cocycle.parent().homology().gens():
            coefficients = HarmonicDifferential._integrate_symbolic(cycle, self._prec)
            self.add_constraint(coefficients, cocycle(cycle), imag=False)

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
            {0: 1.00000000000000 + O(z) at 0, 1: 1.00000000000000 + O(z) at 1}

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
            triangle: PowerSeries(self._surface, triangle, [CC(solution[k][p], solution[k][p + self._prec]) for p in range(self._prec)])
            for (k, triangle) in enumerate(self._surface.label_iterator())
        }
