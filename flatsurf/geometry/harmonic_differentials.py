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
    def _integrate_symbolic(chain, prec):
        r"""
        Return the linear combination of the power series coefficients that
        decsribe the integral of a differential along the homology class
        ``chain``.

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
        surface = chain.surface()

        expression = {}

        for path, multiplicity in chain.voronoi_path().monomial_coefficients().items():

            for S, T in zip(path[1:] + path[:1], path):
                # Integrate from the midpoint of the edge of S to the midpoint of the edge of T
                S = surface.opposite_edge(*S)
                assert S[0] == T[0], f"consecutive elements of a path must be attached to the same face in {path}"

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
            (-0.179500990254010 ... at 0, -0.350210481202233 ... at 1)

        Compute the sum of the constant coefficients::

            sage: η._evaluate({0: [1], 1: [1]})
            -0.529711471456243

        """
        return sum(c*a for face in expression for (c, a) in zip(expression.get(face, []), self._series[face]._coefficients))

    def integrate(self, chain):
        r"""
        Return the integral of this differential along the homology class
        ``chain``.

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

            sage: η.integrate(a)
            (-0.2998018977174907-0.268495653337972j)

            sage: η.integrate(b)
            (0.27770664134510037+0j)

        """
        # TODO: The above outputs are wrong :(
        return self._evaluate(HarmonicDifferential._integrate_symbolic(chain, max(series.precision() for series in self._series.values())))

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
        sage: η.integrate(a)  # not tested
        0
        sage: η.integrate(b)  # not tested
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

    def _element_from_cohomology(self, Φ, /, prec=5):
        # We develop a consistent system of Laurent series at each vertex of the Voronoi diagram
        # to describe a differential.

        # Let η be the differential we are looking for. To describe η we will use ω, the differential
        # corresponding to the flat structure given by our triangulation on this Riemann surface.
        # Then f:=η/ω is a meromorphic function which we can develop locally into a Laurent series.
        # Away from the vertices of the triangulation, ω has no zeros, so f has no poles there and is
        # thus given by a power series.

        # At each vertex of the Voronoi diagram, write f=Σ a_k x^k + O(x^prec). Our task is now to determine
        # the a_k.

        # We use several different constraints to determine these coefficients. We encode each complex
        # coefficient as two real variables, its real and imaginary part.
        constraints = ([], [])

        def add_real_constraint(coefficients, value):
            constraint = [0] * self._surface.num_polygons() * prec * 2
            for triangle in coefficients.keys():
                if triangle == "lagrange":
                    constraint.extend(coefficients[triangle])
                    continue

                constraint[triangle * prec * 2:(triangle + 1) * prec * 2:2] = coefficients[triangle][0]
                constraint[triangle * prec * 2 + 1:(triangle + 1) * prec * 2 + 1:2] = coefficients[triangle][1]

            constraints[0].append(constraint)
            constraints[1].append(value)

        def add_complex_constraint(coefficients, value):
            add_real_constraint({
                triangle: (
                    [c.real() for c in coefficients[triangle]],
                    [-c.imag() for c in coefficients[triangle]],
                ) for triangle in coefficients.keys()
            }, value.real())
            add_real_constraint({
                triangle: (
                    [c.imag() for c in coefficients[triangle]],
                    [c.real() for c in coefficients[triangle]],
                ) for triangle in coefficients.keys()
            }, value.imag())

        def Δ(triangle, edge):
            r"""
            Return the vector from the center of the circumcircle to the center of the edge.
            """
            P = self._surface.polygon(triangle)
            return -P.circumscribing_circle().center() + P.vertex(edge) + P.edge(edge) / 2

        def complex(*x):
            C = self.base_ring().algebraic_closure()
            return C(*x)

        # (1) The radius of convergence of the power series is the distance from the vertex of the Voronoi
        # cell to the closest vertex of the triangulation (since we use a Delaunay triangulation, all vertices
        # are at the same distance in fact.) So the radii of convergence of two neigbhouring cells overlap
        # and the power series must coincide there. Note that this constraint is unrelated to the cohomology
        # class Φ.
        for triangle0, edge0 in self._surface.edge_iterator():
            triangle1, edge1 = self._surface.opposite_edge(triangle0, edge0)

            if triangle1 < triangle0:
                # Add each constraint only once.
                continue

            Δ0 = Δ(triangle0, edge0)
            Δ1 = Δ(triangle1, edge1)

            # TODO: Don't go all the way to prec. The higher degrees are probably not that useful numerically.
            from sage.all import binomial
            for j in range(prec // 2):
                add_complex_constraint({
                    triangle0: [ZZ(0) if k < j else binomial(k, j)*complex(*Δ0)**(k - j) for k in range(prec)],
                    triangle1: [ZZ(0) if k < j else binomial(k, j)*complex(*Δ1)**(k - j) for k in range(prec)]
                }, ZZ(0))

        # (2) Since the area 1 /(2iπ) ∫ η \wedge \overline{η} must be finite [TODO: REFERENCE?] we optimize for
        # this quantity to be minimal. To make our lives easier, we do not optimize this value but instead
        # the sum of the |a_k|^2 radius^k = (Re(a_k)^2 + Im(a_k)^2)·radius^k which are easier to compute.
        # (TODO: Explain why this is reasonable to do instead.)
        # We rewrite this condition using Langrange multipliers into a system of linear equations.

        # If we let L(Re(a), Im(a), λ) = f(Re(a), Im(a)) + λ^t g(Re(a), Im(a)) and denote by f(Re(a), Im(a))
        # the above sum, then we get two constraints for each a_k, one real, one imaginary, namely that the
        # partial derivative for that Re(a_k) or Im(a_k) vanishes.
        # Note that we have one Lagrange multiplier for each linear constraint collected so far.
        lagranges = len(constraints[0])
        for triangle in range(self._surface.num_polygons()):
            R = 1  # TODO: Use the correct radius of convergence of this triangle.
            for k in range(prec):
                add_real_constraint({
                    triangle: (
                        [ZZ(0) if j != k else 2*R**k for j in range(prec)],
                        [ZZ(0) for j in range(prec)]),
                    "lagrange": [constraints[0][j][triangle * prec * 2 + 2*k] for j in range(lagranges)],
                }, ZZ(0))
                add_real_constraint({
                    triangle: (
                        [ZZ(0) for j in range(prec)],
                        [ZZ(0) if j != k else 2*R**k for j in range(prec)]),
                    "lagrange": [constraints[0][j][triangle * prec * 2 + 2*k + 1] for j in range(lagranges)],
                }, ZZ(0))

        # (3) We have that for any cycle γ, Re(∫fω) = Re(∫η) = Φ(γ). We can turn this into constraints
        # on the coefficients as we integrate numerically following the path γ as it intersects the radii of
        # convergence.
        for cycle in Φ.parent().homology().gens():
            value = Φ(cycle)
            constraint = {
                triangle: ([ZZ(0)] * prec, [ZZ(0)] * prec) for triangle in range(self._surface.num_polygons())
            }
            for path, multiplicity in cycle.voronoi_path():
                assert multiplicity == 1
                for c in range(len(path)):
                    # TODO: Zip path and a shifted path instead.
                    triangle_, edge_ = self._surface.opposite_edge(*path[c - 1])
                    triangle, edge = path[c]

                    assert self._surface.singularity(triangle_, edge_) == self._surface.singularity(triangle, edge), "cycle is not closed"

                    while (triangle_, edge_) != (triangle, edge):
                        coefficients = constraint[triangle_][0]

                        P = Δ(triangle_, edge_)
                        for k in range(prec):
                            coefficients[k] -= (complex(*P)**(k + 1)/(k + 1)).real()

                        edge_ = (edge_ + 2) % 3

                        P = Δ(triangle_, edge_)
                        for k in range(prec):
                            coefficients[k] += (complex(*P)**(k + 1)/(k + 1)).real()

                        triangle_, edge_ = self._surface.opposite_edge(triangle_, edge_)

            add_real_constraint(constraint, value)

        # Solve the system of linear constraints to get all the coefficients.
        vars = max(len(constraint) for constraint in constraints[0])

        for constraint in constraints[0]:
            constraint.extend([0] * (vars - len(constraint)))

        import scipy.linalg
        import numpy
        solution, residues, _, _ = scipy.linalg.lstsq(numpy.matrix(constraints[0]), numpy.array(constraints[1]))
        solution = solution[:-lagranges]
        solution = [solution[k*prec:(k+1)*prec] for k in range(self._surface.num_polygons())]
        return self.element_class(self, {
            polygon: PowerSeries(self._surface, polygon, [self._coefficients(c) for c in solution[k]])
            for (k, polygon) in enumerate(self._surface.label_iterator())
            })


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
