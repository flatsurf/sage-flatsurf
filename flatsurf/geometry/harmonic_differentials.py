from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.all import ZZ


class HarmonicDifferential(Element):
    def __init__(self, parent, coefficients):
        super().__init__(parent)
        self._coefficients = coefficients

    def _repr_(self):
        return repr([coefficients[0] for coefficients in self._coefficients])


class HarmonicDifferentials(UniqueRepresentation, Parent):
    Element = HarmonicDifferential

    @staticmethod
    def __classcall__(cls, surface, coefficients=None, category=None):
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
        if x.parent() is cohomology:
            return self._element_from_cohomology(x, *args, **kwargs)

        raise NotImplementedError()

    def _element_from_cohomology(self, Φ, /, prec):
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
        for cycle, value in zip(Φ.parent().homology().basis(), Φ.vector()):
            constraint = {
                triangle: ([ZZ(0)] * prec, [ZZ(0)] * prec) for triangle in range(self._surface.num_polygons())
            }
            for c in range(len(cycle)):
                triangle_, edge_ = self._surface.opposite_edge(*cycle[(c - 1) % len(cycle)])
                triangle, edge = cycle[c]

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
        print(residues)
        solution = solution[:-lagranges]
        solution = [solution[k*prec:(k+1)*prec] for k in range(self._surface.num_polygons())]
        return self.element_class(self, solution)
