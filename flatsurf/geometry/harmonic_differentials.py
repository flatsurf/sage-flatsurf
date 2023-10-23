r"""
TODO: Document this module.
TODO: Rename this module?

EXAMPLES:

We compute harmonic differentials on the square torus::

    sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    sage: T = translation_surfaces.torus((1, 0), (0, 1))
    sage: T.set_immutable()

    sage: H = SimplicialHomology(T)
    sage: a, b = H.gens()

First, the harmonic differentials that sends the horizontal `a` to 1 and the
vertical `b` to zero::

    sage: H = SimplicialCohomology(T)
    sage: f = H({a: 1})
    sage: Ω = HarmonicDifferentials(T)
    sage: ω = Ω(f)  # random output due to deprecation warnings
    sage: ω
    (1.00000000000000 + O(z0^5), 1.00000000000000 + O(z1^5), 1.00000000000000 + O(z2^5), 1.00000000000000 + O(z3^5), 1.00000000000000 + O(z4^5), 1.00000000000000 + O(z5^5), 1.00000000000000 + O(z6^5), 1.00000000000000 + O(z7^5), 1.00000000000000 + O(z8^5))

The harmonic differential that integrates as 0 along `a` but 1 along `b`::

    sage: g = H({b: 1})
    sage: Ω(g)  # tol 1e-9
    (1.00000000000000*I + O(z0^5), 1.00000000000000*I + O(z1^5), 1.00000000000000*I + O(z2^5), 1.00000000000000*I + O(z3^5), 1.00000000000000*I + O(z4^5), 1.00000000000000*I + O(z5^5), 1.00000000000000*I + O(z6^5), 1.00000000000000*I + O(z7^5), 1.00000000000000*I + O(z8^5))

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2023 Julian Rüth
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
from sage.misc.cachefunc import cached_method, cached_function
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement


def integral2(part, α, κ, d, ζd, n, β, λ, dd, ζdd, m, a, b, C, R):
    # Since γ(t) = (1 - t)a + tb, we have |·γ(t)| = |b - a|
    constant = ζd**(κ * (n+1)) / (d+1) * (ζdd**(λ * (m+1)) / (dd+1)).conjugate() * abs(b - a)

    def value(t):
        z = (1 - t) * a + t * b

        if d == 0:
            za = (z-α) ** n
        else:
            za = (z-α).nth_root(d + 1)**(n - d)

        if dd == 0:
            zb = ((z-β) ** m).conjugate()
        else:
            zb = ((z-β).nth_root(dd + 1)**(m - dd)).conjugate()

        value = constant * za * zb

        if part == "Re":
            return float(value.real())
        if part == "Im":
            return float(value.imag())

        raise NotImplementedError

    from scipy.integrate import quad
    integral, error = quad(value, 0, 1)

    return R(integral)


@cached_function
def Cab(C, a):
    return C(*a)


def define_solve():
    import cppyy

    if not hasattr(cppyy.gbl, "solve"):
        cppyy.cppdef(r'''
        #include <cassert>
        #include <vector>
        #include <iostream>
        #include "/home/jule/proj/eskin/sage-flatsurf/mpreal-support.h"
        #include <eigen3/Eigen/LU>
        #include <eigen3/Eigen/Sparse>
        #include <eigen3/Eigen/OrderingMethods>

        using namespace mpfr;
        using namespace Eigen;
        using std::vector;
        typedef SparseMatrix<mpreal> MatrixXmp;
        typedef Matrix<mpreal,Dynamic,1> VectorXmp;

        VectorXmp solve(vector<vector<double>> _A, vector<double> _b)
        {
          // set precision to ? bits (double has only 53 bits)
          mpreal::set_default_prec(256);
          // Declare matrix and vector types with multi-precision scalar type

          const int ROWS = _A.size();
          const int COLS = _A[0].size();
          MatrixXmp A = MatrixXmp(ROWS, COLS);
          VectorXmp b = VectorXmp(ROWS);

          for (int y = 0; y < ROWS; y++) {
            for (int x = 0; x < COLS; x++) {
              if (_A[y][x] != 0) {
                A.insert(y, x) = _A[y][x];
              }
              b[y] = _b[y];
            }
          }

          A.makeCompressed();

          SparseQR<MatrixXmp, NaturalOrdering<int>> QRd;
          QRd.compute(A);

          assert(QRd.info() == Success);

          VectorXmp x = QRd.solve(b);

          assert(QRd.info() == Success);

          return x;
        }
        ''')

    return cppyy.gbl.solve


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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)

            sage: Ω(f) + Ω(f) + Ω(H({a: -2}))
            (O(z0^5), O(z1^5), O(z2^5), O(z3^5), O(z4^5), O(z5^5), O(z6^5), O(z7^5), O(z8^5))

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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
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

        # if kind is None or "midpoint_derivatives" in kind:
        #     C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        #     for derivative in range(self.precision()//3):
        #         for ((label, edge), a, b) in self.parent()._geometry._homology_generators:
        #             opposite_label, opposite_edge = self.parent().surface().opposite_edge(label, edge)
        #             expected = self.evaluate(label, edge, a, C.complex_field()(*self.parent()._geometry.midpoint_on_path_between_centers(label, edge, a, edge, b)), derivative)
        #             other = self.evaluate(opposite_label, opposite_edge, 1 - b, C.complex_field()(*self.parent()._geometry.midpoint_on_path_between_centers(opposite_label, opposite_edge, 1 - b, opposite_edge, 1 - a)), derivative)

        #             abs_error, rel_error = errors(expected, other)

        #             if abs_error > abs_tol or rel_error > rel_tol:
        #                 report = f"Power series defining harmonic differential are not consistent where triangles meet. {derivative}th derivative does not match between {(label, edge, a)} where it is {expected} and {(label, edge, b)} where it is {other}, i.e., there is an absolute error of {abs_error} and a relative error of {rel_error}."
        #                 if verbose:
        #                     print(report)

        #                 error = report
        #                 if not verbose:
        #                     return error

        # if kind is None or "area" in kind:
        #     if verbose:
        #         C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        #         area = self._evaluate(C._area_upper_bound())

        #         report = f"Area (upper bound) is {area}."
        #         print(report)

        if kind is None or "L2" in kind:
            C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
            consistency = C._L2_consistency()
            # print(consistency)
            abs_error = self._evaluate(consistency)

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

        TODO: triangle is not really a triangle anymore.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: η.series(T(0, (1/2, 1/2)))  # abstol 1e-9
            (1.00000000000000 - 6.13512000000000e-17*I) + (6.58662000000000e-33 + 2.58104000000000e-32*I)*z2 + (-2.76745000000000e-31 - 7.82106000000000e-16*I)*z2^2 + (-2.29148000000000e-31 + 1.90965000000000e-31*I)*z2^3 + (1.42943000000000e-15 - 6.22435000000000e-32*I)*z2^4 + O(z2^5)

        """
        return self._series[triangle]

    @cached_method
    def _constraints(self):
        # TODO: This is a hack. Come up with a better class hierarchy!
        return PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)

    def evaluate(self, label, edge, pos, Δ, derivative=0):
        # TODO: This should probably be removed.
        C = self._constraints()
        return self._evaluate(C.evaluate(label, edge, pos, Δ, derivative=derivative))

    def __call__(self, point):
        # TODO: Don't call this __call__. Make it clear in the name that this does not mean anything.
        r"""
        Return the value of this differential at ``point`` on the chart that is
        best suited for evaluating it.

        .. NOTE::

            This operation has no real mathematical meaning. A differential is
            not a function on the surface so it cannot really be evaluated.
            This method exists to compare two differentials, e.g., to plot the
            error between an exact differentials and its numerical
            approximation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True)
            sage: center = S(0, S.polygon(0).circumscribing_circle().center())
            sage: singularity = next(iter(S.vertices()))

            sage: R.<z> = CC[[]]
            sage: ω = Ω({center: R(1, 1), singularity: R(2, 1)})

            sage: ω(center)
            1.00000000000000
            sage: ω(singularity)  # TODO: Make this work at a singularity as well.
            Traceback (most recent call last):
            ...
            ValueError: ...

            sage: ω = Ω({center: R(0, 1), singularity: R(1, 3)})
            sage: point = S(0, S.polygon(0).vertices()[0] + vector((1/1000000, 0)))
            sage: ω(point)  # TODO: Does this make any sense?
            3333.33333333333

        """
        C = self._constraints()
        expression = C.value(point)
        return self._evaluate(expression)

    # TODO: Make this work again.
    # def roots(self):
    #     r"""
    #     Return the roots of this harmonic differential.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()

    #         sage: H = SimplicialHomology(T)
    #         sage: a, b = H.gens()
    #         sage: H = SimplicialCohomology(T)
    #         sage: f = H({a: 1})

    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: η = Ω(f)
    #         sage: η.roots()
    #         []

    #     ::

    #         sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    #         sage: T = translation_surfaces.regular_octagon()
    #         sage: T.set_immutable()

    #         sage: H = SimplicialHomology(T)
    #         sage: a = H.gens()[0]
    #         sage: H = SimplicialCohomology(T)
    #         sage: f = H({a: 1})

    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: η = Ω(f, prec=16, check=False)  # not tested TODO
    #         sage: η.roots()  # not tested TODO

    #     """
    #     roots = []

    #     surface = self.parent().surface()

    #     for triangle in surface.labels():
    #         series = self._series[triangle]
    #         series = series.truncate(series.prec())
    #         for (root, multiplicity) in series.roots():
    #             if multiplicity != 1:
    #                 raise NotImplementedError

    #             root += root.parent()(*self.parent()._geometry.point_on_path_between_centers(triangle))

    #             from sage.all import vector
    #             root = vector(root)

    #             # TODO: Keep roots that are within the circumcircle for sanity checking.
    #             # TODO: Make sure that roots on the edges are picked up by exactly one triangle.

    #             if not surface.polygon(triangle).contains_point(root):
    #                 continue

    #             roots.append((triangle, vector(root)))

    #     # TODO: Deduplicate roots.
    #     # TODO: Compute roots at the vertices.

    #     return roots

    def _evaluate(self, expression):
        r"""
        Evaluate an expression by plugging in the coefficients of the power
        series defining this differential.

        This might not correspond to evaluating the actual power series somewhere.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

        Compute the constant coefficients::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(T, 5, Ω._geometry)
            sage: R = C.symbolic_ring()
            sage: η._evaluate(R(C._gen_nonsingular(0, 0, 0, 0))) # abstol 1e-9
            1.00000000000000 + 1.87854000000000e-72*I

        """
        coefficients = {}

        for variable in expression.variables():
            kind, point, k = variable.describe()

            try:
                coefficient = self._series[point][k]
            except IndexError:
                import warnings
                warnings.warn(f"expected a {k}th coefficient of the power series around {point} but none found")
                coefficients[variable] = 0
                continue

            if kind == "Re":
                coefficients[variable] = coefficient.real()
            else:
                assert kind == "Im"
                coefficients[variable] = coefficient.imag()

        value = expression(coefficients)
        if isinstance(value, SymbolicCoefficientExpression):
            assert value.total_degree() <= 0
            value = value.constant_coefficient()

        return value

    @cached_method
    def precision(self):
        # TODO: This is the number of coefficients of the power series but we use it as bit precision?
        # TODO: There should not be a single global precision. Instead each
        # cell has its own precision. Also these precisions are not comparable,
        # since depending on the degree of the root, we need to scale things.
        precisions = set(series.precision_absolute() for series in self._series.values())
        # assert len(precisions) == 1
        return min(precisions)

    # def cauchy_residue(self, vertex, n, angle=None):
    #     r"""
    #     Return the n-th coefficient of the power series expansion around the
    #     ``vertex`` in local coordinates of that vertex.

    #     Let ``d`` be the degree of the vertex and pick a (d+1)-th root z' of z such
    #     that this differential can be written as a power series Σb_n z'^n. This
    #     method returns the `b_n` by evaluating the Cauchy residue formula,
    #     i.e., 1/2πi ∫ f/z'^{n+1} integrating along a loop around the vertex.

    #     EXAMPLES:

    #     For the square torus, the vertices are no singularities, so harmonic
    #     differentials have no pole at the vertex::

    #         sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()

    #         sage: H = SimplicialHomology(T)
    #         sage: a, b = H.gens()
    #         sage: H = SimplicialCohomology(T)
    #         sage: f = H({a: 1})

    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: η = Ω(f)

    #         sage: vertex = [(0, 0), (0, 1), (0, 2), (0, 3)]
    #         sage: angle = 1

    #         sage: η.cauchy_residue(vertex, 0, 1)  # tol 1e-9
    #         1.0 - 0.0*I
    #         sage: abs(η.cauchy_residue(vertex, -1, 1)) < 1e-9
    #         True
    #         sage: abs(η.cauchy_residue(vertex, -2, 1)) < 1e-9
    #         True
    #         sage: abs(η.cauchy_residue(vertex, 1, 1)) < 1e-9
    #         True
    #         sage: abs(η.cauchy_residue(vertex, 2, 1)) < 1e-9
    #         True

    #     """
    #     surface = self.parent().surface()

    #     if angle is None:
    #         vertex = set(vertex)
    #         for angle, v in surface.angles(return_adjacent_edges=True):
    #             if vertex & set(v):
    #                 assert vertex == set(v)
    #                 vertex = v
    #                 break
    #         else:
    #             raise ValueError("not a vertex in this surface")

    #     parts = []

    #     complex_field = self._series[0].parent().base_ring()

    #     # Integrate real & complex part independently.
    #     for part in ["Re", "Im"]:
    #         integral = 0

    #         # TODO print(f"Building {part} part")

    #         # We integrate along a path homotopic to a (counterclockwise) unit
    #         # circle centered at the vertex.

    #         # We keep track of our last choice of a (d+1)-st root and pick the next
    #         # root that is closest while walking counterclockwise on that circle.
    #         arg = 0

    #         def root(z):
    #             if angle == 1:
    #                 return complex_field(z)

    #             roots = complex_field(z).nth_root(angle, all=True)

    #             positives = [root for root in roots if (root.arg() - arg) % (2*3.14159265358979) > -1e-6]

    #             return min(positives)

    #         for label, edge in vertex:
    #             # We integrate on the line segment from the midpoint of edge to
    #             # the midpoint of the previous edge in triangle, i.e., the next
    #             # edge in counterclockwise order walking around the vertex.
    #             edge_ = (edge - 1) % len(surface.polygon(label).vertices())

    #             # TODO print(f"integrating across {triangle} from the midpoint of {edge} to {edge_}")

    #             P = complex_field(*self.parent()._geometry.midpoint_on_path_between_centers(label, edge))
    #             Q = complex_field(*self.parent()._geometry.midpoint_on_path_between_centers(label, edge_))

    #             # The vector from z=0, the center of the Voronoi cell, to the vertex.
    #             δ = P - complex_field(*surface.polygon(label).edge(edge)) / 2

    #             def at(t):
    #                 z = (1 - t) * P + t * Q
    #                 root_at_z = root(z - δ)
    #                 return z, root_at_z

    #             root_at_P = at(0)[1]
    #             arg = root_at_P.arg()

    #             # TODO print(f"in terms of the Voronoi cell midpoint, this is integrating from {P} to {Q} which lift to {at(0)[1]} and {at(1)[1]} relative to the vertex which is at {δ} from the center of the Voronoi cell")

    #             def integrand(t):
    #                 t = complex_field(t)
    #                 z, root_at_z = at(t)
    #                 denominator = root_at_z ** (n + 1)
    #                 numerator = self.evaluate(label, z)
    #                 integrand = numerator / denominator * (Q - P)
    #                 # TODO print(f"evaluating at {z} resp its root {root_at_z} produced ({numerator}) / ({denominator}) * ({Q - P}) = {integrand}")
    #                 integrand = getattr(integrand, part)()
    #                 return integrand

    #             from sage.all import numerical_integral
    #             integral_on_segment, error = numerical_integral(integrand, 0, 1)
    #             # TODO: What should be do about the error?
    #             integral += integral_on_segment

    #             root_at_Q = at(1)[1]
    #             arg = root_at_Q.arg()

    #         parts.append(integral)

    #     return complex_field(*parts) / complex_field(0, 2*3.14159265358979)

    def integrate(self, cycle, numerical=False):
        # TODO: Generalize to more than just cycles.
        r"""
        Return the integral of this differential along the homology class
        ``cycle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
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
        if numerical:
            raise NotImplementedError

        C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        return self._evaluate(C.integrate(cycle))

    def _repr_(self):
        # TODO: Tune this so we do not loose important information.
        # TODO: Do not print that many digits.
        # TODO: Make it visually clear that we are not printing everything.
        def compress(series):
            def compress_coefficient(coefficient):
                if coefficient.imag().abs() < 1e-9:
                    coefficient = coefficient.parent()(coefficient.real())
                if coefficient.real().abs() < 1e-9:
                    coefficient = coefficient.parent()(coefficient.parent().gen() * coefficient.imag())

                return coefficient

            return series.parent()({exponent: compress_coefficient(coefficient) for (coefficient, exponent) in zip(series.coefficients(), series.exponents())} or 0).add_bigoh(series.precision_absolute())

        return repr(tuple(compress(series) for series in self._series.values()))

    def plot(self, versus=None):
        from sage.all import RealField, I, vector, complex_plot, oo
        S = self.parent().surface()
        GS = S.graphical_surface()

        plot = S.plot(fill=None)

        for label in S.labels():
            P = S.polygon(label)
            PS = GS.graphical_polygon(label)

            if versus:
                if versus.parent().surface().polygon(label).vertices() != P.vertices():
                    raise ValueError

            from sage.all import RR
            PR = P.change_ring(RR)

            def f(z):
                xy = (z.real(), z.imag())
                xy = PS.transform_back(xy)
                if not PR.contains_point(vector(xy)):
                    return oo
                xy = xy - P.circumscribing_circle().center()

                xy = RealField(54)(xy[0]) + I*RealField(54)(xy[1])
                value = self.evaluate(label, edge=None, pos=None, Δ=xy)
                if versus:
                    v = versus.evaluate(label, edge=None, pos=None, Δ=xy)
                    value = (value - v) / v.abs()
                return value

            bbox = PS.bounding_box()
            plot += complex_plot(f, (bbox[0], bbox[2]), (bbox[1], bbox[3]))

        # plot += S.plot(fill=None)

        return plot


class HarmonicDifferentials(UniqueRepresentation, Parent):
    r"""
    The space of harmonic differentials on this surface.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialCohomology
        sage: T = translation_surfaces.torus((1, 0), (0, 1))
        sage: T.set_immutable()

        sage: Ω = HarmonicDifferentials(T); Ω
        Ω(Translation Surface in H_1(0) built from a square)

    ::

        sage: H = SimplicialCohomology(T)
        sage: Ω(H())
        (O(z0^5), O(z1^5), O(z2^5), O(z3^5), O(z4^5), O(z5^5), O(z6^5), O(z7^5), O(z8^5))

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
    def __classcall__(cls, surface, safety=None, singularities=False, category=None):
        r"""
        Normalize parameters when creating the space of harmonic differentials.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: HarmonicDifferentials(T) is HarmonicDifferentials(T)
            True

        """
        return super().__classcall__(cls, surface, HarmonicDifferentials._homology_generators(surface, safety), singularities, category or SetsWithPartialMaps())

    def __init__(self, surface, homology_generators, singularities, category):
        Parent.__init__(self, category=category)

        self._surface = surface
        # TODO: Find a better way to track the L2 circles.
        self._debugs = []

        self._geometry = GeometricPrimitives(surface, homology_generators, singularities=singularities)

    @staticmethod
    def _homology_generators(surface, safety=None):
        # The generators of homology will come from paths crossing from a
        # center of a Delaunay cell to the center of a neighbouring Delaunay
        # cell.
        voronoi_paths = set()
        for label, edge in surface.edges():
            if surface.opposite_edge(label, edge) not in voronoi_paths:
                voronoi_paths.add((label, edge))

        # We now subdivide these generators to attain the required safety.
        # TODO: For now, we just add two more equally spaced points to the path.
        from sage.all import QQ
        gens = []
        if safety is None:
            for path in voronoi_paths:
                # TODO: Hardcoded for the octagon
                gens.append((path, QQ(0), QQ(274)/964))
                gens.append((path, QQ(274)/964, QQ(427)/964))
                gens.append((path, QQ(427)/964, QQ(537)/964))
                gens.append((path, QQ(537)/964, QQ(690)/964))
                gens.append((path, QQ(690)/964, QQ(964)/964))
        elif safety == QQ(1)/3:
            for path in voronoi_paths:
                # TODO: We are not actually doing safety 1/3 here currently but just break paths up evenly in thirds.
                gens.append((path, QQ(0), QQ(1)/3))
                gens.append((path, QQ(1)/3, QQ(2)/3))
                gens.append((path, QQ(2)/3, QQ(1)))
        else:
            # TODO: Currently, we ignore the safety.
            for path in voronoi_paths:
                gens.append((path, QQ(0), QQ(1)))

        gens = tuple(gens)

        return gens

    def plot(self):
        raise NotImplementedError
        # from flatsurf import TranslationSurface
        # S = self._surface
        # G = S.plot()
        # SR = PowerSeriesConstraints(S, prec=20, geometry=self._geometry).symbolic_ring()
        # for (label, edge, pos) in SR._gens:
        #     label, center = self._geometry.point_on_path_between_centers(label, edge, pos, wrap=True)
        #     radius = self._geometry._convergence(label, edge, pos)
        #     P = TranslationSurface(S).surface_point(label, center)
        #     G += P.plot(color="red")

        #     from sage.all import circle
        #     G += circle(P.graphical_surface_point().points()[0], radius, color="green", fill="green", alpha=.05)

        # for (label, edge, pos) in SR._gens:
        #     label, center = self._geometry.point_on_path_between_centers(label, edge, pos, wrap=True)

        #     P = TranslationSurface(S).surface_point(label, center)
        #     p = P.graphical_surface_point().points()[0]

        #     for (lbl, a_edge, a, b_edge, b, Δ0, Δ1, radius) in self._debugs:
        #         if lbl == label and a_edge == edge and a == pos:
        #             Δ = Δ0
        #         elif self._surface.opposite_edge(lbl, a_edge) == (label, edge) and a == 1 - pos:
        #             Δ = Δ0
        #         elif lbl == label and b_edge == edge and b == pos:
        #             Δ = Δ1
        #         elif self._surface.opposite_edge(lbl, b_edge) == (label, edge) and b == 1 - pos:
        #             Δ = Δ1
        #         else:
        #             continue

        #         from sage.all import vector
        #         q = p + vector((Δ.real(), Δ.imag()))

        #         from sage.all import line
        #         G += line((p, q), color="black")

        #         from sage.all import circle
        #         G += circle(q, radius, color="brown", fill="brown", alpha=.1)

        # return G

    def surface(self):
        return self._surface

    def _repr_(self):
        return f"Ω({self._surface})"

    def _element_constructor_(self, x, *args, **kwargs):
        if not x:
            return self.element_class(self, None, *args, **kwargs)

        if isinstance(x, dict):
            return self.element_class(self, x, *args, **kwargs)

        return self._element_from_cohomology(x, *args, **kwargs)

    def _element_from_cohomology(self, cocycle, /, prec=5, algorithm=["L2"], check=True):
        # TODO: In practice we could speed things up a lot with some smarter
        # caching. A lot of the quantities used in the computations only depend
        # on the surface & precision. When computing things for many cocycles
        # we do not need to recompute them. (But currently, we probably do
        # because they live in the constraints instance.)

        # We develop a consistent system of power series at each vertex of the Voronoi diagram
        # to describe a differential.

        # Let η be the differential we are looking for. To describe η we will use ω, the differential
        # corresponding to the flat structure given by our triangulation on this Riemann surface.
        # Then f:=η/ω is a meromorphic function which we can develop locally into a Laurent series.
        # Away from the vertices of the triangulation, ω has no zeros, so f has no poles there and is
        # thus given by a power series.

        # At each vertex of the Voronoi diagram, write f=Σ a_k z^k + O(z^prec). Our task is now to determine
        # the a_k.

        constraints = PowerSeriesConstraints(self.surface(), prec=prec, geometry=self._geometry)
        self._constraints = constraints # TODO: Remove

        # We use a variety of constraints. Which ones to use exactly is
        # determined by the "algorithm" parameter. If algorithm is a dict, it
        # can be used to configure aspects of the constraints.
        def get_parameter(alg, default):
            assert alg in algorithm
            if isinstance(algorithm, dict):
                return algorithm[alg]
            return default

        # (1) The radius of convergence of the power series is the distance from the vertex of the Voronoi
        # cell to the closest vertex of the triangulation (since we use a Delaunay triangulation, all vertices
        # are at the same distance in fact.) So the radii of convergence of two neigbhouring cells overlap
        # and the power series must coincide there. Note that this constraint is unrelated to the cohomology
        # class Φ.
        if "midpoint_derivatives" in algorithm:
            derivatives = get_parameter("midpoint_derivatives", prec//3)
            algorithm = [a for a in algorithm if a != "midpoint_derivatives"]
            constraints.require_midpoint_derivatives(derivatives)

        # (1') TODO: Describe L2 optimization.
        if "L2" in algorithm:
            weight = get_parameter("L2", 1)
            algorithm = [a for a in algorithm if a != "L2"]
            constraints.optimize(weight * constraints._L2_consistency())
            self._debugs = constraints._debugs

        if "squares" in algorithm:
            weight = get_parameter("squares", 1)
            algorithm = [a for a in algorithm if a != "squares"]
            constraints.optimize(weight * constraints._squares())

        # (2) We have that for any cycle γ, Re(∫fω) = Re(∫η) = Φ(γ). We can turn this into constraints
        # on the coefficients as we integrate numerically following the path γ as it intersects the radii of
        # convergence.
        constraints.require_cohomology(cocycle)

        # (3) Since the area ∫ η \wedge \overline{η} must be finite [TODO:
        # REFERENCE?] we optimize for a proxy of this quantity to be minimal.
        if "area_upper_bound" in algorithm:
            weight = get_parameter("area_upper_bound", 1)
            algorithm = [a for a in algorithm if a != "area_upper_bound"]
            constraints.optimize(weight * constraints._area_upper_bound())

        # (3') We can also optimize for the exact quantity to be minimal but
        # this is much slower.
        if "area" in algorithm:
            weight = get_parameter("area", 1)
            algorithm = [a for a in algorithm if a != "area"]
            constraints.optimize(weight * constraints._area())

        if "tykhonov" in algorithm:
            # TODO: Should we still try to do something like this? (Whatever
            # the idea was here?)
            pass

        if algorithm:
            raise ValueError(f"unsupported algorithm {algorithm}")

        solution, residue = constraints.solve()
        η = self.element_class(self, solution, residue=residue, cocycle=cocycle)

        if check:
            if report := η.error():
                raise ValueError(report)

        return η


class GeometricPrimitives(UniqueRepresentation):
    def __init__(self, surface, homology_generators, singularities=False):  # TODO: Make True the default everywhere if this works out.
        # TODO: Require immutable.
        self._surface = surface
        self._homology_generators = homology_generators
        self._singularities = singularities

    @cached_method
    def midpoint_on_path_between_centers(self, label, a_edge, a, b_edge, b):
        r"""
        Return the weighed midpoint of ``center + a * e`` and ``center + b *
        f`` where ``center`` is the center of the circumscribing circle of
        polygon ``label`` and ``e``/``f`` is the straight segment that goes
        from the ``center`` to the midpoint of the polygon across
        ``a_edge``/``b_edge``.

        The midpoint is determined by weighing the two points according to
        their radius of convergence, see :meth:`_convergence`.

        The midpoint is returned as the vector going from the first center,
        i.e., relative to ``center + a * e``.


        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import GeometricPrimitives
            sage: G = GeometricPrimitives(T, None)  # TODO: Should not be None

        Note that the midpoints in this example are not half way between the
        centers because the radius of convergence is not the same (the vertex
        of the square is understood to be a singularity)::

            sage: G.midpoint_on_path_between_centers(0, 0, 0, 0, 1/3)
            (0.000000000000000, -0.142350327708281)
            sage: G.midpoint_on_path_between_centers(0, 2, 2/3, 2, 1)
            (0.000000000000000, 0.190983005625053)
            sage: G.midpoint_on_path_between_centers(0, 2, 2/3, 2, 1) - G.midpoint_on_path_between_centers(0, 0, 0, 0, 1/3)
            (0.000000000000000, 0.333333333333333)

        Here the midpoints are exactly on the edge of the square::

            sage: G.midpoint_on_path_between_centers(0, 0, 1/3, 0, 2/3)
            (0.000000000000000, -0.166666666666667)
            sage: G.midpoint_on_path_between_centers(0, 2, 1/3, 2, 2/3)
            (0.000000000000000, 0.166666666666667)

        Again, the same as in the first example::

            sage: G.midpoint_on_path_between_centers(0, 0, 2/3, 0, 1)
            (0.000000000000000, -0.190983005625053)
            sage: G.midpoint_on_path_between_centers(0, 2, 0, 2, 1/3)
            (0.000000000000000, 0.142350327708281)

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.regular_octagon()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import GeometricPrimitives
            sage: G = GeometricPrimitives(T, None)  # TODO: Should not be None

            sage: G.midpoint_on_path_between_centers(0, 0, 0, 0, 1/2)  # tol 1e-9
            (0.000000000000000, -0.334089318959649)

        A midpoint between two different Voronoi paths::

            sage: G.midpoint_on_path_between_centers(0, 0, 1, 2, 1)
            (1.20710678118655, 1.20710678118655)
            sage: G.midpoint_on_path_between_centers(0, 0, 1, 4, 1)
            (0.000000000000000, 2.41421356237309)

        """
        radii = (
            self._convergence(label, a_edge, a),
            self._convergence(label, b_edge, b),
        )

        centers = (
            self.point_on_path_between_centers(label, a_edge, a)[1],
            self.point_on_path_between_centers(label, b_edge, b)[1],
        )

        return (centers[1] - centers[0]) * radii[1] / (sum(radii))

    @cached_method
    def point_on_path_between_centers(self, label, edge, pos, wrap=False, ring=None):
        r"""
        Return the point at ``center + pos * e`` where ``center`` is the center
        of the circumscribing circle of the polygon ``label`` and ``e`` is the
        straight segment connecting that center to the center of the polygon
        across the ``edge``.

        The point returned might not be inside the polygon ``label`` anymore.

        When ``wrap`` is not set, the point is returned in coordinates of the
        polygon ``label`` even if it is outside of the that polygon.

        Otherwise, the point is returned in coordinates of a polygon that
        contains it.

        OUTPUT:

        A pair ``(label, coordinates)`` where ``coordinates`` are coordinates
        of the point in the polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import GeometricPrimitives
            sage: G = GeometricPrimitives(T, None)  # TODO: Should not be None

            sage: G.point_on_path_between_centers(0, 0, 0)
            (0, (1/2, 1/2))
            sage: G.point_on_path_between_centers(0, 0, 1/2)
            (0, (1/2, 0))
            sage: G.point_on_path_between_centers(0, 0, 1)
            (0, (1/2, -1/2))
            sage: G.point_on_path_between_centers(0, 0, 1, wrap=True)
            (0, (1/2, 1/2))

        """
        # TODO: This method feels a bit hacky. The name is misleading, the return tuple is weird, and the ring argument is a strange performance hack.
        polygon = self._surface.polygon(label)
        polygon_center = polygon.circumscribing_circle().center()

        opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)
        opposite_polygon = self._surface.polygon(opposite_label)
        opposite_polygon = opposite_polygon.translate(-opposite_polygon.vertex((opposite_edge + 1) % len(opposite_polygon.vertices())) + polygon.vertex(edge))

        assert polygon.vertex((edge + 1) % len(polygon.vertices())) == opposite_polygon.vertex(opposite_edge)

        opposite_polygon_center = opposite_polygon.circumscribing_circle().center()

        center = (1 - pos) * polygon_center + pos * opposite_polygon_center

        if not wrap or self._surface.polygon(label).contains_point(center):
            if ring is not None:
                from sage.all import vector
                center = vector((ring(center[0]), ring(center[1])))
            return label, center

        label, edge = self._surface.opposite_edge(label, edge)
        return self.point_on_path_between_centers(label, edge, 1 - pos, wrap=True, ring=ring)

    @cached_method
    def midpoint_on_center_vertex_path(self, label, vertex):
        return (vertex - self._surface.polygon(label).circumscribing_circle().center()) / 2

    @cached_method
    def _convergence(self, label, edge, pos):
        r"""
        Return the radius of convergence at the point at which we develop a
        series around the point that is at ``center + pos * e`` where ``center``
        is the center of the circumscribing circle of the polygon ``label`` and
        ``e`` is the straight segment that goes from that point to the center
        for the polygon across ``edge``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.regular_octagon()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)

            sage: Ω._geometry._convergence(0, 0, 0)
            1.30656296487638
            sage: Ω._geometry._convergence(0, 0, 1/3)
            0.641794946587438
            sage: Ω._geometry._convergence(0, 0, 1/2)
            0.500000000000000
            sage: Ω._geometry._convergence(0, 0, 2/3)
            0.641794946587438
            sage: Ω._geometry._convergence(0, 0, 1)
            1.30656296487638

        """
        label, center = self.point_on_path_between_centers(label, edge, pos)
        polygon = self._surface.polygon(label)

        # TODO: We are assuming that the only relevant singularities are the
        # end points of ``edge``.
        # TODO: Use a ring with more appropriate precision.
        from sage.all import RR
        return min(
            (center - polygon.vertex(edge)).change_ring(RR).norm(),
            (center - polygon.vertex((edge + 1) % len(polygon.vertices()))).change_ring(RR).norm())

    def branch(self, center, Δ):
        center_point = self._surface(*center)
        if not center_point.is_vertex():
            return 0, 0

        # TODO: Hardcoded for octagon
        low = Δ[1] < center[1][1]

        for i, vertex in enumerate(self._surface.polygon(center[0]).vertices()):
            if vertex == center[1]:
                if i == 0:
                    return 0, 2
                if i == 1:
                    return 1, 2
                if i == 2:
                    return (0 if low else 2), 2
                if i == 3:
                    return (1 if low else 0), 2
                if i == 4:
                    return 2, 2
                if i == 5:
                    return 0, 2
                if i == 6:
                    return 1, 2
                if i == 7:
                    return 2, 2

        raise NotImplementedError


class SymbolicCoefficientExpression(CommutativeRingElement):
    # TODO: Make sure that we never have zero coefficients as these would break degree computations.

    def __init__(self, parent, coefficients):
        super().__init__(parent)

        self._coefficients = coefficients

    def _richcmp_(self, other, op):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))

        Due to the internal ordering of the generators, this generator prints
        as ``a2``. The generator ``a0`` is used for the power series developed
        at some other point::

            sage: a
            Im(a2,0)
            sage: b
            Re(a2,0)
            sage: a + b
            Re(a2,0) + Im(a2,0)
            sage: a + b + 1
            Re(a2,0) + Im(a2,0) + 1.00000000000000

        """
        if self.is_constant():
            return repr(self.constant_coefficient())

        def decode(gen):
            if gen < 0:
                return -gen-1,

            for i, g in enumerate(self.parent()._regular_gens()):
                if i == gen:
                    return g

        def variable_name(gen):
            gen = decode(gen)

            if len(gen) == 1:
                return f"λ{gen[0]}"

            if len(gen) == 3:
                kind, point, k = gen
                if k < 0:
                    k = "__minus__" + str(-k)
                index = self.parent()._centers.index(point)
                return f"{kind}__open__a{index}__comma__{k}__close__"

            raise NotImplementedError

        def key(gen):
            gen = decode(gen)

            if len(gen) == 1:
                n = gen[0]
                return 1e9, n

            if len(gen) == 3:
                kind, point, k = gen
                return self.parent()._centers.index(point), k, 0 if kind == "Re" else 1

            kind, label, edge, pos, k = gen
            return label, edge, pos, k, 0 if kind == "Re" else 1

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

        return repr(f).replace('__open__', '(').replace('__close__', ')').replace('__comma__', ',').replace("__minus__", "-")

    def degree(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: x = R.gen(('Im', T(0, (1/2, 1/2)), 0)) + 1; x
            Im(a2,0) + 1.00000000000000
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: -a
            -Im(a2,0)
            sage: -(a + b)
            -Re(a2,0) - Im(a2,0)
            sage: -(a * a)
            -Im(a2,0)^2
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: a + 1
            Im(a2,0) + 1.00000000000000
            sage: a + (-a)
            0.000000000000000
            sage: a + b
            Re(a2,0) + Im(a2,0)
            sage: a * a + b * b
            Re(a2,0)^2 + Im(a2,0)^2

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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: a - 1
            Im(a2,0) - 1.00000000000000
            sage: a - a
            0.000000000000000
            sage: a * a - b * b
            -Re(a2,0)^2 + Im(a2,0)^2

        """
        return self._add_(-other)

    def _mul_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: a * a
            Im(a2,0)^2
            sage: a * b
            Re(a2,0)*Im(a2,0)
            sage: a * R.one()
            Im(a2,0)
            sage: a * R.zero()
            0.000000000000000
            sage: (a + b) * (a - b)
            -Re(a2,0)^2 + Im(a2,0)^2

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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: a * 0
            0.000000000000000
            sage: a * 1
            Im(a2,0)
            sage: a * 2
            2.00000000000000*Im(a2,0)

        """
        return self._lmul_(right)

    def _lmul_(self, left):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: 0 * a
            0.000000000000000
            sage: 1 * a
            Im(a2,0)
            sage: 2 * a
            2.00000000000000*Im(a2,0)

        """
        return type(self)(self.parent(), {key: coefficient for (key, value) in self._coefficients.items() if (coefficient := left * value)})

    def constant_coefficient(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: a.variables()
            {Im(a2,0)}
            sage: (a + b).variables()
            {Im(a2,0), Re(a2,0)}
            sage: (a * a).variables()
            {Im(a2,0)}
            sage: (a * b).variables()
            {Im(a2,0), Re(a2,0)}
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

    def is_variable(self):
        r"""
        Return the label of the polygon affected by this variable.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: a.describe()
            ('Im', Point (1/2, 1/2) of polygon 0, 0)
            sage: b.describe()
            ('Re', Point (1/2, 1/2) of polygon 0, 0)
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

        for i, g in enumerate(self.parent()._regular_gens()):
            if i == variable:
                return g

    def real(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: c = (a + b)**2
            sage: c.real()
            Re(a2,0)^2 + 2.00000000000000*Re(a2,0)*Im(a2,0) + Im(a2,0)^2

        """
        return self.map_coefficients(lambda c: c.real(), self.parent().change_ring(self.parent().real_field()))

    def imag(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
            sage: c = (a + b)**2
            sage: c.imag()
            0.000000000000000

            sage: c = (I*a + b)**2
            sage: c.imag()
            2.00000000000000*Re(a2,0)*Im(a2,0)

        """
        return self.map_coefficients(lambda c: c.imag(), self.parent().change_ring(self.parent().real_field()))

    def __getitem__(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
            2.00000000000000*Im(a2,0)
            sage: c.derivative(b)
            -2.00000000000000*Re(a2,0)

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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: a = R.gen(('Im', T(0, (1/2, 1/2)), 0))
            sage: b = R.gen(('Re', T(0, (1/2, 1/2)), 0))
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
    def __classcall__(cls, surface, base_ring, geometry, category=None):
        from sage.categories.all import CommutativeRings
        return super().__classcall__(cls, surface, base_ring, geometry, category or CommutativeRings())

    def __init__(self, surface, base_ring, geometry, category):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry)
            sage: R.has_coerce_map_from(CC)
            True

            sage: TestSuite(R).run()

        """
        self._surface = surface
        self._base_ring = base_ring
        self._geometry = geometry
        self._homology_generators = geometry._homology_generators

        self._centers = set()
        for (label, edge), a, b in self._homology_generators:
            self._centers.add((label, edge if a else 0, a))
            if b == 1:
                label, edge = self._surface.opposite_edge(label, edge)
                edge = 0
                b = 0
            self._centers.add((label, edge, b))

        self._centers = [self._surface(*self._geometry.point_on_path_between_centers(*gen, wrap=True)) for gen in self._centers]

        if self._geometry._singularities:
            for vertex in self._surface.vertices():
                self._centers.append(vertex)

        CommutativeRing.__init__(self, base_ring, category=category, normalize=False)
        self.register_coercion(base_ring)

    Element = SymbolicCoefficientExpression

    def _repr_(self):
        return f"Ring of Power Series Coefficients over {self.base_ring()}"

    def change_ring(self, ring):
        return SymbolicCoefficientRing(self._surface, ring, self._geometry, category=self.category())

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
        from sage.all import parent, ZZ
        if parent(n) == ZZ:
            n = int(n)

        if isinstance(n, int):
            return self((n,))

        if isinstance(n, tuple):
            if len(n) == 2 and n[0] == "lagrange":
                if n[1] < 0:
                    raise ValueError

                return self.gen(-n[1] - 1)

            if len(n) == 3:
                kind, center, k = n

                for i, gen in enumerate(self._regular_gens()):
                    if gen == n:
                        return self.gen(i)

        raise ValueError(f"symbolic ring has no generator {n}")

    # TODO: Use fixed prec from constructor as default
    @cached_method
    def _regular_gens(self, prec=64):
        def __regular_gens():
            if len(self._surface.angles()) != 1:
                raise NotImplementedError

            d = self._surface.angles()[0] - 1

            # Generate the non-negative coefficients. At singularities we have to create multiple coefficients.
            for i in range(prec):
                for c in self._centers:
                    if c.is_vertex():
                        for j in range(d + 1):
                            yield ("Re", c, i * (d + 1) + j)
                            yield ("Im", c, i * (d + 1) + j)
                    else:
                        yield ("Re", c, i)
                        yield ("Im", c, i)
        return [gen for gen in __regular_gens()]

    def ngens(self):
        raise NotImplementedError


class PowerSeriesConstraints:
    r"""
    A collection of (linear) constraints on the coefficients of power series
    developed at the vertices of the Voronoi cells of a Delaunay triangulation.

    This is used to create harmonic differentials from cohomology classes.
    """

    def __init__(self, surface, prec, geometry, bitprec=None):
        from sage.all import log, ceil, factorial

        self._surface = surface
        self._prec = prec
        # TODO: The default value tends to be gigantic!
        self._bitprec = bitprec or ceil(log(factorial(prec), 2) + 53)
        self._geometry = geometry
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

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: Ω = HarmonicDifferentials(T)

            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)
            sage: C.symbolic_ring()
            Ring of Power Series Coefficients over Complex Field with 54 bits of precision

        """
        # TODO: What's the correct precision here?
        # return SymbolicCoefficientRing(self._surface, base_ring=base_ring or self.complex_field())
        from sage.all import ComplexField
        return SymbolicCoefficientRing(self._surface, base_ring=base_ring or ComplexField(54), geometry=self._geometry)

    @cached_method
    def complex_field(self):
        from sage.all import ComplexField
        return ComplexField(54)
        return ComplexField(self._bitprec)

    @cached_method
    def real_field(self):
        from sage.all import RealField
        return RealField(54)
        return RealField(self._bitprec)

    def gen(self, point, /, conjugate=False):
        raise NotImplementedError

    @cached_method
    def _gen_nonsingular(self, label, edge, pos, k, /, conjugate=False):
        # TODO: Remove this method
        assert conjugate is True or conjugate is False
        real = self._real_nonsingular(label, edge, pos, k)
        imag = self._imag_nonsingular(label, edge, pos, k)

        i = self.symbolic_ring().imaginary_unit()

        if conjugate:
            i = -i

        return real + i*imag

    @cached_method
    def _real_nonsingular(self, label, edge, pos, k):
        r"""
        Return the real part of the kth coefficient of the power series
        developed around ``center + pos * e`` where ``center`` is the center of
        the circumscribing circle of the polygon ``label`` and ``e`` is the
        straight segment connecting that center to the center of the polygon
        across the ``edge``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: Ω = HarmonicDifferentials(T)

            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)
            sage: C._real_nonsingular(0, 0, 0, 0)
            Re(a2,0)
            sage: C._real_nonsingular(0, 0, 0, 1)
            Re(a2,1)
            sage: C._real_nonsingular(0, 0, 1, 0)
            Re(a2,0)
            sage: C._real_nonsingular(0, 1, 0, 0)
            Re(a2,0)
            sage: C._real_nonsingular(0, 2, 137/482, 0)
            Re(a0,0)
            sage: C._real_nonsingular(0, 0,  537/964, 0)
            Re(a3,0)
            sage: C._real_nonsingular(0, 1, 345/482, 0)
            Re(a1,0)
            sage: C._real_nonsingular(0, 1, 537/964, 0)
            Re(a4,0)

        """
        # TODO: Remove this method
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator at this point")

        return self.symbolic_ring().gen(("Re", self._surface(*self._geometry.point_on_path_between_centers(label, edge, pos, wrap=True)), k))

    @cached_method
    def _imag_nonsingular(self, label, edge, pos, k):
        r"""
        Return the imaginary part of the kth generator of the :meth:`symbolic_ring`
        for the polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)
            sage: C._imag_nonsingular(0, 0, 0, 0)
            Im(a2,0)
            sage: C._imag_nonsingular(0, 0, 0, 1)
            Im(a2,1)
            sage: C._imag_nonsingular(0, 0, 0, 2)
            Im(a2,2)

        """
        # TODO: Remove this method
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator at this point")

        return self.symbolic_ring().gen(("Im", self._surface(*self._geometry.point_on_path_between_centers(label, edge, pos, wrap=True)), k))

    @cached_method
    def lagrange(self, k):
        return self.symbolic_ring().gen(("lagrange", k))

    def project(self, x, part):
        r"""
        Return the ``"Re"`` or ``"Im"```inary ``part`` of ``x``.
        """
        if part not in ["Re", "Im"]:
            raise ValueError("part must be one of real or imag")

        if part == "Re":
            part = "real"
        if part == "Im":
            part = "imag"

        # Return the real/imaginary part of a complex number.
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

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)

            sage: C.real_part(1 + I)  # tol 1e-9
            1
            sage: C.real_part(1.)  # tol 1e-9
            1

        ::

            sage: C.real_part(C._gen_nonsingular(0, 0, 0, 0))
            Re(a2,0)
            sage: C.real_part(C._real_nonsingular(0, 0, 0, 0))
            Re(a2,0)
            sage: C.real_part(C._imag_nonsingular(0, 0, 0, 0))
            Im(a2,0)
            sage: C.real_part(2*C._gen_nonsingular(0, 0, 0, 0))  # tol 1e-9
            2*Re(a2,0)
            sage: C.real_part(2*I*C._gen_nonsingular(0, 0, 0, 0))  # tol 1e-9
            -2.0000000000000*Im(a2,0)

        """
        return self.project(x, "Re")

    def imaginary_part(self, x):
        r"""
        Return the imaginary part of ``x``.

        EXAMPLES::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)

            sage: C.imaginary_part(1 + I)  # tol 1e-9
            1
            sage: C.imaginary_part(1.)  # tol 1e-9
            0

        ::

            sage: C.imaginary_part(C._gen_nonsingular(0, 0, 0, 0))
            Im(a2,0)
            sage: C.imaginary_part(C._real_nonsingular(0, 0, 0, 0))  # tol 1e-9
            0
            sage: C.imaginary_part(C._imag_nonsingular(0, 0, 0, 0))  # tol 1e-9
            0
            sage: C.imaginary_part(2*C._gen_nonsingular(0, 0, 0, 0))  # tol 1e-9
            2*Im(a2,0)
            sage: C.imaginary_part(2*I*C._gen_nonsingular(0, 0, 0, 0))  # tol 1e-9
            2*Re(a2,0)

        """
        return self.project(x, "Im")

    def add_constraint(self, expression, rank_check=True):
        total_degree = expression.total_degree()

        if total_degree == -1:
            return

        if total_degree == 0:
            raise ValueError(f"cannot solve for constraint {expression} == 0")

        if total_degree > 1:
            raise NotImplementedError("can only encode linear constraints")

        if expression.parent().base_ring() is self.real_field():
            # TODO: Should we scale?
            # self._constraints.append(expression / expression.norm(1))
            # TODO: This is a very expensive hack to detect dependent conditions.
            if rank_check:
                rank = self.matrix(nowarn=True)[0].rank()
            self._constraints.append(expression)
            if rank_check:
                if self.matrix(nowarn=True)[0].rank() == rank:
                    self._constraints.pop()
        elif expression.parent().base_ring() is self.complex_field():
            self.add_constraint(expression.real(), rank_check=rank_check)
            self.add_constraint(expression.imag(), rank_check=rank_check)
        else:
            raise NotImplementedError("cannot handle expressions over this base ring")

    @cached_method
    def _formal_power_series_nonsingular(self, label, edge, pos, base_ring=None):
        if base_ring is None:
            base_ring = self.symbolic_ring()

        from sage.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, 'z')

        return R([self._gen_nonsingular(label, edge, pos, n) for n in range(self._prec)])

    def develop_nonsingular(self, label, edge, pos, Δ=0, base_ring=None):
        r"""
        Return the power series obtained by developing at z + Δ where z is the
        coordinate at ``center + pos * e`` where ``center`` is the center of
        the circumscribing circle of ``label`` and ``e`` is the straight
        segment connecting that center to the one across the ``edge``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)
            sage: C.develop_nonsingular(0, 0, 0)  # tol 1e-9
            Re(a2,0) + 1.000000000000000*I*Im(a2,0) + (Re(a2,1) + 1.000000000000000*I*Im(a2,1))*z + (Re(a2,2) + 1.000000000000000*I*Im(a2,2))*z^2
            sage: C.develop_nonsingular(0, 0, 0, 1)  # tol 1e-9
            Re(a2,0) + 1.000000000000000*I*Im(a2,0) + Re(a2,1) + 1.000000000000000*I*Im(a2,1) + Re(a2,2) + 1.000000000000000*I*Im(a2,2) + (Re(a2,1) + 1.000000000000000*I*Im(a2,1) + 2.000000000000000*Re(a2,2) + 2.000000000000000*I*Im(a2,2))*z + (Re(a2,2) + 1.000000000000000*I*Im(a2,2))*z^2

        """
        # TODO: Check that Δ is within the radius of convergence.
        f = self._formal_power_series_nonsingular(label, edge, pos, base_ring=base_ring)
        return f(f.parent().gen() + Δ)

    def develop_singular(self, label, vertex, Δ, radius):
        raise NotImplementedError

    def integrate(self, cycle):
        r"""
        Return the linear combination of the power series coefficients that
        describe the integral of a differential along the homology class
        ``cycle``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=1, geometry=Ω._geometry)

            sage: C.integrate(H())
            0.000000000000000

            sage: a, b = H.gens()
            sage: C.integrate(a)  # tol 1e-9
            0.236797907979275*Re(a1,0) + 0.236797907979275*I*Im(a1,0) + 0.247323113451507*Re(a2,0) + 0.247323113451507*I*Im(a2,0) + 0.139540535294972*Re(a4,0) + 0.139540535294972*I*Im(a4,0) + 0.236797907979275*Re(a6,0) + 0.236797907979275*I*Im(a6,0) + 0.139540535294972*Re(a7,0) + 0.139540535294972*I*Im(a7,0)

            sage: C.integrate(b)  # tol 1e-9
            (-0.236797907979275*I)*Re(a0,0) + 0.236797907979275*Im(a0,0) + (-0.247323113451507*I)*Re(a2,0) + 0.247323113451507*Im(a2,0) + (-0.139540535294972*I)*Re(a3,0) + 0.139540535294972*Im(a3,0) + (-0.139540535294972*I)*Re(a5,0) + 0.139540535294972*Im(a5,0) + (-0.236797907979275*I)*Re(a8,0) + 0.236797907979275*Im(a8,0)

            sage: C = PowerSeriesConstraints(T, prec=5, geometry=Ω._geometry)
            sage: C.integrate(a) + C.integrate(-a)  # tol 1e-9
            0.00000000000000000
            sage: C.integrate(b) + C.integrate(-b)  # tol 1e-9
            0.00000000000000000

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=1, geometry=Ω._geometry)

            sage: C.integrate(H())
            0.000000000000000

            sage: a, b, c, d = H.gens()
            sage: C.integrate(b)
            1.41421356237310*Re(a0,0) + 1.41421356237310*I*Im(a0,0) + (-0.389300863573646 + 1.38777878078145e-17*I)*Re(a1,0) + (-1.38777878078145e-17 - 0.389300863573646*I)*Im(a1,0) + (-1.04083408558608e-17 - 0.327551243899061*I)*Re(a1,1) + (0.327551243899061 - 1.04083408558608e-17*I)*Im(a1,1) + 0.270679432377470*Re(a1,2) + 0.270679432377470*I*Im(a1,2)

            sage: C.integrate(d)
            (-1.41421356237310*I)*Re(a0,0) + 1.41421356237310*Im(a0,0) + (5.55111512312578e-17 - 0.389300863573646*I)*Re(a1,0) + (0.389300863573646 + 5.55111512312578e-17*I)*Im(a1,0) + (-2.77555756156289e-17 + 0.327551243899061*I)*Re(a1,1) + (-0.327551243899061 - 2.77555756156289e-17*I)*Im(a1,1) + (-0.270679432377470*I)*Re(a1,2) + 0.270679432377470*Im(a1,2)

            sage: C.integrate(a)
            (1.00000000000000 - 1.00000000000000*I)*Re(a0,0) + (1.00000000000000 + 1.00000000000000*I)*Im(a0,0) + (0.275277280554704 + 0.275277280554704*I)*Re(a1,0) + (-0.275277280554704 + 0.275277280554704*I)*Im(a1,0) + (0.327551243899061 + 6.93889390390723e-18*I)*Re(a1,1) + (-6.93889390390723e-18 + 0.327551243899061*I)*Im(a1,1) + (0.191399262161835 - 0.191399262161835*I)*Re(a1,2) + (0.191399262161835 + 0.191399262161835*I)*Im(a1,2)

        """
        surface = cycle.surface()

        R = self.symbolic_ring()

        expression = R.zero()

        for (label, edge), multiplicity in cycle._chain.monomial_coefficients().items():
            # Integrate along the path crossing the edge into the opposite face.
            pos = 0
            while pos != 1:
                for gen in self._geometry._homology_generators:
                    if gen[0] == (label, edge) and gen[1] == pos:
                        break
                else:
                    assert False, f"cannot continue path from {label}, {edge}, {pos} with generators {self._geometry._homology_generators}"

                if not self._geometry._singularities:
                    # We develop exactly around the endpoints of the generators of homology.
                    a = gen[1]
                    b = gen[2]
                    P = self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, edge, a, edge, b))

                    P_power = P

                    for k in range(self._prec):
                        expression += multiplicity * self._gen_nonsingular(label, edge, a, k) / (k + 1) * P_power
                        P_power *= P

                    opposite_label, opposite_edge = surface.opposite_edge(label, edge)

                    Q = self.complex_field()(*self._geometry.midpoint_on_path_between_centers(opposite_label, opposite_edge, 1 - b, opposite_edge, 1 - a))

                    Q_power = Q

                    for k in range(self._prec):
                        expression -= multiplicity * self._gen_nonsingular(opposite_label, opposite_edge, 1-b, k) / (k + 1) * Q_power

                        Q_power *= Q
                else:
                    if str(surface) != "Translation Surface in H_2(2) built from a regular octagon":
                        raise NotImplementedError

                    assert label == 0

                    if edge == 0:
                        vertex = 0
                    elif edge == 1:
                        vertex = 1
                    elif edge == 2:
                        vertex = 2
                    elif edge == 3:
                        vertex = 4
                    else:
                        raise NotImplementedError

                    octagon = surface.polygon(0)
                    width = octagon.edges()[0][0] + 2 * octagon.edges()[1][0]

                    # Probably not exactly true but close enough.
                    width_on_central_chart = width - octagon.edges()[0][0]
                    width_on_vertex_chart = octagon.edges()[0][0]

                    # First: Integrate along the line crossing over the center of the
                    # octagon from -δ v to +δ v where δ = width_on_central_chart/2.
                    δ = width_on_central_chart / 2
                    v = surface.polygon(label).vertices()[edge] + surface.polygon(label).edges()[edge] / 2 - octagon.circumscribing_circle().center()
                    v_ = v.change_ring(self.real_field())
                    v_ /= v_.norm()

                    # We integrate the summands of the power series \sum a_k z^k.
                    # We have \int_γ a_k z^k = \int_{-δ}^δ a_k γ(t)^k ·γ(t) dt = a_k v^{k + 1} \int_{-δ}^δ t^k dt = a_k v^{k + 1}/(k+1) (δ^{k + 1} - (-δ)^{k + 1}).
                    for k in range(self._prec):
                        expression += multiplicity * self._gen_nonsingular(0, 0, 0, k) * self.complex_field()(*v_)**(k + 1) / (k + 1) * (δ**(k + 1) - (-δ)**(k + 1))

                    # Second: Integrate using the power series at the singularity.
                    # Find segments that describe the path.
                    segments = self._voronoi_diagram().cell(surface(0, 0))
                    from flatsurf.geometry.euclidean import is_parallel
                    segments = [segment for segment in segments if is_parallel(segment.segment()[1].vector(), v)]
                    assert len(segments) == 2, "this is not the octagon"

                    for segment in segments:
                        integrator = self.Integrator(self, segment)
                        # As described in _L2_consistency_voronoi_boundary, we integrate
                        # f = Σ_{n ≥ 0} a_n f_n(z) along the segment γ.
                        # TODO: 3 is hardcoded for the octagon.
                        for n in range(3 * self._prec):
                            expression += multiplicity * integrator.a(n) * integrator.f(n)
                    break

                pos = b

        return expression

    def evaluate(self, label, edge, pos, Δ, derivative=0):
        r"""
        Return the value of the power series evaluated at Δ in terms of
        symbolic variables.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=3, geometry=Ω._geometry)
            sage: C.evaluate(0, 0, 0, 0)  # tol 1e-9
            Re(a2,0) + 1.00000000000000*I*Im(a2,0)
            sage: C.evaluate(0, 0, 0, 2)  # tol 1e-9
            Re(a2,0) + 1.000000000000000*I*Im(a2,0) + 2.000000000000000*Re(a2,1) + 2.000000000000000*I*Im(a2,1) + 4.000000000000000*Re(a2,2) + 4.000000000000000*I*Im(a2,2)

        """
        # TODO: Check that Δ is within the radius of convergence.

        if derivative >= self._prec:
            raise ValueError

        parent = self.symbolic_ring()

        value = parent.zero()

        z = self.complex_field().one()

        from sage.all import factorial
        # factor = self.complex_field()(factorial(derivative))
        from sage.all import ComplexField
        factor = ComplexField(54)(factorial(derivative))

        if edge is None and pos is None:
            raise NotImplementedError
            # TODO: Δ has a different meaning now.
            edge, pos, Δ = self.relativize(label, Δ)

        for k in range(derivative, self._prec):
            value += factor * self._gen_nonsingular(label, edge, pos, k) * z

            factor *= k + 1
            factor /= k - derivative + 1

            z *= Δ

        return value

    def value(self, point):
        # TODO: Make the chart a parameter.
        r"""
        Return f(point) for a differential f dz with z the flat coordinate of
        the polygon containing ``point``.

        This uses the power series that is best suited to describe the value of
        f at the point, e.g., when close to a singularity, it uses the power
        series there.

        TODO: The above statement does not really make sense.

        INPUT:

        - ``point`` -- a non-singular point

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, 2, Ω._geometry)

            sage: center = S(0, S.polygon(0).circumscribing_circle().center())
            sage: C.value(center)
            Re(a0,0) + 1.00000000000000*I*Im(a0,0)

            sage: off_center = S(0, S.polygon(0).circumscribing_circle().center() + vector((1/2, 0)))
            sage: C.value(off_center)
            Re(a0,0) + 1.00000000000000*I*Im(a0,0) + 0.500000000000000*Re(a0,1) + 0.500000000000000*I*Im(a0,1)

        """
        (label, center), Δ = self.relativize(point)

        Δ = self.complex_field()(*Δ)

        # See, _L2_consistency_voronoi_boundary for the notation: we compute
        # the value f(z) =
        # Σ_{n ≥ -d} a_n ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1}
        # at z = center + Δ, i.e., z-α=Δ.
        κ, d = self._geometry.branch((label, center), Δ)

        if Δ == 0 and d:
            raise ValueError("point must be non-singular")

        center_point = self._surface(label, center)

        expression = self.symbolic_ring().zero()
        for n in self._range(center_point, self._prec):
            an = self._gen("Re", center_point, n) + self.complex_field().gen() * self.symbolic_ring(self.complex_field())(self._gen("Im", center_point, n))
            expression += an * self.ζ(d + 1) ** (κ * (n + 1)) / (d + 1) * Δ.nth_root(d + 1)**(n - d)

        return expression

    @cached_method
    def _circumscribing_circle_center(self, label):
        from sage.all import CC
        return CC(*self._surface.polygon(label).circumscribing_circle().center())

    def relativize(self, point):
        r"""
        Determine the power series which is best suited to evaluate a function at ``point``.

        Returns the point at which the power series is developed as a pair
        (label, coordinates) and a vector Δ from that point to the ``point``.

        .. NOTE::

            We do not just return the point at which the power series is
            developed as a surface point since when that center is a
            singularity the information can be ambiguous.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, 1, Ω._geometry)

            sage: C.relativize(S(0, S.polygon(0).circumscribing_circle().center()))
            ((0, (1/2, 1/2*a + 1/2)), (0, 0))


            sage: C.relativize(S(0, S.polygon(0).circumscribing_circle().center() + vector((1/2, 1/2))))
            ((0, (1/2, 1/2*a + 1/2)), (1/2, 1/2))

            sage: C.relativize(S(0, S.polygon(0).circumscribing_circle().center() + vector((2/3, 2/3))))
            ((0, (1, a + 1)), (1/6, -1/2*a + 1/6))

            sage: C.relativize(S(0, S.polygon(0).circumscribing_circle().center() + vector((2/3, 2/3 - 1/1024))))
            ((0, (1/2*a + 1, 1/2*a + 1)), (-1/2*a + 1/6, 509/3072))

        """
        return self._voronoi_diagram().relativize(point)

    def require_midpoint_derivatives(self, derivatives):
        r"""
        Add constraints to verify that the value of the power series and its
        first ``derivatives`` derivatives are compatible.

        Namely, along a path that we use for integration, require that the
        functions coincide at a halfway point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

        We require the power series to be compatible with itself away from the
        midpoint. At this low precision, this just means that the constant
        coefficients need to match::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, prec=1, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) - Re(a6,0), Im(a2,0) - Im(a6,0), Re(a6,0) - Re(a7,0), Im(a6,0) - Im(a7,0), -Re(a4,0) + Re(a7,0), -Im(a4,0) + Im(a7,0), -Re(a1,0) + Re(a4,0), -Im(a1,0) + Im(a4,0), Re(a2,0) - Re(a8,0), Im(a2,0) - Im(a8,0), -Re(a5,0) + Re(a8,0), -Im(a5,0) + Im(a8,0), -Re(a3,0) + Re(a5,0), -Im(a3,0) + Im(a5,0), -Re(a0,0) + Re(a3,0), -Im(a0,0) + Im(a3,0)]


        If we add more coefficients, we get more complicated conditions::

            sage: C = PowerSeriesConstraints(T, prec=2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), Re(a0,0) + 0.160570808419475*Im(a0,1) - Re(a2,0) + 0.123661556725753*Im(a2,1), Im(a0,0) - 0.160570808419475*Re(a0,1) - Im(a2,0) - 0.123661556725753*Re(a2,1)]

        ::

            sage: C = PowerSeriesConstraints(T, prec=2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(2)
            sage: C  # tol 1e-9
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a2,1) - Re(a6,1), Im(a2,1) - Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a6,1) - Re(a7,1), Im(a6,1) - Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a4,1) + Re(a7,1), -Im(a4,1) + Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), -Re(a1,1) + Re(a4,1), -Im(a1,1) + Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a2,1) - Re(a8,1), Im(a2,1) - Im(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a5,1) + Re(a8,1), -Im(a5,1) + Im(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a3,1) + Re(a5,1), -Im(a3,1) + Im(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), -Re(a0,1) + Re(a3,1), -Im(a0,1) + Im(a3,1)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for (label, edge), a, b in self._geometry._homology_generators:
            opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)

            parent = self.symbolic_ring()

            Δ0 = self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, edge, a, edge, b))
            Δ1 = self.complex_field()(*self._geometry.midpoint_on_path_between_centers(opposite_label, opposite_edge, 1-b, opposite_edge, 1-a))

            # Require that the 0th, ..., derivatives-1th derivatives are the same at the midpoint of the edge.
            for derivative in range(derivatives):
                self.add_constraint(
                    parent(self.evaluate(label, edge, a, Δ0, derivative)) - parent(self.evaluate(opposite_label, opposite_edge, 1-b, Δ1, derivative)))

    def _L2_consistency_cost_ball(self, T0, T1, r, debug):
        cost = self.symbolic_ring(self.real_field()).zero()

        b = (T0 - T1).list()

        debug.append(r)

        assert r > 0
        r2 = self.real_field()(r**2)

        r2n = r2
        for n, b_n in enumerate(b):
            # TODO: In the article it says that it should be R^n as a
            # factor but R^{2n+2} is actually more reasonable. See
            # https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/Harmonic.20Differentials/near/308863913
            real = b_n.real()
            imag = b_n.imag()
            cost += (real * real + imag * imag) * r2n

            r2n *= r2

        # TODO: This is not in the article. Does that make sense?
        # Normalize the result with the area of the integral domain.
        cost /= r2

        return cost

    def _L2_consistency_between_nonsingular_points(self, label, a_edge, a, b_edge, b):
        debug = [label, a_edge, a, b_edge, b]

        def Δ(x_edge, x, y_edge, y):
            x_opposite_label, x_opposite_edge = self._surface.opposite_edge(label, x_edge)
            if x_opposite_label != label:
                raise NotImplementedError  # need to shift polygons

            y_opposite_label, y_opposite_edge = self._surface.opposite_edge(label, y_edge)
            if y_opposite_label != label:
                raise NotImplementedError  # need to shift polygons

            Δs = (
                self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, x_edge, x, y_edge, y)),
                self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, x_opposite_edge, 1 - x, y_edge, y)),
                self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, x_edge, x, y_opposite_edge, 1 - y)),
                self.complex_field()(*self._geometry.midpoint_on_path_between_centers(label, x_opposite_edge, 1 - x, y_opposite_edge, 1 - y)),
            )

            return min(Δs, key=lambda v: v.norm())

        # The weighed midpoint of the segment where the power series meet with
        # respect to the centers of the power series.
        # TODO: Assert that the min is attained for the same choice of representatives (or pick the representatives explicitly as the ones with minimal distance.)
        Δ0 = Δ(a_edge, a, b_edge, b)
        Δ1 = Δ(b_edge, b, a_edge, a)

        debug += [Δ0, Δ1]

        # Develop both power series around that midpoint, i.e., Taylor expand them.
        T0 = self.develop_nonsingular(label, a_edge, a, Δ0)
        T1 = self.develop_nonsingular(label, b_edge, b, Δ1)

        # Write b_n for the difference of the n-th coefficient of both power series.
        # We want to minimize the sum of |b_n|^2 r^2n where r is a somewhat
        # randomly chosen small radius around the midpoint.

        a_convergence = self._geometry._convergence(label, a_edge, a)
        b_convergence = self._geometry._convergence(label, b_edge, b)

        convergence = min(a_convergence - Δ0.norm(), b_convergence - Δ1.norm())

        # TODO: What should the divisor be here?
        r = convergence / 2

        cost = self._L2_consistency_cost_ball(T0, T1, r, debug)
        self._debugs.append(debug)

        return cost

    @cached_method
    def _L2_consistency(self):
        r"""
        # TODO: This description is not accurate anymore.
        For each pair of adjacent centers along a homology path we use for
        integrating, let `v` be the weighed midpoint (weighed according to the
        radii of convergence at the adjacent centers.)
        We develop the power series coming from both triangles around that
        midpoint and check them for agreement. Namely, we integrate the square
        of their difference on the circle of maximal radius around `v` as a
        line integral. (If the power series agree, that integral should be
        zero.)

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: consistency = PowerSeriesConstraints(T, prec=η.precision(), geometry=Ω._geometry)._L2_consistency()
            sage: η._evaluate(consistency)  # tol 1e-9
            0

        """
        self._debugs = []

        R = self.symbolic_ring(self.real_field())

        cost = R.zero()

        if not self._geometry._singularities:
            # We develop around the end points of each homology generator.

            for (label, edge), a, b in self._geometry._homology_generators:
                cost += self._L2_consistency_between_nonsingular_points(label, edge, a, edge, b)

            # TODO: Replace these hard-coded conditions with something generic.
            # Maybe, take a Delaunay triangulation of the centers in a polygon and
            # then make sure that we have at least a condition on the four shortest
            # edges of each vertex.
            from sage.all import QQ
            cost += self._L2_consistency_between_nonsingular_points(0, 0, QQ(137)/482, 1, QQ(137)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 1, QQ(137)/482, 2, QQ(137)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 2, QQ(137)/482, 3, QQ(137)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 3, QQ(137)/482, 0, QQ(345)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 0, QQ(345)/482, 1, QQ(345)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 1, QQ(345)/482, 2, QQ(345)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 2, QQ(345)/482, 3, QQ(345)/482)
            cost += self._L2_consistency_between_nonsingular_points(0, 3, QQ(345)/482, 0, QQ(137)/482)

            cost += self._L2_consistency_between_nonsingular_points(0, 0, QQ(427)/964, 1, QQ(427)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 1, QQ(427)/964, 2, QQ(427)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 2, QQ(427)/964, 3, QQ(427)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 3, QQ(427)/964, 0, QQ(537)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 0, QQ(537)/964, 1, QQ(537)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 1, QQ(537)/964, 2, QQ(537)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 2, QQ(537)/964, 3, QQ(537)/964)
            cost += self._L2_consistency_between_nonsingular_points(0, 3, QQ(537)/964, 0, QQ(427)/964)
        else:
            V = self._voronoi_diagram()
            centers = self._voronoi_diagram_centers()

            # We integrate L2 errors along the boundary of Voronoi cells.
            segments = [segment for center in centers for boundary in V.cell(center) for segment in boundary.segments_with_uniform_roots()]
            assert all([s for s in segments if s == -segment] for segment in segments)
            for i, segment in enumerate(segments):
                if -segment in segments[:i]:
                    continue
                print(i, "/", len(segments))
                cost += 2 * self._L2_consistency_voronoi_boundary(segment)

        return cost

    @cached_method
    def ζ(self, d):
        from sage.all import exp, pi, I
        return self.complex_field()(exp(2*pi*I / (d + 1)))

    def _range(self, center, prec):
        if center.is_vertex():
            return range(3*prec)
        return range(prec)

    def _gen(self, kind, center, n):
        if center.is_vertex():
            return self.symbolic_ring(self.real_field()).gen((kind, center, n))
        else:
            return self.symbolic_ring(self.real_field()).gen((kind, center, n))

    class Integrator:
        def __init__(self, constraints, segment):
            # TODO: Verify that roots are consistent along the segment! Always the case for the octagon.
            self._segment = segment
            self._constraints = constraints
            self._surface = constraints._surface
            self.real_field = constraints.real_field()
            self._center, self._label, self._center_coordinates = segment.center()
            self._opposite_center, label, self._opposite_center_coordinates = segment.opposite_center()
            self.complex_field = constraints.complex_field()
            self.I = self.complex_field.gen()
            self.R = constraints.symbolic_ring(self.real_field)
            self.C = constraints.symbolic_ring(self.complex_field)

        def integral(self, α, κ, d, n):
            r"""
            Return

            \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1}
            """
            _, segment = self._segment.segment()
            a, b = segment.endpoints()
            a = self.complex_field(*a)
            b = self.complex_field(*b)

            # Since γ(t) = (1 - t)a + tb, we have ·γ(t) = b - a
            constant = self.ζ(d)**(κ * (n + 1)) / (d + 1) * (b - a)

            def value(part, t):
                z = self.complex_field(*((1 - t) * a + t * b))

                value = constant * (z-α).nth_root(d + 1)**(n - d)

                if part == "Re":
                    return float(value.real())
                if part == "Im":
                    return float(value.imag())

                raise NotImplementedError

            from scipy.integrate import quad
            real, error = quad(lambda t: value("Re", t), 0, 1)
            imag, error = quad(lambda t: value("Im", t), 0, 1)

            return self.complex_field(real, imag)

        def integral2(self, part, α, κ, d, n, β, λ, dd, m):
            r"""
            Return the real/imaginary part of

            \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} \overline{ζ_{dd+1}^{λ (m+1)}/(dd+1) (z-β)^\frac{m-dd}{dd+1}} dz
            """
            C = self.complex_field
            _, segment = self._segment.segment()
            a, b = segment.endpoints()
            # 20s
            a = Cab(C, a)
            # 20s
            b = Cab(C, b)

            # 100s
            return integral2(part, α, κ, d, self.ζ(d), n, β, λ, dd, self.ζ(dd), m, a=a, b=b, C=C, R=self.real_field)

        def ζ(self, d):
            return self._constraints.ζ(d)

        def α(self):
            return self.complex_field(*self._center_coordinates) - self.midpoint()

        def β(self):
            return self.complex_field(*self._opposite_center_coordinates) - self.midpoint()

        def midpoint(self):
            # TODO: Does it matter where we chose the midpoint? Should it be the weighed midpoint?
            return (self.complex_field(*self._center_coordinates) + self.complex_field(*self._opposite_center_coordinates)) / 2

        def κ(self):
            if not self._center.is_vertex():
                return 0

            return self._κλ(self._center_coordinates, self._segment.segment()[1])

        def d(self):
            if not self._center.is_vertex():
                return 0

            return 2

        def λ(self):
            if not self._opposite_center.is_vertex():
                return 0

            return self._κλ(self._opposite_center_coordinates, self._segment.segment()[1])

        def dd(self):
            if not self._opposite_center.is_vertex():
                return 0

            return 2

        def _κλ(self, center, segment):
            # TODO: Mostly duplicated as "branch" in GeometricPrimitives.
            low = min([endpoint[1] for endpoint in segment.endpoints()]) < center[1]

            for i, vertex in enumerate(self._surface.polygon(self._label).vertices()):
                if vertex == center:
                    if i == 0:
                        return 0
                    if i == 1:
                        return 1
                    if i == 2:
                        return 0 if low else 2
                    if i == 3:
                        return 1 if low else 0
                    if i == 4:
                        return 2
                    if i == 5:
                        return 0
                    if i == 6:
                        return 1
                    if i == 7:
                        return 2

            raise NotImplementedError

        def range(self, prec):
            return self._constraints._range(self._center, prec)

        def opposite_range(self, prec):
            return self._constraints._range(self._opposite_center, prec)

        def a(self, n):
            return self.Re_a(n) + self.I * self.C(self.Im_a(n))

        def Re_a(self, n):
            return self._constraints._gen("Re", self._center, n)

        def Im_a(self, n):
            return self._constraints._gen("Im", self._center, n)

        def Re_b(self, n):
            return self._constraints._gen("Re", self._opposite_center, n)

        def Im_b(self, n):
            return self._constraints._gen("Im", self._opposite_center, n)

        def f(self, n):
            return self.integral(self.α(), self.κ(), self.d(), n)

    def _voronoi_diagram_centers(self):
        # TODO: Hardcoded octagon here.
        S = self._surface
        center = S(0, S.polygon(0).centroid())
        centers = S.vertices().union([center])

        return centers

    def _voronoi_diagram(self):
        # TODO: Hardcoded octagon here.
        S = self._surface

        def weight(center):
            if center == S.polygon(0).centroid():
                from sage.all import QQ
                return QQ(center.norm().n())
            if center in S.polygon(0).vertices():
                from sage.all import QQ
                return QQ(1)
            raise NotImplementedError

        from flatsurf.geometry.voronoi import FixedVoronoiWeight, VoronoiDiagram
        return VoronoiDiagram(S, self._voronoi_diagram_centers(), weight=FixedVoronoiWeight(weight))

    def _L2_consistency_voronoi_boundary(self, boundary_segment):
        r"""
        ALGORITHM:

        Two cells meet at the ``boundary_segment``. The harmonic differential
        can on these cells abstractly be described as a power series around
        the center of the cell

        g(y) = Σ_{n ≥ 0} a_n y^n

        where ``d`` is the order of the singularity at the center; this is,
        however, not the representation on any of the charts in which the
        translation surface is given.

        To describe this power series on such a chart, let `z` denote the
        variable on that chart, we are going to have to takes `d+1`-st roots
        of `y`. Note that ``boundary_segment`` is assumed to be such that this
        is consistently possible, namely ``boundary_segment`` does not cross
        the horizontal line on which the singularity lives in the `z`-chart.

        Therefore, we write with the center y=0 being z=α

        y(z) = ζ_{d+1}^κ (z-α)^{1/(d+1)}

        where the last part denotes the principal `d+1`-st root of `z`.

        Hence, we can rewrite `g(y)dy` on the `z`-chart:

        g(y)dy = g(y(z)) dy/dz dz
               = Σ_{n ≥ 0} a_n ζ_{d+1}^{κ n} (z-α)^{n/(d+1)} ζ_{d+1}^κ 1/(d+1) (z-α)^{-d/(d+1)}dz
               = Σ_{n ≥ 0} a_n ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} dz
               =: Σ_{n ≥ 0} a_n f_n(z) dz
               =: f(z) dz

        Note that the formulas above hold when the center is not an actual
        singularity, i.e., d = 0.

        Now, we want to describe the error between two such series when
        integrating along the ``boundary_segment``, namely, for two such
        differentials `f(z)dz` and `g(z)dz` we compute the L2 norm of `f-g` or
        rather the square thereof `\int_γ (f-g)\overline{(f-g)} dz` where γ is
        the ``boundary_segment``.

        TODO: This is not actually the integral but a norm, we get |·γ| in the
        integration and not just the signed derivative.

        As discussed above, we have

        f = Σ_{n ≥ 0} a_n f_n(z),

        g = Σ_{m ≥ 0} b_m g_m(z).

        The `a_n` and `b_n` are symbolic variables since we are going to
        optimize for these later. If we evaluate that L2 norm squared we get a
        polynomial of homogeneous degree two in these variables.

        Namely,

        (f-g)\overline{(f-g)} = f\overline{f} - 2 \Re f\overline{g} + g\overline{g}.

        Note that all terms are real.

        The first term is going to yield polynomials in `a_n a_m`, the second
        term mixed polynomials `a_n b_m`, and the last term polynomials in `b_n
        b_m`.

        In fact, our symbolic variables are going to be `\Re a_n`, `\Im a_n`,
        `\Re b_m`, `\Im b_m`.

        Working through the first of the three components of the L2 norm
        expression, we have

        \int_γ f\overline{f}
        = Σ_{n,m ≥ 0} a_n \overline{a_m} \int_γ f_n(z)\overline{f_m(z)}
        =: Σ_{n,m ≥ 0} a_n \overline{a_m} (f_{\Re, n, m} + i f_{\Im, n, m})

        We can simplify things a bit since we know that the result is real;
        inside the series we have therefore

        \Re (a_n \overline{a_m} (f_{\Re, n, m} + i f_{\Im, n, m}))
        = (\Re a_n \Re a_m + \Im a_n \Im a_m - i \Re a_n \Im a_m + i \Im a_n \Re a_m) (f_{\Re, n, m} + i f_{\Im, n, m})
        = \Re a_n \Re a_m f_{\Re, n, m} + \Im a_n \Im a_m f_{\Re, n, m} + \Re a_n \Im a_m f_{\Im, n, m} - \Im a_n \Re a_m f_{\Im, n, m}.

        We get essentially the same terms for `g\overline{g}`.

        Finally, we need to compute the middle term:

        - \int_γ 2 \Re f\overline{g}
        = - 2 \Re \int_γ f\overline{g}
        = - 2 \Re Σ_{n,m ≥ 0} a_n \overline{b_m} \int_γ f_n{z)\overline{g_m(z)}
        =: - 2 \Re Σ_{n,m ≥ 0} a_n \overline{b_m} ((f,g)_{\Re, n, m} + i (f,g)_{\Im, n, m})

        Again we can simplify and get inside the series

        \Re (a_n \overline{b_m} ((f,g)_{\Re, n, m} + i (f,g)_{\Im, n, m}))
        = \Re a_n \Re b_m (f,g)_{\Re, n, m} + \Im a_n \Im b_m (f,g)_{\Re, n, m} + \Re a_n \Im b_m (f,g)_{\Im, n, m} - \Im a_n \Re b_m (f,g)_{\Im, n, m}

        The integrals `f_{\Re, n, m}`, `f_{\Im, n, m}`, `(f,g)_{\Re, n, m}`,
        and `(f,g)_{\Im, n, m}` are currently all computed numerically. We know
        of but have not implemented a better approach, yet.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True)
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=3, geometry=Ω._geometry)
            sage: V = C._voronoi_diagram()
            sage: centers = C._voronoi_diagram_centers()
            sage: segments = [segment for center in centers for boundary in V.cell(center) for segment in boundary.segments_with_uniform_roots()]

            sage: segment = segments[0]
            sage: segment._center, segment._other
            (Vertex 0 of polygon 0, Vertex 0 of polygon 0)

            sage: E = C._L2_consistency_voronoi_boundary(segment)
            sage: F = C._L2_consistency_voronoi_boundary(-segment)
            sage: (E - F).map_coefficients(lambda c: c if abs(c) > 1e-15 else 0)
            0.000000000000000

            sage: segment = segments[-1]
            sage: segment._center, segment._other
            (Point (1/2, 1/2*a + 1/2) of polygon 0, Vertex 0 of polygon 0)

            sage: E = C._L2_consistency_voronoi_boundary(segment)
            sage: F = C._L2_consistency_voronoi_boundary(-segment)
            sage: (E - F).map_coefficients(lambda c: c if abs(c) > 1e-9 else 0)
            0.000000000000000

        """
        integrator = self.Integrator_(self.Integrator(self, boundary_segment), self._prec)

        ff = integrator.int_f_overline_f()
        gg = integrator.int_g_overline_g()
        fg = integrator.int_Re_f_overline_g()

        return ff - 2*fg + gg

    class Integrator_:
        def __init__(self, integrator, prec):
            self._integrator = integrator
            self._prec = prec

            self._center, self._label, self._center_coordinates = self._integrator._segment.center()
            self._opposite_center, self._label, self._opposite_center_coordinates = self._integrator._segment.opposite_center()

            self._Re = "Re"
            self._Im = "Im"

            self._integrator = integrator

            self._R = integrator.R

            self._Re_a = integrator.Re_a
            self._Im_a = integrator.Im_a

            self._Re_b = integrator.Re_b
            self._Im_b = integrator.Im_b

        def int_f_overline_f(self):
            r"""
            Return the value of `\int_γ f\overline{f}.
            """
            value = self._R.zero()
            for n in self._integrator.range(self._prec):
                for m in self._integrator.range(self._prec):
                    value += (self._Re_a(n) * self._Re_a(m) + self._Im_a(n) * self._Im_a(m)) * self.f_(self._Re, n, m) + (self._Re_a(n) * self._Im_a(m) - self._Im_a(n) * self._Re_a(m)) * self.f_(self._Im, n, m)

            return value

        def int_g_overline_g(self):
            r"""
            Return the value of `\int_γ g\overline{g}.
            """
            value = self._R.zero()
            for n in self._integrator.opposite_range(self._prec):
                for m in self._integrator.opposite_range(self._prec):
                    value += (self._Re_b(n) * self._Re_b(m) + self._Im_b(n) * self._Im_b(m)) * self.g_(self._Re, n, m) + (self._Re_b(n) * self._Im_b(m) - self._Im_b(n) * self._Re_b(m)) * self.g_(self._Im, n, m)

            return value

        def int_Re_f_overline_g(self):
            r"""
            Return the value of `\int_γ f\overline{g}.
            """
            value = self._R.zero()
            for n in self._integrator.range(self._prec):
                for m in self._integrator.opposite_range(self._prec):
                    value += (self._Re_a(n) * self._Re_b(m) + self._Im_a(n) * self._Im_b(m)) * self.fg_(self._Re, n, m) + (self._Re_a(n) * self._Im_b(m) - self._Im_a(n) * self._Re_b(m)) * self.fg_(self._Im, n, m)

            return value

        def f_(self, part, n, m):
            return self._integrator.integral2(part, self._integrator.α(), self._integrator.κ(), self._integrator.d(), n, self._integrator.α(), self._integrator.κ(), self._integrator.d(), m)

        def g_(self, part, n, m):
            return self._integrator.integral2(part, self._integrator.β(), self._integrator.λ(), self._integrator.dd(), n, self._integrator.β(), self._integrator.λ(), self._integrator.dd(), m)

        def fg_(self, part, n, m):
            return self._integrator.integral2(part, self._integrator.α(), self._integrator.κ(), self._integrator.d(), n, self._integrator.β(), self._integrator.λ(), self._integrator.dd(), m)

    @cached_method
    def _elementary_line_integrals(self, label, n, m):
        r"""
        Return the integrals f(z)dx and f(z)dy where f(z) = z^n\overline{z}^m
        along the boundary of the polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 3, geometry=Ω._geometry)

            sage: C._elementary_line_integrals(0, 0, 0)  # tol 1e-9
            (0.0000000000000000, 0.0000000000000000)
            sage: C._elementary_line_integrals(0, 1, 0)  # tol 1e-9
            (0 - 1.0*I, 1.0 + 0.0*I)
            sage: C._elementary_line_integrals(0, 0, 1)  # tol 1e-9
            (0.0 + 1.0*I, 1.0 - 0.0*I)
            sage: C._elementary_line_integrals(0, 1, 1)  # tol 1e-9
            (0, 0)

        """
        ix = self.complex_field().zero()
        iy = self.complex_field().zero()

        polygon = self._surface.polygon(label)
        center = polygon.circumscribing_circle().center()

        for v, e in zip(polygon.vertices(), polygon.edges()):
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
    def _elementary_area_integral(self, label, n, m):
        r"""
        Return the integral of z^n\overline{z}^m on the polygon with ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 3, geometry=Ω._geometry)
            sage: C._elementary_area_integral(0, 0, 0)  # tol 1e-9
            1.0 + 0.0*I

            sage: C._elementary_area_integral(0, 1, 0)  # tol 1e-6
            0.0

        """
        C = self.complex_field()
        # Write f(n, m) for z^n\overline{z}^m.
        # Then 1/(2m + 1) [d/dx f(n, m+1) - d/dy -i f(n, m+1)] = f(n, m).

        # So we can use Green's theorem to compute this integral by integrating
        # on the boundary of the triangle:
        # -i/(2m + 1) f(n, m + 1) dx + 1/(2m + 1) f(n, m + 1) dy

        ix, iy = self._elementary_line_integrals(label, n, m+1)

        from sage.all import I
        return -I/(C(2)*(m + 1)) * ix + C(1)/(2*(m + 1)) * iy

    # def _area(self):
    #     r"""
    #     Return a formula for the area ∫ η \wedge \overline{η}.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()

    #         sage: H = SimplicialCohomology(T)
    #         sage: a, b = H.homology().gens()
    #         sage: f = H({a: 1})

    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: η = Ω(f)

    #         sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: area = PowerSeriesConstraints(T, η.precision(), geometry=Ω._geometry)._area()

    #         sage: η._evaluate(area)  # tol 1e-9
    #         1.0 + 0.0*I

    #     """
    #     R = self.symbolic_ring()

    #     area = R.zero()

    #     # We eveluate the integral of |f|^2 on each triangle where f is the
    #     # power series on the Voronoy cell containing that triangle.
    #     for triangle in self._surface.labels():
    #         # We expand the integrand Σa_nz^n · \overline{Σa_mz^m} naively as
    #         # the sum of a_n \overline{a_m} z^n \overline{z^m}.
    #         for n in range(self._prec):
    #             for m in range(self._prec):
    #                 coefficient = self.gen(triangle, n) * self.gen(triangle, m, conjugate=True)
    #                 # Now we have to integrate z^n \overline{z^m} on the triangle.
    #                 area += coefficient * self._elementary_area_integral(triangle, n, m)

    #     return area

    def _squares(self):
        raise NotImplementedError
        # cost = self.symbolic_ring().zero()
        # for (label, edge, pos) in self.symbolic_ring()._gens:
        #     for k in range(1, self._prec):
        #         cost += self.real(label, edge, pos, k)**2
        #         cost += self.imag(label, edge, pos, k)**2
        # return cost

    # def _area_upper_bound(self):
    #     r"""
    #     Return an upper bound for the area 1/π ∫ η \wedge \overline{η}.

    #     EXAMPLES::

    #         sage: from flatsurf import translation_surfaces, SimplicialCohomology, HarmonicDifferentials
    #         sage: T = translation_surfaces.torus((1, 0), (0, 1))
    #         sage: T.set_immutable()

    #         sage: H = SimplicialCohomology(T)
    #         sage: a, b = H.homology().gens()
    #         sage: f = H({a: 1})

    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: η = Ω(f)

    #         sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
    #         sage: Ω = HarmonicDifferentials(T)
    #         sage: area = PowerSeriesConstraints(T, η.precision(), geometry=Ω._geometry)._area_upper_bound()

    #     The correct area would be 1/π here. However, we approximate the square
    #     with a circle of radius 1/sqrt(2) for a factor π/2::

    #         sage: η._evaluate(area)  # tol 1e-9
    #         0.5

    #     """
    #     # Our upper bound is integrating the quantity over the circumcircles of
    #     # the Delaunay triangles instead, i.e.,
    #     # we integrate the sum of the a_n \overline{a_m} z^n \overline{z}^m.
    #     # Integration in polar coordinates, namely
    #     # ∫ z^n\overline{z}^m = ∫_0^R ∫_0^{2π} ir e^{it(n-m)}
    #     # shows that only when n = m the value does not vanish and is
    #     # π/(n+1) |a_n|^2 R^(2n+2).
    #     area = self.symbolic_ring().zero()

    #     # TODO: This does not match the above description at all anymore.
    #     for (label, edge, pos) in self.symbolic_ring()._gens:
    #         R2 = self.real_field()(self._surface.polygon(label).circumscribing_circle().radius_squared())

    #         for k in range(self._prec):
    #             area += (self.real(label, edge, pos, k)**2 + self.imag(label, edge, pos, k)**2) * R2**(k + 1)

    #     return area

    def optimize(self, f):
        r"""
        Add constraints that optimize the symbolic expression ``f``.

        EXAMPLES:

        TODO: All these examples are a bit pointless::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 1, geometry=Ω._geometry)
            sage: R = C.symbolic_ring()

        We optimize a function in two variables. Since there are no
        constraints, we do not get any Lagrange multipliers in this
        optimization but just for roots of the derivative::

            sage: f = 10*C._real_nonsingular(0, 0, 0, 0)^2 + 16*C._imag_nonsingular(0, 0, 0, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [20.0000000000000*Re(a2,0), 32.0000000000000*Im(a2,0)]

        In this example, the constraints and the optimized values do not
        overlap, so again we do not get Lagrange multipliers::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), Re(a0,0) + 0.160570808419475*Im(a0,1) - Re(a2,0) + 0.123661556725753*Im(a2,1), Im(a0,0) - 0.160570808419475*Re(a0,1) - Im(a2,0) - 0.123661556725753*Re(a2,1)]

            sage: f = 10*C._real_nonsingular(0, 0, 0, 0)^2 + 17*C._imag_nonsingular(0, 0, 0, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), Re(a0,0) + 0.160570808419475*Im(a0,1) - Re(a2,0) + 0.123661556725753*Im(a2,1), Im(a0,0) - 0.160570808419475*Re(a0,1) - Im(a2,0) - 0.123661556725753*Re(a2,1), -λ16 + λ18, -λ17 + λ19, -λ6 + λ8, -λ7 + λ9, 20.0000000000000*Re(a2,0) + λ0 - λ8 + λ10 - λ18, 34.0000000000000*Im(a2,0) + λ1 - λ9 + λ11 - λ19, -λ14 + λ16, -λ15 + λ17, -λ4 + λ6, -λ5 + λ7, -λ12 + λ14, -λ13 + λ15, -λ0 + λ2, -λ1 + λ3, -λ2 + λ4, -λ3 + λ5, -λ10 + λ12, -λ11 + λ13, -0.0762270995598000*λ17 - 0.160570808419475*λ19, 0.0762270995598000*λ16 + 0.160570808419475*λ18, 0.0762270995598000*λ6 + 0.160570808419475*λ8, 0.0762270995598000*λ7 + 0.160570808419475*λ9, 0.123661556725753*λ0 + 0.123661556725753*λ8 - 0.123661556725753*λ11 - 0.123661556725753*λ19, 0.123661556725753*λ1 + 0.123661556725753*λ9 + 0.123661556725753*λ10 + 0.123661556725753*λ18, -0.0570539419087137*λ15 - 0.0824865933862581*λ17, 0.0570539419087137*λ14 + 0.0824865933862581*λ16, 0.0570539419087137*λ4 + 0.0824865933862581*λ6, 0.0570539419087137*λ5 + 0.0824865933862581*λ7, -0.0824865933862581*λ13 - 0.0570539419087137*λ15, 0.0824865933862581*λ12 + 0.0570539419087137*λ14, 0.160570808419475*λ0 + 0.0762270995598000*λ2, 0.160570808419475*λ1 + 0.0762270995598000*λ3, 0.0824865933862581*λ2 + 0.0570539419087137*λ4, 0.0824865933862581*λ3 + 0.0570539419087137*λ5, -0.160570808419475*λ11 - 0.0762270995598000*λ13, 0.160570808419475*λ10 + 0.0762270995598000*λ12]

        ::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), Re(a0,0) + 0.160570808419475*Im(a0,1) - Re(a2,0) + 0.123661556725753*Im(a2,1), Im(a0,0) - 0.160570808419475*Re(a0,1) - Im(a2,0) - 0.123661556725753*Re(a2,1)]

            sage: f = 3*C._real_nonsingular(0, 0, 0, 1)^2 + 5*C._imag_nonsingular(0, 0, 0, 1)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), -Re(a4,0) + 0.0570539419087137*Re(a4,1) + Re(a7,0) + 0.0570539419087137*Re(a7,1), -Im(a4,0) + 0.0570539419087137*Im(a4,1) + Im(a7,0) + 0.0570539419087137*Im(a7,1), -Re(a1,0) + 0.0762270995598000*Re(a1,1) + Re(a4,0) + 0.0824865933862581*Re(a4,1), -Im(a1,0) + 0.0762270995598000*Im(a1,1) + Im(a4,0) + 0.0824865933862581*Im(a4,1), Re(a1,0) + 0.160570808419475*Re(a1,1) - Re(a2,0) + 0.123661556725753*Re(a2,1), Im(a1,0) + 0.160570808419475*Im(a1,1) - Im(a2,0) + 0.123661556725753*Im(a2,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), -Re(a5,0) + 0.0824865933862581*Im(a5,1) + Re(a8,0) + 0.0762270995598000*Im(a8,1), -Im(a5,0) - 0.0824865933862581*Re(a5,1) + Im(a8,0) - 0.0762270995598000*Re(a8,1), -Re(a3,0) + 0.0570539419087137*Im(a3,1) + Re(a5,0) + 0.0570539419087137*Im(a5,1), -Im(a3,0) - 0.0570539419087137*Re(a3,1) + Im(a5,0) - 0.0570539419087137*Re(a5,1), -Re(a0,0) + 0.0762270995598000*Im(a0,1) + Re(a3,0) + 0.0824865933862581*Im(a3,1), -Im(a0,0) - 0.0762270995598000*Re(a0,1) + Im(a3,0) - 0.0824865933862581*Re(a3,1), Re(a0,0) + 0.160570808419475*Im(a0,1) - Re(a2,0) + 0.123661556725753*Im(a2,1), Im(a0,0) - 0.160570808419475*Re(a0,1) - Im(a2,0) - 0.123661556725753*Re(a2,1), -λ16 + λ18, -λ17 + λ19, -λ6 + λ8, -λ7 + λ9, λ0 - λ8 + λ10 - λ18, λ1 - λ9 + λ11 - λ19, -λ14 + λ16, -λ15 + λ17, -λ4 + λ6, -λ5 + λ7, -λ12 + λ14, -λ13 + λ15, -λ0 + λ2, -λ1 + λ3, -λ2 + λ4, -λ3 + λ5, -λ10 + λ12, -λ11 + λ13, -0.0762270995598000*λ17 - 0.160570808419475*λ19, 0.0762270995598000*λ16 + 0.160570808419475*λ18, 0.0762270995598000*λ6 + 0.160570808419475*λ8, 0.0762270995598000*λ7 + 0.160570808419475*λ9, 6.00000000000000*Re(a2,1) + 0.123661556725753*λ0 + 0.123661556725753*λ8 - 0.123661556725753*λ11 - 0.123661556725753*λ19, 10.0000000000000*Im(a2,1) + 0.123661556725753*λ1 + 0.123661556725753*λ9 + 0.123661556725753*λ10 + 0.123661556725753*λ18, -0.0570539419087137*λ15 - 0.0824865933862581*λ17, 0.0570539419087137*λ14 + 0.0824865933862581*λ16, 0.0570539419087137*λ4 + 0.0824865933862581*λ6, 0.0570539419087137*λ5 + 0.0824865933862581*λ7, -0.0824865933862581*λ13 - 0.0570539419087137*λ15, 0.0824865933862581*λ12 + 0.0570539419087137*λ14, 0.160570808419475*λ0 + 0.0762270995598000*λ2, 0.160570808419475*λ1 + 0.0762270995598000*λ3, 0.0824865933862581*λ2 + 0.0570539419087137*λ4, 0.0824865933862581*λ3 + 0.0570539419087137*λ5, -0.160570808419475*λ11 - 0.0762270995598000*λ13, 0.160570808419475*λ10 + 0.0762270995598000*λ12]

        """
        if f:
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

        # TODO: We add lots of lagrange coefficients even if our rank is very
        # small. We should probably determine the rank here and reduce the
        # system first.
        lagranges = len(self._constraints)

        g = self._constraints

        # We form the partial derivative with respect to the variables Re(a_k)
        # and Im(a_k).
        for gen in self.symbolic_ring()._regular_gens(prec=self._prec):
            gen = self.symbolic_ring().gen(gen)
            if self._cost.degree(gen) <= 0:
                # The cost function does not depend on this variable.
                # That's fine, we still need it for the Lagrange multipliers machinery.
                pass

            gen = self._cost.parent()(gen)

            L = self._cost.derivative(gen)

            for i in range(lagranges):
                L += g[i][gen] * self.lagrange(i)

            self.add_constraint(L, rank_check=False)

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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()

        Integrating along the (negative horizontal) cycle `b`, produces
        something with `-Re(a_0)` in the real part.
        Integration along the diagonal `a`, produces essentially `Re(a_0) -
        Im(a_0)`. Note that the two variables ``a0`` and ``a1`` are the same
        because the centers of the Voronoi cells for the two triangles are
        identical::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 1, geometry=Ω._geometry)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.236797907979275*Re(a1,0) + 0.247323113451507*Re(a2,0) + 0.139540535294972*Re(a4,0) + 0.236797907979275*Re(a6,0) + 0.139540535294972*Re(a7,0), 0.236797907979275*Im(a0,0) + 0.247323113451507*Im(a2,0) + 0.139540535294972*Im(a3,0) + 0.139540535294972*Im(a5,0) + 0.236797907979275*Im(a8,0) - 1.00000000000000]

        If we increase precision, we see additional higher imaginary parts.
        These depend on the choice of base point of the integration and will be
        found to be zero by other constraints, not true anymore TODO::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.236797907979275*Re(a1,0) + 0.00998620690459202*Re(a1,1) + 0.247323113451507*Re(a2,0) + 0.139540535294972*Re(a4,0) + 0.00177444290057350*Re(a4,1) + 0.236797907979275*Re(a6,0) - 0.00998620690459202*Re(a6,1) + 0.139540535294972*Re(a7,0) - 0.00177444290057350*Re(a7,1), 0.236797907979275*Im(a0,0) - 0.00998620690459202*Re(a0,1) + 0.247323113451507*Im(a2,0) + 0.139540535294972*Im(a3,0) - 0.00177444290057350*Re(a3,1) + 0.139540535294972*Im(a5,0) + 0.00177444290057350*Re(a5,1) + 0.236797907979275*Im(a8,0) + 0.00998620690459202*Re(a8,1) - 1.00000000000000]

        """
        for cycle in cocycle.parent().homology().gens():
            self.add_constraint(self.real_part(self.integrate(cycle)) - self.real_part(cocycle(cycle)), rank_check=False)

    def lagrange_variables(self):
        return set(variable for constraint in self._constraints for variable in constraint.variables("lagrange"))

    def matrix(self, nowarn=False):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialCohomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialCohomology(T)
            sage: a, b = H.homology().gens()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 6, geometry=Ω._geometry)
            sage: C.require_cohomology(H({a: 1}))
            sage: C.optimize(C._L2_consistency())
            sage: C.matrix()
            ...
            (2 x 54 dense matrix over Real Field with 54 bits of precision,
             (1.00000000000000, 0.000000000000000),
             {1: 0,
              2: 1,
              4: 2,
              5: 3,
              7: 4,
              8: 5,
              11: 6,
              12: 7,
              14: 8,
              17: 9,
              18: 10,
              20: 11,
              24: 12,
              26: 13,
              28: 14,
              30: 15,
              32: 16,
              34: 17,
              37: 18,
              38: 19,
              40: 20,
              41: 21,
              43: 22,
              44: 23,
              47: 24,
              48: 25,
              50: 26,
              53: 27,
              54: 28,
              56: 29,
              60: 30,
              62: 31,
              64: 32,
              66: 33,
              68: 34,
              70: 35,
              73: 36,
              74: 37,
              76: 38,
              77: 39,
              79: 40,
              80: 41,
              83: 42,
              84: 43,
              86: 44,
              89: 45,
              90: 46,
              92: 47,
              96: 48,
              98: 49,
              100: 50,
              102: 51,
              104: 52,
              106: 53},
             set())

        ::

            sage: R = C.symbolic_ring()
            sage: f = 10*C._real_nonsingular(0, 0, 0, 0)^2 + 16*C._imag_nonsingular(0, 0, 0, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C.matrix()  # not tested TODO

        """
        lagranges = {}

        non_lagranges = set()
        for row, constraint in enumerate(self._constraints):
            for monomial, coefficient in constraint._coefficients.items():
                if not monomial:
                    continue

                gen = monomial[0]
                if gen >= 0:
                    non_lagranges.add(gen)

        non_lagranges = {gen: i for (i, gen) in enumerate(non_lagranges)}

        if len(set(self.symbolic_ring()._regular_gens(self._prec))) != len(non_lagranges):
            if not nowarn:
                from warnings import warn
                warn(f"Some power series coefficients are not constrained for this harmonic differential. They will be chosen to be 0 by the solver.")

        from sage.all import matrix, vector
        A = matrix(self.real_field(), len(self._constraints), len(non_lagranges) + len(self.lagrange_variables()))
        b = vector(self.real_field(), len(self._constraints))

        for row, constraint in enumerate(self._constraints):
            for monomial, coefficient in constraint._coefficients.items():
                if not monomial:
                    b[row] -= coefficient
                    continue

                assert len(monomial) == 1

                gen = monomial[0]
                if gen < 0:
                    if gen not in lagranges:
                        lagranges[gen] = len(non_lagranges) + len(lagranges)
                    column = lagranges[gen]
                else:
                    column = non_lagranges[gen]

                A[row, column] += coefficient

        return A, b, non_lagranges, set()

    @cached_method
    def power_series_ring(self, *args):
        r"""
        Return the power series ring to write down the series describing a
        harmonic differential in a Voronoi cell.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: Ω = PowerSeriesConstraints(T, 8, geometry=Ω._geometry)
            sage: Ω.power_series_ring(T(0, (1/2, 1/2)))
            Power Series Ring in z2 over Complex Field with 54 bits of precision

        """
        from sage.all import PowerSeriesRing

        if len(args) != 1:
            raise NotImplementedError

        point = args[0]
        return PowerSeriesRing(self.complex_field(), f"z{self.symbolic_ring()._centers.index(point)}")

    def laurent_series_ring(self, *args):
        return self.power_series_ring(*args).fraction_field()

    def solve(self, algorithm="eigen+mpfr"):
        r"""
        Return a solution for the system of constraints with minimal error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.add_constraint(C._real_nonsingular(0, 0, 0, 0) - C._real_nonsingular(0, 0, 0, 1))
            sage: C.add_constraint(C._real_nonsingular(0, 0, 0, 0) - 1)
            sage: C.solve()  # random output due to random ordering of the dict
            ({(0, 0, 345/482): O(z0_0_0^2),
              (0, 1, 345/482): O(z0_1_0^2),
              (0, 0, 0): 1.00000000000000 + 1.00000000000000*z0_0_1 + O(z0_0_1^2),
              (0, 0, 537/964): O(z0_0_2^2),
              (0, 1, 537/964): O(z0_1_1^2),
              (0, 0, 427/964): O(z0_0_3^2),
              (0, 1, 137/482): O(z0_1_2^2),
              (0, 1, 427/964): O(z0_1_3^2),
              (0, 0, 137/482): O(z0_0_4^2)},
             0.000000000000000)

        """
        self._optimize_cost()

        A, b, decode, _ = self.matrix()

        rows, columns = A.dimensions()
        # TODO
        # print(f"Solving {rows}×{columns} system")
        rank = A.rank()

        from sage.all import RDF, oo
        condition = A.change_ring(RDF).condition()

        if condition == oo:
            print("condition number is not finite")
        else:
            def cond():
                C = A.list()
                # numpy.random.shuffle() is much faster than random.shuffle() on large lists
                import numpy.random
                numpy.random.shuffle(C)
                return A.parent()(C).change_ring(RDF).condition()

            random_condition = min([cond() for i in range(8)])

            if random_condition / condition < .01:
                print(f"condition is {condition:.1e} but could be {random_condition/condition:.1e} of that")

        if rank < columns:
            # TODO: Warn?
            print(f"system underdetermined: {rows}×{columns} matrix of rank {rank}")

        if algorithm == "arb":
            from sage.all import ComplexBallField
            C = ComplexBallField(self.complex_field().prec())
            CA = A.change_ring(C)
            Cb = b.change_ring(C)

            solution = CA.solve_right(Cb)

            solution = solution.change_ring(self.real_field())
        elif algorithm == "scipy":
            CA = A.change_ring(RDF)
            Cb = b.change_ring(RDF)
            solution = CA.solve_right(Cb)
        elif algorithm == "eigen+mpfr":
            import cppyy

            cppyy.load_library("mpfr")

            solution = define_solve()(A, b)

            from sage.all import vector
            solution = vector([self.real_field()(entry) for entry in solution])
        else:
            raise NotImplementedError

        residue = (A*solution - b).norm()

        lagranges = len(self.lagrange_variables())

        if lagranges:
            solution = solution[:-lagranges]

        series = {point: {} for point in self.symbolic_ring()._centers}

        for gen in self.symbolic_ring()._regular_gens(self._prec):
            i = next(iter(self.symbolic_ring().gen(gen)._coefficients))[0]

            if i not in decode:
                continue

            i = decode[i]
            value = solution[i]

            part, center, k = gen
            if k not in series[center]:
                series[center][k] = [None, None]

            if part == "Re":
                part = 0
            elif part == "Im":
                part = 1
            else:
                raise NotImplementedError

            series[center][k][part] = value

        series = {point:
                  sum((self.complex_field()(*entry) * self.laurent_series_ring(point).gen()**k for (k, entry) in series[point].items()), start=self.laurent_series_ring(point).zero()).add_bigoh(max(series[point]) + 1) for point in series if series[point]}

        return series, residue
