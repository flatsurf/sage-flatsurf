r"""
TODO: Document this module.

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
    sage: Ω(f)
    (1.00000000000000 + O(z0_0_0^5), 1.00000000000000 + O(z0_1_0^5), 1.00000000000000 + O(z0_0_1^5), 1.00000000000000 + O(z0_0_2^5), 1.00000000000000 + O(z0_1_1^5), 1.00000000000000 + O(z0_0_3^5), 1.00000000000000 + O(z0_1_2^5), 1.00000000000000 + O(z0_1_3^5), 1.00000000000000 + O(z0_0_4^5))

The harmonic differential that integrates as 0 along `a` but 1 along `b`::

    sage: g = H({b: 1})
    sage: Ω(g)  # tol 1e-9
    (1.00000000000000*I + O(z0_0_0^5), 1.00000000000000*I + O(z0_1_0^5), 1.00000000000000*I + O(z0_0_1^5), 1.00000000000000*I + O(z0_0_2^5), 1.00000000000000*I + O(z0_1_1^5), 1.00000000000000*I + O(z0_0_3^5), 1.00000000000000*I + O(z0_1_2^5), 1.00000000000000*I + O(z0_1_3^5), 1.00000000000000*I + O(z0_0_4^5))


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
from sage.misc.cachefunc import cached_method
from sage.categories.all import SetsWithPartialMaps
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement

import cppyy
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
            (O(z0_0_0^5), O(z0_1_0^5), O(z0_0_1^5), O(z0_0_2^5), O(z0_1_1^5), O(z0_0_3^5), O(z0_1_2^5), O(z0_1_3^5), O(z0_0_4^5))

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

        if kind is None or "midpoint_derivatives" in kind:
            C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
            for derivative in range(self.precision()//3):
                for ((label, edge), a, b) in self.parent()._geometry._homology_generators:
                    opposite_label, opposite_edge = self.parent().surface().opposite_edge(label, edge)
                    expected = self.evaluate(label, edge, a, C.complex_field()(*self.parent()._geometry.midpoint(label, edge, a, edge, b)), derivative)
                    other = self.evaluate(opposite_label, opposite_edge, 1 - b, C.complex_field()(*self.parent()._geometry.midpoint(opposite_label, opposite_edge, 1 - b, opposite_edge, 1 - a)), derivative)

                    abs_error, rel_error = errors(expected, other)

                    if abs_error > abs_tol or rel_error > rel_tol:
                        report = f"Power series defining harmonic differential are not consistent where triangles meet. {derivative}th derivative does not match between {(label, edge, a)} where it is {expected} and {(label, edge, b)} where it is {other}, i.e., there is an absolute error of {abs_error} and a relative error of {rel_error}."
                        if verbose:
                            print(report)

                        error = report
                        if not verbose:
                            return error

        # if kind is None or "area" in kind:
        #     if verbose:
        #         C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
        #         area = self._evaluate(C._area_upper_bound())

        #         report = f"Area (upper bound) is {area}."
        #         print(report)

        if kind is None or "L2" in kind:
            C = PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)
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
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)
            sage: a, b = H.gens()
            sage: H = SimplicialCohomology(T)
            sage: f = H({a: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: η.series((0, 0, 0))  # abstol 1e-9
            1.00000000000000 + 1.87854000000000e-72*I + (-2.25406000000000e-73 + 1.02693000000000e-72*I)*z0_0_1 + (2.43931000000000e-12 + 2.32526000000000e-71*I)*z0_0_1^2 + (3.89566000000000e-72 - 6.49402000000000e-72*I)*z0_0_1^3 + (-6.75204000000000e-13 - 5.90057000000000e-71*I)*z0_0_1^4 + O(z0_0_1^5)

        """
        return self._series[triangle]

    @cached_method
    def _constraints(self):
        # TODO: This is a hack. Come up with a better class hierarchy!
        return PowerSeriesConstraints(self.parent().surface(), self.precision(), geometry=self.parent()._geometry)

    def evaluate(self, label, edge, pos, Δ, derivative=0):
        C = self._constraints()
        return self._evaluate(C.evaluate(label, edge, pos, Δ, derivative=derivative))

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

    #     for triangle in surface.label_iterator():
    #         series = self._series[triangle]
    #         series = series.truncate(series.prec())
    #         for (root, multiplicity) in series.roots():
    #             if multiplicity != 1:
    #                 raise NotImplementedError

    #             root += root.parent()(*self.parent()._geometry.center(triangle))

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
            sage: η._evaluate(R(C.gen(0, 0, 0, 0))) # abstol 1e-9
            1.00000000000000 + 1.87854000000000e-72*I

        """
        coefficients = {}

        for variable in expression.variables():
            kind, label, edge, pos, k = variable.describe()

            coefficient = self._series[(label, edge, pos)][k]

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
        # TODO: This is the number of coefficients of the power series but we use it as bit precision?
        precisions = set(series.precision_absolute() for series in self._series.values())
        assert len(precisions) == 1
        return next(iter(precisions))

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
    #     for part in ["real", "imag"]:
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
    #             edge_ = (edge - 1) % surface.polygon(label).num_edges()

    #             # TODO print(f"integrating across {triangle} from the midpoint of {edge} to {edge_}")

    #             P = complex_field(*self.parent()._geometry.midpoint(label, edge))
    #             Q = complex_field(*self.parent()._geometry.midpoint(label, edge_))

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

    def integrate(self, cycle):
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

            return series.map_coefficients(compress_coefficient)

        return repr(tuple(compress(series) for series in self._series.values()))

    def plot(self, versus=None):
        from sage.all import RealField, I, vector, complex_plot, oo
        S = self.parent().surface()
        GS = S.graphical_surface()

        plot = S.plot(fill=None)

        for label in S.label_iterator():
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
        Ω(TranslationSurface built from 1 polygon)

    ::

        sage: H = SimplicialCohomology(T)
        sage: Ω(H())
        (O(z0_0_0^5), O(z0_1_0^5), O(z0_0_1^5), O(z0_0_2^5), O(z0_1_1^5), O(z0_0_3^5), O(z0_1_2^5), O(z0_1_3^5), O(z0_0_4^5))

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
    def __classcall__(cls, surface, safety=None, category=None):
        r"""
        Normalize parameters when creating the space of harmonic differentials.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: HarmonicDifferentials(T) is HarmonicDifferentials(T)
            True

        """
        return super().__classcall__(cls, surface, HarmonicDifferentials._homology_generators(surface, safety), category or SetsWithPartialMaps())

    def __init__(self, surface, homology_generators, category):
        Parent.__init__(self, category=category)

        self._surface = surface
        # TODO: Find a better way to track the L2 circles.
        self._debugs = []

        self._geometry = GeometricPrimitives(surface, homology_generators)

    @staticmethod
    def _homology_generators(surface, safety=None):
        # safety = float(safety or 7)

        from flatsurf.geometry.homology import SimplicialHomology
        H = SimplicialHomology(surface, generators="voronoi")

        # TODO: Require that the surface is decomposed into Delaunay cells.

        # The generators of homology will come from paths crossing from a
        # center of a Delaunay cell to the center of a neighbouring Delaunay
        # cell.
        voronoi_paths = set()
        for label, edge in surface.edge_iterator():
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
        else:
            for path in voronoi_paths:
                gens.append((path, QQ(0), QQ(1)))

        gens = tuple(gens)

        return gens

    def plot(self):
        S = self._surface
        G = S.plot()
        SR = PowerSeriesConstraints(S, prec=20, geometry=self._geometry).symbolic_ring()
        for (label, edge, pos) in SR._gens:
            label, center = self._geometry.center(label, edge, pos, wrap=True)
            radius = self._geometry._convergence(label, edge, pos)
            from flatsurf import TranslationSurface
            P = TranslationSurface(S).surface_point(label, center)
            G += P.plot(color="red")

            from sage.all import circle
            G += circle(P.graphical_surface_point().points()[0], radius, color="green", fill="green", alpha=.05)

        for (label, edge, pos) in SR._gens:
            label, center = self._geometry.center(label, edge, pos, wrap=True)

            from flatsurf import TranslationSurface
            P = TranslationSurface(S).surface_point(label, center)
            p = P.graphical_surface_point().points()[0]

            for (lbl, a_edge, a, b_edge, b, Δ0, Δ1, radius) in self._debugs:
                if lbl == label and a_edge == edge and a == pos:
                    Δ = Δ0
                elif self._surface.opposite_edge(lbl, a_edge) == (label, edge) and a == 1 - pos:
                    Δ = Δ0
                elif lbl == label and b_edge == edge and b == pos:
                    Δ = Δ1
                elif self._surface.opposite_edge(lbl, b_edge) == (label, edge) and b == 1 - pos:
                    Δ = Δ1
                else:
                    continue

                from sage.all import vector
                q = p + vector((Δ.real(), Δ.imag()))

                from sage.all import line
                G += line((p, q), color="black")

                from sage.all import circle
                G += circle(q, radius, color="brown", fill="brown", alpha=.1)

        return G

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
            constraints.require_midpoint_derivatives(derivatives)

        # (1') TODO: Describe L2 optimization.
        if "L2" in algorithm:
            weight = get_parameter("L2", 1)
            constraints.optimize(weight * constraints._L2_consistency())
            self._debugs = constraints._debugs

        if "squares" in algorithm:
            weight = get_parameter("squares", 1)
            constraints.optimize(weight * constraints._squares())

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

        if "tykhonov" in algorithm:
            pass

        solution, residue = constraints.solve()
        η = self.element_class(self, solution, residue=residue, cocycle=cocycle)

        if check:
            if report := η.error():
                raise ValueError(report)

        return η


class GeometricPrimitives:
    def __init__(self, surface, homology_generators):
        # TODO: Require immutable.
        self._surface = surface
        self._homology_generators = homology_generators

    @cached_method
    def midpoint(self, label, a_edge, a, b_edge, b):
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

            sage: G.midpoint(0, 0, 0, 0, 1/3)
            (0.000000000000000, -0.142350327708281)
            sage: G.midpoint(0, 2, 2/3, 2, 1)
            (0.000000000000000, 0.190983005625053)
            sage: G.midpoint(0, 2, 2/3, 2, 1) - G.midpoint(0, 0, 0, 0, 1/3)
            (0.000000000000000, 0.333333333333333)

        Here the midpoints are exactly on the edge of the square::

            sage: G.midpoint(0, 0, 1/3, 0, 2/3)
            (0.000000000000000, -0.166666666666667)
            sage: G.midpoint(0, 2, 1/3, 2, 2/3)
            (0.000000000000000, 0.166666666666667)

        Again, the same as in the first example::

            sage: G.midpoint(0, 0, 2/3, 0, 1)
            (0.000000000000000, -0.190983005625053)
            sage: G.midpoint(0, 2, 0, 2, 1/3)
            (0.000000000000000, 0.142350327708281)

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
            sage: T = translation_surfaces.regular_octagon()
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import GeometricPrimitives
            sage: G = GeometricPrimitives(T, None)  # TODO: Should not be None

            sage: G.midpoint(0, 0, 0, 0, 1/2)  # tol 1e-9
            (0.000000000000000, -0.334089318959649)

        A midpoint between two different Voronoi paths::

            sage: G.midpoint(0, 0, 1, 2, 1)
            (1.20710678118655, 1.20710678118655)
            sage: G.midpoint(0, 0, 1, 4, 1)
            (0.000000000000000, 2.41421356237309)

        """
        radii = (
            self._convergence(label, a_edge, a),
            self._convergence(label, b_edge, b),
        )

        centers = (
            self.center(label, a_edge, a)[1],
            self.center(label, b_edge, b)[1],
        )

        return (centers[1] - centers[0]) * radii[1] / (sum(radii))

    @cached_method
    def center(self, label, edge, pos, wrap=False, ring=None):
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

            sage: G.center(0, 0, 0)
            (0, (1/2, 1/2))
            sage: G.center(0, 0, 1/2)
            (0, (1/2, 0))
            sage: G.center(0, 0, 1)
            (0, (1/2, -1/2))
            sage: G.center(0, 0, 1, wrap=True)
            (0, (1/2, 1/2))

        """
        # TODO: This method feels a bit hacky. The name is misleading, the return tuple is weird, and the ring argument is a strange performance hack.
        polygon = self._surface.polygon(label)
        polygon_center = polygon.circumscribing_circle().center()

        opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)
        opposite_polygon = self._surface.polygon(opposite_label)
        opposite_polygon = opposite_polygon.translate(-opposite_polygon.vertex((opposite_edge + 1) % opposite_polygon.num_edges()) + polygon.vertex(edge))

        assert polygon.vertex((edge + 1) % polygon.num_edges()) == opposite_polygon.vertex(opposite_edge)

        opposite_polygon_center = opposite_polygon.circumscribing_circle().center()

        center = (1 - pos) * polygon_center + pos * opposite_polygon_center

        if not wrap or self._surface.polygon(label).contains_point(center):
            if ring is not None:
                from sage.all import vector
                center = vector((ring(center[0]), ring(center[1])))
            return label, center

        label, edge = self._surface.opposite_edge(label, edge)
        return self.center(label, edge, 1 - pos, wrap=True, ring=ring)

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
        label, center = self.center(label, edge, pos)
        polygon = self._surface.polygon(label)

        # TODO: We are assuming that the only relevant singularities are the
        # end points of ``edge``.
        # TODO: Use a ring with more appropriate precision.
        from sage.all import RR
        return min(
            (center - polygon.vertex(edge)).change_ring(RR).norm(),
            (center - polygon.vertex((edge + 1) % polygon.num_edges())).change_ring(RR).norm())


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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))

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

            kind = "Im" if gen % 2 else "Re"
            index = gen % (2 * len(self.parent()._gens)) // 2
            label, edge, pos = self.parent()._gens[index]
            k = gen // (2 * len(self.parent()._gens))

            return kind, label, edge, pos, k

        def variable_name(gen):
            gen = decode(gen)

            if len(gen) == 1:
                return f"λ{gen[0]}"

            kind, label, edge, pos, k = gen
            index = self.parent()._gens.index((label, edge, pos))
            return f"{kind}__open__a{index}__comma__{k}__close__"

        def key(gen):
            gen = decode(gen)

            if len(gen) == 1:
                n = gen[0]
                return 1e9, n

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

        return repr(f).replace('__open__', '(').replace('__close__', ')').replace('__comma__', ',')

    def degree(self, gen):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: x = R.gen(('imag', 0, 0, 0, 0)) + 1; x
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
            sage: a.describe()
            ('imag', 0, 0, 0, 0)
            sage: b.describe()
            ('real', 0, 0, 0, 0)
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

        index = variable % (2 * len(self.parent()._gens)) // 2
        label, edge, pos = self.parent()._gens[index]
        k = variable // (2 * len(self.parent()._gens))
        if variable % 2:
            return ("imag", label, edge, pos, k)
        else:
            return ("real", label, edge, pos, k)

    def real(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: a = R.gen(('imag', 0, 0, 0, 0))
            sage: b = R.gen(('real', 0, 0, 0, 0))
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
    def __classcall__(cls, surface, base_ring, homology_generators, category=None):
        from sage.categories.all import CommutativeRings
        return super().__classcall__(cls, surface, base_ring, homology_generators, category or CommutativeRings())

    def __init__(self, surface, base_ring, homology_generators, category):
        r"""
        TESTS::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import SymbolicCoefficientRing, HarmonicDifferentials
            sage: R = SymbolicCoefficientRing(T, CC, HarmonicDifferentials(T)._geometry._homology_generators)
            sage: R.has_coerce_map_from(CC)
            True

            sage: TestSuite(R).run()

        """
        self._surface = surface
        self._base_ring = base_ring
        self._homology_generators = homology_generators

        self._gens = set()
        for (label, edge), a, b in homology_generators:
            self._gens.add((label, edge if a else 0, a))
            if b == 1:
                label, edge = self._surface.opposite_edge(label, edge)
                edge = 0
                b = 0
            self._gens.add((label, edge, b))
        self._gens = list(self._gens)

        CommutativeRing.__init__(self, base_ring, category=category, normalize=False)
        self.register_coercion(base_ring)

    Element = SymbolicCoefficientExpression

    def _repr_(self):
        return f"Ring of Power Series Coefficients over {self.base_ring()}"

    def change_ring(self, ring):
        return SymbolicCoefficientRing(self._surface, ring, self._homology_generators, category=self.category())

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
            if len(n) == 5:
                kind, polygon, edge, pos, k = n

                if (polygon, edge, pos) not in self._gens:
                    pos = 1 - pos
                    polygon, edge = self._surface.opposite_edge(polygon, edge)

                if pos == 1:
                    label, edge = self._surface.opposite_edge(polygon, edge)
                    pos = 0

                if pos == 0:
                    edge = 0

                if (polygon, edge, pos) not in self._gens:
                    raise ValueError(f"{(polygon, edge, pos)} is not a valid generator; valid generators are {self._gens}")

                if kind == "real":
                    kind = 0
                elif kind == "imag":
                    kind = 1
                else:
                    raise NotImplementedError

                polygon = list(self._surface.label_iterator()).index(polygon)
                assert polygon != -1
                n = k * 2 * len(self._gens) + 2 * self._gens.index((polygon, edge, pos)) + kind
            elif len(n) == 2:
                kind, k = n

                if kind != "lagrange":
                    raise ValueError(kind)
                if k < 0:
                    raise ValueError(str(n))

                n = -k-1
            else:
                raise ValueError(str(n))

        from sage.all import parent, ZZ
        if parent(n) is ZZ:
            n = int(n)

        if not isinstance(n, int):
            raise ValueError(str(n))

        return self((n,))

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
        return SymbolicCoefficientRing(self._surface, base_ring=base_ring or ComplexField(54), homology_generators=self._geometry._homology_generators)

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

    @cached_method
    def gen(self, label, edge, pos, k, /, conjugate=False):
        assert conjugate is True or conjugate is False
        real = self.real(label, edge, pos, k)
        imag = self.imag(label, edge, pos, k)

        i = self.symbolic_ring().imaginary_unit()

        if conjugate:
            i = -i

        return real + i*imag

    @cached_method
    def real(self, label, edge, pos, k):
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
            sage: C.real(0, 0, 0, 0)
            Re(a2,0)
            sage: C.real(0, 0, 0, 1)
            Re(a2,1)
            sage: C.real(0, 0, 1, 0)
            Re(a2,0)
            sage: C.real(0, 1, 0, 0)
            Re(a2,0)
            sage: C.real(0, 2, 137/482, 0)
            Re(a0,0)
            sage: C.real(0, 0,  537/964, 0)
            Re(a3,0)
            sage: C.real(0, 1, 345/482, 0)
            Re(a1,0)
            sage: C.real(0, 1, 537/964, 0)
            Re(a4,0)

        """
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator at this point")

        return self.symbolic_ring().gen(("real", label, edge, pos, k))

    @cached_method
    def imag(self, label, edge, pos, k):
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
            sage: C.imag(0, 0, 0, 0)
            Im(a2,0)
            sage: C.imag(0, 0, 0, 1)
            Im(a2,1)
            sage: C.imag(0, 0, 0, 2)
            Im(a2,2)

        """
        if k >= self._prec:
            raise ValueError(f"symbolic ring has no {k}-th generator at this point")

        return self.symbolic_ring().gen(("imag", label, edge, pos, k))

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

            sage: C.real_part(C.gen(0, 0, 0, 0))
            Re(a2,0)
            sage: C.real_part(C.real(0, 0, 0, 0))
            Re(a2,0)
            sage: C.real_part(C.imag(0, 0, 0, 0))
            Im(a2,0)
            sage: C.real_part(2*C.gen(0, 0, 0, 0))  # tol 1e-9
            2*Re(a2,0)
            sage: C.real_part(2*I*C.gen(0, 0, 0, 0))  # tol 1e-9
            -2.0000000000000*Im(a2,0)

        """
        return self.project(x, "real")

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

            sage: C.imaginary_part(C.gen(0, 0, 0, 0))
            Im(a2,0)
            sage: C.imaginary_part(C.real(0, 0, 0, 0))  # tol 1e-9
            0
            sage: C.imaginary_part(C.imag(0, 0, 0, 0))  # tol 1e-9
            0
            sage: C.imaginary_part(2*C.gen(0, 0, 0, 0))  # tol 1e-9
            2*Im(a2,0)
            sage: C.imaginary_part(2*I*C.gen(0, 0, 0, 0))  # tol 1e-9
            2*Re(a2,0)

        """
        return self.project(x, "imag")

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
    def _formal_power_series(self, label, edge, pos, base_ring=None):
        if base_ring is None:
            base_ring = self.symbolic_ring()

        from sage.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, 'z')

        return R([self.gen(label, edge, pos, n) for n in range(self._prec)])

    def develop(self, label, edge, pos, Δ=0, base_ring=None):
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
            sage: C.develop(0, 0, 0)  # tol 1e-9
            Re(a2,0) + 1.000000000000000*I*Im(a2,0) + (Re(a2,1) + 1.000000000000000*I*Im(a2,1))*z + (Re(a2,2) + 1.000000000000000*I*Im(a2,2))*z^2
            sage: C.develop(0, 0, 0, 1)  # tol 1e-9
            Re(a2,0) + 1.000000000000000*I*Im(a2,0) + Re(a2,1) + 1.000000000000000*I*Im(a2,1) + Re(a2,2) + 1.000000000000000*I*Im(a2,2) + (Re(a2,1) + 1.000000000000000*I*Im(a2,1) + 2.000000000000000*Re(a2,2) + 2.000000000000000*I*Im(a2,2))*z + (Re(a2,2) + 1.000000000000000*I*Im(a2,2))*z^2

        """
        # TODO: Check that Δ is within the radius of convergence.
        f = self._formal_power_series(label, edge, pos, base_ring=base_ring)
        return f(f.parent().gen() + Δ)

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
            0.247323113451507*Re(a2,0) + 0.247323113451507*I*Im(a2,0) + 0.236797907979275*Re(a6,0) + 0.236797907979275*I*Im(a6,0) + 0.139540535294972*Re(a7,0) + 0.139540535294972*I*Im(a7,0) + 0.139540535294972*Re(a4,0) + 0.139540535294972*I*Im(a4,0) + 0.236797907979275*Re(a1,0) + 0.236797907979275*I*Im(a1,0)

            sage: C.integrate(b)  # tol 1e-9
            (-0.247323113451507*I)*Re(a2,0) + 0.247323113451507*Im(a2,0) + (-0.236797907979275*I)*Re(a8,0) + 0.236797907979275*Im(a8,0) + (-0.139540535294972*I)*Re(a5,0) + 0.139540535294972*Im(a5,0) + (-0.139540535294972*I)*Re(a3,0) + 0.139540535294972*Im(a3,0) + (-0.236797907979275*I)*Re(a0,0) + 0.236797907979275*Im(a0,0)

            sage: C = PowerSeriesConstraints(T, prec=5, geometry=Ω._geometry)
            sage: C.integrate(a) + C.integrate(-a)  # tol 1e-9
            0.00000000000000000
            sage: C.integrate(b) + C.integrate(-b)  # tol 1e-9
            0.00000000000000000

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

                a = gen[1]
                b = gen[2]
                P = self.complex_field()(*self._geometry.midpoint(label, edge, a, edge, b))

                P_power = P

                for k in range(self._prec):
                    expression += multiplicity * self.gen(label, edge, a, k) / (k + 1) * P_power
                    P_power *= P

                opposite_label, opposite_edge = surface.opposite_edge(label, edge)

                Q = self.complex_field()(*self._geometry.midpoint(opposite_label, opposite_edge, 1 - b, opposite_edge, 1 - a))

                Q_power = Q

                for k in range(self._prec):
                    expression -= multiplicity * self.gen(opposite_label, opposite_edge, 1-b, k) / (k + 1) * Q_power

                    Q_power *= Q

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
            edge, pos, Δ = self.relativize(label, Δ)

        for k in range(derivative, self._prec):
            value += factor * self.gen(label, edge, pos, k) * z

            factor *= k + 1
            factor /= k - derivative + 1

            z *= Δ

        return value

    @cached_method
    def _circumscribing_circle_center(self, label):
        from sage.all import CC
        return CC(*self._surface.polygon(label).circumscribing_circle().center())

    def relativize(self, label, Δ):
        r"""
        Determine the power series which is best suited to evaluate a function at Δ.

        INPUT:

        - ``label`` -- a label of a polygon in the surface

        - ``Δ`` -- a complex number describing coordinates relative to the
          center of the circumscribing circle of that polygon.

        """
        # TODO: For performance reasons we compute in CC. We should probably use complex_field()
        from sage.all import CC

        Δ = CC(Δ)
        Δ += self._circumscribing_circle_center(label)

        # TODO: This sometimes fails due to numerical noise when plotting.
        # if not self._surface.polygon(label).contains_point(Δ):
        #     raise ValueError

        from sage.all import RR

        def center(edge, pos):
            lbl, center = self._geometry.center(label, edge, pos, wrap=True, ring=RR)
            if lbl != label:
                raise NotImplementedError
            return CC(*center)

        edge, pos = min(
            [(edge, pos) for (lbl, edge, pos) in self.symbolic_ring()._gens if lbl == label],
            key=lambda edge_pos: (Δ - center(*edge_pos)).norm())

        # TODO: Once we treat centers that are in a different polygon, we need to take this into account here.
        Δ -= center(edge, pos)

        return edge, pos, self.complex_field()(Δ)

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
            [Re(a2,0) - Re(a6,0), Im(a2,0) - Im(a6,0), Re(a6,0) - Re(a7,0), Im(a6,0) - Im(a7,0), Re(a7,0) - Re(a4,0), Im(a7,0) - Im(a4,0), Re(a4,0) - Re(a1,0), Im(a4,0) - Im(a1,0), Re(a2,0) - Re(a8,0), Im(a2,0) - Im(a8,0), Re(a8,0) - Re(a5,0), Im(a8,0) - Im(a5,0), Re(a5,0) - Re(a3,0), Im(a5,0) - Im(a3,0), Re(a3,0) - Re(a0,0), Im(a3,0) - Im(a0,0)]

        If we add more coefficients, we get more complicated conditions::

            sage: C = PowerSeriesConstraints(T, prec=2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), -Re(a2,0) + 0.123661556725753*Im(a2,1) + Re(a0,0) + 0.160570808419475*Im(a0,1), -Im(a2,0) - 0.123661556725753*Re(a2,1) + Im(a0,0) - 0.160570808419475*Re(a0,1)]

        ::

            sage: C = PowerSeriesConstraints(T, prec=2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(2)
            sage: C  # tol 1e-9
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a2,1) - Re(a6,1), Im(a2,1) - Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a6,1) - Re(a7,1), Im(a6,1) - Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a7,1) - Re(a4,1), Im(a7,1) - Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), Re(a4,1) - Re(a1,1), Im(a4,1) - Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a2,1) - Re(a8,1), Im(a2,1) - Im(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a8,1) - Re(a5,1), Im(a8,1) - Im(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a5,1) - Re(a3,1), Im(a5,1) - Im(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), Re(a3,1) - Re(a0,1), Im(a3,1) - Im(a0,1)]

        """
        if derivatives > self._prec:
            raise ValueError("derivatives must not exceed global precision")

        for (label, edge), a, b in self._geometry._homology_generators:
            opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)

            parent = self.symbolic_ring()

            Δ0 = self.complex_field()(*self._geometry.midpoint(label, edge, a, edge, b))
            Δ1 = self.complex_field()(*self._geometry.midpoint(opposite_label, opposite_edge, 1-b, opposite_edge, 1-a))

            # Require that the 0th, ..., derivatives-1th derivatives are the same at the midpoint of the edge.
            for derivative in range(derivatives):
                self.add_constraint(
                    parent(self.evaluate(label, edge, a, Δ0, derivative)) - parent(self.evaluate(opposite_label, opposite_edge, 1-b, Δ1, derivative)))

    def _L2_consistency_edge(self, label, a_edge, a, b_edge, b):
        cost = self.symbolic_ring(self.real_field()).zero()

        debug = [label, a_edge, a, b_edge, b]

        def Δ(x_edge, x, y_edge, y):
            x_opposite_label, x_opposite_edge = self._surface.opposite_edge(label, x_edge)
            if x_opposite_label != label:
                raise NotImplementedError  # need to shift polygons

            y_opposite_label, y_opposite_edge = self._surface.opposite_edge(label, y_edge)
            if y_opposite_label != label:
                raise NotImplementedError  # need to shift polygons

            Δs = (
                self.complex_field()(*self._geometry.midpoint(label, x_edge, x, y_edge, y)),
                self.complex_field()(*self._geometry.midpoint(label, x_opposite_edge, 1 - x, y_edge, y)),
                self.complex_field()(*self._geometry.midpoint(label, x_edge, x, y_opposite_edge, 1 - y)),
                self.complex_field()(*self._geometry.midpoint(label, x_opposite_edge, 1 - x, y_opposite_edge, 1 - y)),
            )

            return min(Δs, key=lambda v: v.norm())

        # The weighed midpoint of the segment where the power series meet with
        # respect to the centers of the power series.
        # TODO: Assert that the min is attained for the same choice of representatives (or pick the representatives explicitly as the ones with minimal distance.)
        Δ0 = Δ(a_edge, a, b_edge, b)
        Δ1 = Δ(b_edge, b, a_edge, a)

        debug += [Δ0, Δ1]

        # Develop both power series around that midpoint, i.e., Taylor expand them.
        T0 = self.develop(label, a_edge, a, Δ0)
        T1 = self.develop(label, b_edge, b, Δ1)

        # Write b_n for the difference of the n-th coefficient of both power series.
        # We want to minimize the sum of |b_n|^2 r^2n where r is a somewhat
        # randomly chosen small radius around the midpoint.

        a_convergence = self._geometry._convergence(label, a_edge, a)
        b_convergence = self._geometry._convergence(label, b_edge, b)

        convergence = min(a_convergence - Δ0.norm(), b_convergence - Δ1.norm())

        # TODO: What should 4 be here?
        r = convergence / 4
        debug.append(r)

        assert convergence > 0
        r2 = self.real_field()(r**2)

        b = (T0 - T1).list()

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

        self._debugs.append(debug)

        return cost

    def _L2_consistency(self):
        r"""
        For each pair of adjacent centers along a homology path we use for
        integrating, let `v` be the weighed midpoint (weighed according to the
        radii of convergence at the adjacent cents.)
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

        for (label, edge), a, b in self._geometry._homology_generators:
            cost += self._L2_consistency_edge(label, edge, a, edge, b)

        # TODO: Replace these hard-coded conditions with something generic.
        # Maybe, take a Delaunay triangulation of the centers in a polygon and
        # then make sure that we have at least a condition on the four shortest
        # edges of each vertex.
        cost += self._L2_consistency_edge(0, 0, 137/482, 1, 137/482)
        cost += self._L2_consistency_edge(0, 1, 137/482, 2, 137/482)
        cost += self._L2_consistency_edge(0, 2, 137/482, 3, 137/482)
        cost += self._L2_consistency_edge(0, 3, 137/482, 0, 345/482)
        cost += self._L2_consistency_edge(0, 0, 345/482, 1, 345/482)
        cost += self._L2_consistency_edge(0, 1, 345/482, 2, 345/482)
        cost += self._L2_consistency_edge(0, 2, 345/482, 3, 345/482)
        cost += self._L2_consistency_edge(0, 3, 345/482, 0, 137/482)

        cost += self._L2_consistency_edge(0, 0, 427/964, 1, 427/964)
        cost += self._L2_consistency_edge(0, 1, 427/964, 2, 427/964)
        cost += self._L2_consistency_edge(0, 2, 427/964, 3, 427/964)
        cost += self._L2_consistency_edge(0, 3, 427/964, 0, 537/964)
        cost += self._L2_consistency_edge(0, 0, 537/964, 1, 537/964)
        cost += self._L2_consistency_edge(0, 1, 537/964, 2, 537/964)
        cost += self._L2_consistency_edge(0, 2, 537/964, 3, 537/964)
        cost += self._L2_consistency_edge(0, 3, 537/964, 0, 427/964)

        return cost

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
    #     for triangle in self._surface.label_iterator():
    #         # We expand the integrand Σa_nz^n · \overline{Σa_mz^m} naively as
    #         # the sum of a_n \overline{a_m} z^n \overline{z^m}.
    #         for n in range(self._prec):
    #             for m in range(self._prec):
    #                 coefficient = self.gen(triangle, n) * self.gen(triangle, m, conjugate=True)
    #                 # Now we have to integrate z^n \overline{z^m} on the triangle.
    #                 area += coefficient * self._elementary_area_integral(triangle, n, m)

    #     return area

    def _squares(self):
        cost = self.symbolic_ring().zero()

        for (label, edge, pos) in self.symbolic_ring()._gens:
            for k in range(1, self._prec):
                cost += self.real(label, edge, pos, k)**2
                cost += self.imag(label, edge, pos, k)**2

        return cost

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

            sage: f = 10*C.real(0, 0, 0, 0)^2 + 16*C.imag(0, 0, 0, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [32.0000000000000*Im(a2,0), 20.0000000000000*Re(a2,0)]

        In this example, the constraints and the optimized values do not
        overlap, so again we do not get Lagrange multipliers::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), -Re(a2,0) + 0.123661556725753*Im(a2,1) + Re(a0,0) + 0.160570808419475*Im(a0,1), -Im(a2,0) - 0.123661556725753*Re(a2,1) + Im(a0,0) - 0.160570808419475*Re(a0,1)]
            sage: f = 10*C.real(0, 0, 0, 0)^2 + 17*C.imag(0, 0, 0, 0)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), -Re(a2,0) + 0.123661556725753*Im(a2,1) + Re(a0,0) + 0.160570808419475*Im(a0,1), -Im(a2,0) - 0.123661556725753*Re(a2,1) + Im(a0,0) - 0.160570808419475*Re(a0,1), -λ17 + λ19, -λ16 + λ18, 0.0762270995598000*λ16 + 0.160570808419475*λ18, -0.0762270995598000*λ17 - 0.160570808419475*λ19, -λ7 + λ9, -λ6 + λ8, 0.0762270995598000*λ7 + 0.160570808419475*λ9, 0.0762270995598000*λ6 + 0.160570808419475*λ8, 34.0000000000000*Im(a2,0) + λ1 - λ9 + λ11 - λ19, 20.0000000000000*Re(a2,0) + λ0 - λ8 + λ10 - λ18, 0.123661556725753*λ1 + 0.123661556725753*λ9 + 0.123661556725753*λ10 + 0.123661556725753*λ18, 0.123661556725753*λ0 + 0.123661556725753*λ8 - 0.123661556725753*λ11 - 0.123661556725753*λ19, -λ15 + λ17, -λ14 + λ16, 0.0570539419087137*λ14 + 0.0824865933862581*λ16, -0.0570539419087137*λ15 - 0.0824865933862581*λ17, -λ5 + λ7, -λ4 + λ6, 0.0570539419087137*λ5 + 0.0824865933862581*λ7, 0.0570539419087137*λ4 + 0.0824865933862581*λ6, -λ13 + λ15, -λ12 + λ14, 0.0824865933862581*λ12 + 0.0570539419087137*λ14, -0.0824865933862581*λ13 - 0.0570539419087137*λ15, -λ1 + λ3, -λ0 + λ2, 0.160570808419475*λ1 + 0.0762270995598000*λ3, 0.160570808419475*λ0 + 0.0762270995598000*λ2, -λ3 + λ5, -λ2 + λ4, 0.0824865933862581*λ3 + 0.0570539419087137*λ5, 0.0824865933862581*λ2 + 0.0570539419087137*λ4, -λ11 + λ13, -λ10 + λ12, 0.160570808419475*λ10 + 0.0762270995598000*λ12, -0.160570808419475*λ11 - 0.0762270995598000*λ13]

        ::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_midpoint_derivatives(1)
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), -Re(a2,0) + 0.123661556725753*Im(a2,1) + Re(a0,0) + 0.160570808419475*Im(a0,1), -Im(a2,0) - 0.123661556725753*Re(a2,1) + Im(a0,0) - 0.160570808419475*Re(a0,1)]
            sage: f = 3*C.real(0, 0, 0, 1)^2 + 5*C.imag(0, 0, 0, 1)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [Re(a2,0) + 0.123661556725753*Re(a2,1) - Re(a6,0) + 0.160570808419475*Re(a6,1), Im(a2,0) + 0.123661556725753*Im(a2,1) - Im(a6,0) + 0.160570808419475*Im(a6,1), Re(a6,0) + 0.0762270995598000*Re(a6,1) - Re(a7,0) + 0.0824865933862581*Re(a7,1), Im(a6,0) + 0.0762270995598000*Im(a6,1) - Im(a7,0) + 0.0824865933862581*Im(a7,1), Re(a7,0) + 0.0570539419087137*Re(a7,1) - Re(a4,0) + 0.0570539419087137*Re(a4,1), Im(a7,0) + 0.0570539419087137*Im(a7,1) - Im(a4,0) + 0.0570539419087137*Im(a4,1), Re(a4,0) + 0.0824865933862581*Re(a4,1) - Re(a1,0) + 0.0762270995598000*Re(a1,1), Im(a4,0) + 0.0824865933862581*Im(a4,1) - Im(a1,0) + 0.0762270995598000*Im(a1,1), -Re(a2,0) + 0.123661556725753*Re(a2,1) + Re(a1,0) + 0.160570808419475*Re(a1,1), -Im(a2,0) + 0.123661556725753*Im(a2,1) + Im(a1,0) + 0.160570808419475*Im(a1,1), Re(a2,0) + 0.123661556725753*Im(a2,1) - Re(a8,0) + 0.160570808419475*Im(a8,1), Im(a2,0) - 0.123661556725753*Re(a2,1) - Im(a8,0) - 0.160570808419475*Re(a8,1), Re(a8,0) + 0.0762270995598000*Im(a8,1) - Re(a5,0) + 0.0824865933862581*Im(a5,1), Im(a8,0) - 0.0762270995598000*Re(a8,1) - Im(a5,0) - 0.0824865933862581*Re(a5,1), Re(a5,0) + 0.0570539419087137*Im(a5,1) - Re(a3,0) + 0.0570539419087137*Im(a3,1), Im(a5,0) - 0.0570539419087137*Re(a5,1) - Im(a3,0) - 0.0570539419087137*Re(a3,1), Re(a3,0) + 0.0824865933862581*Im(a3,1) - Re(a0,0) + 0.0762270995598000*Im(a0,1), Im(a3,0) - 0.0824865933862581*Re(a3,1) - Im(a0,0) - 0.0762270995598000*Re(a0,1), -Re(a2,0) + 0.123661556725753*Im(a2,1) + Re(a0,0) + 0.160570808419475*Im(a0,1), -Im(a2,0) - 0.123661556725753*Re(a2,1) + Im(a0,0) - 0.160570808419475*Re(a0,1), -λ17 + λ19, -λ16 + λ18, 0.0762270995598000*λ16 + 0.160570808419475*λ18, -0.0762270995598000*λ17 - 0.160570808419475*λ19, -λ7 + λ9, -λ6 + λ8, 0.0762270995598000*λ7 + 0.160570808419475*λ9, 0.0762270995598000*λ6 + 0.160570808419475*λ8, λ1 - λ9 + λ11 - λ19, λ0 - λ8 + λ10 - λ18, 10.0000000000000*Im(a2,1) + 0.123661556725753*λ1 + 0.123661556725753*λ9 + 0.123661556725753*λ10 + 0.123661556725753*λ18, 6.00000000000000*Re(a2,1) + 0.123661556725753*λ0 + 0.123661556725753*λ8 - 0.123661556725753*λ11 - 0.123661556725753*λ19, -λ15 + λ17, -λ14 + λ16, 0.0570539419087137*λ14 + 0.0824865933862581*λ16, -0.0570539419087137*λ15 - 0.0824865933862581*λ17, -λ5 + λ7, -λ4 + λ6, 0.0570539419087137*λ5 + 0.0824865933862581*λ7, 0.0570539419087137*λ4 + 0.0824865933862581*λ6, -λ13 + λ15, -λ12 + λ14, 0.0824865933862581*λ12 + 0.0570539419087137*λ14, -0.0824865933862581*λ13 - 0.0570539419087137*λ15, -λ1 + λ3, -λ0 + λ2, 0.160570808419475*λ1 + 0.0762270995598000*λ3, 0.160570808419475*λ0 + 0.0762270995598000*λ2, -λ3 + λ5, -λ2 + λ4, 0.0824865933862581*λ3 + 0.0570539419087137*λ5, 0.0824865933862581*λ2 + 0.0570539419087137*λ4, -λ11 + λ13, -λ10 + λ12, 0.160570808419475*λ10 + 0.0762270995598000*λ12, -0.160570808419475*λ11 - 0.0762270995598000*λ13]

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
        for (label, edge, pos) in self.symbolic_ring()._gens:
            for k in range(self._prec):
                for kind in ["imag", "real"]:
                    gen = self.symbolic_ring().gen((kind, label, edge, pos, k))

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
            [0.247323113451507*Re(a2,0) + 0.236797907979275*Re(a6,0) + 0.139540535294972*Re(a7,0) + 0.139540535294972*Re(a4,0) + 0.236797907979275*Re(a1,0), 0.247323113451507*Im(a2,0) + 0.236797907979275*Im(a8,0) + 0.139540535294972*Im(a5,0) + 0.139540535294972*Im(a3,0) + 0.236797907979275*Im(a0,0) - 1.00000000000000]

        If we increase precision, we see additional higher imaginary parts.
        These depend on the choice of base point of the integration and will be
        found to be zero by other constraints, not true anymore TODO::

            sage: C = PowerSeriesConstraints(T, 2, geometry=Ω._geometry)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9
            [0.247323113451507*Re(a2,0) + 0.236797907979275*Re(a6,0) - 0.00998620690459202*Re(a6,1) + 0.139540535294972*Re(a7,0) - 0.00177444290057350*Re(a7,1) + 0.139540535294972*Re(a4,0) + 0.00177444290057350*Re(a4,1) + 0.236797907979275*Re(a1,0) + 0.00998620690459202*Re(a1,1), 0.247323113451507*Im(a2,0) + 0.236797907979275*Im(a8,0) + 0.00998620690459202*Re(a8,1) + 0.139540535294972*Im(a5,0) + 0.00177444290057350*Re(a5,1) + 0.139540535294972*Im(a3,0) - 0.00177444290057350*Re(a3,1) + 0.236797907979275*Im(a0,0) - 0.00998620690459202*Re(a0,1) - 1.00000000000000]

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
            sage: f = 10*C.real(0, 0, 0, 0)^2 + 16*C.imag(0, 0, 0, 0)^2
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

        if len(set(self.symbolic_ring().gen((kind, label, edge, pos, k)) for kind in ["real", "imag"] for (label, edge, pos) in self.symbolic_ring()._gens for k in range(self._prec))) != len(non_lagranges):
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
    def power_series_ring(self, label, edge, pos):
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
            sage: Ω.power_series_ring(0, 0, 0)
            Power Series Ring in z0_0_1 over Complex Field with 54 bits of precision

        """
        from sage.all import PowerSeriesRing
        polygon = list(self._surface.label_iterator()).index(label)
        pos = [p for (lbl, e, p) in self.symbolic_ring()._gens if lbl == label and e == edge].index(pos)

        return PowerSeriesRing(self.complex_field(), f"z{polygon}_{edge}_{pos}")

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
            sage: C.add_constraint(C.real(0, 0, 0, 0) - C.real(0, 0, 0, 1))
            sage: C.add_constraint(C.real(0, 0, 0, 0) - 1)
            sage: C.solve()  # rondam output due to random ordering of the dict
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

        from sage.all import RDF
        condition = A.change_ring(RDF).condition()

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

            solution = cppyy.gbl.solve(A, b)

            from sage.all import vector
            solution = vector([self.real_field()(entry) for entry in solution])
        else:
            raise NotImplementedError

        residue = (A*solution - b).norm()

        lagranges = len(self.lagrange_variables())

        if lagranges:
            solution = solution[:-lagranges]

        P = len(self.symbolic_ring()._gens)
        real_solution = [[solution[decode[2*p + 2*P*prec]] if 2*p + 2*P*prec in decode else 0 for prec in range(self._prec)] for p in range(P)]
        imag_solution = [[solution[decode[2*p + 1 + 2*P*prec]] if 2*p + 1 + 2*P*prec in decode else 0 for prec in range(self._prec)] for p in range(P)]

        return {
            (label, edge, pos): self.power_series_ring(label, edge, pos)([self.complex_field()(real_solution[p][k], imag_solution[p][k]) for k in range(self._prec)], self._prec)
            for (p, (label, edge, pos)) in enumerate(self.symbolic_ring()._gens)
        }, residue
