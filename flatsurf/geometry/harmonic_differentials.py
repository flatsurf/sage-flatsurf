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
    sage: ω = Ω(f)
    sage: ω
    (1.0000 + O(z0^5), 1.0000 + O(z1^5))

The harmonic differential that integrates as 0 along `a` but 1 along `b`::

    sage: g = H({b: 1})
    sage: Ω(g)  # tol 1e-9
    (1.0000*I + O(z0^5), 1.0000*I + O(z1^5))

A less trivial example, the regular octagon::

    sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
    sage: S = translation_surfaces.regular_octagon()

    sage: H = SimplicialHomology(S)
    sage: HS = SimplicialCohomology(S, homology=H)
    sage: a, b, c, d = HS.homology().gens()

    sage: f = { a: -1, b: sqrt(2), c: -1, d: 0}

    sage: f = HS(f)
    sage: f._values = {key: RealField(54)(value) for (key, value) in f._values.items()}  # TODO: Why is this hack necessary?

    sage: Omega = HarmonicDifferentials(S)
    sage: omega = Omega(HS(f), prec=3, check=False)
    sage: omega  # TODO: Increase precision once this is faster.
    ((-0.000048417 + 0.000014772*I) + 0.00092363*I*z0 + (1.3146 - 0.000033146*I)*z0^2 + O(z0^3), (-2.0500 + 0.000051690*I) - 0.0014250*I*z1 + (-0.00014525 + 0.000044316*I)*z1^2 + 0.0049924*I*z1^3 - 0.0041365*I*z1^5 + 0.011334*z1^6 - 0.0018723*I*z1^7 + (-3.0794 + 0.000077645*I)*z1^8 + O(z1^9))


The same computation on a triangulation of the octagon::

    sage: from flatsurf import HarmonicDifferentials, SimplicialHomology, SimplicialCohomology, Polygon, translation_surfaces
    sage: S = translation_surfaces.regular_octagon()
    sage: S = S.subdivide().codomain()

    sage: H = SimplicialHomology(S)
    sage: HS = SimplicialCohomology(S, homology=H)
    sage: a, b, c, d = HS.homology().gens()

    sage: f = { a: -1, b: sqrt(2), c: -1, d: 0}

    sage: f = HS(f)
    sage: f._values = {key: RealField(54)(value) for (key, value) in f._values.items()}  # TODO: Why is this hack necessary?

    sage: Omega = HarmonicDifferentials(S, safety=0, singularities=True, centers=False)
    sage: omega = Omega(HS(f), prec=3, check=False)
    sage: omega  # TODO: Increase precision once this is faster.

The same surface but built as the unfolding of a right triangle::

    sage: from flatsurf import similarity_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology, Polygon
    sage: S = similarity_surfaces.billiard(Polygon(angles=[3/8, 1/2, 1/8], lengths=[1/2])).minimal_cover('translation')

    sage: H = SimplicialHomology(S)
    sage: HS = SimplicialCohomology(S, homology=H)
    sage: a, b, c, d = H.gens()

    sage: f = { a: -1, b: sqrt(2), c: -1, d: 0}

    sage: f = HS(f)
    sage: f._values = {key: RealField(54)(value) for (key, value) in f._values.items()}  # TODO: Why is this hack necessary?

    sage: Omega = HarmonicDifferentials(S, safety=0, singularities=True, centers=False)
    sage: omega = Omega(HS(f), prec=3, check=False)
    sage: omega  # TODO: Increase precision once this is faster.

Much more complicated, the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import similarity_surfaces, SimplicialHomology, SimplicialCohomology, HarmonicDifferentials, Polygon

    sage: S = similarity_surfaces.billiard(Polygon(angles=[3, 4, 13])).minimal_cover("translation")
    sage: S = S.erase_marked_points().delaunay_decomposition()

    sage: H = SimplicialHomology(S)
    sage: HS = SimplicialCohomology(S, homology=H)
    sage: gens = H.gens()

    sage: f = { g: 0 for g in gens}
    sage: f[gens[0]] = 1
    sage: f = HS(f)
    sage: f._values = {key: RealField(54)(value) for (key, value) in f._values.items()}  # TODO: Why is this hack necessary?

    sage: Omega = HarmonicDifferentials(S)
    sage: omega = Omega(HS(f), prec=3, check=False)

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

import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    import cppyy

complex = None


def _cppyy():
    global complex
    if complex is None:
        cppyy.include('complex')
        complex = cppyy.gbl.std.complex['double']
    return cppyy


def Ccpp(x):
    _cppyy()
    return complex(float(x.real()), float(x.imag()))


def integral2arb(part, α, κ, d, ζd, n, β, λ, dd, ζdd, m, a, b, C, R):
    r"""
    Return the real/imaginary part of

    \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} \overline{ζ_{dd+1}^{λ (m+1)}/(dd+1) (z-β)^\frac{m-dd}{dd+1}} dz

    where γ(t) = (1-t)a + tb.
    """
    if not hasattr(_cppyy().gbl, "integral2arb"):
        _cppyy().cppdef(r"""
        #include <acb_calc.h>
        #include <string>

        const int ARB_PREC=32;

        struct Arb {
            Arb() {
                arb_init(value);
            }

            Arb(double d) {
                arb_init(value);
                arb_set_d(value, d);
            }

            ~Arb() {
                arb_clear(value);
            }

            operator const arb_t&() const {
                return value;
            }

            operator arb_t&() {
                return value;
            }

            operator double() {
                return arf_get_d(arb_midref(value), ARF_RND_NEAR);
            }

            arb_t value;
        };

        struct Mag {
            Mag(double d) {
                mag_init(value);
                mag_set_d(value, d);
            }

            ~Mag() {
                mag_clear(value);
            }

            operator const mag_t&() const {
                return value;
            }

            mag_t value;
        };

        struct Acb {
            Acb() {
                acb_init(value);
            }

            Acb(const acb_t v) {
                acb_init(value);
                acb_set(value, v);
            }

            Acb(double re, double im) {
                acb_init(value);
                acb_set_d_d(value, re, im);
            }

            Arb real() {
                Arb real;
                acb_get_real(real, value);
                return real;
            }

            Arb imag() {
                Arb imag;
                acb_get_imag(imag, value);
                return imag;
            }

            ~Acb() {
                acb_clear(value);
            }

            operator const acb_t&() const {
                return value;
            }

            operator acb_t&() {
                return value;
            }

            acb_t value;
        };

        struct Args {
            Args(double Re_alpha, double Im_alpha, int kappa, int d, double Re_zeta_d, double Im_zeta_d, int n, double Re_beta, double Im_beta, int lambda, int dd, double Re_zeta_dd, double Im_zeta_dd, int m, double Re_a, double Im_a, double Re_b, double Im_b) :
                alpha(Re_alpha, Im_alpha),
                beta(Re_beta, Im_beta),
                d(d),
                dd(dd),
                n(n),
                m(m) {

                Acb zeta_d_power(Re_zeta_d, Im_zeta_d);
                acb_pow_si(zeta_d_power, zeta_d_power, kappa * (n + 1), ARB_PREC);
                acb_div_si(zeta_d_power, zeta_d_power, d + 1, ARB_PREC);

                Acb zeta_dd_power(Re_zeta_dd, Im_zeta_dd);
                acb_pow_si(zeta_dd_power, zeta_dd_power, lambda * (m + 1), ARB_PREC);
                acb_div_si(zeta_dd_power, zeta_dd_power, dd + 1, ARB_PREC);
                acb_conj(zeta_dd_power, zeta_dd_power);

                acb_sub(ba, Acb(Re_b, Im_b), Acb(Re_a, Im_a), ARB_PREC);

                acb_abs(ba_, ba, ARB_PREC);

                Acb ba__;
                acb_set_round_arb(ba__, ba_, ARB_PREC);

                acb_mul(constant, ba__, zeta_d_power, ARB_PREC);
                acb_mul(constant, constant, zeta_dd_power, ARB_PREC);
            }

            Acb alpha, beta;
            Acb ba;
            Arb ba_;
            int d, dd;
            int n, m;
            Acb constant;
        };

        double integral2arb(std::string part, double Re_alpha, double Im_alpha, int kappa, int d, double Re_zeta_d, double Im_zeta_d, int n, double Re_beta, double Im_beta, int lambda, int dd, double Re_zeta_dd, double Im_zeta_dd, int m, double Re_a, double Im_a, double Re_b, double Im_b) {
            double res_d;

            Acb res;

            Args args(Re_alpha, Im_alpha, kappa, d, Re_zeta_d, Im_zeta_d, n, Re_beta, Im_beta, lambda, dd, Re_zeta_dd, Im_zeta_dd, m, Re_a, Im_a, Re_b, Im_b);

            const auto func = [](acb_ptr out, const acb_t z, void * param, slong order, slong prec) -> int {
                const Args* args = (const Args*)param;

                if (order > 1) {
                    // TODO: When order == 1 we need to verify that we are holomorphic somewhere. But where exactly?
                    throw std::logic_error("derivatives of function not implemented/function not holomorphic");
                }

                Acb za(z);
                acb_sub(za, za, args->alpha, ARB_PREC);
                acb_pow_arb(za, za, Arb(1 / (double)(args->d + 1)), ARB_PREC);
                acb_pow_si(za, za, args->n - args->d, ARB_PREC);

                Acb zb(z);
                acb_sub(zb, zb, args->beta, ARB_PREC);
                acb_pow_arb(zb, zb, Arb(1 / (double)(args->dd + 1)), ARB_PREC);
                acb_pow_si(zb, zb, args->m - args->dd, ARB_PREC);
                acb_conj(zb, zb);

                acb_mul(out, za, zb, ARB_PREC);

                return 0;
            };

            if (acb_calc_integrate(res, func, &args, Acb(Re_a, Im_a), Acb(Re_b, Im_b), ARB_PREC /* rel_goal */, Mag(1e-64) /* abs_tol */, nullptr /* options */, ARB_PREC) != ARB_CALC_SUCCESS) {
                throw std::logic_error("acb_calc_integrate() did not converge");
            }

            acb_mul(res, res, args.constant, ARB_PREC);
            acb_div(res, res, args.ba, ARB_PREC);

            if (part == "Re") {
                return res.real();
            } else if (part == "Im") {
                return res.imag();
            } else {
                throw std::logic_error("unknown part");
            }
        }
        """)

        _cppyy().load_library("arb");

    return R(_cppyy().gbl.integral2arb(part, float(α.real()), float(α.imag()), int(κ), int(d), float(ζd.real()), float(ζd.imag()), int(n), float(β.real()), float(β.imag()), int(λ), int(dd), float(ζdd.real()), float(ζdd.imag()), m, float(a.real()), float(a.imag()), float(b.real()), float(b.imag())))


def integral2cpp(part, α, κ, d, ζd, n, β, λ, dd, ζdd, m, a, b, C, R):
    r"""
    Return the real/imaginary part of

    \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} \overline{ζ_{dd+1}^{λ (m+1)}/(dd+1) (z-β)^\frac{m-dd}{dd+1}} dz

    where γ(t) = (1-t)a + tb.
    """
    # Since γ(t) = (1 - t)a + tb, we have |·γ(t)| = |b - a|
    constant = ζd**(κ * (n+1)) / (d+1) * (ζdd**(λ * (m+1)) / (dd+1)).conjugate() * abs(b - a)

    constant = Ccpp(constant)

    α = Ccpp(α)
    β = Ccpp(β)
    a = Ccpp(a)
    b = Ccpp(b)

    n = int(n)
    m = int(m)
    d = int(d)
    dd = int(dd)

    if not hasattr(_cppyy().gbl, "value"):
        _cppyy().cppdef(r"""
        #include <complex>

        using complex = std::complex<double>;

        complex pow(complex z, double e) {
            if (e == 0) return 1;
            if (e == 1) return z;
            return std::pow(z, e);
        }

        complex value(double t, complex constant, complex a, complex alpha, int d, int n, complex b, complex beta, int dd, int m) {
            complex z = (1 - t) * a + t * b;

            complex za = pow(pow(z - alpha, 1 / (double)(d + 1)), n - d);
            complex zb = std::conj(pow(pow(z - beta, 1 / (double)(dd + 1)), m - dd));

            return constant * za * zb;
        }
        """)

    def pow(z, e):
        z = complex(z)
        e = complex(e)

        if e == 0:
            return complex(1)

        if e == 1:
            return z

        return _cppyy().gbl.std.pow(z, e)

    def value(t):
        value = _cppyy().gbl.value(t, constant, a, α, d, n, b, β, dd, m)

        if part == "Re":
            return value.real
        if part == "Im":
            return value.imag

        raise NotImplementedError

    from scipy.integrate import quad
    integral, error = quad(value, 0, 1)

    return R(integral)


def integral2(part, α, κ, d, ζd, n, β, λ, dd, ζdd, m, a, b, C, R):
    r"""
    Return the real/imaginary part of

    \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} \overline{ζ_{dd+1}^{λ (m+1)}/(dd+1) (z-β)^\frac{m-dd}{dd+1}} dz

    where γ(t) = (1-t)a + tb.
    """
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
    if not hasattr(_cppyy().gbl, "solve"):
        _cppyy().cppdef(r'''
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

    return _cppyy().gbl.solve


# TODO: This works around a problem when pyflatsurf is loaded. If pyflatsurf is loaded first, there are C++ errors.
define_solve()


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

    def _sub_(self, other):
        return self.parent()({
            triangle: self._series[triangle] - other._series[triangle]
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
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True, centers=True)
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

        from flatsurf.geometry.power_series import SymbolicPowerSeries
        if isinstance(value, SymbolicPowerSeries):
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

    def integrate(self, cycle, numerical=False, part=None):  # TODO: Remove debugging hack "part"
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
        return self._evaluate(C.integrate(cycle, part=part))

    def _repr_(self, prec=1e-6):
        # TODO: Tune this so we do not loose important information.
        # TODO: Do not print that many digits.
        # TODO: Make it visually clear that we are not printing everything.
        def compress(series):
            def compress_coefficient(coefficient):
                if coefficient.imag().abs() < prec:
                    coefficient = coefficient.parent()(coefficient.real())
                if coefficient.real().abs() < prec:
                    coefficient = coefficient.parent()(coefficient.parent().gen() * coefficient.imag())

                return coefficient

            from sage.all import ComplexField
            return series.parent()({exponent: compress_coefficient(coefficient) for (coefficient, exponent) in zip(series.coefficients(), series.exponents())} or 0).change_ring(ComplexField(20)).add_bigoh(series.precision_absolute())

        ret = repr(tuple(compress(series) for series in self._series.values()))

        # TODO: Why are zero coefficients not filtered automatically?
        import re
        ret = re.sub(r'[+-] 0\.0*\*z[01](\^\d*|) ', '', ret)

        return ret

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
    def __classcall__(cls, surface, centers=None, homology_generators=None, category=None):
        r"""
        Normalize parameters when creating the space of harmonic differentials.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: HarmonicDifferentials(T) is HarmonicDifferentials(T)
            True

        """
        centers = HarmonicDifferentials._centers(surface, centers)
        homology_generators = HarmonicDifferentials._homology_generators(surface, centers, homology_generators)
        return super().__classcall__(cls, surface, centers, homology_generators, category or SetsWithPartialMaps())

    def __init__(self, surface, centers, homology_generators, category):
        Parent.__init__(self, category=category)

        self._surface = surface

        self._geometry = GeometricPrimitives(surface, centers, homology_generators)

    @staticmethod
    def _centers(surface, algorithm):
        if algorithm is None:
            algorithm = "vertices+centers"

        if algorithm == "vertices":
            return HarmonicDifferentials._centers_vertices(surface)

        if algorithm == "vertices+centers":
            return HarmonicDifferentials._centers_vertices_and_centers(surface)

        from flatsurf.geometry.surface_objects import SurfacePoint
        if all(isinstance(center, SurfacePoint) for center in algorithm):
            return frozenset(algorithm)

        raise NotImplementedError("unsupported algorithm for determining centers of power series")

    @staticmethod
    def _centers_vertices(surface):
        r"""
        Return the vertices of ``surface`` to develop power series around them.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: S = translation_surfaces.regular_octagon()

            sage: HarmonicDifferentials._centers_vertices(S)
            frozenset({Vertex 0 of polygon 0})

        """
        return frozenset(surface.vertices())

    @staticmethod
    def _centers_vertices_and_centers(surface):
        return frozenset(list(surface.vertices()) + [surface(label, surface.polygon(label).centroid()) for label in surface.labels()])

    @staticmethod
    def _homology_generators(surface, centers, algorithm):
        if algorithm is None:
            algorithm = "centers"

        if algorithm == "centers":
            return HarmonicDifferentials._homology_generators_centers(surface)

        if all(isinstance(generator, Path) for generator in algorithm):
            return frozenset(algorithm)

        raise NotImplementedError("unsupported algorithm for determining homology generators")

    @staticmethod
    def _homology_generators_centers(surface):
        r"""
        Return all possible paths (up to orientation) connecting the centroids
        of neighboring polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: S = translation_surfaces.regular_octagon()

            sage: gens = HarmonicDifferentials._homology_generators_centers(S)
            sage: len(gens)
            8

        """
        generators = set()
        for label in surface.labels():
            polygon = surface.polygon(label)
            for edge in range(len(polygon.edges())):
                path = GeodesicPath.across_edge(surface, label, edge)
                if -path not in generators:
                    generators.add(path)

        return frozenset(generators)

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


# TODO: This is overlapping too much with Voronoi cell constructions by now.
class GeometricPrimitives(UniqueRepresentation):
    def __init__(self, surface, centers, homology_generators):
        if surface.is_mutable():
            raise TypeError("surface must be immutable")

        self._surface = surface
        self._homology_generators = homology_generators
        self._centers = list(centers)

    def path_across_edge(self, label, edge):
        r"""
        Return the path from the centroid of the polygon ``label`` to the
        centroid of the polygon on the other side of ``edge`` into a
        :class:`Path` that we can integrate along.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: S = translation_surfaces.regular_octagon()
            sage: Ω = HarmonicDifferentials(S)
            sage: Ω._geometry.path_across_edge(0, 0)
            [(1, Path (0, -a - 1) from (1/2, 1/2*a + 1/2) in polygon 0 to (1/2, 1/2*a + 1/2) in polygon 0)]

        """
        path = GeodesicPath.across_edge(self._surface, label, edge)
        if path in self._homology_generators:
            return [(1, path)]
        if -path in self._homology_generators:
            return [(-1, -path)]

        raise NotImplementedError("cannot rewrite this path as a sum of supported paths yet")

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
        raise NotImplementedError  # use voronoi functionality instead
        center_point = self._surface(*center)
        if not center_point.is_vertex():
            return 0, 0

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

        gens = [f"Re(a{n},?)" for n in range(len(self._geometry._centers))] + [f"Im(a{n},?)" for n in range(len(self._geometry._centers))] + ["λ?"]

        from sage.all import ComplexField
        from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
        return PowerSeriesCoefficientExpressionRing(base_ring or ComplexField(54), tuple(gens))

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
    def lagrange(self, k):
        return self.symbolic_ring().gen(("λ?", k))

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

    def integrate(self, cycle, part=None):  # TODO: Remove (or rework) debugging hack "part"
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
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True, centers=True)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=1, geometry=Ω._geometry)

            sage: C.integrate(H())
            0.000000000000000

            sage: a, b, c, d = H.gens()
            sage: C.integrate(b)  # TODO: It's worthwhile to check these manually. Any slight error here is going to propagate into the differential.
            1.60217526524068*Re(a0,0) + 1.60217526524068*I*Im(a0,0) + (-0.389300863573646 + 1.38777878078145e-17*I)*Re(a1,0) + (-1.38777878078145e-17 - 0.389300863573646*I)*Im(a1,0) + (-1.04083408558608e-17 - 0.327551243899061*I)*Re(a1,1) + (0.327551243899061 - 1.04083408558608e-17*I)*Im(a1,1) + 0.270679432377470*Re(a1,2) + 0.270679432377470*I*Im(a1,2)

            sage: C.integrate(d)
            (-1.60217526524068*I)*Re(a0,0) + 1.60217526524068*Im(a0,0) + (5.55111512312578e-17 - 0.389300863573646*I)*Re(a1,0) + (0.389300863573646 + 5.55111512312578e-17*I)*Im(a1,0) + (-2.77555756156289e-17 + 0.327551243899061*I)*Re(a1,1) + (-0.327551243899061 - 2.77555756156289e-17*I)*Im(a1,1) + (-0.270679432377470*I)*Re(a1,2) + 0.270679432377470*Im(a1,2)

            sage: C.integrate(a)
            (1.13290899470104 - 1.13290899470104*I)*Re(a0,0) + (1.13290899470104 + 1.13290899470104*I)*Im(a0,0) + (0.275277280554704 + 0.275277280554704*I)*Re(a1,0) + (-0.275277280554704 + 0.275277280554704*I)*Im(a1,0) + (0.327551243899061 + 6.93889390390723e-18*I)*Re(a1,1) + (-6.93889390390723e-18 + 0.327551243899061*I)*Im(a1,1) + (0.191399262161835 - 0.191399262161835*I)*Re(a1,2) + (0.191399262161835 + 0.191399262161835*I)*Im(a1,2)

        """
        return sum(multiplicity * sgn * self._integrate_path(path) for (label, edge), multiplicity in cycle._chain.monomial_coefficients().items() for (sgn, path) in self._geometry.path_across_edge(label, edge))

    def _integrate_path(self, path):
        r"""
        Return the linear combination of the power series coefficients that
        describe the integral along the ``path``.

        This is a helper method for :meth:`integrate`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: Ω = HarmonicDifferentials(S)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=1, geometry=Ω._geometry)

            sage: from flatsurf.geometry.harmonic_differentials import GeodesicPath
            sage: path = GeodesicPath.across_edge(S, 0, 0)

            sage: C._integrate_path(path)

        """
        return sum(self._integrate_path_polygon(label, segment) for (label, segment) in path.split())

    def _integrate_path_polygon(self, label, segment):
        r"""
        Return a symbolic expression describing the integral along the
        ``segment`` in the polygon with ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: Ω = HarmonicDifferentials(S)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=1, geometry=Ω._geometry)

            sage: from flatsurf.geometry.euclidean import OrientedSegment

            sage: C._integrate_path_polygon(0, OrientedSegment((0, 0), (1, 0))

        """
        V = self._voronoi_diagram()

        return sum(self._integrate_path_cell(cell, subsegment) for (cell, subsegment) in V.split_segment(label, segment).items())

    def _integrate_path_cell(self, polygon_cell, segment):
        r"""
        Return a symbolic expression describing the integral along the
        ``segment`` in the Voronoi cell ``polygon_cell``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: Ω = HarmonicDifferentials(S)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=1, geometry=Ω._geometry)

            sage: from flatsurf.geometry.euclidean import OrientedSegment

            sage: V = C._voronoi_diagram()
            sage: C._integrate_path_cell(V.polygon_cell(0, (0, 0)), OrientedSegment((0, 0), (1/2, 0)))

        """
        integrator = self.CellIntegrator(self, polygon_cell)
        return sum(integrator.a(n) * integrator.integral(n, segment) for n in range(self._surface(polygon_cell.label(), polygon_cell.center()).angle() * self._prec))

    # TODO: Remove
    def _integrate_across_edge(self, label, edge):
        r"""
        Return a symbolic expression describing the integral from the center of
        the circumscribing circle of the polygon with ``label`` to the center
        of the circumscribing circle of the polygon on the other side of
        ``edge``.
        """
        opposite_label, opposite_edge = self._surface.opposite_edge(label, edge)

        # TODO: Currently, there is no concept of a segment in a surface in
        # sage-flatsurf other than a saddle connection or a flow starting at a
        # vertex. It would be nice to have a proper type for such objects.

        # We integrate from the center of the circumscribing circle of one
        # polygon to the edge and then from the edge to the center of the other
        # circumscribing circle.
        return self._integrate_center_to_edge(label, edge) - self._integrate_center_to_edge(opposite_label, opposite_edge)

    # TODO: Remove
    def _integrate_center_to_edge(self, label, edge):
        r"""
        Return a symbolic expression describing the integral from the center of
        the circumscribing circle of the polygon with ``label`` to the midpoint
        of the ``edge``.
        """
        # Split the segment from the center to the edge into segments that are
        # in a single Voronoi cell.
        V = self._voronoi_diagram()

        polygon = self._surface.polygon(label)

        expression = self.symbolic_ring().zero()

        for cell, segment in V.split_segment(label, polygon.circumscribing_circle().center(), polygon.vertices()[edge] + polygon.edges()[edge]/2).items():
            integrator = self.CellIntegrator(self, cell)

            for n in range(self._surface(label, edge).angle() * self._prec):
                expression += integrator.a(n) * integrator.f(n, segment)

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
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True, centers=True)

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
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True, centers=True)

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
        from sage.all import parallel

        @parallel
        def L2_cost(cells, boundary):
            cell, opposite_cell = list(cells)
            boundary = cell.split_segment_uniform_root_branch(boundary)
            boundary = sum([opposite_cell.split_segment_uniform_root_branch(segment) for segment in boundary], [])
            return sum(self._L2_consistency_voronoi_boundary(cell, segment, opposite_cell) for segment in boundary)

        return sum(cost for _, cost in L2_cost(list(self._voronoi_diagram().boundaries().items())))

    @cached_method
    def ζ(self, d):
        from sage.all import exp, pi, I
        return self.complex_field()(exp(2*pi*I / d))

    def _gen(self, kind, center, n):
        return self.symbolic_ring(self.real_field()).gen((f"{kind}(a{self._geometry._centers.index(center)},?)", n))

    class CellIntegrator:
        def __init__(self, constraints, cell):
            self._constraints = constraints
            self._cell = cell

        @cached_method
        def a(self, n):
            return self.Re_a(n) + self._constraints.complex_field().gen() * self._constraints.symbolic_ring(self._constraints.complex_field())(self.Im_a(n))

        @cached_method
        def Re_a(self, n):
            return self._constraints._gen("Re", self._cell.surface()(self._cell.label(), self._cell.center()), n)

        @cached_method
        def Im_a(self, n):
            return self._constraints._gen("Im", self._cell.surface()(self._cell.label(), self._cell.center()), n)

        def integral(self, n, segment):
            r"""
            Return the sum of

            \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1}

            where `d` is the order of the singularity at the center of the
            Voronoi cell containing ``segment``, γ are segments that partition
            ``segment`` such that the `d+1`st root of `z-α` can be taken
            consistently on the smaller segment and it is precisely the main
            branch of the root up to `ζ^κ` where `ζ` is a `d+1`st root of
            unity, `α` is the center of the Voronoi cell containing the
            segment,
            """
            sum = self._constraints.complex_field().zero()

            C = self._constraints.complex_field()

            d = self._cell.surface()(self._cell.label(), self._cell.center()).angle() - 1
            α = C(*self._cell.center())

            for γ in self._cell.split_segment_uniform_root_branch(segment):
                a = C(*γ.start())
                b = C(*γ.end())

                constant = self._constraints.ζ(d + 1) ** (self._cell.root_branch(γ) * (n + 1)) / (d + 1) * (C(b) - C(a))

                def value(part, t):
                    z = self._constraints.complex_field()(*((1 - t) * a + t * b))

                    value = constant * (z-α).nth_root(d + 1)**(n - d)

                    if part == "Re":
                        return float(value.real())
                    if part == "Im":
                        return float(value.imag())

                    raise NotImplementedError

                from scipy.integrate import quad
                real, error = quad(lambda t: value("Re", t), 0, 1)
                imag, error = quad(lambda t: value("Im", t), 0, 1)

                sum += C(real, imag)

            return sum

        def mixed_integral(self, part, α, κ, d, n, β, λ, dd, m, γ):
            # TODO: In the actual computation we are throwing an absolute value
            # in somewhere, i.e., we are not using γ.(t) but |γ.(t)|.
            # That's fine but the documentation is wrong here and some other places.
            # TODO: It's weird that this one does not split for inuform roots but integral() does.
            r"""
            Return the real/imaginary part of

            \int_γ ζ_{d+1}^{κ (n+1)}/(d+1) (z-α)^\frac{n-d}{d+1} \overline{ζ_{dd+1}^{λ (m+1)}/(dd+1) (z-β)^\frac{m-dd}{dd+1}} dz
            """
            R = self._constraints.real_field()
            C = self._constraints.complex_field()
            a = Cab(C, γ.start())
            b = Cab(C, γ.end())
            α = Cab(C, α)
            β = Cab(C, β)

            return integral2cpp(part, α, κ, d, self._constraints.ζ(d + 1), n, β, λ, dd, self._constraints.ζ(dd + 1), m, a=a, b=b, C=C, R=R)

    class CellBoundaryIntegrator:
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
            a, b = self._segment.endpoints()
            a = self.complex_field(*a)
            b = self.complex_field(*b)

            # Since γ(t) = (1 - t)a + tb, we have ·γ(t) = b - a
            constant = self.ζ(d + 1)**(κ * (n + 1)) / (d + 1) * (b - a)

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

        def α(self):
            return self.complex_field(*self._center_coordinates)

        def β(self):
            return self.complex_field(*self._opposite_center_coordinates)

        def κ(self):
            if not self._center.is_vertex():
                return 0

            return self._κλ(self._center_coordinates, self._segment)

        @cached_method
        def d(self):
            if not self._center.is_vertex():
                return 0

            return 2

        def λ(self):
            if not self._opposite_center.is_vertex():
                return 0

            return self._κλ(self._opposite_center_coordinates, self._segment)

        @cached_method
        def dd(self):
            if not self._opposite_center.is_vertex():
                return 0

            return 2

        def _κλ(self, center, segment):
            raise NotImplementedError  # hardcoded OCTAGON
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

    @cached_method
    def _voronoi_diagram(self):
        from flatsurf.geometry.voronoi import VoronoiDiagram
        return VoronoiDiagram(self._surface, self._geometry._centers, weight="radius_of_convergence")

    def _L2_consistency_voronoi_boundary(self, cell, boundary_segment, opposite_cell):
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
q
        The integrals `f_{\Re, n, m}`, `f_{\Im, n, m}`, `(f,g)_{\Re, n, m}`,
        and `(f,g)_{\Im, n, m}` are currently all computed numerically. We know
        of but have not implemented a better approach, yet.

        TESTS::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()
            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S, safety=0, singularities=True, centers=True)
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(S, prec=3, geometry=Ω._geometry)
            sage: V = C._voronoi_diagram()
            sage: centers = C._voronoi_diagram_centers()
            sage: segments = [segment for center in centers for boundary in V.cell(center) for segment in boundary.segments_with_uniform_root_branch()]

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
        return (self._L2_consistency_voronoi_boundary_f_overline_f(cell, boundary_segment)
                - 2 * self._L2_consistency_voronoi_boundary_Re_f_overline_g(cell, boundary_segment, opposite_cell)
                + self._L2_consistency_voronoi_boundary_g_overline_g(opposite_cell, boundary_segment))

    def _L2_consistency_voronoi_boundary_f_overline_f(self, cell, boundary_segment):
        r"""
        Return the value of `\int_γ f\overline{f}.
        """
        center = self._surface(cell.label(), cell.center())

        int = self.CellIntegrator(self, cell)

        κ = cell.root_branch(boundary_segment)
        α = cell.center()
        d = center.angle() - 1

        def f_(part, n, m):
            return int.mixed_integral(part, α, κ, d, n, α, κ, d, m, boundary_segment)

        return self.symbolic_ring().sum([
            (int.Re_a(n) * int.Re_a(m) + int.Im_a(n) * int.Im_a(m)) * f_("Re", n, m) + (int.Re_a(n) * int.Im_a(m) - int.Im_a(n) * int.Re_a(m)) * f_("Im", n, m)
            for n in range(self._prec * center.angle())
            for m in range(self._prec * center.angle())])

    def _L2_consistency_voronoi_boundary_Re_f_overline_g(self, cell, boundary_segment, opposite_cell):
        r"""
        Return the real part of `\int_γ f\overline{g}.
        """
        center = self._surface(cell.label(), cell.center())
        opposite_center = self._surface(opposite_cell.label(), opposite_cell.center())

        int = self.CellIntegrator(self, cell)
        jnt = self.CellIntegrator(self, opposite_cell)

        κ = cell.root_branch(boundary_segment)
        λ = opposite_cell.root_branch(boundary_segment)
        α = cell.center()
        β = opposite_cell.center()
        d = center.angle() - 1
        dd = opposite_center.angle() - 1

        def fg_(part, n, m):
            return int.mixed_integral(part, α, κ, d, n, β, λ, dd, m, boundary_segment)

        return self.symbolic_ring().sum([
            (int.Re_a(n) * jnt.Re_a(m) + int.Im_a(n) * jnt.Im_a(m)) * fg_("Re", n, m) + (int.Re_a(n) * jnt.Im_a(m) - int.Im_a(n) * jnt.Re_a(m)) * fg_("Im", n, m)
            for n in range(self._prec * center.angle())
            for m in range(self._prec * opposite_center.angle())])

    def _L2_consistency_voronoi_boundary_g_overline_g(self, opposite_cell, boundary_segment):
        r"""
        Return the value of `\int_γ g\overline{g}.
        """
        opposite_center = self._surface(opposite_cell.label(), opposite_cell.center())

        jnt = self.CellIntegrator(self, opposite_cell)

        λ = opposite_cell.root_branch(boundary_segment)
        β = opposite_cell.center()
        dd = opposite_center.angle() - 1

        def g_(part, n, m):
            return jnt.mixed_integral(part, β, λ, dd, n, β, λ, dd, m, boundary_segment)

        return self.symbolic_ring().sum([
            (jnt.Re_a(n) * jnt.Re_a(m) + jnt.Im_a(n) * jnt.Im_a(m)) * g_("Re", n, m) + (jnt.Re_a(n) * jnt.Im_a(m) - jnt.Im_a(n) * jnt.Re_a(m)) * g_("Im", n, m)
            for n in range(self._prec * opposite_center.angle())
            for m in range(self._prec * opposite_center.angle())])

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

    def variables(self):
        terms = list(self._constraints)
        if self._cost is not None:
            terms.append(self._cost)

        return set(variable for constraint in terms for variable in constraint.variables())

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
        for variable in self.variables():
            if variable.describe()[0] == "λ?":
                continue

            L = self._cost.derivative(variable)

            for i in range(lagranges):
                L += g[i][variable] * self.lagrange(i)

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
        return set(variable for variable in self.variables() if variable.describe()[0] == "λ?")

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
            sage: C._optimize_cost()
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
            assert constraint.total_degree() <= 1
            for variable in constraint.variables():
                gen, degree = variable.describe()
                if gen != "λ?":
                    non_lagranges.add(variable)

        non_lagranges = {variable: i for (i, variable) in enumerate(non_lagranges)}

        # TODO: Can we make this warning work again?
        # if len(set(self.symbolic_ring()._regular_gens(self._prec))) != len(non_lagranges):
        #     if not nowarn:
        #         from warnings import warn
        #         warn(f"Some power series coefficients are not constrained for this harmonic differential. They will be chosen to be 0 by the solver.")

        from sage.all import matrix, vector
        A = matrix(self.real_field(), len(self._constraints), len(non_lagranges) + len(self.lagrange_variables()))
        b = vector(self.real_field(), len(self._constraints))

        for row, constraint in enumerate(self._constraints):
            b[row] -= constraint.constant_coefficient()

            for variable in constraint.variables():
                if variable.describe()[0] == "λ?":
                    if variable not in lagranges:
                        lagranges[variable] = len(non_lagranges) + len(lagranges)
                    column = lagranges[variable]
                else:
                    column = non_lagranges[variable]

                A[row, column] += constraint[variable]

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
            Power Series Ring in z0 over Complex Field with 54 bits of precision
            sage: Ω.power_series_ring(T(0, 0))
            Power Series Ring in z1 over Complex Field with 54 bits of precision

        """
        from sage.all import PowerSeriesRing

        if len(args) != 1:
            raise NotImplementedError

        point = args[0]
        return PowerSeriesRing(self.complex_field(), f"z{self._geometry._centers.index(point)}")

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
            ({Point (1/2, 1/2) of polygon 0: 1.00000000000000 + 1.00000000000000*z0 + O(z0^2)},
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
            _cppyy().load_library("mpfr")

            solution = define_solve()(A, b)

            from sage.all import vector
            solution = vector([self.real_field()(entry) for entry in solution])
        else:
            raise NotImplementedError

        residue = (A*solution - b).norm()

        lagranges = len(self.lagrange_variables())

        if lagranges:
            solution = solution[:-lagranges]

        series = {point: {} for point in self._geometry._centers}

        for variable, column in decode.items():
            value = solution[column]

            gen, degree = variable.describe()

            if gen.startswith("Re(a"):
                part = 0
            elif gen.startswith("Im(a"):
                part = 1
            else:
                assert False

            center = self._geometry._centers[int(gen.split(",")[0][4:])]

            if degree not in series[center]:
                series[center][degree] = [None, None]

            series[center][degree][part] = value

        series = {point:
                  sum((self.complex_field()(*entry) * self.laurent_series_ring(point).gen()**k for (k, entry) in series[point].items()), start=self.laurent_series_ring(point).zero()).add_bigoh(max(series[point]) + 1) for point in series if series[point]}

        return series, residue


class Path:
    pass


class GeodesicPath(Path):
    def __init__(self, surface, start, end, holonomy):
        if holonomy.is_mutable():
            raise TypeError("holonomy must be immutable")
        if start[1].is_mutable():
            raise TypeError("start point must be immutable")
        if end[1].is_mutable():
            raise TypeError("end point must be immutable")

        self._surface = surface
        self._start = start
        self._end = end
        self._holonomy = holonomy

    @staticmethod
    def across_edge(surface, label, edge):
        r"""
        Return the :class:`GeodesicPath` that crosses from the center of the
        polygon ``label`` to the center of the polygon across the ``edge`` in
        a straight line.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: from flatsurf.geometry.harmonic_differentials import GeodesicPath
            sage: GeodesicPath.across_edge(S, 0, 0)
            Path (0, -a - 1) from (1/2, 1/2*a + 1/2) in polygon 0 to (1/2, 1/2*a + 1/2) in polygon 0

        """
        polygon = surface.polygon(label)

        if not polygon.is_convex():
            raise NotImplementedError

        opposite_label, opposite_edge = surface.opposite_edge(label, edge)
        opposite_polygon = surface.polygon(opposite_label)

        holonomy = polygon.vertex(edge) - polygon.centroid() + opposite_polygon.centroid() - opposite_polygon.vertex(opposite_edge + 1)
        holonomy.set_immutable()

        return GeodesicPath(surface, (label, polygon.centroid()), (opposite_label, opposite_polygon.centroid()), holonomy)

    def split(self):
        r"""
        Return the path as a sequence of segments in polygons.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: from flatsurf.geometry.harmonic_differentials import GeodesicPath
            sage: path = GeodesicPath.across_edge(S, 0, 0)
            sage: path.split()
            [(0, OrientedSegment((1/2, 1/2*a + 1/2), (1/2, 0))),
             (0, OrientedSegment((1/2, a + 1), (1/2, 1/2*a + 1/2)))]

        """
        polygon = self._surface.polygon(self._start[0])
        if not polygon.is_convex():
            raise NotImplementedError

        from flatsurf.geometry.euclidean import OrientedSegment
        if self._end[0] == self._start[0] and self._start[1] + self._holonomy == self._end[1]:
            return [(self._start[0], OrientedSegment(self._start[1], self._end[1]))]

        path = OrientedSegment(self._start[1], self._start[1] + self._holonomy)
        for v in range(len(polygon.vertices())):
            edge = OrientedSegment(polygon.vertex(v), polygon.vertex(v + 1))
            intersection = edge.intersection(path)

            if intersection is None:
                continue

            if intersection == self._start[1]:
                continue

            assert self._surface(self._start[0], intersection).angle() == 1, "path crosses over a singularity"

            segment = OrientedSegment(self._start[1], intersection)

            prefix = [(self._start[0], segment)]

            holonomy = self._holonomy - (intersection - self._start[1])
            holonomy.set_immutable()

            if not holonomy:
                return prefix

            from flatsurf.geometry.euclidean import is_parallel
            assert is_parallel(holonomy, self._holonomy)

            opposite_label, opposite_edge = self._surface.opposite_edge(self._start[0], v)

            new_start = self._surface.polygon(opposite_label).vertex(opposite_edge + 1) + (intersection - polygon.vertex(v))
            new_start.set_immutable()

            return prefix + GeodesicPath(self._surface, (opposite_label, new_start), self._end, holonomy).split()

        assert False

    def __neg__(self):
        holonomy = -self._holonomy
        holonomy.set_immutable()
        return GeodesicPath(self._surface, self._end, self._start, holonomy)

    def __repr__(self):
        return f"Path {self._holonomy} from {self._start[1]} in polygon {self._start[0]} to {self._end[1]} in polygon {self._end[0]}"

    def __eq__(self, other):
        if not isinstance(other, GeodesicPath):
            return False
        if self._start == other._start and self._holonomy == other._holonomy:
            assert self._end == other._end
            return True

        return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self._start, self._holonomy))
