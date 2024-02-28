r"""
TODO: Document this module.
TODO: Rename this module?

EXAMPLES:

We compute harmonic differentials on the square torus::

    sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialCohomology
    sage: T = translation_surfaces.torus((1, 0), (0, 1))
    sage: T.set_immutable()

    sage: H = SimplicialCohomology(T)
    sage: a, b = H.homology().gens()

First, the harmonic differentials that sends the horizontal `a` to 1 and the
vertical `b` to zero::

    sage: f = H({b: 1})
    sage: Ω = HarmonicDifferentials(T)
    sage: ω = Ω(f)
    sage: ω  # random output
    ((1.00000000000000 + 1.48374000000000e-76*I) + (-1.35679000000000e-77 + 1.67690000000000e-77*I)*z0 + (5.60404000000000e-76 - 3.17779000000000e-76*I)*z0^2 + (-1.38688000000000e-76 + 2.40493000000000e-75*I)*z0^3 + (-1.72352000000000e-74 - 1.88541000000000e-74*I)*z0^4 + O(z0^5), (1.00000000000000 + 2.85447000000000e-77*I) + (4.33241000000000e-78 - 7.66907000000000e-78*I)*z1 + (-3.82322000000000e-77 + 8.71445000000000e-77*I)*z1^2 + (-7.69949000000000e-77 + 6.04329000000000e-76*I)*z1^3 + (4.42836000000000e-75 + 5.93165000000000e-76*I)*z1^4 + O(z1^5))
    sage: ω.simplify()
    (1.00000000000000 + O(z0^5), 1.00000000000000 + O(z1^5))

The harmonic differential that integrates as 0 along `a` but 1 along `b`::

    sage: g = H({a: -1})
    sage: Ω(g).simplify()
    (1.00000000000000*I + O(z0^5), 1.00000000000000*I + O(z1^5))

A less trivial example, the regular octagon::

    sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialCohomology
    sage: S = translation_surfaces.regular_octagon()

    sage: H = SimplicialCohomology(S)
    sage: a, b, c, d = H.homology().gens()

    sage: f = H({ a: sqrt(2) + 1, b: 0, c: -sqrt(2) - 1, d: -sqrt(2) - 2})

    sage: Omega = HarmonicDifferentials(S, ncoefficients=11)
    sage: omega = Omega(f)
    sage: omega.simplify(zero_threshold=1e-4)  # abs-tol 1e-4  # TODO: Why so much tolerance?
    (-2.09019000000000 + (-3.22546000000000)*z0^8 + (-2.61535000000000)*z0^16 + (-2.00435000000000)*z0^24 + (-1.53401000000000)*z0^32 + O(z0^33), 1.31413000000000*z1^2 + (-0.107411000000000)*z1^10 + O(z1^11))

The same computation on a triangulation of the octagon::

    sage: from flatsurf import HarmonicDifferentials, SimplicialCohomology, Polygon, translation_surfaces
    sage: S = translation_surfaces.regular_octagon()
    sage: S = S.subdivide().codomain()

    sage: H = SimplicialCohomology(S)
    sage: a, b, c, d = H.homology().gens()

    sage: f = H({a: -sqrt(2), b: 0, c: -sqrt(2) - 1, d: sqrt(2) + 1})

    sage: Omega = HarmonicDifferentials(S, centers="vertices", ncoefficients=11)
    sage: omega = Omega(f)  # long time
    sage: omega.simplify(zero_threshold=1e-4)  # abs-tol 1e-4  # long time, see above  # TODO: Why so much tolerance?
    (1.31414000000000*z0^2 + (-0.107413000000000)*z0^10 + O(z0^11), -2.09020000000000 + (-3.22547000000000)*z1^8 + (-2.61537000000000)*z1^16 + (-2.00436000000000)*z1^24 + (-1.53401000000000)*z1^32 + O(z1^33))

The same surface but built as the unfolding of a right triangle::

    sage: from flatsurf import similarity_surfaces, HarmonicDifferentials, SimplicialCohomology, Polygon
    sage: S = similarity_surfaces.billiard(Polygon(angles=[3/8, 1/2, 1/8], lengths=[1/2])).minimal_cover('translation')

    sage: H = SimplicialCohomology(S)
    sage: a, b, c, d = H.homology().gens()

    sage: f = H({a: 0, b: sqrt(2) + 2, c: -1, d: -sqrt(2) - 1})

    sage: Omega = HarmonicDifferentials(S, centers="vertices+centers", ncoefficients=11)  # TODO: With just "vertices" this does not terminate. Probably because there are areas in the Voronoi diagram not contained in any cell.
    sage: omega = Omega(f)  # long time  # random output due to precision warnings TODO
    sage: omega.simplify(zero_threshold=1e-4)  # abs-tol 1e-4  # long time, see above  # TODO: We get this output but multiplied with a root of unity.
    (-2.09019000000000 + (-3.22546000000000)*z0^8 + (-2.61535000000000)*z0^16 + (-2.00435000000000)*z0^24 + (-1.53401000000000)*z0^32 + O(z0^33), 1.31413000000000*z1^2 + (-0.107411000000000)*z1^10 + O(z1^11))

Much more complicated, the unfolding of the (3, 4, 13) triangle::

    sage: from flatsurf import similarity_surfaces, SimplicialCohomology, HarmonicDifferentials, Polygon

    sage: S = similarity_surfaces.billiard(Polygon(angles=[3, 4, 13])).minimal_cover("translation")
    sage: S = S.erase_marked_points().codomain().delaunay_triangulation()

    sage: H = SimplicialCohomology(S)
    sage: f = H({H.homology().gens()[0]: 1})
    sage: g = H({H.homology().gens()[1]: 1})
    sage: fg = H({H.homology().gens()[1]: 1, H.homology().gens()[0]: 1})

    sage: Omega = HarmonicDifferentials(S, error=1e-1, centers="vertices")
    sage: omegaf = Omega(f, check=False)  # TODO: Increase precision once this is faster.  # long time
    sage: omegag = Omega(g, check=False)
    sage: omegafg = Omega(fg, check=False)
    sage: omegaf + omegag - omegafg
    sage: omega.simplify()  # long time, see above
    ((0.43474 + 0.32379*I) + O(z0), (0.014497 + 0.099627*I) + O(z1), (-0.050347 + 0.0029207*I) + O(z2), (-0.21127 - 0.090282*I) + O(z3), (0.47454 + 0.087358*I) + O(z4), (-0.19314 + 0.31055*I) + O(z5), (-0.23193 - 0.13398*I) + O(z6), (-0.45636 - 0.42553*I) + O(z7), (-0.21176 + 0.40428*I) + O(z8), (-0.15358 + 0.049085*I) + O(z9), (-0.40648 - 0.48943*I) + O(z10), (-0.0020610 - 0.48399*I) + O(z11), (0.24142 - 0.25153*I) + O(z12), (-0.032511 + 0.55139*I) + O(z13), (-0.32812 + 0.13133*I) + O(z14), (-0.33572 + 0.53022*I) + O(z15), (0.079225 + 0.028393*I) + O(z16), (0.033120 - 0.49532*I) + O(z17), (-0.27160 - 0.34612*I) + O(z18), (-0.21121 + 0.25608*I) + (-0.13894 - 0.56381*I)*z19 + (0.040253 + 0.39319*I)*z19^2 + (-0.088971 - 0.70043*I)*z19^3 + (0.075621 - 0.040790*I)*z19^4 + (-0.76242 - 0.39754*I)*z19^5 + (-0.066434 - 0.12278*I)*z19^6 + (0.0075069 + 0.22176*I)*z19^7 + (-0.40163 - 0.51249*I)*z19^8 + (0.55621 + 0.080413*I)*z19^9 + (-1.0088 - 1.9733*I)*z19^10 + (0.27571 - 0.061058*I)*z19^11 + (-2.4352 - 1.0704*I)*z19^12 + O(z19^13), (-0.14288 + 0.36838*I) + O(z20), (-0.012062 - 0.58947*I) + O(z21), (0.44529 + 0.29461*I) + (0.22745 - 0.67464*I)*z22 + (-0.95440 - 0.53814*I)*z22^2 + O(z22^3), (-2.4328 - 0.60098*I) + O(z23), (-0.71168 - 0.73386*I) + O(z24), (-1.6115 - 1.7234*I) + O(z25))

"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2024 Julian Rüth
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
        import os.path
        _cppyy().include(os.path.join(os.path.dirname(__file__), "..", "..", "mpreal-support.h"))
        _cppyy().cppdef(r'''
        #include <cassert>
        #include <vector>
        #include <iostream>
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

            sage: ω = Ω(f) + Ω(f) + Ω(H({a: -2}))
            sage: ω  # random output due to numerical noise
            ((-2.00000000007046e-82*I) + (-9.99999999996877e-83*I)*z0 + (2.00000000003978e-81 + 3.99999999998751e-82*I)*z0^2 + (3.99999999995683e-81)*z0^3 + (-1.00000000001253e-79*I)*z0^4 + O(z0^5), (3.99999999996833e-83*I) + (-4.00000000006421e-83)*z1 + 0.000000000000000*z1^2 + (-1.99999999997841e-81*I)*z1^3 + (9.99999999983070e-81*I)*z1^4 + O(z1^5))
            sage: ω.simplify()
            (O(z0^5), O(z1^5))

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

    # TODO: Can we increase the default precision? Somehow, we do not get much more precision easily in the octagon. Why?
    def error(self, kind=None, verbose=False, abs_tol=1e-4, rel_tol=1e-4):
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

        if kind is None or "L2" in kind:
            C = PowerSeriesConstraints(self.parent())
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
            sage: f = H({b: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

            sage: η.series(T(0, (1/2, 1/2)))  # abstol 1e-9
            1.00000000000000 + 9.80040000000000e-77*I + (-4.20905000000000e-77 - 8.26568000000000e-79*I)*z0 + (-1.05325000000000e-75 + 6.00322000000000e-78*I)*z0^2 + (-2.11385000000000e-75 - 1.14645000000000e-75*I)*z0^3 + (6.53715000000000e-75 - 1.52468000000000e-74*I)*z0^4 + O(z0^5)

        """
        return self._series[triangle]

    @cached_method
    def _constraints(self):
        # TODO: This is a hack. Come up with a better class hierarchy!
        return PowerSeriesConstraints(differentials=self.parent())

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
            sage: f = H({b: 1})

            sage: Ω = HarmonicDifferentials(T)
            sage: η = Ω(f)

        Compute the constant coefficients::

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(Ω)
            sage: R = C.symbolic_ring()
            sage: gen = C._gen("Re", T(0, (1/2, 1/2)), 0)
            sage: η._evaluate(gen)
            1.00000000000000

        """
        coefficients = {}

        for variable in expression.variables():
            center = self.parent()._gen_center(variable)
            degree = self.parent()._gen_degree(variable)

            try:
                coefficient = self._series[center][degree]
            except (IndexError, KeyError):
                import warnings
                warnings.warn(f"expected a {degree}th coefficient of the power series around {center} but none found")
                coefficients[variable] = 0
                continue

            if self.parent()._gen_is_real(variable):
                coefficients[variable] = coefficient.real()
            elif self.parent()._gen_is_imag(variable):
                coefficients[variable] = coefficient.imag()
            else:
                raise NotImplementedError

        value = expression(coefficients)

        from flatsurf.geometry.power_series import PowerSeriesCoefficientExpression
        if isinstance(value, PowerSeriesCoefficientExpression):
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

        C = PowerSeriesConstraints(self.parent())
        return self._evaluate(C.integrate(cycle))

    def _repr_(self):
        def sparse_parent(series):
            return type(series.parent())(series.parent().base_ring(), series.parent().variable_name(), sparse=True)

        return repr(tuple(sparse_parent(s)(s) for s in self._series.values()))

    def simplify(self, zero_threshold=1e-6):
        def simplify(series):
            def simplify(coefficient):
                if coefficient.real().abs() < zero_threshold:
                    coefficient = coefficient.parent()(coefficient.imag() * coefficient.parent()('I'))
                if coefficient.imag().abs() < zero_threshold:
                    coefficient = coefficient.parent()(coefficient.real())

                return coefficient

            coefficients = {exponent: simplify(coefficient) for (exponent, coefficient) in zip(series.exponents(), series.coefficients())}
            coefficients = {exponent: coefficient for (exponent, coefficient) in coefficients.items() if coefficient}

            return series.parent()(coefficients, prec=series.prec())

        return type(self)(self.parent(), {center: simplify(series) for (center, series) in self._series.items()}, residue=self._residue, cocycle=self._cocycle)

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

        return plot


# TODO: Make these unique for each surface (without using UniqueRepresentation because equal surfaces can be distinct.)
class HarmonicDifferentials(Parent):
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
        (O(z0^5), O(z1^5))

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

    # TODO: Determine ncoefficients automatically
    def __init__(self, surface, error=1e-3, centers=None, category=None):
        Parent.__init__(self, category=category or SetsWithPartialMaps())

        try:
            sorted(surface.labels())
        except Exception:
            raise NotImplementedError("labels on the surface must be sortable so we use label order to make a choice of n-th roots")

        self._surface = surface
        self._error = error
        self._centers = tuple(HarmonicDifferentials._centers(surface, centers))

        # TODO: Why does calling this here fix caching issues in ncoefficients()??
        self.error()

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

    @cached_method
    def ncoefficients(self, center):
        beta = self.beta(self._voronoi_diagram().cell(center))
        if beta >= 1:
            raise ValueError(f"cell at {center} extends beyond radius of convergence")

        from math import log, ceil
        ncoefficients = log(self._error, beta)

        ncoefficients /= center.angle()

        ncoefficients = int(ceil(ncoefficients))

        if ncoefficients <= 0:
            ncoefficients = 1

        # TODO: Should we still multiply here?
        ncoefficients *= center.angle()

        print(f"{ncoefficients} coefficients at {center} with β={beta}")

        return ncoefficients

    def surface(self):
        return self._surface

    @cached_method
    def _voronoi_diagram(self):
        from flatsurf.geometry.voronoi import VoronoiDiagram
        return VoronoiDiagram(self._surface, self._centers, weight="radius_of_convergence")

    @cached_method
    def error(self, cell=None):
        # Returns the a-priori (is it?) error for the value of the differential
        # anywhere in cell (or anywhere in all cells); relative to the maximum
        # of the differential. TODO: Explain algorithm from my notes.
        if cell is None:
            return max(self.error(cell=cell) for cell in self._voronoi_diagram().cells())

        β = self.beta(cell=cell)
        if β >= 1:
            from sage.all import oo
            return oo

        return abs((β**(self.ncoefficients(cell._center) + 1)) / (1 - β))

    def beta(self, cell):
        from math import sqrt
        R = sqrt(float(cell.radius_of_convergence()))
        r = sqrt(float(cell.radius()))

        return r / R

    def error_location(self, cell=None):
        if cell is None:
            cell = self.error_cell()

        return cell.furthest_point()

    def error_cell(self):
        return max(self._voronoi_diagram().cells(), key=lambda cell: self.error(cell=cell))

    def _repr_(self):
        return f"Ω({self._surface})"

    # TODO: Move to some class that is shared between harmonic differentials
    # and constraints that abstracts away details of the symbolic ring.
    def _gen_center(self, gen):
        assert self._gen_is_real(gen) or self._gen_is_imag(gen)

        gen, degree = gen.describe()

        return self._centers[int(gen.split(",")[0][4:])]

    # TODO: Move to some class that is shared between harmonic differentials
    # and constraints that abstracts away details of the symbolic ring.
    def _gen_is_lagrange(self, gen):
        gen, degree = gen.describe()

        return gen == "λ?"

    # TODO: Move to some class that is shared between harmonic differentials
    # and constraints that abstracts away details of the symbolic ring.
    def _gen_is_real(self, gen):
        gen, degree = gen.describe()

        return gen.startswith("Re(a")

    # TODO: Move to some class that is shared between harmonic differentials
    # and constraints that abstracts away details of the symbolic ring.
    def _gen_is_imag(self, gen):
        gen, degree = gen.describe()

        return gen.startswith("Im(a")

    # TODO: Move to some class that is shared between harmonic differentials
    # and constraints that abstracts away details of the symbolic ring.
    def _gen_degree(self, gen):
        gen, degree = gen.describe()
        return degree

    @cached_method(key=lambda self, check: None, do_pickle=True)
    def basis(self, check=True):
        from flatsurf.geometry.homology import SimplicialHomology
        H = SimplicialHomology(self._surface)
        from flatsurf.geometry.cohomology import SimplicialCohomology
        HH = SimplicialCohomology(self._surface)
        from sage.all import ZZ
        return [self(HH({gen: ZZ.one()}), check=check) for gen in H.gens()]

    @cached_method(key=lambda self, check: None, do_pickle=True)
    def period_matrix(self, check=True):
        from flatsurf.geometry.homology import SimplicialHomology
        symplectic_basis = SimplicialHomology(self._surface).symplectic_basis()
        symplectic_basis = symplectic_basis[:len(symplectic_basis)//2]

        from sage.all import matrix
        return matrix([[differential.integrate(path) for differential in self.basis(check=check)] for path in symplectic_basis])

    def _element_constructor_(self, x, *args, **kwargs):
        if not x:
            return self.element_class(self, None, *args, **kwargs)

        if isinstance(x, dict):
            return self.element_class(self, x, *args, **kwargs)

        return self._element_from_cohomology(x, *args, **kwargs)

    def _element_from_cohomology(self, cocycle, /, algorithm=["L2"], check=True):
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

        constraints = PowerSeriesConstraints(self)

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
            derivatives = get_parameter("midpoint_derivatives", self.prec//3)
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


class PowerSeriesConstraints:
    r"""
    A collection of (linear) constraints on the coefficients of power series
    developed at the vertices of the Voronoi cells of a Delaunay triangulation.

    This is used to create harmonic differentials from cohomology classes.
    """

    def __init__(self, differentials):
        self._differentials = differentials
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

            sage: C = PowerSeriesConstraints(Ω)
            sage: C.symbolic_ring()
            Ring of Power Series Coefficients in Re(a0,0),…,Re(a1,0),…,Im(a0,0),…,Im(a1,0),…,λ0,… over Complex Field with 54 bits of precision

        """
        # TODO: What's the correct precision here?

        gens = [f"Re(a{n},?)" for n in range(len(self._differentials._centers))] + [f"Im(a{n},?)" for n in range(len(self._differentials._centers))] + ["λ?"]

        from sage.all import ComplexField
        from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
        return PowerSeriesCoefficientExpressionRing(base_ring or ComplexField(54), tuple(gens))

    @cached_method
    def complex_field(self):
        from sage.all import ComplexField
        # TODO: Why 54?
        return ComplexField(54)

    @cached_method
    def real_field(self):
        from sage.all import RealField
        # TODO: Why 54?
        return RealField(54)

    # TODO: Move to Harmonic Differentials; or maybe some other shared class
    # that abstracts away details of the symbolic ring.
    @cached_method
    def lagrange(self, k):
        return self.symbolic_ring().gen(("λ?", k))

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
            sage: Ω = HarmonicDifferentials(T, centers="vertices", ncoefficients=3)
            sage: C = PowerSeriesConstraints(Ω)

            sage: C.integrate(H())
            0

        Integrating the power series developed around the vertex along the path
        that loops horizontally from the center of the square to itself, we get
        `a_0 - i/2 a_1 - a_2/6`::  # TODO: This is not true anymore.

            sage: a, b = H.gens()
            sage: C.integrate(b)  # TODO: There are many correct answers here. Test for something meaningful.
            Re(a0,0) ...

        :: # TODO: Explain what's the expected output here

            sage: C.integrate(-a)  # TODO: There are many correct answers here. Test for something meaningful.
            (-1.00000000000000*I)*Re(a0,0) ...

        The same integrals but developing the power series at the vertex and at the center of the square::

            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(Ω)

            sage: C.integrate(a)   # not tested # TODO: Check these values
            0.828427124746190*Re(a0,0) + 0.171572875253810*Re(a1,0) + 0.828427124746190*I*Im(a0,0) + 0.171572875253810*I*Im(a1,0) + (-1.38777878078145e-17 - 8.54260702578763e-18*I)*Re(a0,1) + (-3.03576608295941e-18 + 0.0857864376269050*I)*Re(a1,1) + (8.54260702578763e-18 - 1.38777878078145e-17*I)*Im(a0,1) + (-0.0857864376269050 - 3.03576608295941e-18*I)*Im(a1,1) + (0.0473785412436502 + 4.71795158413990e-18*I)*Re(a0,2) + (-0.0424723326565069 - 3.03576608295941e-18*I)*Re(a1,2) + (-4.71795158413990e-18 + 0.0473785412436502*I)*Im(a0,2) + (3.03576608295941e-18 - 0.0424723326565069*I)*Im(a1,2) + (-1.73472347597681e-18 - 2.19851947436667e-18*I)*Re(a0,3) + (2.16840434497101e-18 - 0.0208152801713079*I)*Re(a1,3) + (2.19851947436667e-18 - 1.73472347597681e-18*I)*Im(a0,3) + (0.0208152801713079 + 2.16840434497101e-18*I)*Im(a1,3) + (0.00487732352790257 + 9.71367022318980e-19*I)*Re(a0,4) + (0.0100938339276945 + 1.51788304147971e-18*I)*Re(a1,4) + (-9.71367022318980e-19 + 0.00487732352790257*I)*Im(a0,4) + (-1.51788304147971e-18 + 0.0100938339276945*I)*Im(a1,4)
            sage: C.integrate(b)  # not tested # TODO: Check these values
            (-0.828427124746190*I)*Re(a0,0) + (-0.171572875253810*I)*Re(a1,0) + 0.828427124746190*Im(a0,0) + 0.171572875253810*Im(a1,0) + (-1.38777878078145e-17 + 8.54260702578763e-18*I)*Re(a0,1) + (-3.03576608295941e-18 - 0.0857864376269050*I)*Re(a1,1) + (-8.54260702578763e-18 - 1.38777878078145e-17*I)*Im(a0,1) + (0.0857864376269050 - 3.03576608295941e-18*I)*Im(a1,1) + (7.70371977754894e-34 + 0.0473785412436502*I)*Re(a0,2) + (-2.60208521396521e-18 - 0.0424723326565069*I)*Re(a1,2) + (-0.0473785412436502 + 7.70371977754894e-34*I)*Im(a0,2) + (0.0424723326565069 - 2.60208521396521e-18*I)*Im(a1,2) + (1.73472347597681e-18 - 2.19851947436667e-18*I)*Re(a0,3) + (-2.16840434497101e-18 - 0.0208152801713079*I)*Re(a1,3) + (2.19851947436667e-18 + 1.73472347597681e-18*I)*Im(a0,3) + (0.0208152801713079 - 2.16840434497101e-18*I)*Im(a1,3) + (-9.62964972193618e-35 - 0.00487732352790257*I)*Re(a0,4) + (-1.08420217248550e-18 - 0.0100938339276945*I)*Re(a1,4) + (0.00487732352790257 - 9.62964972193618e-35*I)*Im(a0,4) + (0.0100938339276945 - 1.08420217248550e-18*I)*Im(a1,4)

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology
            sage: S = translation_surfaces.regular_octagon()

            sage: H = SimplicialHomology(S)
            sage: Ω = HarmonicDifferentials(S)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(Ω)

            sage: C.integrate(H())
            0

            sage: a, b, c, d = H.gens()

            sage: C.integrate(a)  # not tested # TODO: Check these values
            (1.13290899470104 - 1.13290899470104*I)*Re(a0,0) + (1.13290899470104 + 1.13290899470104*I)*Im(a0,0) + (0.275277280554704 + 0.275277280554704*I)*Re(a1,0) + (-0.275277280554704 + 0.275277280554704*I)*Im(a1,0) + (0.327551243899061 + 6.93889390390723e-18*I)*Re(a1,1) + (-6.93889390390723e-18 + 0.327551243899061*I)*Im(a1,1) + (0.191399262161835 - 0.191399262161835*I)*Re(a1,2) + (0.191399262161835 + 0.191399262161835*I)*Im(a1,2)

            sage: C.integrate(b)  # not tested # TODO: Check these values
            1.60217526524068*Re(a0,0) + 1.60217526524068*I*Im(a0,0) + (-0.389300863573646 + 1.38777878078145e-17*I)*Re(a1,0) + (-1.38777878078145e-17 - 0.389300863573646*I)*Im(a1,0) + (-1.04083408558608e-17 - 0.327551243899061*I)*Re(a1,1) + (0.327551243899061 - 1.04083408558608e-17*I)*Im(a1,1) + 0.270679432377470*Re(a1,2) + 0.270679432377470*I*Im(a1,2)

            sage: C.integrate(c)  # not tested # TODO: Check these values
            1.60217526524068*Re(a0,0) + 1.60217526524068*I*Im(a0,0) + (-0.389300863573646 + 1.38777878078145e-17*I)*Re(a1,0) + (-1.38777878078145e-17 - 0.389300863573646*I)*Im(a1,0) + (-1.04083408558608e-17 - 0.327551243899061*I)*Re(a1,1) + (0.327551243899061 - 1.04083408558608e-17*I)*Im(a1,1) + 0.270679432377470*Re(a1,2) + 0.270679432377470*I*Im(a1,2)

            sage: C.integrate(d)  # not tested # TODO: Check these values
            (-1.60217526524068*I)*Re(a0,0) + 1.60217526524068*Im(a0,0) + (5.55111512312578e-17 - 0.389300863573646*I)*Re(a1,0) + (0.389300863573646 + 5.55111512312578e-17*I)*Im(a1,0) + (-2.77555756156289e-17 + 0.327551243899061*I)*Re(a1,1) + (-0.327551243899061 - 2.77555756156289e-17*I)*Im(a1,1) + (-0.270679432377470*I)*Re(a1,2) + 0.270679432377470*Im(a1,2)

        """
        return sum((multiplicity * sgn * self._integrate_path(path) for (label, edge), multiplicity in cycle._chain.monomial_coefficients().items() for (sgn, path) in self._integrate_path_along_edge(label, edge)), start=self.symbolic_ring().zero())

    def _integrate_path_along_edge(self, label, edge):
        r"""
        Return the path along the ``edge`` of the polygon with ``label`` as a
        sequence of :class:`Path` that we can integrate along.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces, SimplicialHomology
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: H = SimplicialHomology(T)

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(Ω)
            sage: C._integrate_path_along_edge(0, 1)
            [(1, Path (0, 1) from (1, 0) in polygon 0 to (1, 1) in polygon 0)]

        ::

            sage: from flatsurf import translation_surfaces, HarmonicDifferentials
            sage: S = translation_surfaces.regular_octagon()
            sage: Ω = HarmonicDifferentials(S)
            sage: PowerSeriesConstraints(Ω)._integrate_path_along_edge(0, 0)
            [(1, Path (1, 0) from (0, 0) in polygon 0 to (1, 0) in polygon 0)]

        """
        return [(1, GeodesicPath.along_edge(self._differentials.surface(), label, edge))]

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
            sage: C = PowerSeriesConstraints(Ω)

            sage: from flatsurf.geometry.harmonic_differentials import GeodesicPath
            sage: path = GeodesicPath.along_edge(S, 0, 0)

            sage: C._integrate_path(path)  # not tested # TODO: Check this value
            (-5.55111512312578e-17 - 0.389300863573646*I)*Re(a0,0) + (-1.60217526524068*I)*Re(a1,0) + (0.389300863573646 - 5.55111512312578e-17*I)*Im(a0,0) + 1.60217526524068*Im(a1,0) + (-2.08166817117217e-17 - 0.327551243899060*I)*Re(a0,1) + (1.66533453693773e-16 + 3.19522800018857e-17*I)*Re(a1,1) + (0.327551243899060 - 2.08166817117217e-17*I)*Im(a0,1) + (-3.19522800018857e-17 + 1.66533453693773e-16*I)*Im(a1,1) + (-0.270679432377470*I)*Re(a0,2) + (-1.23259516440783e-32 + 0.342727396656658*I)*Re(a1,2) + 0.270679432377470*Im(a0,2) + (-0.342727396656658 - 1.23259516440783e-32*I)*Im(a1,2) + (1.38777878078145e-17 - 0.219471472765136*I)*Re(a0,3) + (-1.24900090270330e-16 - 3.07576511193400e-17*I)*Re(a1,3) + (0.219471472765136 + 1.38777878078145e-17*I)*Im(a0,3) + (3.07576511193400e-17 - 1.24900090270330e-16*I)*Im(a1,3) + (2.42861286636753e-17 - 0.174329399573979*I)*Re(a0,4) + (1.69481835106077e-32 - 0.131965414609324*I)*Re(a1,4) + (0.174329399573979 + 2.42861286636753e-17*I)*Im(a0,4) + (0.131965414609324 + 1.69481835106077e-32*I)*Im(a1,4) + (3.12250225675825e-17 - 0.135339716188735*I)*Re(a0,5) + (0.135339716188735 + 3.12250225675825e-17*I)*Im(a0,5) + (2.77555756156289e-17 - 0.102340546397005*I)*Re(a0,6) + (0.102340546397005 + 2.77555756156289e-17*I)*Im(a0,6) + (3.46944695195361e-17 - 0.0749845895064374*I)*Re(a0,7) + (0.0749845895064374 + 3.46944695195361e-17*I)*Im(a0,7) + (3.12250225675825e-17 - 0.0527958835241931*I)*Re(a0,8) + (0.0527958835241931 + 3.12250225675825e-17*I)*Im(a0,8) + (2.77555756156289e-17 - 0.0352191503025193*I)*Re(a0,9) + (0.0352191503025193 + 2.77555756156289e-17*I)*Im(a0,9) + (2.08166817117217e-17 - 0.0216611462547145*I)*Re(a0,10) + (0.0216611462547145 + 2.08166817117217e-17*I)*Im(a0,10) + (2.08166817117217e-17 - 0.0115239671919220*I)*Re(a0,11) + (0.0115239671919220 + 2.08166817117217e-17*I)*Im(a0,11) + (1.73472347597681e-17 - 0.00423065874117904*I)*Re(a0,12) + (0.00423065874117904 + 1.73472347597681e-17*I)*Im(a0,12) + (1.04083408558608e-17 + 0.000756226861613830*I)*Re(a0,13) + (-0.000756226861613830 + 1.04083408558608e-17*I)*Im(a0,13) + (8.67361737988404e-18 + 0.00392229868304028*I)*Re(a0,14) + (-0.00392229868304028 + 8.67361737988404e-18*I)*Im(a0,14)

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
            sage: C = PowerSeriesConstraints(Ω)

            sage: from flatsurf.geometry.euclidean import OrientedSegment

            sage: C._integrate_path_polygon(0, OrientedSegment((0, 0), (1, 0)))  # TODO: Check this value
            (1.58740105196839 - 2.64762039020726e-18*I)*Re(a0,0) + (2.64762039020726e-18 + 1.58740105196839*I)*Im(a0,0) + (-1.44328993201270e-15 - 1.11813213514231e-18*I)*Re(a0,1) + (1.11813213514231e-18 - 1.44328993201270e-15*I)*Im(a0,1) + 0.333333333333333*Re(a0,2) + 0.333333333333333*I*Im(a0,2) + (-1.09614712307247e-19*I)*Re(a0,3) + (1.09614712307247e-19)*Im(a0,3) + (0.125992104989487 + 1.78541206316279e-19*I)*Re(a0,4) + (-1.78541206316279e-19 + 0.125992104989487*I)*Im(a0,4) + (5.10269499644730e-18*I)*Re(a0,5) + (-5.10269499644730e-18)*Im(a0,5) + (0.0566928947041920 + 7.90191292835058e-20*I)*Re(a0,6) + (-7.90191292835058e-20 + 0.0566928947041920*I)*Im(a0,6) + (-1.44886888580502e-18*I)*Re(a0,7) + (1.44886888580502e-18)*Im(a0,7) + (0.0277777777777778 - 3.40179666429820e-18*I)*Re(a0,8) + (3.40179666429820e-18 + 0.0277777777777778*I)*Im(a0,8) + (1.73472347597681e-18 + 1.25548457904919e-18*I)*Re(a0,9) + (-1.25548457904919e-18 + 1.73472347597681e-18*I)*Im(a0,9) + (0.0143172846575932 - 2.51141927450074e-19*I)*Re(a0,10) + (2.51141927450074e-19 + 0.0143172846575932*I)*Im(a0,10) + (8.67361737988404e-19 + 1.91351062366774e-18*I)*Re(a0,11) + (-1.91351062366774e-18 + 8.67361737988404e-19*I)*Im(a0,11) + (0.00763173582678622 + 7.56708582916301e-19*I)*Re(a0,12) + (-7.56708582916301e-19 + 0.00763173582678622*I)*Im(a0,12) + (8.67361737988404e-19 + 4.33707190579305e-19*I)*Re(a0,13) + (-4.33707190579305e-19 + 8.67361737988404e-19*I)*Im(a0,13) + (0.00416666666666667 - 1.02053899928946e-18*I)*Re(a0,14) + (1.02053899928946e-18 + 0.00416666666666667*I)*Im(a0,14)

        """
        V = self._differentials._voronoi_diagram()

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
            sage: C = PowerSeriesConstraints(Ω)

            sage: from flatsurf.geometry.euclidean import OrientedSegment

            sage: V = Ω._voronoi_diagram()
            sage: C._integrate_path_cell(V.polygon_cell(0, (0, 0)), OrientedSegment((0, 0), (1/2, 0)))  # TODO: Check this value
            0.793700525984099*Re(a0,0) + 0.793700525984099*I*Im(a0,0) + 0.314980262473718*Re(a0,1) + 0.314980262473718*I*Im(a0,1) + 0.166666666666667*Re(a0,2) + 0.166666666666667*I*Im(a0,2) + 0.0992125657480124*Re(a0,3) + 0.0992125657480124*I*Im(a0,3) + 0.0629960524947437*Re(a0,4) + 0.0629960524947437*I*Im(a0,4) + 0.0416666666666667*Re(a0,5) + 0.0416666666666667*I*Im(a0,5) + 0.0283464473520960*Re(a0,6) + 0.0283464473520960*I*Im(a0,6) + 0.0196862663993065*Re(a0,7) + 0.0196862663993065*I*Im(a0,7) + 0.0138888888888889*Re(a0,8) + 0.0138888888888889*I*Im(a0,8) + 0.00992125657430033*Re(a0,9) + 0.00992125657430033*I*Im(a0,9) + 0.00715864232879659*Re(a0,10) + 0.00715864232879659*I*Im(a0,10) + 0.00520833333333333*Re(a0,11) + 0.00520833333333333*I*Im(a0,11) + 0.00381586791339311*Re(a0,12) + 0.00381586791339311*I*Im(a0,12) + 0.00281232377208854*Re(a0,13) + 0.00281232377208854*I*Im(a0,13) + 0.00208333333333333*Re(a0,14) + 0.00208333333333333*I*Im(a0,14)

        """
        integrator = self.CellIntegrator(self, polygon_cell)
        ncoefficients = self._differentials.ncoefficients(self._differentials.surface()(polygon_cell.label(), polygon_cell.center()))
        return sum(integrator.a(n) * integrator.integral(n, segment) for n in range(ncoefficients))

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
            sage: consistency = PowerSeriesConstraints(Ω)._L2_consistency()
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

        costs = list(cost for _, cost in L2_cost(list(self._differentials._voronoi_diagram().boundaries().items())))

        return sum(costs)

    @cached_method
    def ζ(self, d):
        from sage.all import exp, pi, I
        return self.complex_field()(exp(2*pi*I / d))

    # TODO: Move to HarmonicDifferentials
    def _gen(self, kind, center, n):
        return self.symbolic_ring(self.real_field()).gen((f"{kind}(a{self._differentials._centers.index(center)},?)", n))

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
            sage: Ω = HarmonicDifferentials(S)
            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints
            sage: C = PowerSeriesConstraints(Ω)
            sage: V = Ω._voronoi_diagram()

            sage: cells, segment = next(iter(V.boundaries().items()))
            sage: cell, opposite_cell = list(cells)

            sage: E = C._L2_consistency_voronoi_boundary(cell, segment, opposite_cell)
            sage: F = C._L2_consistency_voronoi_boundary(opposite_cell, -segment, cell)
            sage: (E - F).map_coefficients(lambda c: c if abs(c) > 1e-15 else 0)
            0

        """
        return (self._L2_consistency_voronoi_boundary_f_overline_f(cell, boundary_segment)
                - 2 * self._L2_consistency_voronoi_boundary_Re_f_overline_g(cell, boundary_segment, opposite_cell)
                + self._L2_consistency_voronoi_boundary_g_overline_g(opposite_cell, boundary_segment))

    def _L2_consistency_voronoi_boundary_f_overline_f(self, cell, boundary_segment):
        r"""
        Return the value of `\int_γ f\overline{f}.
        """
        center = self._differentials.surface()(cell.label(), cell.center())

        int = self.CellIntegrator(self, cell)

        κ = cell.root_branch(boundary_segment)
        α = cell.center()
        d = center.angle() - 1

        def f_(part, n, m):
            return int.mixed_integral(part, α, κ, d, n, α, κ, d, m, boundary_segment)

        return self.symbolic_ring().sum([
            (int.Re_a(n) * int.Re_a(m) + int.Im_a(n) * int.Im_a(m)) * f_("Re", n, m) + (int.Re_a(n) * int.Im_a(m) - int.Im_a(n) * int.Re_a(m)) * f_("Im", n, m)
            for n in range(self._differentials.ncoefficients(center))
            for m in range(self._differentials.ncoefficients(center))])

    def _L2_consistency_voronoi_boundary_Re_f_overline_g(self, cell, boundary_segment, opposite_cell):
        r"""
        Return the real part of `\int_γ f\overline{g}.
        """
        center = self._differentials.surface()(cell.label(), cell.center())
        opposite_center = self._differentials.surface()(opposite_cell.label(), opposite_cell.center())

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
            for n in range(self._differentials.ncoefficients(center))
            for m in range(self._differentials.ncoefficients(opposite_center))])

    def _L2_consistency_voronoi_boundary_g_overline_g(self, opposite_cell, boundary_segment):
        r"""
        Return the value of `\int_γ g\overline{g}.
        """
        opposite_center = self._differentials.surface()(opposite_cell.label(), opposite_cell.center())

        jnt = self.CellIntegrator(self, opposite_cell)

        λ = opposite_cell.root_branch(boundary_segment)
        β = opposite_cell.center()
        dd = opposite_center.angle() - 1

        def g_(part, n, m):
            return jnt.mixed_integral(part, β, λ, dd, n, β, λ, dd, m, boundary_segment)

        return self.symbolic_ring().sum([
            (jnt.Re_a(n) * jnt.Re_a(m) + jnt.Im_a(n) * jnt.Im_a(m)) * g_("Re", n, m) + (jnt.Re_a(n) * jnt.Im_a(m) - jnt.Im_a(n) * jnt.Re_a(m)) * g_("Im", n, m)
            for n in range(self._differentials.ncoefficients(opposite_center))
            for m in range(self._differentials.ncoefficients(opposite_center))])

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
            sage: C = PowerSeriesConstraints(Ω)

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

        polygon = self._differentials.surface().polygon(label)
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
            sage: C = PowerSeriesConstraints(Ω)
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
            sage: C = PowerSeriesConstraints(Ω)
            sage: R = C.symbolic_ring()

        We optimize a function in two variables. Since there are no
        constraints, we do not get any Lagrange multipliers in this
        optimization but just for roots of the derivative::

            sage: f = 10*R.gen(0)^2 + 16*R.gen(2)^2
            sage: C.optimize(f)
            sage: C._optimize_cost()
            sage: C
            [20.0000000000000*Re(a0,0), 32.0000000000000*Im(a0,0)]

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
            sage: C = PowerSeriesConstraints(Ω)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # tol 1e-9  # not tested  # TODO: Check this value
            [0.828427124746190*Re(a0,0) + 0.171572875253810*Re(a1,0) - 1.38777878078145e-17*Re(a0,1) - 3.03576608295941e-18*Re(a1,1) + 8.54260702578763e-18*Im(a0,1) + 0.0857864376269050*Im(a1,1) + 0.0473785412436502*Re(a0,2) - 0.0424723326565069*Re(a1,2) - 4.71795158413990e-18*Im(a0,2) - 3.03576608295941e-18*Im(a1,2) - 1.73472347597681e-18*Re(a0,3) + 2.16840434497101e-18*Re(a1,3) + 2.19851947436667e-18*Im(a0,3) - 0.0208152801713079*Im(a1,3) + 0.00487732352790257*Re(a0,4) + 0.0100938339276945*Re(a1,4) - 9.71367022318980e-19*Im(a0,4) + 1.51788304147971e-18*Im(a1,4), 0.828427124746190*Im(a0,0) + 0.171572875253810*Im(a1,0) - 1.38777878078145e-17*Re(a0,1) - 3.03576608295941e-18*Re(a1,1) - 8.54260702578763e-18*Im(a0,1) + 0.0857864376269050*Im(a1,1) + 7.70371977754894e-34*Re(a0,2) - 2.60208521396521e-18*Re(a1,2) - 0.0473785412436502*Im(a0,2) + 0.0424723326565069*Im(a1,2) + 1.73472347597681e-18*Re(a0,3) - 2.16840434497101e-18*Re(a1,3) + 2.19851947436667e-18*Im(a0,3) + 0.0208152801713079*Im(a1,3) - 9.62964972193618e-35*Re(a0,4) - 1.08420217248550e-18*Re(a1,4) + 0.00487732352790257*Im(a0,4) + 0.0100938339276945*Im(a1,4) - 1.00000000000000]


        If we increase precision, we see additional higher imaginary parts.
        These depend on the choice of base point of the integration and will be
        found to be zero by other constraints, not true anymore TODO::

            sage: C = PowerSeriesConstraints(Ω)
            sage: C.require_cohomology(H({b: 1}))
            sage: C  # not tested  # TODO: check this value
            [0.828427124746190*Re(a0,0) + 0.171572875253810*Re(a1,0) - 1.38777878078145e-17*Re(a0,1) - 3.03576608295941e-18*Re(a1,1) + 8.54260702578763e-18*Im(a0,1) + 0.0857864376269050*Im(a1,1) + 0.0473785412436502*Re(a0,2) - 0.0424723326565069*Re(a1,2) - 4.71795158413990e-18*Im(a0,2) - 3.03576608295941e-18*Im(a1,2) - 1.73472347597681e-18*Re(a0,3) + 2.16840434497101e-18*Re(a1,3) + 2.19851947436667e-18*Im(a0,3) - 0.0208152801713079*Im(a1,3) + 0.00487732352790257*Re(a0,4) + 0.0100938339276945*Re(a1,4) - 9.71367022318980e-19*Im(a0,4) + 1.51788304147971e-18*Im(a1,4), 0.828427124746190*Im(a0,0) + 0.171572875253810*Im(a1,0) - 1.38777878078145e-17*Re(a0,1) - 3.03576608295941e-18*Re(a1,1) - 8.54260702578763e-18*Im(a0,1) + 0.0857864376269050*Im(a1,1) + 7.70371977754894e-34*Re(a0,2) - 2.60208521396521e-18*Re(a1,2) - 0.0473785412436502*Im(a0,2) + 0.0424723326565069*Im(a1,2) + 1.73472347597681e-18*Re(a0,3) - 2.16840434497101e-18*Re(a1,3) + 2.19851947436667e-18*Im(a0,3) + 0.0208152801713079*Im(a1,3) - 9.62964972193618e-35*Re(a0,4) - 1.08420217248550e-18*Re(a1,4) + 0.00487732352790257*Im(a0,4) + 0.0100938339276945*Im(a1,4) - 1.00000000000000]

        """
        for cycle in cocycle.parent().homology().gens():
            self.add_constraint(self.integrate(cycle).real() - self.real_field()(cocycle(cycle).real()), rank_check=False)

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
            sage: C = PowerSeriesConstraints(Ω)
            sage: C.require_cohomology(H({a: 1}))
            sage: C.optimize(C._L2_consistency())
            sage: C._optimize_cost()
            sage: C.matrix()  # not tested # TODO: Check this matrix.
            (22 x 22 dense matrix over Real Field with 54 bits of precision,
             (1.00000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000),
             {Im(a0,4): 0,
              Re(a0,2): 1,
              Im(a1,1): 2,
              Re(a1,3): 3,
              Re(a1,1): 4,
              Im(a0,1): 5,
              Im(a0,3): 6,
              Im(a1,3): 7,
              Re(a0,1): 8,
              Im(a1,2): 9,
              Im(a0,0): 10,
              Re(a0,3): 11,
              Im(a1,0): 12,
              Im(a0,2): 13,
              Im(a1,4): 14,
              Re(a1,2): 15,
              Re(a1,4): 16,
              Re(a0,0): 17,
              Re(a0,4): 18,
              Re(a1,0): 19},
            set())

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
            sage: Ω = PowerSeriesConstraints(Ω)
            sage: Ω.power_series_ring(T(0, (1/2, 1/2)))
            Power Series Ring in z0 over Complex Field with 54 bits of precision
            sage: Ω.power_series_ring(T(0, 0))
            Power Series Ring in z1 over Complex Field with 54 bits of precision

        """
        from sage.all import PowerSeriesRing

        if len(args) != 1:
            raise NotImplementedError

        point = args[0]
        return PowerSeriesRing(self.complex_field(), f"z{self._differentials._centers.index(point)}")

    def solve(self, algorithm="eigen+mpfr"):
        r"""
        Return a solution for the system of constraints with minimal error.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: T = translation_surfaces.torus((1, 0), (0, 1))
            sage: T.set_immutable()

            sage: from flatsurf.geometry.harmonic_differentials import PowerSeriesConstraints, HarmonicDifferentials
            sage: Ω = HarmonicDifferentials(T)
            sage: C = PowerSeriesConstraints(Ω)
            sage: R = C.symbolic_ring()
            sage: C.add_constraint(R.gen(0) - R.gen(5))
            sage: C.add_constraint(R.gen(0) - 1)
            sage: C.solve()
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

        series = {point: {} for point in self._differentials._centers}

        for variable, column in decode.items():
            value = solution[column]

            gen, degree = variable.describe()

            # TODO: Use generic functionality to parse variable name.
            if gen.startswith("Re(a"):
                part = 0
            elif gen.startswith("Im(a"):
                part = 1
            else:
                assert False

            center = self._differentials._centers[int(gen.split(",")[0][4:])]

            if degree not in series[center]:
                series[center][degree] = [None, None]

            series[center][degree][part] = value

        series = {point:
                  sum((self.complex_field()(*entry) * self.power_series_ring(point).gen()**k for (k, entry) in series[point].items()), start=self.power_series_ring(point).zero()).add_bigoh(max(series[point]) + 1) for point in series if series[point]}

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
        # TODO: edge() should be immutable.
        holonomy.set_immutable()

        return GeodesicPath(surface, (label, polygon.centroid()), (opposite_label, opposite_polygon.centroid()), holonomy)

    @staticmethod
    def along_edge(surface, label, edge):
        r"""
        Return the :class:`GeodesicPath` along the ``edge`` of the polygon with
        ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()

            sage: from flatsurf.geometry.harmonic_differentials import GeodesicPath
            sage: GeodesicPath.along_edge(S, 0, 0)
            Path (1, 0) from (0, 0) in polygon 0 to (1, 0) in polygon 0

        """
        polygon = surface.polygon(label)
        holonomy = polygon.edge(edge)
        # TODO: edge() should be immutable.
        holonomy.set_immutable()

        return GeodesicPath(surface, (label, polygon.vertex(edge)), (label, polygon.vertex(edge + 1)), holonomy)

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
        # TODO: This only works for convex polygons.
        if self._end[0] == self._start[0] and polygon.get_point_position(self._start[1] + self._holonomy).is_inside():
            return [(self._start[0], OrientedSegment(self._start[1], self._start[1] + self._holonomy))]

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
