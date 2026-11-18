r"""
Checks of the isomorphisms from hyperbolic planes of the quadratic forms
`x^2 + b y^2 - c z^2` to the plane of `x^2 + y^2 - z^2` against hard-coded
matrices, images and refusals.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2026 Vincent Delecroix
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
# ****************************************************************************

import pytest

from sage.all import (
    AA,
    QQ,
    RDF,
    RR,
    NumberField,
    diagonal_matrix,
    identity_matrix,
    matrix,
    polygen,
)

from flatsurf import HyperbolicPlane
from flatsurf.geometry.hyperbolic_isometry import HyperbolicIsometryGroup


TRIANGLE = [(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)]

QFSOLVE_3 = matrix(
    QQ,
    [
        [QQ(7) / 4, QQ(5) / 2, -QQ(15) / 4],
        [0, 3, -3],
        [-QQ(5) / 4, -QQ(7) / 2, QQ(21) / 4],
    ],
)


@pytest.fixture(scope="module")
def K():
    x = polygen(QQ)
    return NumberField(x**2 - 3, "a", embedding=AA(3).sqrt())


@pytest.fixture(scope="module")
def L():
    x = polygen(QQ)
    return NumberField(x**2 - 2, "r", embedding=AA(2).sqrt())


@pytest.fixture(scope="module")
def H3():
    return HyperbolicPlane(QQ, form=[1, 2, -3])


def test_diagonal():
    f = HyperbolicPlane(QQ, form=[1, 4, -9]).standard_form_isometry()
    assert f.matrix() == diagonal_matrix(QQ, [1, 2, 3])
    assert f.codomain() is HyperbolicPlane(QQ)


def test_square_root_of_c():
    f = HyperbolicPlane(QQ, form=[1, 2, -4]).standard_form_isometry()
    assert f.matrix() == matrix(
        QQ,
        [
            [QQ(9) / 8, 0, -QQ(7) / 4],
            [0, 1, 0],
            [-QQ(7) / 8, 0, QQ(9) / 4],
        ],
    )


def test_square_root_of_c_over_b():
    f = HyperbolicPlane(QQ, form=[1, 2, -8]).standard_form_isometry()
    assert f.matrix() == matrix(
        QQ,
        [
            [0, QQ(33) / 8, -QQ(31) / 4],
            [-1, 0, 0],
            [0, -QQ(31) / 8, QQ(33) / 4],
        ],
    )


def test_square_root_of_c_over_number_field(K, H3):
    a = K.gen()
    f = H3.standard_form_isometry(K)
    assert f.domain() is H3.change_ring(K)
    assert f.codomain() is HyperbolicPlane(K)
    assert f.matrix() == matrix(
        K,
        [
            [QQ(7) / 12 * a, 0, -QQ(5) / 4],
            [0, 1, 0],
            [-QQ(5) / 12 * a, 0, QQ(7) / 4],
        ],
    )

    P = H3.polygon(vertices=TRIANGLE)
    Q = f(P)
    assert [v.coordinates(model="klein") for v in Q.vertices()] == [
        (-QQ(5) / 7, 0),
        (QQ(48) / 563 * a - QQ(385) / 563, 0),
        (-QQ(5) / 7, QQ(2) / 7),
    ]
    assert (~f)(Q) == P.change_ring(K)
    assert (~f).matrix().base_ring() is K


def test_qfsolve(H3):
    f = H3.standard_form_isometry()
    assert f.matrix() == QFSOLVE_3
    assert H3.standard_form_isometry(algorithm="qfsolve") == f

    g = HyperbolicPlane(QQ, form=[1, 7, -11]).standard_form_isometry()
    assert g.matrix() == matrix(
        QQ,
        [
            [QQ(51) / 56, QQ(37) / 16, -QQ(407) / 112],
            [0, QQ(11) / 4, -QQ(11) / 4],
            [-QQ(37) / 56, -QQ(51) / 16, QQ(561) / 112],
        ],
    )


def test_anisotropic_over_rationals():
    H = HyperbolicPlane(QQ, form=[1, 1, -3])
    message = "isometry is not defined over Rational Field; adjoin a square root of 3"
    for algorithm in (None, "qfsolve", "conic", "norm_equation"):
        with pytest.raises(ValueError, match=message):
            H.standard_form_isometry(algorithm=algorithm)


def test_isotropic_vector(H3):
    f = H3.standard_form_isometry(isotropic_vector=(1, 1, 1))
    assert f.matrix() == QFSOLVE_3
    with pytest.raises(ValueError, match="is not isotropic for x\\^2 \\+ 2\\*y\\^2"):
        H3.standard_form_isometry(isotropic_vector=(1, 0, 1))
    with pytest.raises(ValueError, match="the vector \\(0, 0, 0\\) is not isotropic"):
        H3.standard_form_isometry(isotropic_vector=(0, 0, 0))
    with pytest.raises(ValueError, match="isotropic_vector and algorithm are exclusive"):
        H3.standard_form_isometry(isotropic_vector=(1, 1, 1), algorithm="conic")


def test_norm_equation_and_conic(L):
    # The isotropic vector found depends on the version of PARI: the matrices
    # are checked by their properties.
    H = HyperbolicPlane(L, form=[1, 1, -3])
    Q = H.quadratic_form()
    J = diagonal_matrix(L, [1, 1, -1])
    f = H.standard_form_isometry()
    assert H.standard_form_isometry(algorithm="norm_equation") == f
    for g in (f, H.standard_form_isometry(algorithm="conic")):
        assert g.domain() is H
        assert g.codomain() is HyperbolicPlane(L)
        A = g.matrix()
        assert A.base_ring() is L
        scale = (A.transpose() * J * A)[0, 0]
        assert scale > 0
        assert A.transpose() * J * A == scale * Q
        assert A.det() > 0
        assert A[2, 2] > 0


def test_definite_at_a_real_place(L):
    r = L.gen()
    H = HyperbolicPlane(L, form=[1, 1, -(r - 1)])
    message = "adjoin a square root of r - 1"
    for algorithm in (None, "conic", "norm_equation"):
        with pytest.raises(ValueError, match=message):
            H.standard_form_isometry(algorithm=algorithm)


def test_large_degree():
    x = polygen(QQ)
    N = NumberField(x**9 - 2, "a", embedding=AA(2).nth_root(9))
    a = N.gen()
    H = HyperbolicPlane(N, form=[1, 1, -(1 + a**2)])
    with pytest.raises(ValueError) as error:
        H.standard_form_isometry()
    assert str(error.value) == (
        "the search for an isotropic vector of x^2 + y^2 - (a^2 + 1)*z^2 has no "
        "time bound over a number field of degree 9; pass "
        'algorithm="norm_equation", "conic" or "extension", an isotropic_vector, '
        "or a ring with a square root of a^2 + 1 or of a^2 + 1"
    )


def test_extension_over_rationals():
    H = HyperbolicPlane(QQ, form=[1, 1, -3])
    f = H.standard_form_isometry(algorithm="extension")
    F = f.domain().base_ring()
    assert F.degree() == 2
    assert F.polynomial() == polygen(QQ) ** 2 - 3
    s = F.gen()
    assert AA(s) == AA(3).sqrt()
    assert f.matrix() == diagonal_matrix(F, [1, 1, s])
    assert f.codomain() is HyperbolicPlane(F)
    p = f(H.point(QQ(1) / 2, QQ(1) / 2, model="klein"))
    assert p.coordinates(model="klein") == (s / 6, s / 6)


def test_extension_over_number_field(L):
    r = L.gen()
    H = HyperbolicPlane(L, form=[1, 1, -(r + 1)])
    f = H.standard_form_isometry(algorithm="extension")
    F = f.domain().base_ring()
    s = F.gen()
    assert F.polynomial() == polygen(QQ) ** 4 - 2 * polygen(QQ) ** 2 - 1
    assert F.has_coerce_map_from(L)
    assert F(r) == s**2 - 1
    assert f.matrix() == diagonal_matrix(F, [1, 1, s])
    Q = f(H.polygon(vertices=TRIANGLE))
    assert Q.parent() is HyperbolicPlane(F)


def test_extension_keeps_the_ring_with_the_roots():
    f = HyperbolicPlane(QQ, form=[1, 4, -9]).standard_form_isometry(
        algorithm="extension"
    )
    assert f.domain().base_ring() is QQ
    assert f.matrix() == diagonal_matrix(QQ, [1, 2, 3])


def test_wrong_ring_and_name(H3):
    with pytest.raises(
        ValueError, match='algorithm "qfsolve" needs the rational field, not'
    ):
        H3.standard_form_isometry(AA, algorithm="qfsolve")
    for algorithm in ("conic", "norm_equation"):
        with pytest.raises(
            ValueError, match="needs the rational field or a number field, not"
        ):
            H3.standard_form_isometry(AA, algorithm=algorithm)
    with pytest.raises(ValueError, match="algorithm must be None"):
        H3.standard_form_isometry(algorithm="pari")


def test_standard_plane():
    H = HyperbolicPlane(QQ)
    assert H.standard_form_isometry().is_identity()
    f = H.standard_form_isometry(AA)
    assert f.domain() is HyperbolicPlane(AA)
    assert f.matrix() == identity_matrix(AA, 3)


def test_real_algebraic_numbers(H3):
    f = H3.standard_form_isometry(AA)
    assert f.matrix() == diagonal_matrix(AA, [1, AA(2).sqrt(), AA(3).sqrt()])


def test_inexact_rings(H3):
    for ring in (RR, RDF):
        with pytest.raises(NotImplementedError, match="needs an exact ring"):
            H3.standard_form_isometry(ring)
        with pytest.raises(NotImplementedError, match="needs an exact ring"):
            HyperbolicPlane(QQ).standard_form_isometry(ring)


def test_isomorphism_needs_a_plane(H3):
    with pytest.raises(ValueError, match="other must be a hyperbolic plane"):
        H3.isomorphism(QQ)


def test_polygon_over_rationals_through_number_field(K, H3):
    f = H3.standard_form_isometry(K)
    P = H3.polygon(vertices=[(0, 0), (QQ(1) / 2, -QQ(1) / 2), (QQ(3) / 4, 0)])
    Q = f(P)
    assert Q.parent() is HyperbolicPlane(K)
    assert Q._cosh2_sides() == P._cosh2_sides()
    assert Q.vertices()[0].coordinates(model="klein") == (-QQ(5) / 7, 0)


def test_psl2_matrix(K):
    a = K.gen()
    G3 = HyperbolicIsometryGroup(QQ, form=[1, 2, -3])
    assert G3.one().psl2_matrix() == identity_matrix(QQ, 2)

    g = G3([[-QQ(1) / 3, -QQ(4) / 3, 0], [QQ(2) / 3, -QQ(1) / 3, 0], [0, 0, 1]])
    assert g.trace() == QQ(1) / 3

    H3 = G3.hyperbolic_plane()
    assert H3.standard_form_isometry()(g).matrix() == matrix(
        QQ,
        [
            [-QQ(9) / 2, -QQ(17) / 3, -QQ(43) / 6],
            [-1, -3, -3],
            [QQ(9) / 2, QQ(19) / 3, QQ(47) / 6],
        ],
    )
    with pytest.raises(
        ValueError,
        match="isometry is not defined over Rational Field; adjoin a square root of 49/12",
    ):
        g.psl2_matrix()

    assert H3.standard_form_isometry(K)(g).matrix() == matrix(
        K,
        [
            [-QQ(31) / 18, -QQ(7) / 9 * a, -QQ(35) / 18],
            [QQ(7) / 9 * a, -QQ(1) / 3, QQ(5) / 9 * a],
            [QQ(35) / 18, QQ(5) / 9 * a, QQ(43) / 18],
        ],
    )
    A = g.psl2_matrix(K)
    assert A == matrix(
        K,
        [
            [QQ(1) / 3 * a + QQ(5) / 6, QQ(7) / 6],
            [-QQ(7) / 6, QQ(1) / 3 * a - QQ(5) / 6],
        ],
    )
    assert A.trace() ** 2 - 1 == g.trace()
