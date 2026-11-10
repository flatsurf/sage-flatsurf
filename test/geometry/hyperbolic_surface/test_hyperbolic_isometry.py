r"""
Checks of the isometries of the hyperbolic plane against hard-coded matrices
of side pairings and of tangent vector maps.
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
    RR,
    I,
    NumberField,
    diagonal_matrix,
    matrix,
    polygen,
    vector,
)

from flatsurf import HyperbolicPlane
from flatsurf.geometry.hyperbolic_isometry import (
    HyperbolicIsometryGroup,
    segment_isometry,
    segment_midpoint,
    tangent_vector_isometry,
)

from hyperbolic_surface_oracles import (
    BOOST,
    CELL_CENTRE,
    CELL_ROTATIONS,
    octagon_vertices,
)


# Klein coordinates of a Lambert quadrilateral (three right angles).
QUADRILATERAL = [(0, 0), (QQ(1) / 2, 0), (QQ(1) / 2, QQ(1) / 2), (0, QQ(1) / 2)]

# Klein coordinates of a non-convex quadrilateral, reflex at the origin.
ARROW = [(0, -QQ(1) / 2), (QQ(1) / 2, QQ(1) / 2), (0, 0), (-QQ(1) / 2, QQ(1) / 2)]

QUARTER_TURN = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]

J = diagonal_matrix(QQ, [1, 1, -1])


@pytest.fixture(scope="module")
def H():
    return HyperbolicPlane(QQ)


@pytest.fixture(scope="module")
def G():
    return HyperbolicIsometryGroup(QQ)


def pairing(P, i, j):
    r"""
    Return the matrix sending the start of the side ``i`` of ``P`` to the end
    of the side ``j`` and the end of ``i`` to the start of ``j``.
    """
    return segment_isometry(P.edge_segment(i), -P.edge_segment(j))


def assert_images(M, images):
    V = M.base_ring() ** 2
    for source, target in images:
        X, Y, Z = M * vector(M.base_ring(), (source[0], source[1], 1))
        assert V((X / Z, Y / Z)) == V(target)


def test_quadrilateral_pairings(H):
    P = H.polygon(vertices=QUADRILATERAL)
    half = QQ(1) / 2

    M = pairing(P, 0, 3)
    assert M == matrix(QQ, QUARTER_TURN)
    assert M.trace() == 1
    assert M.transpose() * J * M == J
    assert_images(M, [((0, 0), (0, 0)), ((half, 0), (0, half))])

    M = pairing(P, 1, 2)
    assert M == matrix(QQ, [[0, 1, 0], [-QQ(5) / 3, 0, QQ(4) / 3], [-QQ(4) / 3, 0, QQ(5) / 3]])
    assert M.trace() == QQ(5) / 3
    assert M.transpose() * J * M == J
    assert_images(M, [((half, 0), (0, half)), ((half, half), (half, half))])


def test_arrow_pairings(H):
    P = H.polygon(vertices=ARROW)
    half = QQ(1) / 2

    M = pairing(P, 0, 3)
    assert M == matrix(QQ, [[13, -16, -8], [16, 11, -4], [-8, 4, 21]]) / 19
    assert M.trace() == QQ(45) / 19
    assert M.transpose() * J * M == J
    assert_images(M, [((0, -half), (0, -half)), ((half, half), (-half, half))])

    M = pairing(P, 1, 2)
    assert M == matrix(QQ, QUARTER_TURN)
    assert M.trace() == 1
    assert_images(M, [((half, half), (-half, half)), ((0, 0), (0, 0))])


def test_self_glued_edge():
    edge = ((0, 0), (QQ(1) / 2, 0))
    reverse = ((QQ(1) / 2, 0), (0, 0))
    with pytest.raises(ValueError):
        segment_isometry(edge, reverse)

    x = polygen(QQ, "x")
    L = NumberField(x**2 - 3, "r", embedding=AA(3).sqrt())
    r = L.gen()
    M = segment_isometry(edge, reverse, ring=L)
    assert M == matrix(L, [[-2 * r / 3, 0, r / 3], [0, -1, 0], [-r / 3, 0, 2 * r / 3]])
    assert M.trace() == -1
    JL = J.change_ring(L)
    assert M.transpose() * JL * M == JL
    assert_images(M, [((0, 0), (QQ(1) / 2, 0)), ((QQ(1) / 2, 0), (0, 0))])


def test_tangent_vectors(H):
    origin = (0, 0)
    M = tangent_vector_isometry((origin, (1, 0)), (origin, (0, 1)))
    assert M == matrix(QQ, QUARTER_TURN)
    assert tangent_vector_isometry((origin, (2, 0)), (origin, (0, 1))) == M
    assert tangent_vector_isometry((H(I), (2, 0)), (H(I), (0, 1))) == M
    with pytest.raises(ValueError):
        tangent_vector_isometry(((QQ(1) / 2, 0), (1, 0)), (origin, (1, 0)))


def test_segment_lengths():
    with pytest.raises(ValueError):
        segment_isometry(((0, 0), (QQ(1) / 2, 0)), ((0, 0), (0, QQ(1) / 3)))


def test_psl2_matrix(H, G):
    A = matrix(QQ, [[2, 3], [1, 2]])
    g = G(H._isometry_gl2_to_sim12(A))
    assert g.matrix() == matrix(QQ, [[7, -4, 8], [4, -1, 4], [8, -4, 9]])
    assert g.trace() == 15
    assert g.matrix().determinant() == 1
    assert g.psl2_matrix() == A
    assert G(A) == g

    boost = G([[QQ(5) / 4, 0, QQ(3) / 4], [0, 1, 0], [QQ(3) / 4, 0, QQ(5) / 4]])
    with pytest.raises(ValueError):
        boost.psl2_matrix()
    B = boost.psl2_matrix(AA)
    assert H._isometry_gl2_to_sim12(B) == boost.matrix()
    assert B.determinant() == 1


def test_action(H, G):
    g = G(QUARTER_TURN)
    v = g * vector(QQ, (QQ(1) / 2, 0))
    assert v == vector(QQ, (0, QQ(1) / 2))
    assert v.is_immutable()
    assert g((QQ(1) / 2, 0)) == v
    assert g * H(2 * I) == H.point(-QQ(3) / 5, QQ(4) / 5, model="half_plane")
    assert g(H(2 * I)) == g * H(2 * I)
    P = H.polygon(vertices=QUADRILATERAL)
    assert g(P) == P.apply_isometry(g) == P.apply_isometry(matrix(QQ, QUARTER_TURN))
    assert g(P).vertices()[1] == H.point(0, QQ(1) / 2, model="klein")


def test_group_laws(G):
    g = G(QUARTER_TURN)
    h = G([[7, -4, 8], [4, -1, 4], [8, -4, 9]])
    v = vector(QQ, (QQ(1) / 3, -QQ(1) / 4))
    assert (g * h) * v == g * (h * v)
    assert (h * g) * v == h * (g * v)
    assert ~g * g == G.one()
    assert ~h * h == G.one()
    assert h * ~h == G.one()


def test_mixed_rings():
    x = polygen(QQ)
    L = NumberField(x**2 - 3, "s3", embedding=AA(3).sqrt())
    M = NumberField(x**2 - 2, "s2", embedding=AA(2).sqrt())
    s3 = L.gen()
    G_L = HyperbolicIsometryGroup(L)
    G_M = HyperbolicIsometryGroup(M)
    g = G_L(matrix(L, [[s3 / 2, -QQ(1) / 2, 0], [QQ(1) / 2, s3 / 2, 0], [0, 0, 1]]))
    H = HyperbolicPlane(QQ)
    P = H.polygon(vertices=[(0, 0), (QQ(1) / 2, 0), (0, QQ(1) / 2)])
    assert g(P) == g(P.change_ring(L))
    assert g(P).parent().base_ring() is L
    assert G_L.one()(P).parent().base_ring() is L
    assert (g * H(2 * I)).parent().base_ring() is L
    assert (G_L.one() * G_M.one()).parent() is HyperbolicIsometryGroup(AA)

    s2 = M.gen()
    P_M = HyperbolicPlane(M).polygon(vertices=[(0, 0), (s2 / 4, 0), (0, s2 / 4)])
    assert g(P_M) == HyperbolicIsometryGroup(AA)(g)(P_M.change_ring(AA))


def proportional(v, w):
    return matrix(QQ, [v, w]).rank() == 1


def cell_vertices(name):
    r"""
    Return the field and the Klein coordinates of the vertices of a regular
    polygon about the origin.
    """
    if name == "square":
        half = QQ(1) / 2
        return QQ, [(half, 0), (0, half), (-half, 0), (0, -half)]
    if name == "octagon":
        return octagon_vertices()
    x = polygen(QQ)
    K = NumberField(x**2 - 3, "s3", embedding=RR(1.7))
    s3 = K.gen()
    return K, [(QQ(1) / 2, 0), (-QQ(1) / 4, s3 / 4), (-QQ(1) / 4, -s3 / 4)]


@pytest.mark.parametrize("name,n,trace", CELL_ROTATIONS)
def test_cell_rotation(name, n, trace):
    K, vertices = cell_vertices(name)
    G = HyperbolicIsometryGroup(K)
    boost = G(matrix(K, BOOST))
    vertices = [boost(v) for v in vertices]

    g = G(segment_isometry((vertices[0], vertices[1]), (vertices[1], vertices[2])))
    assert g.trace() == sum(c * K.gen() ** i for i, c in enumerate(trace))
    assert (g**n).is_one()
    assert not any((g**k).is_one() for k in range(1, n))
    assert g.fixed_point() == HyperbolicPlane(K).point(*CELL_CENTRE, model="klein")


def test_cell_boost():
    K, vertices = cell_vertices("square")
    boost = HyperbolicIsometryGroup(QQ)(matrix(QQ, BOOST))
    assert [tuple(boost(v)) for v in vertices] == [
        (QQ(11) / 13, 0),
        (QQ(3) / 5, QQ(2) / 5),
        (QQ(1) / 7, 0),
        (QQ(3) / 5, -QQ(2) / 5),
    ]


def test_fixed_points(G):
    rho = G(matrix(QQ, [[0, -1], [1, 1]]))
    assert rho.trace() == 0
    assert proportional(rho._fixed_vector(), (1, 0, -2))
    assert rho.fixed_point() == HyperbolicPlane(QQ).point(
        -QQ(1) / 2, 0, model="klein"
    )

    parabolic = G(matrix(QQ, [[1, 1], [0, 1]]))
    assert parabolic.trace() == 3
    assert parabolic.fixed_point() == HyperbolicPlane(QQ).point(0, 1, model="klein")

    boost = G(matrix(QQ, BOOST))
    assert boost.trace() == QQ(7) / 2
    assert proportional(boost._fixed_vector(), (0, 1, 0))
    with pytest.raises(NotImplementedError):
        boost.fixed_point()

    hyperbolic = G(matrix(QQ, [[2, 1], [1, 1]]))
    assert hyperbolic.trace() == 8
    assert proportional(hyperbolic._fixed_vector(), (1, -2, 0))
    with pytest.raises(NotImplementedError):
        hyperbolic.fixed_point()

    with pytest.raises(ValueError):
        G.one().fixed_point()


def test_half_turn_fixed_point(G):
    s = ((QQ(1) / 2, 0), (0, QQ(1) / 2))
    M = segment_isometry(s, (s[1], s[0]))
    assert M == matrix(QQ, [[-8, -1, 4], [-1, -8, 4], [-4, -4, 9]]) / 7
    g = G(M)
    assert proportional(g._fixed_vector(), (1, 1, 4))
    assert g.fixed_point() == HyperbolicPlane(QQ).point(
        QQ(1) / 4, QQ(1) / 4, model="klein"
    )


def test_segment_midpoint():
    half = QQ(1) / 2
    assert segment_midpoint((half, 0), (0, half)) == vector(QQ, (QQ(1) / 4, QQ(1) / 4))
    with pytest.raises(ValueError):
        segment_midpoint((0, 0), (half, 0))

    x = polygen(QQ)
    L = NumberField(x**2 - 3, "s3", embedding=AA(3).sqrt())
    assert segment_midpoint((0, 0), (half, 0), ring=L) == vector(L, (2 - L.gen(), 0))
    m = segment_midpoint((0, 0), (half, 0), ring=AA)
    assert m == vector(AA, (2 - AA(3).sqrt(), 0))


def test_point_on_isometry_right(H, G):
    with pytest.raises(TypeError):
        H(I) * G.one()
    m = matrix(QQ, [[1, 2], [3, 4]])
    assert H(I) * m == H.point(-QQ(7) / 5, QQ(1) / 5, model="half_plane")


def test_constructor_errors(G):
    with pytest.raises(ValueError):
        G([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    with pytest.raises(ValueError):
        G(matrix(QQ, [[1, 0], [0, -1]]))
