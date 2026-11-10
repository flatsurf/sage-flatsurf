r"""
Reference values for hyperbolic cone surfaces and flat orbifolds.

A signature has three groups:

- the genus `g` and the cone points, given by their order `k_i`: the cone
  angle is `2\pi (1 + k_i)` with `k_i > -1` and `k_i` rational (in the flat
  case, `k_i` is the degree of a zero of an Abelian differential);
- the number `c` of cusps, with no parameter (a cusp is both of angle 0 and of
  length 0);
- the geodesic boundaries, each given by the trace `t = 1 + 2 \cosh \ell` of
  its holonomy in `SO(2, 1)`, where `\ell` is its length.

The Euler characteristic with a divisor is
`\chi(S, \beta) = 2 - 2g + \sum k_i - c - b` where `b` is the number of
geodesic boundaries, and the area is `-2\pi \chi(S, \beta)`.

For orbifolds, where every cone angle is `2\pi/m_i`, Katok's signature
``(g; m_1, ..., m_r; s)`` is a second, more readable form: `k_i = 1/m_i - 1`
and `s = c`.

Each table below keeps the notation of its source: Katok's periods `m_i` for
the orbifolds, orders `k_i` for the general cone points.

Units: angles of flat surfaces (``p.angle()``) are multiples of `2\pi`;
``HyperbolicConvexSet.area()`` is the hyperbolic area divided by `2\pi`;
the angle between consecutive edge geodesics of a polygon is the turning
(exterior) angle in multiples of `2\pi`, so that the interior angle is
`1/2` minus it.

This module is not collected by pytest (its name does not start with
``test_``).
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

from sage.all import AA, QQ, RR, NumberField, polygen


def octagon_vertices():
    r"""
    Return the field `K = \QQ(r)` and the Klein coordinates over `K` of the
    vertices of the regular octagon of ``OCTAGON``, counterclockwise from
    `(r, 0)`.
    """
    x = polygen(QQ, "x")
    f = sum(c * x**i for i, c in enumerate(OCTAGON["klein_radius_polynomial"]))
    K = NumberField(f, "r", embedding=RR(OCTAGON["klein_radius_embedding"]))
    r = K.gen()
    s2 = sum(c * r**i for i, c in enumerate(OCTAGON["sqrt2_in_radius"]))
    return K, [
        (r, 0),
        (r / s2, r / s2),
        (0, r),
        (-r / s2, r / s2),
        (-r, 0),
        (-r / s2, -r / s2),
        (0, -r),
        (r / s2, -r / s2),
    ]


# Compact hyperbolic orbifolds, in Katok's orbifold form, rows
# (g, periods, area / pi).
# Source: S. Katok, Fuchsian groups (1992): section 4.3 Example C, p. 95
# (Gamma_8); Exercise 4.9, p. 110 and its hint, p. 162.
COMPACT_SIGNATURES = [
    (2, (), QQ(4)),
    (0, (2, 3, 7), QQ(1) / 21),
    (1, (2,), QQ(1)),
    (0, (2, 2, 2, 2, 2), QQ(1)),
    (0, (2, 2, 2, 3), QQ(1) / 3),
    (0, (3, 3, 4), QQ(1) / 6),
    (0, (2, 4, 5), QQ(1) / 10),
]

# Signatures of area 0, in Katok's orbifold form, rows (g, periods).
# Source: Katok p. 99 (no Fuchsian group of signature (1; )); hint to
# Exercise 4.9, p. 162 (area 0 for (0; 2, 2, 2, 2), (0; 3, 3, 3),
# (0; 2, 3, 6)); Troyanov, Prescribing curvature on compact surfaces with
# conical singularities (1991), Prop. 2, p. 803 (a flat metric exists iff
# chi(S, beta) = 0).
FLAT_SIGNATURES = [
    (1, ()),
    (0, (2, 2, 2, 2)),
    (0, (3, 3, 3)),
    (0, (2, 4, 4)),
    (0, (2, 3, 6)),
]

# Cusped orbifolds, in Katok's orbifold form, rows
# (g, periods, cusps, area / pi). Data for work on polygons with ideal
# vertices; the constructions of compact surfaces do not use it.
# Source: Katok section 3.5 Example A and Exercise 4.10, p. 110 (modular
# group); section 5.5 Example F, pp. 141-142 (Gamma(2)).
CUSPED_SIGNATURES = [
    (0, (2, 3), 1, QQ(1) / 3),
    (0, (), 3, QQ(2)),
]

# Surface of a cube, in the general form: 8 cone points of order -1/4 (angle
# 3/4 in multiples of 2 pi), chi(S) = 2, chi(S, beta) = 0.
# Source: Troyanov 1991, example p. 801.
CUBE_DIVISOR = {
    "genus": 0,
    "orders": [-QQ(1) / 4] * 8,
    "cusps": 0,
    "boundaries": [],
    "euler_characteristic_with_divisor": 0,
}

# Buser's signature (g, p; q): p geodesic boundaries of free length, q cusps.
# Rows (g, p, q, y_pieces, pairs, dimension).
# Source: P. Buser, Geometry and spectra of compact Riemann surfaces, section
# 4.4, p. 111: 2g - 2 + p + q Y-pieces, 3g - 3 + p + q pairs, dimension
# 6g - 6 + 3p + 2q.
BUSER_DIMENSIONS = [
    (0, 3, 0, 1, 0, 3),
    (2, 0, 0, 2, 3, 6),
    (1, 1, 0, 1, 1, 3),
    (0, 0, 3, 1, 0, 0),
]

# Regular octagon with angles pi/4, fundamental domain of Gamma_8.
# Source: Katok section 4.3 Example C, p. 95 (octagon of area 4 pi made of 8
# triangles with angles pi/4, pi/8, pi/8) and Theorem 1.5.2 (iii), p. 16
# (Cosine Rule II gives cosh R = cot^2(pi/8) = 3 + 2 sqrt(2)).
# The Klein radius r is a root of x^4 + 32 x^2 - 32; polynomials are given
# as coefficient lists, constant term first.
OCTAGON = {
    "klein_radius_polynomial": [-32, 0, 32, 0, 1],
    "klein_radius_embedding": 0.98517143100941,
    "sqrt2_in_radius": [QQ(4) / 3, 0, QQ(1) / 12, 0],
    "cosh_circumradius": (3, 2),  # 3 + 2 sqrt(2)
    "area_over_2pi": QQ(2),
    "turning_angle": QQ(3) / 8,
    "interior_angle": QQ(1) / 8,
}

# Side pairings of the octagon above (sides numbered from the vertex of Klein
# coordinates (r, 0), counterclockwise), each giving a closed surface of genus
# 2 whose 8 vertices form one regular point. Rows: pairs of glued sides and the
# trace of the matrix in SO(2, 1)+ of each pairing, the same for all of them,
# as a polynomial in r (coefficient list, constant term first). The sides all
# have cosh^2 of their length equal to OCTAGON_SIDE_COSH2.
# Source: Katok section 4.3 Example C, pp. 95-98 (the gluing "katok", whose
# generators satisfy [A_1, B_1][A_2, B_2] = Id); the traces and the side
# length are exact values in Q(r), computed with sage-flatsurf (isometries
# sending a side onto the other one).
OCTAGON_SIDE_COSH2 = [QQ(331) / 3, 0, QQ(10) / 3, 0]
OCTAGON_GLUINGS = {
    "opposite": ([(0, 4), (1, 5), (2, 6), (3, 7)], [QQ(65) / 3, 0, QQ(2) / 3, 0]),
    "katok": ([(0, 2), (1, 3), (4, 6), (5, 7)], [QQ(31) / 3, 0, QQ(1) / 3, 0]),
}

# Regular polygons about the origin, moved by the boost
# [[5/4, 0, 3/4], [0, 1, 0], [3/4, 0, 5/4]] (so that every centre has Klein
# coordinates (3/5, 0)). Rows: name, order n of the rotation of the polygon,
# trace of that rotation as a polynomial in the generator of the field
# (coefficient list, constant term first). The traces are 1 + 2 cos(2 pi / n)
# with n = 4, 8, 3, written in the field of the vertices (Q, Q(r) of OCTAGON,
# Q(sqrt(3))).
# Source: exact values computed with sage-flatsurf (rotation sending a side
# onto the next one).
CELL_ROTATIONS = [
    ("square", 4, [QQ(1)]),
    ("octagon", 8, [QQ(7) / 3, 0, QQ(1) / 12, 0]),
    ("triangle", 3, [QQ(0)]),
]
CELL_CENTRE = (QQ(3) / 5, QQ(0))
BOOST = [[QQ(5) / 4, 0, QQ(3) / 4], [0, 1, 0], [QQ(3) / 4, 0, QQ(5) / 4]]


def cell_vertices(name):
    r"""
    Return the field and the Klein coordinates of the vertices of the regular
    polygon ``name`` of ``CELL_ROTATIONS`` about the origin, before the boost.
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


def hexagon_vertices():
    r"""
    Return the Klein coordinates over ``AA`` of the vertices of the regular
    right-angled hexagon of ``HEXAGON_SURFACES``, of Klein radius
    `\rho = \sqrt{6}/3`, counterclockwise from `(\rho, 0)`.
    """
    rho = AA(6).sqrt() / 3
    h = AA(2).sqrt() / 2
    return [
        (rho, 0),
        (rho / 2, h),
        (-rho / 2, h),
        (-rho, 0),
        (-rho / 2, -h),
        (rho / 2, -h),
    ]


# Regular right-angled hexagon about the origin, and the surfaces glued from
# two copies of it (polygons 0 and 1). The hexagon is made of six triangles
# with angles pi/3 (at the centre), pi/4, pi/4. Source: Katok, Fuchsian groups
# (1992), Theorem 1.5.2 (iii), p. 16 (Cosine Rule II):
# - the circumradius R has cosh R = (cos(pi/4) + cos(pi/3) cos(pi/4)) /
#   (sin(pi/3) sin(pi/4)) = sqrt(3), so the Klein radius is
#   tanh R = sqrt(6)/3;
# - the side s has cosh s = (cos(pi/3) + cos^2(pi/4)) / sin^2(pi/4) = 2.
# A boundary of length l has the trace t = 1 + 2 cosh l in SO(2, 1).
# - "pants": each of the three boundaries is two sides, cosh 2s = 2 * 2^2 - 1
#   = 7, t = 1 + 2 * 7 = 15;
# - "torus": the one boundary is six sides, cosh 6s = 4 * 7^3 - 3 * 7 = 1351,
#   t = 1 + 2 * 1351 = 2703.
# Rows: the gluings, the genus and the boundary traces.
HEXAGON_SURFACES = {
    "pants": ([((0, 0), (1, 0)), ((0, 2), (1, 4)), ((0, 4), (1, 2))], 0, (15,) * 3),
    "torus": ([((0, 0), (1, 0)), ((0, 2), (1, 2)), ((0, 4), (1, 4))], 1, (2703,)),
}


# Doubles of hyperbolic triangles: a triangle of angles pi/m_i and its mirror
# image, glued along their three edges, a sphere with three cone points whose
# group is the triangle group of periods (m_1, m_2, m_3). Rows: the periods at
# vertices 0, 1, 2 of the triangle; the field of the vertices (defining
# polynomial as a coefficient list, constant term first, and an approximation
# of the embedding of its generator); the Klein coordinates of the vertices
# 1 = (x_1, 0) and 2 = (x_2, y_2), vertex 0 at the origin (coefficient lists in
# the generator); a defining polynomial of the invariant trace field.
# Sources: the coordinates are exact values computed with sage-flatsurf in AA
# from the angles (the sides by Katok, Fuchsian groups (1992), Theorem 1.5.2
# (iii), p. 16), written over their minimal field; the invariant trace field
# k_0 = Q(cos^2(pi/m_i), cos(pi/m_1) cos(pi/m_2) cos(pi/m_3)) is Takeuchi,
# Arithmetic triangle groups (1977), Propositions 4 and 5, pp. 96-99: Q(sqrt 5)
# for (2, 4, 5), Q(cos(2 pi/7)) for (2, 3, 7), Q(sqrt 2) for (3, 3, 4) and
# (4, 4, 4).
TRIANGLE_DOUBLES = [
    (
        (4, 5, 2),
        ([-4, 0, -2, 0, 1], -1.79890743994787),
        ([0, -2, 0, QQ(1) / 2], [0, -1, 0, QQ(1) / 4], [0, -1, 0, QQ(1) / 4]),
        [-1, -1, 1],
    ),
    (
        (4, 2, 5),
        ([-1, 0, -1, 0, 1], 1.27201964951407),
        ([0, 2, 0, -1], [0, 2, 0, -1], [0, 2, 0, -1]),
        [-1, -1, 1],
    ),
    (
        (2, 7, 3),
        ([-1, 0, -2, 0, 1, 0, 1], 1.11668240951376),
        ([0, 2, 0, 0, 0, -1], [0], [0, -1, 0, 1]),
        [-1, -2, 1, 1],
    ),
    (
        (2, 3, 7),
        ([-1, 0, -2, 0, 1, 0, 1], -1.11668240951376),
        ([0, 1, 0, -1], [0], [0, -2, 0, 0, 0, 1]),
        [-1, -2, 1, 1],
    ),
    (
        (4, 3, 3),
        ([-2, 0, 0, 0, 1], 1.18920711500272),
        ([0, 2, 0, -1], [0, -1, 0, 1], [0, -1, 0, 1]),
        [-2, 0, 1],
    ),
    (
        (4, 4, 4),
        ([-1, 0, -2, 0, 1], -1.55377397403004),
        ([0, -3, 0, 1], [0, 2, 0, -1], [0, 2, 0, -1]),
        [-2, 0, 1],
    ),
]


def triangle_double_vertices(row):
    r"""
    Return the field and the Klein coordinates of the vertices of the triangle
    of the row ``row`` of ``TRIANGLE_DOUBLES``, and of its mirror image, whose
    vertex `j` is the image of the vertex `1 - j` under `(x, y) \mapsto (x,
    -y)`, so that edge `i` of the triangle is the image of edge `-i` of the
    mirror.
    """
    _, (polynomial, embedding), (x1, x2, y2), _ = row
    x = polygen(QQ, "x")
    K = NumberField(
        sum(c * x**i for i, c in enumerate(polynomial)), "a", embedding=RR(embedding)
    )
    a = K.gen()
    x1, x2, y2 = (sum(c * a**i for i, c in enumerate(p)) for p in (x1, x2, y2))
    front = [(0, 0), (x1, 0), (x2, y2)]
    return K, front, [(x1, 0), (0, 0), (x2, -y2)]


# Order 3 point of the modular group, fixed by z -> -1/(z + 1), i.e. the
# matrix [[0, -1], [1, 1]].
# Source: Katok section 3.5 Example A, p. 71.
RHO = {
    "klein": (-QQ(1) / 2, QQ(0)),
    "matrix": [[0, -1], [1, 1]],
}

# Compact flat orbifolds built by the surface generators, rows
# (name, signature, angles, euler_characteristic); the signature is in
# Katok's orbifold form and the angles are the sorted cone angles in
# multiples of 2 pi. Their groups are the wallpaper groups p1, p2, p3, p4,
# p6. Source of the signatures: as FLAT_SIGNATURES.
FLAT_ORBIFOLDS = [
    ("square_torus", (1, ()), [QQ(1)], 0),
    ("pillowcase", (0, (2, 2, 2, 2)), [QQ(1) / 2] * 4, 2),
    ("equilateral_double", (0, (3, 3, 3)), [QQ(1) / 3] * 3, 2),
    (
        "isosceles_right_double",
        (0, (2, 4, 4)),
        [QQ(1) / 4, QQ(1) / 4, QQ(1) / 2],
        2,
    ),
    (
        "triangle_236_double",
        (0, (2, 3, 6)),
        [QQ(1) / 6, QQ(1) / 3, QQ(1) / 2],
        2,
    ),
]
