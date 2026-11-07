r"""
Reference values for hyperbolic cone surfaces and flat orbifolds.

A signature has three groups:

- the genus `g` and the cone points, given by their order `k_i`: the cone
  angle is `2\pi (1 + k_i)` with `k_i > -1` and `k_i` rational (in the flat
  case, `k_i` is the degree of a zero of an Abelian differential);
- the number `c` of cusps, with no parameter (a cusp is both of angle 0 and of
  length 0);
- the geodesic boundaries; their parametrisation (trace or `\cosh(\ell/2)`
  rather than the length) is not chosen yet.

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

from sage.all import QQ

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
