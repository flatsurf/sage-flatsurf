r"""
Reference values for hyperbolic cone surfaces and flat orbifolds.

Signatures of a hyperbolic orbifold are written ``(g; m_1, ..., m_r; s)``:
genus ``g``, cone points of angle `2\pi/m_i`, and ``s`` cusps.

The Euler characteristic with a divisor is
`\chi(S, \beta) = \chi(S) + \sum (\theta_i/2\pi - 1)` where the `\theta_i` are
the cone angles and a cusp counts as `\theta = 0`. The area is
`-2\pi \chi(S, \beta)`. Without boundary, `\chi(S) = 2 - 2g` and the area is
`2\pi (2g - 2 + s + \sum (1 - 1/m_i))`.

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

from fractions import Fraction

F = Fraction

# Compact hyperbolic orbifolds, rows (g, periods, area / pi).
# Source: S. Katok, Fuchsian groups (1992): section 4.3 Example C, p. 95
# (Gamma_8); Exercise 4.9, p. 110 and its hint, p. 162.
COMPACT_SIGNATURES = [
    (2, (), F(4)),
    (0, (2, 3, 7), F(1, 21)),
    (1, (2,), F(1)),
    (0, (2, 2, 2, 2, 2), F(1)),
    (0, (2, 2, 2, 3), F(1, 3)),
    (0, (3, 3, 4), F(1, 6)),
    (0, (2, 4, 5), F(1, 10)),
]

# Signatures of area 0, rows (g, periods).
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

# Cusped orbifolds, rows (g, periods, cusps, area / pi). Data for work on
# polygons with ideal vertices; the constructions of compact surfaces do not
# use it.
# Source: Katok section 3.5 Example A and Exercise 4.10, p. 110 (modular
# group); section 5.5 Example F, pp. 141-142 (Gamma(2)).
CUSPED_SIGNATURES = [
    (0, (2, 3), 1, F(1, 3)),
    (0, (), 3, F(2)),
]

# Surface of a cube: 8 cone points of angle 3/4 (multiples of 2 pi),
# chi(S) = 2, chi(S, beta) = 0.
# Source: Troyanov 1991, example p. 801.
CUBE_DIVISOR = {
    "cone_points": 8,
    "angle": F(3, 4),
    "euler_characteristic": 2,
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
    "sqrt2_in_radius": [F(4, 3), 0, F(1, 12), 0],
    "cosh_circumradius": (3, 2),  # 3 + 2 sqrt(2)
    "area_over_2pi": F(2),
    "turning_angle": F(3, 8),
    "interior_angle": F(1, 8),
}

# Order 3 point of the modular group, fixed by z -> -1/(z + 1), i.e. the
# matrix [[0, -1], [1, 1]].
# Source: Katok section 3.5 Example A, p. 71.
RHO = {
    "klein": (F(-1, 2), F(0)),
    "matrix": [[0, -1], [1, 1]],
}

# Compact flat orbifolds built by the surface generators, rows
# (name, signature, angles, euler_characteristic); angles are the sorted cone
# angles in multiples of 2 pi. Their groups are the wallpaper groups p1, p2,
# p3, p4, p6. Source of the signatures: as FLAT_SIGNATURES.
FLAT_ORBIFOLDS = [
    ("square_torus", (1, ()), [F(1)], 0),
    ("pillowcase", (0, (2, 2, 2, 2)), [F(1, 2)] * 4, 2),
    ("equilateral_double", (0, (3, 3, 3)), [F(1, 3)] * 3, 2),
    ("isosceles_right_double", (0, (2, 4, 4)), [F(1, 4), F(1, 4), F(1, 2)], 2),
    ("triangle_236_double", (0, (2, 3, 6)), [F(1, 6), F(1, 3), F(1, 2)], 2),
]
