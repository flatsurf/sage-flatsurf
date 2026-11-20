r"""
The signature of a cone surface.

The signature of a compact cone surface, of any curvature (flat, hyperbolic
or spherical), is topological. It has three groups:

- the genus `g` and the cone points, given by their order `k_i`: the cone
  angle is `2\pi (1 + k_i)` with `k_i > -1`, `k_i \neq 0`, a rational number;
- the number `c` of cusps;
- the boundary components, `b` of them. On a hyperbolic surface with
  geodesic boundary, each is given by the trace of its holonomy in
  `\mathrm{SO}(2, 1)`; on a surface of another geometry the trace is
  ``None``.

The topological Euler characteristic is `\chi(S) = 2 - 2g - c - b`. The
divisor Euler characteristic of Troyanov (*Prescribing curvature on compact
surfaces with conical singularities*, 1991), for the divisor `\beta = \sum k_i
p_i`, is `\chi(S, \beta) = 2 - 2g - c - b + \sum k_i`. By the Gauss-Bonnet
theorem, for any cone angles and with geodesic boundary, a flat surface has
`\chi(S, \beta) = 0`, a hyperbolic surface has `\chi(S, \beta) < 0` and area
`-2 \pi \chi(S, \beta)`, and a spherical surface has `\chi(S, \beta) > 0`.

The cone orders of a flat surface are the data of its stratum: `k_i = m_i`
for a translation surface in `H(m_1, \dots, m_n)`, `k_i = d_i / 2` for a
half-translation surface in `Q(d_1, \dots, d_n)`.

For an orbifold, where every cone angle is `2\pi/m_i`, Katok's signature
`(g; m_1, \dots, m_r; s)` (*Fuchsian groups*, 1992) is a second form, the
orbifold orders, with `k_i = 1/m_i - 1` and `s = c`.

EXAMPLES:

The signature of the hyperbolic orbifold of a triangle of angles `\pi/2`,
`\pi/8`, `\pi/8`::

    sage: from flatsurf.geometry.signature import SurfaceSignature
    sage: s = SurfaceSignature.from_orbifold_orders(0, (2, 8, 8))
    sage: s
    SurfaceSignature(g=0, cone_orders=(-7/8, -7/8, -1/2))
    sage: s.divisor_euler_characteristic()
    -1/4
    sage: s.geometry()
    'hyperbolic'

The signature of a translation surface in `H(2)`, whose topological and
divisor Euler characteristics differ::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon()
    sage: S.signature()
    SurfaceSignature(g=2, cone_orders=(2,))
    sage: S.signature().euler_characteristic()
    -2
    sage: S.signature().divisor_euler_characteristic()
    0

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

from sage.structure.sage_object import SageObject


class SurfaceSignature(SageObject):
    r"""
    The signature of a compact cone surface: its genus, the orders of its
    cone points, its number of cusps and the traces of its boundary
    components.

    The signature is topological and holds for every curvature, see the
    module documentation.

    INPUT:

    - ``g`` -- a non-negative integer, the genus

    - ``cone_orders`` -- a sequence of rationals (default: ``()``); the
      orders `k > -1`, `k \neq 0`, of the cone points, of angle
      `2\pi (1 + k)`

    - ``cusps`` -- a non-negative integer (default: ``0``); the number of
      cusps

    - ``boundary_traces`` -- a sequence (default: ``()``); for each boundary
      component, the trace of its holonomy in `\mathrm{SO}(2, 1)`, see
      :meth:`boundary_traces`, or ``None`` when the surface is not
      hyperbolic

    EXAMPLES::

        sage: from flatsurf.geometry.signature import SurfaceSignature
        sage: s = SurfaceSignature(0, (-1/2, -7/8, -7/8))
        sage: s
        SurfaceSignature(g=0, cone_orders=(-7/8, -7/8, -1/2))
        sage: s == SurfaceSignature.from_orbifold_orders(0, (8, 2, 8))
        True
        sage: TestSuite(s).run()

    The repr is valid input and shows only the non-empty groups::

        sage: s = SurfaceSignature(1, cusps=2, boundary_traces=(15, 3))
        sage: s
        SurfaceSignature(g=1, cusps=2, boundary_traces=(3, 15))
        sage: SurfaceSignature(g=1, cusps=2, boundary_traces=(3, 15)) == s
        True
        sage: SurfaceSignature(g=2)
        SurfaceSignature(g=2)

    TESTS::

        sage: SurfaceSignature(-1)
        Traceback (most recent call last):
        ...
        ValueError: the genus must be non-negative
        sage: SurfaceSignature(0, cusps=-1)
        Traceback (most recent call last):
        ...
        ValueError: the number of cusps must be non-negative
        sage: SurfaceSignature(0, (-1,))
        Traceback (most recent call last):
        ...
        ValueError: the cone order -1 is not larger than -1
        sage: SurfaceSignature(0, (0,))
        Traceback (most recent call last):
        ...
        ValueError: a cone order must be nonzero
        sage: SurfaceSignature(1/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    The boundary traces are all given or all ``None``::

        sage: SurfaceSignature(0, boundary_traces=(None, 15))
        Traceback (most recent call last):
        ...
        ValueError: the boundary traces must be all None or none of them None
        sage: SurfaceSignature(0, boundary_traces=(None, None))
        SurfaceSignature(g=0, boundary_traces=(None, None))

    """

    def __init__(self, g, cone_orders=(), cusps=0, boundary_traces=()):
        from sage.all import QQ, ZZ

        g = ZZ(g)
        if g < 0:
            raise ValueError("the genus must be non-negative")

        cusps = ZZ(cusps)
        if cusps < 0:
            raise ValueError("the number of cusps must be non-negative")

        orders = []
        for k in cone_orders:
            k = QQ(k)
            if k <= -1:
                raise ValueError(f"the cone order {k} is not larger than -1")
            if k == 0:
                raise ValueError("a cone order must be nonzero")
            orders.append(k)

        traces = tuple(boundary_traces)
        unknown = sum(t is None for t in traces)
        if unknown == 0:
            traces = tuple(sorted(traces))
        elif unknown != len(traces):
            raise ValueError(
                "the boundary traces must be all None or none of them None"
            )

        self._genus = g
        self._cone_orders = tuple(sorted(orders))
        self._cusps = cusps
        self._boundary_traces = traces

    @classmethod
    def from_orbifold_orders(cls, g, periods, cusps=0):
        r"""
        Return the signature of an orbifold given by its orbifold orders, i.e.,
        Katok's form `(g; m_1, \dots, m_r; s)`.

        INPUT:

        - ``g`` -- a non-negative integer, the genus

        - ``periods`` -- a sequence of integers `m_i \ge 2`; the cone points
          of angle `2\pi/m_i`

        - ``cusps`` -- a non-negative integer (default: ``0``); the number
          `s` of cusps

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3, 7))
            SurfaceSignature(g=0, cone_orders=(-6/7, -2/3, -1/2))
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3), cusps=1)
            SurfaceSignature(g=0, cone_orders=(-2/3, -1/2), cusps=1)

        The arguments can be Python integers::

            sage: SurfaceSignature.from_orbifold_orders(int(0), (int(2), int(3)))
            SurfaceSignature(g=0, cone_orders=(-2/3, -1/2))

        TESTS::

            sage: SurfaceSignature.from_orbifold_orders(0, (1, 2))
            Traceback (most recent call last):
            ...
            ValueError: the period 1 is not an integer at least 2

        """
        from sage.all import QQ, ZZ

        orders = []
        for m in periods:
            if m not in ZZ or m < 2:
                raise ValueError(f"the period {m} is not an integer at least 2")
            orders.append(1 / QQ(m) - 1)

        return cls(g, orders, cusps=cusps)

    def genus(self):
        r"""
        Return the genus of this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(2).genus()
            2

        """
        return self._genus

    def cone_orders(self):
        r"""
        Return the orders `k_i` of the cone points, sorted: the cone angles
        are `2\pi (1 + k_i)`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 8, 8)).cone_orders()
            (-7/8, -7/8, -1/2)

        """
        return self._cone_orders

    def cone_angles(self):
        r"""
        Return the angles of the cone points, as multiples of `2\pi`, in the
        order of :meth:`cone_orders`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 8, 8)).cone_angles()
            (1/8, 1/8, 1/2)

        """
        return tuple(1 + k for k in self._cone_orders)

    def cusps(self):
        r"""
        Return the number of cusps of this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (), cusps=3).cusps()
            3

        """
        return self._cusps

    def boundary_components(self):
        r"""
        Return the number of boundary components of this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(0, boundary_traces=(15, 15, 15)).boundary_components()
            3

        """
        return len(self._boundary_traces)

    def boundary_traces(self):
        r"""
        Return the traces of the holonomy of the boundary components in
        `\mathrm{SO}(2, 1)`, sorted, or ``None`` for each of them when the
        surface is not hyperbolic (on a flat surface the holonomy has no such
        invariant).

        A hyperbolic isometry of translation length `\ell` has the trace
        `t = 1 + 2 \cosh \ell` in `\mathrm{SO}(2, 1)`, an element of the base
        ring of the surface. The trace of its lift to `\mathrm{SL}(2,
        \mathbb{R})` is `\tau = \sqrt{t + 1}` and the length is `\ell =
        \operatorname{arccosh}((t - 1)/2)`, which need a square root or a
        transcendental function. A cusp, a parabolic isometry, would have
        `t = 3`.

        EXAMPLES:

        The pair of pants glued from two right-angled regular hexagons has
        three boundary components of length `\operatorname{arccosh} 7`::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: s = SurfaceSignature(0, boundary_traces=(15, 15, 15))
            sage: s.boundary_traces()
            (15, 15, 15)
            sage: t = s.boundary_traces()[0]
            sage: arccosh((t - 1) / 2)
            arccosh(7)

        The trace `t = 4` gives the `\mathrm{SL}(2, \mathbb{R})` trace
        `\sqrt{5}`::

            sage: t = 4
            sage: sqrt(t + 1)
            sqrt(5)

        """
        return self._boundary_traces

    def euler_characteristic(self):
        r"""
        Return the topological Euler characteristic `\chi(S) = 2 - 2g - c -
        b` of this signature, where `c` is the number of cusps and `b` the
        number of boundary components.

        It does not depend on the cone orders: it is the
        ``euler_characteristic()`` of a surface of this signature. See
        :meth:`divisor_euler_characteristic` for the one of Troyanov.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3, 7)).euler_characteristic()
            2
            sage: SurfaceSignature(0, boundary_traces=(15, 15, 15)).euler_characteristic()
            -1
            sage: SurfaceSignature(1, cusps=2, boundary_traces=(15, 3)).euler_characteristic()
            -4

        The topological Euler characteristic of a translation surface, where
        it differs from the divisor one::

            sage: from flatsurf import translation_surfaces, hyperbolic_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: S.signature().euler_characteristic()
            -2
            sage: S.signature().euler_characteristic() == S.euler_characteristic()
            True
            sage: S = hyperbolic_surfaces.regular_octagon()
            sage: S.signature().euler_characteristic()
            -2

        """
        return 2 - 2 * self._genus - self._cusps - len(self._boundary_traces)

    def divisor_euler_characteristic(self):
        r"""
        Return the divisor Euler characteristic `\chi(S, \beta) = 2 - 2g - c -
        b + \sum k_i` of Troyanov of this signature, for the divisor `\beta =
        \sum k_i p_i` of the cone points.

        On an orbifold, where every cone angle is `2\pi/m_i`, it is the
        orbifold Euler characteristic `2 - 2g - \sum (1 - 1/m_i) - c` of
        Thurston. By the Gauss-Bonnet theorem, a hyperbolic surface has area
        `-2 \pi \chi(S, \beta)` and a flat one has `\chi(S, \beta) = 0`, for
        every cone angle and with geodesic boundary.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3, 7)).divisor_euler_characteristic()
            -1/42
            sage: SurfaceSignature(0, [-1/4] * 8).divisor_euler_characteristic()
            0
            sage: SurfaceSignature(0, boundary_traces=(15, 15, 15)).divisor_euler_characteristic()
            -1
            sage: SurfaceSignature(1, cusps=2, boundary_traces=(15, 3)).divisor_euler_characteristic()
            -4

        The octagons, flat and hyperbolic::

            sage: from flatsurf import translation_surfaces, hyperbolic_surfaces
            sage: translation_surfaces.regular_octagon().signature().divisor_euler_characteristic()
            0
            sage: hyperbolic_surfaces.regular_octagon().signature().divisor_euler_characteristic()
            -2

        """
        from sage.all import QQ

        return self.euler_characteristic() + sum(self._cone_orders, QQ(0))

    def area(self):
        r"""
        Return the area of a hyperbolic surface of this signature, in units of
        `2\pi`, i.e., `-\chi(S, \beta)`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(2).area()
            2

        """
        return -self.divisor_euler_characteristic()

    def geometry(self):
        r"""
        Return the geometry of the surfaces of this signature, from the sign
        of :meth:`divisor_euler_characteristic`.

        OUTPUT:

        ``"hyperbolic"``, ``"euclidean"`` or ``"spherical"``, for a negative,
        zero or positive divisor Euler characteristic.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(2).geometry()
            'hyperbolic'
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3, 6)).geometry()
            'euclidean'
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3, 5)).geometry()
            'spherical'

        """
        chi = self.divisor_euler_characteristic()
        if chi < 0:
            return "hyperbolic"
        if chi == 0:
            return "euclidean"
        return "spherical"

    def orbifold_orders(self):
        r"""
        Return this signature as orbifold orders, Katok's form `(g; m_1,
        \dots, m_r; s)`.

        OUTPUT:

        A triple ``(g, periods, s)``: the genus, the sorted tuple of the
        orbifold orders `m_i`, with `k_i = 1/m_i - 1`, and the number of
        cusps. A ``ValueError`` is raised if a cone angle is not `2\pi/m` for
        an integer `m \ge 2`, or if there is a boundary.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(0, (-7/8, -1/2, -7/8)).orbifold_orders()
            (0, (2, 8, 8), 0)
            sage: SurfaceSignature.from_orbifold_orders(0, (2, 3), cusps=1).orbifold_orders()
            (0, (2, 3), 1)

        A cone angle that is not `2\pi/m`, or a boundary, has no orbifold
        orders::

            sage: SurfaceSignature(2, (2,)).orbifold_orders()
            Traceback (most recent call last):
            ...
            ValueError: the cone order 2 is not 1/m - 1 for an integer m
            sage: SurfaceSignature(0, boundary_traces=(15, 15, 15)).orbifold_orders()
            Traceback (most recent call last):
            ...
            ValueError: a signature with boundary has no orbifold orders

        """
        if self._boundary_traces:
            raise ValueError("a signature with boundary has no orbifold orders")

        periods = []
        for k in self._cone_orders:
            angle = 1 + k
            if angle.numerator() != 1 or angle.denominator() < 2:
                raise ValueError(f"the cone order {k} is not 1/m - 1 for an integer m")
            periods.append(angle.denominator())

        return (self._genus, tuple(sorted(periods)), self._cusps)

    def _repr_(self):
        r"""
        Return a printable representation of this signature, which is valid
        input and shows only the non-empty groups.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(2)
            SurfaceSignature(g=2)
            sage: SurfaceSignature(0, (2,), boundary_traces=(None, None))
            SurfaceSignature(g=0, cone_orders=(2,), boundary_traces=(None, None))

        """
        text = f"SurfaceSignature(g={self._genus!r}"
        if self._cone_orders:
            text += f", cone_orders={self._cone_orders!r}"
        if self._cusps:
            text += f", cusps={self._cusps!r}"
        if self._boundary_traces:
            text += f", boundary_traces={self._boundary_traces!r}"
        return text + ")"

    def _key(self):
        r"""
        Return the data that determines this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(1, (1,), boundary_traces=(3,))._key()
            (1, (1,), 0, (3,))

        """
        return (self._genus, self._cone_orders, self._cusps, self._boundary_traces)

    def __eq__(self, other):
        r"""
        Return whether this signature is equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(0, (-1/2,) * 4) == SurfaceSignature.from_orbifold_orders(0, (2, 2, 2, 2))
            True
            sage: SurfaceSignature(0, (-1/2,) * 4) == SurfaceSignature(0, (-1/2,) * 3)
            False
            sage: SurfaceSignature(0) == 0
            False

        """
        if not isinstance(other, SurfaceSignature):
            return False
        return self._key() == other._key()

    def __ne__(self, other):
        r"""
        Return whether this signature is not equal to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(1) != SurfaceSignature(2)
            True

        """
        return not self == other

    def __hash__(self):
        r"""
        Return a hash value for this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: hash(SurfaceSignature(0, (-1/2,) * 4)) == hash(SurfaceSignature.from_orbifold_orders(0, (2, 2, 2, 2)))
            True

        """
        # TODO: The hash of a boundary trace depends on its parent, while
        # equality does not: SurfaceSignature(0, boundary_traces=(15, 15, 15))
        # compares equal to the same signature with traces AA(15), but
        # hash(AA(15)) differs from hash(15). Missing: a hash of the traces
        # that does not depend on their parent, for instance through their
        # minimal polynomials, or one that hashes only the topological part of
        # the key.
        return hash(self._key())
