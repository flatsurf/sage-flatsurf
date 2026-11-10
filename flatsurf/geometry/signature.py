r"""
The signature of a cone surface.

The signature of a compact cone surface, flat or hyperbolic, has three
groups:

- the genus `g` and the cone points, given by their order `k_i`: the cone
  angle is `2\pi (1 + k_i)` with `k_i > -1`, `k_i \neq 0`, a rational number;
- the number `c` of cusps;
- the geodesic boundary components, each given by the trace of its holonomy
  in `\mathrm{SO}(2, 1)`.

The Euler characteristic with a divisor (Troyanov, *Prescribing curvature on
compact surfaces with conical singularities*, 1991) is
`\chi(S, \beta) = 2 - 2g + \sum k_i - c - b` where `b` is the number of
boundary components. By the Gauss-Bonnet theorem, it is zero for a flat
surface, and the area of a hyperbolic surface is `-2 \pi \chi(S, \beta)`.

The cone orders of a flat surface are the data of its stratum: `k_i = m_i`
for a translation surface in `H(m_1, \dots, m_n)`, `k_i = d_i / 2` for a
half-translation surface in `Q(d_1, \dots, d_n)`.

For an orbifold, where every cone angle is `2\pi/m_i`, Katok's signature
`(g; m_1, \dots, m_r; s)` (*Fuchsian groups*, 1992) is a second form, with
`k_i = 1/m_i - 1` and `s = c`.

EXAMPLES:

The signature of the hyperbolic orbifold of a triangle of angles `\pi/2`,
`\pi/8`, `\pi/8`::

    sage: from flatsurf.geometry.signature import SurfaceSignature
    sage: s = SurfaceSignature.from_katok(0, (2, 8, 8))
    sage: s
    SurfaceSignature(genus=0, cone_orders=(-7/8, -7/8, -1/2))
    sage: s.euler_characteristic()
    -1/4
    sage: s.geometry()
    'hyperbolic'

The signature of a translation surface in `H(2)`::

    sage: from flatsurf import translation_surfaces
    sage: translation_surfaces.regular_octagon().signature()
    SurfaceSignature(genus=2, cone_orders=(2,))

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

    INPUT:

    - ``genus`` -- a non-negative integer

    - ``cone_orders`` -- a sequence of rationals (default: ``()``); the
      orders `k > -1`, `k \neq 0`, of the cone points, of angle
      `2\pi (1 + k)`

    - ``cusps`` -- a non-negative integer (default: ``0``); the number of
      cusps

    - ``boundary_traces`` -- a sequence (default: ``()``); for each boundary
      component, the trace of its holonomy in `\mathrm{SO}(2, 1)`, see
      :meth:`boundary_traces`, or ``None`` when there is none

    EXAMPLES::

        sage: from flatsurf.geometry.signature import SurfaceSignature
        sage: s = SurfaceSignature(0, (-1/2, -7/8, -7/8))
        sage: s
        SurfaceSignature(genus=0, cone_orders=(-7/8, -7/8, -1/2))
        sage: s == SurfaceSignature.from_katok(0, (8, 2, 8))
        True
        sage: TestSuite(s).run()

    The repr is valid input::

        sage: s = SurfaceSignature(1, cusps=2, boundary_traces=(15, 3))
        sage: s
        SurfaceSignature(genus=1, cone_orders=(), cusps=2, boundary_traces=(3, 15))
        sage: SurfaceSignature(genus=1, cone_orders=(), cusps=2, boundary_traces=(3, 15)) == s
        True

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

    """

    def __init__(self, genus, cone_orders=(), cusps=0, boundary_traces=()):
        from sage.all import QQ, ZZ

        genus = ZZ(genus)
        if genus < 0:
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
        if all(t is not None for t in traces):
            traces = tuple(sorted(traces))

        self._genus = genus
        self._cone_orders = tuple(sorted(orders))
        self._cusps = cusps
        self._boundary_traces = traces

    @classmethod
    def from_katok(cls, g, periods, cusps=0):
        r"""
        Return the signature of an orbifold given in Katok's form `(g; m_1,
        \dots, m_r; s)`.

        INPUT:

        - ``g`` -- a non-negative integer, the genus

        - ``periods`` -- a sequence of integers `m_i \ge 2`; the cone points
          of angle `2\pi/m_i`

        - ``cusps`` -- a non-negative integer (default: ``0``); the number
          `s` of cusps

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(0, (2, 3, 7))
            SurfaceSignature(genus=0, cone_orders=(-6/7, -2/3, -1/2))
            sage: SurfaceSignature.from_katok(0, (2, 3), cusps=1)
            SurfaceSignature(genus=0, cone_orders=(-2/3, -1/2), cusps=1)

        TESTS::

            sage: SurfaceSignature.from_katok(0, (1, 2))
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
            sage: SurfaceSignature.from_katok(2, ()).genus()
            2

        """
        return self._genus

    def cone_orders(self):
        r"""
        Return the orders `k_i` of the cone points, sorted: the cone angles
        are `2\pi (1 + k_i)`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(0, (2, 8, 8)).cone_orders()
            (-7/8, -7/8, -1/2)

        """
        return self._cone_orders

    def cone_angles(self):
        r"""
        Return the angles of the cone points, as multiples of `2\pi`, in the
        order of :meth:`cone_orders`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(0, (2, 8, 8)).cone_angles()
            (1/8, 1/8, 1/2)

        """
        return tuple(1 + k for k in self._cone_orders)

    def cusps(self):
        r"""
        Return the number of cusps of this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(0, (), cusps=3).cusps()
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
        holonomy has no such invariant (on a flat surface).

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
        Return the Euler characteristic with a divisor `\chi(S, \beta) = 2 -
        2g + \sum k_i - c - b` of this signature.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(0, (2, 3, 7)).euler_characteristic()
            -1/42
            sage: SurfaceSignature(0, [-1/4] * 8).euler_characteristic()
            0

        """
        from sage.all import QQ

        return (
            2
            - 2 * self._genus
            + sum(self._cone_orders, QQ(0))
            - self._cusps
            - len(self._boundary_traces)
        )

    def area(self):
        r"""
        Return the area of a hyperbolic surface of this signature, in units of
        `2\pi`, i.e., `-\chi(S, \beta)`.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(2, ()).area()
            2

        """
        return -self.euler_characteristic()

    def geometry(self):
        r"""
        Return the geometry of the surfaces of this signature, from the sign
        of :meth:`euler_characteristic`.

        OUTPUT:

        ``"hyperbolic"``, ``"euclidean"`` or ``"spherical"``, for a negative,
        zero or positive Euler characteristic.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature.from_katok(2, ()).geometry()
            'hyperbolic'
            sage: SurfaceSignature.from_katok(0, (2, 3, 6)).geometry()
            'euclidean'
            sage: SurfaceSignature.from_katok(0, (2, 3, 5)).geometry()
            'spherical'

        """
        chi = self.euler_characteristic()
        if chi < 0:
            return "hyperbolic"
        if chi == 0:
            return "euclidean"
        return "spherical"

    def katok(self):
        r"""
        Return this signature in Katok's orbifold form `(g; m_1, \dots, m_r;
        s)`.

        OUTPUT:

        A triple ``(g, periods, s)``: the genus, the sorted tuple of the
        periods `m_i`, with `k_i = 1/m_i - 1`, and the number of cusps. A
        ``ValueError`` is raised if a cone angle is not `2\pi/m` for an
        integer `m \ge 2`, or if there is a boundary.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(0, (-7/8, -1/2, -7/8)).katok()
            (0, (2, 8, 8), 0)
            sage: SurfaceSignature.from_katok(0, (2, 3), cusps=1).katok()
            (0, (2, 3), 1)

        A cone angle that is not `2\pi/m`, or a boundary, has no Katok
        form::

            sage: SurfaceSignature(2, (2,)).katok()
            Traceback (most recent call last):
            ...
            ValueError: the cone order 2 is not 1/m - 1 for an integer m
            sage: SurfaceSignature(0, boundary_traces=(15, 15, 15)).katok()
            Traceback (most recent call last):
            ...
            ValueError: a signature with boundary has no Katok form

        """
        if self._boundary_traces:
            raise ValueError("a signature with boundary has no Katok form")

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
        input.

        EXAMPLES::

            sage: from flatsurf.geometry.signature import SurfaceSignature
            sage: SurfaceSignature(2)
            SurfaceSignature(genus=2, cone_orders=())
            sage: SurfaceSignature(0, (2,), boundary_traces=(None, None))
            SurfaceSignature(genus=0, cone_orders=(2,), boundary_traces=(None, None))

        """
        text = (
            f"SurfaceSignature(genus={self._genus!r}, cone_orders={self._cone_orders!r}"
        )
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
            sage: SurfaceSignature(0, (-1/2,) * 4) == SurfaceSignature.from_katok(0, (2, 2, 2, 2))
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
            sage: hash(SurfaceSignature(0, (-1/2,) * 4)) == hash(SurfaceSignature.from_katok(0, (2, 2, 2, 2)))
            True

        """
        return hash(self._key())
