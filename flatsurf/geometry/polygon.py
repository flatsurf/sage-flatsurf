r"""
Polygons embedded in the plane `\mathbb{R}^2`.

The emphasis is mostly on convex polygons but there is some limited support for
non-convex polygons.

EXAMPLES::

    sage: from flatsurf.geometry.polygon import Polygon

    sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
    sage: p = Polygon(edges=[(1,0), (-sqrt2,1+sqrt2), (sqrt2-1,-1-sqrt2)])
    sage: p
    Polygon(vertices=[(0, 0), (1, 0), (-sqrt2 + 1, sqrt2 + 1)])

    sage: M = MatrixSpace(K,2)
    sage: m = M([[1,1+sqrt2],[0,1]])
    sage: m * p
    Polygon(vertices=[(0, 0), (1, 0), (sqrt2 + 4, sqrt2 + 1)])
"""

# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020-2023 Julian Rüth
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

from sage.all import (
    cached_method,
    Parent,
    QQ,
    matrix,
    vector,
)

from sage.structure.element import Element
from sage.structure.sequence import Sequence

from flatsurf.geometry.subfield import (
    number_field_elements_from_algebraics,
)

from flatsurf.geometry.categories import EuclideanPolygons
from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
    EuclideanPolygonsWithAngles as EuclideanPolygonsWithAnglesCategory,
)


class EuclideanPolygonPoint(Element):
    r"""
    A point in a polygon.

    EXAMPLES::

        sage: from flatsurf import polygons
        sage: s = polygons.square()

        sage: s.an_element()
        (0, 0)

    TESTS::

        sage: p = s.an_element()

        sage: from flatsurf.geometry.polygon import EuclideanPolygonPoint
        sage: isinstance(p, EuclideanPolygonPoint)
        True

        sage: TestSuite(p).run()

    """

    def __init__(self, parent, xy):
        self._xy = xy

        super().__init__(parent)

    def position(self):
        r"""
        Describe the position of this point in the polygon.

        OUTPUT:

        A :class:`PolygonPosition` object.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()

            sage: p = s.an_element()
            sage: p
            (0, 0)

            sage: p.position()
            point positioned on vertex 0 of polygon

        .. SEEALSO::

            :meth:`~.categories.euclidean_polygons.EuclideanPolygons.ParentMethods.get_point_position`

        """
        return self.parent().get_point_position(self._xy)

    def _repr_(self):
        r"""
        Return a printable representation of this point.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()

            sage: p = s.an_element()
            sage: p
            (0, 0)

        """
        return repr(self._xy)

    def __eq__(self, other):
        r"""
        Return whether this point is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()

            sage: s.an_element() == s.an_element()
            True

            sage: t = polygons.square()

            sage: s.an_element() == t.an_element()
            True

        """
        if not isinstance(other, EuclideanPolygonPoint):
            return False

        return self.parent() == other.parent() and self._xy == other._xy

    def __hash__(self):
        r"""
        Return a hash value of this point that is compatible with
        :meth:`_richcmp_`.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()

            sage: hash(s.an_element()) == hash(s.an_element())
            True

        """
        return hash(self._xy)


class PolygonPosition:
    r"""
    Describes the position of a point within or outside of a polygon.
    """

    # Position Types:
    OUTSIDE = 0
    INTERIOR = 1
    EDGE_INTERIOR = 2
    VERTEX = 3

    def __init__(self, position_type, edge=None, vertex=None):
        self._position_type = position_type
        if self.is_vertex():
            if vertex is None:
                raise ValueError(
                    "Constructed vertex position with no specified vertex."
                )
            self._vertex = vertex
        if self.is_in_edge_interior():
            if edge is None:
                raise ValueError("Constructed edge position with no specified edge.")
            self._edge = edge

    def __repr__(self):
        if self.is_outside():
            return "point positioned outside polygon"
        if self.is_in_interior():
            return "point positioned in interior of polygon"
        if self.is_in_edge_interior():
            return (
                "point positioned on interior of edge "
                + str(self._edge)
                + " of polygon"
            )
        return "point positioned on vertex " + str(self._vertex) + " of polygon"

    def is_outside(self):
        return self._position_type == PolygonPosition.OUTSIDE

    def is_inside(self):
        r"""
        Return true if the position is not outside the closure of the polygon
        """
        return bool(self._position_type)

    def is_in_interior(self):
        return self._position_type == PolygonPosition.INTERIOR

    def is_in_boundary(self):
        r"""
        Return true if the position is in the boundary of the polygon
        (either the interior of an edge or a vertex).
        """
        return (
            self._position_type == PolygonPosition.EDGE_INTERIOR
            or self._position_type == PolygonPosition.VERTEX
        )

    def is_in_edge_interior(self):
        return self._position_type == PolygonPosition.EDGE_INTERIOR

    def is_vertex(self):
        return self._position_type == PolygonPosition.VERTEX

    def get_position_type(self):
        return self._position_type

    def get_edge(self):
        if not self.is_in_edge_interior():
            raise ValueError("Asked for edge when not in edge interior.")
        return self._edge

    def get_vertex(self):
        if not self.is_vertex():
            raise ValueError("Asked for vertex when not a vertex.")
        return self._vertex


class EuclideanPolygon(Parent):
    r"""
    A (possibly non-convex) simple polygon in the plane `\mathbb{R}^2`.

    EXAMPLES::

        sage: from flatsurf import polygons, Polygon
        sage: s = polygons.square()
        sage: s
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

    TESTS::

        sage: from flatsurf.geometry.polygon import EuclideanPolygon
        sage: isinstance(s, EuclideanPolygon)
        True

        sage: TestSuite(s).run()

    """

    Element = EuclideanPolygonPoint

    def __init__(self, base_ring, vertices, category=None):
        V = base_ring**2
        self._v = tuple(map(V, vertices))
        for vv in self._v:
            vv.set_immutable()
        if category is None:
            category = EuclideanPolygons(base_ring)

        category &= EuclideanPolygons(base_ring)

        super().__init__(base_ring, category=category)

        if "Convex" not in category.axioms() and self.is_convex():
            self._refine_category_(category.Convex())

        # The category is not refined automatically to the WithAngles()
        # subcategory since computation of angles can be very costly.
        # The category gets further refined when angles() is invoked.

    def _an_element_(self):
        r"""
        Return a point of this polygon.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()
            sage: s.an_element()
            (0, 0)

        """
        return self(self.vertices()[0])

    def parent(self):
        r"""
        Return the category this polygon belongs to.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()
            sage: s.parent()
            doctest:warning
            ...
            UserWarning: parent() of a polygon has been deprecated and will be removed in a future version of sage-flatsurf; use category() instead
            Category of convex simple euclidean rectangles over Rational Field

        Note that the parent may change during the lifetime of a polygon as
        more of its features are discovered::

            sage: from flatsurf import Polygon
            sage: p = Polygon(vertices=[(0, 0), (1, 0), (1, 1)])
            sage: p.parent()
            Category of convex simple euclidean polygons over Rational Field
            sage: p.angles()
            (1/8, 1/4, 1/8)
            sage: p.parent()
            Category of convex simple euclidean triangles with angles (1/8, 1/4, 1/8) over Rational Field

        """
        import warnings

        warnings.warn(
            "parent() of a polygon has been deprecated and will be removed in a future version of sage-flatsurf; use category() instead"
        )

        return self.category()

    @cached_method
    def __hash__(self):
        return hash(self._v)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf import polygons, Polygon
            sage: p1 = polygons.square()
            sage: p2 = Polygon(edges=[(1,0),(0,1),(-1,0),(0,-1)], base_ring=QQbar)
            sage: p1 == p2
            True

            sage: p3 = Polygon(edges=[(2,0),(-1,1),(-1,-1)])
            sage: p1 == p3
            False

        TESTS::

            sage: from flatsurf import Polygon, polygons
            sage: p1 = polygons.square()
            sage: p2 = Polygon(edges=[(1,0),(0,1),(-1,0),(0,-1)], base_ring=QQbar)
            sage: p1 != p2
            False

            sage: p3 = Polygon(edges=[(2,0),(-1,1),(-1,-1)])
            sage: p1 != p3
            True
        """
        if not isinstance(other, EuclideanPolygon):
            return False

        return self._v == other._v

    def cmp(self, other):
        r"""
        Implement a total order on polygons
        """
        if not isinstance(other, EuclideanPolygon):
            raise TypeError("__cmp__ only implemented for ConvexPolygons")
        if not self.base_ring() == other.base_ring():
            raise ValueError(
                "__cmp__ only implemented for ConvexPolygons defined over the same base_ring"
            )
        sign = len(self.vertices()) - len(other.vertices())
        if sign > 0:
            return 1
        if sign < 0:
            return -1
        sign = self.area() - other.area()
        if sign > self.base_ring().zero():
            return 1
        if sign < self.base_ring().zero():
            return -1
        for v in range(1, len(self.vertices())):
            p = self.vertex(v)
            q = other.vertex(v)
            sign = p[0] - q[0]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
            sign = p[1] - q[1]
            if sign > self.base_ring().zero():
                return 1
            if sign < self.base_ring().zero():
                return -1
        return 0

    def translate(self, u):
        r"""
        Return a copy of this polygon that has been translated by ``u``.

        TESTS::

            sage: from flatsurf import Polygon
            sage: Polygon(vertices=[(0,0), (2,0), (1,1)]).translate((3,-2))
            Polygon(vertices=[(3, -2), (5, -2), (4, -1)])

        """
        u = (self.base_ring() ** 2)(u)
        return Polygon(
            base_ring=self.base_ring(),
            vertices=[u + v for v in self._v],
            check=False,
            category=self.category(),
        )

    def change_ring(self, ring):
        r"""
        Return a copy of this polygon whose vertices have coordinates over the
        base ring ``ring``.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2)**(1/2))
            sage: S.change_ring(K)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
            sage: S.change_ring(K).base_ring()
            Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.4142...

        """
        if ring is self.base_ring():
            return self

        return Polygon(
            base_ring=ring, vertices=self._v, category=self.category().change_ring(ring)
        )

    def is_strictly_convex(self):
        r"""
        Return whether this polygon is strictly convex.

        EXAMPLES::

            sage: from flatsurf import Polygon
            sage: Polygon(vertices=[(0,0), (1,0), (1,1)]).is_strictly_convex()
            doctest:warning
            ...
            UserWarning: is_strictly_convex() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex(strict=True) instead
            True
            sage: Polygon(vertices=[(0,0), (1,0), (2,0), (1,1)]).is_strictly_convex()
            False

        TESTS::

            sage: Polygon(vertices=[(0, 0), (1, 1/2), (2, 0), (1, 1)]).is_strictly_convex()
            False

        """
        import warnings

        warnings.warn(
            "is_strictly_convex() has been deprecated and will be removed in a future version of sage-flatsurf; use is_convex(strict=True) instead"
        )

        return self.is_convex(strict=True)

    def num_edges(self):
        r"""
        Return the number of edges of this polygon.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: S.num_edges()
            doctest:warning
            ...
            UserWarning: num_edges() has been deprecated and will be removed in a future version of sage-flatsurf; use len(vertices()) instead
            4

        """
        import warnings

        warnings.warn(
            "num_edges() has been deprecated and will be removed in a future version of sage-flatsurf; use len(vertices()) instead"
        )

        return len(self.vertices())

    def _repr_(self):
        r"""
        Return a printable representation of this polygon.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: S = polygons.square()
            sage: S
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        """
        return f"Polygon(vertices={repr(list(self.vertices()))})"

    def vertices(self, translation=None, marked_vertices=True):
        r"""
        Return the vertices of this polygon in counterclockwise order as
        vectors in the real plane.

        INPUT:

        - ``marked_vertices`` -- a boolean (default: ``True``); whether to
          include vertices with a π angle in the output.

        EXAMPLES::

            sage: from flatsurf import polygons
            sage: s = polygons.square()
            sage: s.vertices()
            ((0, 0), (1, 0), (1, 1), (0, 1))

        """
        if translation is not None:
            import warnings

            warnings.warn(
                "the translation keyword of vertices() has been deprecated and will be removed in a future version of sage-flatsurf; use translate().vertices() instead"
            )

            return self.translate(translation).vertices(marked_vertices=marked_vertices)

        if not marked_vertices:
            return tuple(
                vertex
                for (vertex, slope) in zip(self._v, self.slopes(relative=True))
                if slope[1] != 0
            )

        return self._v

    def __iter__(self):
        import warnings

        warnings.warn(
            "iterating over the vertices of a polygon implicitly has been deprecated, this functionality will be removed in a future version of sage-flatsurf; iterate over vertices() instead"
        )
        return iter(self.vertices())


class PolygonsConstructor:
    def square(self, side=1, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons.square()
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
            sage: polygons.square(base_ring=QQbar).category()
            Category of convex simple euclidean rectangles over Algebraic Field

        """
        return self.rectangle(side, side, **kwds)

    def rectangle(self, width, height, **kwds):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: polygons.rectangle(1,2)
            Polygon(vertices=[(0, 0), (1, 0), (1, 2), (0, 2)])

            sage: K.<sqrt2> = QuadraticField(2)
            sage: polygons.rectangle(1,sqrt2)
            Polygon(vertices=[(0, 0), (1, 0), (1, sqrt2), (0, sqrt2)])
            sage: _.category()
            Category of convex simple euclidean rectangles over Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?

        """
        if width <= 0:
            raise ValueError("width must be positive")

        if height <= 0:
            raise ValueError("height must be positive")

        if not kwds:
            # No need to verify that the edges and the angles are consistent.
            kwds = {"check": False}

        return Polygon(
            edges=[(width, 0), (0, height), (-width, 0), (0, -height)],
            angles=(1, 1, 1, 1),
            **kwds,
        )

    def triangle(self, a, b, c):
        """
        Return the triangle with angles a*pi/N,b*pi/N,c*pi/N where N=a+b+c.

        INPUT:

        - ``a``, ``b``, ``c`` -- integers

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons
            sage: T = polygons.triangle(3,4,5)
            sage: T
            Polygon(vertices=[(0, 0), (1, 0), (-1/2*c0 + 3/2, -1/2*c0 + 3/2)])
            sage: T.base_ring()
            Number Field in c0 with defining polynomial x^2 - 3 with c0 = 1.732050807568878?

            sage: polygons.triangle(1,2,3).angles()
            (1/12, 1/6, 1/4)

        Some fairly complicated examples::

            sage: polygons.triangle(1, 15, 21)  # long time (2s)
            Polygon(vertices=[(0, 0),
                              (1, 0),
                              (1/2*c^34 - 17*c^32 + 264*c^30 - 2480*c^28 + 15732*c^26 - 142481/2*c^24 + 237372*c^22 - 1182269/2*c^20 +
                               1106380*c^18 - 1552100*c^16 + 3229985/2*c^14 - 2445665/2*c^12 + 654017*c^10 - 472615/2*c^8 + 107809/2*c^6 - 13923/2*c^4 + 416*c^2 - 6,
                               -1/2*c^27 + 27/2*c^25 - 323/2*c^23 + 1127*c^21 - 10165/2*c^19 + 31009/2*c^17 - 65093/2*c^15 + 46911*c^13 - 91091/2*c^11 + 57355/2*c^9 - 10994*c^7 + 4621/2*c^5 - 439/2*c^3 + 6*c)])

            sage: polygons.triangle(2, 13, 26)  # long time (3s)
            Polygon(vertices=[(0, 0),
                              (1, 0),
                              (1/2*c^30 - 15*c^28 + 405/2*c^26 - 1625*c^24 + 8625*c^22 - 31878*c^20 + 168245/2*c^18 - 159885*c^16 + 218025*c^14 - 209950*c^12 + 138567*c^10 - 59670*c^8 + 15470*c^6 - 2100*c^4 + 225/2*c^2 - 1/2,
                               -1/2*c^39 + 19*c^37 - 333*c^35 + 3571*c^33 - 26212*c^31 + 139593*c^29 - 557844*c^27 + 1706678*c^25 - 8085237/2*c^23 + 7449332*c^21 -
                               10671265*c^19 + 11812681*c^17 - 9983946*c^15 + 6317339*c^13 - 5805345/2*c^11 + 1848183/2*c^9 - 378929/2*c^7 + 44543/2*c^5 - 2487/2*c^3 + 43/2*c)])
        """
        return Polygon(angles=[a, b, c], check=False)

    @staticmethod
    def regular_ngon(n, field=None):
        r"""
        Return a regular n-gon with unit length edges, first edge horizontal, and other vertices lying above this edge.

        Assuming field is None (by default) the polygon is defined over a NumberField (the minimal number field determined by n).
        Otherwise you can set field equal to AA to define the polygon over the Algebraic Reals. Other values for the field
        parameter will result in a ValueError.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import polygons

            sage: p = polygons.regular_ngon(17)
            sage: p
            Polygon(vertices=[(0, 0), (1, 0), ..., (-1/2*a^14 + 15/2*a^12 - 45*a^10 + 275/2*a^8 - 225*a^6 + 189*a^4 - 70*a^2 + 15/2, 1/2*a)])

            sage: polygons.regular_ngon(3,field=AA)
            Polygon(vertices=[(0, 0), (1, 0), (1/2, 0.866025403784439?)])
        """
        # The code below crashes for n=4!
        if n == 4:
            return polygons.square(QQ(1), base_ring=field)

        from sage.rings.qqbar import QQbar

        c = QQbar.zeta(n).real()
        s = QQbar.zeta(n).imag()

        if field is None:
            field, (c, s) = number_field_elements_from_algebraics((c, s))
        cn = field.one()
        sn = field.zero()
        edges = [(cn, sn)]
        for _ in range(n - 1):
            cn, sn = c * cn - s * sn, c * sn + s * cn
            edges.append((cn, sn))

        ngon = Polygon(base_ring=field, edges=edges)
        ngon._refine_category_(ngon.category().WithAngles([1] * n))
        return ngon

    @staticmethod
    def right_triangle(angle, leg0=None, leg1=None, hypotenuse=None):
        r"""
        Return a right triangle in a number field with an angle of pi*angle.

        You can specify the length of the first leg (``leg0``), the second leg (``leg1``),
        or the ``hypotenuse``.

        EXAMPLES::

            sage: from flatsurf import polygons

            sage: P = polygons.right_triangle(1/3, 1)
            sage: P
            Polygon(vertices=[(0, 0), (1, 0), (1, a)])
            sage: P.base_ring()
            Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?

            sage: polygons.right_triangle(1/4,1)
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])
            sage: polygons.right_triangle(1/4,1).base_ring()
            Rational Field
        """
        from sage.rings.qqbar import QQbar

        angle = QQ(angle)
        if angle <= 0 or angle > QQ((1, 2)):
            raise ValueError("angle must be in ]0,1/2]")

        z = QQbar.zeta(2 * angle.denom()) ** angle.numerator()
        c = z.real()
        s = z.imag()

        nargs = (leg0 is not None) + (leg1 is not None) + (hypotenuse is not None)

        if nargs == 0:
            leg0 = 1
        elif nargs > 1:
            raise ValueError("only one length can be specified")

        if leg0 is not None:
            c, s = leg0 * c / c, leg0 * s / c
        elif leg1 is not None:
            c, s = leg1 * c / s, leg1 * s / s
        elif hypotenuse is not None:
            c, s = hypotenuse * c, hypotenuse * s

        field, (c, s) = number_field_elements_from_algebraics((c, s))

        return Polygon(
            base_ring=field, edges=[(c, field.zero()), (field.zero(), s), (-c, -s)]
        )

    def __call__(self, *args, **kwargs):
        r"""
        EXAMPLES::

            sage: from flatsurf import polygons

            sage: polygons((1,0),(0,1),(-1,0),(0,-1))
            doctest:warning
            ...
            UserWarning: calling Polygon() with positional arguments has been deprecated and will not be supported in a future version of sage-flatsurf; use edges= or vertices= explicitly instead
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
            sage: polygons((1,0),(0,1),(-1,0),(0,-1), ring=QQbar)
            doctest:warning
            ...
            UserWarning: ring has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use base_ring instead
            Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])
            sage: _.category()
            Category of convex simple euclidean polygons over Algebraic Field

            sage: polygons(vertices=[(0,0), (1,0), (0,1)])
            Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

            sage: polygons(edges=[(2,0),(-1,1),(-1,-1)], base_point=(3,3))
            doctest:warning
            ...
            UserWarning: base_point has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use .translate() on the resulting polygon instead
            Polygon(vertices=[(3, 3), (5, 3), (4, 4)])
            sage: polygons(vertices=[(0,0),(2,0),(1,1)], base_point=(3,3))
            Polygon(vertices=[(3, 3), (5, 3), (4, 4)])

            sage: polygons(angles=[1,1,1,2], length=1)
            doctest:warning
            ...
            UserWarning: length has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use lengths instead
            Polygon(vertices=[(0, 0), (1, 0), (-1/2*c^2 + 5/2, 1/2*c), (-1/2*c^2 + 2, 1/2*c^3 - 3/2*c)])
            sage: polygons(angles=[1,1,1,2], length=2)
            Polygon(vertices=[(0, 0), (2, 0), (-c^2 + 5, c), (-c^2 + 4, c^3 - 3*c)])
            sage: polygons(angles=[1,1,1,2], length=AA(2)**(1/2))  # tol 1e-9
            Polygon(vertices=[(0, 0), (1.414213562373095?, 0), (0.9771975379242739?, 1.344997023927915?), (0.270090756737727?, 0.831253875554907?)])

            sage: polygons(angles=[1]*5).angles()
            (3/10, 3/10, 3/10, 3/10, 3/10)
            sage: polygons(angles=[1]*8).angles()
            (3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8, 3/8)

            sage: P = polygons(angles=[1,1,3,3], lengths=[3,1])
            sage: P.angles()
            (1/8, 1/8, 3/8, 3/8)
            sage: e0 = P.edge(0); assert e0[0]**2 + e0[1]**2 == 3**2
            sage: e1 = P.edge(1); assert e1[0]**2 + e1[1]**2 == 1

            sage: polygons(angles=[1, 1, 1, 2])
            Polygon(vertices=[(0, 0), (1/10*c^3 + c^2 - 1/5*c - 3, 0), (1/20*c^3 + 1/2*c^2 - 1/20*c - 3/2, 1/20*c^2 + 1/2*c), (1/2*c^2 - 3/2, 1/2*c)])

            sage: polygons(angles=[1,1,1,8])
            Polygon(vertices=[(0, 0), (c^6 - 6*c^4 + 8*c^2 + 3, 0), (1/2*c^4 - 3*c^2 + 9/2, 1/2*c^9 - 9/2*c^7 + 13*c^5 - 11*c^3 - 3*c), (1/2*c^6 - 7/2*c^4 + 7*c^2 - 3, 1/2*c^9 - 5*c^7 + 35/2*c^5 - 49/2*c^3 + 21/2*c)])
            sage: polygons(angles=[1,1,1,8], lengths=[1, 1])
            Polygon(vertices=[(0, 0), (1, 0), (-1/2*c^4 + 2*c^2, 1/2*c^7 - 7/2*c^5 + 7*c^3 - 7/2*c), (1/2*c^6 - 7/2*c^4 + 13/2*c^2 - 3/2, 1/2*c^9 - 9/2*c^7 + 27/2*c^5 - 29/2*c^3 + 5/2*c)])

        TESTS::

            sage: from itertools import product
            sage: for a,b,c in product(range(1,5), repeat=3):  # long time (1.5s)
            ....:     if gcd([a,b,c]) != 1:
            ....:         continue
            ....:     T = polygons(angles=[a,b,c])
            ....:     D = 2*(a+b+c)
            ....:     assert T.angles() == (a/D, b/D, c/D)

        """
        import warnings

        warnings.warn(
            "calling polygons() has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead"
        )
        return Polygon(*args, **kwargs)


polygons = PolygonsConstructor()


def ConvexPolygons(base_ring):
    r"""
    EXAMPLES::

        sage: from flatsurf import ConvexPolygons
        sage: P = ConvexPolygons(QQ)
        doctest:warning
        ...
        UserWarning: ConvexPolygons() has been deprecated and will be removed from a future version of sage-flatsurf; use Polygon() to create polygons.
        If you really need the category of convex polygons over a ring use EuclideanPolygons(ring).Simple().Convex() instead.
        sage: P(vertices=[(0, 0), (1, 0), (0, 1)])
        doctest:warning
        ...
        UserWarning: ConvexPolygons(…)(…) has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon() instead
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

    """
    import warnings

    warnings.warn(
        "ConvexPolygons() has been deprecated and will be removed from a future version of sage-flatsurf; use Polygon() to create polygons. "
        "If you really need the category of convex polygons over a ring use EuclideanPolygons(ring).Simple().Convex() instead."
    )
    return EuclideanPolygons(base_ring).Simple().Convex()


def Polygon(
    *args,
    vertices=None,
    edges=None,
    angles=None,
    lengths=None,
    base_ring=None,
    category=None,
    check=True,
    **kwds,
):
    r"""
    Return a polygon from the given ``vertices``, ``edges``, or ``angles``.

    INPUT:

    - ``vertices`` -- a sequence of vertices or ``None`` (default: ``None``); the
      vertices of the polygon as points in the real plane

    - ``edges`` -- a sequence of vectors or ``None`` (default: ``None``); the
      vectors connecting the vertices of the polygon

    - ``angles`` -- a sequence of numbers that prescribe the inner angles of
      the polygon or ``None`` (default: ``None``); the angles are rescaled so
      that their sum matches the sum of the angles in an ngon.

    - ``lengths`` -- a sequence of numbers that prescribe the lengths of the
      edges of the polygon or ``None`` (default: ``None``)

    - ``base_ring`` -- a ring or ``None`` (default: ``None``); the ring over
      which the polygon will be defined

    - ``category`` -- a category or ``None`` (default: ``None``); the category
      the polygon will be in (further refined from the features of the polygon
      that are found during the construction).

    - ``check`` -- a boolean (default: ``True``); whether to check the
      consistency of the parameters or blindly trust them. Setting this to
      ``False`` allows creation of degenerate polygons in some cases. While
      they might be somewhat functional, no guarantees are made about such
      polygons.

    EXAMPLES:

    A right triangle::

        sage: from flatsurf import Polygon
        sage: Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

    A right triangle that is not based at the origin::

        sage: Polygon(vertices=[(1, 0), (2, 0), (1, 1)])
        Polygon(vertices=[(1, 0), (2, 0), (1, 1)])

    A right triangle at the origin, specified by giving the edge vectors::

        sage: Polygon(edges=[(1, 0), (-1, 1), (0, -1)])
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

    When redundant information is given, it is checked for consistency::

        sage: Polygon(vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])
        sage: Polygon(vertices=[(1, 0), (2, 0), (1, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
        Polygon(vertices=[(1, 0), (2, 0), (1, 1)])
        sage: Polygon(vertices=[(0, 0), (2, 0), (1, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
        Traceback (most recent call last):
        ...
        ValueError: vertices and edges are not compatible

    Polygons given by edges must be closed (in particular we do not add an edge
    automatically to close things up since this is often not what the user
    wanted)::

        sage: Polygon(edges=[(1, 0), (0, 1), (1, 1)])
        Traceback (most recent call last):
        ...
        ValueError: polygon not closed

    A polygon with prescribed angles::

        sage: Polygon(angles=[2, 1, 1])
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

    Again, if vertices and edges are also specified, they must be compatible
    with the angles::

        sage: Polygon(angles=[2, 1, 1], vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
        Polygon(vertices=[(0, 0), (1, 0), (0, 1)])

        sage: Polygon(angles=[1, 2, 3], vertices=[(0, 0), (1, 0), (0, 1)], edges=[(1, 0), (-1, 1), (0, -1)])
        Traceback (most recent call last):
        ...
        ValueError: polygon does not have the prescribed angles

    When angles are specified, side lengths can also be prescribed::

        sage: Polygon(angles=[1, 1, 1], lengths=[1, 1, 1])
        Polygon(vertices=[(0, 0), (1, 0), (1/2, 1/2*c)])

    The function will deduce lengths if one or two are missing::

        sage: Polygon(angles=[1, 1, 1, 1], lengths=[1, 1, 1])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(angles=[1, 1, 1, 1], lengths=[1, 1])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(angles=[1, 1, 1, 1], lengths=[1])
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot construct a quadrilateral from 4 angles and 2 vertices

    Equally, we deduce vertices or edges::

        sage: Polygon(angles=[1, 1, 1, 1], vertices=[(0, 0), (1, 0), (1, 1)])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(angles=[1, 1, 1, 1], edges=[(1, 0), (0, 1)])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

    When the angles are incompatible with the data, an error is reported (that
    might be somewhat cryptic at times)::

        sage: Polygon(angles=[1, 1, 1, 1], edges=[(1, 0), (0, 1), (1, 2)])
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot recover a rational angle from these numerical results

    When lengths are given in addition to vertices or edges, they are checked for consistency::

        sage: Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)], lengths=[1, 1, 1, 1])
        Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)])

        sage: Polygon(vertices=[(0, 0), (1, 0), (0, 1)], lengths=[1, 1, 1])
        Traceback (most recent call last):
        ...
        ValueError: polygon does not have the prescribed lengths

    Currently, we cannot create a polygon from just lengths::

        sage: Polygon(lengths=[1, 1, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: one of vertices, edges, or angles must be set

    Polygons do not have to be convex::

        sage: p = Polygon(vertices=[(0, 0), (1, 1), (2, 0), (2, 4), (0, 4)])
        sage: p.describe_polygon()
        ('a', 'non-convex pentagon', 'non-convex pentagons')

    Polygons must be positively oriented::

        sage: Polygon(vertices=[(0, 0), (0, 1), (1, 0)])
        Traceback (most recent call last):
        ...
        ValueError: polygon has negative area; probably the vertices are not in counter-clockwise order

    Polygons must have at least three sides::

        sage: Polygon(vertices=[(0, 0), (1, 0)])
        Traceback (most recent call last):
        ...
        ValueError: polygon must have at least three sides

        sage: Polygon(vertices=[(0, 0), (1, 0), (2, 0)])
        Traceback (most recent call last):
        ...
        ValueError: polygon has zero area

    Currently, polygons must not self-intersect::

        sage: p = Polygon(vertices=[(0, 0), (2, 0), (0, 1), (1, -1), (2, 1)])
        Traceback (most recent call last):
        ...
        NotImplementedError: polygon self-intersects

    Currently, all angles must be less than 2π::

        sage: p = Polygon(angles=[14, 1, 1, 1, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: each angle must be in (0, 2π)

    """
    if "base_point" in kwds:
        base_point = kwds.pop("base_point")
        import warnings

        warnings.warn(
            "base_point has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use .translate() on the resulting polygon instead"
        )
        return Polygon(
            *args,
            vertices=vertices,
            edges=edges,
            angles=angles,
            lengths=lengths,
            base_ring=base_ring,
            category=category,
            **kwds,
        ).translate(base_point)

    if "ring" in kwds:
        import warnings

        warnings.warn(
            "ring has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use base_ring instead"
        )
        base_ring = kwds.pop("ring")

    if "field" in kwds:
        import warnings

        warnings.warn(
            "field has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use base_ring instead"
        )
        base_ring = kwds.pop("field")

    convex = None
    if "convex" in kwds:
        convex = kwds.pop("convex")
        import warnings

        if convex:
            warnings.warn(
                "convex has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; it has no effect other than checking the input for convexity so you may just drop it"
            )
        else:
            warnings.warn(
                "convex has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; it has no effect anymore, polygons are always allowed to be non-convex"
            )

    if args:
        import warnings

        warnings.warn(
            "calling Polygon() with positional arguments has been deprecated and will not be supported in a future version of sage-flatsurf; use edges= or vertices= explicitly instead"
        )

        edges = args

    if angles:
        if "length" in kwds:
            import warnings

            warnings.warn(
                "length has been deprecated as a keyword argument to Polygon() and will be removed in a future version of sage-flatsurf; use lengths instead"
            )

            lengths = [kwds.pop("length")] * (len(angles) - 2)

    if kwds:
        raise ValueError("keyword argument not supported by Polygon()")

    # Determine the number of sides of this polygon.
    if angles:
        n = len(angles)
    elif edges:
        n = len(edges)
    elif vertices:
        n = len(vertices)
    else:
        raise NotImplementedError("one of vertices, edges, or angles must be set")

    if n < 3:
        raise ValueError("polygon must have at least three sides")

    # Determine the base ring of the polygon
    if base_ring is None:
        base_ring = _Polygon_base_ring(vertices, edges, angles, lengths)

    if category is None:
        from flatsurf.geometry.categories import EuclideanPolygons

        # Currently, all polygons are assumed to be without self-intersection, i.e., simple.
        category = EuclideanPolygons(base_ring).Simple()
        if angles:
            category = category.WithAngles(angles)

        if n == 3:
            category = category.Convex()

    # We now rewrite the given data into vertices. Whenever there is
    # redundancy, we check that things are compatible. Note that much of the
    # complication of the below comes from the "angles" keyword. When angles
    # are given, some of the vertex coordinates can be deduced automatically.

    choice, vertices, edges, angles, lengths = _Polygon_normalize_arguments(
        category, n, vertices, edges, angles, lengths
    )

    assert vertices

    vertices = [vector(base_ring, vertex) for vertex in vertices]

    # Deduce missing vertices for prescribed angles
    if angles and len(vertices) != n:
        vertices = _Polygon_complete_vertices(n, vertices, angles, choice=choice)
        angles = None

    assert (
        len(vertices) == n
    ), f"expected to build {n}-gon from {n} vertices but found {vertices}"

    p = EuclideanPolygon(base_ring=base_ring, vertices=vertices, category=category)

    if check:
        _Polygon_check(p, vertices, edges, angles, lengths, convex)

    return p


def _Polygon_base_ring(vertices, edges, angles, lengths):
    r"""
    Return the base ring a polygon can be defined over.

    This is a helper function for :func:`Polygon`.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import _Polygon_base_ring
        sage: _Polygon_base_ring(vertices=[(0, 0), (1, 0), (0, 1)], edges=None, angles=None, lengths=None)
        Rational Field
        sage: _Polygon_base_ring(vertices=None, edges=[(1, 0), (-1, 1), (0, -1)], angles=None, lengths=None)
        Rational Field
        sage: _Polygon_base_ring(vertices=None, edges=None, angles=[1, 1, 1], lengths=None)
        Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?
        sage: _Polygon_base_ring(vertices=None, edges=None, angles=[1, 1, 1], lengths=[AA(2).sqrt(), 1])
        Algebraic Real Field

    """
    from sage.categories.pushout import pushout

    base_ring = QQ

    if angles:
        from flatsurf import EuclideanPolygonsWithAngles

        base_ring = pushout(base_ring, EuclideanPolygonsWithAngles(angles).base_ring())

    if vertices:
        base_ring = pushout(
            base_ring,
            Sequence([v[0] for v in vertices] + [v[1] for v in vertices]).universe(),
        )

    if edges:
        base_ring = pushout(
            base_ring,
            Sequence([e[0] for e in edges] + [e[1] for e in edges]).universe(),
        )

    if lengths:
        base_ring = pushout(base_ring, Sequence(lengths).universe())

        if angles and not edges:
            with_angles = (
                EuclideanPolygonsWithAngles(angles)
                ._without_axiom("Simple")
                ._without_axiom("Convex")
            )
            for slope, length in zip(with_angles.slopes(), lengths):
                scale = base_ring(length**2 / (slope[0] ** 2 + slope[1] ** 2))
                try:
                    is_square = scale.is_square()
                except NotImplementedError:
                    import warnings

                    warnings.warn(
                        "Due to https://github.com/flatsurf/exact-real/issues/173, we cannot compute the minimal base ring over which this polygon is defined. The polygon could possibly have been defined over a smaller ring."
                    )
                    is_square = False

                if not is_square:
                    # Note that this ring might not be minimal.
                    base_ring = pushout(base_ring, with_angles._cosines_ring())

    return base_ring


def _Polygon_normalize_arguments(category, n, vertices, edges, angles, lengths):
    r"""
    Return the normalized arguments defining a polygon. Additionally, a flag is
    returned that indicates whether we made a choice in normalizing these
    arguments.

    This is a helper function for :func:`Polygon`.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import _Polygon_normalize_arguments
        sage: from flatsurf.geometry.categories import EuclideanPolygons
        sage: category = EuclideanPolygons(AA)
        sage: _Polygon_normalize_arguments(category=category, n=3, vertices=[(0, 0), (1, 0), (0, 1)], edges=None, angles=None, lengths=None)
        (False, [(0, 0), (1, 0), (0, 1)], None, None, None)
        sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=[(1, 0), (-1, 1), (0, -1)], angles=None, lengths=None)
        (False, [(0, 0), (1, 0), (0, 1)], None, None, None)

        sage: category = category.WithAngles([1, 1, 1])
        sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=None, angles=[1, 1, 1], lengths=None)
        (True, [(0, 0), (1, 0), (1/2, 0.866025403784439?)], None, None, None)
        sage: _Polygon_normalize_arguments(category=category, n=3, vertices=None, edges=None, angles=[1, 1, 1], lengths=[AA(2).sqrt(), 1])
        (False,
         [(0, 0), (1.414213562373095?, 0), (0.9142135623730951?, 0.866025403784439?)],
         None,
         [1, 1, 1],
         None)

    """
    base_ring = category.base_ring()

    # Track whether we made a choice that possibly is the reason that we fail
    # to find a polygon with the given data.
    choice = False

    # Rewrite angles and lengths as angles and edges.
    if angles and lengths and not edges:
        edges = []
        for slope, length in zip(category.slopes(), lengths):
            scale = base_ring((length**2 / (slope[0] ** 2 + slope[1] ** 2)).sqrt())
            edges.append(scale * slope)

        if len(edges) == n:
            angles = 0

        lengths = None

    # Deduce edges if only angles are given
    if angles and not edges and not vertices:
        assert not lengths

        choice = True

        # We pick the edges such that they form a closed polygon with the
        # prescribed angles. However, there might be self-intersection which
        # currently leads to an error.
        edges = [
            length * slope
            for (length, slope) in zip(
                sum(r.vector() for r in category.lengths_polytope().rays()),
                category.slopes(),
            )
        ]

        angles = None

    # Rewrite edges as vertices.
    if edges and not vertices:
        vertices = [vector(base_ring, (0, 0))]
        for edge in edges:
            vertices.append(vertices[-1] + vector(base_ring, edge))

        if len(vertices) == n + 1:
            if vertices[-1]:
                raise ValueError("polygon not closed")
            vertices.pop()

        edges = None

    return choice, vertices, edges, angles, lengths


def _Polygon_complete_vertices(n, vertices, angles, choice):
    r"""
    Return vertices that define a polygon by completing the ``vertices`` to an
    ``n``-gon with ``angles``.

    This is a helper function for :func:`Polygon`.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import _Polygon_complete_vertices
        sage: _Polygon_complete_vertices(3, [vector((0, 0)), vector((1, 0))], [1, 1, 1], choice=False)
        [(0, 0), (1, 0), (1/2, 1/2*c)]

    """
    if len(vertices) == n - 1:
        # We do not use category.slopes() since the matrix formed by such
        # slopes might not be invertible (because exact-reals do not have a
        # fraction field implemented).
        slopes = EuclideanPolygonsWithAngles(angles).slopes()

        # We do not use solve_left() because the vertices might not live in
        # a ring that has a fraction field implemented (such as an
        # exact-real ring).
        s, t = (vertices[0] - vertices[n - 2]) * matrix(
            [slopes[-1], slopes[n - 2]]
        ).inverse()
        assert vertices[0] - s * slopes[-1] == vertices[n - 2] + t * slopes[n - 2]

        if s <= 0 or t <= 0:
            raise (NotImplementedError if choice else ValueError)(
                "cannot determine polygon with these angles from the given data"
            )

        vertices.append(vertices[0] - s * slopes[-1])

    if len(vertices) != n:
        from flatsurf.geometry.categories import Polygons

        raise NotImplementedError(
            f"cannot construct {' '.join(Polygons._describe_polygon(n)[:2])} from {n} angles and {len(vertices)} vertices"
        )

    return vertices


def _Polygon_check(p, vertices, edges, angles, lengths, convex):
    r"""
    Verify that ``p`` is a valid polygon and that it satisfies the constraints
    given.

    This is a helper function for :func:`Polygon`.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import _Polygon_check, Polygon
        sage: p = Polygon(angles=[1, 1, 1])
        sage: _Polygon_check(p, vertices=None, edges=None, angles=[1, 1, 1], lengths=None, convex=None)

    """
    # Check that the polygon satisfies the assumptions of EuclideanPolygon
    area = p.area()

    if area < 0:
        raise ValueError(
            "polygon has negative area; probably the vertices are not in counter-clockwise order"
        )

    if area == 0:
        raise ValueError("polygon has zero area")

    if any(edge == 0 for edge in p.edges()):
        raise ValueError("polygon has zero edge")

    for i in range(len(p.vertices())):
        from flatsurf.geometry.euclidean import is_anti_parallel

        if is_anti_parallel(p.edge(i), p.edge(i + 1)):
            raise ValueError("polygon has anti-parallel edges")

    from flatsurf.geometry.categories import EuclideanPolygons

    if not EuclideanPolygons.ParentMethods.is_simple(p):
        raise NotImplementedError("polygon self-intersects")

    # Check that any redundant data is compatible
    if edges:
        # Check compatibility of vertices and edges
        edges = [vector(p.base_ring(), edge) for edge in edges]
        if len(edges) != len(vertices):
            raise ValueError("vertices and edges must have the same length")

        for i in range(len(p.vertices())):
            if vertices[i - 1] + edges[i - 1] != vertices[i]:
                raise ValueError("vertices and edges are not compatible")

    if angles:
        # Check that the polygon has the prescribed angles
        from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
            EuclideanPolygonsWithAngles,
        )
        from flatsurf.geometry.categories.euclidean_polygons import (
            EuclideanPolygons,
        )

        # Use EuclideanPolygon's angle() so we do not use the precomputed angles set by the category.
        if EuclideanPolygonsWithAnglesCategory._normalize_angles(angles) != tuple(
            EuclideanPolygons.ParentMethods.angle(p, i)
            for i in range(len(p.vertices()))
        ):
            raise ValueError("polygon does not have the prescribed angles")

    if lengths:
        for edge, length in zip(p.edges(), lengths):
            if edge.norm() != length:
                raise ValueError("polygon does not have the prescribed lengths")

    if convex and not p.is_convex():
        raise ValueError("polygon is not convex")


def EuclideanPolygonsWithAngles(*angles):
    r"""
    Return the category of Euclidean polygons with prescribed ``angles``
    over a (minimal) number field.

    This method is a convenience to interact with that category. To create
    polygons with prescribed angles over such a field, one should just use
    ``Polygon()`` directly, see below.

    INPUT:

    - ``angles`` -- a sequence of integers or rationals describing the angles
      of the polygon (the number get normalized so that they sum to (n-2)π
      automatically.

    TESTS::

        sage: from flatsurf import EuclideanPolygonsWithAngles

    The polygons with inner angles `\pi/4`, `\pi/2`, `5\pi/4`::

        sage: P = EuclideanPolygonsWithAngles(1, 2, 5)
        sage: P
        Category of simple euclidean triangles with angles (1/16, 1/8, 5/16) over Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    Internally, polygons are given by their vertices' coordinates over some
    number field, in this case a quadratic field::

        sage: P.base_ring()
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    Polygons with these angles can be created by providing a single length,
    however this feature is deprecated::

        sage: P(1)
        doctest:warning
        ...
        UserWarning: calling EuclideanPolygonsWithAngles() has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon(angles=[...], lengths=[...]) instead.
        To make the resulting polygon non-normalized, i.e., the lengths are not actual edge lengths but the multiple of slope vectors, use Polygon(edges=[length * slope for (length, slope) in zip(lengths, EuclideanPolygonsWithAngles(angles).slopes())]).
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])

    Instead, one should use :func:`Polygon`::

        sage: from flatsurf import Polygon
        sage: Polygon(angles=[1, 2, 5], lengths=[1])
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])

    It is actually faster not to specify lengths since normalization can be
    costly (only relevant for polygons living in big number fields)::

        sage: Polygon(angles=[1, 2, 5])
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])

    Polygons can also be defined over other number field implementations::

        sage: from pyeantic import RealEmbeddedNumberField # optional: pyeantic  # random output due to matplotlib warnings with some combinations of setuptools and matplotlib
        sage: K = RealEmbeddedNumberField(P.base_ring()) # optional: pyeantic
        sage: p = P(K(1)) # optional: pyeantic  # random output due to deprecation warnings
        sage: p  # optional: pyeantic
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])
        sage: p.base_ring() # optional: pyeantic
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

    However, specific instances of such polygons might be defined over another ring::

        sage: P(1).base_ring()
        Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?

        sage: P(AA(1))
        Polygon(vertices=[(0, 0), (1, 0), (0.7071067811865475?, 0.2928932188134525?)])
        sage: _.base_ring()
        Algebraic Real Field

    Polygons can also be defined over a module containing transcendent parameters::

        sage: from pyexactreal import ExactReals # optional: pyexactreal  # random output due to deprecation warnings with some versions of pkg_resources
        sage: R = ExactReals(P.base_ring()) # optional: pyexactreal
        sage: P(R(1)) # optional: pyexactreal
        Polygon(vertices=[(0, 0), (1, 0), ((1/2*c0 ~ 0.70710678), (-1/2*c0+1 ~ 0.29289322))])
        sage: P(R(R.random_element([0.2, 0.3]))) # random output including some deprecation warnings, optional: pyexactreal
        Polygon(vertices=[(0, 0),])
                 (ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54), 0),
                 (((12*c0+17 ~ 33.970563)*ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54))/((17*c0+24 ~ 48.041631)),
                 ((5*c0+7 ~ 14.071068)*ℝ(0.287373=2588422249976937p-53 + ℝ(0.120809…)p-54))/((17*c0+24 ~ 48.041631)))
        sage: _.base_ring() # optional: pyexactreal
        Real Numbers as (Real Embedded Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?)-Module

    ::

        sage: L = P.lengths_polytope()    # polytope of admissible lengths for edges
        sage: L
        A 1-dimensional polyhedron in (Number Field in c0 with defining polynomial x^2 - 2 with c0 = 1.414213562373095?)^3 defined as the convex hull of 1 vertex and 1 ray
        sage: lengths = L.rays()[0].vector()
        sage: lengths
        (1, -1/2*c0 + 1, -1/2*c0 + 1)
        sage: p = P(*lengths)    # build one polygon with the given lengths
        sage: p
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c0, -1/2*c0 + 1)])
        sage: p.angles()
        (1/16, 1/8, 5/16)
        sage: P.angles(integral=False)
        (1/16, 1/8, 5/16)
        sage: P.angles(integral=True)
        (1, 2, 5)

        sage: P = EuclideanPolygonsWithAngles(1, 2, 1, 2, 2, 1)
        sage: L = P.lengths_polytope()
        sage: L
        A 4-dimensional polyhedron in (Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?)^6 defined as the convex hull of 1 vertex and 6 rays
        sage: rays = [r.vector() for r in L.rays()]
        sage: rays
        [(1, 0, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (0, 1, 0, 0, c^2 - 3, c^2 - 2),
         (1/3*c^4 - 2*c^2 + 3, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (-c^4 + 4*c^2, 0, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c, 0, -1/6*c^5 + 5/6*c^3 - 2/3*c),
         (0, 1/3*c^4 - 2*c^2 + 3, c^2 - 3, 0, 0, 1/3*c^4 - c^2),
         (0, -c^4 + 4*c^2, 0, c^2 - 3, 0, -c^4 + 5*c^2 - 3)]
        sage: lengths = 3*rays[0] + rays[2] + 2*rays[3] + rays[4]
        sage: p = P(*lengths)
        sage: p
        Polygon(vertices=[(0, 0),
                          (-5/3*c^4 + 6*c^2 + 6, 0),
                          (3*c^5 - 5/3*c^4 - 16*c^3 + 6*c^2 + 18*c + 6, c^4 - 6*c^2 + 9),
                          (2*c^5 - 2*c^4 - 10*c^3 + 15/2*c^2 + 9*c + 5, -1/2*c^5 + c^4 + 5/2*c^3 - 3*c^2 - 2*c),
                          (2*c^5 - 10*c^3 - 3/2*c^2 + 9*c + 9, -3/2*c^5 + c^4 + 15/2*c^3 - 3*c^2 - 6*c),
                          (2*c^5 - 10*c^3 - 3*c^2 + 9*c + 12, -3*c^5 + c^4 + 15*c^3 - 3*c^2 - 12*c)])

        sage: p.angles()
        (2/9, 4/9, 2/9, 4/9, 4/9, 2/9)

        sage: EuclideanPolygonsWithAngles(1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 2, 1)
        Category of simple euclidean pentadecagons with angles (13/46, 13/23, 13/46, 13/23, 13/46, 13/23, 13/46, 13/23, 13/23, 13/23, 13/23, 13/46, 13/46, 13/23, 13/46) over Number Field in c with defining polynomial ...

    A regular pentagon::

        sage: E = EuclideanPolygonsWithAngles(1, 1, 1, 1, 1)
        sage: E(1, 1, 1, 1, 1, normalized=True)
        doctest:warning
        ...
        UserWarning: calling EuclideanPolygonsWithAngles() has been deprecated and will be removed in a future version of sage-flatsurf; use Polygon(angles=[...], lengths=[...]) instead.
        Polygon(vertices=[(0, 0), (1, 0), (1/2*c^2 - 1/2, 1/2*c), (1/2, 1/2*c^3 - c), (-1/2*c^2 + 3/2, 1/2*c)])

    """
    if len(angles) == 1 and isinstance(angles[0], (tuple, list)):
        angles = angles[0]

    angles = EuclideanPolygonsWithAnglesCategory._normalize_angles(angles)

    from flatsurf.geometry.categories.euclidean_polygons_with_angles import (
        _base_ring,
    )

    base_ring = _base_ring(angles)

    return EuclideanPolygons(base_ring).WithAngles(angles).Simple()


def EquiangularPolygons(*angles, **kwds):
    r"""
    EXAMPLES::

        sage: from flatsurf import EquiangularPolygons
        sage: EquiangularPolygons(1, 1, 1)
        doctest:warning
        ...
        UserWarning: EquiangularPolygons() has been deprecated and will be removed in a future version of sage-flatsurf; use EuclideanPolygonsWithAngles() instead
        Category of simple euclidean equilateral triangles over Number Field in c with defining polynomial x^2 - 3 with c = 1.732050807568878?

    """
    import warnings

    warnings.warn(
        "EquiangularPolygons() has been deprecated and will be removed in a future version of sage-flatsurf; use EuclideanPolygonsWithAngles() instead"
    )

    if "number_field" in kwds:
        from warnings import warn

        warn(
            "The number_field parameter has been removed in this release of sage-flatsurf. "
            "To create an equiangular polygon over a number field, do not pass this parameter; to create an equiangular polygon over the algebraic numbers, do not pass this parameter but call the returned object with algebraic lengths."
        )
        kwds.pop("number_field")

    if kwds:
        raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

    return EuclideanPolygonsWithAngles(*angles)
