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
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.all import FreeModule
from sage.misc.cachefunc import cached_method

from flatsurf.geometry.categories.polygons import Polygons


class RealProjectivePolygons(Category_over_base_ring):
    def super_categories(self):
        return [Polygons(self.base_ring())]

    class ParentMethods:
        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.vector_space()
                Vector space of dimension 2 over Rational Field

            """
            return self.category().vector_space()

        def module(self):
            r"""
            Return the free module of rank 2 in which this polygon embeds.

            EXAMPLES::

                sage: from flatsurf import polygons
                sage: S = polygons.square()
                sage: S.module()
                Vector space of dimension 2 over Rational Field

            """
            return self.category().module()

        def field(self):
            import warnings
            warnings.warn("field() has been deprecated and will be removed from a future version of sage-flatsurf; use base_ring() instead")

            return self.base_ring()

    class SubcategoryMethods:
        @cached_method
        def module(self):
            r"""
            Return the free module of rank 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.module()
                Vector space of dimension 2 over Rational Field

            """
            return FreeModule(self.base_ring(), 2)

        @cached_method
        def vector_space(self):
            r"""
            Return the vector space of dimension 2 in which these polygons embed.

            EXAMPLES::

                sage: from flatsurf import Polygons
                sage: C = Polygons(QQ)
                sage: C.vector_space()
                Vector space of dimension 2 over Rational Field

            """
            from sage.all import VectorSpace
            return VectorSpace(self.base_ring().fraction_field(), 2)

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf import Polygons, ConvexPolygons

            sage: C = Polygons(QQ)
            sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
            sage: p
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: C(p) is p
            True
            sage: C((1,0), (0,1), (-1, 1))
            Traceback (most recent call last):
            ...
            ValueError: the polygon does not close up

            sage: D = ConvexPolygons(QQbar)
            sage: D(p)
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(vertices=p.vertices())
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
            sage: D(edges=p.edges())
            polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
        """
        # TODO: We should deprecate this. It will conflict with more complex categories eventually.
        check = kwds.pop("check", True)

        from flatsurf.geometry.polygon import Polygon
        if len(args) == 1 and isinstance(args[0], Polygon):
            if args[0].category() is self:
                return args[0]
            vertices = map(self.vector_space(), args[0].vertices())
            args = ()

        else:
            vertices = kwds.pop("vertices", None)
            edges = kwds.pop("edges", None)
            base_point = kwds.pop("base_point", (0, 0))

            if (vertices is None) and (edges is None):
                if len(args) == 1:
                    edges = args[0]
                elif args:
                    edges = args
                else:
                    raise ValueError(
                        "exactly one of 'vertices' or 'edges' must be provided"
                    )
            if kwds:
                raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

            if edges is not None:
                v = self.vector_space()(base_point)
                vertices = []
                for e in map(self.vector_space(), edges):
                    vertices.append(v)
                    v += e
                if v != vertices[0]:
                    raise ValueError("the polygon does not close up")

        from flatsurf.geometry.polygon import EuclideanPolygon
        return EuclideanPolygon(ring=self.base(), vertices=vertices, category=self, check=check)

    class Convex(CategoryWithAxiom_over_base_ring):
        def __call__(self, *args, **kwds):
            r"""
            TESTS::

                sage: from flatsurf import ConvexPolygons

                sage: C = ConvexPolygons(QQ)
                sage: p = C(vertices=[(0,0),(1,0),(2,0),(1,1)])
                sage: p
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: C(p) is p
                True
                sage: C((1,0), (0,1), (-1, 1))
                Traceback (most recent call last):
                ...
                ValueError: the polygon does not close up

                sage: D = ConvexPolygons(QQbar)
                sage: D(p)
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: D(vertices=p.vertices())
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])
                sage: D(edges=p.edges())
                polygon(vertices=[(0, 0), (1, 0), (2, 0), (1, 1)])

            """
            # TODO: We should deprecate this. It will conflict with more complex categories eventually.
            check = kwds.pop("check", True)

            from flatsurf.geometry.polygon import Polygon
            if len(args) == 1 and isinstance(args[0], Polygon):
                if args[0].category() is self:
                    return args[0]

                vertices = map(self.vector_space(), args[0].vertices())
                args = ()

            else:
                vertices = kwds.pop("vertices", None)
                edges = kwds.pop("edges", None)
                base_point = kwds.pop("base_point", (0, 0))

                if (vertices is None) and (edges is None):
                    if len(args) == 1:
                        edges = args[0]
                    elif args:
                        edges = args
                    else:
                        raise ValueError(
                            "exactly one of 'vertices' or 'edges' must be provided"
                        )
                if kwds:
                    raise ValueError("invalid keyword {!r}".format(next(iter(kwds))))

                if edges is not None:
                    v = self.module()(base_point)
                    vertices = []
                    for e in map(self.module(), edges):
                        vertices.append(v)
                        v += e
                    if v != vertices[0]:
                        raise ValueError("the polygon does not close up")

            from flatsurf.geometry.polygon import EuclideanPolygon
            return EuclideanPolygon(ring=self.base(), vertices=vertices, category=self, check=check)
