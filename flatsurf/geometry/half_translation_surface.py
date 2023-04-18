# ****************************************************************************
#       Copyright (C) 2013-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2013-2019 W. Patrick Hooper <wphooper@gmail.com>
#                          2023 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .polygon import wedge_product
from .half_dilation_surface import HalfDilationSurface
from .rational_cone_surface import RationalConeSurface

from sage.rings.all import QQ, AA
from sage.matrix.constructor import matrix


class HalfTranslationSurface(HalfDilationSurface, RationalConeSurface):
    r"""
    A half translation surface has gluings between polygons whose monodromy is +I or -I.
    """

    def __init__(self, surface, category=None):
        from flatsurf.geometry.categories.half_translation_surfaces import HalfTranslationSurfaces
        super().__init__(surface, category or surface.category() & HalfTranslationSurfaces().Oriented())

    def holonomy_field(self):
        r"""
        Return the relative holonomy field of this translation or half-translation surface.

        EXAMPLES::

            sage: from flatsurf import *

            sage: S = translation_surfaces.veech_2n_gon(5)
            sage: S.holonomy_field()
            Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = 1.618033988749895?
            sage: S.base_ring()
            Number Field in a with defining polynomial y^4 - 5*y^2 + 5 with a = 1.175570504584947?

            sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
            sage: T.holonomy_field()
            Rational Field

            sage: T = polygons.triangle(1,6,11)
            sage: S = similarity_surfaces.billiard(T)
            sage: S = S.minimal_cover("translation")
            sage: S.base_ring()
            Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?
            sage: S.holonomy_field()
            Number Field in c0 with defining polynomial x^3 - 3*x - 1 with c0 = 1.879385241571817?
        """
        return self.normalized_coordinates()[0].base_ring()

    def normalized_coordinates(self):
        r"""
        Return a pair ``(new_surface, matrix)`` where ``new_surface`` is defined over the
        holonomy field and ``matrix`` is the transition matrix that maps this surface to
        ``new_surface``.

        EXAMPLES::

            sage: from flatsurf import *

            sage: S = translation_surfaces.veech_2n_gon(5)
            sage: U, mat = S.normalized_coordinates()
            sage: U.base_ring()
            Number Field in a0 with defining polynomial x^2 - x - 1 with a0 = 1.618033988749895?
            sage: mat
            [             0 -2/5*a^3 + 2*a]
            [            -1 -3/5*a^3 + 2*a]

            sage: T = translation_surfaces.torus((1, AA(2).sqrt()), (AA(3).sqrt(), 3))
            sage: U, mat = T.normalized_coordinates()
            sage: U.base_ring()
            Rational Field
            sage: U.holonomy_field()
            Rational Field
            sage: mat
            [-2.568914100752347?  1.816496580927726?]
            [-5.449489742783178?  3.146264369941973?]
            sage: TestSuite(U).run()

            sage: T = polygons.triangle(1,6,11)
            sage: S = similarity_surfaces.billiard(T)
            sage: S = S.minimal_cover("translation")
            sage: U, _ = S.normalized_coordinates()
            sage: U.base_ring()
            Number Field in c0 with defining polynomial x^3 - 3*x - 1 with c0 = 1.879385241571817?
            sage: U.holonomy_field() == U.base_ring()
            True
            sage: S.base_ring()
            Number Field in c with defining polynomial x^6 - 6*x^4 + 9*x^2 - 3 with c = 1.969615506024417?
            sage: TestSuite(U).run()

            sage: from flatsurf import EquiangularPolygons
            sage: E = EquiangularPolygons(1, 3, 1, 1)
            sage: r1, r2 = [r.vector() for r in E.lengths_polytope().rays()]
            sage: p = E(r1 + r2)
            sage: B = similarity_surfaces.billiard(p)
            sage: B.minimal_cover("translation")
            TranslationSurface built from 6 polygons
            sage: S = B.minimal_cover("translation")
            sage: S, _ = S.normalized_coordinates()
            sage: S
            TranslationSurface built from 6 polygons
        """
        if not self.is_finite():
            raise ValueError("the surface must be finite")
        if self.base_ring() is QQ:
            return (self, matrix(QQ, 2, 2, 1))

        lab = next(self.label_iterator())
        p = self.polygon(lab)
        u = p.edge(1)
        v = -p.edge(0)
        i = 1
        while wedge_product(u, v) == 0:
            i += 1
            u = p.edge(i)
            v = -p.edge(i - 1)
        M = matrix(2, [u, v]).transpose().inverse()
        assert M.det() > 0
        hols = []
        for lab in self.label_iterator():
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                w = M * p.edge(e)
                hols.append(w[0])
                hols.append(w[1])
        if self.base_ring() is AA:
            from .subfield import number_field_elements_from_algebraics

            K, new_hols = number_field_elements_from_algebraics(hols)
        else:
            from .subfield import subfield_from_elements

            K, new_hols, _ = subfield_from_elements(self.base_ring(), hols)

        from .polygon import ConvexPolygons
        from .surface import Surface_list

        S = Surface_list(K)
        C = ConvexPolygons(K)
        relabelling = {}
        k = 0
        for lab in self.label_iterator():
            m = self.polygon(lab).num_edges()
            relabelling[lab] = S.add_polygon(
                C(
                    edges=[
                        (new_hols[k + 2 * i], new_hols[k + 2 * i + 1]) for i in range(m)
                    ]
                )
            )
            k += 2 * m

        for (p1, e1), (p2, e2) in self.edge_iterator(gluings=True):
            S.set_edge_pairing(relabelling[p1], e1, relabelling[p2], e2)

        return (type(self)(S, category=self.category()), M)
