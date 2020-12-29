from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from .similarity_surface import SimilaritySurface
from .matrix_2x2 import is_cosine_sine_of_rational

class RationalSimilaritySurface(SimilaritySurface):
    r"""
    A similarity surface such that the monodromy around any loop is similarity
    whose rotational part has finite order.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s = Surface_list(AA)
        sage: CP = ConvexPolygons(AA)
        sage: s.add_polygon(CP(vertices=[(0,0),(2,0),(1,1)]))
        0
        sage: s.add_polygon(CP(vertices=[(0,0),(1,0),(1,sqrt(3))]))
        1
        sage: for i in range(3):
        ....:     s.change_edge_gluing(0,i,1,i)
        sage: s.change_base_label(0)
        sage: s.set_immutable()
        sage: TestSuite(s).run()
        sage: from flatsurf.geometry.rational_similarity_surface import RationalSimilaritySurface
        sage: ss = RationalSimilaritySurface(s)
        sage: TestSuite(ss).run()

    Example of a similarity surface which is not rational::

        sage: from flatsurf import *
        sage: s = Surface_list(QQ)
        sage: CP = ConvexPolygons(QQ)
        sage: s.add_polygon(CP(vertices=[(0,0),(2,0),(1,1)]))
        0
        sage: s.add_polygon(CP(vertices=[(0,0),(1,0),(1,2)]))
        1
        sage: for i in range(3):
        ....:     s.change_edge_gluing(0,i,1,i)
        sage: s.change_base_label(0)
        sage: s.set_immutable()
        sage: TestSuite(s).run()
        sage: from flatsurf.geometry.rational_similarity_surface import RationalSimilaritySurface
        sage: ss = RationalSimilaritySurface(s)
        sage: TestSuite(ss).run()
        ...
        The following tests failed: _test_edge_matrix
    """
    def _test_edge_matrix(self, **options):
        r"""
        Check the compatibility condition
        """
        tester = self._tester(**options)

        from .similarity_surface import SimilaritySurface
        from sage.rings.qqbar import AA

        if self.is_finite():
            it = self.label_iterator()
        else:
            from itertools import islice
            it = islice(self.label_iterator(), 30)

        for lab in it:
            p = self.polygon(lab)
            for e in range(p.num_edges()):
                # Warning: check the matrices computed from the edges,
                # rather the ones overriden by TranslationSurface.
                m = SimilaritySurface.edge_matrix(self,lab,e)
                a = AA(m[0,0])
                b = AA(m[1,0])
                q = (a**2 + b**2).sqrt()
                a /= q
                b /= q
                tester.assertTrue(is_cosine_sine_of_rational(a, b))
