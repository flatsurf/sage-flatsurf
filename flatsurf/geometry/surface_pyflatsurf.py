r"""
Wrapper over pyflatsurf flat triangulations.

This module implements the protocol that defines a surface in sage-flatsurf
from the data of a flat triangulation in pyflatsurf.
"""

import cppyy
import gmpxxyy
from pyeantic import RealEmbeddedNumberField

import pyflatsurf
import pyflatsurf.vector

from .surface import Surface
from .pyflatsurf_conversion import sage_base_ring
from .polygon import ConvexPolygons

class Surface_pyflatsurf(Surface):
    r"""
    A surface built from an underlying pyflatsurf triangulation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf # optional: pyflatsurf
        sage: from flatsurf.geometry.surface_pyflatsurf import Surface_pyflatsurf # optional: pyflatsurf
        sage: T = to_pyflatsurf(translation_surfaces.veech_double_n_gon(5)) # optional: pyflatsurf
        sage: S = Surface_pyflatsurf(T) # optional: pyflatsurf
        sage: S.num_polygons() # optional: pyflatsurf
        6
        sage: S.num_edges() # optional: pyflatsurf
        9
        sage: S.base_ring() # optional: pyflatsurf
        Number Field in a with defining polynomial x^4 - 5*x^2 + 5 with a = 1.902113032590308?
        sage: S.opposite_edge(0, 0) # optional: pyflatsurf
        (2, 2)
        sage: S.opposite_edge(2, 2) # optional: pyflatsurf
        (0, 0)
        sage: TestSuite(S).run() # optional: pyflatsurf
    """
    def __init__(self, flat_triangulation):
        self._flat_triangulation = flat_triangulation
        self._base_ring, self._to_base_ring = sage_base_ring(flat_triangulation)
        self._base_label = 0
        self._mutable = False
        self._finite = True

        self._polygons = [] # list of triangles (as sage-flatsurf ConvexPolygon)
        self._half_edge_to_label = {} # pyflatsurf half edge -> (polygon_label, edge_number)

        P = ConvexPolygons(self._base_ring)
        V = P.module()
        for face in flat_triangulation.faces():
            a, b, c = map(cppyy.gbl.flatsurf.HalfEdge, face)
            vectors = [flat_triangulation.fromHalfEdge(he) for he in face]
            vectors = [V([self._to_base_ring(v.x()), self._to_base_ring(v.y())]) for v in vectors]
            triangle = P(vectors)

            face_id = len(self._polygons)
            self._polygons.append(triangle)
            assert a not in self._half_edge_to_label
            self._half_edge_to_label[a] = (face_id, 0)
            assert b not in self._half_edge_to_label
            self._half_edge_to_label[b] = (face_id, 1)
            assert c not in self._half_edge_to_label
            self._half_edge_to_label[c] = (face_id, 2)

    def __reduce__(self):
        return Surface_pyflatsurf, (self._flat_triangulation,)

    def is_triangulated(self):
        return True

    def polygon(self, label):
        r"""
        Return the polygon with the given ``label``.
        """
        return self._polygons[label]

    def opposite_edge(self, l, e):
        h = self._flat_triangulation.faces()[int(l)][int(e)]
        return self._half_edge_to_label[-h]

    def num_polygons(self):
        return len(self._polygons)

    def num_edges(self):
        return 3 * len(self._polygons) // 2

    def area(self):
        raise NotImplementedError

    def label_iterator(self):
        return iter(range(self.num_polygons()))

    def edge_iterator(self):
        return ((l, e) for l in self.label_iterator() for e in range(3))
