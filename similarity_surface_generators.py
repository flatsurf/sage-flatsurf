from sage.rings.integer_ring import ZZ

ZZ_2 = ZZ(2)

class SimilaritySurfaceGenerators:
    r"""
    Examples of similarity surfaces.
    """
    @staticmethod
    def example():
        r"""
        Construct a SimilaritySurface from a pair of triangles.
        """
        from similarity_surface import SimilaritySurface_polygons_and_gluings
        from polygon import PolygonCreator
        pc=PolygonCreator()
        pc.add_vertex((0,0))
        pc.add_vertex((2,-2))
        pc.add_vertex((2,0))
        p0=pc.get_polygon()
        pc=PolygonCreator()
        pc.add_vertex((0,0))
        pc.add_vertex((2,0))
        pc.add_vertex((1,3))
        p1=pc.get_polygon()
        ps=(p0,p1)
        glue={ (0,2):(1,0), (0,0):(1,1), (0,1):(1,2), (1,0):(0,2), (1,1):(0,0), (1,2):(0,1) }
        return SimilaritySurface_polygons_and_gluings(ps,glue)


class TranslationSurfaceGenerators:
    r"""
    Common and less common translation surfaces.
    """
    @staticmethod
    def regular_octagon():
        r"""
        Return the translation surface built from the regular octagon by
        identifying opposite sides.

        EXAMPLES::

            sage: T = translation_surfaces.regular_octagon()
            sage: T
            Translation surface built from the regular octagon
            sage: T.stratum()
            H_2(2)
        """
        from polygon import regular_octagon
        from similarity_surface import TranslationSurface_polygons_and_gluings
        polygons = [regular_octagon()]
        identifications = {}
        identifications.update(dict(((0,i),(0,i+4)) for i in xrange(4)))
        return TranslationSurface_polygons_and_gluings(polygons=polygons, identifications=identifications)

    @staticmethod
    def octagon_and_squares():
        from polygon import square, regular_octagon
        from sage.matrix.matrix_space import MatrixSpace
        from similarity_surface import TranslationSurface_polygons_and_gluings

        o = regular_octagon()
        K = o.parent().field()
        sqrt2 = K.gen()
        rot = MatrixSpace(K,2)([[sqrt2/ZZ_2,-sqrt2/ZZ_2],[sqrt2/ZZ_2,sqrt2/ZZ_2]])
        polygons = [regular_octagon(), ZZ_2*square(K), ZZ_2*rot*square(K)]
        identifications = {
            (0,0): (1,3),
            (0,1): (2,3),
            (0,2): (1,0),
            (0,3): (2,0),
            (0,4): (1,1),
            (0,5): (2,1),
            (0,6): (1,2),
            (0,7): (2,2),
            }
        return TranslationSurface_polygons_and_gluings(polygons=polygons, identifications=identifications)

    @staticmethod
    def origami(r,u,rr=None,uu=None,domain=None):
        r"""
        Return the origami defined by the permutations ``r`` and ``u``.

        EXAMPLES::

            sage: S = SymmetricGroup(3)
            sage: r = S('(1,2)')
            sage: u = S('(1,3)')
            sage: o = translation_surfaces.origami(r,u)
            sage: o
            Origami defined by r=(1,2) and u=(1,3)
            sage: o.stratum()
            H_2(2)
        """
        from similarity_surface import Origami
        return Origami(r,u,rr,uu,domain)

    @staticmethod
    def infinite_origami_example():
        from similarity_surface import Origami

        return Origami(
            lambda x: x+1,
            lambda x: x-1,
            lambda x: x-1,
            lambda x: x+1,
            ZZ)
