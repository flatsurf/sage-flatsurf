from sage.matrix.constructor import matrix

from flatsurf.geometry.surface import OrientedSimilaritySurface


def _is_finite(surface):
    r"""
    Return whether ``surface`` is a finite rational cone surface.
    """
    if not surface.is_finite():
        return False

    from flatsurf.geometry.categories import ConeSurfaces

    if surface in ConeSurfaces().Rational():
        return True

    surface = surface.reposition_polygons(relabel=True)

    for label in surface.label_iterator():
        polygon = surface.polygon(label)

        for e in range(polygon.num_edges()):
            m = surface.edge_matrix(label, e)

            from flatsurf.geometry.matrix_2x2 import is_cosine_sine_of_rational

            if not is_cosine_sine_of_rational(m[0][0], m[0][1]):
                return False

    return True


class MinimalTranslationCover(OrientedSimilaritySurface):
    r"""
    We label copy by cartesian product (polygon from bot, matrix).

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.minimal_cover import MinimalTranslationCover
        sage: s = Surface_list(QQ)
        sage: P = ConvexPolygons(QQ)
        sage: s.add_polygon(P(vertices=[(0,0),(5,0),(0,5)]))
        0
        sage: s.add_polygon(P(vertices=[(0,0),(3,4),(-4,3)]))
        1
        sage: s.change_polygon_gluings(0,[(1,2),(1,1),(1,0)])
        sage: s.set_immutable()
        sage: ss=MinimalTranslationCover(s)
        sage: ss.is_finite()
        True
        sage: len(ss.polygons())
        8
        sage: TestSuite(ss).run()

    The following is to test that unfolding is reasonably fast on the instances reported
    in https://github.com/flatsurf/sage-flatsurf/issues/47::

        sage: T = polygons.triangle(2, 13, 26)
        sage: S = similarity_surfaces.billiard(T, rational=True)
        sage: S = S.minimal_cover("translation")
        sage: S
        MinimalTranslationCover(Surface built from 2 polygons)

    TESTS::

        sage: from flatsurf.geometry.categories import TranslationSurfaces
        sage: S in TranslationSurfaces()
        True

    """

    def __init__(self, similarity_surface, category=None):
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                self._ss = MutableOrientedSimilaritySurface.from_surface(similarity_surface)
            else:
                raise ValueError(
                    "Can not construct MinimalTranslationCover of a surface that is mutable and infinite."
                )
        else:
            self._ss = similarity_surface

        if category is None:
            from flatsurf.geometry.categories import TranslationSurfaces
            category = TranslationSurfaces()

            if _is_finite(self._ss):
                category = category.FiniteType()
            else:
                category = category.InfiniteType()

        OrientedSimilaritySurface.__init__(self, self._ss.base_ring(), category=category)

    def base_label(self):
        self._F = self._ss.base_ring()
        return (self._ss.base_label(), self._F.one(), self._F.zero())

    def is_mutable(self):
        return False

    def polygon(self, lab):
        if not isinstance(lab, tuple) or len(lab) != 3:
            raise ValueError("invalid label {!r}".format(lab))
        return matrix([[lab[1], -lab[2]], [lab[2], lab[1]]]) * self._ss.polygon(lab[0])

    def opposite_edge(self, p, e):
        pp, a, b = p  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, e)
        m = self._ss.edge_matrix(p2, e2)
        aa = a * m[0][0] - b * m[1][0]
        bb = b * m[0][0] + a * m[1][0]
        return ((p2, aa, bb), e2)

    def _repr_(self):
        return f"MinimalTranslationCover({repr(self._ss)})"

    def __hash__(self):
        return super().__hash__()

    def _cache_key(self):
        return (MinimalTranslationCover, self._ss, self.category())

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import polygons, similarity_surfaces
            sage: T = polygons.triangle(2, 13, 26)
            sage: S = similarity_surfaces.billiard(T, rational=True)
            sage: S = S.minimal_cover("translation")

            sage: S == S
            True

        ::

            sage: TT = polygons.triangle(2, 15, 26)
            sage: SS = similarity_surfaces.billiard(TT, rational=True)
            sage: SS = SS.minimal_cover("translation")

            sage: S == SS
            False

        """
        if not isinstance(other, MinimalTranslationCover):
            return False

        if self.category() != other.category():
            return False

        if self._ss == other._ss:
            return True

        return self._eq_oriented_similarity_surfaces(other)


class MinimalHalfTranslationCover(OrientedSimilaritySurface):
    r"""
    We label copy by cartesian product (polygon from bot, matrix).

    EXAMPLES::

        sage: from flatsurf import *
        sage: from flatsurf.geometry.minimal_cover import MinimalHalfTranslationCover
        sage: s = Surface_list(QQ)
        sage: P = ConvexPolygons(QQ)
        sage: s.add_polygon(P(vertices=[(0,0),(5,0),(0,5)]))
        0
        sage: s.add_polygon(P(vertices=[(0,0),(3,4),(-4,3)]))
        1
        sage: s.change_polygon_gluings(0,[(1,2),(1,1),(1,0)])
        sage: s.set_immutable()
        sage: ss=MinimalHalfTranslationCover(s)
        sage: ss.is_finite()
        True
        sage: len(ss.polygons())
        4
        sage: TestSuite(ss).run()

    The following is to test that unfolding is reasonably fast on the instances reported
    in https://github.com/flatsurf/sage-flatsurf/issues/47::

        sage: T = polygons.triangle(2, 13, 26)
        sage: S = similarity_surfaces.billiard(T, rational=True)
        sage: S = S.minimal_cover("half-translation")
        sage: S
        MinimalHalfTranslationCover(Surface built from 2 polygons)

    TESTS::

        sage: from flatsurf.geometry.categories import DilationSurfaces
        sage: S in DilationSurfaces()
        True

    """

    def __init__(self, similarity_surface, category=None):
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                self._ss = MutableOrientedSimilaritySurface.from_surface(similarity_surface)
            else:
                raise ValueError(
                    "Can not construct MinimalTranslationCover of a surface that is mutable and infinite."
                )
        else:
            self._ss = similarity_surface

        if category is None:
            from flatsurf.geometry.categories import HalfTranslationSurfaces
            category = HalfTranslationSurfaces()

            if _is_finite(self._ss):
                category = category.FiniteType()
            else:
                category = category.InfiniteType()

        OrientedSimilaritySurface.__init__(self, self._ss.base_ring(), category=category)

    def base_label(self):
        self._F = self._ss.base_ring()
        return (self._ss.base_label(), self._F.one(), self._F.zero())

    def is_mutable(self):
        return False

    def _repr_(self):
        return f"MinimalHalfTranslationCover({repr(self._ss)})"

    def polygon(self, lab):
        if not isinstance(lab, tuple) or len(lab) != 3:
            raise ValueError("invalid label {!r}".format(lab))
        return matrix([[lab[1], -lab[2]], [lab[2], lab[1]]]) * self._ss.polygon(lab[0])

    def opposite_edge(self, p, e):
        pp, a, b = p  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, e)
        m = self._ss.edge_matrix(pp, e)
        aa = a * m[0][0] + b * m[1][0]
        bb = b * m[0][0] - a * m[1][0]
        if aa > 0 or (aa == 0 and bb > 0):
            return ((p2, aa, bb), e2)
        else:
            return ((p2, -aa, -bb), e2)

    def __hash__(self):
        return super().__hash__()

    def _cache_key(self):
        return (MinimalHalfTranslationCover, self._ss, self.category())

    def __eq__(self, other):
        if not isinstance(other, MinimalHalfTranslationCover):
            return False

        if self.category() != other.category():
            return False

        if self._ss == other._ss:
            return True

        return self._eq_oriented_similarity_surfaces(other)


class MinimalPlanarCover(OrientedSimilaritySurface):
    r"""
    The minimal planar cover of a surface S is the smallest cover C so that the
    developing map from the universal cover U to the plane induces a well
    defined map from C to the plane. This is a translation surface.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s = translation_surfaces.square_torus()
        sage: from flatsurf.geometry.minimal_cover import MinimalPlanarCover
        sage: pc = MinimalPlanarCover(s)
        sage: pc.is_finite()
        False
        sage: sing = pc.singularity(pc.base_label(),0,limit=4)
        doctest:warning
        ...
        UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
        sage: len(sing.vertex_set())
        doctest:warning
        ...
        UserWarning: vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
        4
        sage: TestSuite(s).run()
    """

    def __init__(self, similarity_surface, base_label=None, category=None):
        if similarity_surface.is_mutable():
            if similarity_surface.is_finite():
                from flatsurf.geometry.surface import MutableOrientedSimilaritySurface
                self._ss = MutableOrientedSimilaritySurface.from_surface(similarity_surface)
            else:
                raise ValueError(
                    "Can not construct MinimalPlanarCover of a surface that is mutable and infinite."
                )
        else:
            self._ss = similarity_surface

        if base_label is None:
            base_label = self._ss.base_label()

        # The similarity group containing edge identifications.
        self._sg = self._ss.edge_transformation(self._ss.base_label(), 0).parent()
        self._base_label = (self._ss.base_label(), self._sg.one())

        if category is None:
            from flatsurf.geometry.categories import TranslationSurfaces
            category = TranslationSurfaces().InfiniteType()

        OrientedSimilaritySurface.__init__(self, self._ss.base_ring(), category=category)

    def _repr_(self):
        return f"MinimalPlanarCover({repr(self._ss)})"

    def base_label(self):
        return self._base_label

    def is_mutable(self):
        return False

    def polygon(self, lab):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: C = translation_surfaces.chamanara(1/2)
            sage: C.polygon('a')
            Traceback (most recent call last):
            ...
            ValueError: invalid label 'a'
        """
        if not isinstance(lab, tuple) or len(lab) != 2:
            raise ValueError("invalid label {!r}".format(lab))
        return lab[1](self._ss.polygon(lab[0]))

    def opposite_edge(self, p, e):
        pp, m = p  # this is the polygon m * ss.polygon(p)
        p2, e2 = self._ss.opposite_edge(pp, e)
        me = self._ss.edge_transformation(pp, e)
        mm = m * ~me
        return ((p2, mm), e2)

    def __hash__(self):
        return super().__hash__()

    def _cache_key(self):
        return (MinimalPlanarCover, self._ss, self._base_label, self.category())

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.

        Note that this is not implemented in most non-trivial cases.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.square_torus()
            sage: from flatsurf.geometry.minimal_cover import MinimalPlanarCover
            sage: pc = MinimalPlanarCover(s)
            sage: pc == pc
            True

        """
        if not isinstance(other, MinimalPlanarCover):
            return False

        if self._base_label != other._base_label:
            return False
        if self.category() != other.category():
            return False
        if self._ss == other._ss:
            return True

        return self._eq_oriented_similarity_surfaces(other)
