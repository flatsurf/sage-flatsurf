from flatsurf.geometry.surface import OrientedSimilaritySurface


class Surface_pyflatsurf(OrientedSimilaritySurface):
    r"""
    EXAMPLES::

        sage: from flatsurf import Polygon, MutableOrientedSimilaritySurface

        sage: S = MutableOrientedSimilaritySurface(QQ)
        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 0), (1, 1)]), label=0)
        0
        sage: S.add_polygon(Polygon(vertices=[(0, 0), (1, 1), (0, 1)]), label=1)
        1

        sage: S.glue((0, 0), (1, 1))
        sage: S.glue((0, 1), (1, 2))
        sage: S.glue((0, 2), (1, 0))

        sage: S.set_immutable()

        sage: T = S.pyflatsurf().codomain()

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

        sage: isinstance(T, Surface_pyflatsurf)
        True
        sage: TestSuite(T).run()

    """

    def __init__(self, flat_triangulation):
        # TODO: Check that _flat_triangulation is never used by _flat_triangulation instead.
        self._flat_triangulation = flat_triangulation

        # TODO: We have to be smarter about the ring bridge here.
        from flatsurf.geometry.pyflatsurf_conversion import RingConversion

        base_ring = RingConversion.from_pyflatsurf_from_flat_triangulation(flat_triangulation).domain()

        from flatsurf.geometry.categories import TranslationSurfaces
        # TODO: This is assuming that the surface is connected. Currently that's the case for all surfaces in libflatsurf?
        category = TranslationSurfaces().FiniteType().Connected()
        if flat_triangulation.hasBoundary():
            category = category.WithBoundary()
        else:
            category = category.WithoutBoundary()
        super().__init__(base=base_ring, category=category)

    def __eq__(self, other):
        if not isinstance(other, Surface_pyflatsurf):
            return False

        return self._flat_triangulation == other._flat_triangulation

    def __hash__(self):
        # TODO: This is not working. Flat triangulations that are equal do not produce the same hashes.
        return hash(self._flat_triangulation)

    def flat_triangulation(self):
        return self._flat_triangulation

    def is_mutable(self):
        return False

    def roots(self):
        # TODO: This is assuming that the surface is connected. Currently that's the case for all surfaces in libflatsurf?
        return [self._normalize_label(self._flat_triangulation.face(1))]

    def pyflatsurf(self):
        from flatsurf.geometry.deformation import IdentityDeformation
        return IdentityDeformation(self)

    @classmethod
    def _from_flatsurf(cls, surface):
        r"""
        Return an isomorphism to a :class:`Surface_pyflatsurf` built from
        ``surface``, i.e., represent ``surface`` in pyflatsurf wrapped as a
        :class:`Surface` for sage-flatsurf.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: Surface_pyflatsurf._from_flatsurf(S)
            Generic morphism:
              From: Triangulation of Translation Surface in H_1(0) built from a square
              To:   FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        if isinstance(surface, Surface_pyflatsurf):
            return surface

        if not surface.is_triangulated():
            triangulation = surface.triangulate()
            to_pyflatsurf = cls._from_flatsurf(triangulation.codomain())
            return to_pyflatsurf * triangulation

        from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion

        to_pyflatsurf = FlatTriangulationConversion.to_pyflatsurf(surface)

        surface_pyflatsurf = Surface_pyflatsurf(to_pyflatsurf.codomain())

        from flatsurf.geometry.pyflatsurf.morphism import Morphism_to_pyflatsurf

        return Morphism_to_pyflatsurf._create_morphism(surface, surface_pyflatsurf, to_pyflatsurf)

    def __repr__(self):
        r"""
        Return a printable representation of this surface, namely, print this
        surface as pyflatsurf would.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: S.pyflatsurf().codomain()
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 0), 2: (0, 1), 3: (-1, -1)}

        """
        return repr(self._flat_triangulation)

    def apply_matrix(self, m):
        from sage.all import matrix
        m = matrix(m, ring=self.base_ring())

        from flatsurf.geometry.pyflatsurf_conversion import RingConversion

        to_pyflatsurf = RingConversion.from_pyflatsurf_from_flat_triangulation(self._flat_triangulation)

        m = [to_pyflatsurf(x) for x in m.list()]

        deformation = self._flat_triangulation.applyMatrix(*m)
        codomain = Surface_pyflatsurf(deformation.codomain())

        from flatsurf.geometry.pyflatsurf.deformation import Deformation_pyflatsurf
        return Deformation_pyflatsurf(self, codomain, deformation)

    @classmethod
    def _normalize_label(cls, label):
        label = tuple(edge.id() if hasattr(edge, "id") else int(edge) for edge in label)

        shift = label.index(min(*label))

        label = label[shift:] + label[:shift]
        return label

    def polygon(self, label):
        r"""
        Return the polygon with ``label`` in this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()
            sage: S = S.pyflatsurf().codomain()

            sage: S.polygon((1, 2, 3))
            Polygon(vertices=[(0, 0), (1, 0), (1, 1)])

        """
        label = Surface_pyflatsurf._normalize_label(label)

        from pyflatsurf import flatsurf
        half_edges = (flatsurf.HalfEdge(half_edge) for half_edge in label)
        vectors = [self._flat_triangulation.fromHalfEdge(half_edge) for half_edge in half_edges]

        from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
        vector_space_conversion = VectorSpaceConversion.from_pyflatsurf_from_elements(vectors)
        vectors = [vector_space_conversion.section(vector) for vector in vectors]

        from flatsurf.geometry.polygon import Polygon
        return Polygon(edges=vectors)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon and edge that is across from ``edge`` of ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().codomain()
            sage: S = S.pyflatsurf().codomain()

            sage: S.opposite_edge((1, 2, 3), 0)
            ((-3, -1, -2), 1)
            sage: S.opposite_edge((1, 2, 3), 1)
            ((-3, -1, -2), 2)
            sage: S.opposite_edge((1, 2, 3), 2)
            ((-3, -1, -2), 0)

        """
        label = Surface_pyflatsurf._normalize_label(label)

        from pyflatsurf import flatsurf
        half_edge = flatsurf.HalfEdge(label[edge])

        opposite_half_edge = -half_edge

        opposite_label = self._flat_triangulation.face(opposite_half_edge)
        opposite_label = Surface_pyflatsurf._normalize_label(opposite_label)

        return opposite_label, opposite_label.index(opposite_half_edge.id())
