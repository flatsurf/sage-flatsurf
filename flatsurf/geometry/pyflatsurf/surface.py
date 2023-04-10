from flatsurf.geometry.surface import Surface


class Surface_pyflatsurf(Surface):
    r"""
    EXAMPLES::

        sage: from flatsurf import polygons
        sage: from flatsurf.geometry.surface import Surface_dict

        sage: S = Surface_dict(QQ)
        sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 0), (1, 1)]), label=0)
        0
        sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 1), (0, 1)]), label=1)
        1

        sage: S.set_edge_pairing(0, 0, 1, 1)
        sage: S.set_edge_pairing(0, 1, 1, 2)
        sage: S.set_edge_pairing(0, 2, 1, 0)

        sage: T = S.pyflatsurf().codomain()  # random output due to deprecation warnings

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

        sage: isinstance(T, Surface_pyflatsurf)
        True
        sage: TestSuite(T).run()  # TODO: not tested

    """

    def __init__(self, flat_triangulation):
        self._flat_triangulation = flat_triangulation

        # TODO: We have to be smarter about the ring bridge here.
        from flatsurf.geometry.pyflatsurf_conversion import RingConversion

        base_ring = RingConversion.from_pyflatsurf_from_flat_triangulation(flat_triangulation).domain()

        base_face = next(iter(flat_triangulation.faces()))
        base_label = tuple(sorted(half_edge.id() for half_edge in base_face))

        super().__init__(
            base_ring=base_ring, base_label=base_label, finite=True, mutable=False
        )

    def pyflatsurf(self):
        from flatsurf.geometry.deformation import IdentityDeformation
        return IdentityDeformation(self, self)

    @classmethod
    def _from_flatsurf(cls, surface):
        r"""
        Return an isomorphism to a :class:`Surface_pyflatsurf` built from
        ``surface``, i.e., represent ``surface`` in pyflatsurf wrapped as a
        :class:`Surface` for sage-flatsurf.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().underlying_surface()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: Surface_pyflatsurf._from_flatsurf(S)
            Deformation from mutable domain to FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 1), 2: (-1, 0), 3: (0, -1)}

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

        from flatsurf.geometry.pyflatsurf.deformation import Deformation_to_pyflatsurf

        return Deformation_to_pyflatsurf(surface, surface_pyflatsurf, to_pyflatsurf)

    def __repr__(self):
        r"""
        Return a printable representation of this surface, namely, print this
        surface as pyflatsurf would.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().underlying_surface()

            sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf
            sage: S.pyflatsurf().codomain()
            FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 1), 2: (-1, 0), 3: (0, -1)}

        """
        return repr(self._flat_triangulation)

    def apply_matrix(self, m):
        deformation = self._flat_triangulation.applyMatrix(*m.list())
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
            sage: S = translation_surfaces.square_torus().triangulate().underlying_surface()
            sage: S = S.pyflatsurf().codomain()

            sage: S.polygon((1, 2, 3))
            Polygon: (0, 0), (1, 1), (0, 1)

        """
        label = Surface_pyflatsurf._normalize_label(label)

        from pyflatsurf import flatsurf
        half_edges = (flatsurf.HalfEdge(half_edge) for half_edge in label)
        vectors = [self._flat_triangulation.fromHalfEdge(half_edge) for half_edge in half_edges]

        from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
        vector_space_conversion = VectorSpaceConversion.from_pyflatsurf_from_elements(vectors)
        vectors = [vector_space_conversion.section(vector) for vector in vectors]

        from flatsurf.geometry.polygon import polygons
        return polygons(edges=vectors)

    def opposite_edge(self, label, edge):
        r"""
        Return the polygon and edge that is across from ``edge`` of ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.square_torus().triangulate().underlying_surface()
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
