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

        base_label = map(int, next(iter(flat_triangulation.faces())))

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
            sage: Surface_pyflatsurf._from_flatsurf(S)
            Deformation from mutable domain to FlatTriangulationCombinatorial(vertices = (1, -3, 2, -1, 3, -2), faces = (1, 2, 3)(-1, -2, -3)) with vectors {1: (1, 1), 2: (-1, 0), 3: (0, -1)}

        """
        return repr(self._flat_triangulation)
