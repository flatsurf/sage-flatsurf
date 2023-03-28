r"""
Deformations between Surfaces

.. NOTE::

    We call this a deformation because it might not actually be a map between
    surfaces.

    For all practical purposes, one should think of these as maps. However,
    they might not be maps on the points of the surface or really on anything
    in some extreme cases.

EXAMPLES:

We can use deformations to follow a surface through a retriangulation process::

    sage: from flatsurf import translation_surfaces
    sage: S = translation_surfaces.regular_octagon().underlying_surface()
    sage: deformation = S.subdivide_edges(3)
    sage: deformation = deformation.codomain().subdivide() * deformation
    sage: T = deformation.codomain()

    sage: deformation
    Deformation from <flatsurf.geometry.surface.Surface_list object at 0x...> to <flatsurf.geometry.surface.Surface_dict object at 0x...>

We can then map points through the deformation::

    sage: from flatsurf.geometry.surface_objects import SurfacePoint
    sage: p = SurfacePoint(S, 0, (0, 0))
    sage: p
    Surface point with 8 coordinate representations

    sage: q = deformation(p)
    sage: q
    Surface point with 16 coordinate representations

A non-singular point::

    sage: p = SurfacePoint(S, 0, (1, 1))
    sage: p
    Surface point located at (1, 1) in polygon 0

    sage: q = deformation(p)
    sage: q
    Surface point with 2 coordinate representations

"""
# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2023 Julian Rüth
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
# ********************************************************************


class Deformation:
    # TODO: docstring
    def __init__(self, domain, codomain):
        # TODO: docstring

        self._domain = None
        if not domain.is_mutable():
            self._domain = domain

        self._codomain = None
        if not codomain.is_mutable():
            self._codomain = codomain

    def domain(self):
        # TODO: docstring
        if self._domain is None:
            raise ValueError("cannot determine the domain of a deformation if the domain is mutable")

        return self._domain

    def codomain(self):
        # TODO: docstring
        if self._codomain is None:
            raise ValueError("cannot determine the codomain of a deformation if the codomain is mutable")

        return self._codomain

    def __call__(self, x):
        # TODO: docstring
        from flatsurf.geometry.surface_objects import SurfacePoint
        if isinstance(x, SurfacePoint):
            return self._image_point(x)

        from flatsurf.geometry.homology import SimplicialHomologyClass
        if isinstance(x, SimplicialHomologyClass):
            return self._image_homology(x)

        from flatsurf.geometry.cohomology import SimplicialCohomologyClass
        if isinstance(x, SimplicialCohomologyClass):
            return self._image_cohomology(x)

        raise NotImplementedError(f"cannot map a {type(x)} through a deformation yet")

    def _image_point(self, p):
        # TODO: docstring
        raise NotImplementedError(f"a {type(self)} cannot compute the image of a point yet")

    def _image_homology(self, γ):
        # TODO: docstring
        raise NotImplementedError(f"a {type(self)} cannot compute the image of a homology class yet")

    def _image_cohomology(self, f):
        # TODO: docstring
        from sage.all import ZZ, matrix

        from flatsurf.geometry.homology import SimplicialHomology

        # TODO: Make homology/cohomology work without having an explicit TranslationSurface (but only a Surface.)
        from flatsurf import TranslationSurface

        domain_homology = SimplicialHomology(TranslationSurface(self.domain()))
        codomain_homology = SimplicialHomology(TranslationSurface(self.codomain()))

        M = matrix(ZZ, [
            [gamma.coefficient(gen) for gen in codomain_homology.gens()]
            for gamma in [self(gen) for gen in domain_homology.gens()]
        ])

        values = M.solve_right([f(gen) for gen in codomain_homology.gens()])
        return self.codomain({gen: value for (gen, value) in zip(codomain_homology.gens(), values)})

    def __mul__(self, other):
        # TODO: docstring
        return CompositionDeformation(self, other)

    def __repr__(self):
        # TODO: docstring
        return f"Deformation from {self.domain()} to {self.codomain()}"


class CompositionDeformation(Deformation):
    # TODO: docstring
    def __init__(self, lhs, rhs):
        # TODO: docstring
        super().__init__(rhs.domain(), lhs.codomain())

        self._lhs = lhs
        self._rhs = rhs

    def __call__(self, x):
        # TODO: docstring
        return self._lhs(self._rhs(x))


class SubdivideDeformation(Deformation):
    # TODO: docstring
    def _image_point(self, p):
        # TODO: docstring
        from flatsurf.geometry.pyflatsurf_conversion import FlatTriangulationConversion
        to_pyflatsurf = FlatTriangulationConversion.to_pyflatsurf(domain=self.codomain())

        label = next(iter(p.labels()))
        coordinates = next(iter(p.coordinates(label))) - self.domain().polygon(label).vertex(0)

        from flatsurf.geometry.pyflatsurf_conversion import VectorSpaceConversion
        to_pyflatsurf_vector = VectorSpaceConversion.to_pyflatsurf(coordinates.parent())

        def ccw(v, w):
            r"""
            Return whether v->w describe a non-clockwise turn.
            """
            return v[0] * w[1] >= w[0] * v[1]

        initial_edge = ((label, 0), 0)
        mid_edge = ((label, self.codomain().num_polygons() - 1), 1)

        face = mid_edge if ccw(self.codomain().polygon(mid_edge[0]).edge(mid_edge[1]), coordinates) else initial_edge
        face = to_pyflatsurf(face)
        coordinates = to_pyflatsurf_vector(coordinates)

        import pyflatsurf

        p = pyflatsurf.flatsurf.Point[type(to_pyflatsurf.codomain())](to_pyflatsurf.codomain(), face, coordinates)

        return to_pyflatsurf.section(p)


class SubdivideEdgesDeformation(Deformation):
    # TODO: docstring
    def _image_point(self, p):
        r"""
        Return the image of ``p`` in the codomain of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon().underlying_surface()
            sage: deformation = S.subdivide_edges(2)

            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: p = SurfacePoint(S, 0, (1, 1))

            sage: deformation._image_point(p)
            Surface point located at (1, 1) in polygon 0

        """
        from flatsurf.geometry.surface_objects import SurfacePoint

        label = next(iter(p.labels()))
        return SurfacePoint(self.codomain(), label, next(iter(p.coordinates(label))))


class DelaunayDeformation(Deformation):
    # TODO: docstring
    def __init__(self, domain, codomain, flip_sequence):
        # TODO: docstring
        super().__init__(domain, codomain)
        self._flip_sequence = flip_sequence
