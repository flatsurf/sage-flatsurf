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
    sage: S = translation_surfaces.regular_octagon()
    sage: deformation = S.subdivide_edges(3)
    sage: deformation = deformation.codomain().subdivide() * deformation
    sage: T = deformation.codomain()

    sage: deformation
    Deformation from Translation Surface in H_2(2) built from a regular octagon to Translation Surface in H_2(2, 0^9) built from 8 isosceles triangles and 16 triangles

We can then map points through the deformation::

    sage: from flatsurf.geometry.surface_objects import SurfacePoint
    sage: p = SurfacePoint(S, 0, (0, 0))
    sage: p
    Vertex 0 of polygon 0

    sage: q = deformation(p)
    sage: q
    Vertex 0 of polygon (0, 0)

A non-singular point::

    sage: p = SurfacePoint(S, 0, (1, 1))
    sage: p
    Point (1, 1) of polygon 0

    sage: q = deformation(p)
    sage: q
    Point (1, 1) of polygon (0, 5)

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

from sage.misc.cachefunc import cached_method


class Deformation:
    # TODO: docstring
    def __init__(self, domain, codomain):
        # TODO: docstring

        self._domain = domain
        if domain is None or domain.is_mutable():
            self._domain = None

        self._codomain = codomain
        if codomain is None or codomain.is_mutable():
            self._codomain = None

    def domain(self):
        # TODO: docstring
        if self._domain is None:
            raise ValueError("cannot determine the domain of this deformation")

        return self._domain

    def codomain(self):
        # TODO: docstring
        if self._codomain is None:
            raise ValueError("cannot determine the codomain of this deformation")

        return self._codomain

    def __getattr__(self, name):
        if name in ["__cached_methods"]:
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        try:
            attr = getattr(self._codomain, name)
        except AttributeError:
            raise AttributeError(f"'{type(self)}' has no attribute '{name}'")

        import warnings
        warnings.warn(f"This methods returns a deformation instead of a surface. Use .codomain().{name} to access the surface instead of the deformation.")

        return attr

    def section(self):
        return SectionDeformation(self)

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

        from flatsurf.geometry.saddle_connection import SaddleConnection
        if isinstance(x, SaddleConnection):
            return self._image_saddle_connection(x)

        raise NotImplementedError(f"cannot map a {type(x)} through a deformation yet")

    def _image_point(self, p):
        # TODO: docstring
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of a point yet")

    def _image_saddle_connection(self, c):
        # TODO: docstring
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of a saddle connection yet")

    def _image_homology(self, γ):
        # TODO: docstring

        from flatsurf.geometry.homology import SimplicialHomology
        codomain_homology = SimplicialHomology(self.codomain())

        from sage.all import vector
        image = self._image_homology_matrix() * vector(γ._homology())

        homology, to_chain, to_homology = codomain_homology._homology()

        image = sum(coefficient * gen for (coefficient, gen) in zip(image, homology.gens()))

        return codomain_homology(to_chain(image))

    def _image_edge(self, label, edge):
        raise NotImplementedError(f"a {type(self).__name__} cannot compute the image of an edge yet")

    @cached_method
    def _image_homology_matrix(self):
        # TODO: docstring
        from flatsurf.geometry.homology import SimplicialHomology
        domain_homology = SimplicialHomology(self.domain())
        codomain_homology = SimplicialHomology(self.codomain())

        domain_gens = domain_homology.gens()
        codomain_gens = codomain_homology.gens()

        from sage.all import matrix, ZZ
        M = matrix(ZZ, len(codomain_gens), len(domain_gens), sparse=True)

        for x, domain_gen in enumerate(domain_gens):
            image = self._image_homology_gen(domain_gen)
            for y, codomain_gen in enumerate(codomain_gens):
                M[y, x] = image.coefficient(codomain_gen)

        M.set_immutable()
        return M

    def _image_homology_gen(self, gen):
        from flatsurf.geometry.homology import SimplicialHomology
        codomain_homology = SimplicialHomology(self.codomain())

        chain = gen._chain
        image = codomain_homology.zero()
        for label, edge in chain.support():
            coefficient = chain[(label, edge)]
            assert coefficient
            for step in self._image_edge(label, edge):
                image += coefficient * image.parent()(step)

        return codomain_homology(image)

    def _image_cohomology(self, f):
        # TODO: docstring
        from sage.all import ZZ, matrix

        from flatsurf.geometry.homology import SimplicialHomology

        domain_homology = SimplicialHomology(self.domain())
        codomain_homology = SimplicialHomology(self.codomain())

        M = matrix(ZZ, [
            [domain_homology_gen_image.coefficient(codomain_homology_gen) for codomain_homology_gen in codomain_homology.gens()]
            for domain_homology_gen_image in [self(domain_homology_gen) for domain_homology_gen in domain_homology.gens()]
        ])

        from sage.all import vector
        values = M.solve_right(vector([f(gen) for gen in domain_homology.gens()]))

        from flatsurf.geometry.cohomology import SimplicialCohomology
        codomain_cohomology = SimplicialCohomology(self.codomain())
        return codomain_cohomology({gen: value for (gen, value) in zip(codomain_homology.gens(), values)})

    def __mul__(self, other):
        # TODO: docstring
        return CompositionDeformation(self, other)

    def __repr__(self):
        # TODO: docstring
        return f"Deformation from {self.domain() if self._domain is not None else '?'} to {self.codomain() if self._codomain is not None else '?'}"


class IdentityDeformation(Deformation):
    # TODO: docstring
    def __init__(self, domain):
        super().__init__(domain, domain)

    def __call__(self, x):
        return x

    def _image_homology(self, γ):
        return γ

    def _image_edge(self, label, edge):
        return [(label, edge)]


class SectionDeformation(Deformation):
    def __init__(self, deformation):
        self._deformation = deformation
        super().__init__(deformation.codomain(), deformation.domain())

    def _image_homology_matrix(self):
        M = self._deformation._image_homology_matrix()
        return M.parent()(M.inverse())

    def _image_edge(self, label, edge):
        for (l, e) in self.codomain().edges():
            if self._deformation._image_edge(l, e) == (label, edge):
                return (l, e)

        raise NotImplementedError


class CompositionDeformation(Deformation):
    # TODO: docstring
    def __init__(self, lhs, rhs):
        # TODO: docstring
        super().__init__(rhs._domain, lhs._codomain)

        self._lhs = lhs
        self._rhs = rhs

    def __call__(self, x):
        # TODO: docstring
        return self._lhs(self._rhs(x))

    def _image_homology_matrix(self):
        return self._lhs._image_homology_matrix() * self._rhs._image_homology_matrix()

    def _image_edge(self, label, edge):
        return self._lhs._image_edge(*self._rhs._image_edge(label, edge))


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
        mid_edge = ((label, len(self.codomain().polygons()) // len(self.domain().polygons()) - 1), 1)

        face = mid_edge if ccw(self.codomain().polygon(mid_edge[0]).edge(mid_edge[1]), coordinates) else initial_edge
        face = to_pyflatsurf(face)
        coordinates = to_pyflatsurf_vector(coordinates)

        import pyflatsurf

        p = pyflatsurf.flatsurf.Point[type(to_pyflatsurf.codomain())](to_pyflatsurf.codomain(), face, coordinates)

        return to_pyflatsurf.section(p)

    def _image_homology(self, γ):
        # TODO: homology should allow to write this code more naturally somehow
        # TODO: docstring
        from flatsurf.geometry.homology import SimplicialHomology

        H = SimplicialHomology(self.codomain())
        C = H.chain_module(dimension=1)

        image = H()

        for gen in γ.parent().gens():
            for step in gen._path():
                codomain_gen = (step, 0)
                sgn = 1
                if codomain_gen not in H.simplices():
                    sgn = -1
                    codomain_gen = self.codomain().opposite_edge(*codomain_gen)
                image += sgn * γ.coefficient(gen) * H(C(codomain_gen))

        return image


class SubdivideEdgesDeformation(Deformation):
    # TODO: docstring

    def __init__(self, domain, codomain, parts):
        super().__init__(domain, codomain)
        self._parts = parts

    def _image_point(self, p):
        r"""
        Return the image of ``p`` in the codomain of this surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.regular_octagon()
            sage: deformation = S.subdivide_edges(2)

            sage: from flatsurf.geometry.surface_objects import SurfacePoint
            sage: p = SurfacePoint(S, 0, (1, 1))

            sage: deformation._image_point(p)
            Point (1, 1) of polygon 0

        """
        from flatsurf.geometry.surface_objects import SurfacePoint

        if not p.surface() == self.domain():
            raise ValueError

        label = next(iter(p.labels()))
        return SurfacePoint(self.codomain(), label, next(iter(p.coordinates(label))))

    # def _image_homology(self, γ):
    #     # TODO: homology should allow to write this code more naturally somehow
    #     # TODO: docstring
    #     from flatsurf.geometry.homology import SimplicialHomology

    #     H = SimplicialHomology(self.codomain())
    #     C = H.chain_module(dimension=1)

    #     image = H()

    #     for gen in γ.parent().gens():
    #         for step in gen._path():
    #             for p in range(self._parts):
    #                 codomain_gen = (step[0], step[1] * self._parts + p)
    #                 sgn = 1
    #                 if codomain_gen not in H.simplices():
    #                     sgn = -1
    #                     codomain_gen = self.codomain().opposite_edge(*codomain_gen)
    #                 image += sgn * γ.coefficient(gen) * H(C(codomain_gen))

    #     return image


class TrianglesFlipDeformation(Deformation):
    # TODO: docstring
    def __init__(self, domain, codomain, flip_sequence):
        # TODO: docstring
        super().__init__(domain, codomain)
        assert not self.domain().is_mutable()
        self._flip_sequence = flip_sequence

    @cached_method
    def _flip_deformation(self):
        domains = [self.domain()]

        for flip in self._flip_sequence:
            from flatsurf.geometry.surface import Surface_dict
            codomain = Surface_dict(surface=domains[-1], mutable=True)
            codomain.flip(*flip)
            codomain.set_immutable()
            domains.append(codomain)

        assert domains[-1] == self.codomain()
        domains.pop()
        domains.append(self.codomain())

        # TODO: Get rid of this trivial step.
        deformation = IdentityDeformation(self.domain())
        for i, flip in enumerate(self._flip_sequence):
            deformation = TriangleFlipDeformation(domains[i], domains[i + 1], flip) * deformation

        return deformation

    def _image_homology(self,  γ):
        return self._flip_deformation()._image_homology(γ)

    def _image_homology_matrix(self):
        return self._flip_deformation()._image_homology_matrix()

    def _image_point(self, p):
        return self._flip_deformation()(p)


class TriangleFlipDeformation(Deformation):
    # TODO: docstring
    def __init__(self, domain, codomain, flip):
        # TODO: docstring
        super().__init__(domain, codomain)
        self._flip = flip

    def _image_edge(self, label, edge):
        opposite_flip = self.domain().opposite_edge(*self._flip)

        def find_in_codomain(vector):
            candidates = [(l, e) for l in [self._flip[0], opposite_flip[0]] for e in [0, 1, 2]]
            candidates = [(l, e) for (l, e) in candidates if self.codomain().polygon(l).edge(e) == vector]
            assert candidates
            if len(candidates) > 1:
                raise NotImplementedError("cannot break ties on such a non-translation surface yet")

            return candidates[0]

        if label == self._flip[0]:
            if edge == self._flip[1]:
                return self._image_edge(*self.domain().opposite_edge(label, (edge + 1) % 3)) + self._image_edge(*self.domain().opposite_edge(label, (edge + 2) % 3))

            return [find_in_codomain(self.domain().polygon(label).edge(edge))]
        if label == opposite_flip[0]:
            if edge == opposite_flip[1]:
                return self._image_edge(*self.domain().opposite_edge(label, (edge + 1) % 3)) + self._image_edge(*self.domain().opposite_edge(label, (edge + 2) % 3))

            return [find_in_codomain(self.domain().polygon(label).edge(edge))]

        return [(label, edge)]

    def _image_point(self, p):
        # TODO: This is a dumb way to do this (and wrong for non-translation surfaces?)

        from flatsurf.geometry.surface_objects import SurfacePoint
        affected = {self._flip[0], self.domain().opposite_edge(*self._flip)[0]}

        for label in p.labels():
            if label not in affected:
                return SurfacePoint(self.codomain(), label, next(iter(p.coordinates(label))))

        for label in p.labels():
            polygon = self.domain().polygon(label)
            for edge in range(polygon.num_edges()):
                cross_label, cross_edge = self.domain().opposite_edge(label, edge)
                if cross_label in affected:
                    continue

                Δ = next(iter(p.coordinates(label))) - polygon.vertex(edge)
                recross_label, recross_edge = self.codomain().opposite_edge(cross_label, cross_edge)
                recross_polygon = self.codomain().polygon(recross_label)
                recross_position = recross_polygon.vertex(recross_edge) + Δ
                if not recross_polygon.get_point_position(recross_position).is_inside():
                    continue

                return SurfacePoint(self.codomain(), recross_label, recross_position)

        raise NotImplementedError


class TriangulationDeformation(Deformation):
    def _image_edge(self, label, edge):
        return [self._codomain._triangulation(label)[1][edge]]

    def _image_saddle_connection(self, connection):
        (label, edge) = connection.start_data()

        (label, edge) = self._image_edge(label, edge)[0]

        from flatsurf.geometry.euclidean import ccw
        while ccw(connection.direction(), -self._codomain.polygon(label).edges()[(edge - 1) % len(self._codomain.polygon(label).edges())]) <= 0:
            (label, edge) = self._codomain.opposite_edge(label, (edge - 1) % len(self._codomain.polygon(label).edges()))

        # TODO: This is extremely slow.
        from flatsurf.geometry.saddle_connection import SaddleConnection
        return SaddleConnection(self._codomain, (label, edge), direction=connection.direction())


class DelaunayDecompositionDeformation(Deformation):
    def __repr__(self):
        # TODO: docstring
        return f"Delaunay cell decomposition of {self._domain or '?'}"


class GL2RDeformation(Deformation):
    def __init__(self, domain, codomain, m):
        super().__init__(domain, codomain)

        from sage.all import matrix
        self._matrix = matrix(m, immutable=True)
