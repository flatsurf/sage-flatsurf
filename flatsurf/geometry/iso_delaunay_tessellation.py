r"""
EXAMPLES::

    sage: from flatsurf import translation_surfaces
    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: s = translation_surfaces.veech_double_n_gon(5)
    sage: IsoDelaunayTessellation(s)
    IsoDelaunay Tessellation of TranslationSurface built from 2 polygons

TODO: Write something smart here.
TODO: Implement this surface in sage-flatsurf.

::

    sage: def prym3(w, h, t, e):
    ....:     from flatsurf import Surface_dict, TranslationSurface, polygons
    ....:     k.<rtD> = QuadraticField(e**2 + 8 * w * h)
    ....:     lmbd = (e + rtD)/2
    ....:
    ....:     square = polygons.square(lmbd, field=k)
    ....:     parallelogram_top = polygons((w - lmbd, 0), (lmbd, 0), (t, h), (-lmbd, 0), (lmbd - w, 0), (-t, -h), field=k)
    ....:     parallelogram_bottom = polygons((lmbd, 0), (w - lmbd, 0), (t, h), (lmbd - w, 0), (-lmbd, 0), (-t, -h), field=k)
    ....:
    ....:     surface = Surface_dict(base_ring=k)
    ....:
    ....:     surface.add_polygon(parallelogram_bottom, label=0)
    ....:     surface.add_polygon(square, label=1)
    ....:     surface.add_polygon(parallelogram_top, label=2)
    ....:
    ....:     surface.change_base_label(1)
    ....:
    ....:     surface.change_edge_gluing(0, 0, 2, 3)
    ....:     surface.change_edge_gluing(0, 1, 0, 3)
    ....:     surface.change_edge_gluing(0, 2, 0, 5)
    ....:     surface.change_edge_gluing(0, 4, 1, 0)
    ....:     surface.change_edge_gluing(1, 1, 1, 3)
    ....:     surface.change_edge_gluing(1, 2, 2, 1)
    ....:     surface.change_edge_gluing(2, 0, 2, 4)
    ....:     surface.change_edge_gluing(2, 2, 2, 5)
    ....:
    ....:     return TranslationSurface(surface)
    ....:
    sage: s = prym3(5, 1, 0, 0)
    sage: t = s.delaunay_triangulation(in_place=False)

    sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
    sage: idt = IsoDelaunayTessellation(t)
    sage: idt.explore()
    sage: idt.plot()
    Graphics object consisting of 218 graphics primitives

REFERENCES:

TODO: Write a text citing some references.

.. [JB2009] \J. Bowman, "Flat Structures and Complex Structures in Teichmüller
            Theory", PhD Thesis, https://hdl.handle.net/1813/13979

"""
# *********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022 Sam Freedman
#                      2022 Julian Rüth
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
# *********************************************************************

from sage.structure.parent import Parent
from flatsurf.geometry.hyperbolic import HyperbolicPlane
from sage.misc.cachefunc import cached_method


class IsoDelaunayTessellation(Parent):
    def __init__(self, surface):
        r"""
        TESTS:

        TODO: There's something wrong here. The IDT has a single polygon with
        three edges but only two of them are glued to each other::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.square_torus()
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()

        ::

            sage: S = SymmetricGroup(4)
            sage: r = S('(1)(2,3)(4)')
            sage: u = S('(1,2)(3,4)')
            sage: s = translation_surfaces.origami(r, u)
            sage: idt.explore()

        """
        from sage.all import Graph

        self._surface_original = surface

        self._hyperbolic_plane = HyperbolicPlane(surface.base_ring())

        self._surface = self._nondegenerate_delaunay_triangulation(surface)
        self._surface.set_immutable()

        self._dual_graph = Graph(multiedges=True, loops=True)
        self._surface_classes = {}

        self._ensure_dual_graph_vertex(
            self.root(), self._surface, self.root().edges()[0]
        )

    def _repr_(self):
        return f"IsoDelaunay Tessellation of {self._surface_original}"

    def explore(self, limit=None, tessellation_face=None):
        r"""
        Explore the dual graph of the IsoDelaunay tessellation up to the combinatorial ``limit`` where you start from ``vertex`` and then first cross ``edge``.
        When ``vertex`` is ``None``, start exploring from the vertex bounded by ``edge``.
        When ``edge`` is ``None``, explore all edges of the IsoDelaunay regions bounding ``vertex``.
        When both ``vertex`` and ``edge`` are ``None``, start exploring from a vertex that contains the point ``i``.


        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()

        """
        from sage.all import oo

        limit = oo if limit is None else limit
        if limit <= 0:
            return

        if tessellation_face is None:
            tessellation_face = self.root()

        if tessellation_face not in self._dual_graph:
            raise ValueError(
                "tessellation_face must be a polygon of the explored tessellation"
            )

        queued_tessellation_faces = set()

        from collections import deque

        queue = deque([(0, tessellation_face)])
        queued_tessellation_faces.add(tessellation_face)

        while queue:
            distance, tessellation_face = queue.popleft()

            for tessellation_edge in tessellation_face.edges():
                other_tessellation_face = self._explore(
                    tessellation_face, tessellation_edge
                )
                if (
                    distance + 1 < limit
                    and other_tessellation_face not in queued_tessellation_faces
                ):
                    queued_tessellation_faces.add(other_tessellation_face)
                    queue.append((distance + 1, other_tessellation_face))

    def _explore(self, tessellation_face, tessellation_edge):
        assert tessellation_face.dimension() == 2
        assert tessellation_edge.dimension() == 1

        cross_tessellation_face = self._cross(tessellation_face, tessellation_edge)

        # nothing to do if we have already explored across this polygon edge
        if cross_tessellation_face is not None:
            return cross_tessellation_face

        cross_tessellation_face, target_triangulation = self._develop(
            tessellation_face, tessellation_edge
        )

        target_triangulation.set_immutable()
        (
            cross_tessellation_face,
            cross_tessellation_edge,
            isometry,
        ) = self._ensure_dual_graph_vertex(
            cross_tessellation_face, target_triangulation, -tessellation_edge
        )

        self._dual_graph.add_edge(
            tessellation_face,
            cross_tessellation_face,
            label=(
                {tessellation_edge, cross_tessellation_edge},
                (isometry, tessellation_edge),
            ),
        )

        return cross_tessellation_face

    def _develop(self, tessellation_face, tessellation_edge, source_triangulation=None):
        if source_triangulation is None:
            _, source_triangulation = self._dual_graph.get_vertex(tessellation_face)

        target_triangulation = source_triangulation.copy()

        while True:
            for triangulation_edge in target_triangulation.edge_iterator():
                half_plane = self._half_plane(target_triangulation, triangulation_edge)
                if half_plane is None:
                    continue
                if half_plane.boundary() == tessellation_edge.geodesic():
                    target_triangulation = target_triangulation.triangle_flip(
                        *triangulation_edge
                    )
                    break

                # glue triangles across triangulation_edge that flips to get mock Delaunay cells
                # solve for self-isomorphism on level of mock cells?
            else:
                break

        cross_tessellation_face = self._hyperbolic_plane.polygon(
            self._iso_delaunay_region(target_triangulation)
        )

        assert (
            -tessellation_edge in cross_tessellation_face.edges()
        ), f"edge {-tessellation_edge} is not in the polygon {cross_tessellation_face} after crossing {tessellation_edge} from {tessellation_face}"

        return cross_tessellation_face, target_triangulation

    def _cross(self, tessellation_face, tessellation_edge):
        r"""
        Return the hyperbolic polygon obtained by crossing from
        ``tessellation_face`` over ``tessellation_edge``.

        This method returns a polygon that forms part of the
        :meth:`fundamental_domain`. If we have not explored what is on the
        other side of ``tesselation_edge`` yet, returns ``None``.
        """
        mod, _ = self._dual_graph.get_vertex(tessellation_face)

        if mod is not None:
            edges = list(tessellation_face.edges())
            edge = set(edges[edges.index(tessellation_edge) % mod :: mod])
        else:
            edge = {tessellation_edge}

        for v, w, (edges, isometry) in self._dual_graph.edges(
            tessellation_face, labels=True, sort=False
        ):
            if any(e in edges for e in edge):
                if v == tessellation_face:
                    return w
                if w == tessellation_face:
                    return v
                assert False

        return None

    def insert_orbifold_points(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.mcmullen_genus2_prototype(1, 1, 0, -1)
            sage: t = s.delaunay_triangulation(in_place=False)
            sage: idt = IsoDelaunayTessellation(t)
            sage: z = idt._hyperbolic_plane(i)
            sage: idt.explore()

            sage: idt._dual_graph.vertices(sort=True)
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0}]
            sage: idt.insert_orbifold_points()
            sage: idt._dual_graph.vertices(sort=True)
            [{2*l*(x^2 + y^2) + (-8*l + 4)*x ≥ 0} ∩ {(8*l - 4)*x - 4*l + 2 ≤ 0} ∩ {x ≥ 0} ∪ {I}]

        """
        # TODO: Should this mutate the tessellation or create a copy instead?
        for tessellation_face in self._dual_graph.vertices(sort=False):
            for (
                source_tessellation_face,
                target_tessellation_face,
                (tessellation_edges, isometry),
            ) in list(
                self._dual_graph.edges(tessellation_face, labels=True, sort=False)
            ):
                # crossing edge of the polygon cycles back to the very edge in the
                # polygon, so there is an orbifold point on that edge.
                # We patch the polygon by inserting a marked point.
                if len(tessellation_edges) == 1:
                    assert source_tessellation_face == target_tessellation_face
                    tessellation_face, _ = self._insert_orbifold_point(
                        tessellation_face, next(iter(tessellation_edges)).midpoint()
                    )

        # TODO: Insert orbifold points in the interior of a polygon, i.e., the ones detected with isomorphism()

    def _insert_orbifold_point(self, tessellation_face, point):
        r"""
        Insert ``point`` as a marked point on an edge of the polygon ``vertex``
        (and update the graph representing the explored part of the fundamental
        domain.)

        Return the new polygon and the edges adjacent to point.
        """
        for tessellation_edge in tessellation_face.edges():
            if point in tessellation_edge:
                break
        else:
            assert False

        tessellation_face_with_marked_vertices = tessellation_face.parent().polygon(
            tessellation_face.half_spaces(),
            marked_vertices=tuple(tessellation_face.vertices()) + (point,),
        )

        assert tessellation_face_with_marked_vertices != tessellation_face

        # TODO: This breaks self._surface_classes currently.
        self._dual_graph.add_vertex(tessellation_face_with_marked_vertices)
        self._dual_graph.set_vertex(
            tessellation_face_with_marked_vertices,
            self._dual_graph.get_vertex(tessellation_face),
        )
        for (
            source_tessellation_face,
            target_tessellation_face,
            (tessellation_edges, isometry),
        ) in list(self._dual_graph.edges(tessellation_face, labels=True, sort=False)):
            self._dual_graph.delete_edge(
                source_tessellation_face, target_tessellation_face, tessellation_edges
            )
            if source_tessellation_face == tessellation_face:
                source_tessellation_face = tessellation_face_with_marked_vertices
            if target_tessellation_face == tessellation_face:
                target_tessellation_face = tessellation_face_with_marked_vertices

            self._dual_graph.add_edge(
                source_tessellation_face,
                target_tessellation_face,
                label=(tessellation_edges, isometry),
            )

        self._dual_graph.delete_vertex(tessellation_face)

        for e in tessellation_face_with_marked_vertices.edges():
            if e.start() == tessellation_edge.start():
                edge_before_marked_point = e
                break
        else:
            assert False

        for e in tessellation_face_with_marked_vertices.edges():
            if e.end() == edge_before_marked_point.end():
                edge_after_marked_point = e
                break
        else:
            assert False

        return tessellation_face_with_marked_vertices, (
            edge_before_marked_point,
            edge_after_marked_point,
        )

    @cached_method
    def _to_pyflatsurf(self, triangulation):
        # TODO: Clear this cache when explore() is finished.
        from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf

        return to_pyflatsurf(triangulation)

    @cached_method
    def _automorphisms_quotient(self):
        r"""
        Return the number of automorphisms that correspond to the matrix
        identity and -identity, i.e., the subgroup of the automorphisms that
        consists just of relabelings or deformations that do not affect the
        hyperbolic picture.
        """
        surface = self._to_pyflatsurf(self._surface)

        from pyflatsurf import flatsurf

        S = type(surface)
        unlabeled_equivalence = flatsurf.Equivalence[S].unlabeled()
        order = flatsurf.EquivalenceClass[S](
            surface, unlabeled_equivalence
        ).automorphisms()

        if unlabeled_equivalence.isomorphic(
            surface, surface.applyMatrix(-1, 0, 0, -1).codomain()
        ):
            order *= 2

        return order

    def _ensure_dual_graph_vertex(self, tessellation_face, surface, tessellation_edge):
        r"""
        Return vertex and edge of hyperbolic polygon TODO
        """
        for tessellation_face_ in self._dual_graph:
            if tessellation_face_ == tessellation_face:
                return tessellation_face_, tessellation_edge, False

        surface_ = self._to_pyflatsurf(surface)

        from pyflatsurf import flatsurf

        S = type(surface_)
        linear_equivalence = flatsurf.Equivalence[S].linear()
        clazz = flatsurf.EquivalenceClass[S](surface_, linear_equivalence)

        if clazz in self._surface_classes:
            tessellation_face_ = self._surface_classes[clazz]

            # TODO: We don't need to run the isomorphism() machinery to recover the isomorphism here.
            isomorphism = None

            def capture_matrix(a, b, c, d):
                if a * d - b * c != 1:
                    return False
                if hasattr(a, "parent"):
                    from pyeantic import RealEmbeddedNumberField

                    k = RealEmbeddedNumberField(a.parent())
                    a, b, c, d = k(a), k(b), k(c), k(d)
                    k = k.number_field
                    a, b, c, d = k(a), k(b), k(c), k(d)
                else:
                    from sage.all import QQ

                    a, b, c, d = QQ(a), QQ(b), QQ(c), QQ(d)

                nonlocal isomorphism
                isomorphism = (a, b, c, d)
                return True

            # TODO: This is not correct after insert_orbifold_points() has been called.
            assert len(tessellation_face_.edges()) == len(tessellation_face.edges())

            surface__ = self._to_pyflatsurf(
                self._dual_graph.get_vertex(tessellation_face_)[1]
            )

            assert surface_.isomorphism(
                surface__, filter_matrix=capture_matrix
            ).has_value()
            assert isomorphism is not None

            a, b, c, d = isomorphism
            from sage.all import matrix

            mob = matrix(self.base_ring(), 2, [a, -b, -c, d])
            image_edge = tessellation_edge.apply_isometry(mob, model="half_plane")

            assert image_edge in tessellation_face_.edges()
            return tessellation_face_, image_edge, mob

        # We have to add a new vertex.

        # First, we check if the new face has self-symmetries.
        assert clazz.automorphisms() % self._automorphisms_quotient() == 0
        order = clazz.automorphisms() // self._automorphisms_quotient()

        if order != 1:
            assert len(tessellation_face.edges()) % order == 0
            mod = len(tessellation_face.edges()) // order

            # Add the new vertex.
            self._dual_graph.add_vertex(tessellation_face)
            assert surface is not None
            self._dual_graph.set_vertex(tessellation_face, (mod, surface))
            self._surface_classes[clazz] = tessellation_face
            return tessellation_face, tessellation_edge, None

        # Add the new vertex.
        self._dual_graph.add_vertex(tessellation_face)
        assert surface is not None
        self._surface_classes[clazz] = tessellation_face
        self._dual_graph.set_vertex(tessellation_face, (None, surface))
        return tessellation_face, tessellation_edge, None

    def is_vertex(self, translation_surface):
        r"""
        Return whether this is a vertex.
        """
        raise NotImplementedError

    def root(self):
        r"""
        Return the tessellation face from which we started to build the IsoDelaunay tessellation.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.root()
            {(6*a + 8)*(x^2 + y^2) + (-20*a - 28)*x + 14*a + 20 ≥ 0} ∩ {(2*a + 3)*(x^2 + y^2) + (-4*a - 6)*x - 2*a - 3 ≤ 0} ∩ {(4*a + 6)*(x^2 + y^2) - 4*a - 6 ≥ 0}
        """
        from sage.all import I

        return self.face(I)

    def point(self, surface):
        r"""
        Return the point in the :class:`HyperbolicPlane` corresponding to ``surface``.

        INPUT:

        - ``surface`` -- A surface in the SL(2, R)-orbit of the defining surface
        """
        raise NotImplementedError

    def edge(self, translation_surface):
        r"""
        Return the unoriented tessellation edge this translation surface is on.
        """
        raise NotImplementedError

    @classmethod
    def _face(cls, surface, point):
        A = cls._point_to_matrix(point)
        A_T = surface.apply_matrix(A, in_place=False).delaunay_triangulation(
            in_place=False
        )
        T = A_T.apply_matrix(~A, in_place=False)
        T.set_immutable()
        half_planes = cls._iso_delaunay_region(T)
        iso_delaunay_region = point.parent().polygon(half_planes)

        if iso_delaunay_region.dimension() < 2:
            return None, T

        return iso_delaunay_region, T

    def face(self, point):
        r"""
        Return a tessellation face containing the hyperbolic ``point``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.face(i)
            {(6*a + 8)*(x^2 + y^2) + (-20*a - 28)*x + 14*a + 20 ≥ 0} ∩ {(2*a + 3)*(x^2 + y^2) + (-4*a - 6)*x - 2*a - 3 ≤ 0} ∩ {(4*a + 6)*(x^2 + y^2) - 4*a - 6 ≥ 0}

        """
        point = self._hyperbolic_plane(point)
        x, y = point.coordinates()

        from sage.all import QQ

        shift = QQ(1) / 2

        while True:
            shifted = point.parent().point(x + shift, y, model="half_plane")
            face, pulled_back_triangulation = self._face(self._surface, shifted)
            if point in face:
                self._surface = pulled_back_triangulation
                return face
            shift /= 2

    def surface(self, point):
        r"""
        Return a translation surface corresponding to this ``point`` in the
        hyperbolic plane.
        """
        # This should be shared with the IDR code.
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return self._surface.apply_matrix(matrix([[1, x], [0, y]]), in_place=False)

    def fundamental_domain(self):
        r"""
        Return the fundamental domain as a polygon with edge pairings.
        """
        # The result would be a fundamental domain of the group generated by the
        # symmetries discovered so far.
        raise NotImplementedError

    def plot(self):
        # TODO: Why is plotting so slow?
        return sum(idr.plot() for idr in self._dual_graph)

    def polygon(self, vertex_or_edge):
        r"""
        Return the polygon obtained as the union of the
        triangles bounded by this edge and its reverse /
        by the edges adjacent to this vertex.
        """
        raise NotImplementedError

    def geodesics(self, vertex):
        r"""
        Return the geodesics through ``vertex``.

        These are the geodesics underlying the :meth:`segments`, with end points removed.
        """

        raise NotImplementedError

    @classmethod
    def _point_to_matrix(cls, point):
        from sage.all import matrix

        x, y = point.coordinates(model="half_plane")
        return matrix(2, [1, x, 0, y])

    @classmethod
    def _iso_delaunay_region(cls, triangulation):
        return [
            half_plane
            for edge in triangulation.edge_iterator()
            if (half_plane := cls._half_plane(triangulation, edge)) is not None
        ]

    @classmethod
    def _half_plane(cls, surface, edge):
        r"""
        Build the halfplane associated to ``edge``.
        If the hinge is not convex, return ``None``.
        """
        # check if hinge is convex
        if not surface.triangle_flip(*edge, test=True):
            return None

        index_triangle, index_edge = edge
        index_opposite_triangle, index_opposite_edge = surface.opposite_edge(edge)
        v0 = surface.polygon(index_opposite_triangle).edge(
            (index_opposite_edge + 1) % 3
        )
        v1 = surface.polygon(index_triangle).edge(index_edge)
        v2 = -surface.polygon(index_triangle).edge((index_edge - 1) % 3)

        x0, y0 = v0
        x1, y1 = v1
        x2, y2 = v2

        a = (
            x0 * y1 * y2 * (y2 - y1)
            + x1 * y0 * y2 * (y0 - y2)
            + x2 * y0 * y1 * (y1 - y0)
        )
        b = (
            x0 * y0 * (x1 * y2 - x2 * y1)
            + x1 * y1 * (x2 * y0 - x0 * y2)
            + x2 * y2 * (x0 * y1 - x1 * y0)
        )
        c = (
            x0 * x1 * y2 * (x0 - x1)
            + x1 * x2 * y0 * (x1 - x2)
            + x0 * x2 * y1 * (x2 - x0)
        )

        H = HyperbolicPlane(surface.base_ring())
        return H.geodesic(a, 2 * b, c, model="half_plane").left_half_space()

    @classmethod
    def _nondegenerate_delaunay_triangulation(cls, surface):
        r"""
        Return a Delaunay triangulation for `_surface` whose IDR is 2-dimensional
        """
        # Let M_i be starting surf with arbitrary DT producing potentially degen IDR
        # Consider T_e * M_i until it has a nondegenerate IDR
        # let M_{i+e} be the Delaunay triangulated surface with a nondegen IDR
        # then T_{-e} * M_{i + e} is a new DT of M_i with a nondegen IDR

        # If M_i has a degenerate IDR, do flips until not the case
        # e.g. only 8 / 250 of the triangulations of the regular octagon lead to a nondegen IDR

        surface = surface.delaunay_triangulation(in_place=False)

        H = HyperbolicPlane(surface.base_ring())

        shift = 0
        while True:
            shifted = H.point(shift, 1, model="half_plane")
            face, _ = cls._face(surface, shifted)

            from sage.all import I

            if face is not None and I in face:
                break

            if shift == 0:
                from sage.all import QQ

                shift = QQ(1)

            shift /= 2

        perturbation = cls._point_to_matrix(shifted)
        return (
            surface.apply_matrix(perturbation, in_place=False)
            .delaunay_triangulation(in_place=False)
            .apply_matrix(~perturbation, in_place=False)
        )

    def vertices(self):
        r"""
        Return the vertices of the completed fundamental domain.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()
            sage: len(idt.vertices())
            3

        """
        # TODO: Check that the fundamental domain has been computed.

        half_edges = set(
            (tessellation_face, tessellation_edge)
            for tessellation_face in self._dual_graph.vertices(sort=False)
            for tessellation_edge in tessellation_face.edges()
        )

        vertices = []

        while half_edges:
            vertex = [next(iter(half_edges))]

            while True:
                (tessellation_face, tessellation_edge) = vertex[-1]

                assert tessellation_edge in tessellation_face.edges()

                previous_tessellation_edge = tessellation_face.edges()[
                    tessellation_face.edges().index(tessellation_edge) - 1
                ]

                # TODO: Merge this with the code in _cross maybe.
                for (
                    source_tessellation_face,
                    target_tessellation_face,
                    (tessellation_edges, isometry),
                ) in self._dual_graph.edges(tessellation_face, labels=True, sort=False):

                    if previous_tessellation_edge in tessellation_edges:
                        if len(tessellation_edges) == 1:
                            cross_tessellation_edge = previous_tessellation_edge
                        else:
                            tessellation_edges = list(tessellation_edges)
                            if previous_tessellation_edge == tessellation_edges[0]:
                                cross_tessellation_edge = tessellation_edges[1]
                            else:
                                cross_tessellation_edge = tessellation_edges[0]
                    else:
                        continue

                    cross_tessellation_edge = (
                        target_tessellation_face
                        if source_tessellation_face == tessellation_face
                        else source_tessellation_face,
                        cross_tessellation_edge,
                    )

                    vertex.append(cross_tessellation_edge)
                    half_edges.remove(cross_tessellation_edge)
                    break
                else:
                    assert (
                        False
                    ), f"{previous_tessellation_edge} of {tessellation_face} not glued to another edge"

                if vertex[0] == vertex[-1]:
                    vertex.pop()
                    vertices.append(vertex)
                    break

        return vertices

    def angle(self, vertex, algorithm="develop"):
        r"""
        Return the total angle at ``vertex`` divided by 2π.

        The returned angle will be of the form 1/n or 0.

        INPUT:

        - ``vertex`` -- a vertex description as provided by :meth:`vertices`.

        ALGORITHM:

        We arrange the polygons in ``vertex`` by gluing them at their that
        touch the vertex. We copy this widget, rotate it and glue it to itself.
        We repeat the process until we have completed a full circle around the
        vertex. If we needed n copies, then the angle is 1/n.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()
            sage: vertices = idt.vertices()
            sage: angles = sorted([idt.angle(vertex) for vertex in vertices])
            sage: angles
            [0, 0, 1/4]

        """
        R = self._surface.base_ring()

        if vertex[0][1].start().is_ideal():
            return 0

        def next_edge(polygon, edge):
            edges = polygon.edges()
            return -edges[edges.index(edge) - 1]

        if algorithm == "isometry":
            # TODO: This is wrong in the above example.
            # Rotate the polygons so that they are joined at their edges.
            widget = []
            for polygon, edge in vertex:
                if not widget:
                    widget.append((polygon, edge))
                    continue

                previous_polygon, previous_edge = widget[-1]

                # Determine the edge in the previous polygon that is glued to edge
                # so we know by how much we need to rotate polygon to make it line
                # up with previous_polygon.
                previous_edge = next_edge(previous_polygon, previous_edge)

                from sage.all import AA

                isometry = (
                    edge.change_ring(AA)
                    .isometry(previous_edge.change_ring(AA))
                    .change_ring(R)
                )
                polygon = polygon.apply_isometry(isometry)
                edge = previous_edge

                widget.append((polygon, edge))

            # Now all the polygons are glued correctly into one "widget". We only
            # need the first and the last half edge of this glued widget.
            start = widget[0][1]
            end = next_edge(*widget[-1])

            # We determine how many copies of widget we need to go around vertex in a full circle.
            copies = 1

            from sage.all import AA

            rotation = (
                start.change_ring(AA).isometry(end.change_ring(AA)).change_ring(R)
            )
            while start != end:
                copies += 1
                end = end.apply_isometry(rotation)

            from sage.all import QQ

            return QQ(1) / copies

        elif algorithm == "develop":
            start = vertex[0]
            current = start
            triangulation = None
            polygons = [start[0]]
            while True:
                cross_edge = next_edge(*current)

                if cross_edge == start[1]:
                    break

                cross_polygon, triangulation = self._develop(
                    current[0], -cross_edge, triangulation
                )
                polygons.append(cross_polygon)
                current = cross_polygon, cross_edge

            assert len(polygons) % len(vertex) == 0

            from sage.all import QQ

            return QQ(1) / (len(polygons) // len(vertex))
        else:
            raise ValueError

    def orbifold_euler_characteristic(self):
        r"""
        Following Farb-Margalit "Primer on Mapping Class Groups",
        Section 7.2.2.
        Y is a 2-dimensional hyperbolic orbifold
        having signature (g; p_1, p_2, ..., p_m)
        chi_orb = (2 - 2g) - m + sum 1/p_i 
        """
        raise NotImplementedError

    def genus(self):
        r"""
        2 - 2g = vertices - edges + faces
        """

        def step(edge, face):
            r"""
            Walk CCW about source vertex of edge, stepping into next face
            """            
            idx_edge = list(face.edges()).index(edge) 
            edge_turn = face.edges()[idx_edge - 1]
            for (source_tessellation_face, target_tessellation_face, (tessellation_edges, isometry)) in list(self._dual_graph.edges(face, labels=True, sort=False)):

                if edge_turn in tessellation_edges:
                    face = source_tessellation_face if face == target_tessellation_face else target_tessellation_face
                    if len(tessellation_edges) == 1:
                        return edge_turn, face

                    tessellation_edges = list(tessellation_edges)
                    return (tessellation_edges[0] if edge_turn != tessellation_edges[0] else tessellation_edges[1]), face
            assert False, f"{edge_turn} not found"

        nfaces = len(self._dual_graph)
        nedges = self._dual_graph.num_edges()
        nvertices = 0

        edges_seen = set()
        for face in self._dual_graph.vertices(sort=False):
            for edge in face.edges():
                if edge in edges_seen:
                    continue
                nvertices += 1
                edges_seen.add(edge)

                next_edge, next_face = step(edge, face)
                while next_edge != edge:
                    assert next_edge not in edges_seen
                    edges_seen.add(next_edge)
                    next_edge, next_face = step(next_edge, next_face)

        print(nvertices,nedges,nfaces)
        chi = nvertices - nedges + nfaces
        return (2 - chi)//2

    def cusps(self):
        r"""
        Return the punctures of the quotient H^2/Gamma

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()
            sage: idt.cusps()

        """
        return [
            vertex_equivalence_class
            for vertex_equivalence_class in self.vertices()
            if vertex_equivalence_class[0][1].start().is_ideal()
        ]

    def cusp_matrix(self, cusp):
        r"""
        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: from flatsurf.geometry.iso_delaunay_tessellation import IsoDelaunayTessellation
            sage: s = translation_surfaces.veech_2n_gon(4)
            sage: idt = IsoDelaunayTessellation(s)
            sage: idt.explore()
            sage: cusps = idt.cusps()
            sage: idt.cusp_matrix(cusps[0])
            [    -a  a + 1]
            [-a - 1  a + 2]
            sage: idt.cusp_matrix(cusps[1])
            [     1/2*a -3/2*a - 2]
            [     1/2*a -1/2*a - 2]
        """
        from sage.all import MatrixSpace

        matrix = MatrixSpace(self.base_ring(), 2, 2).one()

        for tessellation_face, tessellation_edge in cusp:
            dual_graph_edges = self._dual_graph.edges(
                sort=False, vertices=[tessellation_face], labels=True
            )
            for (
                _,
                __,
                (tessellation_edges, (isometry, source_edge)),
            ) in dual_graph_edges:
                if isometry is None:
                    continue
                if tessellation_edge in tessellation_edges:
                    if source_edge != tessellation_edge:
                        isometry = ~isometry
                    matrix *= isometry

        return matrix

    def cusp_direction(self, cusp):
        m = self.cusp_matrix(cusp)
        evns = m.eigenvectors_right()
        assert len(evns) == 1
        e, v, n = evns[0]
        assert len(v) == 1
        return v[0]

    def multitwist(self, cusp):
        v = self.cusp_direction(cusp)
        from flatsurf import GL2ROrbitClosure

        O = GL2ROrbitClosure(self._surface)
        decomposition = O.decomposition(v)
        components = decomposition.components()
        assert all(c.cylinder() for c in components)
        holonomies = [c.circumferenceHolonomy() for c in components]
        areas = [c.area() / 2 for c in components]
        inv_moduli = [
            (v.x() * v.x() + v.y() * v.y()) / area for v, area in zip(holonomies, areas)
        ]

        import pyeantic

        inv_moduli = [
            pyeantic.RealEmbeddedNumberField(x.parent())(x) for x in inv_moduli
        ]

        inv_moduli = [x.parent().number_field(x) for x in inv_moduli]

        from sage.all import QQ, lcm, matrix

        multiplicities = [QQ(x // inv_moduli[0]) for x in inv_moduli]
        denominator = lcm(rat.denominator() for rat in multiplicities)
        multiplicities = [denominator * rat for rat in multiplicities]
        c = lcm(multiplicities) * inv_moduli[0]

        # want to carry v = [x y] to span(e_1)
        # can use [x -y y x]^-1
        COB = ~matrix(2, [v[0], -v[1], v[1], v[0]])
        return ~COB * matrix(2, [1, c, 0, 1]) * COB

    def base_ring(self):
        return self._surface.base_ring()

    def orbifold_points(self, order=None):
        r"""
        Return the set of orbifold points, i.e., the fixed points of finite-order rotations in Gamma.

        When ``order = k``, return only the orbifold points with total angle ``2pi/k``.
        """
        # TODO: implement by developing in a circle around vertex
        pass

    def gens(self):
        r"""
        Return the generators of the Veech group that we see while exploring the isoDelaunay tessellation.
        """

        def twist(m):
            from sage.all import matrix

            ((a, b), (c, d)) = m
            return matrix([[a, -b], [-c, d]])

        for (
            _,
            _,
            (source_tessellation_edge, target_tessellation_edge),
        ) in self._dual_graph.edges(labels=True, sort=False):
            isometry = target_tessellation_edge[0]
            if isometry is not None:
                yield twist(isometry)

        for tessellation_face in self._dual_graph.vertices(sort=False):
            mod, _ = self._dual_graph.get_vertex(tessellation_face)
            edges = tessellation_face.edges()
            if mod is not None:
                yield twist(edges[0].isometry(edges[mod]))
