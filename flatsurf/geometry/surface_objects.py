r"""
Geometric objects on surfaces.

This includes singularities, saddle connections and cylinders.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2017-2020 W. Patrick Hooper
#                      2017-2020 Vincent Delecroix
#                           2023 Julian Rüth
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
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.plot.graphics import Graphics
from sage.plot.polygon import polygon2d
from sage.rings.qqbar import AA
from sage.structure.sage_object import SageObject
from sage.structure.element import Element

from flatsurf.geometry.similarity import SimilarityGroup


def Singularity(similarity_surface, label, v, limit=None):
    r"""
    Return the point of ``similarity_surface`` at the ``v``-th vertex of the
    polygon ``label``.

    If the surface is infinite, the ``limit`` can be set. In this case the
    construction of the singularity is successful if the sequence of vertices
    hit by passing through edges closes up in ``limit`` or less steps.

    EXAMPLES::

        sage: from flatsurf.geometry.similarity_surface_generators import TranslationSurfaceGenerators
        sage: s=TranslationSurfaceGenerators.veech_2n_gon(5)
        sage: from flatsurf.geometry.surface_objects import Singularity
        sage: sing=Singularity(s, 0, 1)
        doctest:warning
        ...
        UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
        sage: print(sing)
        Vertex 1 of polygon 0
        sage: TestSuite(sing).run()

    """
    import warnings

    warnings.warn(
        "Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead."
    )
    return similarity_surface.point(
        label, similarity_surface.polygon(label).vertex(v), limit=limit
    )


class SurfacePoint_base(Element):
    pass


class SurfacePoint(SurfacePoint_base):
    r"""
    A point on ``surface``.

    INPUT:

    - ``surface`` -- a similarity surface

    - ``label`` -- a polygon label for the polygon with respect to which the
      ``point`` coordinates can be made sense of

    - ``point`` -- coordinates of a point in the polygon ``label`` or the index
      of the vertex of the polygon with ``label``

    - ``ring`` -- a SageMath ring or ``None`` (default: ``None``); the
      coordinate ring for ``point``

    - ``limit`` -- an integer or ``None`` (default: ``None`` for an unlimited
      number of steps); if this is a singularity of the surface, then this
      limits the number of edges that are crossed to determine all the edges
      adjacent to that singularity. An error is raised if the limit is
      insufficient.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces

        sage: permutation = SymmetricGroup(2)('(1, 2)')
        sage: S = translation_surfaces.origami(permutation, permutation)

    A point can have a single representation with coordinates when it is in
    interior of a polygon::

        sage: S.point(1, (1/2, 1/2))
        Point (1/2, 1/2) of polygon 1

    A point can have two representations when it is in interior of an edge::

        sage: p = S.point(1, (1/2, 0))
        sage: q = S.point(2, (1/2, 1))
        sage: p == q
        True
        sage: p.coordinates(2)
        ((1/2, 1),)

    A point can have even more representations when it is a vertex::

        sage: S.point(1, (0, 0))
        Vertex 0 of polygon 1

    TESTS:

    Test that points can be created on disconnected surfaces::

        sage: from flatsurf import MutableOrientedSimilaritySurface, polygons
        sage: S = MutableOrientedSimilaritySurface(QQ)
        sage: S.add_polygon(polygons.square())
        0
        sage: S(0, 0)
        Vertex 0 of polygon 0

    """

    def __init__(self, surface, label, point, ring=None, limit=None):
        self._surface = surface

        if ring is None:
            ring = surface.base_ring()

        if ring is not surface.base_ring():
            import warnings

            warnings.warn(
                "the ring parameter is deprecated and will be removed in a future version of sage-flatsurf; define the surface over a larger ring instead so that this points' coordinates live in the base ring"
            )

        polygon = surface.polygon(label)

        from sage.all import ZZ

        if point in ZZ:
            point = surface.polygon(label).vertex(point)

        point = (ring**2)(point)
        point.set_immutable()

        position = polygon.get_point_position(point)

        if not position.is_inside():
            raise NotImplementedError(
                "point must be positioned within the polygon with the given label"
            )

        if position.is_in_interior():
            self._representatives = {(label, point)}
        elif position.is_in_edge_interior():
            self._representatives = {(label, point)}

            opposite_edge = surface.opposite_edge(label, position.get_edge())
            if opposite_edge is not None:
                cross_label, cross_edge = opposite_edge
                cross_point = surface.edge_transformation(label, position.get_edge())(point)
                cross_point.set_immutable()

                self._representatives.add((cross_label, cross_point))
        elif position.is_vertex():
            self._representatives = set()

            source_edge = position.get_vertex()

            def collect_representatives(label, source_edge, direction, limit):
                def rotate(label, source_edge, direction):
                    if direction == -1:
                        source_edge = (source_edge - 1) % len(surface.polygon(label).vertices())
                    opposite_edge = surface.opposite_edge(label, source_edge)
                    if opposite_edge is None:
                        return None
                    label, source_edge = opposite_edge
                    if direction == 1:
                        source_edge = (source_edge + 1) % len(surface.polygon(label).vertices())
                    return label, source_edge

                while True:
                    self._representatives.add((label, source_edge))

                    # Rotate to the next edge that is leaving at the vertex
                    rotated = rotate(label, source_edge, direction)
                    if rotated is None:
                        # Surface is disconnected
                        return

                    label, source_edge = rotated

                    if limit is not None:
                        limit -= 1
                        if limit < 0:
                            raise ValueError("number of edges at singularity exceeds limit")

                    if (label, source_edge) in self._representatives:
                        break

            # Collect respresentatives of edge walking clockwise and counterclockwise
            collect_representatives(label, source_edge, 1, limit)
            collect_representatives(label, source_edge, -1, limit)

            self._representatives = {
                (label, surface.polygon(label).vertex(vertex))
                for (label, vertex) in self._representatives
            }
        else:
            raise NotImplementedError

        self._representatives = frozenset(self._representatives)

        super().__init__(surface)

    def surface(self):
        r"""
        Return the surface containing this point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (1/2, 1/2))
            sage: p.surface() is S
            True

        """
        return self._surface

    def is_vertex(self):
        r"""
        Return whether this point is a singularity of the surface.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))

            sage: p.is_vertex()
            True

        """
        label, coordinates = self.representative()
        position = self.surface().polygon(label).get_point_position(coordinates)
        return position.is_vertex()

    def one_vertex(self):
        r"""
        Return a pair (l, v) from the equivalence class of this singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.one_vertex()  # random output: depends on the Python version
            doctest:warning
            ...
            UserWarning: one_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use (label, coordinates) = point.representative(); vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
            (2, 1)

        """
        import warnings

        warnings.warn(
            "one_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use (label, coordinates) = point.representative(); vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead"
        )
        label, coordinates = self.representative()
        vertex = (
            self.surface().polygon(label).get_point_position(coordinates).get_vertex()
        )
        return label, vertex

    def representatives(self):
        r"""
        Return the representatives of this point as pairs of polygon labels and
        coordinates.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.representatives()
            frozenset({(1, (0, 0)), (1, (1, 1)), (2, (0, 1)), (2, (1, 0))})

        """
        return self._representatives

    def representative(self):
        r"""
        Return a representative of this point, i.e., the first of
        :meth:`representatives`.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.representative()  # random output: depends on the Python version
            (2, (1, 0))

        """
        return next(iter(self.representatives()))

    def vertex_set(self):
        r"""
        Return the list of pairs (l, v) in the equivalence class of this singularity.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: list(p.vertex_set())  # random output: ordering depends on the Python version
            doctest:warning
            ...
            UserWarning: vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead
            [(2, 1), (1, 2), (1, 0), (2, 3)]

        """
        import warnings

        warnings.warn(
            "vertex_set() is deprecated and will be removed in a future version of sage-flatsurf; use representatives() and then vertex = surface.polygon(label).get_point_position(coordinates).get_vertex() instead"
        )

        return [
            (
                label,
                self.surface()
                .polygon(label)
                .get_point_position(coordinates)
                .get_vertex(),
            )
            for label, coordinates in self.representatives()
        ]

    def contains_vertex(self, label, v=None):
        r"""
        Checks if the pair ``(label, v)`` is in the equivalence class returning
        true or false. If ``v`` is ``None``, the both the pair ``(label, v)``
        is passed as a single parameter in ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.contains_vertex((1, 0))
            doctest:warning
            ...
            UserWarning: contains_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use the == operator instead
            doctest:warning
            ...
            UserWarning: Singularity() is deprecated and will be removed in a future version of sage-flatsurf. Use surface.point() instead.
            True
            sage: p.contains_vertex(label=1, v=0)
            True

        """
        import warnings

        warnings.warn(
            "contains_vertex() is deprecated and will be removed in a future version of sage-flatsurf; use the == operator instead"
        )

        if v is None:
            label, v = label

        return Singularity(self.surface(), label, v) == self

    def num_coordinates(self):
        r"""
        Return the number of different coordinate representations of the point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces

            sage: permutation = SymmetricGroup(2)('(1, 2)')
            sage: S = translation_surfaces.origami(permutation, permutation)
            sage: p = S.point(1, (0, 0))
            sage: p.num_coordinates()
            doctest:warning
            ...
            UserWarning: num_coordinates() is deprecated and will be removed in a future version of sage-flatsurf; use len(representatives()) instead.
            4

        """
        import warnings

        warnings.warn(
            "num_coordinates() is deprecated and will be removed in a future version of sage-flatsurf; use len(representatives()) instead."
        )

        return len(self._representatives)

    def labels(self):
        r"""
        Return the labels of polygons containing the point.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

        For a point in the interior of polygon, there is exactly one label::

            sage: p = S.point(0, (1/2, 1/2))
            sage: p.labels()
            {0}

        For a point in the interior of an edge of a polygon, there can be up to
        two labels::

            sage: p = S.point(0, (0, 1/2))
            sage: p.labels()
            {0, 2}

        For a point at a vertex, there can be more labels::

            sage: p = S.point(0, (0, 0))
            sage: p.labels()
            {0, 1, 2}

        """
        return {label for (label, _) in self._representatives}

    def coordinates(self, label):
        r"""
        Return coordinates for the point in the in the polygon ``label``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.mcmullen_L(1, 1, 1, 1)

            sage: p = S.point(0, (0, 0))
            sage: p.coordinates(0)  # random output: order depends on the Python version
            ((0, 0), (1, 0), (0, 1), (1, 1))

        """
        return tuple(
            coordinates for (l, coordinates) in self._representatives if l == label
        )

    def graphical_surface_point(self, graphical_surface=None):
        r"""
        Return a
        :class:`flatsurf.graphical.surface_point.GraphicalSurfacePoint` to
        represent this point graphically.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: G = p.graphical_surface_point()

        """
        from flatsurf.graphical.surface_point import GraphicalSurfacePoint

        return GraphicalSurfacePoint(self, graphical_surface=graphical_surface)

    def plot(self, *args, **kwargs):
        r"""
        Return a plot of this point.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (0, 0))
            sage: p.plot()
            ...Graphics object consisting of 1 graphics primitive

            sage: p = S.point(0, (0, 25/12))
            sage: p.plot()
            ...Graphics object consisting of 1 graphics primitive

        """
        graphical_surface = None
        if args:
            graphical_surface = args[0]
            args = args[1:]

        return self.graphical_surface_point(graphical_surface=graphical_surface).plot(
            *args, **kwargs
        )

    def __repr__(self):
        r"""
        Return a printable representation of this point.

        TESTS::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: S.point(0, (1/2, 1/2))
            Point (1/2, 1/2) of polygon 0

        """

        def render(label, coordinates):
            if self.is_vertex():
                vertex = (
                    self.surface()
                    .polygon(label)
                    .get_point_position(coordinates)
                    .get_vertex()
                )
                return "Vertex {} of polygon {}".format(vertex, label)

            return "Point {} of polygon {}".format(coordinates, label)

        # We pick a specific representative to make our lives easier when doctesting
        return min(
            render(label, coordinates)
            for (label, coordinates) in self.representatives()
        )

    def __eq__(self, other):
        r"""
        Return whether this point is indistinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p == p
            True

            sage: q = S.point(0, (1/2, 1/3))
            sage: p == q
            False

        TESTS:

        Verify that points can be compared to non-points so they can be put into sets and dicts with other objects::

            sage: p == 42
            False

        """
        if self is other:
            return True
        if not isinstance(other, SurfacePoint):
            return False
        if self._surface != other._surface:
            return False
        return self._representatives == other._representatives

    def _test_category(self, **options):
        r"""
        Check that this point inherits from the element class of its surface's
        category.

        Overridden to disable these tests when this is a point of a mutable
        surface since the category might then change as the surface becomes
        immutable.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p._test_category()

        """
        if self.surface().is_mutable():
            return

        super()._test_category(**options)

    def __hash__(self):
        r"""
        Return a hash value of this point.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (0, 0))
            sage: q = S.point(4, (0, 0))
            sage: hash(p) == hash(q)
            True

        """
        return hash(self._representatives)

    def __ne__(self, other):
        r"""
        Return whether this point is distinguishable from ``other``.

        EXAMPLES::

            sage: from flatsurf import half_translation_surfaces
            sage: S = half_translation_surfaces.step_billiard([1, 1, 1, 1], [1, 1/2, 1/3, 1/4])
            sage: p = S.point(0, (1/2, 1/2))
            sage: p != p
            False

            sage: q = S.point(0, (1/2, 1/3))
            sage: p != q
            True

        """
        return not (self == other)


class Cylinder(SageObject):
    r"""
    Represents a cylinder in a SimilaritySurface. A cylinder for these purposes is a
    topological annulus in a surface bounded by a finite collection of saddle connections
    meeting at 180 degree angles.

    To Do
    -----
    * Support cylinders whose monodromy is a dilation.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s = translation_surfaces.octagon_and_squares()
        sage: from flatsurf.geometry.surface_objects import Cylinder
        sage: cyl = Cylinder(s, 0, [2, 3, 3, 3, 2, 0, 1, 3, 2, 0])
        sage: cyl.initial_label()
        0
        sage: cyl.edges()
        (2, 3, 3, 3, 2, 0, 1, 3, 2, 0)
        sage: # a = sqrt(2) below.
        sage: cyl.area()
        2*a + 4
        sage: cyl.circumference().minpoly()
        x^4 - 680*x^2 + 400
        sage: cyl.holonomy()
        (8*a + 12, 4*a + 6)
    """

    def __init__(self, s, label0, edges):
        r"""
        Construct a cylinder on the surface `s` from an initial label and a
        sequence of edges crossed.

        Parameters
        ----------
        s: A SimilaritySurface
            the surface containing the cylinder
        label0: An initial label
            representing a polygon the cylinder passes through.
        edges: a list
            giving the sequence of edges the cylinder crosses until it closes.
        """
        self._surface = s
        self._label0 = label0
        self._edges = tuple(edges)
        ss = s.minimal_cover(cover_type="planar")
        SG = SimilarityGroup(s.base_ring())
        labels = [(label0, SG.one())]  # labels of polygons on the cover ss.
        for e in edges:
            labels.append(ss.opposite_edge(labels[-1], e)[0])
        if labels[0][0] != labels[-1][0]:
            raise ValueError("Combinatorial path does not close.")
        trans = labels[-1][1]
        if not trans.is_translation():
            raise NotImplementedError(
                "Only cylinders with translational monodromy are currently supported"
            )
        m = trans.matrix()
        v = vector(s.base_ring(), (m[0][2], m[1][2]))  # translation vector
        from flatsurf.geometry.euclidean import ccw

        p = ss.polygon(labels[0])
        e = edges[0]
        min_y = ccw(v, p.vertex(e))
        max_y = ccw(v, p.vertex((e + 1) % len(p.vertices())))
        if min_y >= max_y:
            raise ValueError("Combinatorial data does not represent a cylinder")

        # Stores the vertices where saddle connections starts:
        min_list = [0]
        max_list = [0]

        for i in range(1, len(edges)):
            e = edges[i]
            p = ss.polygon(labels[i])
            y = ccw(v, p.vertex(e))
            if y == min_y:
                min_list.append(i)
            elif y > min_y:
                min_list = [i]
                min_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")
            y = ccw(v, p.vertex((e + 1) % len(p.vertices())))
            if y == max_y:
                max_list.append(i)
            elif y < max_y:
                max_list = [i]
                max_y = y
                if min_y >= max_y:
                    raise ValueError("Combinatorial data does not represent a cylinder")

        # Extract the saddle connections on the right side:
        from flatsurf.geometry.saddle_connection import SaddleConnection

        sc_set_right = set()
        vertices = []
        for i in min_list:
            label = labels[i]
            p = ss.polygon(label)
            vertices.append((i, p.vertex(edges[i])))
        i, vert_i = vertices[-1]
        vert_i = vert_i - v
        j, vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v) * li[1])
            lio = ss.opposite_edge(li, edges[i])
            lj = labels[j]
            sc = SaddleConnection(
                s,
                (lio[0][0], (lio[1] + 1) % len(ss.polygon(lio[0]).vertices())),
                (~lio[0][1])(vert_j) - (~lio[0][1])(vert_i),
            )
            sc_set_right.add(sc)
        i = j
        vert_i = vert_j
        for j, vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                li = (li[0], SG(-v) * li[1])
                lio = ss.opposite_edge(li, edges[i])
                lj = labels[j]
                sc = SaddleConnection(
                    s,
                    (lio[0][0], (lio[1] + 1) % len(ss.polygon(lio[0]).vertices())),
                    (~lio[0][1])(vert_j) - (~lio[0][1])(vert_i),
                    limit=j - i,
                )
                sc_set_right.add(sc)
            i = j
            vert_i = vert_j

        # Extract the saddle connections on the left side:
        sc_set_left = set()
        vertices = []
        for i in max_list:
            label = labels[i]
            p = ss.polygon(label)
            vertices.append((i, p.vertex((edges[i] + 1) % len(p.vertices()))))
        i, vert_i = vertices[-1]
        vert_i = vert_i - v
        j, vert_j = vertices[0]
        if vert_i != vert_j:
            li = labels[i]
            li = (li[0], SG(-v) * li[1])
            lio = ss.opposite_edge(li, edges[i])
            lj = labels[j]
            sc = SaddleConnection(
                s,
                (lj[0], (edges[j] + 1) % len(ss.polygon(lj).vertices())),
                (~lj[1])(vert_i) - (~lj[1])(vert_j),
            )
            sc_set_left.add(sc)
        i = j
        vert_i = vert_j
        for j, vert_j in vertices[1:]:
            if vert_i != vert_j:
                li = labels[i]
                lio = ss.opposite_edge(li, edges[i])
                lj = labels[j]
                sc = SaddleConnection(
                    s,
                    (lj[0], (edges[j] + 1) % len(ss.polygon(lj).vertices())),
                    (~lj[1])(vert_i) - (~lj[1])(vert_j),
                )
                sc_set_left.add(sc)
            i = j
            vert_i = vert_j
        self._boundary1 = frozenset(sc_set_right)
        self._boundary2 = frozenset(sc_set_left)
        self._boundary = frozenset(self._boundary1.union(self._boundary2))

        edge_intersections = []
        i = min_list[0]
        label = labels[i]
        p = ss.polygon(label)
        right_point = p.vertex(edges[i])  # point on the right boundary
        i = max_list[0]
        label = labels[i]
        p = ss.polygon(label)
        left_point = p.vertex((edges[i] + 1) % len(p.vertices()))
        from flatsurf.geometry.euclidean import solve

        for i in range(len(edges)):
            label = labels[i]
            p = ss.polygon(label)
            e = edges[i]
            v1 = p.vertex(e)
            v2 = p.vertex((e + 1) % len(p.vertices()))
            a, b = solve(left_point, v, v1, v2 - v1)
            w1 = (~(label[1]))(v1 + b * (v2 - v1))
            a, b = solve(right_point, v, v1, v2 - v1)
            w2 = (~(label[1]))(v1 + b * (v2 - v1))
            edge_intersections.append((w1, w2))

        polygons = []
        pair1 = edge_intersections[-1]
        l1 = labels[-2][0]
        e1 = edges[-1]
        for i in range(len(edges)):
            l2 = labels[i][0]
            pair2 = edge_intersections[i]
            e2 = edges[i]
            trans = s.edge_transformation(l1, e1)
            pair1p = (trans(pair1[0]), trans(pair1[1]))
            polygon_verts = [pair1p[0], pair1p[1]]
            if pair2[1] != pair1p[1]:
                polygon_verts.append(pair2[1])
            if pair2[0] != pair1p[0]:
                polygon_verts.append(pair2[0])

            from flatsurf import Polygon

            polygons.append(
                (l2, Polygon(vertices=polygon_verts, base_ring=s.base_ring()))
            )
            l1 = l2
            pair1 = pair2
            e1 = e2
        self._polygons = tuple(polygons)

    def surface(self):
        return self._surface

    def initial_label(self):
        r"""
        Return one label on the surface that the cylinder passes through.
        """
        return self._label0

    def edges(self):
        r"""
        Return a tuple of edge numbers representing the edges crossed
        when the cylinder leaves the polygon with `initial_label` until
        it returns by closing.
        """
        return self._edges

    def boundary(self):
        r"""
        Return the set of saddle connections in the boundary, oriented so that
        the surface is on the left.
        """
        return self._boundary

    def polygons(self):
        r"""
        Return a list of pairs each consisting of a label and a polygon.
        Each polygon represents a sub-polygon of the polygon on the surface
        with the given label. The union of these sub-polygons form the
        cylinder. The subpolygons are listed in cyclic order.
        """
        return self._polygons

    @cached_method
    def area(self):
        r"""
        Return the area of this cylinder if it is contained in a ConeSurface.
        """
        from flatsurf.geometry.categories import ConeSurfaces

        if self._surface not in ConeSurfaces():
            raise NotImplementedError("area only makes sense for cone surfaces")

        area = 0
        for label, p in self.polygons():
            area += p.area()
        return area

    def plot(self, **options):
        r"""
        Plot this cylinder in coordinates used by a graphical surface. This
        plots this cylinder as a union of subpolygons. Only the intersections
        with polygons visible in the graphical surface are shown.

        Parameters other than `graphical_surface` are passed to `polygon2d`
        which is called to render the polygons.

        Parameters
        ----------
        graphical_surface : a GraphicalSurface
            If not provided or `None`, the plot method uses the default graphical
            surface for the surface.
        """
        if "graphical_surface" in options and options["graphical_surface"] is not None:
            gs = options["graphical_surface"]
            if gs.get_surface() != self._surface:
                raise ValueError("graphical surface for the wrong surface")
            del options["graphical_surface"]
        else:
            gs = self._surface.graphical_surface()
        plt = Graphics()
        for label, p in self.polygons():
            if gs.is_visible(label):
                gp = gs.graphical_polygon(label)
                t = gp.transformation()
                pp = t(p)
                poly = polygon2d(pp.vertices(), **options)
                plt += poly.plot()
        return plt

    @cached_method
    def labels(self):
        r"""
        Return the set of labels that this cylinder passes through.
        """
        polygons = self.polygons()
        return frozenset([label for label, p in polygons])

    def boundary_components(self):
        r"""
        Return a set of two elements: the set of saddle connections on
        the right and left sides. Saddle connections are oriented so that
        the surface is on the left.
        """
        return frozenset([self._boundary1, self._boundary2])

    def next(self, sc):
        r"""
        Return the next saddle connection as you move around the cylinder boundary
        moving from sc in the direction of its orientation.
        """
        if sc not in self._boundary:
            raise ValueError

        v = sc.end_tangent_vector()
        v = v.clockwise_to(-v.vector())
        from flatsurf.geometry.euclidean import is_parallel

        for sc2 in self._boundary:
            if sc2.start_data() == (
                v.polygon_label(),
                v.vertex(),
            ) and is_parallel(sc2.direction(), v.vector()):
                return sc2
        raise ValueError("Failed to find next saddle connection in boundary set.")

    def previous(self, sc):
        r"""
        Return the previous saddle connection as you move around the cylinder boundary
        moving from sc in the direction opposite its orientation.
        """
        if sc not in self._boundary:
            raise ValueError
        v = sc.start_tangent_vector()
        v = v.counterclockwise_to(-v.vector())
        from flatsurf.geometry.euclidean import is_parallel

        for sc2 in self._boundary:
            if sc2.end_data() == (v.polygon_label(), v.vertex()) and is_parallel(
                sc2.end_direction(), v.vector()
            ):
                return sc2
        raise ValueError("Failed to find previous saddle connection in boundary set.")

    @cached_method
    def holonomy(self):
        r"""
        In a translation surface, return one of the two holonomy vectors of the cylinder,
        which differ by a sign.
        """
        from flatsurf.geometry.categories import TranslationSurfaces

        if self._surface not in TranslationSurfaces():
            raise NotImplementedError(
                "holonomy currently only computable for translation surfaces"
            )

        V = self._surface.base_ring() ** 2
        total = V.zero()
        for sc in self._boundary1:
            total += sc.holonomy()

        # Debugging:
        total2 = V.zero()
        for sc in self._boundary2:
            total2 += sc.holonomy()
        assert (
            total + total2 == V.zero()
        ), "Holonomy of the two boundary components should sum to zero."

        return total

    @cached_method
    def circumference(self):
        r"""
        In a cone surface, return the circumference, i.e., the length
        of a geodesic loop running around the cylinder. Since this may
        not lie in the field of definition of the surface, it is returned
        as an element of the Algebraic Real Field.
        """
        from flatsurf.geometry.categories import ConeSurfaces

        if self._surface not in ConeSurfaces():
            raise NotImplementedError(
                "circumference only makes sense for cone surfaces"
            )

        total = 0
        for sc in self._boundary1:
            total += sc.length()
        return total
