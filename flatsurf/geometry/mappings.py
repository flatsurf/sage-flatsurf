r"""Mappings between translation surfaces."""

# *********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2022 W. Patrick Hooper
#                      2016-2022 Vincent Delecroix
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
# *********************************************************************


class SurfaceMapping:
    r"""Abstract class for any mapping between surfaces."""

    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    def domain(self):
        r"""
        Return the domain of the mapping.
        """
        return self._domain

    def codomain(self):
        r"""
        Return the range of the mapping.
        """
        return self._codomain

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        raise NotImplementedError

    def pull_vector_back(self, tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        raise NotImplementedError

    def __mul__(self, other):
        # Compose SurfaceMappings
        return SurfaceMappingComposition(other, self)

    def __rmul__(self, other):
        return SurfaceMappingComposition(self, other)


class SurfaceMappingComposition(SurfaceMapping):
    r"""
    Composition of two mappings between surfaces.
    """

    def __init__(self, mapping1, mapping2):
        r"""
        Represent the mapping of mapping1 followed by mapping2.
        """
        if mapping1.codomain() != mapping2.domain():
            raise ValueError(
                "Codomain of mapping1 must be equal to the domain of mapping2"
            )
        self._m1 = mapping1
        self._m2 = mapping2
        SurfaceMapping.__init__(self, self._m1.domain(), self._m2.codomain())

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        return self._m2.push_vector_forward(
            self._m1.push_vector_forward(tangent_vector)
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        return self._m1.pull_vector_back(self._m2.pull_vector_back(tangent_vector))

    def factors(self):
        r"""
        Return the two factors of this surface mapping as a pair (f,g),
        where the original map is f o g.
        """
        return self._m2, self._m1


class IdentityMapping(SurfaceMapping):
    r"""
    Construct an identity map between two "equal" surfaces.
    """

    def __init__(self, domain, codomain):
        SurfaceMapping.__init__(self, domain, codomain)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self._codomain.tangent_vector(
            tangent_vector.polygon_label(),
            tangent_vector.point(),
            tangent_vector.vector(),
            ring=ring,
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self._domain.tangent_vector(
            tangent_vector.polygon_label(),
            tangent_vector.point(),
            tangent_vector.vector(),
            ring=ring,
        )


class SimilarityJoinPolygonsMapping(SurfaceMapping):
    r"""
    Return a SurfaceMapping joining two polygons together along the edge provided to the constructor.

    EXAMPLES::

        sage: from flatsurf import MutableOrientedSimilaritySurface, Polygon
        sage: s = MutableOrientedSimilaritySurface(QQ)
        sage: s.add_polygon(Polygon(edges=[(1,0),(0,1),(-1,-1)]))
        0
        sage: s.add_polygon(Polygon(edges=[(-1,0),(0,-1),(1,1)]))
        1
        sage: s.glue((0, 0), (1, 0))
        sage: s.glue((0, 1), (1, 1))
        sage: s.glue((0, 2), (1, 2))
        sage: s.set_immutable()

        sage: from flatsurf.geometry.mappings import SimilarityJoinPolygonsMapping
        sage: m=SimilarityJoinPolygonsMapping(s, 0, 2)
        sage: s2=m.codomain()
        sage: s2.labels()
        (0,)
        sage: s2.polygons()
        (Polygon(vertices=[(0, 0), (1, 0), (1, 1), (0, 1)]),)
        sage: s2.gluings()
        (((0, 0), (0, 2)), ((0, 1), (0, 3)), ((0, 2), (0, 0)), ((0, 3), (0, 1)))

    """

    def __init__(self, s, p1, e1):
        r"""
        Join polygon with label p1 of s to polygon sharing edge e1.
        """
        if s.is_mutable():
            raise ValueError(
                "Can only construct SimilarityJoinPolygonsMapping for immutable surfaces."
            )

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        ss2 = MutableOrientedSimilaritySurface.from_surface(s)
        s2 = ss2

        poly1 = s.polygon(p1)
        p2, e2 = s.opposite_edge(p1, e1)
        poly2 = s.polygon(p2)
        t = s.edge_transformation(p2, e2)
        dt = t.derivative()
        vs = []  # actually stores new edges...
        edge_map = {}  # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(vs)] = (p1, i)
            vs.append(poly1.edge(i))
        ne = len(poly2.vertices())
        for i in range(1, ne):
            ee = (e2 + i) % ne
            edge_map[len(vs)] = (p2, ee)
            vs.append(dt * poly2.edge(ee))
        for i in range(e1 + 1, len(poly1.vertices())):
            edge_map[len(vs)] = (p1, i)
            vs.append(poly1.edge(i))

        inv_edge_map = {}
        for key, value in edge_map.items():
            inv_edge_map[value] = (p1, key)

        if p2 in s.roots():
            # The polygon with the base label is being removed.
            s2.set_roots(tuple(p1 if label == p2 else label for label in s.roots()))

        s2.remove_polygon(p1)
        from flatsurf import Polygon

        s2.add_polygon(Polygon(edges=vs, base_ring=s.base_ring()), label=p1)

        for i in range(len(vs)):
            p3, e3 = edge_map[i]
            p4, e4 = s.opposite_edge(p3, e3)
            if p4 == p1 or p4 == p2:
                pp, ee = inv_edge_map[(p4, e4)]
                s2.glue((p1, i), (pp, ee))
            else:
                s2.glue((p1, i), (p4, e4))

        s2.remove_polygon(p2)
        s2.set_immutable()

        self._saved_label = p1
        self._removed_label = p2
        self._remove_map = t
        self._remove_map_derivative = dt
        self._glued_edge = e1
        SurfaceMapping.__init__(self, s, ss2)

    def removed_label(self):
        r"""
        Return the label that was removed in the joining process.
        """
        return self._removed_label

    def glued_vertices(self):
        r"""
        Return the vertices of the newly glued polygon which bound the diagonal formed by the glue.
        """
        return (
            self._glued_edge,
            self._glued_edge
            + len(self._domain.polygon(self._removed_label).vertices()),
        )

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._removed_label:
            return self._codomain.tangent_vector(
                self._saved_label,
                self._remove_map(tangent_vector.point()),
                self._remove_map_derivative * tangent_vector.vector(),
                ring=ring,
            )
        else:
            return self._codomain.tangent_vector(
                tangent_vector.polygon_label(),
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring,
            )

    def pull_vector_back(self, tangent_vector):
        r"""
        Applies the inverse of the mapping to the provided vector.
        """
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._saved_label:
            p = tangent_vector.point()
            v = self._domain.polygon(self._saved_label).vertex(self._glued_edge)
            e = self._domain.polygon(self._saved_label).edge(self._glued_edge)
            from flatsurf.geometry.euclidean import ccw

            wp = ccw(p - v, e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector(
                    self._removed_label,
                    (~self._remove_map)(tangent_vector.point()),
                    (~self._remove_map_derivative) * tangent_vector.vector(),
                    ring=ring,
                )
            if wp < 0:
                # in polygon with the removed label
                return self.domain().tangent_vector(
                    self._saved_label,
                    tangent_vector.point(),
                    tangent_vector.vector(),
                    ring=ring,
                )
            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = ccw(w, e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector(
                    self._removed_label,
                    (~self._remove_map)(tangent_vector.point()),
                    (~self._remove_map_derivative) * tangent_vector.vector(),
                    ring=ring,
                )
            return self.domain().tangent_vector(
                self._saved_label,
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring,
            )
        else:
            return self._domain.tangent_vector(
                tangent_vector.polygon_label(),
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring,
            )


class SplitPolygonsMapping(SurfaceMapping):
    r"""
    Class for cutting a polygon along a diagonal.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.mappings import SplitPolygonsMapping
        sage: m = SplitPolygonsMapping(s,0,0,2)
        sage: s2=m.codomain()
        sage: TestSuite(s2).run()
        sage: s2.labels()
        (0, 1)
        sage: s2.polygons()
        (Polygon(vertices=[(0, 0), (1/2*a + 1, 1/2*a), (1/2*a + 1, 1/2*a + 1), (1, a + 1), (0, a + 1), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)]), Polygon(vertices=[(0, 0), (-1/2*a - 1, -1/2*a), (-1/2*a, -1/2*a)]))
        sage: s2.gluings()
        (((0, 0), (1, 0)), ((0, 1), (0, 5)), ((0, 2), (0, 6)), ((0, 3), (1, 1)), ((0, 4), (1, 2)), ((0, 5), (0, 1)), ((0, 6), (0, 2)), ((1, 0), (0, 0)), ((1, 1), (0, 3)), ((1, 2), (0, 4)))

    """

    def __init__(self, s, p, v1, v2, new_label=None):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.

        Warning: We do not ensure that new_label is not already in the list of labels unless it is None (as by default).
        """
        if s.is_mutable():
            raise ValueError("The surface should be immutable.")

        poly = s.polygon(p)
        ne = len(poly.vertices())
        if v1 < 0 or v2 < 0 or v1 >= ne or v2 >= ne:
            raise ValueError("Provided vertices out of bounds.")
        if abs(v1 - v2) <= 1 or abs(v1 - v2) >= ne - 1:
            raise ValueError("Provided diagonal is not a diagonal.")
        if v2 < v1:
            temp = v1
            v1 = v2
            v2 = temp

        newedges1 = [poly.vertex(v2) - poly.vertex(v1)]
        for i in range(v2, v1 + ne):
            newedges1.append(poly.edge(i))

        from flatsurf import Polygon

        newpoly1 = Polygon(edges=newedges1, base_ring=s.base_ring())

        newedges2 = [poly.vertex(v1) - poly.vertex(v2)]
        for i in range(v1, v2):
            newedges2.append(poly.edge(i))
        newpoly2 = Polygon(edges=newedges2, base_ring=s.base_ring())

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        ss2 = MutableOrientedSimilaritySurface.from_surface(s)
        s2 = ss2
        s2.remove_polygon(p)
        s2.add_polygon(newpoly1, label=p)
        new_label = s2.add_polygon(newpoly2, label=new_label)

        old_to_new_labels = {}
        for i in range(ne):
            if i < v1:
                old_to_new_labels[i] = (p, i + ne - v2 + 1)
            elif i < v2:
                old_to_new_labels[i] = (new_label, i - v1 + 1)
            else:  # i>=v2
                old_to_new_labels[i] = (p, i - v2 + 1)
        new_to_old_labels = {}
        for i, pair in old_to_new_labels.items():
            new_to_old_labels[pair] = i

        # This glues the split polygons together.
        s2.glue((p, 0), (new_label, 0))
        for e in range(ne):
            ll, ee = old_to_new_labels[e]
            lll, eee = s.opposite_edge(p, e)
            if lll == p:
                gl, ge = old_to_new_labels[eee]
                s2.glue((ll, ee), (gl, ge))
            else:
                s2.glue((ll, ee), (lll, eee))

        s2.set_immutable()

        self._p = p
        self._v1 = v1
        self._v2 = v2
        self._new_label = new_label
        from flatsurf.geometry.similarity import SimilarityGroup

        TG = SimilarityGroup(s.base_ring())
        self._tp = TG(-s.polygon(p).vertex(v1))
        self._tnew_label = TG(-s.polygon(p).vertex(v2))
        SurfaceMapping.__init__(self, s, ss2)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            point = tangent_vector.point()
            vertex1 = self._domain.polygon(self._p).vertex(self._v1)
            vertex2 = self._domain.polygon(self._p).vertex(self._v2)

            from flatsurf.geometry.euclidean import ccw

            wp = ccw(vertex2 - vertex1, point - vertex1)

            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector(
                    self._p,
                    self._tp(tangent_vector.point()),
                    tangent_vector.vector(),
                    ring=ring,
                )
            if wp < 0:
                # in new polygon 2
                return self.codomain().tangent_vector(
                    self._new_label,
                    self._tnew_label(tangent_vector.point()),
                    tangent_vector.vector(),
                    ring=ring,
                )

            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = ccw(vertex2 - vertex1, w)
            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector(
                    self._p,
                    self._tp(tangent_vector.point()),
                    tangent_vector.vector(),
                    ring=ring,
                )
            # in new polygon 2
            return self.codomain().tangent_vector(
                self._new_label,
                self._tnew_label(tangent_vector.point()),
                tangent_vector.vector(),
                ring=ring,
            )
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._codomain.tangent_vector(
                tangent_vector.polygon_label(),
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring,
            )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            return self._domain.tangent_vector(
                self._p,
                (~self._tp)(tangent_vector.point()),
                tangent_vector.vector(),
                ring=ring,
            )
        elif tangent_vector.polygon_label() == self._new_label:
            return self._domain.tangent_vector(
                self._p,
                (~self._tnew_label)(tangent_vector.point()),
                tangent_vector.vector(),
                ring=ring,
            )
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._domain.tangent_vector(
                tangent_vector.polygon_label(),
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring,
            )


def subdivide_a_polygon(s):
    r"""
    Return a SurfaceMapping which cuts one polygon along a diagonal or None if the surface is triangulated.
    """
    from flatsurf.geometry.euclidean import ccw

    for label, poly in zip(s.labels(), s.polygons()):
        n = len(poly.vertices())
        if n > 3:
            for i in range(n):
                e1 = poly.edge(i)
                e2 = poly.edge((i + 1) % n)
                if ccw(e1, e2) != 0:
                    return SplitPolygonsMapping(s, label, i, (i + 2) % n)
            raise ValueError(
                "Unable to triangulate polygon with label "
                + str(label)
                + ": "
                + str(poly)
            )
    return None


def triangulation_mapping(s):
    r"""
    Return a SurfaceMapping triangulating ``s``.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.mappings import triangulation_mapping
        sage: m=triangulation_mapping(s)
        sage: s2=m.codomain()
        sage: TestSuite(s2).run()
        sage: s2.polygons()
        (Polygon(vertices=[(0, 0), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)]),
         Polygon(vertices=[(0, 0), (1/2*a, -1/2*a - 1), (1/2*a, 1/2*a)]),
         Polygon(vertices=[(0, 0), (-1/2*a - 1, -1/2*a - 1), (0, -1)]),
         Polygon(vertices=[(0, 0), (-1, -a - 1), (1/2*a, -1/2*a)]),
         Polygon(vertices=[(0, 0), (0, -a - 1), (1, 0)]),
         Polygon(vertices=[(0, 0), (-1/2*a - 1, -1/2*a), (-1/2*a, -1/2*a)]))

    """
    if not s.is_finite_type():
        raise NotImplementedError

    m = subdivide_a_polygon(s)
    if m is None:
        return None
    s1 = m.codomain()
    while True:
        m2 = subdivide_a_polygon(s1)
        if m2 is None:
            return m
        s1 = m2.codomain()
        m = SurfaceMappingComposition(m, m2)
    return m


def flip_edge_mapping(s, p1, e1):
    r"""
    Return a mapping whose domain is s which flips the provided edge.
    """
    m1 = SimilarityJoinPolygonsMapping(s, p1, e1)
    v1, v2 = m1.glued_vertices()
    removed_label = m1.removed_label()
    m2 = SplitPolygonsMapping(
        m1.codomain(), p1, (v1 + 1) % 4, (v1 + 3) % 4, new_label=removed_label
    )
    return SurfaceMappingComposition(m1, m2)


def one_delaunay_flip_mapping(s):
    r"""
    Returns one delaunay flip, or none if no flips are needed.
    """
    for p, poly in zip(s.labels(), s.polygons()):
        for e in range(len(poly.vertices())):
            if s._delaunay_edge_needs_flip(p, e):
                return flip_edge_mapping(s, p, e)
    return None


def delaunay_triangulation_mapping(s):
    r"""
    Returns a mapping to a Delaunay triangulation or None if the surface already is Delaunay triangulated.
    """
    if not s.is_finite_type():
        raise NotImplementedError

    m = triangulation_mapping(s)
    if m is None:
        s1 = s
    else:
        s1 = m.codomain()
    m1 = one_delaunay_flip_mapping(s1)
    if m1 is None:
        return m
    if m is None:
        m = m1
    else:
        m = SurfaceMappingComposition(m, m1)
    s1 = m1.codomain()
    while True:
        m1 = one_delaunay_flip_mapping(s1)
        if m1 is None:
            return m
        s1 = m1.codomain()
        m = SurfaceMappingComposition(m, m1)


def delaunay_decomposition_mapping(s):
    r"""
    Returns a mapping to a Delaunay decomposition or possibly None if the surface already is Delaunay.
    """
    m = delaunay_triangulation_mapping(s)
    if m is None:
        s1 = s
    else:
        s1 = m.codomain()

    joins = set()
    edge_vectors = []

    for p, poly in zip(s1.labels(), s1.polygons()):
        for e in range(len(poly.vertices())):
            pp, ee = s1.opposite_edge(p, e)
            if (pp, ee) in joins:
                continue
            if s1._delaunay_edge_needs_join(p, e):
                joins.add((p, e))
                edge_vectors.append(s1.tangent_vector(p, poly.vertex(e), poly.edge(e)))

    if len(edge_vectors) > 0:
        ev = edge_vectors.pop()
        p, e = ev.edge_pointing_along()
        m1 = SimilarityJoinPolygonsMapping(s1, p, e)
        s2 = m1.codomain()
        while len(edge_vectors) > 0:
            ev = edge_vectors.pop()
            ev2 = m1.push_vector_forward(ev)
            p, e = ev2.edge_pointing_along()
            mtemp = SimilarityJoinPolygonsMapping(s2, p, e)
            m1 = SurfaceMappingComposition(m1, mtemp)
            s2 = m1.codomain()
        if m is None:
            return m1
        else:
            return SurfaceMappingComposition(m, m1)
    return m


def canonical_first_vertex(polygon):
    r"""
    Return the index of the vertex with smallest y-coordinate.
    If two vertices have the same y-coordinate, then the one with least x-coordinate is returned.
    """
    best = 0
    best_pt = polygon.vertex(best)
    for v in range(1, len(polygon.vertices())):
        pt = polygon.vertex(v)
        if pt[1] < best_pt[1]:
            best = v
            best_pt = pt
    if best == 0:
        if pt[1] == best_pt[1]:
            return v
    return best


class CanonicalizePolygonsMapping(SurfaceMapping):
    r"""
    This is a mapping to a surface with the polygon vertices canonically determined.
    A canonical labeling is when the canonocal_first_vertex is the zero vertex.
    """

    def __init__(self, s):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.
        """
        if not s.is_finite_type():
            raise ValueError("Currently only works with finite surfaces.")
        ring = s.base_ring()
        from flatsurf.geometry.similarity import SimilarityGroup

        T = SimilarityGroup(ring)
        cv = {}  # dictionary for canonical vertices
        translations = {}  # translations bringing the canonical vertex to the origin.
        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        s2 = MutableOrientedSimilaritySurface(ring)
        for label, polygon in zip(s.labels(), s.polygons()):
            cv[label] = cvcur = canonical_first_vertex(polygon)
            newedges = []
            for i in range(len(polygon.vertices())):
                newedges.append(polygon.edge((i + cvcur) % len(polygon.vertices())))

            from flatsurf import Polygon

            s2.add_polygon(Polygon(edges=newedges, base_ring=ring), label=label)
            translations[label] = T(-polygon.vertex(cvcur))
        for l1, polygon in zip(s.labels(), s.polygons()):
            for e1 in range(len(polygon.vertices())):
                l2, e2 = s.opposite_edge(l1, e1)
                ee1 = (e1 - cv[l1] + len(polygon.vertices())) % len(polygon.vertices())
                polygon2 = s.polygon(l2)
                ee2 = (e2 - cv[l2] + len(polygon2.vertices())) % len(
                    polygon2.vertices()
                )
                # newgluing.append( ( (l1,ee1),(l2,ee2) ) )
                s2.glue((l1, ee1), (l2, ee2))
        s2.set_roots(s.roots())
        s2.set_immutable()
        ss2 = s2

        self._cv = cv
        self._translations = translations

        SurfaceMapping.__init__(self, s, ss2)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        label = tangent_vector.polygon_label()
        return self.codomain().tangent_vector(
            label,
            self._translations[label](tangent_vector.point()),
            tangent_vector.vector(),
            ring=ring,
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        label = tangent_vector.polygon_label()
        return self.domain().tangent_vector(
            label,
            (~self._translations[label])(tangent_vector.point()),
            tangent_vector.vector(),
            ring=ring,
        )


class ReindexMapping(SurfaceMapping):
    r"""
    Apply a dictionary to relabel the polygons.
    """

    def __init__(self, s, relabler, new_base_label=None):
        r"""
        The parameters should be a surface and a dictionary which takes as input a label and produces a new label.
        """
        if not s.is_finite_type():
            raise ValueError("Currently only works with finite surfaces." "")
        f = {}  # map for labels going forward.
        b = {}  # map for labels going backward.
        for label in s.labels():
            if label in relabler:
                l2 = relabler[label]
                f[label] = l2
                if l2 in b:
                    raise ValueError(
                        "Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change."
                    )
                b[l2] = label
            else:
                # If no key then don't change the label
                f[label] = label
                if label in b:
                    raise ValueError(
                        "Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change."
                    )
                b[label] = label

        self._f = f
        self._b = b

        if new_base_label is None:
            if s.root() in f:
                new_base_label = f[s.root()]
            else:
                new_base_label = s.root()
        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        s2 = MutableOrientedSimilaritySurface.from_surface(s)
        s2.relabel(relabler, in_place=True)
        s2.set_roots([new_base_label])
        s2.set_immutable()

        SurfaceMapping.__init__(self, s, s2)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        # There is no change- we just move it to the new surface.
        ring = tangent_vector.bundle().base_ring()
        return self.codomain().tangent_vector(
            self._f[tangent_vector.polygon_label()],
            tangent_vector.point(),
            tangent_vector.vector(),
            ring=ring,
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self.domain().tangent_vector(
            self._b[tangent_vector.polygon_label()],
            tangent_vector.point(),
            tangent_vector.vector(),
            ring=ring,
        )


def my_sgn(val):
    if val > 0:
        return 1
    elif val < 0:
        return -1
    else:
        return 0


def polygon_compare(poly1, poly2):
    r"""
    Compare two polygons first by area, then by number of sides,
    then by lexicographical ordering on edge vectors."""
    # This should not be used is broken!!
    # from sage.functions.generalized import sgn
    res = my_sgn(-poly1.area() + poly2.area())
    if res != 0:
        return res
    res = my_sgn(len(poly1.vertices()) - len(poly2.vertices()))
    if res != 0:
        return res
    ne = len(poly1.vertices())
    for i in range(0, ne - 1):
        edge_diff = poly1.edge(i) - poly2.edge(i)
        res = my_sgn(edge_diff[0])
        if res != 0:
            return res
        res = my_sgn(edge_diff[1])
        if res != 0:
            return res
    return 0


def canonicalize_translation_surface_mapping(s):
    r"""
    Return the translation surface in a canonical form.

    EXAMPLES::

        sage: from flatsurf import translation_surfaces
        sage: s = translation_surfaces.octagon_and_squares().canonicalize()

        sage: TestSuite(s).run()

        sage: a = s.base_ring().gen()  # a is the square root of 2.

        sage: from flatsurf.geometry.mappings import GL2RMapping
        sage: from flatsurf.geometry.mappings import canonicalize_translation_surface_mapping
        sage: mat=Matrix([[1,2+a],[0,1]])
        sage: from flatsurf.geometry.mappings import GL2RMapping
        sage: m1=GL2RMapping(s, mat)
        sage: m2=canonicalize_translation_surface_mapping(m1.codomain())
        sage: m=m2*m1
        sage: m.domain().cmp(m.codomain())
        0
        sage: TestSuite(m.codomain()).run()
        sage: s=m.domain()
        sage: v=s.tangent_vector(0,(0,0),(1,1))
        sage: w=m.push_vector_forward(v)
        sage: print(w)
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (a + 3, 1)
    """
    from flatsurf.geometry.categories import TranslationSurfaces

    if not s.is_finite_type():
        raise NotImplementedError
    if s not in TranslationSurfaces():
        raise ValueError("Only defined for TranslationSurfaces")
    m1 = delaunay_decomposition_mapping(s)
    if m1 is None:
        s2 = s
    else:
        s2 = m1.codomain()
    m2 = CanonicalizePolygonsMapping(s2)
    if m1 is None:
        m = m2
    else:
        m = SurfaceMappingComposition(m1, m2)
    s2 = m.codomain()

    # This is essentially copy & paste from canonicalize() from TranslationSurfaces()
    from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

    s2copy = MutableOrientedSimilaritySurface.from_surface(s2)
    ss = MutableOrientedSimilaritySurface.from_surface(s2)
    labels = set(s2.labels())
    for label in labels:
        ss.set_roots([label])
        if ss.cmp(s2copy) > 0:
            s2copy.set_roots([label])

    s2copy.set_immutable()

    # We now have the base_label correct.
    # We will use the label walk to generate the canonical labeling of polygons.
    labels = {label: i for (i, label) in enumerate(s2copy.labels())}

    m3 = ReindexMapping(s2, labels, 0)
    return SurfaceMappingComposition(m, m3)


class GL2RMapping(SurfaceMapping):
    r"""
    This class pushes a surface forward under a matrix.

    Note that for matrices of negative determinant we need to relabel edges (because
    edges must have a counterclockwise cyclic order). For each n-gon in the surface,
    we relabel edges according to the involution `e \mapsto n-1-e`.

    EXAMPLE::

        sage: from flatsurf import translation_surfaces
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.mappings import GL2RMapping
        sage: mat=Matrix([[2,1],[1,1]])
        sage: m=GL2RMapping(s,mat)
        sage: TestSuite(m.codomain()).run()
    """

    def __init__(self, s, m, category=None):
        r"""
        Hit the surface s with the 2x2 matrix m which should have positive determinant.
        """
        from flatsurf.geometry.lazy import GL2RImageSurface

        codomain = GL2RImageSurface(s, m, category=category or s.category())
        self._m = m
        self._im = ~m
        SurfaceMapping.__init__(self, s, codomain)

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        return self.codomain().tangent_vector(
            tangent_vector.polygon_label(),
            self._m * tangent_vector.point(),
            self._m * tangent_vector.vector(),
        )

    def pull_vector_back(self, tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        return self.domain().tangent_vector(
            tangent_vector.polygon_label(),
            self._im * tangent_vector.point(),
            self._im * tangent_vector.vector(),
        )
