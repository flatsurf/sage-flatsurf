# -*- coding: utf-8 -*-
r"""
Similarity surfaces.
"""
#*********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2020 Vincent Delecroix
#                      2020      Julian Rüth
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
#*********************************************************************

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.sage_unittest import TestSuite

from sage.structure.sage_object import SageObject

from sage.rings.infinity import Infinity

from sage.rings.all import ZZ, QQ, AA, RIF, RR, NumberField

from sage.modules.free_module_element import vector

from sage.matrix.constructor import matrix, identity_matrix
from sage.modules.free_module import VectorSpace

from sage.all import FreeModule

from .matrix_2x2 import (is_similarity,
                    homothety_rotation_decomposition,
                    similarity_from_vectors,
                    rotation_matrix_angle,
                    is_cosine_sine_of_rational)

from .similarity import SimilarityGroup
from .polygon import ConvexPolygons, wedge_product, triangulate, build_faces

from .surface import Surface, Surface_dict, Surface_list, LabelComparator
from .surface_objects import Singularity, SaddleConnection, SurfacePoint
from .circle import Circle

ZZ_1 = ZZ.one()
ZZ_2 = ZZ_1 + ZZ_1


class SimilaritySurface(SageObject):
    r"""
    An oriented surface built from a set of polygons and edges identified with
    similarities (i.e. composition of homothety, rotations and translations).

    Each polygon is identified with a unique key (its label). The choice of the
    label of the polygons is done at startup. If the set is finite then by
    default the labels are the first non-negative integers 0,1,...

    The edges are labeled by a pair ``(polygon label, edge number)``.

    EXAMPLES:

    The easiest way to construct a similarity surface is to use the pre-built
    constructions from
    :class:`flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators`::

        sage: from flatsurf import polygons, similarity_surfaces
        sage: P = polygons(vertices=[(0,0), (2,0), (1,4), (0,5)])
        sage: similarity_surfaces.self_glued_polygon(P)
        HalfTranslationSurface built from 1 polygon

    The second way is to build a surface (using e.g. :class:`flatsurf.geometry.surface.Surface_list`)
    and then use this surface as an argument for class:`SimilaritySurface`)::

        sage: from flatsurf.geometry.similarity_surface import SimilaritySurface
        sage: from flatsurf.geometry.surface import Surface_list
        sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
        sage: Stop = Surface_list(QQ)
        sage: Stop.add_polygon(P)
        0
        sage: Stop.add_polygon(2*P)
        1
        sage: Stop.add_polygon(3*P)
        2
        sage: Stop.set_edge_pairing(0, 1, 1, 3)
        sage: Stop.set_edge_pairing(0, 0, 2, 2)
        sage: Stop.set_edge_pairing(0, 2, 2, 0)
        sage: Stop.set_edge_pairing(0, 3, 1, 1)
        sage: Stop.set_edge_pairing(1, 2, 2, 1)
        sage: Stop.set_edge_pairing(1, 0, 2, 3)
        sage: S = SimilaritySurface(Stop)
        sage: S
        SimilaritySurface built from 3 polygons

    To perform a sanity check on the obtained surface, you can run its test
    suite::

        sage: TestSuite(S).run()

    In the following example, we build two broken surfaces and
    check that the test suite fails as expected::

        sage: P = polygons(vertices=[(0,0), (1,0), (1,1), (0,1)])
        sage: Stop = Surface_list(QQ)
        sage: Stop.add_polygon(P)
        0
        sage: S = SimilaritySurface(Stop)
        sage: TestSuite(S).run()
        ...
          AssertionError: edge (0, 0) is not glued
          ------------------------------------------------------------
          The following tests failed: _test_gluings
        Failure in _test_underlying_surface
        The following tests failed: _test_underlying_surface

        sage: Stop.set_edge_pairing(0, 0, 0, 3)
        sage: Stop.set_edge_pairing(0, 1, 0, 3)
        sage: Stop.set_edge_pairing(0, 2, 0, 3)
        sage: S = SimilaritySurface(Stop)
        sage: TestSuite(S).run()
        ...
          AssertionError: edge gluing is not a pairing:
          (0, 0) -> (0, 3) -> (0, 2)
          ------------------------------------------------------------
          The following tests failed: _test_gluings
        Failure in _test_underlying_surface
        The following tests failed: _test_underlying_surface

    Finally, you can also implement a similarity surface by inheriting from
    :class:`SimilaritySurface` and implement the methods:

    - ``base_ring(self)``: the base ring in which coordinates lives

    - ``polygon(self, lab)``: the polygon associated to the label ``lab``

    - ``base_label(self)``: which label to use as the base one

    - ``opposite_edge(self, lab, edge)``: a pair (``other_label``,
      ``other_edge``) representing the edge being glued

    - ``is_finite(self)``: whether the surface is built from finitely many polygons
    """
    def __init__(self, surface):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity_surface import SimilaritySurface
            sage: SimilaritySurface(3)
            Traceback (most recent call last):
            ...
            TypeError: invalid argument surface=3 to build a half-translation surface
        """
        if isinstance(surface, SimilaritySurface):
            self._s = surface.underlying_surface()
        elif isinstance(surface, Surface):
            self._s = surface
        else:
            raise TypeError("invalid argument surface={} to build a half-translation surface".format(surface))

    def underlying_surface(self):
        r"""
        Return the surface underlying this SimilaritySurface.
        """
        return self._s

    def _test_underlying_surface(self, **options):
        is_sub_testsuite = 'tester' in options
        tester = self._tester(**options)
        tester.info("")

        # TODO: this nested subsuite is very fragile since we have no
        #       way of forwarding the doctests to be skipped... Since
        #       for now, the unique usage of this is for pickling of
        #       infinite surface we provide that manually
        if not self._s.is_finite():
            skip = ['_test_pickling']
        else:
            skip = []

        TestSuite(self._s).run(verbose = tester._verbose,
                                prefix = tester._prefix + "  ",
                                raise_on_failure=is_sub_testsuite,
                                skip=skip)
        tester.info(tester._prefix + " ", newline=False)

    def base_ring(self):
        r"""
        The field on which the coordinates of ``self`` live.

        This method must be overriden in subclasses!
        """
        return self._s.base_ring()

    def polygon(self, lab):
        r"""
        Return the polygon with label ``lab``.
        """
        return self._s.polygon(lab)

    def base_label(self):
        r"""
        Always returns the same label.
        """
        return self._s.base_label()

    def opposite_edge(self, l, e=None):
        r"""
        Given the label ``l`` of a polygon and an edge ``e`` in that polygon
        returns the pair (``ll``, ``ee``) to which this edge is glued.
        If e is not provided, then it expects the only parameter to be
        the pair (``l``,``e``) and will again return a the pair (``ll``,``ee``).
        """
        if e is None:
            return self._s.opposite_edge(l[0],l[1])
        return self._s.opposite_edge(l,e)

    def is_finite(self):
        r"""
        Return whether or not the surface is finite.
        """
        return self._s.is_finite()

    def is_mutable(self):
        r"""
        Return if the surface is mutable.
        """
        return self._s.is_mutable()

    def set_immutable(self):
        r"""
        Mark the surface as immutable.
        """
        self._s.set_immutable()

    def is_triangulated(self, limit=None):
        return self._s.is_triangulated(limit=limit)

    #
    # generic methods
    #

    #def compute_surface_type_from_gluings(self,limit=None):
    #    r"""
    #    Compute the surface type by looking at the edge gluings.
    #    If limit is defined, we try to guess the type by looking at limit many edges.
    #    """
    #    if limit is None:
    #        if not self.is_finite():
    #            raise ValueError("Need a limit when working with an infinite surface.")
    #        it = self.edge_iterator()
    #        label,edge = it.next()
    #        # Use honest matrices!
    #        m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #        surface_type = surface_type_from_matrix(m)
    #        for label,edge in it:
    #            # Use honest matrices!
    #            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #            surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
    #        return surface_type
    #    else:
    #        count=0
    #        it = self.edge_iterator()
    #        label,edge = it.next()
    #        # Use honest matrices!
    #        m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #        surface_type = surface_type_from_matrix(m)
    #        for label,edge in it:
    #            # Use honest matrices!
    #            m = SimilaritySurface_generic.edge_matrix(self,label,edge)
    #            surface_type = combine_surface_types(surface_type, surface_type_from_matrix(m))
    #            count=count+1
    #            if count >= limit:
    #                return surface_type
    #        return surface_type

    def walker(self):
        return self._s.walker()

    def label_iterator(self, polygons=False):
        r"""
        Iterator over all polygon labels.

        If the keyword polygons is True then we return pairs (label, polygon)
        instead of just labels.
        """
        if polygons:
            return self._s.label_polygon_iterator()
        else:
            return self._s.label_iterator()

    def edge_iterator(self, gluings=False):
        r"""
        Iterate over the edges of polygons, which are pairs (l,e) where l is a polygon label, 0 <= e < N and N is the number of edges of the polygon with label l.

        If the keyword gluings is set to true, then we iterate over ordered
        pairs of edges ((l,e),(ll,ee)) where edge (l,e) is glued to (ll,ee).

        EXAMPLES::

            sage: from flatsurf import ConvexPolygons
            sage: P = ConvexPolygons(QQ)
            sage: tri0=P([(1,0),(0,1),(-1,-1)])
            sage: tri1=P([(-1,0),(0,-1),(1,1)])
            sage: gluings=[((0,0),(1,0)),((0,1),(1,1)),((0,2),(1,2))]
            sage: from flatsurf.geometry.surface import surface_list_from_polygons_and_gluings
            sage: from flatsurf.geometry.translation_surface import TranslationSurface
            sage: s=TranslationSurface(surface_list_from_polygons_and_gluings([tri0,tri1], gluings))
            sage: for edge in s.edge_iterator():
            ....:     print(edge)
            (0, 0)
            (0, 1)
            (0, 2)
            (1, 0)
            (1, 1)
            (1, 2)
        """
        if gluings:
            return self._s.edge_gluing_iterator()
        else:
            return self._s.edge_iterator()

    def num_polygons(self):
        r"""
        Return the number of polygons.
        """
        return self._s.num_polygons()

    def num_edges(self):
        r"""
        Return the total number of edges of all polygons used.
        """
        return self._s.num_edges()

    def num_singularities(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *

            sage: translation_surfaces.regular_octagon().num_singularities()
            1

            sage: S = SymmetricGroup(4)
            sage: r = S('(1,2)(3,4)')
            sage: u = S('(2,3)')
            sage: translation_surfaces.origami(r,u).num_singularities()
            2

            sage: S = SymmetricGroup(8)
            sage: r = S('(1,2,3,4,5,6,7,8)')
            sage: u = S('(1,8,5,4)(2,3)(6,7)')
            sage: translation_surfaces.origami(r,u).num_singularities()
            4
        """
        if not self.is_finite():
            raise ValueError("the method only work for finite surfaces")

        # NOTE:
        # the very same code is implemented in the method angles (translation
        # surfaces). we should factor out the code
        edges = set((p,e) for p in self.label_iterator() for e in range(self.polygon(p).num_edges()))

        n = ZZ(0)
        while edges:
            p,e = edges.pop()
            n += 1
            ee = (e-1) % self.polygon(p).num_edges()
            p,e = self.opposite_edge(p,ee)
            while (p,e) in edges:
                edges.remove((p,e))
                ee = (e-1) % self.polygon(p).num_edges()
                p,e = self.opposite_edge(p,ee)
        return n

    def _repr_(self):
        if self.num_polygons() == Infinity:
            num = 'infinitely many'
        else:
            num = str(self.num_polygons())

        if self.num_polygons() == 1:
            end = ""
        else:
            end = "s"

        return "{} built from {} polygon{}".format(self.__class__.__name__, num, end)

    def edge_matrix(self, p, e=None):
        r"""
        Returns the 2x2 matrix representing a similarity which when applied to the polygon with label `p`
        makes it so the edge `e` can be glued to its opposite edge by translation.

        If `e` is not provided, then `p` should be a pair consisting of a polygon label and an edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: print(s.polygon(0))
            Polygon: (0, 0), (2, -2), (2, 0)
            sage: print(s.polygon(1))
            Polygon: (0, 0), (2, 0), (1, 3)
            sage: s.opposite_edge(0,0)
            (1, 1)
            sage: m = s.edge_matrix(0, 0)
            sage: m
            [   1  1/2]
            [-1/2    1]
            sage: m * vector((2,-2)) == -vector((-1, 3))
            True
        """
        if e is None:
            p,e = p
        u = self.polygon(p).edge(e)
        pp,ee = self.opposite_edge(p,e)
        v = self.polygon(pp).edge(ee)

        # be careful, because of the orientation, it is -v and not v
        return similarity_from_vectors(u,-v)

    def edge_transformation(self, p, e):
        r"""
        Return the similarity bringing the provided edge to the opposite edge.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s = SimilaritySurfaceGenerators.example()
            sage: print(s.polygon(0))
            Polygon: (0, 0), (2, -2), (2, 0)
            sage: print(s.polygon(1))
            Polygon: (0, 0), (2, 0), (1, 3)
            sage: print(s.opposite_edge(0,0))
            (1, 1)
            sage: g = s.edge_transformation(0,0)
            sage: g((0,0))
            (1, 3)
            sage: g((2,-2))
            (2, 0)
        """
        G=SimilarityGroup(self.base_ring())
        q=self.polygon(p)
        a=q.vertex(e)
        b=q.vertex(e+1)
        # This is the similarity carrying the origin to a and (1,0) to b:
        g=G(b[0]-a[0],b[1]-a[1],a[0],a[1])

        pp,ee = self.opposite_edge(p,e)
        qq=self.polygon(pp)
        # Be careful here: opposite vertices are identified
        aa=qq.vertex(ee+1)
        bb=qq.vertex(ee)
        # This is the similarity carrying the origin to aa and (1,0) to bb:
        gg=G(bb[0]-aa[0],bb[1]-aa[1],aa[0],aa[1])

        # This is the similarity carrying (a,b) to (aa,bb):
        return gg/g

    def set_vertex_zero(self, label, v, in_place=False):
        r"""
        Applies a combinatorial rotation to the polygon with the provided label.

        This makes what is currently vertex v of this polygon vertex 0. In other words,
        what is currently vertex (or edge) e will now become vertex (e-v)%n where
        n is the number of sides of the polygon.

        For the updated polygons, the polygons will be translated so that vertex
        0 is the origin.

        EXAMPLES:

        Example with polygon glued to another polygon::

            sage: from flatsurf import *
            sage: s = translation_surfaces.veech_double_n_gon(4)
            sage: s.polygon(0)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: [s.opposite_edge(0,i) for i in range(4)]
            [(1, 0), (1, 1), (1, 2), (1, 3)]
            sage: ss = s.set_vertex_zero(0,1)
            sage: ss.polygon(0)
            Polygon: (0, 0), (0, 1), (-1, 1), (-1, 0)
            sage: [ss.opposite_edge(0,i) for i in range(4)]
            [(1, 1), (1, 2), (1, 3), (1, 0)]
            sage: TestSuite(ss).run()

        Example with polygon glued to self::

            sage: s = translation_surfaces.veech_2n_gon(2)
            sage: s.polygon(0)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: [s.opposite_edge(0,i) for i in range(4)]
            [(0, 2), (0, 3), (0, 0), (0, 1)]
            sage: ss = s.set_vertex_zero(0,3)
            sage: ss.polygon(0)
            Polygon: (0, 0), (0, -1), (1, -1), (1, 0)
            sage: [ss.opposite_edge(0,i) for i in range(4)]
            [(0, 2), (0, 3), (0, 0), (0, 1)]
            sage: TestSuite(ss).run()
        """
        if in_place:
            us = self.underlying_surface()
            if not us.is_mutable():
                raise ValueError("set_vertex_zero can only be done in_place for a mutable surface.")
            p = us.polygon(label)
            n=p.num_edges()
            assert 0<=v and v<n
            glue=[]
            P = ConvexPolygons(us.base_ring())
            pp = P(edges=[p.edge((i+v)%n) for i in range(n)])

            for i in range(n):
                e=(v+i)%n
                ll,ee = us.opposite_edge(label,e)
                if ll==label:
                    ee = (ee+n-v)%n
                glue.append((ll,ee))

            us.change_polygon(label,pp,gluing_list=glue)
            return self
        else:
            return self.copy(mutable=True).set_vertex_zero(label,v,in_place=True)

    def _label_comparator(self):
        r"""
        Return a LabelComparator, which provides a fixed total ordering on the polygon labels.
        """
        try:
            return self._lc
        except AttributeError:
            self._lc = LabelComparator()
            return self._lc

    def relabel(self, relabeling_map, in_place=False):
        r"""
        Attempt to relabel the polygons according to a relabeling_map, which takes as input
        a current label and outputs a new label for the same polygon. The method returns a pair
        (surface,success) where surface is the relabeled surface, and success is a boolean value
        indicating the success of the operation. The operation will fail if the implementation of the
        underlying surface does not support labels used in the image of the relabeling map. In this case,
        other (arbitrary) labels will be used to replace the labels of the surface, and the resulting
        surface should still be okay.

        Currently, the relabeling_map must be a dictionary.

        If in_place is True then the relabeling is done to the current surface, otherwise a
        mutable copy is made before relabeling.

        ToDo:
          - Allow relabeling_map to be a function rather than just a dictionary.
            This will allow it to work for infinite surfaces.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.veech_double_n_gon(5)
            sage: ss,valid=s.relabel({0:1,1:2})
            sage: valid
            True
            sage: ss.base_label()
            1
            sage: ss.opposite_edge(1,0)
            (2, 0)
            sage: ss.num_polygons()
            2
            sage: TestSuite(ss).run()
        """
        if in_place:
            us = self.underlying_surface()
            if not us.is_mutable():
                raise ValueError("Your surface is not mutable, so can not be relabeled in place.")
            if not isinstance(relabeling_map,dict):
                raise NotImplementedError("Currently relabeling is only implemented via a dictionary.")
            domain=set()
            codomain=set()
            data={}
            for l1,l2 in iteritems(relabeling_map):
                p=us.polygon(l1)
                glue = []
                for e in range(p.num_edges()):
                    ll,ee = us.opposite_edge(l1,e)
                    try:
                        lll=relabeling_map[ll]
                    except KeyError:
                        lll=ll
                    glue.append((lll,ee))
                data[l2]=(p,glue)
                domain.add(l1)
                codomain.add(l2)
            if len(domain)!=len(codomain):
                raise ValueError("The relabeling_map must be injective. Received "+str(relabeling_map))
            changed_labels = domain.intersection(codomain)
            added_labels=codomain.difference(domain)
            removed_labels=domain.difference(codomain)
            # Pass to add_polygons
            relabel_errors={}
            for l2 in added_labels:
                p,glue=data[l2]
                l3 = us.add_polygon(p, label=l2)
                if not l2==l3:
                    # This means the label l2 could not be added for some reason.
                    # Perhaps the implementation does not support this type of label.
                    # Or perhaps there is already a polygon with this label.
                    relabel_errors[l2]=l3
            # Pass to change polygons
            for l2 in changed_labels:
                p,glue=data[l2]
                # This should always work since the domain of the relabeling map should be labels for polygons.
                us.change_polygon(l2,p)
            # Deal with the base_label
            base_label = us.base_label()
            if base_label in relabeling_map:
                base_label = relabeling_map[base_label]
                if base_label in relabel_errors:
                    base_label = relabel_errors[base_label]
                us.change_base_label(base_label)
            # Pass to remove polygons:
            for l1 in removed_labels:
                us.remove_polygon(l1)
            # Pass to update the edge gluings
            if len(relabel_errors)==0:
                # No problems. Update the gluings.
                for l2 in codomain:
                    p,glue=data[l2]
                    us.change_polygon_gluings(l2, glue)
            else:
                # Use the gluings provided by relabel_errors when necessary
                for l2 in codomain:
                    p,glue=data[l2]
                    for e in range(p.num_edges()):
                        ll,ee=glue[e]
                        try:
                            # First try the error dictionary
                            us.change_edge_gluing(l2, e, relabel_errors[ll],ee)
                        except KeyError:
                            us.change_edge_gluing(l2, e, ll,ee)
            return self, len(relabel_errors)==0
        else:
            return self.copy(mutable=True).relabel(relabeling_map, in_place=True)

    def copy(self, relabel=False, mutable=False, lazy=None, new_field=None, optimal_number_field=False):
        r"""
        Returns a copy of this surface. The method takes several flags to modify how the copy is taken.

        If relabel is True, then instead of returning an exact copy, it returns a copy indexed by the
        non-negative integers. This uses the Surface_list implementation. If relabel is False (default),
        then we return an exact copy. The returned surface uses the Surface_dict implementation.

        The mutability flag returns if the resulting surface should be mutable or not. By default, the
        resulting surface will not be mutable.

        If lazy is True, then the surface is copied by reference. This is the only type of copy
        possible for infinite surfaces. The parameter defaults to False for finite surfaces, and
        defaults to True for infinite surfaces.

        The new_field parameter can be used to place the vertices in a larger field than the basefield
        for the original surface.

        The optimal_number_field option can be used to find a best NumberField containing the
        (necessarily finite) surface.

        EXAMPLES::

            sage: from flatsurf import *
            sage: ss=translation_surfaces.ward(3)
            sage: print(ss.is_mutable())
            False
            sage: s=ss.copy(mutable=True)
            sage: print(s.is_mutable())
            True
            sage: TestSuite(s).run()
            sage: print(s==ss)
            True

            sage: # Changing the base field
            sage: from flatsurf import *
            sage: s=translation_surfaces.veech_double_n_gon(5)
            sage: ss=s.copy(mutable=False,new_field=AA)
            sage: TestSuite(ss).run()
            sage: ss.base_ring()
            Algebraic Real Field

            sage: # Optimization of number field
            sage: from flatsurf import *
            sage: s = translation_surfaces.arnoux_yoccoz(3)
            sage: ss = s.copy(new_field=AA).copy(optimal_number_field=True)
            sage: TestSuite(ss).run()
            sage: ss.base_ring().discriminant()
            -44
        """
        s = None  # This will be the surface we copy. (Likely we will set s=self below.)
        if new_field is not None and optimal_number_field:
            raise ValueError("You can not set a new_field and also set optimal_number_field=True.")
        if optimal_number_field == True:
            assert self.is_finite(), "Can only optimize_number_field for a finite surface."
            assert not lazy, "Lazy copying is unavailable when optimize_number_field=True."
            coordinates_AA = []
            for l,p in self.label_iterator(polygons = True):
                for e in p.edges():
                    coordinates_AA.append(AA(e[0]))
                    coordinates_AA.append(AA(e[1]))
            from sage.rings.qqbar import number_field_elements_from_algebraics
            field,coordinates_NF,hom = number_field_elements_from_algebraics(coordinates_AA, minimal = True)
            if field is QQ:
                new_field = QQ
                # We pretend new_field = QQ was passed as a parameter.
                # It will now get picked up by the "if new_field is not None:" line below.
            else:
                # Unfortunately field doesn't come with an real embedding (which is given by hom!)
                # So, we make a copy of the field, and add the embedding.
                field2 = NumberField(field.polynomial(), name = "a", embedding = hom(field.gen()))
                # The following converts from field to field2:
                hom2 = field.hom(im_gens = [field2.gen()])

                ss = Surface_dict(base_ring = field2)
                index = 0
                P = ConvexPolygons(field2)
                for l,p in self.label_iterator(polygons = True):
                    new_edges = []
                    for i in range(p.num_edges()):
                        new_edges.append( (hom2(coordinates_NF[index]), hom2(coordinates_NF[index+1]) ) )
                        index += 2
                    pp = P(edges = new_edges)
                    ss.add_polygon(pp, label = l)
                ss.change_base_label(self.base_label())
                for (l1,e1),(l2,e2) in self.edge_iterator(gluings = True):
                    ss.change_edge_gluing(l1, e1, l2, e2)
                s = self.__class__(ss)
                if not relabel:
                    if not mutable:
                        s.set_immutable()
                    return s
                # Otherwise we are supposed to relabel. We will make a relabeled copy of s below.
        if new_field is not None:
            from flatsurf.geometry.surface import BaseRingChangedSurface
            s = BaseRingChangedSurface(self,new_field)
        if s is None:
            s = self
        if s.is_finite():
            if relabel:
                return self.__class__(Surface_list(surface=s, copy=not lazy, mutable=mutable))
            else:
                return self.__class__(Surface_dict(surface=s, copy=not lazy, mutable=mutable))
        else:
            if lazy==False:
                raise ValueError("Only lazy copying available for infinite surfaces.")
            if self.underlying_surface().is_mutable():
                raise ValueError("An infinite surface can only be copied if it is immutable.")
            if relabel:
                return self.__class__(Surface_list(surface=s, copy=False, mutable=mutable))
            else:
                return self.__class__(Surface_dict(surface=s, copy=False, mutable=mutable))

    def triangle_flip(self, l1, e1, in_place=False, test=False, direction=None):
        r"""
        Flips the diagonal of the quadrilateral formed by two triangles
        glued together along the provided edge (l1,e1). This can be broken
        into two steps: join along the edge to form a convex quadilateral,
        then cut along the other diagonal. Raises a ValueError if this
        quadrilateral would be non-convex. In this case no changes to the
        surface are made.

        The direction parameter defaults to (0,1). This is used to decide how
        the triangles being glued in are labeled. Let p1 be the triangle
        associated to label l1, and p2 be the triangle associated to l2
        but moved by a similarity to share the edge (l1,e1). Each triangle
        has a exactly one separatrix leaving a vertex which travels in the
        provided direction or its opposite. (For edges we only count as sepatrices
        traveling counter-clockwise around the triangle.) This holds for p1
        and p2 and the separatrices must point in opposite directions.

        The above description gives two new triangles t1 and t2 which must be
        glued in (obtained by flipping the diagonal of the quadrilateral).
        Up to swapping t1 and t2 we can assume the separatrix in t1 in the
        provided direction (or its opposite) points in the same direction as
        that of p1. Further up to cyclic permutation of vertex labels we can
        assume that the separatrices in p1 and t1 start at the vertex with the
        same index (an element of {0,1,2}). The same can be done for p2 and t2.
        We apply the label l1 to t1 and the label l2 to t2. This precisely
        determines how t1 and t2 should be used to replace p1 and p2.

        INPUT:

        - ``l1`` - label of polygon

        - ``e1`` - (integer) edge of the polygon

        - ``in_place`` (boolean) - If True do the flip to the current surface
          which must be mutable. In this case the updated surface will be
          returned.  Otherwise a mutable copy is made and then an edge is
          flipped, which is then returned.

        - ``test`` (boolean) - If True we don't actually flip, and we return
          True or False depending on whether or not the flip would be
          successful.

        - ``direction`` (2-dimensional vector) - Defaults to (0,1). The choice
          of this vector determines how the newly added triangles are labeled.

        EXAMPLES::

            sage: from flatsurf import *

            sage: s = similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, 0), (0, 1)
            sage: s.triangle_flip(0, 0, test=True)
            False
            sage: s.triangle_flip(0, 1, test=True)
            True
            sage: s.triangle_flip(0, 2, test=True)
            False

            sage: s = similarity_surfaces.right_angle_triangle(ZZ(1),ZZ(1))
            sage: from flatsurf.geometry.surface import Surface_list
            sage: s = s.__class__(Surface_list(surface=s, mutable=True))
            sage: try:
            ....:     s.triangle_flip(0,0,in_place=True)
            ....: except ValueError as e:
            ....:     print(e)
            Gluing triangles along this edge yields a non-convex quadrilateral.
            sage: s.triangle_flip(0,1,in_place=True)
            ConeSurface built from 2 polygons
            sage: s.polygon(0)
            Polygon: (0, 0), (1, 1), (0, 1)
            sage: s.polygon(1)
            Polygon: (0, 0), (-1, -1), (0, -1)
            sage: for p in s.edge_iterator(gluings=True):
            ....:     print(p)
            ((0, 0), (1, 0))
            ((0, 1), (0, 2))
            ((0, 2), (0, 1))
            ((1, 0), (0, 0))
            ((1, 1), (1, 2))
            ((1, 2), (1, 1))
            sage: try:
            ....:     s.triangle_flip(0,2,in_place=True)
            ....: except ValueError as e:
            ....:     print(e)
            ....:
            Gluing triangles along this edge yields a non-convex quadrilateral.

            sage: p = polygons((2,0),(-1,3),(-1,-3))
            sage: s = similarity_surfaces.self_glued_polygon(p)
            sage: from flatsurf.geometry.surface import Surface_list
            sage: s = s.__class__(Surface_list(surface=s,mutable=True))
            sage: s.triangle_flip(0,1,in_place=True)
            HalfTranslationSurface built from 1 polygon
            sage: for x in s.label_iterator(polygons=True):
            ....:     print(x)
            (0, Polygon: (0, 0), (-3, -3), (-1, -3))
            sage: for x in s.edge_iterator(gluings=True):
            ....:     print(x)
            ((0, 0), (0, 0))
            ((0, 1), (0, 1))
            ((0, 2), (0, 2))
            sage: TestSuite(s).run()
        """
        if test:
            # Just test if the flip would be successful
            p1=self.polygon(l1)
            if not p1.num_edges()==3:
                return false
            l2,e2 = self.opposite_edge(l1,e1)
            p2 = self.polygon(l2)
            if not p2.num_edges()==3:
                return false
            sim = self.edge_transformation(l2,e2)
            hol = sim( p2.vertex( (e2+2)%3 ) - p1.vertex((e1+2)%3) )
            from flatsurf.geometry.polygon import wedge_product
            return wedge_product(p1.edge((e1+2)%3), hol) > 0 and \
                wedge_product(p1.edge((e1+1)%3), hol) > 0

        if in_place:
            s=self
            assert s.is_mutable(), "Surface must be mutable for in place triangle_flip."
        else:
            s=self.copy(mutable=True)

        p1=s.polygon(l1)
        if not p1.num_edges()==3:
            raise ValueError("The polygon with the provided label is not a triangle.")
        l2,e2 = s.opposite_edge(l1,e1)

        sim = s.edge_transformation(l2,e2)
        m = sim.derivative()
        p2=s.polygon(l2)
        if not p2.num_edges()==3:
            raise ValueError("The polygon opposite the provided edge is not a triangle.")
        P=p1.parent()
        p2=P(vertices=[sim(v) for v in p2.vertices()])

        if direction is None:
            direction=s.vector_space()((0,1))
        # Get vertices corresponding to separatices in the provided direction.
        v1=p1.find_separatrix(direction=direction)[0]
        v2=p2.find_separatrix(direction=direction)[0]
        # Our quadrilateral has vertices labeled:
        # * 0=p1.vertex(e1+1)=p2.vertex(e2)
        # * 1=p1.vertex(e1+2)
        # * 2=p1.vertex(e1)=p2.vertex(e2+1)
        # * 3=p2.vertex(e2+2)
        # Record the corresponding vertices of this quadrilateral.
        q1 = (3+v1-e1-1)%3
        q2 = (2+(3+v2-e2-1)%3)%4

        new_diagonal=p2.vertex((e2+2)%3)-p1.vertex((e1+2)%3)
        # This list will store the new triangles which are being glued in.
        # (Unfortunately, they may not be cyclically labeled in the correct way.)
        new_triangle=[]
        try:
            new_triangle.append(P(edges=[p1.edge((e1+2)%3),p2.edge((e2+1)%3),-new_diagonal]))
            new_triangle.append(P(edges=[p2.edge((e2+2)%3),p1.edge((e1+1)%3),new_diagonal]))
            # The above triangles would be glued along edge 2 to form the diagonal of the quadrilateral being removed.
        except ValueError:
            raise ValueError("Gluing triangles along this edge yields a non-convex quadrilateral.")

        # Find the separatrices of the two new triangles, and in particular which way they point.
        new_sep=[]
        new_sep.append(new_triangle[0].find_separatrix(direction=direction)[0])
        new_sep.append(new_triangle[1].find_separatrix(direction=direction)[0])
        # The quadrilateral vertices corresponding to these separatrices are
        # new_sep[0]+1 and (new_sep[1]+3)%4 respectively.

        # i=0 if the new_triangle[0] should be labeled l1 and new_triangle[1] should be labeled l2.
        # i=1 indicates the opposite labeling.
        if new_sep[0]+1==q1:
            # For debugging:
            assert (new_sep[1]+3)%4==q2, \
                "Bug: new_sep[1]="+str(new_sep[1])+" and q2="+str(q2)
            i=0
        else:
            # For debugging:
            assert (new_sep[1]+3)%4==q1
            assert new_sep[0]+1==q2
            i=1

        # These quantities represent the cyclic relabeling of triangles needed.
        cycle1 = (new_sep[i]-v1+3)%3
        cycle2 = (new_sep[1-i]-v2+3)%3

        # This will be the new triangle with label l1:
        tri1=P(edges=[new_triangle[i].edge(cycle1), \
                      new_triangle[i].edge((cycle1+1)%3), \
                      new_triangle[i].edge((cycle1+2)%3)])
        # This will be the new triangle with label l2:
        tri2=P(edges=[new_triangle[1-i].edge(cycle2), \
                      new_triangle[1-i].edge((cycle2+1)%3), \
                      new_triangle[1-i].edge((cycle2+2)%3)])
        # In the above, edge 2-cycle1 of tri1 would be glued to edge 2-cycle2 of tri2
        diagonal_glue_e1=2-cycle1
        diagonal_glue_e2=2-cycle2

        # FOR CATCHING BUGS:
        assert p1.find_separatrix(direction=direction)==tri1.find_separatrix(direction=direction)
        assert p2.find_separatrix(direction=direction)==tri2.find_separatrix(direction=direction)

        # Two opposite edges will not change their labels (label,edge) under our regluing operation.
        # The other two opposite ones will change and in fact they change labels.
        # The following finds them (there are two cases).
        # At the end of the if statement, the following will be true:
        # * new_glue_e1 and new_glue_e2 will be the edges of the new triangle with label l1 and l2 which need regluing.
        # * old_e1 and old_e2 will be the corresponding edges of the old triangles.
        # (Note that labels are swapped between the pair. The appending 1 or 2 refers to the label used for the triangle.)
        if p1.edge(v1)==tri1.edge(v1):
            # We don't have to worry about changing gluings on edge v1 of the triangles with label l1
            # We do have to worry about the following edge:
            new_glue_e1=3-diagonal_glue_e1-v1 # returns the edge which is neither diagonal_glue_e1 nor v1.
            # This corresponded to the following old edge:
            old_e1 = 3 - e1 - v1 # Again this finds the edge which is neither e1 nor v1
        else:
            temp = (v1+2)%3
            # FOR CATCHING BUGS:
            assert p1.edge(temp)==tri1.edge(temp)
            # We don't have to worry about changing gluings on edge (v1+2)%3 of the triangles with label l1
            # We do have to worry about the following edge:
            new_glue_e1=3-diagonal_glue_e1-temp # returns the edge which is neither diagonal_glue_e1 nor temp.
            # This corresponded to the following old edge:
            old_e1 = 3 - e1 - temp # Again this finds the edge which is neither e1 nor temp
        if p2.edge(v2)==tri2.edge(v2):
            # We don't have to worry about changing gluings on edge v2 of the triangles with label l2
            # We do have to worry about the following edge:
            new_glue_e2=3-diagonal_glue_e2-v2 # returns the edge which is neither diagonal_glue_e2 nor v2.
            # This corresponded to the following old edge:
            old_e2 = 3 - e2 - v2 # Again this finds the edge which is neither e2 nor v2
        else:
            temp = (v2+2)%3
            # FOR CATCHING BUGS:
            assert p2.edge(temp)==tri2.edge(temp)
            # We don't have to worry about changing gluings on edge (v2+2)%3 of the triangles with label l2
            # We do have to worry about the following edge:
            new_glue_e2=3-diagonal_glue_e2-temp # returns the edge which is neither diagonal_glue_e2 nor temp.
            # This corresponded to the following old edge:
            old_e2 = 3 - e2 - temp # Again this finds the edge which is neither e2 nor temp

        # remember the old gluings.
        old_opposite1 = s.opposite_edge(l1, old_e1)
        old_opposite2 = s.opposite_edge(l2, old_e2)

        # We make changes to the underlying surface
        us=s.underlying_surface()

        # Replace the triangles.
        us.change_polygon(l1,tri1)
        us.change_polygon(l2,tri2)
        # Glue along the new diagonal of the quadrilateral
        us.change_edge_gluing(l1,diagonal_glue_e1,
                             l2,diagonal_glue_e2)
        # Now we deal with that pair of opposite edges of the quadrilateral that need regluing.
        # There are some special cases:
        if old_opposite1==(l2,old_e2):
            # These opposite edges were glued to each other.
            # Do the same in the new surface:
            us.change_edge_gluing(l1,new_glue_e1,
                                 l2,new_glue_e2)
        else:
            if old_opposite1==(l1,old_e1):
                # That edge was "self-glued".
                us.change_edge_gluing(l2,new_glue_e2,
                                     l2,new_glue_e2)
            else:
                # The edge (l1,old_e1) was glued in a standard way.
                # That edge now corresponds to (l2,new_glue_e2):
                us.change_edge_gluing(l2,new_glue_e2,
                                     old_opposite1[0],old_opposite1[1])
            if old_opposite2==(l2,old_e2):
                # That edge was "self-glued".
                us.change_edge_gluing(l1,new_glue_e1,
                                     l1,new_glue_e1)
            else:
                # The edge (l2,old_e2) was glued in a standard way.
                # That edge now corresponds to (l1,new_glue_e1):
                us.change_edge_gluing(l1,new_glue_e1,
                                     old_opposite2[0],old_opposite2[1])
        return s

    def join_polygons(self, p1, e1, test=False, in_place=False):
        r"""
        Join polygons across the provided edge (p1,e1). By default,
        it returns the surface obtained by joining the two polygons
        together. It raises a ValueError if gluing the two polygons
        together results in a non-convex polygon. This is done to the
        current surface if in_place is True, and otherwise a mutable
        copy is made and then modified.

        If test is True then instead of changing the surface, it just
        checks to see if the change would be successful and returns
        True if successful or False if not.

        EXAMPLES::

            sage: from flatsurf import *
            sage: ss=translation_surfaces.ward(3)
            sage: s=ss.copy(mutable=True)
            sage: s.join_polygons(0,0, in_place=True)
            TranslationSurface built from 2 polygons
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (0, 2*a), (-1, a)
            sage: s.join_polygons(0,4, in_place=True)
            TranslationSurface built from 1 polygon
            sage: print(s.polygon(0))
            Polygon: (0, 0), (1, -a), (2, 0), (3, a), (2, 2*a), (1, 3*a), (0, 2*a), (-1, a)
        """
        poly1=self.polygon(p1)
        p2,e2 = self.opposite_edge(p1,e1)
        poly2=self.polygon(p2)
        if p1==p2:
            if test:
                return False
            else:
                raise ValueError("Can't glue polygon to itself.")
        t=self.edge_transformation(p2, e2)
        dt=t.derivative()
        vs = []
        edge_map={} # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))
        ne=poly2.num_edges()
        for i in range(1,ne):
            ee=(e2+i)%ne
            edge_map[len(vs)]=(p2,ee)
            vs.append(dt * poly2.edge( ee ))
        for i in range(e1+1, poly1.num_edges()):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))

        try:
            new_polygon = ConvexPolygons(self.base_ring())(vs)
        except (ValueError, TypeError):
            if test:
                return False
            else:
                raise ValueError("Joining polygons along this edge results in a non-convex polygon.")

        if test:
            # Gluing would be successful
            return True

        # Now no longer testing. Do the gluing.
        if in_place:
            ss=self
        else:
            ss=self.copy(mutable=True)
        s=ss.underlying_surface()

        inv_edge_map={}
        for key, value in iteritems(edge_map):
            inv_edge_map[value]=(p1,key)

        glue_list=[]
        for i in range(len(vs)):
            p3,e3 = edge_map[i]
            p4,e4 = self.opposite_edge(p3,e3)
            if p4 == p1 or p4 == p2:
                glue_list.append(inv_edge_map[(p4,e4)])
            else:
                glue_list.append((p4,e4))

        if s.base_label()==p2:
             s.change_base_label(p1)

        s.remove_polygon(p2)

        s.change_polygon(p1, new_polygon, glue_list)

        return ss

    def subdivide_polygon(self, p, v1, v2, test=False, new_label=None):
        r"""
        Cut the polygon with label p along the diagonal joining vertex
        v1 to vertex v2. This cuts p into two polygons, one will keep the same
        label. The other will get a new label, which can be provided
        via new_label. Otherwise a default new label will be provided.
        If test=False, then the surface will be changed (in place). If
        test=True, then it just checks to see if the change would be successful

        The convention is that the resulting subdivided polygon which has an oriented
        edge going from the original vertex v1 to vertex v2 will keep the label p.
        The other polygon will get a new label.

        The change will be done in place.
        """
        poly=self.polygon(p)
        ne=poly.num_edges()
        if v1<0 or v2<0 or v1>=ne or v2>=ne:
            if test:
                return False
            else:
                raise ValueError('Provided vertices out of bounds.')
        if abs(v1-v2)<=1 or abs(v1-v2)>=ne-1:
            if test:
                return False
            else:
                raise ValueError('Provided diagonal is not actually a diagonal.')

        if v2<v1:
            v2=v2+ne

        newedges1=[poly.vertex(v2)-poly.vertex(v1)]
        for i in range(v2, v1+ne):
            newedges1.append(poly.edge(i))
        newpoly1 = ConvexPolygons(self.base_ring())(newedges1)

        newedges2=[poly.vertex(v1)-poly.vertex(v2)]
        for i in range(v1,v2):
            newedges2.append(poly.edge(i))
        newpoly2 = ConvexPolygons(self.base_ring())(newedges2)

        # Store the old gluings
        old_gluings = {(p,i): self.opposite_edge(p,i) for i in range(ne)}

        # Update the polygon with label p, add a new polygon.
        self.underlying_surface().change_polygon(p, newpoly1)
        if new_label is None:
            new_label = self.underlying_surface().add_polygon(newpoly2)
        else:
            new_label = self.underlying_surface().add_polygon(newpoly2, label=new_label)
        # This gluing is the diagonal we used.
        self.underlying_surface().change_edge_gluing(p, 0, new_label, 0)

        # Setup conversion from old to new labels.
        old_to_new_labels={}
        for i in range(v1, v2):
            old_to_new_labels[(p,i%ne)]=(new_label,i-v1+1)
        for i in range(v2, ne+v1):
            old_to_new_labels[(p,i%ne)]=(p,i-v2+1)

        for e in range(1, newpoly1.num_edges()):
            pair = old_gluings[(p,(v2+e-1)%ne)]
            if pair in old_to_new_labels:
                pair = old_to_new_labels[pair]
            self.underlying_surface().change_edge_gluing(p, e, pair[0], pair[1])

        for e in range(1, newpoly2.num_edges()):
            pair = old_gluings[(p,(v1+e-1)%ne)]
            if pair in old_to_new_labels:
                pair = old_to_new_labels[pair]
            self.underlying_surface().change_edge_gluing(new_label, e, pair[0], pair[1])

    def singularity(self, l, v, limit=None):
        r"""
        Represents the Singularity associated to the v-th vertex of the polygon with
        label l.

        If the surface is infinite, the limit needs to be set. In this case the construction
        of the singularity is successful if the sequence of vertices hit by passing through
        edges closes up in limit or less steps.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.square_torus()
            sage: pc = s.minimal_cover(cover_type="planar")
            sage: pc.singularity(pc.base_label(),0)
            Traceback (most recent call last):
            ...
            ValueError: need a limit when working with an infinite surface
            sage: pc.singularity(pc.base_label(),0,limit=4)
            singularity with vertex equivalence class frozenset(...)
        """
        return Singularity(self,l,v,limit)

    def point(self, label, point, ring=None, limit=None):
        r"""
        Return a point in this surface.

        INPUT:

        - ``label`` - label of the polygon

        - ``point`` - coordinates of the point inside the polygon

        - ``ring`` (optional) - a ring for the coordinates

        - ``limit`` (optional) - undocumented (only necessary if the point corresponds
          to a singularity in an infinite surface)

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = translation_surfaces.square_torus()
            sage: pc = s.minimal_cover(cover_type="planar")
            sage: pc.surface_point(pc.base_label(),(0,0))
            Traceback (most recent call last):
            ...
            ValueError: need a limit when working with an infinite surface
            sage: pc.surface_point(pc.base_label(),(1,0),limit=4)
            Surface point with 4 coordinate representations
            sage: z = pc.surface_point(pc.base_label(),(sqrt(2)-1,sqrt(3)-1),ring=AA)
            sage: next(iter(z.coordinates(z.labels()[0]))).parent()
            Vector space of dimension 2 over Algebraic Real Field
        """
        return SurfacePoint(self, label, point, ring=ring, limit=limit)

    # TODO: deprecate
    surface_point = point

    def ramified_cover(self, degree, data):
        r"""
        Build a ramified cover of this surface with given ``degree`` and ramification ``data``.

        INPUT:

        - ``degree`` (integer) -- the degree of the cover

        - ``data`` -- dictionary that associates a pair ``(polygon_label, edge_number)`` a permutation
          of ``{1, 2, ..., d}``

        EXAMPLES:

        The L-shape origami::

            sage: import flatsurf
            sage: T = flatsurf.translation_surfaces.square_torus()
            sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
            TranslationSurface built from 3 polygons
            sage: O = T.ramified_cover(3, {(0,0): '(1,2)', (0,1): '(1,3)'})
            sage: O.stratum()
            H_2(2)

        TESTS::

            sage: import flatsurf
            sage: T = flatsurf.translation_surfaces.square_torus()
            sage: T.ramified_cover(3, {(0,0): '(1,2)', (0,2): '(1,3)'})
            Traceback (most recent call last):
            ...
            ValueError: inconsistent covering data
        """
        if not self.is_finite():
            raise ValueError("this method is only available for finite surfaces")
        return type(self)(self._s.ramified_cover(degree, data))

    def minimal_cover(self, cover_type = "translation"):
        r"""
        Return the minimal translation or half-translation cover of the surface.

        Cover type may be either "translation", "half-translation" or "planar".

        The minimal planar cover of a surface S is the smallest cover C so that
        the developing map from the universal cover U to the plane induces a
        well defined map from C to the plane. This is an infinite translation
        surface that is naturally a branched cover of the plane.

        EXAMPLES::

            sage: from flatsurf.geometry.surface import Surface_list
            sage: s = Surface_list(QQ)
            sage: from flatsurf.geometry.polygon import polygons
            sage: square = polygons.square(field=QQ)
            sage: s.add_polygon(square)
            0
            sage: s.change_edge_gluing(0,0,0,1)
            sage: s.change_edge_gluing(0,2,0,3)
            sage: from flatsurf.geometry.cone_surface import ConeSurface
            sage: cs = ConeSurface(s)
            sage: ts = cs.minimal_cover(cover_type="translation")
            sage: ts
            TranslationSurface built from 4 polygons
            sage: hts = cs.minimal_cover(cover_type="half-translation")
            sage: hts
            HalfTranslationSurface built from 2 polygons
            sage: TestSuite(hts).run()
            sage: ps = cs.minimal_cover(cover_type="planar")
            sage: ps
            TranslationSurface built from infinitely many polygons
            sage: TestSuite(ps).run(skip="_test_pickling")

            sage: from flatsurf import similarity_surfaces
            sage: S = similarity_surfaces.example()
            sage: T = S.minimal_cover(cover_type="translation")
            sage: T
            TranslationSurface built from infinitely many polygons
            sage: T.polygon(T.base_label())
            Polygon: (0, 0), (2, -2), (2, 0)
        """
        if cover_type == "translation":
            from flatsurf.geometry.translation_surface import TranslationSurface
            from flatsurf.geometry.minimal_cover import MinimalTranslationCover
            return TranslationSurface(MinimalTranslationCover(self))
        if cover_type == "half-translation":
            from flatsurf.geometry.half_translation_surface import HalfTranslationSurface
            from flatsurf.geometry.minimal_cover import MinimalHalfTranslationCover
            return HalfTranslationSurface(MinimalHalfTranslationCover(self))
        if cover_type == "planar":
            from flatsurf.geometry.translation_surface import TranslationSurface
            from flatsurf.geometry.minimal_cover import MinimalPlanarCover
            return TranslationSurface(MinimalPlanarCover(self))
        raise ValueError("Provided cover_type is not supported.")

    def minimal_translation_cover(self):
        r"""
        Return the minimal translation cover.

        "Be careful that if the surface is not built from one polygon, this is
        not the smallest translation cover of the surface." - Vincent

        "I disagree with the prior statement. Can you provide an example?" -Pat
        """
        from sage.misc.superseded import deprecation
        deprecation(13109, "minimal_translation_cover is deprecated. Use minimal_cover(cover_type = \"translation\") instead.")
        from flatsurf.geometry.translation_surface import MinimalTranslationCover, TranslationSurface
        return TranslationSurface(MinimalTranslationCover(self))

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), 2)

    def fundamental_group(self, base_label=None):
        r"""
        Return the fundamental group of this surface.
        """
        if not self.is_finite():
            raise ValueError("the method only work for finite surfaces")
        if base_label is None:
            base_label = self.base_label()
        from .fundamental_group import FundamentalGroup
        return FundamentalGroup(self, base_label)

    def tangent_bundle(self, ring=None):
        r"""
        Return the tangent bundle

        INPUT:

        - ``ring`` -- an optional field (defaults to the coordinate field of the
          surface)
        """
        if ring is None:
            ring = self.base_ring()

        try:
            return self._tangent_bundle_cache[ring]
        except AttributeError:
            self._tangent_bundle_cache = {}
        except KeyError:
            pass

        from .tangent_bundle import SimilaritySurfaceTangentBundle
        self._tangent_bundle_cache[ring] = SimilaritySurfaceTangentBundle(self, ring)
        return self._tangent_bundle_cache[ring]

    def tangent_vector(self, lab, p, v, ring=None):
        r"""
        Return a tangent vector.

        INPUT:

        - ``lab`` -- label of a polygon

        - ``p`` -- coordinates of a point in the polygon

        - ``v`` -- coordinates of a vector in R^2

        EXAMPLES::

            sage: from flatsurf.geometry.chamanara import chamanara_surface
            sage: S = chamanara_surface(1/2)
            sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,1))
            SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, 1)
            sage: K.<sqrt2> = QuadraticField(2)
            sage: S.tangent_vector(S.base_label(), (1/2,1/2), (1,sqrt2), ring=K)
            SimilaritySurfaceTangentVector in polygon (1, -1, 0) based at (1/2, -3/2) with vector (1, sqrt2)
        """
        p = vector(p)
        v = vector(v)

        if p.parent().dimension() != 2 or v.parent().dimension() != 2:
            raise ValueError("p (={!r}) and v (={!v}) should have two coordinates")

        if ring is None:
            ring = self.base_ring()
            try:
                return self.tangent_bundle(ring)(lab, p, v)
            except TypeError:
                raise TypeError("Use the ring=??? option to construct tangent vectors in other field different from the base_ring().")
            # Old version seemed to be to accepting of inputs (eg, from Symbolic Ring)
            #R = p.base_ring()
            #if R != v.base_ring():
            #    from sage.structure.element import get_coercion_model
            #    cm = get_coercion_model()
            #    R = cm.common_parent(R, v.base_ring())
            #    p = p.change_ring(R)
            #    v = v.change_ring(R)

            #R2 = self.base_ring()
            #if R != R2:
            #    if R2.has_coerce_map_from(R):
            #        p = p.change_ring(R2)
            #        v = v.change_ring(R2)
            #        R = R2
            #    elif not R.has_coerce_map_from(R2):
            #        raise ValueError("not able to find a common ring for arguments")
            #return self.tangent_bundle(R)(lab, p, v)
        else:
            return self.tangent_bundle(ring)(lab, p, v)

    def reposition_polygons(self, in_place=False, relabel=False):
        r"""
        We choose a maximal tree in the dual graph of the decomposition into
        polygons, and ensure that the gluings between two polygons joined by
        an edge in this tree is by translation.

        This guarantees that the group generated by the edge identifications is
        minimal among representions of the surface. In particular, if for instance
        you have a translation surface which is anot representable as a translation
        surface (because polygons are presented with rotations) then after this
        change it will be representable as a translation surface.
        """
        if not self.is_finite():
            raise NotImplementedError("Only implemented for finite surfaces.")
        if in_place:
            if not self.is_mutable():
                raise ValueError("reposition_polygons in_place is only available "+\
                    "for mutable surfaces.")
            s=self
        else:
            s=self.copy(relabel=relabel, mutable=True)
        w=s.walker()
        from flatsurf.geometry.similarity import SimilarityGroup
        S=SimilarityGroup(self.base_ring())
        identity=S.one()
        it = iter(w)
        label = next(it)
        changes = {label:identity}
        for label in it:
            edge = w.edge_back(label)
            label2,edge2 = s.opposite_edge(label, edge)
            changes[label] = changes[label2] * s.edge_transformation(label,edge)
        it = iter(w)
        # Skip the base label:
        label = next(it)
        for label in it:
            p = s.polygon(label)
            p = changes[label].derivative()*p
            s.underlying_surface().change_polygon(label,p)
        return s

    def triangulation_mapping(self):
        r"""
        Return a SurfaceMapping triangulating the suface or None if the surface is already triangulated.
        """
        from flatsurf.geometry.mappings import triangulation_mapping
        return triangulation_mapping(self)

    def triangulate(self, in_place=False, label = None, relabel=False):
        r"""
        Return a triangulated version of this surface. (This may be mutable
        or not depending on the input.)

        If label=None (as default) all polygons are triangulated. Otherwise,
        label should be a polygon label. In this case, just this polygon
        is split into triangles.

        This is done in place if in_place is True (defaults to False).

        If we are not doing triangulation in_place, then we must make a copy.
        This can be a relabeled copy (indexed by the non-negative ints)
        or a label preserving copy. The copy is relabeled if relabel=True
        (default False).

        EXAMPLES::

            sage: from flatsurf import *
            sage: s=translation_surfaces.mcmullen_L(1,1,1,1)
            sage: ss=s.triangulate()
            sage: gs=ss.graphical_surface()
            sage: gs.make_all_visible()
            sage: print(gs)
            Graphical version of Similarity Surface TranslationSurface built from 6 polygons

        A non-strictly convex example that caused trouble:

            sage: from flatsurf import *
            sage: s=similarity_surfaces.self_glued_polygon(polygons(edges=[(1,1),(-3,-1),(1,0),(1,0)]))
            sage: s=s.triangulate()
            sage: s.polygon(0).num_edges()
            3
        """
        if label is None:
            # We triangulate the whole surface
            if self.is_finite():
                # Store the current labels.
                labels = [label for label in self.label_iterator()]
                if in_place:
                    s=self
                else:
                    s=self.copy(mutable=True)
                # Subdivide each polygon in turn.
                for l in labels:
                    s = s.triangulate(in_place=True, label=l)
                return s
            else:
                if in_place:
                    raise ValueError("You can't triangulate an infinite surface in place.")
                from flatsurf.geometry.delaunay import LazyTriangulatedSurface
                return self.__class__(LazyTriangulatedSurface(self))
        else:
            poly = self.polygon(label)
            n=poly.num_edges()
            if n>3:
                if in_place:
                    s=self
                else:
                    s=self.copy(mutable=True)
            else:
                # This polygon is already a triangle.
                return self
            from flatsurf.geometry.polygon import wedge_product
            for i in range(n-3):
                poly = s.polygon(label)
                n=poly.num_edges()
                for i in range(n):
                    e1=poly.edge(i)
                    e2=poly.edge((i+1)%n)
                    if wedge_product(e1,e2) != 0:
                        # This is in case the polygon is a triangle with subdivided edge.
                        e3=poly.edge((i+2)%n)
                        if wedge_product(e1+e2,e3) != 0:
                            s.subdivide_polygon(label,i,(i+2)%n)
                            break
            return s
        raise RuntimeError("Failed to return anything!")

    def _edge_needs_flip(self,p1,e1):
        r"""
        Returns -1 if the the provided edge incident to two triangles which
        should be flipped to get closer to the Delaunay decomposition.
        Returns 0 if the quadrilateral formed by the triangles is inscribed
        in a circle, and returns 1 otherwise.

        A ValueError is raised if the edge is not indident to two triangles.
        """
        p2,e2=self.opposite_edge(p1,e1)
        poly1=self.polygon(p1)
        poly2=self.polygon(p2)
        if poly1.num_edges()!=3 or poly2.num_edges()!=3:
            raise ValueError("Edge must be adjacent to two triangles.")
        from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
        sim1=similarity_from_vectors(poly1.edge(e1+2),-poly1.edge(e1+1))
        sim2=similarity_from_vectors(poly2.edge(e2+2),-poly2.edge(e2+1))
        sim=sim1*sim2
        return sim[1][0]<0

    def _edge_needs_join(self,p1,e1):
        r"""
        Returns -1 if the the provided edge incident to two triangles which
        should be flipped to get closer to the Delaunay decomposition.
        Returns 0 if the quadrilateral formed by the triangles is inscribed
        in a circle, and returns 1 otherwise.

        A ValueError is raised if the edge is not indident to two triangles.
        """
        p2,e2=self.opposite_edge(p1,e1)
        poly1=self.polygon(p1)
        poly2=self.polygon(p2)
        from flatsurf.geometry.matrix_2x2 import similarity_from_vectors
        sim1=similarity_from_vectors(poly1.vertex(e1) - poly1.vertex(e1+2),\
            -poly1.edge(e1+1))
        sim2=similarity_from_vectors(poly2.vertex(e2) - poly2.vertex(e2+2),\
            -poly2.edge(e2+1))
        sim=sim1*sim2
        from sage.functions.generalized import sgn
        return sim[1][0]==0

    def delaunay_single_flip(self):
        r"""
        Does a single in place flip of a triangulated mutable surface.
        """
        if not self.is_finite():
            raise NotImplementedError("Not implemented for infinite surfaces.")
        lc = self._label_comparator()
        for (l1,e1),(l2,e2) in self.edge_iterator(gluings=True):
            if (lc.lt(l1,l2) or (l1==l2 and e1<=e2)) and self._edge_needs_flip(l1,e1):
                self.triangle_flip(l1, e1, in_place=True)
                return True
        return False

    def is_delaunay_triangulated(self, limit=None):
        r"""
        Return if the surface is triangulated and the triangulation is Delaunay.
        If limit is set, then it checks this only limit many edges.
        Limit must be set for infinite surfaces.
        """
        if limit is None:
            if not self.is_finite():
                raise NotImplementedError("A limit must be set for infinite surfaces.")
            limit = self.num_edges()
        count = 0
        for (l1,e1),(l2,e2) in self.edge_iterator(gluings=True):
            if count >= limit:
                break
            count =  count+1
            if self.polygon(l1).num_edges()!=3:
                print("Polygon with label "+str(l1)+" is not a triangle.")
                return False
            if self.polygon(l2).num_edges()!=3:
                print("Polygon with label "+str(l2)+" is not a triangle.")
                return False
            if self._edge_needs_flip(l1,e1):
                print("Edge "+str((l1,e1))+" needs to be flipped.")
                print("This edge is glued to "+str((l2,e2))+".")
                return False
        return True

    def is_delaunay_decomposed(self, limit=None):
        r"""
        Return if the decomposition of the surface into polygons is Delaunay.
        If limit is set, then it checks this only limit many polygons.
        Limit must be set for infinite surfaces.
        """
        if limit is None:
            if not self.is_finite():
                raise NotImplementedError("A limit must be set for infinite surfaces.")
            limit = self.num_polygons()
        count = 0
        for (l1,p1) in self.label_iterator(polygons=True):
            try:
                c1=p1.circumscribing_circle()
            except ValueError:
                # p1 is not circumscribed
                return False
            for e1 in range(p1.num_edges()):
                c2=self.edge_transformation(l1,e1)*c1
                l2,e2=self.opposite_edge(l1,e1)
                if c2.point_position(self.polygon(l2).vertex(e2+2))!=-1:
                    # The circumscribed circle developed into the adjacent polygon
                    # contains a vertex in its interior or boundary.
                    return False
            return True

    def delaunay_triangulation(self, triangulated=False, in_place=False, limit=None, direction=None, relabel=False):
        r"""
        Returns a Delaunay triangulation of a surface, or make some
        triangle flips to get closer to the Delaunay decomposition.

        INPUT:

        - ``triangulated`` (boolean) - If true, the algorithm assumes the
          surface is already triangulated. It does this without verification.

        - ``in_place`` (boolean) - If true, the triangulating and the
          triangle flips are done in place.  Otherwise, a mutable copy of the
          surface is made.

        - ``limit`` (None or Integer) - If None, this will return a
          Delaunay triangulation. If limit is an integer 1 or larger, then at
          most limit many diagonal flips will be done.

        - ``direction`` (None or Vector) - with two entries in the base field
            Used to determine labels when a pair of triangles is flipped. Each triangle
            has a unique separatrix which points in the provided direction or its
            negation. As such a vector determines a sign for each triangle.
            A pair of adjacent triangles have opposite signs. Labels are chosen
            so that this sign is preserved (as a function of labels).

        - ``relabel`` (boolean) - If in_place is False, then a copy must be
          made. By default relabel is False and labels will be respected by
          this copy. If relabel is True then polygons will be reindexed in an
          arbitrary way by the non-negative integers.

        EXAMPLES::

            sage: from flatsurf import *
            sage: from flatsurf.geometry.delaunay import *

            sage: m = matrix([[2,1],[1,1]])
            sage: s = m*translation_surfaces.infinite_staircase()
            sage: ss = s.delaunay_triangulation(relabel=True)
            sage: ss.base_label()
            0
            sage: ss.polygon(0)
            Polygon: (0, 0), (1, 1), (0, 1)
            sage: TestSuite(ss).run(skip="_test_pickling")
            sage: ss.is_delaunay_triangulated(limit=10)
            True
        """
        if not self.is_finite() and limit is None:
            if in_place:
                raise ValueError("in_place delaunay triangulation is not possible for infinite surfaces unless a limit is set.")
            if self.underlying_surface().is_mutable():
                raise ValueError("delaunay_triangulation only works on infinite "+\
                    "surfaces if they are immutable or if a limit is set.")
            from flatsurf.geometry.delaunay import LazyDelaunayTriangulatedSurface
            return self.__class__(LazyDelaunayTriangulatedSurface( \
                self,direction=direction, relabel=relabel))
        if in_place and not self.is_mutable():
            raise ValueError("in_place delaunay_triangulation only defined for mutable surfaces")
        if triangulated:
            if in_place:
                s=self
            else:
                s=self.copy(mutable=True, relabel=False)
        else:
            if in_place:
                s=self
                self.triangulate(in_place=True)
            else:
                s=self.copy(relabel=True,mutable=True)
                s.triangulate(in_place=True)
        loop=True
        if direction is None:
            base_ring = self.base_ring()
            direction = self.vector_space()( (base_ring.zero(), base_ring.one()) )
        else:
            assert not direction.is_zero()
        if s.is_finite() and limit is None:
            from collections import deque
            unchecked_labels=deque(label for label in s.label_iterator())
            checked_labels = set()
            while unchecked_labels:
                label = unchecked_labels.popleft()
                flipped=False
                for edge in range(3):
                    if s._edge_needs_flip(label,edge):
                        # Record the current opposite edge:
                        label2,edge2=s.opposite_edge(label,edge)
                        # Perform the flip.
                        s.triangle_flip(label, edge, in_place=True, direction=direction)
                        # Move the opposite polygon to the list of labels we need to check.
                        if label2 != label:
                            try:
                                checked_labels.remove(label2)
                                unchecked_labels.append(label2)
                            except KeyError:
                                # Occurs if label2 is not in checked_labels
                                pass
                        flipped=True
                        break
                if flipped:
                    unchecked_labels.append(label)
                else:
                    checked_labels.add(label)
            return s
        else:
            # Old method for infinite surfaces, or limits.
            count=0
            lc = self._label_comparator()
            while loop:
                loop=False
                for (l1,e1),(l2,e2) in s.edge_iterator(gluings=True):
                    if (lc.lt(l1,l2) or (l1==l2 and e1<=e2)) and s._edge_needs_flip(l1,e1):
                        s.triangle_flip(l1, e1, in_place=True, direction=direction)
                        count += 1
                        if not limit is None and count>=limit:
                            return s
                        loop=True
                        break
            return s

    def delaunay_single_join(self):
        if not self.is_finite():
            raise NotImplementedError("Not implemented for infinite surfaces.")
        lc = self._label_comparator()
        for (l1,e1),(l2,e2) in self.edge_iterator(gluings=True):
            if (lc.lt(l1,l2) or (l1==l2 and e1<=e2)) and self._edge_needs_join(l1,e1):
                self.join_polygons(l1, e1, in_place=True)
                return True
        return False


    def delaunay_decomposition(self, triangulated=False, \
            delaunay_triangulated=False, in_place=False, direction=None,\
            relabel=False):
        r"""
        Return the Delaunay Decomposition of this surface.

        INPUT:

        - ``triangulated`` (boolean) - If true, the algorithm assumes the
          surface is already triangulated. It does this without verification.

        - ``delaunay_triangulated`` (boolean) - If true, the algorithm assumes
          the surface is already delaunay_triangulated. It does this without
          verification.

        - ``in_place`` (boolean) - If true, the triangulating and the triangle
          flips are done in place. Otherwise, a mutable copy of the surface is
          made.

        - ``relabel`` (None or Integer) - If in_place is False, then a copy
          must be made of the surface.  If relabel is False (as default), the
          copy has the same labels as the original surface. Note that in this
          case, labels will be added if it is necessary to subdivide polygons
          into triangles.  If relabel is True, the new surface will have
          polygons labeled by the non-negative integers in an arbitrary way.

        - ``direction`` - (None or Vector with two entries in the base field) -
          Used to determine labels when a pair of triangles is flipped. Each triangle
          has a unique separatrix which points in the provided direction or its
          negation. As such a vector determines a sign for each triangle.
          A pair of adjacent triangles have opposite signs. Labels are chosen
          so that this sign is preserved (as a function of labels).

        EXAMPLES::

            sage: from flatsurf import *
            sage: s0 = translation_surfaces.octagon_and_squares()
            sage: a = s0.base_ring().gens()[0]
            sage: m = Matrix([[1,2+a],[0,1]])
            sage: s = m*s0
            sage: s = s.triangulate()
            sage: ss = s.delaunay_decomposition(triangulated=True)
            sage: ss.num_polygons()
            3

            sage: p = polygons((4,0),(-2,1),(-2,-1))
            sage: s0 = similarity_surfaces.self_glued_polygon(p)
            sage: s = s0.delaunay_decomposition()
            sage: TestSuite(s).run()

            sage: m = matrix([[2,1],[1,1]])
            sage: s = m*translation_surfaces.infinite_staircase()
            sage: ss = s.delaunay_decomposition()
            sage: ss.base_label()
            0
            sage: ss.polygon(0)
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: TestSuite(ss).run(skip="_test_pickling")
            sage: ss.is_delaunay_decomposed(limit=10)
            True
        """
        if not self.is_finite():
            if in_place:
                raise ValueError("in_place delaunay_decomposition is not possible for infinite surfaces.")
            if self.underlying_surface().is_mutable():
                raise ValueError("delaunay_decomposition only works on infinite "+\
                    "surfaces if they are immutable.")
            from flatsurf.geometry.delaunay import LazyDelaunaySurface
            return self.__class__(LazyDelaunaySurface( \
                self,direction=direction, relabel=relabel))
        if in_place:
            s=self
        else:
            s=self.copy(mutable=True, relabel=relabel)
        if not delaunay_triangulated:
            s.delaunay_triangulation(triangulated=triangulated, in_place=True, \
                direction=direction)
        # Now s is Delaunay Triangulated
        loop=True
        lc = self._label_comparator()
        while loop:
            loop=False
            for (l1,e1),(l2,e2) in s.edge_iterator(gluings=True):
                if (lc.lt(l1,l2) or (l1==l2 and e1<=e2)) and s._edge_needs_join(l1,e1):
                    s.join_polygons(l1, e1, in_place=True)
                    loop=True
                    break
        return s

    def saddle_connections(self, squared_length_bound, initial_label=None, initial_vertex=None, sc_list=None, check=False):
        r"""
        Returns a list of saddle connections on the surface whose length squared is less than or equal to squared_length_bound.
        The length of a saddle connection is measured using holonomy from polygon in which the trajectory starts.

        If initial_label and initial_vertex are not provided, we return all saddle connections satisfying the bound condition.

        If initial_label and initial_vertex are provided, it only provides saddle connections emanating from the corresponding
        vertex of a polygon. If only initial_label is provided, the added saddle connections will only emanate from the
        corresponding polygon.

        If sc_list is provided the found saddle connections are appended to this list and the resulting list is returned.

        If check==True it uses the checks in the SaddleConnection class to sanity check our results.

        EXAMPLES::
            sage: from flatsurf import *
            sage: s = translation_surfaces.square_torus()
            sage: sc_list = s.saddle_connections(13, check=True)
            sage: len(sc_list)
            32
        """
        assert squared_length_bound > 0
        if sc_list is None:
            sc_list = []
        if initial_label is None:
            assert self.is_finite()
            assert initial_vertex is None, "If initial_label is not provided, then initial_vertex must not be provided either."
            for label in self.label_iterator():
                self.saddle_connections(squared_length_bound, initial_label=label, sc_list=sc_list)
            return sc_list
        if initial_vertex is None:
            for vertex in range( self.polygon(initial_label).num_edges() ):
                self.saddle_connections(squared_length_bound, initial_label=initial_label, initial_vertex=vertex, sc_list=sc_list)
            return sc_list

        # Now we have a specified initial_label and initial_vertex
        SG = SimilarityGroup(self.base_ring())
        start_data = (initial_label, initial_vertex)
        circle = Circle(self.vector_space().zero(), squared_length_bound, base_ring =   self.base_ring())
        p = self.polygon(initial_label)
        v = p.vertex(initial_vertex)
        last_sim = SG(-v[0],-v[1])

        # First check the edge eminating rightward from the start_vertex.
        e = p.edge(initial_vertex)
        if e[0]**2 + e[1]**2 <= squared_length_bound:
            sc_list.append( SaddleConnection(self, start_data, e) )

        # Represents the bounds of the beam of trajectories we are sending out.
        wedge = ( last_sim( p.vertex((initial_vertex+1)%p.num_edges()) ),
                  last_sim( p.vertex((initial_vertex+p.num_edges()-1)%p.num_edges()) ))

        # This will collect the data we need for a depth first search.
        chain = [(last_sim, initial_label, wedge, [(initial_vertex+p.num_edges()-i)%p.num_edges() for i in range(2,p.num_edges())])]

        while len(chain)>0:
            # Should verts really be edges?
            sim, label, wedge, verts = chain[-1]
            if len(verts) == 0:
                chain.pop()
                continue
            vert = verts.pop()
            #print("Inspecting "+str(vert))
            p = self.polygon(label)
            # First check the vertex
            vert_position = sim(p.vertex(vert))
            #print(wedge[1].n())
            if wedge_product(wedge[0], vert_position) > 0 and \
               wedge_product(vert_position, wedge[1]) > 0 and \
               vert_position[0]**2 + vert_position[1]**2 <= squared_length_bound:
                    sc_list.append( SaddleConnection(self, start_data, vert_position,
                                                   end_data = (label,vert),
                                                   end_direction = ~sim.derivative()*-vert_position,
                                                   holonomy = vert_position,
                                                   end_holonomy = ~sim.derivative()*-vert_position,
                                                   check = check) )
            # Now check if we should develop across the edge
            vert_position2 = sim(p.vertex( (vert+1)%p.num_edges() ))
            if wedge_product(vert_position,vert_position2)>0 and \
               wedge_product(wedge[0],vert_position2)>0 and \
               wedge_product(vert_position,wedge[1])>0 and \
               circle.line_segment_position(vert_position, vert_position2)==1:
                if wedge_product(wedge[0], vert_position) > 0:
                    # First in new_wedge should be vert_position
                    if wedge_product(vert_position2, wedge[1]) > 0:
                        new_wedge = (vert_position, vert_position2)
                    else:
                        new_wedge = (vert_position, wedge[1])
                else:
                    if wedge_product(vert_position2, wedge[1]) > 0:
                        new_wedge = (wedge[0], vert_position2)
                    else:
                        new_wedge=wedge
                new_label, new_edge = self.opposite_edge(label, vert)
                new_sim = sim*~self.edge_transformation(label,vert)
                p = self.polygon(new_label)
                chain.append( (new_sim, new_label, new_wedge, [(new_edge+p.num_edges()-i)%p.num_edges() for i in range(1,p.num_edges())]) )
        return sc_list

    def set_default_graphical_surface(self, graphical_surface):
        r"""
        Replace the default graphical surface with the provided GraphicalSurface.
        """
        from flatsurf.graphical.surface import GraphicalSurface
        if not isinstance(graphical_surface, GraphicalSurface):
            raise ValueError("graphical_surface must be a GraphicalSurface")
        if self != graphical_surface.get_surface():
            raise ValueError("The provided graphical_surface renders a different surface!")
        self._gs =  graphical_surface

    def graphical_surface(self, *args, **kwds):
        r"""
        Return a GraphicalSurface representing this surface.

        By default this returns a cached version of the GraphicalSurface. If
        ``cached=False`` is provided as a keyword option then a new
        GraphicalSurface is returned. Other keyword options:

        INPUT:

        - ``cached`` -- a boolean (default ``True``). If true return a cached
          GraphicalSurface. Otherwise we make a new one.

        - ``polygon_labels`` -- a boolean (default ``True``) whether the label
          of polygons are displayed

        - ``edge_labels`` -- option to control the display of edge labels. It
          can be one of

            - ``False`` or ``None`` for no labels

            - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
              name of the polygon to which it is glued

            - ``'number'`` -- to put on each side of each edge the number of the
              edge

            - ``'gluings and numbers'`` -- full information

            - ``'letter'`` -- add matching letters to glued edges in an arbitrary way

        - ``default_position_function`` -- a function mapping polygon labels to
          similarities describing the position of the corresponding polygon.

        EXAMPLES:

        Test the difference between the cached graphical_surface and the uncached version::

            sage: from flatsurf import *
            sage: s = translation_surfaces.octagon_and_squares()
            sage: s.plot()     # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 32 graphics primitives
            sage: s.graphical_surface(cached=False,adjacencies=[]).plot()   # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 18 graphics primitives
        """
        from flatsurf.graphical.surface import GraphicalSurface
        if "cached" in kwds:
            if not kwds["cached"]:
                # cached=False: return a new surface.
                kwds.pop("cached",None)
                return GraphicalSurface(self, *args, **kwds)
            kwds.pop("cached",None)
        if hasattr(self, '_gs'):
            self._gs.process_options(*args, **kwds)
        else:
            self._gs = GraphicalSurface(self, *args, **kwds)
        return self._gs

    def plot(self, *args, **kwds):
        r"""
        Returns a plot of the surface.

        There may be zero or one argument. If provided the single argument
        should be a GraphicalSurface whick will be used in the plot.

        INPUT:

        - ``polygon_labels`` -- a boolean (default ``True``) whether the label
          of polygons are displayed

        - ``edge_labels`` -- option to control the display of edge labels. It
          can be one of

            - ``False`` or ``None`` for no labels

            - ``'gluings'`` -- to put on each side of each non-adjacent edge, the
              name of the polygon to which it is glued

            - ``'number'`` -- to put on each side of each edge the number of the
              edge

            - ``'gluings and number'`` -- full information

        - ``adjacencies`` -- a list of pairs ``(p,e)`` to be used to set
          adjacencies of polygons.

        - ``default_position_function`` -- a function mapping polygon labels to
          similarities describing the position of the corresponding polygon.
        """
        if len(args) > 1:
            raise ValueError("SimilaritySurface.plot() can take at most one non-keyword argument.")
        if len(args)==1:
            from flatsurf.graphical.surface import GraphicalSurface
            if not isinstance(args[0], GraphicalSurface):
                raise ValueError("If an argument is provided, it must be a GraphicalSurface.")
            gs = args[0]
            gs.process_options(**kwds)
        else:
            gs = self.graphical_surface(**kwds)
        return gs.plot()

    def plot_polygon(self, label, graphical_surface = None,
                     plot_polygon = True, plot_edges = True, plot_edge_labels = True,
                     edge_labels = None,
                     polygon_options = {"axes":True}, edge_options = None, edge_label_options = None):
        r"""
        Returns a plot of the polygon with the provided label.

        Note that this method plots the polygon in its coordinates as opposed to
        graphical coordinates that the :func:``plot`` method uses. This makes it useful
        for visualizing the natural coordinates of the polygon.

        INPUT:

            - ``graphical_surface`` -- (default ``None``) If provided this function pulls graphical options
              from the graphical surface. If not provided, we use the default graphical surface.

            - ``plot_polygon`` -- (default ``True``) If True, we plot the solid polygon.

            - ``polygon_options`` -- (default ``{"axes":True}``) Options for the rendering of the polygon.
              These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_polygon`.
              This should be either None or a dictionary.

            - ``plot_edges`` -- (default ``True``) If True, we plot the edges of the polygon as segments.

            - ``edge_options`` -- (default ``None``) Options for the rendering of the polygon edges.
              These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_edge`.
              This should be either None or a dictionary.

            - ``plot_edge_labels`` -- (default ``True``) If True, we plot labels on the edges.

            - ``edge_label_options`` -- (default ``None``) Options for the rendering of the edge labels.
              These options will be passed to :func:`~flatsurf.graphical.polygon.GraphicalPolygon.plot_edge_label`.
              This should be either None or a dictionary.

            - ``edge_labels`` -- (default ``None``) If None and plot_edge_labels is True, we write the edge
              number on each edge. Otherwise edge_labels should be a list of strings of length equal to the
              number of edges of the polygon. The strings will be printed on each edge.

        EXAMPLES::

            sage: from flatsurf import *
            sage: s = similarity_surfaces.example()
            sage: s.plot() # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 13 graphics primitives
            s.plot_polygon(1) # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 7 graphics primitives

            sage: labels = []
            sage: p = s.polygon(1)
            sage: for e in range(p.num_edges()): \
                labels.append(str(p.edge(e)))
            sage: s.plot_polygon(1, polygon_options=None, plot_edges=False, \
                edge_labels=labels, edge_label_options={"color":"red"}) # not tested (problem with matplotlib font caches on Travis)
            Graphics object consisting of 4 graphics primitives
        """
        if graphical_surface is None:
            graphical_surface = self.graphical_surface()
        p = self.polygon(label)
        from flatsurf.graphical.polygon import GraphicalPolygon
        gp = GraphicalPolygon(p)

        if plot_polygon:
            if polygon_options is None:
                o = graphical_surface.polygon_options
            else:
                o = graphical_surface.polygon_options.copy()
                o.update(polygon_options)
            plt = gp.plot_polygon(**o)

        if plot_edges:
            if edge_options is None:
                o = graphical_surface.non_adjacent_edge_options
            else:
                o = graphical_surface.non_adjacent_edge_options.copy()
                o.update(edge_options)
            for e in range(p.num_edges()):
                plt += gp.plot_edge(e, **o)

        if plot_edge_labels:
            if edge_label_options is None:
                o = graphical_surface.edge_label_options
            else:
                o = graphical_surface.edge_label_options.copy()
                o.update(edge_label_options)
            for e in range(p.num_edges()):
                if edge_labels is None:
                    el = str(e)
                else:
                    el = edge_labels[e]
                plt += gp.plot_edge_label(e, el, **o)
        return plt

# I'm not sure we want to support this...
#
#    def minimize_monodromy_mapping(self):
#        r"""
#        Return a mapping from this surface to a similarity surface
#        with a minimal monodromy group.
#        Note that this may be slow for infinite surfaces.
#
#        EXAMPLES::
#            sage: from flatsurf.geometry.polygon import ConvexPolygons
#            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
#            sage: octagon = ConvexPolygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
#            sage: square = ConvexPolygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
#            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
#            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
#            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
#            sage: print s
#            Rational cone surface built from 3 polygons
#            sage: m=s.minimize_monodromy_mapping()
#            sage: s2=m.codomain()
#            sage: print s2
#            Translation surface built from 3 polygons
#            sage: v=s.tangent_vector(2,(0,0),(1,0))
#            sage: print m.push_vector_forward(v)
#            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (-1/2*sqrt2, -1/2*sqrt2)
#            sage: w=s2.tangent_vector(2,(0,0),(0,-1))
#            sage: print m.pull_vector_back(w)
#            SimilaritySurfaceTangentVector in polygon 2 based at (0, 0) with vector (1/2*sqrt2, 1/2*sqrt2)
#        """
#        lw = self.walker()
#        class MatrixFunction:
#            def __init__(self, lw):
#                self._lw=lw
#                from sage.matrix.constructor import identity_matrix
#                self._d = {lw.surface().base_label():
#                    identity_matrix(lw.surface().base_ring(), n=2)}
#            def __call__(self, label):
#                try:
#                    return self._d[label]
#                except KeyError:
#                    e = self._lw.edge_back(label)
#                    label2,e2 = self._lw.surface().opposite_edge(label,e)
#                    m=self._lw.surface().edge_matrix(label,e) * self(label2)
#                    self._d[label]=m
#                    return m
#        mf = MatrixFunction(lw)
#        from flatsurf.geometry.mappings import (
#            MatrixListDeformedSurfaceMapping,
#            IdentityMapping)
#        mapping = MatrixListDeformedSurfaceMapping(self, mf)
#        surface_type = mapping.codomain().compute_surface_type_from_gluings(limit=100)
#        new_codomain = convert_to_type(mapping.codomain(),surface_type)
#        identity = IdentityMapping(mapping.codomain(), new_codomain)
#        return identity * mapping
#
#    def minimal_monodromy_surface(self):
#        r"""
#        Return an equivalent similarity surface with minimal monodromy.
#        Note that this may be slow for infinite surfaces.
#
#        EXAMPLES::
#            sage: from flatsurf.geometry.polygon import ConvexPolygons
#            sage: K.<sqrt2> = NumberField(x**2 - 2, embedding=1.414)
#            sage: octagon = ConvexPolygons(K)([(1,0),(sqrt2/2, sqrt2/2),(0, 1),(-sqrt2/2, sqrt2/2),(-1,0),(-sqrt2/2, -sqrt2/2),(0, -1),(sqrt2/2, -sqrt2/2)])
#            sage: square = ConvexPolygons(K)([(1,0),(0,1),(-1,0),(0,-1)])
#            sage: gluings = [((0,i),(1+(i%2),i//2)) for i in range(8)]
#            sage: from flatsurf.geometry.surface import surface_from_polygons_and_gluings
#            sage: s=surface_from_polygons_and_gluings([octagon,square,square],gluings)
#            sage: print s
#            Rational cone surface built from 3 polygons
#            sage: s2=s.minimal_monodromy_surface()
#            sage: print s2
#            Translation surface built from 3 polygons
#        """
#        return self.minimize_monodromy_mapping().codomain()

    def __eq__(self, other):
        r"""
        Implements a naive notion of equality where two finite surfaces are equal if:
        - their base labels are equal,
        - their polygons are equal and labeled and glued in the same way.
        For infinite surfaces we use reference equality.
        Raises a value error if the surfaces are defined over different rings.
        """
        if not self.is_finite():
            return self is other
        if self is other:
            return True
        if not isinstance(other, SimilaritySurface):
            raise TypeError
        if not other.is_finite():
            raise ValueError("Can not compare infinite surfaces.")
        if self.base_ring() != other.base_ring():
            raise ValueError("Refusing to compare surfaces with different base rings.")
        if not self.is_mutable() and not other.is_mutable():
            hash1 = hash(self)
            hash2 = hash(other)
            if hash1 != hash2:
                return False
        if self.base_label() != other.base_label():
            return False
        if self.num_polygons() != other.num_polygons():
            return False
        for label,polygon in self.label_iterator(polygons=True):
            try:
                polygon2 = other.polygon(label)
            except ValueError:
                return False
            if polygon != polygon2:
                return False
            for edge in range(polygon.num_edges()):
                if self.opposite_edge(label,edge) != other.opposite_edge(label,edge):
                    return False
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        r"""
        Hash compatible with equals.
        """
        if self._s.is_mutable():
            raise ValueError("Attempting to hash with mutable underlying surface.")
        if hasattr(self, '_hash'):
            # Return the cached hash.
            return self._hash
        # Compute the hash
        h = 17*hash(self.base_ring())+23*hash(self.base_label())
        for pair in self.label_iterator(polygons=True):
            h = h + 7*hash(pair)
        for edgepair in self.edge_iterator(gluings=True):
            h = h + 3*hash(edgepair)
        self._hash=h
        return h

    def erase_marked_points(self):
        r"""
        Return an isometric or similar surface with a minimal number of regular
        vertices of angle 2π.

        EXAMPLES::

            sage: import flatsurf

            sage: G = SymmetricGroup(4)
            sage: S = flatsurf.translation_surfaces.origami(G('(1,2,3,4)'), G('(1,4,2,3)'))
            sage: S.stratum()
            H_2(2, 0)
            sage: S.erase_marked_points().stratum() # optional: pyflatsurf
            H_2(2)

            sage: for (a,b,c) in [(1,4,11), (1,4,15), (3,4,13)]: # optional: pyflatsurf
            ....:     T = flatsurf.polygons.triangle(a,b,c)
            ....:     S = flatsurf.similarity_surfaces.billiard(T)
            ....:     S = S.minimal_cover("translation")
            ....:     print(S.erase_marked_points().stratum())
            H_6(10)
            H_6(2^5)
            H_8(12, 2)

        TESTS:

        Verify that https://github.com/flatsurf/flatsurf/issues/263 has been resolved::

            sage: from flatsurf import EquiangularPolygons, similarity_surfaces
            sage: E = EquiangularPolygons((10, 8, 3, 1, 1, 1))
            sage: P = E((1, 1, 2, 4), normalized=True)
            sage: B = similarity_surfaces.billiard(P, rational=True)
            sage: S = B.minimal_cover(cover_type="translation")
            sage: S = S.erase_marked_points() # optional: pyflatsurf

        ::

            sage: from flatsurf import EquiangularPolygons, similarity_surfaces
            sage: E = EquiangularPolygons((10, 7, 2, 2, 2, 1))
            sage: P = E((1, 1, 2, 3), normalized=True)
            sage: B = similarity_surfaces.billiard(P, rational=True)
            sage: S_mp = B.minimal_cover(cover_type="translation")
            sage: S = S_mp.erase_marked_points() # optional: pyflatsurf

        """
        from .pyflatsurf_conversion import from_pyflatsurf, to_pyflatsurf
        S = to_pyflatsurf(self)
        S = S.eliminateMarkedPoints().surface()
        S.delaunay()
        return from_pyflatsurf(S)
