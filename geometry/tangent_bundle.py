from geometry.polygon import *

class SimilaritySurfaceTangentVector:
    def __init__(self, tangent_bundle, polygon_label, point, vector):
        self._bundle=tangent_bundle
        p=self.surface().polygon(polygon_label)
        pos=p.get_point_position(point)
        if vector == self._bundle.vector_space().zero():
            raise ValueError("Provided vector is zero. (Temporarily not supported.)")
        if pos.is_in_interior():
            self._polygon_label=polygon_label
            self._point=point
            self._vector=vector
            self._position=pos
        elif pos.is_in_edge_interior():
            e=pos.get_edge()
            edge_v=p.edge(e)
            if wedge_product(edge_v,vector)<0 or is_opposite_direction(edge_v,vector):
                # Need to move point and vector to opposite edge.
                label2,e2 = self.surface().opposite_edge(polygon_label,e)
                similarity = self.surface().edge_transformation(polygon_label,e)
                point2=similarity(point)
                vector2=similarity.derivative()*vector
                self._polygon_label=label2
                self._point=point2
                self._vector=vector2
                self._position=self.surface().polygon(label2).get_point_position(point2)
            else:
                self._polygon_label=polygon_label
                self._point=point
                self._vector=vector
                self._position=pos
        elif pos.is_vertex():
            v=pos.get_vertex()
            # subsequent edge:
            edge1 = self.surface().polygon(polygon_label).edge(v)
            # prior edge:
            edge0 = self.surface().polygon(polygon_label).edge(v-1)
            wp1 = wedge_product(edge1,vector)
            wp0 = wedge_product(edge0,vector)
            if wp1<0 or wp0<0:
                raise ValueError("Singular point with vector pointing away from polygon")
            if wp0 == 0:
                # vector points backward along edge 0
                label2,e2 = self.surface().opposite_edge(polygon_label,v-1)
                similarity = self.surface().edge_transformation(polygon_label,v-1)
                point2=similarity(point)
                vector2=similarity.derivative()*vector
                self._polygon_label=label2
                self._point=point2
                self._vector=vector2
                self._position=self.surface().polygon(label2).get_point_position(point2)
            else:
                # vector points along edge1 in that directior or points into polygons interior
                self._polygon_label=polygon_label
                self._point=point
                self._vector=vector
                self._position=pos
        else:
            raise ValueError("Provided point lies outside the indexed polygon")
    
    def __repr__(self):
        return "SimilaritySurfaceTangentVector in polygon "+repr(self._polygon_label)+\
            " based at "+repr(self._point)+" with vector "+repr(self._vector)
    
    def surface(self):
        r"""Return the underlying surface."""
        return self._bundle.surface()

    def is_based_at_singularity(self):
        r"""
        Return the truth value of the statement 'the base point for this vector is a singularity.'
        """
        return self._position.is_vertex()
    
    def is_in_boundary_of_polygon(self):
        r"""
        Return the truth value of the statement 
        'the base point for this vector lies on the boundary of 
        one of the polygons making up the surface.'
        """
        return self._position.is_in_boundary()
    
    def bundle(self):
        r""" Return the tangent bundle containing this vector. """
        return self._bundle
    
    def polygon_label(self):
        return self._polygon_label

    def polygon(self):
        return self.surface().polygon(self.polygon_label())
        
    def point(self):
        r""" Return the coordinates of the basepoint of the vector within the assigned polygon. """
        return self._point
    
    def vector(self):
        r""" Return the coordinates of this vector within the assigned polygon. """
        return self._vector
        
    def differs_by_scaling(self, another_tangent_vector):
        r"""
        Returns true if the other vector just differs by scaling. This means they should lie
        in the same polygon, be based at the same point, and point in the same direction.
        """
        return self.polygon_label()==another_tangent_vector.polygon_label() and \
            self.point()==another_tangent_vector.point() and \
            is_same_direction(self.vector(),another_tangent_vector.vector())
        
    def invert(self):
        r"""
        Returns the negation of this tangent vector. 
        Raises a ValueError if the vector is based at a singularity.'
        """
        if self.is_based_at_singularity():
            raise ValueError("Can't invert tangent vector based at a singularity.")
        return SimilaritySurfaceTangentVector(
            self.bundle(), 
            self.polygon_label(), 
            self.point(),
            -self.vector())

    def forward_to_polygon_boundary(self):
        r"""
        Flows forward (in the direction of the tangent vector) until the end
        of the polygon is reached.
        Returns the tangent vector based at the endpoint which point backward along the trajectory.
        
        NOTES:: 
        
            We return the backward trajectory, because continuing forward does not make sense if a
            singularity is reached. You can obtain the forward vector by subsequently applying invert().
        
        EXAMPLES::
        
            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print("Polygon 0 is "+str(s.polygon(0)))
            Polygon 0 is Polygon: (0, 0), (2, -2), (2, 0)
            sage: print("Polygon 1 is "+str(s.polygon(1)))
            Polygon 1 is Polygon: (0, 0), (2, 0), (1, 3)
            sage: from geometry.tangent_bundle import SimilaritySurfaceTangentVector
            sage: V=tb.surface().vector_space()
            sage: v=SimilaritySurfaceTangentVector(tb, 0, V((0,0)), V((3,-1)))
            sage: print(v)
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (3, -1)
            sage: v2=v.forward_to_polygon_boundary()
            sage: print(v2)
            SimilaritySurfaceTangentVector in polygon 0 based at (2, -2/3) with vector (-3, 1)
            sage: print(v2.invert())
            SimilaritySurfaceTangentVector in polygon 1 based at (2/3, 2) with vector (4, -3)
        """
        p=self.polygon()
        point2,pos2 = p.flow_to_exit(self.point(), self.vector())
        #diff=point2-point
        new_vector = SimilaritySurfaceTangentVector(
            self.bundle(), 
            self.polygon_label(), 
            point2,
            -self.vector())
        return new_vector
    
class SimilaritySurfaceTangentBundle:
    def __init__(self, similarity_surface):
        self._s=similarity_surface

    def base_ring(self):
        return self._s.base_ring()
    
    field=base_ring

    def vector_space(self):
        r"""
        Return the vector space over the field of the bundle.
        """
        return self._s.vector_space()

    def surface(self):
        r"""Return the surface this bundle is over."""
        return self._s
        
    def edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels 
        counterclockwise around the boundary of the polygon along the edge with the provided index. 
        The length of the vector matches the length of the indexed edge.

        EXAMPLES::

            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print(s.polygon(0))
            Polygon: (0, 0), (2, -2), (2, 0)
            sage: print(tb.edge(0,0))
            SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (2, -2)
        """
        polygon=self.surface().polygon(polygon_label)
        point=polygon.vertex(edge_index)
        vector=polygon.edge(edge_index)
        return SimilaritySurfaceTangentVector(self, polygon_label, point, vector)

    def clockwise_edge(self, polygon_label, edge_index):
        r"""Return the vector leaving a vertex of the polygon which under straight-line flow travels 
        *clockwise* around the boundary of the polygon along the edge with the provided index. 
        The length of the vector matches the length of the indexed edge.
        Note that the point will be based in the polgon opposite the provided edge.
        
        EXAMPLES::

            sage: from geometry.similarity_surface_generators import SimilaritySurfaceGenerators
            sage: s=SimilaritySurfaceGenerators.example()
            sage: from geometry.tangent_bundle import SimilaritySurfaceTangentBundle
            sage: tb = SimilaritySurfaceTangentBundle(s)
            sage: print("Polygon 0 is "+str(s.polygon(0)))
            Polygon 0 is Polygon: (0, 0), (2, -2), (2, 0)
            sage: print("Polygon 1 is "+str(s.polygon(1)))
            Polygon 1 is Polygon: (0, 0), (2, 0), (1, 3)
            sage: print("Opposite edge to (0,0) is "+repr(s.opposite_edge(0,0)))
            Opposite edge to (0,0) is (1, 1)
            sage: print(tb.clockwise_edge(0,0))
            SimilaritySurfaceTangentVector in polygon 1 based at (2, 0) with vector (-1, 3)
        """
        polygon=self.surface().polygon(polygon_label)
        point=polygon.vertex(edge_index+1)
        vector=-polygon.edge(edge_index)
        return SimilaritySurfaceTangentVector(self, polygon_label, point, vector)

