r"""Mappings between translation surfaces."""

class SimilaritySurfaceMapping:
    r"""Abstract class for any mapping between similarity surfaces."""
    
    def __init__(self, domain, codomain):
        self._domain=domain
        self._codomain=codomain
    
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

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        raise NotImplementedError

    def pull_vector_back(self,tangent_vector):
        r"""Applies the inverse of the mapping to the provided vector."""
        raise NotImplementedError

class SimilarityJoinPolygonsMapping(SimilaritySurfaceMapping):
    def __init__(self, s, p1, e1):
        poly1=s.polygon(p1)
        p2,e2 = s.opposite_edge(p1,e1)
        poly2=s.polygon(p2)
        t=s.edge_transformation(p2, e2)
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
        
        inv_edge_map={}
        for key, value in edge_map.iteritems():
            inv_edge_map[value]=(p1,key)
        
        from flatsurf.geometry.polygon import Polygons
        newpoly = Polygons(s.base_ring())(vs)
        
        if s.is_finite():
            # Finite surface. Use a dictionary for polygons
            polygons={}
            for label in s.polygon_labels():
                if label==p1:
                    polygons[label]=newpoly
                elif label != p2:
                    polygons[label]=s.polygon(label)
            gluings=[]
            # Iterate over the new edges:
            for ll1,pp1 in polygons.iteritems():
                for ee1 in range(pp1.num_edges()):
                    if ll1 == p1:
                        ll2,ee2 = edge_map[ee1]
                    else:
                        ll2,ee2 = (ll1,ee1)
                    ll3,ee3=s.opposite_edge(ll2,ee2)
                    if ll3 == p1 or ll3 == p2:
                        gluings.append( ( (ll1,ee1), inv_edge_map[(ll3,ee3)] ) )
                    else:
                        gluings.append( ( (ll1,ee1), (ll3,ee3) ) )
            print gluings
            from flatsurf.geometry.similarity_surface import SimilaritySurface_polygons_and_gluings
            s2=SimilaritySurface_polygons_and_gluings(polygons,gluings)
        else:
            # There is no reason this could not be implemented.
            raise NotImplementedError
        self._saved_label=p1
        self._removed_label=p2
        self._remove_map = t
        self._remove_map_derivative = dt
        self._glued_edge=e1
        SimilaritySurfaceMapping.__init__(self, s, s2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._removed_label:
            return self._codomain.tangent_vector( \
                self._saved_label, \
                self._remove_map(tangent_vector.point()), \
                self._remove_map_derivative*tangent_vector.vector(), \
                ring = ring)
        else:
            return self._codomain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""
        Applies the inverse of the mapping to the provided vector.
        """
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._saved_label:
            p=tangent_vector.point()
            v=self._domain.polygon(self._saved_label).vertex(self._glued_edge)
            e=self._domain.polygon(self._saved_label).edge(self._glued_edge)
            from flatsurf.geometry.polygon import wedge_product
            wp = wedge_product(p-v,e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._removed_label, \
                    (~ self._remove_map)(tangent_vector.point()), \
                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
                    ring = ring)
            if wp < 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._saved_label, \
                    tangent_vector.point(), \
                    tangent_vector.vector(), \
                    ring = ring)
            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = wedge_product(w,e)
            if wp > 0:
                # in polygon with the removed label
                return self.domain().tangent_vector( \
                    self._removed_label, \
                    (~ self._remove_map)(tangent_vector.point()), \
                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
                    ring = ring)
            return self.domain().tangent_vector( \
                self._saved_label, \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            return self._domain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

#class AffineMappingJoin(self):
#    def __init__(self, s, p1, e1, p2, e2):
        
