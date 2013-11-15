class SurfacePoint:

    def __init__(self, similarity_surface, polygon_label, point):
        self._ss=similarity_surface
        self._l=polygon_label
        self._p=point
        # check that point is within the polygon
        p=self._ss.polygon(self._l)
        if not p.contains_point(self._p):
            raise ValueError("point not within specified polygon")

    def __repr__(self):
        return "point within polygon "+str(self._l)+" located at "+str(self._p)

    def _flow_step(self,holonomy):
        p=self._ss.polygon(self._l)
        q,hol,pos=p.flow(self._p,holonomy)
        if pos.is_in_edge_interior():
            # This is annoyingly slow because of vertices(). Maybe we should cache translations too...
            e=pos.get_edge()
            ll,ee=self._ss.opposite_edge(self._l,e)
            pp=self._ss.polygon(ll)
            v0=p.vertices()[e]
            v1=pp.vertices()[(ee+1)%pp.num_edges()]
            m=self._ss.edge_matrix(self._l,e)
            p_res=SurfacePoint(self._ss,ll,m*(q-v0)+v1)
            hol_res=m*hol
            return q,p_res,hol_res,pos
        else:
            # Either hit a vertex, or still within interior of polygon
            return q,SurfacePoint(self._ss,self._l,q),hol,pos

    def flow_segments(self,holonomy):
        pair_list=[]
        hol=holonomy
        V=self._ss.vector_space()
        pt0=self
        while hol != V.zero():
            q,pt1,hol,pos=pt0._flow_step(hol)
            if pos.is_vertex() and hol != V.zero():
                raise ValueError("Straight line flow hit a vertex")
            pair_list.append((pt0,q))
            pt0=pt1
        print pt0
        return pair_list

    def flow_segments_wait(self,holonomy):
        pair_list=self.flow_segments(holonomy)
        pair_list2=[]
        for i in range(len(pair_list)/2,len(pair_list)):
            pair_list2.append(pair_list[i])
        return pair_list2

    def flow(self,holonomy):
        hol=holonomy
        V=self._ss.vector_space()
        pt=self
        while hol != V.zero():
            q,pt,hol,pos=pt._flow_step(hol)
            if pos.is_vertex() and hol != V.zero():
                raise ValueError("Straight line flow hit a vertex")
        return pt

    def get_label(self):
        return self._l

    def get_point(self):
        return self._p
