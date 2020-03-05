from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.qqbar import AA, number_field_elements_from_algebraics
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.functions.other import sqrt

from flatsurf.geometry.polygon import ConvexPolygons
from flatsurf.geometry.surface import surface_list_from_polygons_and_gluings
from flatsurf.geometry.cone_surface import ConeSurface
from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectory, SegmentInPolygon 
from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector

class ConeSurfaceToPolyhedronMap(SageObject):
    r"""
    A map sending objects defined on a ConeSurface built from a polyhedron to the polyhedron. Currently, this 
    works to send a trajectory to a list of points.
    
    This class should not be called directly. You get an object of this type from polyhedron_to_cone_surface.
    """
    def __init__(self, cone_surface, polyhedron, mapping_data):
        self._s=cone_surface
        self._p=polyhedron
        self._md=mapping_data
    
    def __call__(self,o):
        r"""
        This method is used to convert from an object on the cone surface to an object on the polyhedron.
        
        Currently works with 
        - StraightLineTrajectory -- returns the corresponding list of points on the polyhedron
        - SegmentInPolygon -- returns the corresponding pair of points on the polyhedron
        - SimilaritySurfaceTangentVector -- returns a pair of points corresponding to the image point and image of the tangent vector.
        """
        if isinstance(o, StraightLineTrajectory):
            points=[]
            it = iter(o.segments())
            s = next(it)
            label = s.polygon_label()
            points.append(self._md[label][0]+self._md[label][1]*s.start().point())
            points.append(self._md[label][0]+self._md[label][1]*s.end().point())
            for s in it:
                label = s.polygon_label()
                points.append(self._md[label][0]+self._md[label][1]*s.end().point())
            return points    
        if isinstance(o, SegmentInPolygon):
            # Return the pair of images of the endpoints.
            label = o.polygon_label()
            return ( self._md[label][0]+self._md[label][1]*o.start().point(), \
                     self._md[label][0]+self._md[label][1]*o.end().point() )
        if isinstance(o,SimilaritySurfaceTangentVector):
            # Map to a pair of vectors conisting of the image of the basepoint and the image of the vector.
            label = o.polygon_label()
            point = o.point()
            vector = o.vector()
            return (self._md[label][0]+self._md[label][1]*point, \
                    self._md[label][1]*vector)
        raise ValueError("Failed to recognize type of passed object")
        

def polyhedron_to_cone_surface(polyhedron, use_AA=False, scaling_factor=ZZ(1)):
    r"""Construct the Euclidean Cone Surface associated to the surface of a polyhedron and a map
    from the cone surface to the polyhedron.
    
    INPUT:

    - ``polyhedron`` -- A 3-dimensional polyhedron, which should be define over something that coerces into AA

    - ``use_AA`` -- If True, the surface returned will be defined over AA. If false, the algorithm will find the smallest NumberField and write the field there.
    
    - ``scaling_factor`` -- The surface returned will have a metric scaled by multiplication by this factor (compared with the original polyhendron). This can be used to produce a surface defined over a smaller NumberField.
    
    OUTPUT:
    
    A pair consisting of a ConeSurface and a ConeSurfaceToPolyhedronMap.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import *
    sage: vertices=[]
    sage: for i in range(3):
    ....:     temp=vector([1 if k==i else 0 for k in range(3)])
    ....:     for j in range(-1,3,2):
    ....:         vertices.append(j*temp)
    sage: octahedron=Polyhedron(vertices=vertices)
    sage: surface,surface_to_octahedron = \
    ....:     polyhedron_to_cone_surface(octahedron,scaling_factor=AA(1/sqrt(2)))
    sage: TestSuite(surface).run()
    sage: TestSuite(surface_to_octahedron).run(skip="_test_pickling")
    sage: surface.num_polygons()
    8
    sage: surface.base_ring()
    Number Field in a with defining polynomial y^2 - 3 with a = 1.732050807568878?
    sage: sqrt3=surface.base_ring().gen()
    sage: tangent_bundle=surface.tangent_bundle()
    sage: v=tangent_bundle(0,(0,0),(sqrt3,2))
    sage: traj=v.straight_line_trajectory()
    sage: traj.flow(10)
    sage: traj.is_saddle_connection()
    True
    sage: traj.combinatorial_length()
    8
    sage: path3d = surface_to_octahedron(traj)
    sage: len(path3d)
    9
    sage: # We will show that the length of the path is sqrt(42):
    sage: total_length = 0
    sage: for i in range(8):
    ....:     start = path3d[i]
    ....:     end = path3d[i+1]
    ....:     total_length += (vector(end)-vector(start)).norm()
    sage: ZZ(total_length**2)
    42
    """
    assert polyhedron.dim()==3
    c=polyhedron.center()
    vertices=polyhedron.vertices()
    vertex_order={}
    for i,v in enumerate(vertices):
        vertex_order[v]=i
    faces=polyhedron.faces(2)
    face_order={}
    face_edges=[]
    face_vertices=[]
    face_map_data=[]
    for i,f in enumerate(faces):
        face_order[f]=i
        edges=f.as_polyhedron().faces(1)
        face_edges_temp = set()
        for edge in edges:
            edge_temp=set()
            for vertex in edge.vertices():
                v=vertex.vector()
                v.set_immutable()
                edge_temp.add(v)
            face_edges_temp.add(frozenset(edge_temp))

            
        last_edge = next(iter(face_edges_temp))
        v = next(iter(last_edge))
        face_vertices_temp=[v]
        for j in range(len(face_edges_temp)-1):
            for edge in face_edges_temp:
                if v in edge and edge!=last_edge:
                    # bingo
                    last_edge=edge
                    for vv in edge:
                        if vv!=v:
                            v=vv
                            face_vertices_temp.append(vv)
                            break
                    break


        v0=face_vertices_temp[0]
        v1=face_vertices_temp[1]
        v2=face_vertices_temp[2]
        n=(v1-v0).cross_product(v2-v0)
        if (v0-c).dot_product(n)<0:
            n=-n
            face_vertices_temp.reverse()
            v0=face_vertices_temp[0]
            v1=face_vertices_temp[1]
            v2=face_vertices_temp[2]

        face_vertices.append(face_vertices_temp)    
        n=n/AA(n.norm())
        w=v1-v0
        w=w/AA(w.norm())
        m=1/scaling_factor*matrix(AA,[w,n.cross_product(w),n]).transpose()
        mi=~m
        mis=mi.submatrix(0,0,2,3)
        face_map_data.append((
            v0, # translation to bring origin in plane to v0
            m.submatrix(0,0,3,2),
            -mis*v0,
            mis
        ))

        it=iter(face_vertices_temp)    
        v_last=next(it)
        face_edge_dict={}
        j=0
        for v in it:
            edge=frozenset([v_last,v])
            face_edge_dict[edge]=j
            j+=1
            v_last=v
        v=next(iter(face_vertices_temp))
        edge=frozenset([v_last,v])
        face_edge_dict[edge]=j
        face_edges.append(face_edge_dict)
        
    gluings={}
    for p1,face_edge_dict1 in enumerate(face_edges):
        for edge, e1 in face_edge_dict1.items():
            found=False
            for p2, face_edge_dict2 in enumerate(face_edges):
                if p1!=p2 and edge in face_edge_dict2:
                    e2=face_edge_dict2[edge]
                    gluings[(p1,e1)]=(p2,e2)
                    found=True
                    break
            if not found:
                print(p1)
                print(e1)
                print(edge)
                raise RuntimeError("Failed to find glued edge")
    polygon_vertices_AA=[]
    for p,vs in enumerate(face_vertices):
        trans=face_map_data[p][2]
        m=face_map_data[p][3]
        polygon_vertices_AA.append([trans+m*v for v in vs])
        
    
    if use_AA==True:
        Polys=ConvexPolygons(AA)
        polygons=[]
        for vs in polygon_vertices_AA:
            polygons.append(Polys(vertices=vs))
        S=ConeSurface(surface_list_from_polygons_and_gluings(polygons,gluings))
        return S, \
            ConeSurfaceToPolyhedronMap(S,polyhedron,face_map_data)
    else:
        elts=[]
        for vs in polygon_vertices_AA:
            for v in vs:
                elts.append(v[0])
                elts.append(v[1])
                
        # Find the best number field:
        field,elts2,hom = number_field_elements_from_algebraics(elts,minimal=True)
        if field==QQ:
            # Defined over the rationals!
            polygon_vertices_field2=[]
            j=0
            for vs in polygon_vertices_AA:
                vs2=[]
                for v in vs:
                    vs2.append(vector(field,[elts2[j],elts2[j+1]]))
                    j=j+2
                polygon_vertices_field2.append(vs2)
            Polys=ConvexPolygons(field)
            polygons=[]
            for vs in polygon_vertices_field2:
                polygons.append(Polys(vertices=vs))
            S=ConeSurface(surface_list_from_polygons_and_gluings(polygons,gluings))
            return S, \
                ConeSurfaceToPolyhedronMap(S,polyhedron,face_map_data)

        else:        
            # Unfortunately field doesn't come with an real embedding (which is given by hom!)
            # So, we make a copy of the field, and add the embedding.
            field2=NumberField(field.polynomial(),name="a",embedding=hom(field.gen()))
            # The following converts from field to field2:
            hom2=field.hom(im_gens=[field2.gen()])

            polygon_vertices_field2=[]
            j=0
            for vs in polygon_vertices_AA:
                vs2=[]
                for v in vs:
                    vs2.append(vector(field2,[hom2(elts2[j]),hom2(elts2[j+1])]))
                    j=j+2
                polygon_vertices_field2.append(vs2)
            Polys=ConvexPolygons(field2)
            polygons=[]
            for vs in polygon_vertices_field2:
                polygons.append(Polys(vertices=vs))
            S=ConeSurface(surface_list_from_polygons_and_gluings(polygons,gluings))
            return S, \
                ConeSurfaceToPolyhedronMap(S,polyhedron,face_map_data)

def platonic_tetrahedron():
    r"""Produce a triple consisting of a polyhedral version of the platonic tetrahedron,
    the associated cone surface, and a ConeSurfaceToPolyhedronMap from the surface
    to the polyhedron.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import platonic_tetrahedron
    sage: polyhedron,surface,surface_to_polyhedron = platonic_tetrahedron()
    sage: TestSuite(surface).run()
    r"""
    vertices=[]
    for x in range(-1,3,2):
        for y in range(-1,3,2):
                vertices.append(vector(QQ,(x,y,x*y)))
    p=Polyhedron(vertices=vertices)
    s,m = polyhedron_to_cone_surface(p,scaling_factor=AA(1/sqrt(2)))
    return p,s,m

def platonic_cube():
    r"""Produce a triple consisting of a polyhedral version of the platonic cube,
    the associated cone surface, and a ConeSurfaceToPolyhedronMap from the surface
    to the polyhedron.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import platonic_cube
    sage: polyhedron,surface,surface_to_polyhedron = platonic_cube()
    sage: TestSuite(surface).run()
    r"""
    vertices=[]
    for x in range(-1,3,2):
        for y in range(-1,3,2):
            for z in range(-1,3,2):
                vertices.append(vector(QQ,(x,y,z)))
    p=Polyhedron(vertices=vertices)
    s,m = polyhedron_to_cone_surface(p,scaling_factor=QQ(1)/2)
    return p,s,m

def platonic_octahedron():
    r"""Produce a triple consisting of a polyhedral version of the platonic octahedron,
    the associated cone surface, and a ConeSurfaceToPolyhedronMap from the surface
    to the polyhedron.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import platonic_octahedron
    sage: polyhedron,surface,surface_to_polyhedron = platonic_octahedron()
    sage: TestSuite(surface).run()
    r"""
    vertices=[]
    for i in range(3):
        temp=vector(QQ,[1 if k==i else 0 for k in range(3)])
        for j in range(-1,3,2):
            vertices.append(j*temp)
    octahedron=Polyhedron(vertices=vertices)
    surface,surface_to_octahedron = \
        polyhedron_to_cone_surface(octahedron,scaling_factor=AA(sqrt(2)))
    return octahedron,surface,surface_to_octahedron

def platonic_dodecahedron():
    r"""Produce a triple consisting of a polyhedral version of the platonic dodecahedron,
    the associated cone surface, and a ConeSurfaceToPolyhedronMap from the surface
    to the polyhedron.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import platonic_dodecahedron
    sage: polyhedron,surface,surface_to_polyhedron = platonic_dodecahedron()
    sage: TestSuite(surface).run()
    r"""
    vertices=[]
    phi=AA(1+sqrt(5))/2
    F=NumberField(phi.minpoly(),"phi",embedding=phi)
    phi=F.gen()
    for x in range(-1,3,2):
        for y in range(-1,3,2):
            for z in range(-1,3,2):
                vertices.append(vector(F,(x,y,z)))
    for x in range(-1,3,2):
        for y in range(-1,3,2):
            vertices.append(vector(F,(0,x*phi,y/phi)))
            vertices.append(vector(F,(y/phi,0,x*phi)))
            vertices.append(vector(F,(x*phi,y/phi,0)))
    scale=AA(2/sqrt(1+(phi-1)**2+(1/phi-1)**2))
    p=Polyhedron(vertices=vertices)
    s,m = polyhedron_to_cone_surface(p,scaling_factor=scale)
    return p,s,m

def platonic_icosahedron():
    r"""Produce a triple consisting of a polyhedral version of the platonic icosahedron,
    the associated cone surface, and a ConeSurfaceToPolyhedronMap from the surface
    to the polyhedron.

    EXAMPLES::

    sage: from flatsurf.geometry.polyhedra import platonic_icosahedron
    sage: polyhedron,surface,surface_to_polyhedron = platonic_icosahedron()
    sage: TestSuite(surface).run()
    r"""
    vertices=[]
    phi=AA(1+sqrt(5))/2
    F=NumberField(phi.minpoly(),"phi",embedding=phi)
    phi=F.gen()
    for i in range(3):
        for s1 in range(-1,3,2):
            for s2 in range(-1,3,2):
                p=3*[None]
                p[i]=s1*phi
                p[(i+1)%3]=s2
                p[(i+2)%3]=0
                vertices.append(vector(F,p))
    p=Polyhedron(vertices=vertices)
    
    s,m = polyhedron_to_cone_surface(p)
    return p,s,m
