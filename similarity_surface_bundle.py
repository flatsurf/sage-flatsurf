from editor_actor import *
from editor_renderer import *
from edge_gluings import *
from similarity_surface import *
from surface_bundle import *

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector

from Tkinter import *
import tkSimpleDialog
import tkMessageBox

from collections import defaultdict

class SimilaritySurfaceBundle(SurfaceBundle, EditorRenderer):

    def __init__(self, name, editor, similarity_surface):
        r"""
        INPUT:

        - ``name`` - a name for the surface (a string)

        - ``editor`` - the identification of the edges. A list of pairs
          ((p0,e0),(p1,e1)) or
        """
        SurfaceBundle.__init__(self,name, editor, field=similarity_surface.base_ring() )
        EditorRenderer.__init__(self,editor)
        self._ss=similarity_surface
        # matrices which act affinely on the polygons:
        self._gl=defaultdict(self._default_gl)
        # translation vectors for the polygons:
        self._t=defaultdict(self._default_t)
        # List of visible polygons
        self._visible=set()
        # cache for transformed vertices
        self._polygon_cache={}
        # handles for the polygons in the canvas
        self._polygon_to_handles={}
        # stores labels:
        self._edge_labels=SurfaceLabels(self._ss)
        
    def _default_gl(self):
        return identity_matrix( 2, self.field() )

    def _default_t(self):
        return self.vector_space().zero()

    def get_transformed_vertices(self, i):
        r"""
        Get the transformed vertices for a polygon.
        """
        try:
            return self._polygon_cache[i]
        except KeyError:
            return self.reset_transformed_vertices(i)

    def redraw_all(self):
        self._editor.get_canvas().delete("all")
        self._visible=set()
        self._polygon_cache={}
        self._polygon_to_handles={}
        self.initial_render()

    def initial_render(self):
        if self._ss.polygons().is_finite():
            for i in self._ss.polygons().keys():
                self._visible.add(i)
        self._render_all_polygons()
        self._render_all_edge_labels()

    def is_visible(self,polygon_index):
        r"""
        Return true if the polygon with the given index is marked as visible.
        """
        return polygon_index in self._visible

    def is_adjacent(self,p1,e1):
        r"""
        Return true if the given edge, e of the polygon p, is adjacent to its neighbor.
        """
        p2,e2 = self._ss.opposite_edge(p1,e1)
        if self.is_visible(p1) and self.is_visible(p2):
            vs1=self.get_transformed_vertices(p1)
            vs2=self.get_transformed_vertices(p2)
            return vs1[e1]==vs2[(e2+1)%len(vs2)] and vs1[(e1+1)%len(vs1)]==vs2[e2]
        return 0

    def reset_transformed_vertices(self, i):
        r"""
        Reset the cache storing the transformed vertices of polygon i.
        """
        t=self._t[i]
        gl=self._gl[i]
        imgs=[]
        vs=self._ss.polygon(i).vertices()
        for i in range(self._ss.polygon(i).num_edges()):
            imgs.append(gl*vs[i]+t)
        res=tuple(imgs)
        self._polygon_cache[i]=res
        return res

    def _render_polygon_fill(self,i):
        vs=self.get_transformed_vertices(i)
        imgs=[]
        for v in vs:
            img = self.math_to_screen_coordinates( v )
            imgs.append(img[0])
            imgs.append(img[1])
        handle=self._editor.get_canvas().create_polygon(imgs,
            outline="",fill="white",tags=("SimilaritySurfaceBundle","polygon"))
        self._polygon_to_handles[i]=handle

    def _render_edge(self,p,e):
        vs=self.get_transformed_vertices(p)
        img0=self.math_to_screen_coordinates( vs[e] )
        img1=self.math_to_screen_coordinates( vs[(e+1)%len(vs)] )
        if self.is_adjacent(p,e):
            handle=self._editor.get_canvas().create_line(
                img0[0], img0[1],img1[0],img1[1],
                fill="#f00",tags=("SimilaritySurfaceBundle","edge"),dash=(1, 4))
        else:
            handle=self._editor.get_canvas().create_line(
                img0[0], img0[1],img1[0],img1[1],
                fill="#f00",tags=("SimilaritySurfaceBundle","edge"))

    def _render_all_edge_labels(self):
        print "rendering all edge labels"
        for p in self._visible:
            vs=self.get_transformed_vertices(p)
            for e in range(len(vs)):
                if not self.is_adjacent(p,e):
                    p2,e2=self._ss.opposite_edge(p,e)
                    if p2 in self._visible:
                        self._render_edge_label(p,e)

    def _render_all_polygons(self):
        for p in self._visible:
            self._render_polygon_fill(p)
        for p in self._visible:
            self._render_polygon_outline(p)

    def _render_edge_label(self, p, e):
        canvas = self._editor.get_canvas()
        vs=self.get_transformed_vertices(p)
        mid=self.math_to_screen_coordinates( (vs[e] + vs[(e+1)%len(vs)])/2 )
        dx, dy = vs[(e+1)%len(vs)] - vs[e]
        dy = -dy
        label=self._edge_labels.get_label(p,e)
        print "label="+str(label)
        if dx>0:
            if dy>0:
                anchor="ne"
            elif dy<0:
                anchor="nw"
            else:
                anchor="n"
        elif dx<0:
            if dy>0:
                anchor="se"
            elif dy<0:
                anchor="sw"
            else:
                anchor="s"
        else:
            if dy>0:
                anchor="e"
            elif dy<0:
                anchor="w"
            else:
                anchor="center"
        print "mid="+str(mid)
        label_handle=canvas.create_text(mid[0], mid[1], text=label, 
            fill="#7efffc", font=("Helvetica","12"), anchor=anchor, 
            tags=("SimilaritySurfaceBundle","label") )

    def _render_polygon_outline(self,p):
        vs=self.get_transformed_vertices(p)
        for e in range(len(vs)):
            self._render_edge(p,e)

    def set_polygon_translation(self, polygon_index, x, y):
        F=self.field()
        self._t[polygon_index]=vector([F(x),F(y)])
        self.reset_transformed_vertices(polygon_index)

