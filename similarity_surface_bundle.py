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

    # STATIC VARIABLES
    count = 0

    def __init__(self, similarity_surface, editor = None, name = None):
        r"""
        INPUT:

        - ``name`` - a name for the surface (a string)

        - ``editor`` - the identification of the edges. A list of pairs
          ((p0,e0),(p1,e1)) or
        """
        if name is None:
            SimilaritySurfaceBundle.count = SimilaritySurfaceBundle.count + 1
            name = "Similarity Surface #"+ str(SimilaritySurfaceBundle.count)
        if editor is None:
            from surface_manipulator import SurfaceManipulator
            editor = SurfaceManipulator.launch()
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
        self._polygon_to_handle={}
        self._handle_to_polygon={}
        # stores labels:
        self._edge_labels=SurfaceLabels(self._ss)

    def after_zoom_change(self):
        self._render_all_edge_labels()

    def before_zoom_change(self):
        # remove all labels
        self._editor.get_canvas().delete("label")
        
    def _default_gl(self):
        return identity_matrix( self.field(), n=2)

    def _default_t(self):
        return self.vector_space().zero()

    def get_surface(self):
        return self._ss

    def get_transformed_vertices(self, i):
        r"""
        Get the transformed vertices for a polygon.
        """
        try:
            return self._polygon_cache[i]
        except KeyError:
            return self.reset_transformed_vertices(i)

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

    def make_create_menu(self,menu):
        #menu.add_command(label="Polygon", command=self._on_draw_polygon)
        #menu.add_command(label="To Similarity Surface", command=self._on_to_similarity_surface)
        #menu.add_command(label="Produce Similarity Surface", command=self._on_glue)
        pass

    def make_action_menu(self,menu):
        menu.add_separator()
        menu.add_command(label="Make adjacent", command=self._on_make_adjacent)

    def make_visible(self, polygon_index):
        if not self.is_visible(polygon_index):
            self._visible.add(polygon_index)
            self._render_polygon_fill(polygon_index)
            self._render_polygon_outline(polygon_index)
            self._render_polygon_edge_labels(polygon_index)

    def _on_make_adjacent(self):
        ps=EdgeSelector(self._editor,self._on_make_adjacent_callback)
        self._editor.set_actor(ps)

    def _on_make_adjacent_callback(self, polygon_handle, e1):
        p1=self._handle_to_polygon[polygon_handle]
        p2,e2 = self._ss.opposite_edge(p1,e1)
        m1=self._gl[p1]
        mc=self._ss.edge_matrix(p2,e2)
        m2=m1*mc
        vs=self.get_transformed_vertices(p1)
        pt1=vs[e1]
        vs=self._ss.polygon(p2).vertices()
        pt2=m2*vs[(e2+1)%self._ss.polygon(p2).num_edges()]
        t2=pt1-pt2
        self.set_polygon_view(p2,m2,t2[0],t2[1])
        self._visible.add(p2)
        self.redraw_all()
        # Do it again:
        ps=EdgeSelector(self._editor,self._on_make_adjacent_callback)
        self._editor.set_actor(ps)

    def redraw_all(self):
        r"""
        Remove and redraw everything on the canvas.
        """
        self._editor.get_canvas().delete("all")
        self._visible=set()
        self._polygon_cache={}
        self._polygon_to_handle={}
        self._handle_to_polygon={}
        self.initial_render()

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
        self._polygon_to_handle[i]=handle
        self._handle_to_polygon[handle]=i

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
        for p in self._visible:
            self._render_polygon_edge_labels(p)

    def _render_polygon_edge_labels(self, p):
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
        offset=(0,0)
        if dx>0:
            if dy>0:
                anchor="ne"
                offset=(1,-1)
            elif dy<0:
                anchor="nw"
                offset=(-1,-1)
            else:
                anchor="n"
                offset=(0,-2)
        elif dx<0:
            if dy>0:
                anchor="se"
                offset=(1,1)
            elif dy<0:
                anchor="sw"
                offset=(-1,1)
            else:
                anchor="s"
                offset=(0,2)
        else:
            if dy>0:
                anchor="e"
                offset=(2,0)
            elif dy<0:
                anchor="w"
                offset=(-2,0)
            else:
                anchor="center"
                offset=(0,0)
        label_handle=canvas.create_text(mid[0]-offset[0], mid[1]-offset[1],
            text=label, 
            fill="#7efffc", font=("Helvetica","12"), anchor=anchor, 
            tags=("SimilaritySurfaceBundle","label") )

    def _render_polygon_outline(self,p):
        vs=self.get_transformed_vertices(p)
        for e in range(len(vs)):
            self._render_edge(p,e)

    def set_polygon_view(self, polygon_index, m, tx, ty):
        F=self.field()
        self._gl[polygon_index]=m
        self._t[polygon_index]=vector([F(tx),F(ty)])
        self.reset_transformed_vertices(polygon_index)

    def set_polygon_translation(self, polygon_index, x, y):
        F=self.field()
        self._t[polygon_index]=vector([F(x),F(y)])
        self.reset_transformed_vertices(polygon_index)

