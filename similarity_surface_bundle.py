from editor_actor import *
from editor_renderer import *
from edge_gluings import *
from geometry.similarity_surface import *
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
import threading
import time




class SimilaritySurfaceBundle(SurfaceBundle, EditorRenderer):

    # STATIC VARIABLES
    count = 0

    def __init__(self, similarity_surface, editor=None, name=None):
        r"""
        INPUT:

        - ``name`` - a name for the surface (a string)

        - ``editor`` - the identification of the edges. A list of pairs
          ((p0,e0),(p1,e1)) or
        """
        if name is None:
            #SimilaritySurfaceBundle.count = SimilaritySurfaceBundle.count + 1
            #name = "Similarity Surface #"+ str(SimilaritySurfaceBundle.count)
            name=repr(similarity_surface)
        if editor is None:
            from surface_manipulator import SurfaceManipulator
            editor = SurfaceManipulator.launch()
        SurfaceBundle.__init__(self, name, editor, field=similarity_surface.base_ring() )
        EditorRenderer.__init__(self, editor)
        self._ss = similarity_surface
        # matrices which act affinely on the polygons:
        self._gl = defaultdict(self._default_gl)
        # translation vectors for the polygons:
        self._t = defaultdict(self._default_t)
        # List of visible polygons
        self._visible = set()
        self._visible.add(self._ss.base_label())
        # cache for transformed vertices
        self._polygon_cache = {}
        # handles for the polygons in the canvas
        self._polygon_to_handle = {}
        self._handle_to_polygon = {}
        # stores labels:
        self._edge_labels = SurfaceLabels(self._ss)

    def __repr__(self):
        s = "Similarity surface bundle"
        s += "  _handle_to_polygon : {!r}\n".format(self._handle_to_polygon)
        s += "  _polygon_to_handle : {!r}\n".format(self._polygon_to_handle)
        s += "  _visible           : {!r}\n".format(self._visible)
        return s

    def after_zoom_change(self):
        self._render_all_edge_labels()

    def before_zoom_change(self):
        # remove all labels
        self.remove_all_labels()
        
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
        if self._ss.is_finite():
            for i in self._ss.polygon_labels():
                self._visible.add(i)
        self._render_all_polygons()
        self._render_all_edge_labels()

    def zoom_fit_nice_boundary(self):
        self.zoom_fit(boundary=QQ(1)/20)

    def zoom_fit(self,boundary=0):
        x1,y1,x2,y2=self.get_math_bbox()
        #print "fit math box="+str((x1,y1,x2,y2))
        if boundary:
            boundary=QQ(boundary)
            xc=(x1+x2)/2
            yc=(y1+y2)/2
            width=x2-x1
            height=y2-y1
            x1=xc-(1+boundary)*width/2
            x2=xc+(1+boundary)*width/2
            y1=yc-(1+boundary)*height/2
            y2=yc+(1+boundary)*height/2
        self.zoom_math_box(x1,y1,x2,y2)
        #self.remove_all_labels()
        #bbox=self._editor.get_canvas().bbox(ALL)
        #x1,y1,x2,y2=bbox
        #print "bbox="+str(bbox)
        #xc=(QQ(x1)+x2)/2
        #yc=(QQ(y1)+y2)/2
        #width=QQ(x2)-x1
        #height=QQ(y2)-y1
        #self.zoom_screen_box(xc-11*width/20,yc-11*height/20,xc+11*width/20,yc+11*height/20)


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

    def make_menus(self,menubar):
        # Setup Surface Bundle's Action Menu
        SurfaceBundle.make_menus(self,menubar)
        self._action_menu.add_command(label="Zoom fit", command=self.zoom_fit_nice_boundary,accelerator="Ctrl+F")
        self._editor.bind_all("<Control-f>", lambda event: self.zoom_fit_nice_boundary() )
        self._action_menu.add_separator()
        self._action_menu.add_command(label="Make adjacent", command=self._on_make_adjacent)
        self._action_menu.add_command(label="Move show", command=self._on_move_show)
        menubar.add_cascade(label="Create", underline=0, menu=self._create_menu)
        self._create_menu.add_command(label="Straight-line trajectory", command=self.construct_trajectory)


    def make_visible(self, polygon_index):
        if not self.is_visible(polygon_index):
            self._visible.add(polygon_index)
            self._render_polygon_fill(polygon_index)
            self._render_polygon_outline(polygon_index)
            self._render_polygon_edge_labels(polygon_index)

    def _on_move_show(self):
        ps=PolygonEdgeDragger(self._editor,self._on_move_show_callback)
        self._editor.set_actor(ps)

    def _on_move_show_callback(self, polygon_handle, e1):
        p1 = self._handle_to_polygon[polygon_handle]
        p2,e2 = self._ss.opposite_edge(p1,e1)
        #print "p1="+str(p1)+" e1="+str(e1)+" p2="+str(p2)+" e2="+str(e2)
        if not (p2 in self._visible):
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

    def _on_make_adjacent(self):
        ps = EdgeSelector(self._editor, self._on_make_adjacent_callback)
        self._editor.set_actor(ps)

    def _on_make_adjacent_callback(self, polygon_handle, e1):
        p1 = self._handle_to_polygon[polygon_handle]
        p2,e2 = self._ss.opposite_edge(p1,e1)
        #print "p1="+str(p1)+" e1="+str(e1)+" p2="+str(p2)+" e2="+str(e2)
        m1 = self._gl[p1]
        mc = self._ss.edge_matrix(p2,e2)
        m2 = m1*mc
        vs = self.get_transformed_vertices(p1)
        pt1 = vs[e1]
        vs = self._ss.polygon(p2).vertices()
        pt2 = m2*vs[(e2+1)%self._ss.polygon(p2).num_edges()]
        t2 = pt1-pt2
        self.set_polygon_view(p2,m2,t2[0],t2[1])
        self._visible.add(p2)
        self.redraw_all()
        # Do it again:
        self._on_make_adjacent()

    def pick_edge(self):
        self.done_picking=IntVar()
        self.done_picking.set(0)
        ps=EdgeSelector(self._editor,self._pick_edge_callback)
        self._editor.set_actor(ps)
        self._editor.wait_variable(self.done_picking)
        return self.picked_polygon_handle, self.picked_edge

    def _pick_edge_callback(self,polygon_handle, e):
        self.picked_polygon_handle=polygon_handle
        self.picked_edge=e
        self.done_picking.set(1)

    def to_math_coordinates(self, surface_point):
        r"""Convert a surface point to the math version of screen coordinates."""
        l=surface_point.get_label()
        p=surface_point.get_point()
        t=self._t[l]
        gl=self._gl[l]
        return gl*p+t
            
    def pick_point(self,callback):
        r"""Has the user select a point. 
        Returns the point as a SurfacePoint via callback."""
        #self.done_picking=IntVar()
        #self.done_picking.set(0)

        def _pick_point_callback(polygon_handle, x,y):
            from geometry.surface_point import SurfacePoint
            #print "x="+str(x)+" and y="+str(y)
            v=self.screen_to_math_coordinates(x,y)
            #print "v="+str(v)
            i=self._handle_to_polygon[polygon_handle]
            t=self._t[i]
            gl=self._gl[i]
            p=gl.inverse()*(v-t)
            try:
                #self._picked_point=SurfacePoint(self._ss,i,p)
                callback(SurfacePoint(self._ss,i,p))
            except ValueError:
                #self._picked_point=None
                callback(None)
            #self.done_picking.set(1)

        ps=PointSelector(self._editor,_pick_point_callback)
        self._editor.set_actor(ps)
        #self._editor.wait_variable(self.done_picking)


    def remove_all_labels(self):
        self._editor.get_canvas().delete("label")

    def draw_flow(self, holonomy):
        self._holonomy=self._ss.vector_space()((holonomy[0],holonomy[1]))
        self.done_picking=IntVar()
        self.done_picking.set(0)
        ps=PointSelector(self._editor,self._draw_flow_callback)
        self._editor.set_actor(ps)
        

    def _draw_flow_callback(self,polygon_handle, x,y):
        from geometry.surface_point import SurfacePoint
        #print "x="+str(x)+" and y="+str(y)
        v=self.screen_to_math_coordinates(x,y)
        #print "v="+str(v)
        i=self._handle_to_polygon[polygon_handle]
        t=self._t[i]
        gl=self._gl[i]
        p=gl.inverse()*(v-t)
        #try:
        pt=SurfacePoint(self._ss,i,p)
        segments=pt.flow_segments(self._holonomy)
        self.render_segments(segments)
        #except ValueError:
        #    print
        #    pass

    def draw_flow_with_skip(self, holonomy):
        self._holonomy=self._ss.vector_space()((holonomy[0],holonomy[1]))
        self.done_picking=IntVar()
        self.done_picking.set(0)
        ps=PointSelector(self._editor,self._draw_flow_with_skip_callback)
        self._editor.set_actor(ps)
        

    def _draw_flow_with_skip_callback(self,polygon_handle, x,y):
        from geometry.surface_point import SurfacePoint
        #print "x="+str(x)+" and y="+str(y)
        v=self.screen_to_math_coordinates(x,y)
        #print "v="+str(v)
        i=self._handle_to_polygon[polygon_handle]
        t=self._t[i]
        gl=self._gl[i]
        p=gl.inverse()*(v-t)
        #try:
        pt=SurfacePoint(self._ss,i,p)
        segments=pt.flow_segments_wait(self._holonomy)
        self.render_segments(segments)
        #except ValueError:
        #    print
        #    pass


    def redraw_all(self):
        r"""
        Remove and redraw everything on the canvas.
        """
        self._editor.get_canvas().delete("all")
        #self._visible=set()
        self._polygon_cache={}
        self._polygon_to_handle={}
        self._handle_to_polygon={}
        self._render_all_polygons()
        self._render_all_edge_labels()

    def reset_transformed_vertices(self, i):
        r"""
        Reset the cache storing the transformed vertices of polygon i.
        """
        #print "Updating i="+str(i)
        t=self._t[i]
        gl=self._gl[i]
        imgs=[]
        vs=self._ss.polygon(i).vertices()
        for j in range(self._ss.polygon(i).num_edges()):
            imgs.append(gl*vs[j]+t)
        res=tuple(imgs)
        self._polygon_cache[i]=res
        return res

    def get_math_bbox(self):
        first = True
        xmin=0
        xmax=0
        ymin=0
        ymax=0
        for i in self._visible:
            vs=self.get_transformed_vertices(i)
            for v in vs:
                #print "v="+str(v)
                if first:
                    #print "running first"
                    xmin=v[0]
                    xmax=v[0]
                    ymin=v[1]
                    ymax=v[1]
                    first=False
                else:
                    if v[0]<xmin:
                        xmin=v[0]
                    #print "comparing v0="+str(v[0])+" with xmax="+str(xmax)
                    if v[0]>xmax:
                        #print "v0 is larger"
                        xmax=v[0]
                    #else:
                    #    print "v0 is smaller"
                    if v[1]<ymin:
                        ymin=v[1]
                    if v[1]>ymax:
                        ymax=v[1]
        #print "bbox="+str((xmin,ymin,xmax,ymax))
        return xmin,ymin,xmax,ymax

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

    def render_segments(self,pair_list):
        for pair in pair_list:
            a,b=pair
            i=a.get_label()
            if i in self._visible:
                t=self._t[i]
                gl=self._gl[i]
                img0=self.math_to_screen_coordinates( gl*a.get_point()+t )
                img1=self.math_to_screen_coordinates( gl*b+t )
                handle=self._editor.get_canvas().create_line(
                    img0[0], img0[1],img1[0],img1[1],
                    fill="#000",tags=("segment"))


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
            fill="#7efffc", 
            font=("Helvetica","12"), anchor=anchor, 
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

    def construct_trajectory(self):
        #self._ct_start=None
        #self._ct_start_screen=None
        def _callback1(pt,dictionary={}):
            if pt is None:
                return
            dictionary["start"]=pt
            dictionary["start_screen"]=self.to_math_coordinates(pt)
            vs=VectorSelector(self._editor,self.math_to_screen_coordinates(dictionary["start_screen"]),_callback2,dictionary=dictionary)
            self._editor.set_actor(vs)
        def _callback2(x,y,dictionary={}):
            dictionary["end"]=self.screen_to_math_coordinates(x,y)
            _LengthDialog(self,_callback3,dictionary=dictionary)
            
        class _LengthDialog(tkSimpleDialog.Dialog):
            def __init__(self,bundle,callback3,dictionary={}):
                self._callback3=callback3
                self._dog=dictionary
                tkSimpleDialog.Dialog.__init__(self,bundle.get_editor().get_parent(),title="Length of trajectory")

            def body(self, master):
                Label(master, text="Length multiplier:").grid(row=0)
                self.entry = Entry(master)
                self.entry.insert(0, "1000")
                self.entry.grid(row=0, column=1)
                return self.entry # initial focus

            def apply(self):
                self._callback3(self.entry.get(), dictionary=self._dog)

        def _callback3(length,dictionary={}):
            hol=10000*(dictionary["end"]-dictionary["start_screen"])
            segments=dictionary["start"].flow_segments(hol)
            self.render_segments(segments)
            self._editor.set_actor(None)

        pt=self.pick_point(_callback1)


