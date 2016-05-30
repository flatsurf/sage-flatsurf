from editor_actor import *
from editor_renderer import *
from graphical.edge_gluings import *
from similarity_surface_bundle import *

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ

from Tkinter import *
import tkSimpleDialog
import tkMessageBox

class CreateSimilaritySurfaceBundle(SurfaceBundle, EditorRenderer):

    def __init__(self, n, editor):
        SurfaceBundle.__init__(self,"New Similarity Surface "+str(n), editor)
        EditorRenderer.__init__(self,editor)
        self._polygons=[]
        self._translations=[]
        # List of polygons created 
        self._polygon_to_handles={}
        self._handles_to_polygon={}
        self._glue=EditableEdgeGluing()
        self._handles_to_labels={}
        self._edge_to_label_handle={}

    def initial_render(self):
        self._polygon_to_handles={}
        self._handles_to_polygon={}
        for i in range(len(self._polygons)):
            self.render_polygon(i)
        self.render_edge_labels()
            
    def render_polygon(self,i):
        coords=self.polygon_to_screen_coordinates(self._polygons[i],self._translations[i])
        handle=self._editor.get_canvas().create_polygon(coords,
            outline="#f00",fill="white",tags=("CreateSimilaritySurfaceBundle","polygon"))
        self._polygon_to_handles[i]=handle
        self._handles_to_polygon[handle]=i

    def render_edge_labels(self):
        canvas = self._editor.get_canvas()
        polygons = canvas.find_withtag("polygon")
        for handle in polygons:
            index=self._handles_to_polygon[handle]
            polygon=self._polygons[index]
            n=polygon.num_edges()
            coords=canvas.coords(handle)
            for edge in range(n):
                if (index,edge) in self._edge_to_label_handle:
                    label_handle = self._edge_to_label_handle[(index,edge)]
                    canvas.delete(label_handle)
                    del self._edge_to_label_handle[(index,edge)]
                    del self._handles_to_labels[label_handle]
                label = self._glue.get_label(index,edge)
                if label is not None:
                    midx=(coords[2*edge]+coords[(2*edge+2)%(2*n)])/2
                    dx=coords[(2*edge+2)%(2*n)]-coords[2*edge]
                    midy=(coords[2*edge+1]+coords[(2*edge+3)%(2*n)])/2
                    dy=coords[(2*edge+3)%(2*n)]-coords[2*edge+1]
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
                    label_handle=canvas.create_text(midx, midy, text=label, 
                        fill="#7efffc", font=("Helvetica","12"), anchor=anchor, 
                        tags="CreateSimilaritySurfaceBundle")
                    self._edge_to_label_handle[(index,edge)]=label_handle
                    self._handles_to_labels[label_handle]=(index,edge)

    def make_create_menu(self,menu):
        menu.add_command(label="Polygon", command=self._on_draw_polygon)
        menu.add_command(label="To Similarity Surface", command=self._on_to_similarity_surface)
        #menu.add_command(label="Produce Similarity Surface", command=self._on_glue)

    def make_action_menu(self,menu):
        menu.add_separator()
        menu.add_command(label="Rename surface", command=self.launch_rename_dialog)
        menu.add_command(label="Glue pair of edges", command=self._on_glue)


    def _on_draw_polygon(self):
        pd=PolygonDrawer(self._editor,self,self._receive_drawn_polygon)
        self._editor.set_actor(pd)

    def _on_glue(self):
        ps=EdgePairSelector(self._editor,self._glue_receive)
        self._editor.set_actor(ps)

    def _on_to_similarity_surface(self):
        from geometry.similarity_surface_generators import SimilaritySurface_polygons_and_gluings
        s=SimilaritySurface_polygons_and_gluings(self._polygons,self._glue.get_edge_pair_list())
        sb=SimilaritySurfaceBundle(s, editor=self._editor, name=self._name+" [SS]")
        for i in range(len(self._translations)):
            v=self._translations[i]
            sb.set_polygon_translation(i,v[0],v[1])
        a1,a2,a3 = self.get_transform()
        sb.set_transform(a1,a2,a3)
        self._editor.set_surface(sb)

    def polygon_to_screen_coordinates(self,polygon,translation=None):
        r"""
        Return a list of vertices for a polygon to be drawn to the screen.

        INPUT:

        - ``polygon`` -- a polygon defined over the field

        - ``translation`` -- translation of polygon measured in math coordinates
        """
        ret=[]
        if translation is None:
            for vertex in polygon:
                p=self.math_to_screen_coordinates(vertex)
                ret.append(p[0])
                ret.append(p[1])
        else:
            for vertex in polygon:
                p=self.math_to_screen_coordinates(vertex+translation)
                ret.append(p[0])
                ret.append(p[1])
        return ret

    def _glue_receive(self, ph1, e1, ph2, e2):
        try:
            p1=self._handles_to_polygon[ph1]
            p2=self._handles_to_polygon[ph2]
            self._glue.glue_and_label(p1,e1,p2,e2)
            self.render_edge_labels()
        except ValueError:
            tkMessageBox.showerror("Bad gluing","You can not glue an edge to itself!")
        self._on_glue()

    def on_demo(self):
        da=DemoActor(self._editor)
        self._editor.set_actor(da)

    def redraw(self):
        self._editor.get_canvas().delete("CreateSimilaritySurfaceBundle")
        self.initial_render()

    def _receive_drawn_polygon(self,polygon,v0x,v0y):
        r"""
        Return a list of vertices for a polygon to be drawn to the screen.

        INPUT:

        - ``polygon`` -- a polygon defined over the field

        - ``v0x`` -- screen x-coordinate of first drawn vertex

        - ``v0y`` -- screen y-coordinate of first drawn vertex
        """
        self._polygons.append(polygon)
        self._translations.append(self.screen_to_math_coordinates(v0x,v0y))
        self.render_polygon(len(self._polygons)-1)
