from editor_actor import *
from editor_renderer import *
from edge_gluings import *


from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ

from Tkinter import *
import tkSimpleDialog
import tkMessageBox

class RenameDialog(tkSimpleDialog.Dialog):
    r"""
    Rename dialog for a SurfaceBundle.
    """
    def __init__(self,bundle):
        self._bundle=bundle
        tkSimpleDialog.Dialog.__init__(self,bundle.get_editor().get_parent(),title="Rename surface")

    def body(self, master):
        Label(master, text="Name:").grid(row=0)
        self.entry = Entry(master)
        self.entry.insert(0, self._bundle.get_name())
        self.entry.grid(row=0, column=1)
        return self.entry # initial focus

    def apply(self):
        self._bundle.set_name(self.entry.get())


class SurfaceBundle:

    def __init__(self, name, editor, field=QQ):
        self._name = name
        self._editor = editor
        self._field = field
        # scaling constants
        self._sx = self._field(1)
        self._sy = self._field(-1)
        # translation constants
        self._tx = self._field.zero()
        self._ty = self._field.zero()

    def get_canvas(self):
        r"""Returns the canvas of the editor."""
        return self._editor.get_canvas()

    def get_editor(self):
        r"""Return the editor manipulating this surface"""
        return self._editor

    def make_create_menu(self,menu):
        pass

    def make_action_menu(self,menu):
        pass

    def get_name(self):
        return self._name

    def launch_rename_dialog(self):
        RenameDialog(self)

    def math_to_screen_coordinates(self,v):
        return ( RDF(self._sx*v[0]+self._tx), RDF(self._sy*v[1]+self._ty) )

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

    def recenter_screen(self,x0,y0):
        c=self._editor.get_center()
        self._editor.get_canvas().move(ALL,c[0]-x0,c[1]-y0)
        self._tx=self._tx+c[0]-x0
        self._ty=self._ty+c[1]-y0

    def screen_to_math_coordinates(self, x, y):
        return self.vector_space()(( (self._field(x)-self._tx)/self._sx, (self._field(y)-self._ty)/self._sy ))

    def set_name(self, newname):
        self._name=newname
        self._editor.surface_renamed()

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        return VectorSpace(self._field, 2)

    def zoom(self,factor,xc,yc):
        r"""
        Scale picture by a factor fixing the point (xc,yc) in screen coordinates
        - ``factor`` -- rational scaling factor

        - ``xc`` -- rational or integer point in screen coordinates
        - ``yc`` -- rational or integer point in screen coordinates
        """
        x0=xc-factor*xc
        y0=yc-factor*yc
        print "scale by "+str(RDF(factor))+" xc="+str(xc)+" yc="+str(yc)
        self._editor.get_canvas().scale(ALL,xc,yc,RDF(factor),RDF(factor))
        self._sx = self._sx * factor
        self._sy = self._sy * factor
        self._tx = factor*self._tx + x0
        self._ty = factor*self._ty + y0



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
                if self._edge_to_label_handle.has_key((index,edge)):
                    label_handle = self._edge_to_label_handle[(index,edge)]
                    canvas.delete(label_handle)
                    del self._edge_to_label_handle[(index,edge)]
                    del self._handles_to_labels[label_handle]
                label = self._glue.get_label(index,edge)
                if label is not None:
                    midx=(coords[2*edge]+coords[2*edge+2])/2
                    dx=coords[2*edge+2]-coords[2*edge]
                    midy=(coords[2*edge+1]+coords[2*edge+3])/2
                    dy=coords[2*edge+3]-coords[2*edge+1]
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
        menu.add_command(label="Polygon", command=self.on_draw_polygon)
        #menu.add_command(label="Produce Similarity Surface", command=self._on_glue)

    def make_action_menu(self,menu):
        menu.add_separator()
        menu.add_command(label="Rename surface", command=self.launch_rename_dialog)
        menu.add_command(label="Glue pair of edges", command=self._on_glue)


    def on_draw_polygon(self):
        pd=PolygonDrawer(self._editor,self,self._receive_drawn_polygon)
        self._editor.set_actor(pd)

    def _on_glue(self):
        ps=EdgePairSelector(self._editor,self._glue_receive)
        self._editor.set_actor(ps)

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
