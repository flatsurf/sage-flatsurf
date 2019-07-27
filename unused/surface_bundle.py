from editor_actor import *
from editor_renderer import *
from graphical.edge_gluings import *

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
        self._s = self._field(1)
        # translation constants
        self._tx = self._field.zero()
        self._ty = self._field.zero()

    def after_zoom_change(self):
        pass

    def before_zoom_change(self):
        pass

    def field(self):
        return self._field

    def get_canvas(self):
        r"""Returns the canvas of the editor."""
        return self._editor.get_canvas()

    def get_editor(self):
        r"""Return the editor manipulating this surface"""
        return self._editor

    def get_surface(self):
        return None

    def get_transform(self):
        r"""
        Get the parts of the transformation which convert to screen coordinates.
        """
        return (self._s, self._tx, self._ty)


    def make_menus(self,menubar):
        self._action_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Action", underline=0, menu=self._action_menu)
        self._action_menu.add_command(label="Rename", underline=2, command=self.launch_rename_dialog)
        self._action_menu.add_command(label="Recenter", underline=2, command=self._editor._on_recenter)
        self._action_menu.add_command(label="Zoom", underline=0, command=self._editor._on_zoom,accelerator="Ctrl+Z")
        self._editor.bind_all("<Control-z>", lambda event: self._editor._on_zoom() )
        self._action_menu.add_command(label="Zoom Box", command=self._editor._on_zoom_box)
        self._action_menu.add_command(label="Redraw All", underline=0, command=self._editor._on_redraw_all)
        self._create_menu=Menu(menubar, tearoff=0)
        

    def get_name(self):
        return self._name

    def launch_rename_dialog(self):
        RenameDialog(self)

    def math_to_screen_coordinates(self,v):
        return ( RDF(self._s*v[0]+self._tx), RDF(-self._s*v[1]+self._ty) )

    def redraw_all(self):
        pass

    def recenter_screen(self,x0,y0):
        c=self._editor.get_center()
        self.set_transform(self._s, self._tx+c[0]-x0, self._ty+c[1]-y0)

    def screen_to_math_coordinates(self, x, y):
        x=QQ(x)
        y=QQ(y)
        return self.vector_space()(( (self._field(x)-self._tx)/self._s, (-self._field(y)+self._ty)/self._s ))

    def set_name(self, newname):
        self._name=newname
        self._editor.surface_renamed()

    def set_transform(self, s, tx, ty):
        r"""
        Set the parts of the transformation which convert to screen coordinates.
        """
        s=QQ(s)
        tx=QQ(tx)
        ty=QQ(ty)
        ratio=self._s/s
        if (ratio > QQ(999)/1000) and (ratio < QQ(1001)/1000):
            # ignore negligible change in scale!
            self._editor.get_canvas().move(ALL,RDF(tx-self._tx),RDF(ty-self._ty))
            self._tx=self._field(tx)
            self._ty=self._field(ty)
        else:
            self.before_zoom_change()
            scale=1/ratio
            offset_x=((self._s*tx)-(self._tx*s))/(self._s-s)
            offset_y=((self._s*ty)-(self._ty*s))/(self._s-s)
            self._editor.get_canvas().scale(ALL,RDF(offset_x),RDF(offset_y),RDF(scale),RDF(scale))
            self._s=self._field(s)
            self._tx=self._field(tx)
            self._ty=self._field(ty)
            self.after_zoom_change()

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
        self.set_transform(self._s * factor, factor*self._tx + x0, factor*self._ty + y0)

    def zoom_math_box(self,x1,y1,x2,y2):
        x1p,y2p=self.math_to_screen_coordinates(self.vector_space()((x1,y1)))
        x2p,y1p=self.math_to_screen_coordinates(self.vector_space()((x2,y2)))
        self.zoom_screen_box(x1p,y1p,x2p,y2p)

    def editor_ready():
        return self._editor.get_width()!=1

    def zoom_screen_box(self,x1,y1,x2,y2):
        r"""
        Scale picture by a factor fixing the point (xc,yc) in screen coordinates
        - ``factor`` -- rational scaling factor

        - ``xc`` -- rational or integer point in screen coordinates
        - ``yc`` -- rational or integer point in screen coordinates
        """
        width=x2-x1
        height=y2-y1
        screen_width=self._editor.get_width()
        screen_height=self._editor.get_height()
        width_change=QQ(screen_width)/width
        height_change=QQ(screen_height)/height
        # Commenting this line out seems to have solved the label loss issue:
        #self.before_zoom_change()

        # proposed scale change:
        scale_change=min(width_change,height_change)
        new_scale=self._s*scale_change

        screen_center=self._editor.get_center()
        vertex1=self.screen_to_math_coordinates(x1,y1)
        vertex2=self.screen_to_math_coordinates(x2,y2)
        math_center_x=(vertex1[0]+vertex2[0])/2
        math_center_y=(vertex1[1]+vertex2[1])/2
        tx=screen_center[0]-(new_scale*math_center_x)
        ty=screen_center[1]+(new_scale*math_center_y)
        self.set_transform(new_scale,tx,ty)


