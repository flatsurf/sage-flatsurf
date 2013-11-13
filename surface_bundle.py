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

    def field(self):
        return self._field

    def get_canvas(self):
        r"""Returns the canvas of the editor."""
        return self._editor.get_canvas()

    def get_editor(self):
        r"""Return the editor manipulating this surface"""
        return self._editor

    def get_transform(self):
        r"""
        Get the parts of the transformation which convert to screen coordinates.
        """
        return (self._sx, self._sy, self._tx, self._ty)

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

    def set_transform(self, sx, sy, tx, ty):
        self._sx=sx
        self._sy=sy
        self._tx=tx
        self._ty=ty
        r"""
        Set the parts of the transformation which convert to screen coordinates.
        """

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
        self._editor.get_canvas().scale(ALL,xc,yc,RDF(factor),RDF(factor))
        self._sx = self._sx * factor
        self._sy = self._sy * factor
        self._tx = factor*self._tx + x0
        self._ty = factor*self._ty + y0



