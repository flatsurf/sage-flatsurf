from editor_actor import *
from editor_renderer import *

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import VectorSpace
from sage.rings.real_double import RDF
from sage.rings.rational_field import QQ

from Tkinter import ALL

class SurfaceBundle:
    def __init__(self, name, editor, field=QQ):
        self.name = name
        self._editor=editor
        self._field=field
        """scaling constants"""
        self._sx=self._field(1)
        self._sy=self._field(-1)
        """translation constants"""
        self._tx=self._field.zero()
        self._ty=self._field.zero()
    def make_create_menu(self,menu):
        pass
    def recenter_screen(self,x0,y0):
        c=self._editor.get_center()
        self._editor.get_canvas().move(ALL,c[0]-x0,c[1]-y0)
        self._tx=self._tx+c[0]-x0
        self._ty=self._ty+c[1]-y0
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
    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        return VectorSpace(self._field, 2)
    def screen_to_math_coordinates(self, x, y):
        return self.vector_space()(( (self._field(x)-self._tx)/self._sx, (self._field(y)-self._ty)/self._sy ))
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

class CreateSimilaritySurfaceBundle(SurfaceBundle, EditorRenderer):
    def __init__(self, n, editor):
        SurfaceBundle.__init__(self,"New Similarity Surface "+str(n), editor)
        EditorRenderer.__init__(self,editor)
        self._polygons=[]
        self._translations=[]

    def initial_render(self):
        for i in range(len(self._polygons)):
            self.render_polygon(i)
            
    def render_polygon(self,i):
        coords=self.polygon_to_screen_coordinates(self._polygons[i],self._translations[i])
        self.editor.get_canvas().create_polygon(coords,
            outline="#f00",fill="",tags="CreateSimilaritySurfaceBundle")

    def make_create_menu(self,menu):
        menu.add_command(label="Draw Polygon", command=self.on_draw_polygon)
        menu.add_command(label="Demo Actor", command=self.on_demo)
        menu.add_command(label="Location Actor", command=self.on_location)
        menu.add_command(label="Redraw Polygons", command=self.redraw)

    def on_draw_polygon(self):
        pd=PolygonDrawer(self.editor,self,self.receive_drawn_polygon)
        self.editor.set_actor(pd)

    def on_demo(self):
        da=DemoActor(self.editor)
        self.editor.set_actor(da)

    def redraw(self):
        self.editor.get_canvas().delete("CreateSimilaritySurfaceBundle")
        self.initial_render()

    def on_location(self):
        la=LocationActor(self.editor)
        self.editor.set_actor(la)

    def receive_drawn_polygon(self,polygon,v0x,v0y):
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
