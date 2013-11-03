from Tkinter import *

class EditorActor:
    def __init__(self, editor):
        self._editor = editor        
    def on_activate(self):
        pass
    def on_deactivate(self):
        pass
    def single_left_click(self, event):
        pass

    def double_left_click(self, event):
        pass

    def triple_left_click(self, event):
        pass

    def single_middle_click(self, event):
        pass    

    def double_middle_click(self, event):
        pass

    def triple_middle_click(self, event):
        pass

    def single_right_click(self, event):
        pass

    def double_right_click(self, event):
        pass

    def triple_right_click(self, event):
        pass

    def shift_click(self, event):
        pass

    def mouse_moved(self, event):
        pass

    def focus_in(self, event):
        pass
    def focus_out(self, event):
        pass
    def key_press(self, event):
        pass
    def key_release(self, event):
        pass

class RecenterActor(EditorActor):
    def __init__(self, editor):
        EditorActor.__init__(self, editor)

    def on_activate(self):
        self._editor.set_text("Click to recenter the screen.")

    def recenter(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        bundle=self._editor.get_surface_bundle()
        if bundle is None:
            pass
        else:
            bundle.recenter_screen(x,y)

    def single_left_click(self, event):
        self.recenter(event)

    def single_middle_click(self, event):
        self.recenter(event)

    def single_right_click(self, event):
        self.recenter(event)

class ZoomActor(EditorActor):
    def __init__(self, editor):
        EditorActor.__init__(self, editor)

    def on_activate(self):
        self._editor.set_text("Left click to zoom in. Right click to zoom out.")

    def recenter(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        bundle=self._editor.get_surface_bundle()
        if bundle is None:
            pass
        else:
            bundle.recenter_screen(x,y)

    def single_left_click(self, event):
        from sage.rings.rational_field import QQ
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        bundle=self._editor.get_surface_bundle()
        if bundle is None:
            pass
        else:
            bundle.zoom(QQ(8)/QQ(7),x,y)

    def double_left_click(self, event):
        self.single_left_click(event)

    def triple_left_click(self, event):
        self.single_left_click(event)

    def single_right_click(self, event):
        from sage.rings.rational_field import QQ
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        bundle=self._editor.get_surface_bundle()
        if bundle is None:
            pass
        else:
            bundle.zoom(QQ(7)/QQ(8),x,y)

    def double_right_click(self, event):
        self.single_right_click(event)

    def triple_right_click(self, event):
        self.single_right_click(event)



from polygon import PolygonCreator
from sage.rings.rational_field import QQ

class PolygonDrawer(EditorActor):
    def __init__(self, editor, bundle, polygon_receiver, field = QQ):
        EditorActor.__init__(self, editor)        
        self.polygon_receiver=polygon_receiver
        self._pc=PolygonCreator(field)
        self._v=[]
        self._bundle=bundle
    def on_activate(self):
        self._editor.set_text("Left click to start drawing a polygon.")
    def on_deactivate(self):
        self._editor.get_canvas().delete("PolygonDrawer")

    def single_left_click(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        newv=self._bundle.screen_to_math_coordinates(x,y)
        if self._pc.add_vertex(newv):
            self._v.append(x)
            self._v.append(y)
            self._editor.get_canvas().delete("PolygonDrawer")
            if (len(self._v)==2):
                self._editor.get_canvas().create_oval(x-2,y-2,x+2,y+2,outline="",fill="#f70",tags="PolygonDrawer")
            else:
                self._editor.get_canvas().create_line(self._v,fill="#f70",tags="PolygonDrawer")
            if (len(self._v)<=4):
                self._editor.set_text("Left click again to choose a new vertex.")
            else:
                self._editor.set_text("Click again, or press return to stop.")
        else:
            self._editor.set_text("That point would make the polygon non-convex.")

    def key_press(self, event):
        if event.keycode==36:
            self.polygon_receiver(self._pc.get_polygon(),self._v[0],self._v[1])
            self._editor.get_canvas().delete("PolygonDrawer")
            self._editor.set_actor(None)

class LocationActor(EditorActor):
    def __init__(self, editor):
        EditorActor.__init__(self, editor)
    def on_activate(self):
        print "Print coordinates of click"
    def single_left_click(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        bundle=self._editor.get_surface_bundle()
        if bundle is None:
            print "single left click at (%s,%s)"%(x,y)
        else:
            print "math coordinates of click: "+str(bundle.screen_to_math_coordinates(x,y))

class DemoActor(EditorActor):
    def __init__(self, editor):
        EditorActor.__init__(self, editor)

    def on_activate(self):
        print "DemoActor activated"

    def single_left_click(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        print "single left click at (%s,%s)"%(x,y)
        self._editor.get_canvas().create_rectangle(x-5,y-5,x+5,y+5, 
            fill="#eeff00", tags="junk")

    def single_middle_click(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        print "single middle click at (%s,%s)"%(x,y)
        self._editor.get_canvas().create_rectangle(x-5,y-5,x+5,y+5, 
            fill='#00aa0f', tags="junk")

    def single_right_click(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        print "single right click at (%s,%s)"%(x,y)
        self._editor.get_canvas().create_rectangle(x-5,y-5,x+5,y+5, 
            fill='#ff00ff', tags="junk")

    def double_left_click(self, event):
        print "double left click"

    def shift_click(self, event):
        print "shift click"

    def mouse_moved(self, event):
        x = self._editor.get_canvas().canvasx(event.x)
        y = self._editor.get_canvas().canvasy(event.y)
        print "mouse moved to (%s,%s)"%(x,y)

    def focus_in(self, event):
        print "focus in"

    def focus_out(self, event):
        print "focus out"

    def key_press(self, event):
        print "key %s press with code %s"%(event.char,event.keycode)

    def key_release(self, event):
        print "key %s release with code %s"%(event.char,event.keycode)

