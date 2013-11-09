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

class EditorRedirectingActor(EditorActor):
    r"""
    An actor which redirects input to another actor, which can be changed.
    The purpose of this class is to allow the construction of more complicated
    user actions by stringing together lots of simple actions.
    """
    def __init__(self, editor):
        r"""
        Create an actor which redirects input to a subactor.
        """
        EditorActor.__init__(self, editor)
        self._subactor = None
        self._is_active=0

    def setActor(self, newactor):
        r"""
        Set a new actor to handle input
        """
        if (self._subactor is not None) and (self._is_active): 
            self._subactor.on_deactivate()
        self._subactor=newactor
        if (self._subactor is not None) and (self._is_active): 
            self._subactor.on_activate()

    def getActor(self):
        return self._subactor

    def on_activate(self):
        self._is_active=1
        if self._subactor is not None: 
            self._subactor.on_activate()

    def on_deactivate(self):
        self._is_active=0
        if self._subactor is not None: 
            self._subactor.on_deactivate()

    def single_left_click(self, event):
        if self._subactor is not None: 
            self._subactor.single_left_click(event)

    def double_left_click(self, event):
        if self._subactor is not None: 
            self._subactor.double_left_click(event)

    def triple_left_click(self, event):
        if self._subactor is not None: 
            self._subactor.triple_left_click(event)

    def single_middle_click(self, event):
        if self._subactor is not None: 
            self._subactor.single_middle_click(event)

    def double_middle_click(self, event):
        if self._subactor is not None: 
            self._subactor.double_middle_click(event)

    def triple_middle_click(self, event):
        if self._subactor is not None: 
            self._subactor.triple_middle_click(event)

    def single_right_click(self, event):
        if self._subactor is not None: 
            self._subactor.single_right_click(event)

    def double_right_click(self, event):
        if self._subactor is not None: 
            self._subactor.double_right_click(event)

    def triple_right_click(self, event):
        if self._subactor is not None: 
            self._subactor.triple_right_click(event)

    def shift_click(self, event):
        if self._subactor is not None: 
            self._subactor.shift_click(event)

    def mouse_moved(self, event):
        if self._subactor is not None: 
            self._subactor.mouse_moved(event)

    def focus_in(self, event):
        if self._subactor is not None: 
            self._subactor.focus_in(event)

    def focus_out(self, event):
        if self._subactor is not None: 
            self._subactor.focus_out(event)

    def key_press(self, event):
        if self._subactor is not None: 
            self._subactor.key_press(event)

    def key_release(self, event):
        if self._subactor is not None: 
            self._subactor.key_release(event)


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

class PolygonSelector(EditorActor):
    r"""
    A class for selecting a polygon. 
    The polygons must be tagged with the tag "polygon". 
    When a polygon is clicked the class calls handle_reciever with the handle of the poilygon clicked.
    """
    def __init__(self, editor, handle_receiver, msg="Select a polygon."):
        EditorActor.__init__(self, editor)
        self._over_handle=None
        self._old_fill=None
        self._handle_receiver=handle_receiver
        self._msg=msg

    def on_activate(self):
        self._editor.set_text(self._msg)

    def on_deactivate(self):
        self._revert_higlight()

    def _highlight(self,handle):
        self._over_handle=handle
        self._old_fill=self._editor.get_canvas().itemcget(handle,"fill")
        self._editor.get_canvas().itemconfig(handle,fill="#ffaa66")

    def _revert_higlight(self):
        if self._over_handle is not None:
            self._editor.get_canvas().itemconfig(self._over_handle,fill=self._old_fill)
        self._over_handle=None
        self._old_fill=None

    def mouse_moved(self, event):
        canvas=self._editor.get_canvas()
        handles=canvas.find_withtag(CURRENT)
        if (len(handles)>0):
            handle=handles[0]
        else :
            handle=0
        if self._over_handle is not None:
            if self._over_handle != handle:
                self._revert_higlight()
            else: 
                # still over same polygon
                return
        if handle:
            tags=canvas.gettags(CURRENT)
            if "polygon" in tags:
                # Over polygon
                self._highlight(handle)

    def single_left_click(self, event):
        if self._over_handle is not None:
            self._handle_receiver(self._over_handle)
        # Old Method:
        ## We check to see if a polygon has been clicked and 
        #canvas=self._editor.get_canvas()
        #handle=canvas.find_withtag(CURRENT)
        #if handle:
        #    tags=canvas.gettags(CURRENT)
        #    if "polygon" in tags:
        #        # A polygon was clicked.
        #        # canvas.itemconfig(CURRENT, fill="blue")
        #        # Pass the polygon handle:
        #        self._handle_receiver(handle)

class PolygonEdgeSelector(EditorActor):
    r"""
    A class for selecting an edge from a given polygon.
    """
    def __init__(self, editor, polygon_handle, edge_receiver, msg="Select an edge from the polygon."):
        EditorActor.__init__(self, editor)
        self._polygon_handle=polygon_handle
        self._over_handle=None
        self._polygon_fill=self._editor.get_canvas().itemcget(polygon_handle,"fill")
        self._edge_receiver=edge_receiver
        self._msg=msg

    def on_activate(self):
        handle=self._polygon_handle
        canvas=self._editor.get_canvas()
        self._polygon_fill=canvas.itemcget(self._polygon_handle,"fill")
        canvas.itemconfig(handle,fill="#ffaa66")
        coords=canvas.coords(handle)
        self._edge_handles=[]
        n=len(coords)/2-1
        for i in range(n):
            eh=canvas.create_line(coords[2*i], coords[2*i+1], coords[2*i+2], coords[2*i+3], 
                fill="#dd5500", width=5.0, tags="PolygonEdge")
            self._edge_handles.append(eh)
        self._editor.set_text(self._msg)

    def on_deactivate(self):
        self._editor.get_canvas().delete("PolygonEdge")
        self._editor.get_canvas().itemconfig(self._polygon_handle,fill=self._polygon_fill)

    def _highlight(self,handle):
        self._over_handle=handle
        self._old_fill=self._editor.get_canvas().itemcget(handle,"fill")
        self._editor.get_canvas().itemconfig(handle,fill="#ff7700")

    def _revert_higlight(self):
        if self._over_handle is not None:
            self._editor.get_canvas().itemconfig(self._over_handle,fill=self._old_fill)
        self._over_handle=None
        self._old_fill=None

    def mouse_moved(self, event):
        canvas=self._editor.get_canvas()
        handles=canvas.find_withtag(CURRENT)
        if (len(handles)>0):
            handle=handles[0]
        else :
            handle=0
        if self._over_handle is not None:
            if self._over_handle != handle:
                self._revert_higlight()
            else: 
                # still over same polygon
                return
        if handle:
            tags=canvas.gettags(handle)
            if "PolygonEdge" in tags:
                # Over polygon
                self._highlight(handle)

    def single_left_click(self, event):
        if self._over_handle is not None:
            for i in range(len(self._edge_handles)):
                if self._over_handle == self._edge_handles[i]:
                    self._edge_receiver(i)
                    return

class EdgeSelector(EditorRedirectingActor):
    r"""
    A class for selecting an edge from a collection of polygons polygon. 
    The polygons must be tagged with the tag "polygon". 
    """
    def __init__(self, editor, polygon_handle_and_edge_receiver, 
            polygon_msg="Select a polygon.", 
            edge_msg="Select an edge from the polygon."):
        EditorRedirectingActor.__init__(self, editor)
        self.setActor(PolygonSelector(editor,self._receive_polygon_handle,msg=polygon_msg))
        self._emsg=edge_msg
        self._polygon_handle=None
        self._polygon_handle_and_edge_receiver=polygon_handle_and_edge_receiver

    def _receive_polygon_handle(self,handle):
        self._polygon_handle=handle
        self.setActor(PolygonEdgeSelector(self._editor,handle,self._receive_edge,msg=self._emsg))

    def _receive_edge(self,edge):
        self._polygon_handle_and_edge_receiver(self._polygon_handle,edge)

class EdgePairSelector(EditorRedirectingActor):
    r"""
    A class for selecting a pair of edges from a collection of polygons.
    """
    def __init__(self, editor, pair_receiver):
        EditorRedirectingActor.__init__(self, editor)
        self._pair_receiver=pair_receiver
        self.setActor(EdgeSelector(editor,self._receive_edge_1,
            polygon_msg="Select a polygon.", 
            edge_msg="Select an edge from the polygon."))

    def _receive_edge_1(self,polygon_handle,edge):
        self._ph1=polygon_handle
        self._e1=edge
        self.setActor(EdgeSelector(self._editor,self._receive_edge_2,
            polygon_msg="Select another polygon.", 
            edge_msg="Select an edge from the second polygon."))

    def _receive_edge_2(self,polygon_handle,edge):
        self._ph2=polygon_handle
        self._e2=edge
        self._pair_receiver(self._ph1, self._e1, self._ph2, self._e2)


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

