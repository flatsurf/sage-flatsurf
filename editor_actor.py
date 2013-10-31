from Tkinter import Tk, Frame, Menu, Canvas, Label
import tkFileDialog

class EditorActor:
    def __init__(self, editor):
        self.editor = editor        

    def single_left_click(self, event):
        pass

    def single_middle_click(self, event):
        pass    

    def single_right_click(self, event):
        pass

    def double_click(self, event):
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

class DemoActor(EditorActor):
    def __init__(self, editor):
        EditorActor.__init__(self, editor)        

    def single_left_click(self, event):
        x = self.editor.canvas.canvasx(event.x)
        y = self.editor.canvas.canvasy(event.y)
        print "single left click at (%s,%s)"%(x,y)
        self.editor.canvas.create_rectangle(x-5,y-5,x+5,y+5, 
            fill="#eeff00", tags="junk")

    def single_middle_click(self, event):
        x = self.editor.canvas.canvasx(event.x)
        y = self.editor.canvas.canvasy(event.y)
        print "single middle click at (%s,%s)"%(x,y)
        self.editor.canvas.create_rectangle(x-5,y-5,x+5,y+5, 
            fill='#00aa0f', tags="junk")

    def single_right_click(self, event):
        x = self.editor.canvas.canvasx(event.x)
        y = self.editor.canvas.canvasy(event.y)
        print "single right click at (%s,%s)"%(x,y)
        self.editor.canvas.create_rectangle(x-5,y-5,x+5,y+5, 
            fill='#ff00ff', tags="junk")

    def double_click(self, event):
        print "double click"

    def shift_click(self, event):
        print "shift click"

    def mouse_moved(self, event):
        x = self.editor.canvas.canvasx(event.x)
        y = self.editor.canvas.canvasy(event.y)
        print "mouse moved to (%s,%s)"%(x,y)

    def focus_in(self, event):
        print "focus in"

    def focus_out(self, event):
        print "focus out"

    def key_press(self, event):
        print "key %s press with code %s"%(event.char,event.keycode)

    def key_release(self, event):
        print "key %s release with code %s"%(event.char,event.keycode)

