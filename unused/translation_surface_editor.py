import Tkinter

class TranslationSurfaceEditor:
    r"""
    A translation surface editor in tk.
    """
    def __init__(self, surface=None):
        r"""
        INPUT:

        - ``surface`` -- a translation surface that the editor may modify
        """
        self.window = Tkinter.Tk()

        if surface is None:
            surface = FlatSurface()
        self.surface = surface

        # put a frame and a canvas
        self.frame = Tkinter.Frame(self.window,
                borderwidth=0,
                relief=Tkinter.FLAT,
                background='#dcecff')
        self.canvas = Tkinter.Canvas(self.frame,
                bg='#dcecff',
                width=500,
                height=500,
                highlightthickness=0)
        self.frame.pack(padx=0, pady=0, fill=Tkinter.BOTH, expand=Tkinter.YES)
        self.canvas.pack(padx=0, pady=0, fill=Tkinter.BOTH, expand=Tkinter.YES)

        # Event bindings
        self.canvas.bind('<Button-1>', self.single_left_click)
        self.canvas.bind('<Button-2>', self.single_middle_click)
        self.canvas.bind('<Button-3>', self.single_right_click)
        self.canvas.bind('<Double-Button-1>', self.double_click)
        self.canvas.bind('<Shift-Button-1>', self.shift_click)
        self.canvas.bind('<Motion>', self.mouse_moved)
        self.window.bind('<FocusIn>', self.focus_in)
        self.window.bind('<FocusOut>', self.focus_out)
        self.window.bind('<Key>', self.key_press)
        self.window.bind('<KeyRelease>', self.key_release)
        self.window.protocol("WM_DELETE_WINDOW", self.done)

    def single_left_click(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        print("single left click at (%s,%s)"%(x,y))
        self.canvas.create_rectangle(x-5,y-5,x+5,y+5, fill="#eeff00")

    def single_middle_click(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        print("single middle click at (%s,%s)"%(x,y))
        self.canvas.create_rectangle(x-5,y-5,x+5,y+5, fill='#00aa0f')

    def single_right_click(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        print("single right click at (%s,%s)"%(x,y))
        self.canvas.create_rectangle(x-5,y-5,x+5,y+5, fill='#ff00ff')

    def double_click(self, event):
        print("double click")

    def shift_click(self, event):
        print("shift click")

    def mouse_moved(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        print("mouse moved to (%s,%s)"%(x,y))

    def focus_in(self, event):
        print("focus in")

    def focus_out(self, event):
        print("focus out")

    def key_press(self, event):
        print("key %s press with code %s"%(event.char,event.keycode))

    def key_release(self, event):
        print("key %s release with code %s"%(event.char,event.keycode))

    def done(self, event=None):
        self.surface.data = 'nicer data'
        self.window.destroy()



