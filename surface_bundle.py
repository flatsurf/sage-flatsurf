from editor_actor import *
from editor_renderer import *


class SurfaceBundle:
    def __init__(self, name, editor):
        self.name = name
        self.editor=editor
    def makeCreateMenu(self,menu):
        pass

class BilliardTableBundle(SurfaceBundle,EditorRenderer):
    def __init__(self, n, editor):
        SurfaceBundle.__init__(self,"New Billiard Table "+str(n), editor)
        EditorRenderer.__init__(self,editor)
    def makeCreateMenu(self,menu):
        pass
    def initial_render(self):
        r"""
        Called when the renderer initially becomes active.
        """
        self.editor.canvas.create_rectangle(50, 25, 150, 75, fill="red")

