from Tkinter import Tk, Frame, Menu, Canvas, Label
import tkFileDialog

class EditorRenderer:
    r"""
    Generic class which can draw to a TranslationSurfaceEditor
    THIS IS PRETTY ROUGH
    """
    def __init__(self, editor):
        r"""
        INPUT:
        - ``editor`` -- the TranslationSurfaceEditor we will draw on
        """
        self.editor = editor        

    def initial_render(self):
        r"""
        Called when the renderer initially becomes active.
        """
        pass

    def pause_render(self, event):
        r"""
        Called when the renderer temporarily becomes inactive.
        """
        pass    

    def restart_render(self, event):
        r"""
        Called when the renderer reactivates.
        """
        pass

    def end_render(self, event):
        r"""
        Called when the renderer permenantly stops.
        """
        pass
