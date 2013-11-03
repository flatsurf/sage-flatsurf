#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Current status of window for manipulating translation surfaces.
This file can be run from bash via "sage surface_manipulator.py"
"""

from Tkinter import Tk, Frame, Menu, Canvas, Label, Radiobutton, IntVar
import tkFileDialog
from editor_actor import *
from editor_renderer import *
from surface_bundle import *

class SurfaceManipulator(Frame):
    r"""
    A translation surface editor in tk.
    """
  
    def __init__(self, parent, surface=None, surfaces=[]):
        r"""
        INPUT:
        - ``surfaces`` -- a list of surfaces that the editor may modify
        - ``surface`` -- surface selected by default
        - ``parent`` -- parent Tk window
        """
        Frame.__init__(self, parent)   
        self.parent = parent
        self.pack(fill="both", expand=1)

        r"""Surface currently being manipulated"""
        self._surface=None
        r"""List of surfaces in editor"""
        self._surfaces=[]
        r"""More variables to initialize"""
        self._currentActor=None
        r"""Initialization of GUI"""
        self._init_menu()
        self._init_gui()
        r"""Setup surface list"""
        for s in surfaces:
            self.add_surface(s)
        r"""Setup initial surface"""
        if (surface != None):
            self.add_surface(surface)
        self.set_surface(surface)

    def _init_menu(self):
        
        menubar = Menu(self.parent)
        self.parent.config(menu=menubar)
        
        #new_menu = Menu(menubar, tearoff=0)
        #new_menu.add_command(label="Billiard Table", command=self.on_new_similarity_surface)
        
        file_menu = Menu(menubar, tearoff=0)
        #file_menu.add_cascade(label="New", menu=new_menu)
        file_menu.add_command(label="New Similarity Surface", command=self.on_new_similarity_surface)
        file_menu.add_separator()
        file_menu.add_command(label="About", command=self.on_about)
        file_menu.add_command(label="Export PostScript", command=self.on_export)
        file_menu.add_command(label="Exit", command=self.on_exit,accelerator="Alt+F4")
        menubar.add_cascade(label="File", underline=0, menu=file_menu)

        self._surface_menu = Menu(menubar, tearoff=0)
        self._selected_surface = IntVar()
        self._selected_surface.set(-1)
        self._surface_menu.add_radiobutton(label="None", 
            command=self.menu_select_surface, variable=self._selected_surface, 
            value=-1)
        menubar.add_cascade(label="Surface", underline=0, menu=self._surface_menu)

        self._create_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Create", underline=0, menu=self._create_menu)

        self._action_menu = Menu(menubar, tearoff=0)
        self._reset_action_menu()
        menubar.add_cascade(label="Action", underline=0, menu=self._action_menu)
    
        help_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)

    def _init_gui(self):
        self.parent.title("FlatSurf Editor")

        self._canvas = Canvas(self, bg="white", width=300, height=200)
        self._canvas.pack(fill="both", expand=1)
        
        self.bottom_text = Label(self, text="Welcome to FlatSurf.",
            anchor="w")
        self.bottom_text.pack(fill="x", expand=0)

        self.set_actor(None)

    def add_surface(self,newsurface):
        r"""
        Add a surface to the display list for the window.
        Returns the index of the new surface in the surface list.
        """
        if (newsurface==None):
            return -1
        i=0
        for s in self._surfaces:            
            if (s==newsurface):
                return i
            i=i+1
        self._surfaces.append(newsurface)
        self._surface_menu.add_radiobutton(label=newsurface.name,
            command=self.menu_select_surface, variable=self._selected_surface,
            value=len(self._surfaces)-1)
        return len(self._surfaces)-1

    def get_canvas(self):
        return self._canvas

    def get_center(self):
        r"""
        Return the center of the canvas as a pair of integers.
        """
        return ( self.get_width()/2, self.get_height()/2 )

    def get_height(self):
        r"""
        Return the height of the canvas (an integer).
        """
        return self.get_canvas().winfo_height()

    def get_surface_bundle(self):
        r"""
        Get the current surface bundle, or None if there is none.
        """
        return self._surface

    def get_width(self):
        r"""
        Return the width of the canvas (an integer).
        """
        return self.get_canvas().winfo_width()

    def menu_select_surface(self):
        r"""
        Called when a surface is selected from a menu.
        """
        i=self._selected_surface.get()
        if (i==-1):
            self.set_surface(None)
        else:
            self.set_surface(self._surfaces[i])

    def on_about(self):
        self.set_text("Written by Vincent Delecroix and Pat Hooper.")

    def on_delete_junk(self):
        self._canvas.delete("junk")

    def on_exit(self):
        self.quit()

    def on_export(self):
        r"""
        Export image as postscript file.
        """
        myFormats = [('PostScript','*.ps')]
        fileName = tkFileDialog.asksaveasfilename(parent=self,
            filetypes=myFormats , title="Save image as...")
        if len(fileName ) > 0:
            self._canvas.update() 
            self._canvas.postscript(file = fileName) 
            self.set_text("Wrote image to "+fileName)

    def on_new_similarity_surface(self):  
        s=CreateSimilaritySurfaceBundle(len(self._surfaces),self)
        if (s!=None):
            i=self.set_surface(s)
            self.set_text("Created new surface `"+self._surfaces[i].name+"'.")

    def _on_no_surface(self):
        self._canvas.delete("all")

    def _on_zoom(self):
        self.set_actor(ZoomActor(self))

    def _on_recenter(self):
        self.set_actor(RecenterActor(self))

    def _reset_action_menu(self):
        for i in range(100):
            self._action_menu.delete(0)
        self._action_menu.add_command(label="Recenter", underline=2, command=self._on_recenter)
        self._action_menu.add_command(label="Zoom", underline=0, command=self._on_zoom,accelerator="Alt+Z")
        self._action_menu.add_command(label="Delete Junk", command=self.on_delete_junk)

    def _reset_create_menu(self):
        for i in range(100):
            self._create_menu.delete(0)
        if self._surface!= None:
            self._surface.make_create_menu(self._create_menu)

    def set_text(self, text):
        self.bottom_text["text"]=text

    def set_actor(self, actor):
        r"""
        Set the current mode of user interaction.
        """
        if (actor != self._currentActor):
            if self._currentActor != None:
                self._currentActor.on_deactivate()
            if (actor==None):
                self.set_text("Nothing going on.")
                # Event bindings
                self._canvas.unbind('<Button-1>')
                self._canvas.unbind('<Button-2>')
                self._canvas.unbind('<Button-3>')
                self._canvas.unbind('<Double-Button-1>')
                self._canvas.unbind('<Shift-Button-1>')
                self._canvas.unbind('<Motion>')
                self.unbind('<FocusIn>')
                self.unbind('<FocusOut>')
                self.parent.unbind('<Key>')
                self.parent.unbind('<KeyRelease>')
            else:
                # Event bindings
                self._canvas.bind('<Button-1>', actor.single_left_click)
                self._canvas.bind('<Double-Button-1>', actor.double_left_click)
                self._canvas.bind('<Triple-Button-1>', actor.double_left_click)
                self._canvas.bind('<Button-2>', actor.single_middle_click)
                self._canvas.bind('<Double-Button-2>', actor.double_middle_click)
                self._canvas.bind('<Triple-Button-2>', actor.double_middle_click)
                self._canvas.bind('<Button-3>', actor.single_right_click)
                self._canvas.bind('<Double-Button-3>', actor.double_right_click)
                self._canvas.bind('<Triple-Button-3>', actor.double_right_click)
                self._canvas.bind('<Shift-Button-1>', actor.shift_click)
                self._canvas.bind('<Motion>', actor.mouse_moved)
                self.bind('<FocusIn>', actor.focus_in)
                self.bind('<FocusOut>', actor.focus_out)
                self.parent.bind('<Key>', actor.key_press)
                self.parent.bind('<KeyRelease>', actor.key_release)
                self._currentActor=actor
                self._currentActor.on_activate()

    def set_surface(self,surface):
        i=self.add_surface(surface)
        if (surface != self._surface):
            self._canvas.delete("all")
            self._surface=surface
            self._surface_menu.invoke(i+1)
            if (i>=0):
                self.set_text("Switched to `"+self._surface.name+"'.")
                self.parent.title(self._surface.name)
                self._reset_create_menu()
                if (isinstance(self._surface, EditorRenderer)):
                    self._surface.initial_render()
            else:
                self.set_text("No surface selected.")
                self.parent.title("FlatSurf Editor")
        return i

def main():
    root = Tk()
    root.geometry("400x300+300+300")
    app = SurfaceManipulator(root)
    root.mainloop()  

if __name__ == '__main__':
    main()  
