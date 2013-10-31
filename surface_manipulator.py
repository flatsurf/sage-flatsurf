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
        self.surface=None
        r"""List of surfaces in editor"""
        self.surfaces=[]
        r"""Initialization of Interfaces"""
        self.initMenu()
        self.initUI()
        r"""Setup surface list"""
        for s in surfaces:
            self.addSurface(s)
        r"""Setup initial surface"""
        if (surface != None):
            self.addSurface(surface)
        self.setSurface(surface)
    
    def initMenu(self):
        
        menubar = Menu(self.parent)
        self.parent.config(menu=menubar)
        
        newMenu = Menu(menubar, tearoff=0)
        newMenu.add_command(label="Billiard Table", command=self.onNewBilliardTable)
        
        fileMenu = Menu(menubar, tearoff=0)
        fileMenu.add_cascade(label="New", menu=newMenu)
        fileMenu.add_command(label="About", command=self.onAbout)    
        fileMenu.add_command(label="Export PostScript", command=self.onExport)    
        fileMenu.add_command(label="Exit", command=self.onExit)
        menubar.add_cascade(label="File", menu=fileMenu)

        self.surfaceMenu = Menu(menubar, tearoff=0)
        self.selectedSurface = IntVar()
        self.selectedSurface.set(-1)
        self.surfaceMenu.add_radiobutton(label="None", 
            command=self.menuSelectSurface, variable=self.selectedSurface, 
            value=-1)
        menubar.add_cascade(label="Surface", menu=self.surfaceMenu)

        self.createMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Create", menu=self.createMenu)

        self.actionMenu = Menu(menubar, tearoff=0)
        self.actionMenu.add_command(label="Delete Junk", command=self.onDeleteJunk)
        menubar.add_cascade(label="Action", menu=self.actionMenu)
    
        helpMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=helpMenu)

    def initUI(self):
        self.parent.title("FlatSurf Editor")

        self.canvas = Canvas(self, bg="white", width=300, height=200)
        self.canvas.pack(fill="both", expand=1)

        self.canvas.create_line(0, 0, 200, 100)
        self.canvas.create_line(0, 100, 200, 0, fill="red", dash=(4, 4))
        self.canvas.create_rectangle(50, 25, 150, 75, fill="blue")
        
        self.bottom_text = Label(self, text="Welcome to FlatSurf.",
            anchor="w")
        self.bottom_text.pack(fill="x", expand=0)

        self.setCurrentActor(DemoActor(self))

    def addSurface(self,newsurface):
        r"""
        Add a surface to the display list for the window.
        Returns the index of the new surface in the surface list.
        """
        if (newsurface==None):
            return -1
        i=0
        for s in self.surfaces:            
            if (s==newsurface):
                return i
            i=i+1
        self.surfaces.append(newsurface)
        self.surfaceMenu.add_radiobutton(label=newsurface.name,
            command=self.menuSelectSurface, variable=self.selectedSurface,
            value=len(self.surfaces)-1)
        return len(self.surfaces)-1

    def menuSelectSurface(self):
        r"""
        Called when a surface is selected from a menu.
        """
        i=self.selectedSurface.get()
        if (i==-1):
            self.setSurface(None)
        else:
            self.setSurface(self.surfaces[i])

    def onAbout(self):
        self.setBottomText("Written by Vincent Delecroix and Pat Hooper.")

    def onDeleteJunk(self):
        self.canvas.delete("junk")

    def onExit(self):
        self.quit()

    def onExport(self):
        r"""
        Export image as postscript file.
        """
        myFormats = [('PostScript','*.ps')]
        fileName = tkFileDialog.asksaveasfilename(parent=self,
            filetypes=myFormats , title="Save image as...")
        if len(fileName ) > 0:
            self.canvas.update() 
            self.canvas.postscript(file = fileName) 
            self.setBottomText("Wrote image to "+fileName)

    def onNewBilliardTable(self):  
        s=BilliardTableBundle(len(self.surfaces),self)
        if (s!=None):
            i=self.setSurface(s)
            self.setBottomText("Created new surface `"+self.surfaces[i].name+"'.")

    def onNoSurface(self):
        self.canvas.delete("all")

    def setBottomText(self, text):
        self.bottom_text["text"]=text

    def setCurrentActor(self, actor):
        r"""
        Set the current mode of user interaction.
        """
        self.currentAction=actor
        # Event bindings
        self.canvas.bind('<Button-1>', actor.single_left_click)
        self.canvas.bind('<Button-2>', actor.single_middle_click)
        self.canvas.bind('<Button-3>', actor.single_right_click)
        self.canvas.bind('<Double-Button-1>', actor.double_click)
        self.canvas.bind('<Shift-Button-1>', actor.shift_click)
        self.canvas.bind('<Motion>', actor.mouse_moved)
        self.bind('<FocusIn>', actor.focus_in)
        self.bind('<FocusOut>', actor.focus_out)
        self.bind('<Key>', actor.key_press)
        self.bind('<KeyRelease>', actor.key_release)

    def setSurface(self,surface):
        i=self.addSurface(surface)
        if (surface != self.surface):
            self.canvas.delete("all")
            self.surface=surface
            self.surfaceMenu.invoke(i+1)
            if (i>=0):
                self.setBottomText("Switched to `"+self.surface.name+"'.")
                self.parent.title(self.surface.name)
                if (isinstance(self.surface, EditorRenderer)):
                    self.surface.initial_render()
            else:
                self.setBottomText("No surface selected.")
                self.parent.title("FlatSurf Editor")
        return i

def main():
    root = Tk()
    root.geometry("400x300+300+300")
    app = SurfaceManipulator(root)
    root.mainloop()  

if __name__ == '__main__':
    main()  
