r"""
Translation surface in Sage.

EXAMPLES:

Right now we can do from sage::

    sage: %runfile translation_surface.py
    sage: t = TranslationSurface()
    sage: t.edit()
"""

class TranslationSurfaceDisplay:
    r"""
    Specify position of polygons in the plane.

    Internally it is a mapping from the index set of polygon to K^2 zhere K is a
    field embedded in RR^2.

    TODO:

        Think about a class for plotting saddle connections, cylinders, circles.
    """


class TranslationSurface:
    r"""
    A surface with a flat metric and conical singularities (not necessarily
    multiple angle of pi or 2pi).

    - polygon = polygon + vertex (or equivalently, canonical ordering)

    A translation surface is:

    - field embedded in R
    - index set for the (convex) polygons + favorite polygon
    - edges: ((t1,e1),(t2,e2))

    For finite case:

    - canonical labelings of polygons
    - Delaunay triangulation
    """
    data = None  # some data that can be modified in the .edit option

    def edit(self):
        r"""
        Launch the tk editor to interactively modify ``self``.
        """
        from translation_surface_editor import TranslationSurfaceEditor
        fse = TranslationSurfaceEditor(self)
        fse.window.mainloop()
