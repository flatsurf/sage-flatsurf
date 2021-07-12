---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: SageMath 9.2
  language: sage
  name: sagemath
---

# Defining Surfaces in FlatSurf

Initial version by Pat Hooper <whooper@ccny.cuny.edu>, Dec 16, 2017.

```{code-cell} ipython3
from flatsurf import *
```

## Built in surfaces

Veech's double n-gon surfaces:

```{code-cell} ipython3
s = translation_surfaces.veech_double_n_gon(5)
s.plot()
```

The Arnoux-Yoccoz surface of arbitrary genus is built in:

```{code-cell} ipython3
s=translation_surfaces.arnoux_yoccoz(3)
s.plot()
```

Chamanara's infinite translation surface:

```{code-cell} ipython3
s=translation_surfaces.chamanara(1/2)
```

```{code-cell} ipython3
s.plot(polygon_labels=False,edge_labels=False)
```

```{code-cell} ipython3
s=translation_surfaces.infinite_staircase()
```

```{code-cell} ipython3
s.plot()
```

## Billiard tables

```{code-cell} ipython3
s=similarity_surfaces.billiard(polygons(vertices=[(0,0), (3,0), (0,4)]))
```

```{code-cell} ipython3
s.plot()
```

## Minimal translation surface covers

Continuing the billiard example above, we get an infinite translation surface below:

```{code-cell} ipython3
ss = s.minimal_cover(cover_type="translation").copy(relabel=True)
```

```{code-cell} ipython3
gs = ss.graphical_surface()
```

```{code-cell} ipython3
gs.make_all_visible(limit=12)
```

```{code-cell} ipython3
gs.plot()
```

## Building surfaces from polygons

This defines a regular 12-gon with algebraic real coordinates (AA) with first vector given by (1,0):

```{code-cell} ipython3
p0 = polygons.regular_ngon(12,field=AA)
p1 = polygons.regular_ngon(3,field=AA)
```

```{code-cell} ipython3
p0.plot()+p1.plot()
```

The vertices of n-gons are numbered by $\{0,...,n-1\}$, with the $0$-th vertex at the origin. Edge $i$ joins vertex $i$ to vertex $i+1 \pmod{n}$.

We can act on polygon with $2 \times 2$ matrices. We define the rotation by $\frac{\pi}{6}$ below:

```{code-cell} ipython3
R = matrix(AA,[[cos(pi/6),-sin(pi/6)],[sin(pi/6),cos(pi/6)]])
show(R)
```

```{code-cell} ipython3
R*p1
```

Define a surface over the field <code>AA</code> of algebraic reals.

```{code-cell} ipython3
surface = Surface_dict(base_ring=AA)
```

Add two polygons to the surface with labels 0 and 1:

```{code-cell} ipython3
surface.add_polygon(p0,label=0)
```

```{code-cell} ipython3
surface.add_polygon(p1,label=1)
```

Set the "base label" for the surface. This is just a choice of a favorite polygon label.

```{code-cell} ipython3
surface.change_base_label(0)
```

Glue the edges of polygon 0 to the parallel edges of polygon 1.

```{code-cell} ipython3
surface.change_edge_gluing(0,6,1,0)
surface.change_edge_gluing(0,10,1,1)
surface.change_edge_gluing(0,2,1,2)
```

Add three more rotated triangles and glue them appropriately.

```{code-cell} ipython3
for i in range(1,4):
    surface.add_polygon((R**i)*p1,label=i+1)
    surface.change_edge_gluing(0,6+i,i+1,0)
    surface.change_edge_gluing(0,(10+i)%12,i+1,1)
    surface.change_edge_gluing(0,2+i,i+1,2)
```

Now we have a closed surface. In fact we have defined a Translation Surface. The package also supports
SimilaritySurface, ConeSurface, HalfDilationSurface, DilationSurface, and HalfTranslationSurface.

```{code-cell} ipython3
s=TranslationSurface(surface)
```

Test to insure that we created the translation surface correctly. (Errors would be printed if you did not glue parallel edges, or have some unglued edges, etc.)

```{code-cell} ipython3
TestSuite(s).run(verbose=True)
```

We can plot the surface. Edges are labeled according to the polygon they are glued to.

```{code-cell} ipython3
s.plot()
```

The field containing the vertices:

```{code-cell} ipython3
s.base_ring()
```

Computations in the Algebraic Real Field (AA) are slow. It is better to use a NumberField. The following finds the smallest embedding into a NumberField:

```{code-cell} ipython3
ss=s.copy(optimal_number_field=True)
```

```{code-cell} ipython3
ss.base_ring()
```

## Getting a surface from Flipper

<span style="color: red; font-weight: bold; ">This does not work as of SageMath 9.0. Code is commented out below.</span>

<a href="http://flipper.readthedocs.io/en/latest/">Flipper</a> is a program written by Mark Bell which understands mapping classes and can compute the flat structure associated to a pseudo-Anosov mapping class. FlatSurf can import this structure.

This code below requires flipper to be installed. You can do this by running the shell within sage:
<code>sage --sh</code>
Then within the shell execute:
<code>python -m pip install flipper --user --upgrade</code>
More information including pitfalls are described in <a href="http://flipper.readthedocs.io/en/latest/start.html#installation">Flipper's installation instructions</a>.

```{code-cell} ipython3
# import flipper
```

```{code-cell} ipython3
# T = flipper.load('SB_4')
```

```{code-cell} ipython3
# h = T.mapping_class('s_0S_1s_2S_3s_1S_2') 
```

```{code-cell} ipython3
# h.is_pseudo_anosov()
```

```{code-cell} ipython3
# s = translation_surfaces.from_flipper(h)
```

The surface s is actually a half translation surface

```{code-cell} ipython3
# type(s)
```

```{code-cell} ipython3
# s.plot()
```

## From polyhedra

```{code-cell} ipython3
from flatsurf.geometry.polyhedra import *
```

```{code-cell} ipython3
polyhedron,s,mapping = platonic_dodecahedron()
```

The surface $s$ is a Euclidean cone surface.

```{code-cell} ipython3
type(s)
```

```{code-cell} ipython3
s.plot()
```

Sage has a built in polyhedron class. You can build a polyhedron as a convex hull of a list of vertices.

```{code-cell} ipython3
polyhedron=Polyhedron([(0,0,0),(1,0,0),(0,1,0),(0,0,1)])
```

```{code-cell} ipython3
polyhedron.plot()
```

The following computes the boundary surface as a Euclidean cone surface. It also provides a map from the surface to the polyhedron.

```{code-cell} ipython3
s,mapping = polyhedron_to_cone_surface(polyhedron)
```

```{code-cell} ipython3
s.plot()
```

## Defining an infinite surface from scratch

The following demonstrates the implementation of a TranslationSurface. Each geometric structure has an underlying "Surface". The following defines a surface and then uses it to construct a translation surface.

```{code-cell} ipython3
from flatsurf.geometry.surface import Surface

class ParabolaSurface(Surface):
    def __init__(self):
        # The space of polygons with vertices in the rationals:
        self._P = Polygons(QQ)
        self._inv = matrix(QQ,[[-1,0],[0,-1]])
        
        # Set the base field to QQ, the base label to be 1, and note that the surface is infinite.
        Surface.__init__(self, QQ, ZZ(1), finite=False)
    
    def polygon(self, label):
        if label not in ZZ:
            raise ValueError("invalid label {!r}".format(lab))
        assert label != 0, "Label should not be zero."
        if label >= 0:
            if label==1:
                return self._P(vertices=[(0,0),(1,1),(-1,1)])
            else:
                return self._P( vertices=[
                    (label-1, (label-1)**2),
                    (label, label**2),
                    (-label, label**2),
                    (-label+1, (label-1)**2) ] )
        else:
            return self._inv*self.polygon(-label)

    def opposite_edge(self, label, e):
        if label not in ZZ:
            raise ValueError("invalid label {!r}".format(lab))
        assert label != 0, "Label should not be zero."

        if label==1 or label==-1:
            if e==1:
                return 2*label,3
            else:
                return -label,e
        else:
            if e==0 or e==2:
                return -label,e
            if e==1:
                if label>0:
                    return label+1,3
                else:
                    return label-1,3
            if e==3:
                if label>0:
                    return label-1,1
                else:
                    return label+1,1
```

```{code-cell} ipython3
s = TranslationSurface(ParabolaSurface())
```

```{code-cell} ipython3
TestSuite(s).run(verbose=True, skip="_test_pickling")
```

A graphical surface controls the display of graphical data. For an infinite surface you need to configure the display manually.

```{code-cell} ipython3
gs=s.graphical_surface()
```

We make six polygons nearest to the polygon with the base label visible.

```{code-cell} ipython3
gs.make_all_visible(limit=6)
```

```{code-cell} ipython3
s.plot()
```
