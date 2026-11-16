---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.1
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# Defining Surfaces

sage-flatsurf handles flat and hyperbolic surfaces. Both are built the same
way, either from collections of pre-built surfaces or by gluing polygons, and
they share many methods. A surface with infinitely many polygons is written as
a class that computes its polygons and their gluings. This is explained, for a
flat and for a hyperbolic example, in
[Building finite and infinite surfaces as classes](./surfaces_as_classes).

## Built in surfaces

There are four built-in collections of flat surface
- `similarity_surfaces`
- `dilation_surfaces`
- `half_translation_surfaces`
- `translation_surfaces`

Veech's double n-gon surfaces:

```{code-cell}
from flatsurf import translation_surfaces

s = translation_surfaces.veech_double_n_gon(5)
s.plot()
```

The Arnoux-Yoccoz surface of arbitrary genus:

```{code-cell}
s = translation_surfaces.arnoux_yoccoz(3)
s.plot()
```

A dilation square:

```{code-cell}
from flatsurf import dilation_surfaces

ds = dilation_surfaces.genus_two_square(1 / 2, 1 / 3, 1 / 4, 1 / 5)
ds.plot(polygon_labels=False, edge_labels="letter")
```

Chamanara's infinite translation surface:

```{code-cell}
s = translation_surfaces.chamanara(1 / 2)
```

```{code-cell}
s.plot(polygon_labels=False, edge_labels=False)
```

The infinite stairacase (which is an infinite ramified cover of the torus):

```{code-cell}
s = translation_surfaces.infinite_staircase()
```

```{code-cell}
s.plot()
```

Hyperbolic surfaces have their own collection, `hyperbolic_surfaces`. The
regular octagon with angles $\pi/4$, its opposite sides glued, is a closed
surface of genus 2:

```{code-cell}
from flatsurf import hyperbolic_surfaces

S = hyperbolic_surfaces.regular_octagon()
S
```

A pair of pants has three geodesic boundary components. Each is given by the
trace $t = 1 + 2 \cosh \ell$ of its holonomy in $\mathrm{SO}(2, 1)$, where
$\ell$ is its length:

```{code-cell}
S = hyperbolic_surfaces.pants(boundary_traces=(15, 15, 15))
S
```

What is specific to the hyperbolic setting is described in
[Hyperbolic Surfaces](./hyperbolic_surfaces).

## Billiard tables and pillowcases

The pillowcase of a polygon is the polygon glued to its mirror image along its
edges. The same call builds it in both geometries. A flat triangle with angles
proportional to $1, 2, 3$, that is $\pi/6, \pi/3, \pi/2$, gives a flat sphere;
`euclidean_polygons.triangle(1, 2, 3)` takes the angles proportional to
`1, 2, 3`:

```{code-cell}
from flatsurf import similarity_surfaces, euclidean_polygons

F = similarity_surfaces.pillowcase(euclidean_polygons.triangle(1, 2, 3))
F
```

A hyperbolic triangle with angles $\pi/3, \pi/4, \pi/4$ gives a hyperbolic
sphere. The hyperbolic triangle is built the same way, from its angles, which
are multiples of $2\pi$:

```{code-cell}
from sage.all import QQ

from flatsurf import hyperbolic_polygons, hyperbolic_surfaces

T = hyperbolic_polygons.triangle(QQ(1) / 6, QQ(1) / 8, QQ(1) / 8)
X = hyperbolic_surfaces.pillowcase(T)
X
```

## Minimal translation surface covers

The minimal translation cover of the flat pillowcase `F` built above is a
translation surface:

```{code-cell}
ss = F.minimal_cover(cover_type="translation")
```

```{code-cell}
ss
```

```{code-cell}
gs = ss.graphical_surface()
```

```{code-cell}
gs.make_all_visible(limit=12)
```

```{code-cell}
gs.plot()
```

## Building surfaces from polygons

This defines a regular 12-gon with algebraic real coordinates (AA) with first
vector given by (1,0):

```{code-cell}
from flatsurf import polygons

p0 = polygons.regular_ngon(12, field=AA)
p1 = polygons.regular_ngon(3, field=AA)
```

```{code-cell}
p0.plot() + p1.plot()
```

The vertices of n-gons are numbered by $\{0,...,n-1\}$, with the $0$-th vertex
at the origin. Edge $i$ joins vertex $i$ to vertex $i+1 \pmod{n}$.

We can act on polygon with $2 \times 2$ matrices. We define the rotation by
$\frac{\pi}{6}$ below:

```{code-cell}
R = matrix(AA, [[cos(pi / 6), -sin(pi / 6)], [sin(pi / 6), cos(pi / 6)]])
show(R)
```

```{code-cell}
R * p1
```

Define a surface over the field <code>AA</code> of algebraic reals.

```{code-cell}
from flatsurf import MutableOrientedSimilaritySurface

surface = MutableOrientedSimilaritySurface(AA)
```

Add two polygons to the surface with labels 0 and 1:

```{code-cell}
surface.add_polygon(p0, label=0)
```

```{code-cell}
surface.add_polygon(p1, label=1)
```

Glue the edges of polygon 0 to the parallel edges of polygon 1.

```{code-cell}
surface.glue((0, 6), (1, 0))
surface.glue((0, 10), (1, 1))
surface.glue((0, 2), (1, 2))
```

Add three more rotated triangles and glue them appropriately.

```{code-cell}
for i in range(1, 4):
    surface.add_polygon((R**i) * p1, label=i + 1)
    surface.glue((0, 6 + i), (i + 1, 0))
    surface.glue((0, (10 + i) % 12), (i + 1, 1))
    surface.glue((0, 2 + i), (i + 1, 2))
```

Now we have a closed surface. In fact this is a translation surface:

```{code-cell}
surface
```

Once we are done building the surface, it is recommended to make the surface
immutable. This lets sage-flatsurf figure out of which nature this surface is,
e.g., that it is a translation surface. This speeds up many operations on the
surface and makes it possible to compute things that are only defined or
implemented for some types of surfaces:

```{code-cell}
surface.set_immutable()
surface
```

If you want to compute things, such as the stratum without making a surface
immutable, please refer to the details in the documentation of the
``flatsurf.geometry.categories`` module in the module reference.

We can plot the surface. Edges are labeled according to the polygon they are
glued to.

```{code-cell}
surface.plot()
```

The field containing the vertices:

```{code-cell}
surface.base_ring()
```

Computations in the Algebraic Real Field (AA) are slow. It is better to use a
NumberField. The following finds a smaller number field::

```{code-cell}
points = [surface.polygon(p).vertex(v) for (p, v) in surface.edges()]
coordinates = [point[0] for point in points] + [point[1] for point in points]
base_ring = Sequence(
    [coordinate.as_number_field_element()[1] for coordinate in coordinates]
).universe()
ss = surface.change_ring(base_ring)
```

```{code-cell}
ss.base_ring()
```

Hyperbolic surfaces are glued from polygons of the hyperbolic plane, whose
points are given here by their coordinates in the Klein model, the unit disk
in which hyperbolic geodesics are straight segments. The quadrilateral whose
vertices are those of a Euclidean square in the Klein model has two pairs of
sides of equal lengths, as the hyperbolic cosines squared of its sides show:

```{code-cell}
from sage.all import QQ

from flatsurf import HyperbolicPlane

H = HyperbolicPlane(QQ)
Q = H.polygon(
    vertices=[
        H.point(0, 0, model="klein"),
        H.point(QQ(1) / 2, 0, model="klein"),
        H.point(QQ(1) / 2, QQ(1) / 2, model="klein"),
        H.point(0, QQ(1) / 2, model="klein"),
    ]
)
[Q.edge_cosh_squared(i) for i in range(4)]
```

The surface is built with the same `add_polygon` and `glue` as the flat
surfaces above; every gluing is checked to be an isometry defined over the
base ring of the surface:

```{code-cell}
from flatsurf import MutableOrientedHyperbolicSurface

S = MutableOrientedHyperbolicSurface(QQ)
S.add_polygon(Q)
S.glue((0, 0), (0, 3))
```

Edges of different lengths cannot be glued: executing
`S.glue((0, 0), (0, 1))` would result in an error.

Gluing the two remaining edges and making the surface immutable with
`set_immutable`, as for the flat surfaces, gives a sphere with three cone
points:

```{code-cell}
S.glue((0, 1), (0, 2))
S.set_immutable()
S
```

```{code-cell}
[S(0, i).angle(numerical=True) for i in range(3)]
```

## Basic geometric invariants

Flat and hyperbolic surfaces share many methods. Here they are compared on the
flat pillowcase `F` and the hyperbolic pillowcase `X` built above. Angles are
multiples of $2\pi$ in both geometries, so that a right angle is `1/4`:

```{code-cell}
F.angles(), X.angles()
```

The flat area is the plain Euclidean area, while the hyperbolic area is in
units of $2\pi$:

```{code-cell}
F.area(), X.area()
```

The signature records the genus and the orders $k = a - 1$ of the cone points
of angles $2\pi a$:

```{code-cell}
F.signature(), X.signature()
```

Each cone point of these two spheres has an angle $2\pi/m$, so that both
surfaces are orbifolds, with orbifold points of orders $m$:

```{code-cell}
F.is_orbifold(), X.is_orbifold()
```

```{code-cell}
F.orbifold_points()
```

```{code-cell}
X.orbifold_points()
```

The holonomy representation sends a closed path to the isometry obtained by
developing the polygons along it, in the group of similarities of the plane
for a flat surface. The image for `F` is generated by a rotation of order 6
and a half-turn: it is the wallpaper group $(0; 2, 3, 6)$, which uniformises
the orbifold `F`:

```{code-cell}
G = F.holonomy_representation().image()
g, h = G.gens()
[(g**k).is_one() for k in range(1, 7)], [(h**k).is_one() for k in range(1, 3)]
```

For a hyperbolic surface, the holonomy lands in $\mathrm{SO}(2, 1)^+$. The
image for `X` is generated by rotations of orders 3 and 4: it is the triangle
group $(3, 4, 4)$, which uniformises the orbifold `X`:

```{code-cell}
G = X.holonomy_representation().image()
g, h = G.gens()
[(g**k).is_one() for k in range(1, 4)], [(h**k).is_one() for k in range(1, 5)]
```

Signatures, and the construction of a hyperbolic surface from its signature
with `hyperbolic_surfaces.from_signature`, are described in
[Hyperbolic Surfaces](./hyperbolic_surfaces).

## Getting a surface from Flipper

The constructions that follow are for flat surfaces only for now.

<a href="http://flipper.readthedocs.io/en/latest/">Flipper</a> is a program
written by Mark Bell which understands mapping classes and can compute the
flat structure associated to a pseudo-Anosov mapping class. FlatSurf can
import this structure.

This code below requires flipper to be installed. You can do this by running
the shell within sage:
<code>sage --sh</code>
Then within the shell execute:
<code>python -m pip install flipper --user --upgrade</code>
More information including pitfalls are described in
<a href="http://flipper.readthedocs.io/en/latest/start.html#installation">Flipper's installation instructions</a>.

```{code-cell}
import flipper
```

```{code-cell}
T = flipper.load("SB_4")
```

```{code-cell}
h = T.mapping_class("s_0S_1s_2S_3s_1S_2")
```

```{code-cell}
h.is_pseudo_anosov()
```

```{code-cell}
s = translation_surfaces.from_flipper(h)
```

The surface s is actually a half translation surface

```{code-cell}
s
```

```{code-cell}
s.plot()
```

## Flat surface of a polyhedron

```{code-cell}
from flatsurf.geometry.polyhedra import platonic_dodecahedron

polyhedron, s, mapping = platonic_dodecahedron()
```

The surface $s$ is a Euclidean cone surface.

```{code-cell}
s
```

```{code-cell}
s.plot()
```

Sage has a built in polyhedron class. You can build a polyhedron as a convex
hull of a list of vertices.

```{code-cell}
polyhedron = Polyhedron([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
```

```{code-cell}
polyhedron.plot()
```

The following computes the boundary surface as a Euclidean cone surface. It
also provides a map from the surface to the polyhedron.

```{code-cell}
from flatsurf.geometry.polyhedra import polyhedron_to_cone_surface

s, mapping = polyhedron_to_cone_surface(polyhedron)
s
```

```{code-cell}
s.plot()
```
