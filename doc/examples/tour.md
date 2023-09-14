---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.0
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# A Tour of the flatsurf Software Suite

+++

## Defining Surfaces
see also [this section](https://flatsurf.github.io/sage-flatsurf/examples/defining_surfaces.html) in the documentation.

+++

### Predefined Surfaces

```{code-cell} ipython3
from flatsurf import translation_surfaces, similarity_surfaces, dilation_surfaces

translation_surfaces.cathedral(1, 4).plot().show()
translation_surfaces.infinite_staircase().plot().show()
```

### Building Surfaces from Polygons

+++

Surfaces can be built by specifying a (finite) set of polygons and gluings between the sides of the polygons.

Once you call `set_immutable()`, the type of the surface is determined, here a translation surface:

```{code-cell} ipython3
from flatsurf import MutableOrientedSimilaritySurface, Polygon

hexagon = Polygon(vertices=((0, 0), (3, 0), (3, 1), (3, 2), (0, 2), (0, 1)))
# hexagon.plot().show()

square = Polygon(vertices=((0, 0), (1, 0), (1, 1), (0, 1)))
# square.plot().show()

S = MutableOrientedSimilaritySurface(QQ)
S.add_polygon(hexagon)
S.add_polygon(square)

S.glue((0, 0), (0, 3))
S.glue((0, 1), (1, 3))
S.glue((0, 2), (0, 4))
S.glue((0, 5), (1, 1))
S.glue((1, 0), (1, 2))
S.set_immutable()
print(S)
S.plot()
```

We can also create a half-translation surface:

```{code-cell} ipython3
T = MutableOrientedSimilaritySurface.from_surface(S)
T.glue((1, 1), (0, 2))
T.glue((0, 4), (0, 5))
T.set_immutable()
print(T)
T.plot().show()
```

Or anything that can be built by gluing with similarities… (though the labeling is not overly helpful currently)

```{code-cell} ipython3
T = MutableOrientedSimilaritySurface.from_surface(S)
T.glue((0, 0), (1, 1))
T.glue((0, 5), (0, 3))
T.set_immutable()
print(T)
T.plot().show()
```

Or a surface with boundary, with self-gluings, …

+++

### Billiards

+++

see below at "Base Rings".

+++

### Building Surfaces that are none of the Above

+++

There is a relatively small contract that a surface needs to implement, things such as "what is on the other side of this edge?", "is this surface compact?", … so it is not hard to implement surfaces from scratch. See the documentation.

+++

## Computing with Trajectories & Saddle Connections

```{code-cell} ipython3
from flatsurf import translation_surfaces

S = translation_surfaces.infinite_staircase()
v = S.tangent_bundle()(0, (1/47, 1/49), (1, 1/31)).straight_line_trajectory()
v.flow(100)
v.plot(color='red') + S.plot()
```

```{code-cell} ipython3
from flatsurf import translation_surfaces
S = translation_surfaces.octagon_and_squares()

connections = S.saddle_connections(squared_length_bound=100)
connections.sort(key=lambda c: c.length())

lengths = sorted(set(c.length() for c in connections))

def color(c):
    return colormaps.Accent(lengths.index(c.length()) / len(lengths))[:3]

S.plot(polygon_labels=False, edge_labels=False) + sum(sc.plot(color=color(sc)) for sc in connections)
```

## Base Rings

Most **exact** subrings of the reals from SageMath are supported, in particular $\mathbb{Q}$, NumberField and "exact reals".

```{code-cell} ipython3
from flatsurf import similarity_surfaces, Polygon
P = Polygon(angles=(3, 4, 13), vertices=[(0, 0), (1, 0)])
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")
print(f"{P.base_ring() = }")
```

The plotting in sage-flatsurf does not produce a very appealing picture of this surface yet. The ipyvue-flatsurf widget, if installed, has a more sophisticated layout algorithm for translation surfaces. You can see its output by running:

```python
from ipyvue_flatsurf import Widget
Widget(S).show()
```

```{code-cell} ipython3
from pyexactreal import ExactReals
R = ExactReals(P.base_ring())
almost_one = R.random_element(1)
print(f"{almost_one = }")

from flatsurf import similarity_surfaces, Polygon
Q = Polygon(angles=(3, 4, 13), vertices=[(0, 0), (almost_one, 0)])
S = similarity_surfaces.billiard(Q).minimal_cover(cover_type="translation")
```

## Flow Decompositions

+++

We compute flow decompositions on the unfolding of the (2, 2, 5) triangle. We only run 4 iterations of decomposition algorithm initially which does not always manage to find all the cylinders there are:

```{code-cell} ipython3
from flatsurf import similarity_surfaces, Polygon, GL2ROrbitClosure
P = Polygon(angles=(2, 2, 5))
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")

O = GL2ROrbitClosure(S)
decompositions = iter(O.decompositions(bound=64, limit=4))
```

If you run this once and then run the cell below, you get a nice picture with 4 cylinders.

Run it again and you see one cylinder and one undetermined component. That component is going to decompose into three cylinders eventually. However, these components cannot be plotted since they are extremely long; the Widget loads all the intersections of the cylinders' boundary with half edges…

Run it a few more times to see a situation with two minimal components. Unfortunately, the widget struggles a bit to display minimal components and produces some gap artifacts.

```{code-cell} ipython3
decomposition = next(decompositions)
print(f"{decomposition = }")
```

If you have the (optional) ipyvue-flatsurf widget installed, you can run the following cell to visualize the decomposition:

```python
from ipyvue_flatsurf import Widget
Widget(decomposition)
```

```{code-cell} ipython3
decomposition.decompose()
print(f"{decomposition = }")
```

Here is a somewhat mysterious case. On the same surface, trying to decompose into a direction of a saddle connection, we can never certify that the components is minimal. (We ran this for a very long time and there really don't seem to be any cylinders out there.)

```{code-cell} ipython3
from flatsurf import similarity_surfaces, Polygon, GL2ROrbitClosure
P = Polygon(angles=(2, 2, 5))
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")

c = P.base_ring().gen()
direction = (12*c^4 - 60*c^2 + 56, 16*c^5 - 104*c^3 + 132*c)

O = GL2ROrbitClosure(S)
D = O.decomposition(direction, limit=1024)
print(f"{D = }")

D = O.decomposition(direction, limit=0)
iet = D.components()[0R].intervalExchangeTransformation()
print(f"{iet = }")
```

## $SL_2(\mathbb{Z})$ orbits of square-tiled Surfaces

Square-tiled surfaces (also called origamis) are translation surfaces built from unit squares. Many properties of their $GL_2(\mathbb{R}$-orbit closures can be computed easily.

Note that the `Origami` object that we manipulate below are different from all other translation surfaces that were built above (that were using `sage-flatsurf`)

```{code-cell} ipython3
from surface_dynamics import Origami, OrigamiDatabase
```

```{code-cell} ipython3
o = Origami('(1,2)', '(1,3)')
```

```{code-cell} ipython3
o.plot()
```

```{code-cell} ipython3
T = o.teichmueller_curve()
```

```{code-cell} ipython3
T.sum_of_lyapunov_exponents()
```

```{code-cell} ipython3
o.lyapunov_exponents_approx()
```

```{code-cell} ipython3

```

```{code-cell} ipython3
# as with sage-flatsurf translation surfaces we can get the stratum
o.stratum()
```

```{code-cell} ipython3
# and we can even compute the stratum component
o.stratum_component()
```

```{code-cell} ipython3

```

```{code-cell} ipython3
# There is a relatively big database of pre-computed arithmetic Teichmüller curves that can be querried
D = OrigamiDatabase()
```

```{code-cell} ipython3
# get the list of properties that are stored for each Teichmüller curves
D.cols()
```

```{code-cell} ipython3
# get of summary of the content of the database
D.info()
```

```{code-cell} ipython3
# get the list of Teichmüller curves in genus gsuch that in any direction we have >= g cylinders 
for g in range(2, 7):
    q = D.query(('genus', '=', g), ('min_nb_of_cyls', '>=', g))
    print('g={}: got {} examples'.format(g, q.number_of()))
```

```{code-cell} ipython3
# get some more information than counting by setting columns in the query
g = 3
q = D.query(('genus', '=', g), ('min_nb_of_cyls', '>=', g))
q.cols('stratum')
print(set(q))
```

```{code-cell} ipython3
# the actual origami representatives of the Teichmüller curves are in the column "representative"
q.cols('representative')
o = choice(q.list())
print(o)
```

```{code-cell} ipython3
g = 2
q = D.query(('genus', '=', g), ('min_nb_of_cyls', '>=', g))
q.cols('stratum')
print(set(q))
```

```{code-cell} ipython3

```

## Combinatorial Graphs up to Isomorphism

```{code-cell} ipython3
from surface_dynamics import FatGraphs
```

```{code-cell} ipython3
# making the list of graphs in genus 2, 2 faces and vertex degree at least 3
fg = FatGraphs(g=2, nf=2, vertex_min_degree=3)
L = fg.list()
```

```{code-cell} ipython3
print(len(L))
```

```{code-cell} ipython3
fg = L[0]
```

```{code-cell} ipython3
# the display (and the encoding) uses the standard representation with permutations
print(fg)
```

```{code-cell} ipython3

```

## Hyperbolic geometry

```{code-cell} ipython3
from flatsurf import HyperbolicPlane
```

```{code-cell} ipython3
# the hyperbolic plane whose coordinates belong to QQ[2^(1/3)]
x = polygen(QQ)
K = NumberField(x^3 - 2, 'cbrt3', embedding=AA(2)**(1/3))
cbrt3 = K.gen()
H2 = HyperbolicPlane(K)
```

```{code-cell} ipython3
g0 = H2.geodesic(0, 1)
g1 = H2.geodesic(1/2 - cbrt3/5, 2)
```

```{code-cell} ipython3
g0.plot(color='blue') + g1.plot(color='red')
```

```{code-cell} ipython3
g0.plot('klein', color='blue') + g1.plot('klein', color='red')
```

```{code-cell} ipython3
p = g0.intersection(g1)
```

```{code-cell} ipython3
p.coordinates('klein')
```

```{code-cell} ipython3
# computing coordinates in the upper half plane requires square-root so this would fail
# p.coordinates('half_plane')

p.change_ring(AA).coordinates('half_plane')
```

```{code-cell} ipython3

```

```{code-cell} ipython3
g2 = H2.geodesic(p, 1/2)
```

```{code-cell} ipython3
g0.plot(color='blue') + g1.plot(color='red') + g2.plot(color='orange')
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

## Veech Surfaces & Iso-Delaunay Tessellations

(for Iso-Delaunay tessellation, see [sage-flatsurf#163](https://github.com/flatsurf/sage-flatsurf/pull/163))

```{code-cell} ipython3
from flatsurf import polygons, similarity_surfaces, GL2ROrbitClosure
```

```{code-cell} ipython3
T = polygons.triangle(1, 4, 7)
S = similarity_surfaces.billiard(T).minimal_cover('translation')
S = S.erase_marked_points()
O = GL2ROrbitClosure(S)
```

```{code-cell} ipython3
O.decompositions(10)
```

```{code-cell} ipython3
d = next(O.decompositions(10))
```

```{code-cell} ipython3
# the direction is completely periodic
d
```

```{code-cell} ipython3
# one can check for parabolicity... though we get a tribool
d.parabolic()
```

```{code-cell} ipython3
bool(d.parabolic())
```

## Strata & Orbit Closures

```{code-cell} ipython3
from flatsurf import similarity_surfaces, polygons
from surface_dynamics import AbelianStratum
```

```{code-cell} ipython3
T = polygons.triangle(3, 4, 13)
S = similarity_surfaces.billiard(T).minimal_cover('translation')
H = S.stratum()
```

```{code-cell} ipython3
print(H)
```

```{code-cell} ipython3
S0 = S.erase_marked_points()
H0 = S0.stratum()
print(H0)
```

```{code-cell} ipython3
print(H.dimension(), H0.dimension())
```

```{code-cell} ipython3

```

## Veering Triangulations

- [M. Bell, V. Delecroix, V. Gadre, R. Gutiérrez-Romo, Saul Schleimer  arXiv:1909.00890](https://arxiv.org/abs/1909.00890)
- [B. Zykoski arXiv:2206.04143](https://arxiv.org/abs/2206.04143)

```{code-cell} ipython3
# not by default in the flatsurf stack... but soon
%pip install git+https://github.com/flatsurf/veerer
```

```{code-cell} ipython3
from veerer import VeeringTriangulation
from surface_dynamics import AbelianStratum
```

```{code-cell} ipython3
H2 = AbelianStratum(2)
vt = VeeringTriangulation.from_stratum(H2)
```

```{code-cell} ipython3
vt
```

```{code-cell} ipython3
vt.flat_structure_middle().plot()
```

```{code-cell} ipython3

```
