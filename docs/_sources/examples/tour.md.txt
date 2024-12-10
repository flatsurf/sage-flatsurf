---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: SageMath 10.2
  language: sage
  name: sagemath
---

# A Tour of the flatsurf Suite

The [flatsurf software suite](https://flatsurf.github.io) is a collection of mathematical software libraries to study translation surfaces and related objects.

The [SageMath](https://sagemath.org) library [``sage-flatsurf``](https://github.com/flatsurf/sage-flatsurf) provides the most convenient interface to the flatsurf suite. Here, we showcase some of its major features.

+++

## Defining Surfaces

Many surfaces that have been studied in the literature are readily available in sage-flatsurf.

You can find predefined translation surfaces in the collection {class}`translation_surfaces <flatsurf.geometry.similarity_surface_generators.TranslationSurfaceGenerators>`.

```{code-cell}
from flatsurf import translation_surfaces

S = translation_surfaces.cathedral(1, 2)
S.plot()
```

There are also infinite type surfaces (built from an infinite number of polygons) in that collection.

```{code-cell}
from flatsurf import translation_surfaces

S = translation_surfaces.infinite_staircase()
S.plot()
```

Some more general surfaces where the gluings are dilations are defined in the collection {class}`dilation_surfaces <flatsurf.geometry.similarity_surface_generators.DilationSurfaceGenerators>`.

```{code-cell}
from flatsurf import dilation_surfaces

S = dilation_surfaces.genus_two_square(1/2, 1/3, 1/4, 1/5)
S.plot()
```

Even more generality can be found in the collection {class}`similarity_surfaces <flatsurf.geometry.similarity_surface_generators.SimilaritySurfaceGenerators>`.

```{code-cell}
from flatsurf import Polygon, similarity_surfaces

P = Polygon(edges=[(2, 0), (-1, 3), (-1, -3)])
S = similarity_surfaces.self_glued_polygon(P)
S.plot()
```

### Building Surfaces from Polygons

Surfaces can also be built from scratch by specifying a (finite) set of polygons and gluings between the sides of the polygons.

Once you call `set_immutable()`, the type of the surface is determined, here a translation surface:

```{code-cell}
from flatsurf import MutableOrientedSimilaritySurface, Polygon

hexagon = Polygon(vertices=((0, 0), (3, 0), (3, 1), (3, 2), (0, 2), (0, 1)))

square = Polygon(vertices=((0, 0), (1, 0), (1, 1), (0, 1)))

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

```{code-cell}
T = MutableOrientedSimilaritySurface.from_surface(S)

T.glue((1, 1), (0, 2))
T.glue((0, 4), (0, 5))
T.set_immutable()

print(T)
T.plot()
```

Or anything that can be built by gluing with similarities…

```{code-cell}
T = MutableOrientedSimilaritySurface.from_surface(S)

T.glue((0, 0), (1, 1))
T.glue((0, 5), (0, 3))
T.set_immutable()

print(T)
T.plot()
```

There is a relatively small contract that a surface needs to implement, things such as "what is on the other side of this edge?", "is this surface compact?", … so it is not too hard to implement infinite type surfaces from scratch.

### Further Reading

* [Defining Surfaces](./defining_surfaces)
* {mod}`flatsurf.geometry.similarity_surface_generators`

+++

## Trajectories & Saddle Connections

Starting from a tangent vector, we can shoot a trajectory along that tangent vector. The trajectory will stop when it hits a singularity of the surface or after some (combinatorial) limit has been reached.

```{code-cell}
from flatsurf import translation_surfaces

S = translation_surfaces.infinite_staircase()
T = S.tangent_vector(0, (1/47, 1/49), v=(1, 1/31))

S.plot() + T.plot(color="red")
```

```{code-cell}
trajectory = T.straight_line_trajectory()
trajectory.flow(100)

S.plot() + trajectory.plot(color="red")
```

Note that on a finite type surface, we can compute a `FlowDecomposition` (see below) to determine how a surface decomposes in a direction into areas with periodic and dense trajectories.

+++

We can also determine all the saddle connections up to a certain length bound (on finite type surfaces).

```{code-cell}
from flatsurf import translation_surfaces
S = translation_surfaces.octagon_and_squares()

connections = S.saddle_connections(squared_length_bound=50)

# To get a more interesting picture, we color the saddle connections according to their length.
lengths = sorted(set(c.length() for c in connections))
color = lambda connection: colormaps.Accent(lengths.index(connection.length()) / len(lengths))[:3]

S.plot(polygon_labels=False, edge_labels=False) + sum(sc.plot(color=color(sc)) for sc in connections)
```

### Further Reading

* [Straight Line Flows](./straight_line_flow)
* [Saddle Connections](./saddle_connections)
* {mod}`flatsurf.geometry.straight_line_trajectory`

+++

## Supported Base Rings

Surfaces can be defined over most **exact** subrings of the reals from SageMath.

Here is an example defined over a number field:

```{code-cell}
from flatsurf import similarity_surfaces, Polygon

P = Polygon(angles=(1, 1, 4), vertices=[(0, 0), (1, 0)])
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")

print(S.base_ring())
S.plot(polygon_labels=False, edge_labels=False)
```

Here is a surface with a transcendental coordinate, supported through [exact-real](https://github.com/flatsurf/exact-real).

```{code-cell}
from pyexactreal import ExactReals

R = ExactReals(QuadraticField(3))
almost_one = R.random_element(1)

almost_one
```

```{code-cell}
from flatsurf import similarity_surfaces, Polygon

P = Polygon(angles=(1, 1, 4), vertices=[(0, 0), (almost_one, 0)])
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")

print(S.base_ring())
S.plot(polygon_labels=False, edge_labels=False)
```

## Flow Decompositions

We compute flow decompositions on the unfolding of the (2, 2, 5) triangle. We pick some direction coming from a saddle connection and decompose the surface into cylinders and minimal components in that direction.

```{code-cell}
from flatsurf import similarity_surfaces, Polygon
P = Polygon(angles=(2, 2, 5))
S = similarity_surfaces.billiard(P).minimal_cover(cover_type="translation")

G = S.graphical_surface(polygon_labels=False, edge_labels=False)
```

```{code-cell}
connection = next(iter(S.saddle_connections(squared_length_bound=4)))

G.plot() + connection.plot(color="red")
```

```{code-cell}
from flatsurf import GL2ROrbitClosure

O = GL2ROrbitClosure(S)
decomposition = O.decomposition(connection.holonomy())

decomposition
```

In this direction, the surface decomposes into 4 cylinders. Currently, sage-flatsurf cannot plot these cylinders. The optional ipyvue-flatsurf widgets can be used to visualize these cylinders in a Jupyter notebook.

```python
from ipyvue_flatsurf import Widget
Widget(decomposition)
```

+++

Here is a somewhat mysterious case. On the same surface, trying to decompose into a direction of a particular saddle connection, we can never certify that the components is minimal. (We ran this for a very long time and there really don't seem to be any cylinders out there.)

```{code-cell}
c = P.base_ring().gen()
direction = (12*c^4 - 60*c^2 + 56, 16*c^5 - 104*c^3 + 132*c)

O = GL2ROrbitClosure(S)
decomposition = O.decomposition(direction, limit=1024)

decomposition
```

We can also see the underlying Interval Exchange Transformation, to try to understand what is going on here.

```{code-cell}
O.decomposition(direction, limit=0).components()[0R].intervalExchangeTransformation()
```

### Further Reading

* [Siegel-Veech constants](./siegel_veech)
* the {mod}`module reference <flatsurf.geometry.gl2r_orbit_closure>` for computations related to the $GL_2(\mathbb{R})$-orbit closure

+++

## $SL_2(\mathbb{Z})$ orbits of square-tiled Surfaces

Square-tiled surfaces (also called origamis) are translation surfaces built from unit squares. Many properties of their $GL_2(\mathbb{R})$-orbit closures can be computed easily.

Note that the `Origami` object that we manipulate below are different from all other translation surfaces that were built above (that were using `sage-flatsurf`).

```{code-cell}
from surface_dynamics import Origami

O = Origami('(1,2)', '(1,3)')

O.plot()
```

```{code-cell}
T = O.teichmueller_curve()
```

```{code-cell}
T.sum_of_lyapunov_exponents()
```

```{code-cell}
O.lyapunov_exponents_approx()
```

The stratum of this surface and the stratum component.

```{code-cell}
O.stratum(), O.stratum_component()
```

### A Database of Arithmetic Teichmüller Curves

There is a relatively big database of pre-computed arithmetic Teichmüller curves that can be queried.

```{code-cell}
from surface_dynamics import OrigamiDatabase

D = OrigamiDatabase()

D.info()
```

The properties that are stored for each curve.

```{code-cell}
D.cols()
```

#### Some sample queries

The Teichmüller curves in genus $g$ such that in any direction we have at least $g$ cylinders.

```{code-cell}
for g in range(2, 7):
    q = D.query(('genus', '=', g), ('min_nb_of_cyls', '>=', g))
    print(f"g={g}: {q.number_of()}")
```

Some more information than just the count:

```{code-cell}
g = 3

q = D.query(('genus', '=', g), ('min_nb_of_cyls', '>=', g))

q.cols('stratum')

print(set(q))
```

The actual Origami representations from the "representative" column.

```{code-cell}
q.cols('representative')

o = choice(q.list())

print(o)
```

### Further Reading

* the [surface-dynamics documentation](https://flatsurf.github.io/surface-dynamics/)

+++

## Combinatorial Graphs up to Isomorphism

We list the graphs of genus 2 with 2 faces and vertex degree at least 3.

```{code-cell}
from surface_dynamics import FatGraphs

graphs = FatGraphs(g=2, nf=2, vertex_min_degree=3)

len(graphs.list())
```

One such graph in the list:

```{code-cell}
graphs.list()[0]
```

### Further Reading

* the [surface-dynamics documentation](https://flatsurf.github.io/surface-dynamics/)

+++

## Hyperbolic Geometry

sage-flatsurf provides an implementation of hyperbolic geometry that unlike the one in SageMath does not rely on the symbolic ring.

We define the hyperbolic plane over a number field containing a third root of 2.

```{code-cell}
from flatsurf import HyperbolicPlane

K.<a> = NumberField(x^3 - 2, embedding=1)

H = HyperbolicPlane(K)
```

We plot some geodesics in the upper half plane and in the Klein model.

```{code-cell}
g = H.geodesic(0, 1)
h = H.geodesic(1/2 - a/5, 2)

g.plot(color="red") + h.plot(color="blue")
```

```{code-cell}
g.plot(model="klein", color="red") + h.plot(model="klein", color="blue")
```

We determine the point of intersection of these geodesics. Note that that point has no coordinates in the upper half plane (without going to a quadratic extension).

```{code-cell}
P = g.intersection(h)
P.coordinates(model="klein")
```

```{code-cell}
P.change_ring(AA).coordinates(model="half_plane")
```

We use that point to define another geodesic to the infinite point ½.

```{code-cell}
g.plot(color="red") + h.plot(color="blue") + H.geodesic(P, 1/2).plot(color="orange")
```

### Further Reading

* the {mod}`module reference <flatsurf.geometry.hyperbolic>` for hyperbolic geometry

+++

## Veech Surfaces

We can explore how a Veech surfaces decomposes into cylinders.

```{code-cell}
from flatsurf import polygons, similarity_surfaces, GL2ROrbitClosure

T = polygons.triangle(1, 4, 7)
S = similarity_surfaces.billiard(T).minimal_cover('translation')
S = S.erase_marked_points()

S.plot()
```

```{code-cell}
O = GL2ROrbitClosure(S)
decomposition = O.decomposition((1, 0))

decomposition
```

```{code-cell}
bool(decomposition.parabolic())  # the underlying value is a tribool, since it could be undetermined
```

### Further Reading

* the {mod}`module reference <flatsurf.geometry.gl2r_orbit_closure>` related to the $GL_2(\mathbb{R})$-orbit closure

+++

## Strata & Orbit Closures

We can query the stratum a surface belongs to, and then (often) determine whether the surface has dense orbit closure.

```{code-cell}
from flatsurf import similarity_surfaces, polygons

T = polygons.triangle(2, 3, 8)
S = similarity_surfaces.billiard(T).minimal_cover('translation')

S.stratum()
```

```{code-cell}
from flatsurf import GL2ROrbitClosure

O = GL2ROrbitClosure(S)
O
```

```{code-cell}
for decomposition in O.decompositions(10, limit=20):
    if O.dimension() == O.ambient_stratum().dimension():
        break
    O.update_tangent_space_from_flow_decomposition(decomposition)

O
```

### Further Reading

* the {mod}`module reference <flatsurf.geometry.gl2r_orbit_closure>` related to the $GL_2(\mathbb{R})$-orbit closure
* [Exploring Orbit Closures](./apisa_wright)

+++

## Veering Triangulations

We create a Veering triangulation from the stratum $H_2(2)$.

```{code-cell}
from veerer import VeeringTriangulation, FlatVeeringTriangulation
from surface_dynamics import Stratum

H2 = Stratum([2], 1)
VT = VeeringTriangulation.from_stratum(H2)

VT
```

If you aim to study a specific surface, you might need to input the veering triangulation manually.

```{code-cell}
triangles = "(0,~7,6)(1,~5,~2)(2,4,~3)(3,8,~4)(5,7,~6)(~8,~1,~0)"
edge_slopes = "RBBRBRBRB"
VT2 = VeeringTriangulation(triangles, edge_slopes)

VT2
```

```{code-cell}
VT2.stratum()
```

We play with flat structures on a Veering triangulation. There are several pre-built constructions.

```{code-cell}
F0 = VT.flat_structure_middle()
F0.plot()
```

```{code-cell}
F1 = VT.flat_structure_geometric_middle()
F1.plot()
```

```{code-cell}
F2 = VT.flat_structure_min()
F2.plot()
```

We can work with $L^{\infty}$-Delaunay flat structures. These are also called *geometry* flat structures in `veerer`.

We check that the Veering triangulation `VT` corresponds to an open cell of the $L^\infty$-Delaunay decomposition of $H_2(2)$.

```{code-cell}
VT.is_geometric()
```

We compute the cone of $L^\infty$-Delaunay data for the given Veering triangulation `VT`.

Each point in the cone corresponds to a geometric flat structure given as $(x_0, \ldots, x_8, y_0, \ldots, y_8)$ where $(x_i, y_i)$ is the holonomy of the $i$-th edge.

The geometric structure is a polytope. Here the ambient dimension 18 corresponds to the fact that we have 9 edges (each edge has an $x$ and a $y$ coordinate). The dimension 8 is the real dimension of the stratum ($\dim_\mathbb{C} H_2(2) = 4$).

```{code-cell}
geometric_structures = VT.geometric_polytope()

geometric_structures
```

The rays allow to build any vector in the cone via linear combination (with non-negative coefficients).

```{code-cell}
rays = list(map(QQ**18, geometric_structures.rays()))
```

If all entries are positive, we have a valid geometric structure on `VT`.

```{code-cell}
xy = rays[0] + rays[12] + rays[5] + 10 * rays[6] + rays[10] + rays[16]

xy
```

Construct the associated flat structure (note that x and y are inverted).

```{code-cell}
flat_veering = VT._flat_structure_from_train_track_lengths(xy[9:], xy[:9])

flat_veering.plot()
```

We now explore the $L^\infty$-Delaunay of a given linear family.

We do it for the ambient stratum, i.e., the tangent space is everything.

The object one needs to start from is a pair of a Veering triangulation and a linear subspace.

```{code-cell}
L = VT.as_linear_family()
```

The *geometric automaton* is the set of such pairs that one obtains by moving around in the moduli space.

Initially, there is only the pair we provided.

```{code-cell}
A = L.geometric_automaton(run=False)

A
```

```{code-cell}
A.run(10)

A
```

We could compute everything (until there is nothing more to be explored).

If the computation terminates, it proves that `L` was indeed the tangent space to a $GL_2(\mathbb{R})$-orbit closure.

+++

### Further Reading

- [M. Bell, V. Delecroix, V. Gadre, R. Gutiérrez-Romo, Saul Schleimer  arXiv:1909.00890](https://arxiv.org/abs/1909.00890)
- [B. Zykoski arXiv:2206.04143](https://arxiv.org/abs/2206.04143)
