---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: SageMath 9.7
    language: sage
    name: sagemath
---

# Harmonic Differentials


First, some differentials on a square torus.

(TODO: Unfortunately, we need to explicitly Delaunay triangulate the torus for this to work. We should remove this limitation and hide the triangulation as an implementation detail…)

```sage
from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
T = translation_surfaces.torus((1, 0), (0, 1)).delaunay_triangulation()
T.set_immutable()
T.plot()
```

We create differentials with prescribed values on generators of homology. The differential is, modulo some numerical noise, a constant:

```sage
H = SimplicialHomology(T)
a, b = H.gens()
```

```sage
H = SimplicialCohomology(T)
f = H({a: 1})
```

```sage
Omega = HarmonicDifferentials(T)
omega = Omega(f, prec=2)
omega
```

The power series is developed around the centers of the circumcircle of the triangulation. In this example the centers for the triangles are the same, so the coefficients are (essentially) forced to be identical.

We can recover the series for each triangle of the triangulation:

```sage
omega.series(0)
```

```sage
omega.series(1)
```

We can ask the differential how well it solves the constraints that were used to created it:

```sage
omega.error(verbose=True)
```

If we create the differential with much more precision, we see some numerical noise here:

```sage
Omega(f, prec=40).error(verbose=True)
```

We can also use other strategies to determine the differential. The supported strategies are `L2` (the default), `midpoint_derivatives` (forcing derivatives up to some point to match at the midpoints of the triangulation edges,) `area_upper_bound` (minimize an approximation of the area,) `area` (minimize the area).

```sage
Omega(f, prec=2, algorithm=["midpoint_derivatives"])
```

These strategies can also be mixed and weighted differently (in the case of `midpoint_derivatives`, this controls up to which derivative we force derivatives to match).

There are checks for obvious errors in the computation, e.g., when the error in the L2 norm gets too big:

```sage
Omega(f, prec=10, algorithm={"midpoint_derivatives": 2, "area_upper_bound": 10, "L2": 0})
```

These checks can be disabled though:

```sage
Omega(f, prec=10, algorithm={"midpoint_derivatives": 2, "area_upper_bound": 10, "L2": 0}, check=False)
```

There are some other basic operations supported. We can, e.g., ask for the roots of a differential (TODO: This does not include roots at the vertices of the triangulation yet):

```sage
omega.roots()
```

At the vertices we can ask for the coefficients of the power series developed around that vertex:

```sage
vertex = T.angles(return_adjacent_edges=True)[0][1]
```

```sage
omega.cauchy_residue(vertex, 0)
```

```sage
omega.cauchy_residue(vertex, -1)
```

## A Less Trivial Example, the Regular Octagon

```sage
from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
T = translation_surfaces.regular_octagon().delaunay_triangulation()
T.set_immutable()
T.plot()
```

```sage
H = SimplicialHomology(T)
a, b, c, d = H.gens()
```

We create a differential whose integral along `a` is 1 and 0 on the other generators of homology. Note that `a` can be written as the edge 0 on the polygon 0, i.e., the right edge of that polygon:

```sage
a
```

```sage
H = SimplicialCohomology(T)
f = H({a: 1})
```

**TODO**: Unfortunately this fails. We don't find solution here.

```sage
Omega = HarmonicDifferentials(T)
omega = Omega(f, prec=20)
```

### An Explicit Series for the Octagon
We can also provide the series for the Voronoi cells explicitly if we don't want to solve for a cohomology class.

```sage
from flatsurf import translation_surfaces, HarmonicDifferentials, SimplicialHomology, SimplicialCohomology
T = translation_surfaces.regular_octagon().delaunay_triangulation()
T = T.apply_matrix(diagonal_matrix([2.32718 / 2, 2.32718 / 2]))
T.set_immutable()

Omega = HarmonicDifferentials(T)
```

```sage
R.<z> = CC[[]]
f = z^2 - 1/9*z^10 + 20/1377*z^18 - 14/6885*z^26 + 2044/6952473*z^34 -111097/2565462537*z^42 + O(z^43)
```

```sage
omega = Omega({triangle: f for triangle in T.label_iterator()})
```

```sage
omega.error(verbose=True)
```
