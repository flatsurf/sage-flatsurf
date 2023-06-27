---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.6
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# Playing with the monotile

The monotile is a polygon that tiles the plane by rotation and reflections. Interestingly
one can build two translation surfaces out of it.

The monotile is the following polygon.
```{code-cell}
from flatsurf import Polygon, MutableOrientedSimilaritySurface

K = QuadraticField(3)
a = K.gen()

# build the vectors
l = [(1, 0), (1, 2), (a, 11), (a, 1), (1, 4), (1, 6), (a, 3),
     (a, 5), (1, 8), (1, 6), (a, 9), (a, 7), (1, 10), (1, 0)]
vecs = []
for m, e in l:
    v = vector(K, [m * cos(2*pi*e/12), m * sin(2*pi*e/12)])
    vecs.append(v)
p = Polygon(edges=vecs)
```

One can build translation surfaces by gluing parallel edges. There is an
ambiguity in doing so because of the horizontal segments. We a surface
`Sbase` where non-ambiguous gluings are performed.
```{code-cell}
from collections import defaultdict
d = defaultdict(list)
for i, e in enumerate(p.edges()):
    e.set_immutable()
    d[e].append(i)
```

```{code-cell}
Sbase = MutableOrientedSimilaritySurface(K)
_ = Sbase.add_polygon(p)
for v in list(d):
    if v in d:
        indices = d[v]
        v_op = -v
        v_op.set_immutable()
        opposite_indices = d[v_op]
        assert len(indices) == len(opposite_indices), (len(indices), len(opposite_indices))
        if len(indices) == 1:
            del d[v]
            del d[v_op]
            Sbase.glue((0, indices[0]), (0, opposite_indices[0]))
```

Next we recover the ambiguous edges and build the two possible remaining gluings.
```{code-cell}
assert len(d) == 2
(i0, j0), (i1, j1) = d.values()
```

```{code-cell}
S1 = MutableOrientedSimilaritySurface.from_surface(Sbase)
S1.glue((0, i0), (0, i1))
S1.glue((0, j0), (0, j1))
S1.set_immutable()
```

```{code-cell}
S2 = MutableOrientedSimilaritySurface.from_surface(Sbase)
S2.glue((0, i0), (0, j1))
S2.glue((0, j0), (0, i1))
S2.set_immutable()
```

We indeed obtain translation surfaces
```{code-cell}
print(S1.category())
print(S2.category())
```

And one can compute their genera
```{code-cell}
print(S1.genus(), S2.genus())
```

and strata
```{code-cell}
print(S1.stratum(), S2.stratum())
```

```{code-cell}
S1.triangulate()
```{code-cell}
S2.triangulate()
```
