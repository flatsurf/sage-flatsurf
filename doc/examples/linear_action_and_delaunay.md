---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# The GL(2,R) Action, the Veech Group, Delaunay Decomposition

## Acting on surfaces by matrices.

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
from flatsurf import translation_surfaces

s = translation_surfaces.veech_double_n_gon(5)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s.plot()
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
m = matrix([[2, 1], [1, 1]])
```

You can act on surfaces with the $GL(2,R)$ action

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
ss = m * s
ss
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
ss.plot()
```

To "renormalize" you can improve the presentation using the Delaunay decomposition.

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
sss = ss.delaunay_decomposition()
sss
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
sss.plot()
```

## The Veech group

Set $s$ to be the double pentagon again.

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
s = translation_surfaces.veech_double_n_gon(5)
```

The surface has a horizontal cylinder decomposition all of whose moduli are given as below

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
p = s.polygon(0)
modulus = (p.vertex(3)[1] - p.vertex(2)[1]) / (p.vertex(2)[0] - p.vertex(4)[0])
AA(modulus)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
m = matrix(s.base_ring(), [[1, 1 / modulus], [0, 1]])
show(m)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
show(matrix(AA, m))
```

The following can be used to check that $m$ is in the Veech group of $s$.

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s.canonicalize() == (m * s).canonicalize()
```

## Infinite surfaces

Infinite surfaces support multiplication by matrices and computing the Delaunay decomposition. (Computation is done "lazily")

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s = translation_surfaces.chamanara(1 / 2)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s.plot(edge_labels=False, polygon_labels=False)
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
ss = s.delaunay_decomposition()
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
gs = ss.graphical_surface(edge_labels=False, polygon_labels=False)
gs.make_all_visible(limit=20)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
gs.plot()
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
m = matrix([[2, 0], [0, 1 / 2]])
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
ms = m * s
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
gs = ms.graphical_surface(edge_labels=False, polygon_labels=False)
gs.make_all_visible(limit=20)
gs.plot()
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
mss = ms.delaunay_decomposition()
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
gs = mss.graphical_surface(edge_labels=False, polygon_labels=False)
gs.make_all_visible(limit=20)
gs.plot()
```

You can tell from the above picture that $m$ is in the Veech group.
