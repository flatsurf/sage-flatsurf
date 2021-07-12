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

+++ {"deletable": true, "editable": true}

# The GL(2,R) action, the Veech group, Delaunay decomposition

Initial version by Pat Hooper <whooper@ccny.cuny.edu>, Dec 16, 2017.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
from flatsurf import *
```

+++ {"deletable": true, "editable": true}

## Acting on surfaces by matrices.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
s = translation_surfaces.veech_double_n_gon(5)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
m=matrix([[2,1],[1,1]])
```

+++ {"deletable": true, "editable": true}

You can act on surfaces with the $GL(2,R)$ action

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss = m*s
ss
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss.plot()
```

+++ {"deletable": true, "editable": true}

To "renormalize" you can improve the presentation using the Delaunay decomposition.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sss = ss.delaunay_decomposition().copy(relabel=True)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sss.plot()
```

+++ {"deletable": true, "editable": true}

## The Veech group

Set $s$ to be the double pentagon again.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
s = translation_surfaces.veech_double_n_gon(5)
```

+++ {"deletable": true, "editable": true}

It is best to work in the field in which the surfact is defined.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
p=s.polygon(0)
p
```

+++ {"deletable": true, "editable": true}

The surface has a horizontal cylinder decomposition all of whose moduli are given as below

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
modulus = (p.vertex(3)[1]-p.vertex(2)[1])/(p.vertex(2)[0]-p.vertex(4)[0])
AA(modulus)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
m = matrix(s.base_ring(),[[1, 1/modulus],[0,1]])
show(m)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
show(matrix(AA,m))
```

+++ {"deletable": true, "editable": true}

The following can be used to check that $m$ is in the Veech group of $s$.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.canonicalize() == (m*s).canonicalize()
```

+++ {"deletable": true, "editable": true}

## Infinite surfaces

Infinite surfaces support multiplication by matrices and computing the Delaunay decomposition. (Computation is done "lazily")

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s=translation_surfaces.chamanara(1/2)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot(edge_labels=False,polygon_labels=False)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
ss=s.delaunay_decomposition()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
ss.graphical_surface().make_all_visible(limit=20)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss.plot(edge_labels=False,polygon_labels=False)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
m = matrix([[2,0],[0,1/2]])
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
ms = m*s
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ms.graphical_surface().make_all_visible(limit=20)
ms.plot(edge_labels=False,polygon_labels=False)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
mss = ms.delaunay_decomposition()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
mss.graphical_surface().make_all_visible(limit=20)
mss.plot(edge_labels=False,polygon_labels=False)
```

+++ {"deletable": true, "editable": true}

You can tell from the above picture that $m$ is in the Veech group.
