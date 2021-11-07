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

# Relative Period Deformations

Initial version by Pat Hooper <whooper@ccny.cuny.edu>, Dec. 16, 2017.

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

## The Arnoux-Yoccoz surface

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
s = translation_surfaces.arnoux_yoccoz(3).canonicalize()
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
field=s.base_ring()
field
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
alpha = field.gen()
AA(alpha)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
m=matrix(field,[[alpha,0],[0,1/alpha]])
show(m)
```

+++ {"deletable": true, "editable": true}

Check that $m$ is the derivative of a pseudo-Anosov of $s$.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
(m*s).canonicalize()==s
```

+++ {"deletable": true, "editable": true}

## Rel deformation

A singularity of the surface is an equivalence class of vertices of the polygons making up the surface.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.singularity(0,0)
```

+++ {"deletable": true, "editable": true}

We'll move this singularity to the right by two different amounts:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s1=s.rel_deformation({s.singularity(0,0):vector(field,(alpha/(1-alpha),0))}).canonicalize()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
s2=s.rel_deformation({s.singularity(0,0):vector(field,(1/(1-alpha),0))}).canonicalize()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
# Note that by the action of the derivative of the pseudo-Anosov we have:
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s1==m*s2
```

+++ {"deletable": true, "editable": true}

By a Theorem of Barak Weiss and the author of this notebook, these surfaces are all periodic in the vertical direction. You can see the vertical cylinders:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s1.plot()
```
