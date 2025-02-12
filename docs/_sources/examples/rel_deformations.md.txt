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

# Relative Period Deformations

## The Arnoux-Yoccoz surface

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
from flatsurf import translation_surfaces

s = translation_surfaces.arnoux_yoccoz(3).canonicalize()
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
field = s.base_ring()
field
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
alpha = field.gen()
AA(alpha)
```

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
m = matrix(field, [[alpha, 0], [0, 1 / alpha]])
show(m)
```

Check that $m$ is the derivative of a pseudo-Anosov of $s$.

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
(m * s).canonicalize() == s
```

## Rel deformation

A singularity of the surface is an equivalence class of vertices of the polygons making up the surface.

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s.point(0, 0)
```

We'll move this singularity to the right by two different amounts:

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s1 = s.rel_deformation(
    {s.point(0, 0): vector(field, (alpha / (1 - alpha), 0))}
).canonicalize()
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
---
s2 = s.rel_deformation(
    {s.point(0, 0): vector(field, (1 / (1 - alpha), 0))}
).canonicalize()
```

+++ {"jupyter": {"outputs_hidden": true}}

Note that by the action of the derivative of the pseudo-Anosov we have:

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s1 == (m * s2).canonicalize()
```

By a Theorem of Barak Weiss and the author of this notebook, these surfaces are all periodic in the vertical direction. You can see the vertical cylinders:

```{code-cell}
---
jupyter:
  outputs_hidden: false
---
s1.plot()
```
