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

# Straight-line flow

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

+++ {"deletable": true, "editable": true}

Defines the tangent_bundle on the surface defined over the base_ring of s.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
TB = s.tangent_bundle()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
baricenter = sum(s.polygon(0).vertices())/5
```

+++ {"deletable": true, "editable": true}

Define the tangent vector based at the baricenter of polygon 0 aimed downward.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
v = TB(0, baricenter, (0,-1))
```

+++ {"deletable": true, "editable": true}

Convert to a straight-line trajectory. Trajectories are unions of segments in polygons.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
traj = v.straight_line_trajectory()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot()+traj.plot()
```

+++ {"deletable": true, "editable": true}

Flow into the next $100$ polygons or until the trajectory hits a vertex.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
traj.flow(100)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot()+traj.plot()
```

+++ {"deletable": true, "editable": true}

We can tell its type.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
traj.is_saddle_connection()
```

+++ {"deletable": true, "editable": true}

You can also test if a straight-line trajectory is closed or a forward/backward separatrix.

Lets do it again but in the slope one direction.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
v = TB(0, baricenter, (1,1))
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
traj=v.straight_line_trajectory()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
traj.flow(100)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot()+traj.plot()
```

+++ {"deletable": true, "editable": true}

We remark that it follows from work of Veech that the slope one direction is ergodic for the straight-line flow.
