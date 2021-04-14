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

# Saddle Connection Demo

Initial version by Pat Hooper <whooper@ccny.cuny.edu> on Dec. 16, 2017.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
from flatsurf import *
```

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

We get a list of all saddle connections of length less than $\sqrt{100}$.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sc_list = s.saddle_connections(100)
len(sc_list)
```

+++ {"deletable": true, "editable": true}

The following removes dublicate saddle connections which appear with opposite orientations.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sc_set = set()
for sc in sc_list:
    if sc.invert() not in sc_set:
        sc_set.add(sc)
sc_list2 = [sc for sc in sc_set]
len(sc_list2)
```

+++ {"deletable": true, "editable": true}

We pick two saddle connections:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sc1 = sc_list2[99]
sc2 = sc_list2[55]
```

+++ {"deletable": true, "editable": true}

We can find their holonomies and other information about them using methods.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
print("Holonomy of sc1 is"+str(sc1.holonomy())+" = "+str(sc1.holonomy().n()))
print("Holonomy of sc2 is"+str(sc2.holonomy())+" = "+str(sc2.holonomy().n()))
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot() + sc1.plot() + sc2.plot(color="green")
```

+++ {"deletable": true, "editable": true}

We can test that they intersect. By default the singularity does not count.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sc1.intersects(sc2)
```

+++ {"deletable": true, "editable": true}

We can get an iterator over the set of intersection points:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
for p in sc1.intersections(sc2):
    print p
```

+++ {"deletable": true, "editable": true}

It is a good idea to store the intersections in a list if you want to reuse them:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
intersections = [p for p in sc1.intersections(sc2)]
```

+++ {"deletable": true, "editable": true}

We can plot the intersection points:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
plot = s.plot() + sc1.plot() + sc2.plot(color="green")
for p in intersections:
    plot += p.plot(color="red", zorder=3)
plot
```

+++ {"deletable": true, "editable": true}

We can plot all the saddle connections:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
plot=s.plot(edge_labels=False, polygon_labels=False)
for sc in sc_list:
    plot += sc.plot(thickness=0.05)
plot
```

We will build a subset of the saddle connection graph where vertices are saddle connections and two vertices are joined by an edge if and only if the saddle connections do not intersect (on their interiors).

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
# Build intersection graph
d={}
for i in xrange(len(sc_list2)):
    for j in xrange(i+1,len(sc_list2)):
        if not sc_list2[i].intersects(sc_list2[j]):
            if i not in d:
                d[i]=[j]
            else:
                d[i].append(j)
            if j not in d:
                d[j]=[i]
            else:
                d[j].append(i)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
g=Graph(d)
```

We place the vertex of a saddle connection with holonomy $z \in {\mathbb C}$ at the point $z^2$.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
pos={}
for i in xrange(len(sc_list2)):
    sc = sc_list2[i]
    val = sc.holonomy().n()
    z = val[0]+I*val[1]
    w = z**2/z.abs()
    pos[i]=(w.real(),w.imag())
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
g.plot(pos=pos,vertex_labels=False,vertex_size=0)
```
