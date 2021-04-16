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

# Graphics Configuration

## Rearranging Polygons

```{code-cell} ipython3
from flatsurf import *
from flatsurf.geometry.polyhedra import *
```

```{code-cell} ipython3
s = platonic_dodecahedron()[1]
```

The default way plotting looks:

```{code-cell} ipython3
s.plot()
```

Labels in the center of each polygon indicate the label of the polygon. Edge labels above indicate which polygon the edge is glued to.

Plotting the surface is controlled by a GraphicalSurface object. You can get the surface as follows:

```{code-cell} ipython3
gs = s.graphical_surface()
```

The graphical surface controls where polygons are drawn. You can glue a polygon across an edge using `gs.make_adjacent(label,edge)`. A difficulty is that you need to know which edge is which. You can enable `zero_flags` to see the zero vertex of each polygon.

```{code-cell} ipython3
gs.will_plot_zero_flags = True
```

```{code-cell} ipython3
gs.plot()
```

FlatSurf just uses some simple algorithm to layout polygons. Sometimes they overlap for example. I'm troubled by that the picture is asymmetric: Polygon 4 should have a pentagon 10 sticking to it. We can see by counting counterclockwise from the red flag in polygon 4 that the 10 appears on edge 3. We can verify that with:

```{code-cell} ipython3
s.opposite_edge(4,3)
```

We can move polygon 10 so that it is adjacent to polygon 4 with the command:

```{code-cell} ipython3
gs.make_adjacent(4,3)
```

Lets check that it worked:

```{code-cell} ipython3
s.plot()
```

Oops. Lets also move polygon 11.

```{code-cell} ipython3
gs.make_adjacent(7,4)
```

```{code-cell} ipython3
gs.plot()
```

There are other options which can change what is displayed on an edge. The option can also be passed to `s.graphical_surface()`.

```{code-cell} ipython3
gs.process_options(edge_labels="gluings and number")
```

```{code-cell} ipython3
gs.plot()
```

```{code-cell} ipython3
gs.process_options(edge_labels="number")
```

```{code-cell} ipython3
gs.plot()
```

## Moving between coordinate systems

The Euclidean Cone Surface `s` works in a different coordinate system then the graphical surface `gs`. So, when we moved the polygon above, we had no affect on `s`. In fact, the polygons of `s` are all the same:

```{code-cell} ipython3
s.polygon(0)
```

```{code-cell} ipython3
s.polygon(1)
```

So really `s` is a disjoint union of twelve copies of a standard pentagon with some edge gluings.

Lets now look at "graphical coordinates" i.e., the coordinates in which `gs` works.

```{code-cell} ipython3
show(gs.plot(), axes=True)
```

I can tell that the point `(3,5)` is in the unfolding, but I can't tell if it is in polygon 3 or 7. The GraphicalSurface `gs` is made out of GraphicalPolygons which we can use to deal with this sort of thing.

```{code-cell} ipython3
gs.graphical_polygon(3).contains((3,5))
```

```{code-cell} ipython3
gs.graphical_polygon(7).contains((3,5))
```

Great. Now we can get the position of the point on the surface!

```{code-cell} ipython3
gp = gs.graphical_polygon(7)
pt = gp.transform_back((3,5))
pt
```

Here we plot polygon 6 in its geometric coordinates with `pt`.

```{code-cell} ipython3
s.polygon(6).plot()+point2d([pt],zorder=100)
```

Lets convert it to a surface point and plot it!

```{code-cell} ipython3
spt = s.surface_point(7,pt)
spt
```

```{code-cell} ipython3
s.plot() + spt.plot(color="red", size=20)
```

Now I want to plot an upward trajectory through this point. Again, I have to deal with the fact that the coordinates might not match. You can get access to the transformation (a similarity) from geometric coordinates
to graphical coordinates:

```{code-cell} ipython3
transformation = gs.graphical_polygon(6).transformation()
transformation
```

Really we want the inverse:

```{code-cell} ipython3
inverse_transformation = ~transformation
inverse_transformation
```

We just want the derivative of this similarity to transform the vertical direction. The derivative is a $2 \times 2$ matrix.

```{code-cell} ipython3
show(inverse_transformation.derivative())
```

```{code-cell} ipython3
direction = inverse_transformation.derivative()*vector((0,1))
direction
```

We can use the point and the direction to get a tangent vector, which we convert to a trajectory, flow and plot.

```{code-cell} ipython3
tangent_vector = s.tangent_vector(7, pt, direction)
tangent_vector
```

```{code-cell} ipython3
traj = tangent_vector.straight_line_trajectory()
traj.flow(100)
traj.is_closed()
```

```{code-cell} ipython3
s.plot()+spt.plot(color="red")+traj.plot()
```

## Multiple graphical surfaces

It is possible to have more than one graphical surface. Maybe you want to have one where things look better.
To get a new suface, you can call `s.graphical_surface()` again but with a `cached=False` parameter.

```{code-cell} ipython3
pretty_gs = s.graphical_surface(polygon_labels=False, edge_labels=False, cached=False)
```

```{code-cell} ipython3
pretty_gs.plot()
```

Current polygon printing options:

```{code-cell} ipython3
pretty_gs.polygon_options
```

```{code-cell} ipython3
del pretty_gs.polygon_options["color"]
pretty_gs.polygon_options["rgbcolor"]="#ffeeee"
```

```{code-cell} ipython3
pretty_gs.non_adjacent_edge_options["thickness"] = 0.5
pretty_gs.non_adjacent_edge_options["color"] = "lightblue"
pretty_gs.will_plot_adjacent_edges = False
```

```{code-cell} ipython3
pretty_gs.plot()
```

To use a non-default graphical surface you need to pass the graphical surface as a parameter.

```{code-cell} ipython3
pretty_gs.plot()+spt.plot(pretty_gs, color="red")+traj.plot(pretty_gs)
```

Lets make it prettier by drawing some stars on the faces!

Find all saddle connections of length at most $\sqrt{16}$:

```{code-cell} ipython3
saddle_connections = s.saddle_connections(16)
```

The edges have length two so we will keep anything that has a different length.

```{code-cell} ipython3
saddle_connections2 = []
for sc in saddle_connections:
    h = sc.holonomy()
    if h[0]**2 + h[1]**2 != 4:
        saddle_connections2.append(sc)
len(saddle_connections2)
```

```{code-cell} ipython3
plot = pretty_gs.plot()
for sc in saddle_connections2:
    plot += sc.plot(pretty_gs, color="red")
plot
```

Plot using the original graphical surface.

```{code-cell} ipython3
plot = s.plot()
for sc in saddle_connections2:
    plot += sc.plot(color="red")
plot
```

## Manipulating edge labels

```{code-cell} ipython3
s = translation_surfaces.arnoux_yoccoz(4)
```

```{code-cell} ipython3
s.plot()
```

Here is an example with the edge labels centered on the edge.

```{code-cell} ipython3
gs = s.graphical_surface(cached=False)
```

```{code-cell} ipython3
del gs.polygon_options["color"]
gs.polygon_options["rgbcolor"]="#eee"
gs.edge_label_options["position"]="edge"
gs.edge_label_options["t"]=0.5
gs.edge_label_options["push_off"]=0
gs.edge_label_options["color"]="green"
gs.adjacent_edge_options["thickness"]=0.5
gs.will_plot_non_adjacent_edges=False
```

```{code-cell} ipython3
gs.plot()
```

```{code-cell} ipython3
gs = s.graphical_surface(cached=False)
```

```{code-cell} ipython3
del gs.polygon_options["color"]
gs.polygon_options["rgbcolor"]="#eef"
gs.edge_label_options["position"]="outside"
gs.edge_label_options["t"]=0.5
gs.edge_label_options["push_off"]=0.02
gs.edge_label_options["color"]="green"
gs.adjacent_edge_options["thickness"]=0.5
gs.non_adjacent_edge_options["thickness"]=0.25
```

```{code-cell} ipython3
gs.plot()
```

## Hacking GraphicalSurface

```{code-cell} ipython3
s = translation_surfaces.infinite_staircase()
```

```{code-cell} ipython3
gs = s.graphical_surface(polygon_labels=False, edge_labels=False)
```

```{code-cell} ipython3
gs.plot()
```

The methods `plot_polygon`, `plot_polygon_label`, `plot_edge`, `plot_edge_label`, `plot_zero_flag` are fairly
simple. They just call methods in a GraphicalPolygon which was passed as a parameter. We can replace any of these functions to customize their behaviour.

```{code-cell} ipython3
# Define a replacement method.
def plot_polygon(self, label, graphical_polygon, upside_down):
    if label%2==0:
        return graphical_polygon.plot_polygon(color="lightgreen")
    else:
        return graphical_polygon.plot_polygon(color="yellow")        
```

```{code-cell} ipython3
# Replace the method in gs.
from types import MethodType
gs.plot_polygon = MethodType(plot_polygon, gs)
```

```{code-cell} ipython3
gs.plot()
```
