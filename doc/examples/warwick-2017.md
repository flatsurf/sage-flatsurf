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

# FlatSurf Demo

## Warwick EPSRC Symposium on "Computation in geometric topology"

By Pat Hooper <whooper@ccny.cuny.edu>

This material is based upon work supported by the National Science Foundation
under Grant 1500965.

Flatsurf is a Sage module to work with "flat surfaces"
written so far mostly by Vincent Delecroix and Pat Hooper.
It is developing in random directions depending on our interests,
but we would like it to be as helpful as possible for everyone...
We welcome other contributors!

To get flatsurf you first need <a href="https://www.sagemath.org/">Sage</a>.
Then you can get flatsurf by following the directions here on the
<a href="https://github.com/flatsurf/sage-flatsurf#installing-the-module">flatsurf GitHub page</a>.
(Currently we recommend running

<code>sage -pip install git+https://github.com/flatsurf/sage-flatsurf</code>

to install.)

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

## Veech group elements (affine symmetries)

Veech's double n-gon surfaces:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s = translation_surfaces.veech_double_n_gon(5).canonicalize()
s.plot()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
p = s.polygon(0)
modulus = (p.vertex(3)[1] - p.vertex(2)[1]) / (p.vertex(2)[0] - p.vertex(4)[0])
AA(modulus)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
m = matrix(s.base_ring(), [[1, 2], [0, 1]])
show(matrix(AA, m))
ss = m*s
ss.plot()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss.delaunay_decomposition().plot()
```

+++ {"deletable": true, "editable": true}

The following checks that the matrix m stabilizes s:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss.canonicalize() == s
```

+++ {"deletable": true, "editable": true}

## Geodesics

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s = translation_surfaces.veech_double_n_gon(5)
```

+++ {"deletable": true, "editable": true}

The tangent bundle of the surface:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
TB = s.tangent_bundle()
```

+++ {"deletable": true, "editable": true}

Define a tangent vector in polygon $0$ starting at $(\frac{1}{2}, 0)$ and pointed in some direction:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
direction = s.polygon(0).vertex(2) + 3*s.polygon(0).vertex(3)
v = TB(0, (1/2, 0), direction)
```

+++ {"deletable": true, "editable": true}

Convert the vector to a straight-line trajectory.

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
s.plot() + traj.plot()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
traj.flow(1000)
print(traj.is_closed())
print(traj.combinatorial_length())
s.plot() + traj.plot()
```

+++ {"deletable": true, "editable": true}

## Cone surfaces from polyhedra

Polyhedra are built into Sage and you can use them to build a translation surface.
In this demo we only use a built-in function for a Platonic Solid.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
from flatsurf.geometry.polyhedra import *
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
polyhedron, s, mapping = platonic_dodecahedron()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
polyhedron.plot(frame=False)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot(polygon_labels=False, edge_labels=False)
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
TB = s.tangent_bundle()
direction = s.polygon(0).vertex(2) + 2*s.polygon(0).vertex(3)
v = TB(0, (1/2, 0), direction)
traj = v.straight_line_trajectory()
traj.flow(100)
print(traj.is_closed())
print(traj.combinatorial_length())
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s.plot() + traj.plot()
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
G = polyhedron.plot(frame=False, point=False, line=False, wireframe=None)
G += line3d(mapping(traj), radius=0.02, frame=False)
G
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
TB = s.tangent_bundle()
direction = s.polygon(0).vertex(2) + 3*s.polygon(0).vertex(3)
v = TB(0, (1/2, 0), direction)
traj = v.straight_line_trajectory()
traj.flow(1000)
print(traj.is_closed())
print(traj.combinatorial_length())
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
show(s.plot() + traj.plot())
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
p = polyhedron.plot(frame=False, point=False, line=False, wireframe=None)
p += line3d(mapping(traj), radius=0.02, frame=False)
p.show(viewer='tachyon', frame=False)
```

+++ {"deletable": true, "editable": true}

## Relative period deformations

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s = translation_surfaces.veech_2n_gon(5)
s.plot(edge_labels=False, polygon_labels=False)
```

+++ {"deletable": true, "editable": true}

Currently we have to triangulate to do a rel deformation.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s = s.triangulate().copy(relabel=True, mutable=True)
```

+++ {"deletable": true, "editable": true}

A singularity is an equivalence class of vertices of polygons.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
sing = s.singularity(0, 0)
sing
```

+++ {"deletable": true, "editable": true}

We can now deform by moving one singularity relative to the others.
Here is a small deformation in the slope one direction.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss = s.rel_deformation({sing: vector(s.base_ring(), (1/20, 1/20))})
ss.plot()
```

+++ {"deletable": true, "editable": true}

A larger deformation:

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
ss = s.rel_deformation({sing:vector(s.base_ring(), (100, 100))})
ss.plot()
```

+++ {"deletable": true, "editable": true}

## The Necker Cube Surface

I'm demonstrating a result (in progress) of Pavel Javornik,
an undergraduate at City College of New York.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
from flatsurf.geometry.straight_line_trajectory import StraightLineTrajectory

class SurfaceToSpaceMapping(SageObject):

    def __init__(self, similarity_surface, tranformation):
        self._s = similarity_surface
        from types import FunctionType
        if isinstance(transformation, FunctionType):
            self.transformation = transformation

    def transformation(self, label):
        r"""
        Return a pair (m, t) where m is a 3x2 matrix and t is a vector with 3 entries.

        The associated tranformation from the polygon with the given label
        is v mapsto m*v + t where v is a point in the polygon.
        """
        return self._t[label]
    
    def image_polygon(self, label):
        r"""
        Return a 2-dimensional polyhedron in 3-space representing
        the image of the polygon with the given label.
        """
        p = self._s.polygon(label)
        m, t = self.transformation(label)
        vertices = [m*v + t for v in p.vertices()]
        return Polyhedron(vertices=vertices)

    def plot(self, labels, point=False, line=False, polygon=None,
             wireframe=None, frame=False, label_to_color=None):
        r"""
        Return a 3d plot of the polygonal images in 3-space
        corresponding to the collection of labels.

        The other parameters are passed to a Polyhedron.plot method
        and affect the rendering.
        """
        it = iter(labels)
        label = next(it)
        if label_to_color is None:
            p = self.image_polygon(label).plot(point=point,
                                               line=line,
                                               polygon=polygon,
                                               wireframe=wireframe,
                                               frame=frame,
                                               color="pink")
        else:
            p = self.image_polygon(label).plot(point=point,
                                               line=line,
                                               polygon=polygon,
                                               wireframe=wireframe,
                                               frame=frame,
                                               color=label_to_color(label))
        for label in it:
            if label_to_color is None:
                p += self.image_polygon(label).plot(point=point,
                                                    line=line,
                                                    polygon=polygon,
                                                    wireframe=wireframe,
                                                    frame=frame,
                                                    color="pink")
            else:
                p += self.image_polygon(label).plot(point=point,
                                                    line=line,
                                                    polygon=polygon,
                                                    wireframe=wireframe,
                                                    frame=frame,
                                                    color=label_to_color(label))
        from sage.modules.free_module_element import vector
        p.frame_aspect_ratio(tuple(vector(p.bounding_box()[1])
                                   - vector(p.bounding_box()[0])))
        return p

    def __call__(self, o):
        r"""
        This method is used to convert from an object on the surface to an object in space.

        Currently works with

        - ``StraightLineTrajectory`` -- returns the corresponding list of points in space
        - ``SegmentInPolygon`` -- returns the corresponding pair of points in space
        - ``SimilaritySurfaceTangentVector`` -- returns a pair of points corresponding
          to the image point and image of the tangent vector.
        """
        if isinstance(o, StraightLineTrajectory):
            points = []
            it = iter(o.segments())
            s = next(it)
            label = s.polygon_label()
            m, t = self.transformation(label)
            points.append(t + m*s.start().point())
            points.append(t + m*s.end().point())
            for s in it:
                label = s.polygon_label()
                m, t = self.transformation(label)
                points.append(t + m*s.end().point())
            return points
        if isinstance(o, SegmentInPolygon):
            # Return the pair of images of the endpoints.
            label = o.polygon_label()
            m, t = self.transformation(label)
            return (t + m*o.start().point(), t + m*o.end().point())
        if isinstance(o, SimilaritySurfaceTangentVector):
            # Map to a pair of vectors consisting of the image
            # of the basepoint and the image of the vector.
            label = o.polygon_label()
            m, t = self.transformation(label)
            point = o.point()
            vector = o.vector()
            return (t + m*point, m*vector)
        raise ValueError("Failed to recognize type of passed object")
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
from flatsurf.geometry.surface import Surface
from flatsurf.geometry.polygon import ConvexPolygons
from flatsurf.geometry.similarity import SimilarityGroup
class CubeSurf(Surface):
    def __init__(self, F):
        ZZ3 = IntegerModRing(3)
        P = ConvexPolygons(F)
        self._faceA = P(vertices = [(0, 0), (1, 0), (1, 1), (0, 1)])
        self._faceB = P(vertices = [(0, 0), (1, 0), (1, 1), (0, 1)])
        self._faceC = P(vertices = [(0, 0), (1, 0), (1, 1), (0, 1)])
        Surface.__init__(self, F, (ZZ(0), ZZ(0), ZZ3(0)), finite=False, mutable=False)

    def polygon(self, label):
        x, y, l = label
        if l == 0:
            return self._faceA
        if l == 1:
            return self._faceB
        if l == 2:
            return self._faceC

    def opposite_edge(self, label, edge):
        x, y, l = label
        # l(0) = A, l(1) = B, l(2) = C
        if l == 0:
            if edge == 0:
                return((x, y - 1, l + 2), 2)
            if edge == 1:
                return((x, y, l + 1), 3)
            if edge == 2:
                return((x, y, l + 2), 0)
            if edge == 3:
                return((x - 1, y, l + 1), 1)
        if l == 1:
            if edge == 0:
                return((x + 1, y - 1, l + 1), 3)
            if edge == 1:
                return((x + 1, y, l + 2), 3)
            if edge == 2:
                return((x, y, l + 1), 1)
            if edge == 3:
                return((x, y, l + 2), 1)
        if l == 2:
            if edge == 0:
                return((x, y, l + 1), 2)
            if edge == 1:
                return((x, y, l + 2), 2)
            if edge == 2:
                return((x, y + 1, l + 1), 0)
            if edge == 3:
                return((x - 1 , y + 1, l + 2), 0)

SG = SimilarityGroup(QQ)
def default_position(label):
    x, y, l = label
    if(ZZ(l) == 0):
        return SG(2*x, 2*y)  # (b + c) x, (a + c) y
    if(ZZ(l) == 1):
        return SG(2*x + 1, 2*y)  # (b + c) x + c, (a + c) y
    if(ZZ(l) == 2):
        return SG(2*x, 2*y + 1)  # (b + c) x, (a + c) y + c
    # Reminder to parameterize a, b, c here for positions.
    # Rework this to work for surfaces of different sizes.
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
s = SimilaritySurface(CubeSurf(QQ))
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
MM = matrix(QQ,[[0, 1, 0],
                [-1, 0, 0],
                [0, 0, 1]
])
def transformation(label):
    M = MatrixSpace(QQ, 3, 2)
    V = VectorSpace(QQ, 3)
    x, y, l = label
    if l == 0:
        return MM*M([[1, 0], [0, 1], [0, 0]]), MM*V([x, y, -x - y])
    elif l == 1:
        return MM*M([[0, 0], [0, 1], [-1, 0]]), MM*V([x + 1, y, -x - y])
    else:  # l == 2
        return MM*M([[1, 0], [0, 0], [0, -1]]), MM*V([x, y + 1, -x - y])
m = SurfaceToSpaceMapping(s, transformation)
def label_to_color(label):
    if label[2] == 0:
        return "pink"
    if label[2] == 1:
        return "yellow"
    if label[2] == 2:
        return "beige"
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
it = s.label_iterator()
m.plot({next(it) for i in range(30)}, label_to_color=label_to_color)
```

+++ {"deletable": true, "editable": true}

<b>Theorem (Pavel Javornik).</b>
A trajectory of rational slope (measured on one of the squares interpreted to
have horizontal and vertical sides) on the Necker Cube Surface closes up
if and only if the slope can be expressed as the ratio of two odd integers.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: true
---
B = s.tangent_bundle()
```

+++ {"deletable": true, "editable": true}

The following builds a trajectory starting in the base polygon at the point
$(\frac{1}{4}, \frac{1}{4})$ and traveling in a direction of slope one.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
v = B(s.base_label(), (1/4, 1/4), (-1, 1))
traj = v.straight_line_trajectory()
traj.flow(100)
if traj.is_closed():
    print("The trajectory closed up.")
labels = [seg.polygon_label() for seg in traj.segments()]
m.plot(labels, label_to_color=label_to_color) + line3d(m(traj), radius=0.02)
```

+++ {"deletable": true, "editable": true}

A trajectory of slope $5/4$.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
v = B(s.base_label(), (1/3, 1/4), (4, 5))
traj = v.straight_line_trajectory()
traj.flow(50)
labels = [seg.polygon_label() for seg in traj.segments()]
p = m.plot(labels, label_to_color=label_to_color) + line3d(m(traj),
           radius=0.04, label_to_color=label_to_color)
p.frame_aspect_ratio(tuple(vector(p.bounding_box()[1])
                           - vector(p.bounding_box()[0])))
p
```

+++ {"deletable": true, "editable": true}

A trajectory of slope $11/9$

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
v = B(s.base_label(), (1/3, 1/4), (9, 11))
traj = v.straight_line_trajectory()
traj.flow(1000)
while not traj.is_closed():
    traj.flow(1000)
labels = [seg.polygon_label() for seg in traj.segments()]
p = m.plot(labels, label_to_color=label_to_color)
p += line3d(m(traj), radius=0.04)
# p
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
show(p, frame=False, viewer="tachyon")
```

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
# show(p, frame=False)
```
