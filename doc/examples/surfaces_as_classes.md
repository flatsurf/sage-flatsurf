---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.1
kernelspec:
  display_name: SageMath 9.7
  language: sage
  name: sagemath
---

# Building finite and infinite surfaces as classes

The ready-made surfaces and the gluing of polygons one by one are described in
[Defining Surfaces](./defining_surfaces). A surface with infinitely many
polygons cannot be built that way, and neither can a surface whose polygons
are better computed than stored. Such a surface is written as a class: instead
of storing its polygons, the class says how to compute the polygon with a
given label and what lies on the other side of each of its edges. This page is
more technical than the others; it explains what such a class must provide and
walks through two examples, a flat one and a hyperbolic one.

## What a class must implement

A flat surface is a subclass of `OrientedSimilaritySurface`, a hyperbolic one a
subclass of `OrientedHyperbolicSurface`. In both cases the constructor takes
the base ring of the surface and its category:
`super().__init__(base_ring, category=...)`.

The polygons of the surface are indexed by labels, which can be any hashable
objects: integers, tuples, strings, ... The following methods describe the
surface, the mandatory ones first.

- **`roots()`** returns a tuple with one label for each connected component of
  the surface. The polygons of a component are found by walking from its
  root through the gluings, so that the roots are where every iteration over
  the labels starts.
- **`polygon(label)`** returns the polygon with this label, a Euclidean polygon
  for a flat surface and a hyperbolic polygon for a hyperbolic one. Its edges
  are numbered counterclockwise, the edge `i` going from the vertex `i` to the
  vertex `i + 1`.
- **`opposite_edge(label, edge)`** returns the pair `(label, edge)` of the edge
  glued to the given one, or `None` for an edge without a neighbour, a
  boundary edge. Gluing is symmetric: the opposite of the opposite edge is the
  edge we started from.
- **`is_mutable()`** returns whether the surface can be modified. A surface
  that computes its polygons is not mutable: it returns `False`.
- **`is_compact()`** returns whether the surface is compact. A surface with
  infinitely many polygons is not: it returns `False`.
- **`__eq__(other)` and `__hash__()`** make two surfaces built by the same
  class and the same data equal. Without them a surface is only equal to
  itself. As soon as `__eq__` is defined, `__hash__` must be defined as well,
  or the surface cannot be hashed.

The other methods are optional:

- **`__repr__()`** returns the text printed for the surface.
- **`graphical_surface()`** returns the object used to plot a flat surface;
  an infinite surface can make a few of its polygons visible by default.

Nothing else is needed: the other methods of a surface (`labels()`, `edges()`,
`angles()`, ...) are provided by the category of the surface, and they use
only the methods above. The testing framework of Sage, `TestSuite`, checks
that the mandatory methods exist and that the surface is consistent; it is
the first thing to run on a new class.

## The role of the categories

The category that is passed to `super().__init__` states what the surface is,
a translation surface or a hyperbolic cone surface for instance, and which
axioms it satisfies:

- `InfiniteType` says that the surface has infinitely many polygons,
- `WithoutBoundary` that every edge has a neighbour,
- `Connected` that the surface has a single component.

The methods that a surface offers come from its category. A category that is
too weak gives fewer methods: only a translation surface has the methods of
translation surfaces. A category that is wrong gives wrong answers, since the
methods trust it. For a finite surface, the category can be computed from the
polygons and the gluings, which is what happens when a mutable surface is made
immutable. For an infinite surface nothing can be computed, since only finitely
many polygons can ever be looked at, so the category is stated explicitly. The
examples below show the category of the surface with `category()`.

## A flat example: the parabola surface

The following translation surface is made of infinitely many polygons that
stack up inside the parabola $y = x^2$. The polygon with label $1$ is the
triangle with vertices $(0, 0)$, $(1, 1)$ and $(-1, 1)$. For a label $n > 1$ it
is the quadrilateral of the points inside the parabola between the heights
$(n-1)^2$ and $n^2$. A negative label $-n$ is the polygon $n$ rotated by a
half-turn. The class implements the mandatory methods above, and
`graphical_surface()` so that the first six polygons are shown when plotting.
We state the category explicitly: a translation surface of infinite type,
without boundary and connected.

```{code-cell}
from sage.all import QQ, ZZ, matrix, TestSuite

from flatsurf import Polygon
from flatsurf.geometry.surface import OrientedSimilaritySurface
from flatsurf.geometry.categories import TranslationSurfaces


class ParabolaSurface(OrientedSimilaritySurface):
    def __init__(self):
        # For finite surfaces, the category can be determined automatically
        # but for infinite surfaces, we need to make an explicit choice here.
        super().__init__(
            QQ,
            category=TranslationSurfaces().InfiniteType().WithoutBoundary().Connected(),
        )

    def __repr__(self):
        r"""
        Return a printable representation of this surface.
        """
        return "ParabolaSurface()"

    def roots(self):
        r"""
        Return a label for each connected component of the surface.

        Iterating the polygons of the connected component starts at these
        labels.
        """
        return (1,)

    def is_mutable(self):
        r"""
        Return whether this surface can be modified by the user.
        """
        return False

    def is_compact(self):
        r"""
        Return whether this surface is a compact space.
        """
        return False

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.
        """
        return isinstance(other, ParabolaSurface)

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with ``__eq__``.
        """
        return hash(type(self))

    def graphical_surface(self, **kwds):
        r"""
        Return a plottable representation of this surface.
        """
        graphical_surface = super().graphical_surface(**kwds)
        # Make the first six polygons of the surface visible by default when
        # plotting.
        graphical_surface.make_all_visible(limit=6)
        return graphical_surface

    def polygon(self, label):
        r"""
        Return the polygon making up this surface labeled ``label``.
        """
        if label not in ZZ or label == 0:
            raise ValueError(f"invalid label {label}")

        if label < 0:
            return matrix(QQ, [[-1, 0], [0, -1]]) * self.polygon(-label)

        if label == 1:
            return Polygon(vertices=[(0, 0), (1, 1), (-1, 1)], base_ring=QQ)

        return Polygon(
            vertices=[
                (label - 1, (label - 1) ** 2),
                (label, label**2),
                (-label, label**2),
                (-label + 1, (label - 1) ** 2),
            ],
            base_ring=QQ,
        )

    def opposite_edge(self, label, e):
        if label not in ZZ or label == 0:
            raise ValueError(f"invalid label {label}")
        if label in [-1, 1] and e not in [0, 1, 2]:
            raise ValueError("no such edge")
        if e not in [0, 1, 2, 3]:
            raise ValueError("no such edge")

        if label in [-1, 1] and e == 1:
            return 2 * label, 3

        if e in [0, 2]:
            return -label, e

        if e == 1:
            return label + label.sign(), 3

        return label - label.sign(), 1
```

```{code-cell}
s = ParabolaSurface()
s
```

```{code-cell}
s.category()
```

```{code-cell}
s.plot()
```

We can run a test suite to ensure that we have implemented everything that is
needed to make this a fully functional surface.

```{code-cell}
TestSuite(s).run(verbose=True)
```

## A hyperbolic example: a tree of pants

The hyperbolic example is an infinite surface, as above but in the hyperbolic
plane: pairs of pants glued along their boundary components so that the gluing
pattern is the infinite 3-regular tree. A pair of pants of
`hyperbolic_surfaces.pants` is two right-angled hexagons, with labels 0 (the
front one) and 1 (its mirror image). They are glued along three edges, the
*seams*; the three other edges of each hexagon are not glued. The unglued
edges form the three boundary geodesics of the pants, the *cuffs*, each made of
one edge of each hexagon. Here the three cuffs have the same trace $t = 21/4$:

```{code-cell}
from sage.all import QQ

from flatsurf import hyperbolic_surfaces

P = hyperbolic_surfaces.pants(boundary_traces=(QQ(21) / 4,) * 3)
P.gluings(), P.boundary_components()
```

The gluings are the seams: the edge `(0, 0)` of the front hexagon is glued to
the edge `(1, 0)` of the mirror one, and similarly for the edges 2 and 4. The
edges 1, 3 and 5 of each hexagon are unglued, and `boundary_components()`
lists the cuffs. The first cuff is made of the edge `(0, 1)` and the edge
`(1, 5)`, so that the edges of the mirror hexagon come in the reverse order of
those of the front one.

To get the tree, we take one copy of the pants for each vertex of the
3-regular tree, and we glue the cuff number `i` of a copy to the cuff number
`i` of the neighbouring copy in direction `i`. The labels of the surface are
the pairs `(word, side)`: `word` is a vertex of the tree, a tuple over
$\{0, 1, 2\}$ without two equal consecutive letters (the neighbour across the
cuff `i` appends `i` to the word, or removes it when it is its last letter),
and `side` is 0 or 1, the hexagon of the pants. The polygon of a label is the
hexagon of its side. An edge of a hexagon is either a seam, whose opposite is
in the same pants, or a cuff, whose opposite is the same edge of the same
hexagon of the neighbouring pants. It is the same edge and not another edge of
the cuff: the gluing is the half-turn about the midpoint of the edge, so that
there is no twist between the two pants. The base class is
`OrientedHyperbolicSurface`, and the category is the one of a hyperbolic cone
surface of infinite type, without boundary and connected.

```{code-cell}
from flatsurf.geometry.categories import HyperbolicConeSurfaces
from flatsurf.geometry.hyperbolic_surface import OrientedHyperbolicSurface


class TreeOfPants(OrientedHyperbolicSurface):
    def __init__(self):
        self._pants = hyperbolic_surfaces.pants(
            boundary_traces=(QQ(21) / 4, QQ(21) / 4, QQ(21) / 4)
        )
        # The cuff of each unglued edge of the pants: 0, 1 or 2.
        self._cuff = {
            side_edge: cuff
            for cuff, component in enumerate(self._pants.boundary_components())
            for side_edge in component
        }
        super().__init__(
            self._pants.hyperbolic_plane(),
            category=HyperbolicConeSurfaces()
            .Oriented()
            .InfiniteType()
            .WithoutBoundary()
            .Connected(),
        )

    def __repr__(self):
        r"""
        Return a printable representation of this surface.
        """
        return "TreeOfPants()"

    def roots(self):
        r"""
        Return a label for each connected component of the surface.

        The surface is connected: its root is the pants at the root of the
        tree.
        """
        return (((), 0),)

    def is_mutable(self):
        r"""
        Return whether this surface can be modified by the user.
        """
        return False

    def is_compact(self):
        r"""
        Return whether this surface is a compact space.
        """
        return False

    def __eq__(self, other):
        r"""
        Return whether this surface is indistinguishable from ``other``.
        """
        return isinstance(other, TreeOfPants)

    def __hash__(self):
        r"""
        Return a hash value for this surface that is compatible with ``__eq__``.
        """
        return hash(type(self))

    def polygon(self, label):
        r"""
        Return the polygon making up this surface labeled ``label``.
        """
        word, side = label
        return self._pants.polygon(side)

    def opposite_edge(self, label, edge):
        r"""
        Return the label and the edge glued to the given edge of ``label``.
        """
        word, side = label
        opposite = self._pants.opposite_edge(side, edge)
        if opposite is not None:
            # a seam: the other hexagon of the same pants
            return (word, opposite[0]), opposite[1]
        # a cuff: the same edge of the neighbouring pants
        cuff = self._cuff[side, edge]
        if word and word[-1] == cuff:
            return (word[:-1], side), edge
        return (word + (cuff,), side), edge
```

```{code-cell}
S = TreeOfPants()
S
```

The category states what we told the class, and the hyperbolic plane is that
of the pants, the plane of a quadratic form of its own over its base ring:

```{code-cell}
S.category()
```

```{code-cell}
S.hyperbolic_plane()
```

The surface has infinitely many labels, which are enumerated starting from the
root. Here are the first eight:

```{code-cell}
import itertools

list(itertools.islice(S.labels(), 8))
```

The edge 0 of the front hexagon at the root is a seam, glued to the mirror
hexagon of the same pants:

```{code-cell}
S.opposite_edge(((), 0), 0)
```

The edge 1 is a cuff, glued to the front hexagon of the pants in direction 0:

```{code-cell}
S.opposite_edge(((), 0), 1)
```

Crossing a cuff twice leads back to the starting pants. From the pants `(1,)`,
the edge 3 leads back to the root, since the word loses its last letter:

```{code-cell}
S.opposite_edge(((1,), 0), 3)
```

The gluing across a cuff is the half-turn about the midpoint of the edge; the
trace of its matrix in $\mathrm{SO}(2, 1)$ is $-1$:

```{code-cell}
S.edge_transformation(((), 0), 1).trace()
```

The vertices of the hexagons are right angles, and four of them meet at a
vertex of the surface, which is not a cone point: its angle is `1`, that is
$2\pi$:

```{code-cell}
S(((), 0), 0).angle()
```

Finally the test suite checks the consistency of the surface:

```{code-cell}
TestSuite(S).run()
```

## When the gluings are wrong

Nothing of what a class returns is checked at construction. The gluing must be
an involution, `opposite_edge` applied twice gives the edge we started from,
and the walk around a vertex of the surface, which jumps from a polygon to the
next one through the gluings, must close. On a finite surface, the walk
visits finitely many polygons. On an infinite surface, nothing bounds it, so
that a mismatch makes the construction of a vertex run forever. For instance,
matching the cuffs of the mirror hexagon in the order of the front one, the
edges 1, 3, 5 on both sides, instead of the reversed order returned by
`boundary_components()`, produces a surface where the walk around a vertex
never closes, and asking for the angle `S(label, 0)` does not return. Running
the `TestSuite` on a new class is the first check.
