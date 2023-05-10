r"""
Translation Surfaces.
"""

from sage.rings.rational_field import QQ
from flatsurf.geometry.surface import OrientedSimilaritySurface


class AbstractOrigami(OrientedSimilaritySurface):
    r"""
    Abstract base class for (connected) origamis.

    Realization needs just to define a _domain and four cardinal directions.
    """

    def __init__(self, domain, root=None, base_label=None, category=None):
        self._domain = domain

        if base_label is not None:
            import warnings
            warnings.warn("base_label has been deprecated as a keyword argument for AbstractOrigami and will be removed in a future version of sage-flatsurf; use root instead")
            root = base_label
            base_label = None

        if root is None:
            root = domain.an_element()
        self._root = root

        from flatsurf.geometry.categories import TranslationSurfaces

        if category is None:
            category = TranslationSurfaces()

        category &= TranslationSurfaces().WithoutBoundary().Connected()

        finite = domain.is_finite()
        if finite:
            category &= category.FiniteType()
        else:
            category &= category.InfiniteType()

        from flatsurf.geometry.polygon import polygons
        self._square = polygons.square()

        super().__init__(QQ, category=category)

    def roots(self):
        return (self._root,)

    def labels(self):
        from flatsurf.geometry.surface import LabelsView

        return LabelsView(self, self._domain)

    def is_mutable(self):
        return False

    def up(self, label):
        raise NotImplementedError

    def down(self, label):
        raise NotImplementedError

    def right(self, label):
        raise NotImplementedError

    def left(self, label):
        raise NotImplementedError

    def _repr_(self):
        return "Some AbstractOrigami"

    def polygon_labels(self):
        return self._domain

    def polygon(self, lab):
        if lab not in self._domain:
            # Updated to print a possibly useful error message
            raise ValueError("Label " + str(lab) + " is not in the domain")

        return self._square

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError
        if e == 0:
            return self.down(p), 2
        if e == 1:
            return self.right(p), 3
        if e == 2:
            return self.up(p), 0
        if e == 3:
            return self.left(p), 1
        raise ValueError


class Origami(AbstractOrigami):
    def __init__(
        self, r, u, rr=None, uu=None, domain=None, root=None, base_label=None, category=None
    ):
        if domain is None:
            domain = r.parent().domain()

        self._r = r
        self._u = u
        if rr is None:
            rr = ~r
        else:
            for a in domain.some_elements():
                if r(rr(a)) != a:
                    raise ValueError("r o rr is not identity on %s" % a)
                if rr(r(a)) != a:
                    raise ValueError("rr o r is not identity on %s" % a)
        if uu is None:
            uu = ~u
        else:
            for a in domain.some_elements():
                if u(uu(a)) != a:
                    raise ValueError("u o uu is not identity on %s" % a)
                if uu(u(a)) != a:
                    raise ValueError("uu o u is not identity on %s" % a)

        self._perms = [uu, r, u, rr]  # down,right,up,left
        AbstractOrigami.__init__(self, domain, root=root, base_label=base_label, category=category)

    def opposite_edge(self, p, e):
        if p not in self._domain:
            raise ValueError(
                "Polygon label p=" + str(p) + " is not in domain=" + str(self._domain)
            )
        if e < 0 or e > 3:
            raise ValueError("Edge value e=" + str(e) + " does not satisfy 0<=e<4.")
        return self._perms[e](p), (e + 2) % 4

    def up(self, label):
        return self.opposite_edge(label, 2)[0]

    def down(self, label):
        return self.opposite_edge(label, 0)[0]

    def right(self, label):
        return self.opposite_edge(label, 1)[0]

    def left(self, label):
        return self.opposite_edge(label, 3)[0]

    def _repr_(self):
        return "Origami defined by r=%s and u=%s" % (self._r, self._u)

    def __eq__(self, other):
        if not isinstance(other, Origami):
            return False

        return (
            self._perms == other._perms
            and self._domain is other._domain
            and self.roots() == other.roots()
        )

    def __hash__(self):
        return hash((Origami, tuple(self._perms), self._domain, self.roots()))


class LazyStandardizedPolygonSurface(OrientedSimilaritySurface):
    r"""
    This class handles standardizing polygons for infinite translation surfaces.
    See the TranslationSurface.standardize_polygons method.

    This class should not be instantiated directly.
    Instead use TranslationSurface.standardize_polygons.
    """

    def __init__(self, surface, relabel=None, category=None):
        if relabel is not None:
            if relabel:
                raise NotImplementedError(
                    "the relabel keyword has been removed from LazyStandardizedPolygonSurface; use relabel({old: new for (new, old) in enumerate(surface.labels())}) to use integer labels instead"
                )
            else:
                import warnings

                warnings.warn(
                    "the relabel keyword will be removed in a future version of sage-flatsurf; do not pass it explicitly anymore to LazyStandardizedPolygonSurface"
                )

        from flatsurf.geometry.surface import MutableOrientedSimilaritySurface

        self._s = MutableOrientedSimilaritySurface.from_surface(surface)
        self._labels = set()

        super().__init__(self._s.base_ring(), category=category or self._s.category())

    def roots(self):
        return self._s.roots()

    def is_mutable(self):
        return False

    def standardize(self, label):
        best = 0
        polygon = self._s.polygon(label)
        best_pt = polygon.vertex(best)
        for v in range(1, polygon.num_edges()):
            pt = polygon.vertex(v)
            if (pt[1] < best_pt[1]) or (pt[1] == best_pt[1] and pt[0] < best_pt[0]):
                best = v
                best_pt = pt
        if best != 0:
            self._s.set_vertex_zero(label, best, in_place=True)
        self._labels.add(label)

    def polygon(self, label):
        r"""
        Return the polygon with the provided label.

        This method must be overridden in subclasses.
        """
        if label in self._labels:
            return self._s.polygon(label)
        else:
            self.standardize(label)
            return self._s.polygon(label)

    def opposite_edge(self, label, e):
        r"""
        Given the label ``label`` of a polygon and an edge ``e`` in that
        polygon returns the pair (``ll``, ``ee``) to which this edge is glued.
        """
        if label not in self._labels:
            self.standardize(label)
        ll, ee = self._s.opposite_edge(label, e)
        if ll in self._labels:
            return (ll, ee)
        self.standardize(ll)
        return self._s.opposite_edge(label, e)
