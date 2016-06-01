from sage.misc.cachefunc import cached_method

from sage.structure.element import parent

from sage.categories.groups import Groups
from sage.groups.group import Group
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation


class Path(MultiplicativeGroupElement):
# activating the following somehow break the discovery of the Python _mul_
# method below...
#    __slots__ = ['_polys', '_edges', '_edges_rev']

    def __init__(self, parent, polys, edge, edge_rev):
        self._polys = tuple(polys)
        self._edges = tuple(edge)
        self._edges_rev = tuple(edge_rev)
        MultiplicativeGroupElement.__init__(self, parent)
        self._check()

    def __hash__(self):
        return hash(self._polys) ^ hash(self._edges)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group()
            sage: a,b = F.gens()
            sage: a == b
            False
            sage: a*b == b*a
            False
            sage: a*b == a*b
            True
        """
        return parent(self) is parent(other) and \
               self._polys == other._polys and  \
               self._edges == other._edges 

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: t = translation_surfaces.square_torus()
            sage: F = t.fundamental_group()
            sage: a,b = F.gens()
            sage: a != b
            True
            sage: a*b != b*a
            True
            sage: a*b != a*b
            False
        """
        return parent(self) is not parent(other) or \
               self._polys != other._polys or \
               self._edges != other._edges

    def _check(self):
        if not(len(self._polys)-1 == len(self._edges) == len(self._edges_rev)):
            raise ValueError("polys = {}\nedges = {}\nedges_rev={}".format(
                self._polys, self._edges, self._edges_rev))
        assert self._polys[0] == self.parent()._b == self._polys[-1]

    def is_one(self):
        return not self._edges

    def _repr_(self):
        return " ".join("{} --{}--".format(p,e)
              for p,e in zip(self._polys, self._edges)) + \
               "{}".format(self._polys[-1]) 

    def _mul_(self, other):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: t = translation_surfaces.square_torus()
            sage: a,b = t.fundamental_group().gens()
            sage: a*b
            0 --0-- 0 --1--0
        """
        sp = self._polys[:]
        se = self._edges[:]
        ser = self._edges_rev[:]

        op = other._polys[:]
        oe = other._edges[:]
        oer = other._edges_rev[:]

        if sp[-1] != op[0]:
            return None

        i = 0
        while i < len(se) and i < len(oe) and se[-i-1] == oer[i]:
            i += 1

        P = self.parent()
        return P.element_class(
                P,
                sp[:len(sp)-i] + op[i+1:],
                se[:len(se)-i]+ oe[i:],
                ser[:len(ser)-i] + oer[i:])

    def __invert__(self):
        r"""
        TESTS::

            sage: from flatsurf import *
            sage: o = translation_surfaces.octagon_and_squares()
            sage: F = o.fundamental_group()
            sage: a1,a2,a3,a4,a5,a6 = F.gens()
            sage: (a1 * a2 * ~a2 * ~a1).is_one()
            True
            sage: (a1 * ~a2 * a2 * a1) == a1 * a1
            True
        """
        P = self.parent()
        return P.element_class(
                P,
                self._polys[::-1],
                self._edges_rev[::-1],
                self._edges[::-1])

class FundamentalGroup(UniqueRepresentation, Group):
    r"""
    The fundamental group of a punctured surface

    EXAMPLES::

        sage: from flatsurf import *
        sage: t = translation_surfaces.square_torus()
        sage: TestSuite(t.fundamental_group()).run()

    """
    Element = Path

    def __init__(self, surface, base):
        if not surface.is_finite():
            raise ValueError("the method only work for finite surfaces")
        self._s = surface
        self._b = base

        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        return "Fundamental group of {} based at polygon {}".format(
                self._s,
                self._b)

    @cached_method
    def one(self):
        return self.element_class(self, [self._b], [], [])

    @cached_method
    def gens(self):
        r"""
        EXAMPLES::

            sage: from flatsurf import *
            sage: S = SymmetricGroup(8)
            sage: r = S('(1,2,3,4,5,6,7,8)')
            sage: u = S('(1,8,5,4)(2,3)(6,7)')
            sage: o = translation_surfaces.origami(r,u)
            sage: len(o.fundamental_group().gens())
            9
        """
        p = self._b
        s = self._s
        tree = {}   # a tree whose root is base_label
        basis = []

        tree[p] = (None,None,None)

        wait = [] # list of edges of the dual graph, ie p1 -- (e1,e2) --> p2
        for e in xrange(s.polygon(p).num_edges()):
            pp,ee = s.opposite_edge(p,e)
            wait.append((pp,ee,p,e))
        while wait:
            p1,e1,p2,e2 = wait.pop()
            assert p2 in tree
            if p1 in tree: # new cycle?
                if (p1,e1) > (p2,e2):
                    continue
                polys = [p1]
                edges = []
                edges_rev = []

                p1,e,e_back = tree[p1]
                while p1 is not None:
                    edges.append(e_back)
                    edges_rev.append(e)
                    polys.append(p1)
                    p1,e,e_back = tree[p1]
                polys.reverse()
                edges.reverse()
                edges_rev.reverse()

                polys.append(p2)
                edges.append(e1)
                edges_rev.append(e2)

                p2,e,e_back = tree[p2]
                while p2 is not None:
                    edges.append(e)
                    edges_rev.append(e_back)
                    polys.append(p2)
                    p2,e,e_back = tree[p2]

                basis.append((polys, edges, edges_rev))

            else: # new branch
                tree[p1] = (p2,e1,e2)
                for e in xrange(s.polygon(p1).num_edges()):
                    if e != e1:
                        pp,ee = s.opposite_edge(p1,e)
                        wait.append((pp,ee,p1,e))

        basis.sort()
        return tuple([self.element_class(self,p,e,er) for p,e,er in basis])

