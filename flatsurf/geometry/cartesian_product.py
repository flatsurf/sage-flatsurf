from sage.categories.cartesian_product import cartesian_product
CartesianProduct = cartesian_product

# for now we keep the code below, but we should get rid of it!!

from sage.rings.integer_ring import ZZ
from sage.structure.parent import Parent

ZZ_0 = ZZ(0)

class OLD_CartesianProduct(Parent):
    def __init__(self, sets, category=None):
        from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
        from sage.categories.enumerated_sets import EnumeratedSets

        self._sets = sets
        if category is None:
            if all(s in FiniteEnumeratedSets() for s in sets):
                category = FiniteEnumeratedSets()
            else:
                category = EnumeratedSets()
        Parent.__init__(self, category=category)

    def cardinality(self):
        from sage.misc.misc_c import prod
        if any(s.cardinality() == 0 for s in self._sets):
            return ZZ_0
        return prod(s.cardinality() for s in self._sets)

    def _repr_(self):
        return "Cartesian products of " + ", ".join(map(str,self._sets))

    def __iter__(self):
        r"""
        TESTS::

            sage: from flatsurf.geometry.cartesian_product import CartesianProduct
            sage: A2 = IntegerRange(2)
            sage: A3 = FiniteEnumeratedSet('abc')
            sage: A4 = Alphabet([-1,-2,-3,-4])
            sage: len(CartesianProduct([A2,A3]).list())
            6
            sage: len(CartesianProduct([A3,A2]).list())
            6

            sage: len(CartesianProduct([A2,A3,A4]).list())
            24
            sage: len(CartesianProduct([A3,A2,A4]).list())
            24
            sage: len(CartesianProduct([A4,A2,A3]).list())
            24
        """
        from itertools import product
        from sage.combinat.set_partition_ordered import OrderedSetPartitions

        n = len(self._sets)
        iterators = map(iter, self._sets)
        not_yet_terminated = [True] * len(iterators)
        values = [[] for _ in xrange(n)]
        last_values = [None] * len(iterators)
        s = [(list(s2),list(s1)) for s1,s2 in OrderedSetPartitions(range(n),2)]
        s.append((range(n),[]))
        l = [None] * n
        r = []

        while True:
            has_changed = False
            for i in xrange(n):
                if not_yet_terminated[i]:
                    try:
                        last_values[i] = iterators[i].next()
                    except StopIteration:
                        not_yet_terminated[i] = False
                        has_changed = True
            if has_changed:
                r = [i for i in xrange(n) if not not_yet_terminated[i]]
                rr = [i for i in xrange(n) if not_yet_terminated[i]]
                # first atom: the last values
                s = [(list(s2),list(s1)+r) for s1,s2 in OrderedSetPartitions(rr,2)]
                s.append((rr,r))

            if len(r) == n:
                return

            for s1,s2 in s:
                for j in s1:
                    l[j] = last_values[j]
                for p in product(*tuple(values[i] for i in s2)):
                    for i,j in enumerate(s2):
                        l[j] = p[i]
                    yield tuple(l)

            for i in xrange(n):
                if not_yet_terminated[i]:
                    values[i].append(last_values[i])

    def an_element(self):
        return tuple(s.an_element() for s in self._sets)

    def __contains__(self, elt):
        if not isinstance(elt,tuple):
            return False
        if len(elt) != len(self._sets):
            return False
        return all(elt[i] in self._sets[i] for i in xrange(len(self._sets)))
