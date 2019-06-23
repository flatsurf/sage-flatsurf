from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.structure.sage_object import SageObject

class FlowPolygonMap(SageObject):
    r"""
    The map obtained as the return map of the flow on the sides of a (convex)
    polygon.

    Formally, this can be defined as follows: one start with two partition into
    (finitely many) intervals of a given interval. The map corresponds to
    changing from the first partition to the second. In other words, points are
    identified as pairs ``(i, x)`` where ``i`` is an atom and ``x`` is the
    relative position in this atom.

    Contrarily to an interval exchange transformation, things here are going
    from bottom to top.

    Note that this could also be used for homothetic surfaces. And to some
    extent to half translation surface.

    EXAMPLES::

        sage: from flatsurf.geometry.interval_exchange_transformation import FlowPolygonMap
        sage: T = FlowPolygonMap(QQ, [0,1,2], [2,3,1], [2,1,0], [1,3,2])
        sage: [T.forward_image(0,x) for x in range(3)]
        [(2, 0), (1, 0), (1, 1)]
        sage: [T.forward_image(1,x) for x in range(4)]
        [(1, 1), (1, 2), (0, 0), (0, 1)]
        sage: [T.forward_image(2,x) for x in range(1)]
        [(0, 1)]
    """
    def __init__(self, ring, bot_labels, bot_lengths, top_labels, top_lengths):
        r"""
        INPUT:

        - ``ring`` -- the base ring for the lengths of the interval

        - ``bot_labels`` -- labels for the bottom partition

        - ``bot_lengths`` -- lengths for the bottom partition

        - ``top_labels`` -- labels for the top partition

        - ``top_lengths`` -- lengths for top partition
        """
        assert all(x > ring.zero() for x in bot_lengths)
        assert all(x > ring.zero() for x in top_lengths)
        assert len(bot_labels) == len(bot_lengths)
        assert len(top_labels) == len(top_lengths)
        assert sum(top_lengths)  == sum(bot_lengths)

        self._ring = ring

        self._bot_labels = bot_labels
        self._top_labels = top_labels
        self._bot_labels_to_index = {j:i for i,j in enumerate(bot_labels)}
        self._top_labels_to_index = {j:i for i,j in enumerate(top_labels)}
        if len(self._bot_labels) != len(self._bot_labels_to_index):
            raise ValueError("non unique labels for bot: {}".format(bot_labels))
        if len(self._top_labels) != len(self._top_labels_to_index):
            raise ValueError("non unique labels in top: {}".format(top_labels))

        self._bot_lengths = list(map(ring,bot_lengths))
        self._top_lengths = list(map(ring,top_lengths))


        # forward image of intervals
        it = 0
        lt = x1 = top_lengths[it]
        self._forward_images = []
        for ib in range(len(bot_lengths)):
            lenb = bot_lengths[ib]
            self._forward_images.append((it,lt-x1))

            while lenb and lenb >= x1:
                lenb -= x1
                it += 1
                if it < len(top_lengths):
                    lt = x1 = top_lengths[it]
                else:
                    lt = ring.zero()
            if lenb:
                x1 -= lenb

        # backward image of intervals
        ib = 0
        lb = x1 = bot_lengths[ib]
        self._backward_images = []
        for it in range(len(top_lengths)):
            lent = top_lengths[it]
            self._backward_images.append((ib,lb-x1))

            while lent and lent >= x1:
                lent -= x1
                ib += 1
                if ib < len(bot_lengths):
                    lb = x1 = bot_lengths[ib]
                else:
                    lb = ring.zero()
            if lent:
                x1 -= lent

    def length_bot(self, i):
        i = self._bot_labels_to_index[i]
        return self._bot_lengths[i]

    def length_top(self, i):
        i = self._top_labels_to_index[i]
        return self._top_lengths[i]

    def _repr_(self):
        s = ["Flow polygon map:"]
        s.append(" " + " ".join(str(x) for x in self._top_labels))
        s.append(" " + " ".join(str(x) for x in self._bot_labels))
        s.append("top lengths: {}".format(self._top_lengths))
        s.append("bot lengths: {}".format(self._bot_lengths))
        return "\n".join(s)

    def forward_image(self, i, x):
        r"""
        Return the forward image.

        EXAMPLES::

            sage: from flatsurf.geometry.interval_exchange_transformation import FlowPolygonMap

        Singularities are always sent to a ``(i,0)``::

            sage: T = FlowPolygonMap(QQ, [0,1], [2,1], [2,3,4], [1,1,1])
            sage: T.forward_image(0, 0)
            (2, 0)
            sage: T.forward_image(0, 1)  # could have equally been (2, 1)
            (3, 0)
        """
        i = self._bot_labels_to_index[i]
        if x < self._ring.zero() or x > self._bot_lengths[i]:
            raise ValueError("x = {} is out of the interval".format(x))
        j,y = self._forward_images[i]
        if x+y < self._top_lengths[j]:
            return (self._top_labels[j], x+y)
        x -= self._top_lengths[j]-y
        j += 1
        while x > self._top_lengths[j]:
            x -= self._top_lengths[j]
            j += 1
        return (self._top_labels[j],x)

    def backward_image(self, i, x):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.interval_exchange_transformation import \
            ....:     FlowPolygonMap

            sage: x = polygen(ZZ)
            sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())

            sage: T = FlowPolygonMap(K, [0,1,2],
            ....:       [sqrt2,1+sqrt2,1], [2,0,1], [1,sqrt2,1+sqrt2])
            sage: T.backward_image(*T.forward_image(1, 1))
            (1, 1)
            sage: T.backward_image(*T.forward_image(0, 1))
            (0, 1)
            sage: T.backward_image(*T.forward_image(0, sqrt2-1))
            (0, sqrt2 - 1)
            sage: T.backward_image(*T.forward_image(1, sqrt2-1))
            (1, sqrt2 - 1)

        Singularities are always sent to a ``(i,0)``::

            sage: T = FlowPolygonMap(QQ, [2,3,4], [1,1,1], [0,1], [2,1])
            sage: T.backward_image(0, 0)
            (2, 0)
            sage: T.backward_image(0, 1)  # could have equally been (2, 1)
            (3, 0)

        TESTS::

            sage: T = FlowPolygonMap(K, [0,1,2],
            ....:       [5*sqrt2,1,1], [2,1,0], [1,1,5*sqrt2])
            sage: for x in range(1,8):
            ....:     p0 = (0,x)
            ....:     for n in range(5):
            ....:         p1 = T.forward_image(*p0)
            ....:         assert T.backward_image(*p1) == p0, "p0 = {}, p1 = {}".format(p0,p1)
            ....:         p0 = p1

        """
        i = self._top_labels_to_index[i]
        if x < self._ring.zero() or x > self._top_lengths[i]:
            raise ValueError("x = {} is out of the interval".format(x))
        j,y = self._backward_images[i]
        if x+y < self._bot_lengths[j]:
            return (self._bot_labels[j], x+y)
        x -= self._bot_lengths[j]-y
        j += 1
        while x > self._bot_lengths[j]:
            x -= self._bot_lengths[j]
            j += 1
        return (self._bot_labels[j],x)
