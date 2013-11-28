r"""
Class to store and manipulate edge identifications between polygons.

TODO:

- The current structure only label external edges but it is sometimes convenient
  to label all of them. Especially when one consider cutting sequences in
  Bouw-Moeller lattice surfaces.

- Sometimes we want labels to be fixed (for example on an origami where the
  squares are labeld 1,2,... it is convenient to use a1, b1, a2, b2,...).
"""
from similarity_surface import *

lower_and_upper_case = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
class InfiniteAlphabet(Parent, UniqueRepresentation):
    r"""
    EXAMPLES::

        sage: I = InfiniteAlphabet(); I
        {a, b, c, ..., A, B, ..., aa, ab, ... }
        sage: it = iter(I)
        sage: ' '.join(it.next() for _ in xrange(26))
        'a b c d e f g h i j k l m n o p q r s t u v w x y z'
        sage: ' '.join(it.next() for _ in xrange(26))
        'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'
        sage: ' '.join(it.next() for _ in xrange(26))
        'aa ab ac ad ae af ag ah ai aj ak al am an ao ap aq ar as at au av aw ax ay az'
        sage: for _ in xrange(2**16): it.next();
        sage: it.next()
        'xmQ'
    """
    def __init__(self):
        from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        return '{a, b, c, ..., A, B, ..., aa, ab, ... }'

    def __contains__(self, elt):
        if not isinstance(elt, str):
            return False
        for i in elt:
            if ord(i) < 65 or ord(i) > 122:
                return False
            if ord(i) > 90 and ord(i) < 97:
                return False
        return True

    def _rank_one_letter(self, letter):
        r"""
        TESTS::

            sage: I = InfiniteAlphabet()
            sage: I._rank_one_letter('a')
            0
            sage: I._rank_one_letter('A')
            26
            sage: I._rank_one_letter('Z')
            51
        """
        # ord('a') = 97
        # ord('z') = 122
        # ord('A') = 65
        # ord('Z') = 90
        n = ord(letter)
        if n > 96:
            if n > 122:
                raise ValueError("%s not in self"%letter)
            return n-97
        if n < 65 or n > 90:
            raise ValueError("%s not in self"%letter)
        return 26+n-65

    def __iter__(self):
        from itertools import product,count

        for n in count(1):
            for word in product(lower_and_upper_case, repeat=n):
                yield ''.join(word)

class EdgeLabels:
    r"""
    Abstract class for storing labels for edges in a similarity surface (or similar)

    The method "get_opposite" must be overriden.
    """
    def __init__(self, label_set=None):
        r"""
        INPUT:

        - ``label_set`` - an optional set to be used for labels. It can be any
          iterable (a list, a string, ...).
        """
        self._labels={}         # dictionnary: (p1,e1) -> label and (p2,e2) -> label
        self._label_to_edge={}  # dictionnary: label -> (p1,e1)
        if label_set is None:
            label_set = InfiniteAlphabet()
        self._label_set = label_set

    def get_edge(self,label):
        return self._label_to_edge.get(label)

    def get_opposite(self, p1, e1):
        raise NotImplementedError("get_opposite must be overriden in subclasses")

    def remove_label(self, label):
        if self._label_to_edge.has_key(label):
            pair1=self.get_edge(label)
            pair2=self.get_opposite(pair1[0], pair1[1])
            del self._labels[pair1]
            del self._labels[pair2]
            del self._label_to_edge[label]

    def get_edge_pair_list(self):
        res=[]
        for pair1, pair2 in self._gluings.iteritems():
            res.append( (pair1,pair2) )
        return res

    def get_label(self,p1,e1):
        return self._labels.get( (p1,e1) )

    def get_unused_label(self):
        for label in self._label_set:
            if label not in self._label_to_edge:
                return label
        raise ValueError("All available labels already used!")

    def add_label(self, p1, e1, label=None):
        pair2=self.get_opposite(p1,e1)
        if pair2 is not None:
            if label is None:
               label=self.get_unused_label()
            self._labels[(p1,e1)]=label
            self._labels[pair2]=label
            self._label_to_edge[label]=(p1,e1)

class SurfaceLabels(EdgeLabels):
    r"""
    Class for storing labels for edges in a similarity surface.
    """
    def __init__(self, similarity_surface):
        AbstractLabels.__init__(self)
        self._ss=similarity_surface

    def get_opposite(self, p1, e1):
        return self._ss.opposite_edge(p1,e1)

    def get_label(self,p1,e1):
        res=AbstractLabels.get_label(self,p1,e1)
        if res is None:
            self.add_label(p1,e1)
            res=AbstractLabels.get_label(self,p1,e1)
        return res

class EditableEdgeGluing(EdgeLabels):
    r"""
    Class for manipulating gluings and labels of a similiarity surface.
    """
    def __init__(self):
        AbstractLabels.__init__(self)
        self._gluings={}

    def get_opposite(self, p1, e1):
        return self._gluings.get( (p1,e1) )

    def remove_glue(self, p1, e1):
        pair2=self.get_opposite(p1,e1)
        if pair2 is None:
            return
        label=self.get_label(p1,e1)
        if label is not None:
            self.remove_label(label)
        del self._gluings[(p1,e1)]
        del self._gluings[pair2]

    def get_edge_pair_list(self):
        res=[]
        for pair1, pair2 in self._gluings.iteritems():
            res.append( (pair1,pair2) )
        return res

    def glue(self, p1, e1, p2, e2):
        if (p1 != p2) or (e1 != e2):
            self.remove_glue(p1,e1)
            self.remove_glue(p2,e2)
            self._gluings[(p1,e1)]=(p2,e2)
            self._gluings[(p2,e2)]=(p1,e1)
        else:
            raise ValueError("The edges being identified must be distinct!")

    def glue_and_label(self, p1, e1, p2, e2, label=None):
        self.glue(p1, e1, p2, e2)
        self.add_label(p1, e1, label=label)

