r"""
Class for alphabet with inverses.

The main class is AlphabetWithInverses. See examples in there.
"""

#*****************************************************************************
#       Copyright (C) 2013 Thierry Coulbois <thierry.coulbois@univ-amu.fr>
#                     2013 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.combinat.words.alphabet import build_alphabet

def is_lower_character(a):
    r"""
    EXAMPLES::

        sage: filter(is_lower_character, map(chr,range(256)))
    """
    return isinstance(a,str) and len(a) == 1 and ord(a) > 96 and ord(a) < 123

def is_positive_integer(a):
    return isinstance(a,(int,long,Integer)) and a > 0

def lower_upper(a):
    r"""
    EXAMPLES::

        sage: lower_upper('B')
        'b'
        sage: lower_upper('t')
        'T'
    """
    if ord(a) < 91: return a.lower()
    return a.upper()

def neg_operation(a):
    r"""
    EXAMPLES::

        sage: neg_operation(1)
        -1
        sage: neg_operation(12)
        -12
    """
    return -a

class CallableDict(dict):
    def __call__(self, a):
        return self[a]

class InverseFromRank(object):
    def __init__(self, pos, neg):
        self.pos = pos
        self.neg = neg

    def __call__(self, a):
        if a in self.pos:
            return self.neg.unrank(self.pos.rank(a))
        if a in self.neg:
            return self.pos.unrank(self.neg.rank(a))
        raise ValueError("a does not belong to the alphabet")

def fast_inverse_dictionnary(pos,neg):
    inverse = {}
    inverse.update(dict((p,n) for p,n in zip(pos,neg)))
    inverse.update(dict((n,p) for p,n in zip(pos,neg)))
    return inverse

def build_alphabet_with_inverse(data, names=None, name=None, neg=None, inverse=None):
    r"""
    Returns a triple (positive_letters, negative_letters, involution on pos u
    neg).

    EXAMPLES::

        sage: pos,neg,inverse = build_alphabet_with_inverse('abc')
        sage: pos
        {'a', 'b', 'c'}
        sage: neg
        {'A', 'B', 'C'}
        sage: inverse('b')
        'B'
        sage: inverse('C')
        'c'

        sage: pos,neg,inverse=build_alphabet_with_inverse(3,'a')
        sage: pos
        {'a0', 'a1', 'a2'}
        sage: neg
        {'A0', 'A1', 'A2'}
        sage: inverse('a0')
        'A0'
        sage: inverse('A0')
        'a0'
    """
    if data in Sets():
        if neg in Sets():
            assert data.cardinality() == neg.cardinality()
            if inverse is not None:
                return data,neg,inverse
            else:
                return data,neg,InverseFromRank(data,neg)

    if isinstance(data, (int,long,Integer)):
        if names is None:
            from sage.sets.integer_range import IntegerRange
            return IntegerRange(1,data+1), IntegerRange(-data,0), neg_operation
        elif isinstance(names, str):
            if is_lower_character(names):
                pos = TotallyOrderedFiniteSet([names + '%d'%i for i in xrange(data)])
                neg = TotallyOrderedFiniteSet([lower_upper(names) + '%d'%i for i in xrange(data)])
                return pos,neg,fast_inverse_dictionnary(pos,neg)
        raise ValueError("not possible")

    if data is not None:
        data = build_alphabet(data)
    if neg is not None:
        neg = build_alphabet(data)
        assert data.cardinality() == neg.cardinality()
        if data.is_finite():
            return data,neg,fast_inverse_dictionnary(data,neg)
        return data,neg,InverseFromRank(data,neg)

    # now we want to build an inverse
    if not data.is_finite():
        raise ValueError("if ``pos`` is infinite you should provide ``neg``")
    if all(is_lower_character(a) for a in data):
        neg = build_alphabet([a.upper() for a in data])
        inverse = lower_upper
    elif all(is_positive_integer(a) for a in data):
        neg = build_alphabet([-a for a in data])
        inverse = neg_operation
    else:
        raise ValueError("does not know how to build ``neg`` from the given ``pos``")

    return data,neg,inverse



class AlphabetWithInverses(Parent):
    """
    Class for alphabet with inverse letters.

    Intended to be used by FreeGroup.  Builds a finite ordered
    alphabet with an inverse for each letter. There must be no
    duplicate. Inverse letters are either given or assumed to be
    capitalized letters.

    EXAMPLES::

        sage: AlphabetWithInverse(['a','b','c'],['A','B','C'])
        Alphabet with inverse on ['a', 'b', 'c']

        sage: AlphabetWithInverses(3,'a')
        Alphabet with inverses {'a0', 'a1', 'a2'}, {'A0', 'A1', 'A2'}

        sage: AlphabetWithInverses('abcde')
        Alphabet with inverses {'a', 'b', 'c', 'd', 'e'}, {'A', 'B', 'C', 'D', 'E'}

        sage: pos = IntegerRange(1,Infinity)
        sage: neg = IntegerRange(-1,-Infinity,-1)
        sage: AlphabetWithInverses(pos, neg=neg, inverse= lambda x: -x)
        Alphabet with inverses {1, 2, ..}, {-1, -2, ..}
    """
    def __init__(self, data, names=None, name=None, neg=None, inverse=None):
        """
        INPUT:

        - ``pos`` -- the set of positive letters

        - ``neg`` -- the set of negative letters

        - ``inverse`` -- a function that invert a letter
        """
        pos,neg,inverse = build_alphabet_with_inverse(data, names=names, name=name, neg=neg, inverse=inverse)

        self._pos = pos
        self._neg = neg
        self.invert = inverse

    def _repr_(self):
        """
        String representation of self.
        """
        return "Alphabet with inverses %s, %s"%(self._pos,self._neg)

    def __iter__(self):
        """
        Iterator through the letters of the alphabet.

        WARNING:

        The iterator is on all the letters of the alphabet (both
        positive and negative). This is NOT consistent with ```len()``.
        """
        if self._pos.is_finite():
            for p in self._pos: yield pos
            for n in self._neg: yield neg
        itp = iter(self._pos)
        itn = iter(self._neg)
        while True:
            yield next(itp)
            yield next(itn)

    def an_element(self):
        return self._pos.an_element()

    def some_elements(self):
        return tuple(self._pos.some_elements()) + tuple(self._neg.some_elements())

    def cardinality(self):
        """
        The cardinality of the positive letters.
        """
        return self._pos.cardinality() + self._neg.cardinality()

    def __contains__(self,letter):
        """
        Test whether the letter is contained in self
        """
        return letter in self._pos or letter in self._neg

    def __len__(self):
        r"""
        Raise an error if infinite.
        """
        return len(self.cardinality())

    def rank(self,letter):
        """
        Return the rank of the letter

        from 0 to card(self)-1: positive letters
        from card(self) to 2card(self)-1: negative letters
        """
        if letter in self._pos:
            return self._pos.index(letter)
        else:
            return self.cardinality()+self._neg.index(letter)

    def __getitem__(self,n):
        """
        Return the letter with rank n.

        from 0 to card(self)-1: positive letters
        from card(self) to 2card(self)-1: negative letters
        """
        if n < self.cardinality():
            return self._pos[n]
        else:
            return self._neg[n-self.cardinality()]

    def are_inverse(self,a,b):
        """
        Test if the two letters are inverse of each other.
        """
        return self._inverse(a) == b

    def is_positive_letter(self,letter):
        """
        Test if the letter is a positive letter.
        """
        return letter in self._pos

    def is_negative_letter(self,letter):
        """
        Test if the letter is a negative letter.
        """
        return letter in self._neg

    def to_positive_letter(self,letter):
        """
        Given letter a or a^-1 returns a.

        EXAMPLES::

            sage: A = AlphabetWithInverse(['a','b','c'],['A','B','C'])
            sage: A.to_positive_letter('b')
            'b'
            sage: A.to_positive_letter('B')
            'b'
        """
        if letter in self._pos:
            return letter
        elif letter in self._neg:
            return self._inverse(letter)
        else:
           raise ValueError("The letter %s is not in the alphabet %s"%(letter,self))

    def positive_letters(self):
        """
        The set of positive letters.

        EXAMPLES::

            sage: AlphabetWithInverses('abc').positive_letters()
            {'a', 'b', 'c'}
        """
        return self._pos

    def negative_letters(self):
        """
        The set of negative letters.

        EXAMPLES::

            sage: AlphabetWithInverses('abc').negative_letters()
            {'A', 'B', 'C'}
        """
        return self._neg

    def random_element(self, exclude=None):
        """
        A random letter, different from the letters in ``exclude``.
        """
        if exclude is None:
            exclude = []
        while True:
            if random() < .5:
                test = self._pos.random_element()
            else:
                test = self._neg.random_element()
            if test not in exclude:
                return test
