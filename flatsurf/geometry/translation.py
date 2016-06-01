from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.groups.group import Group
from sage.categories.groups import Groups
from sage.structure.element import MultiplicativeGroupElement
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer import Integer


ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_2 = Integer(2)
ZZ_3 = Integer(3)
ZZ_4 = Integer(4)


class Translation(MultiplicativeGroupElement):
    r"""Class for a translation of the plane. 
    This is an element of a group written multiplicatively (for composition, and
    compatibility with Similarity)."""
    
    def __init__(self, parent, *args):
        r'''
        Construct the translation (x,y) mapsto (x+s,y+t).
        
        Arguments: 
        The first argument must be the parent. 
        If there is one additional argument (call it a) then s=a[0] and t=a[1].
        If there are two elements, these perscribe s and t.
        '''
        if parent is None:
            raise ValueError("The parent must be provided")
        if len(args)==1:
            self._s=args[0][0]
            self._t=args[0][1]
        if len(args)==2:
            self._s=args[0]
            self._t=args[1]
        self._parent=parent
        MultiplicativeGroupElement.__init__(self,parent)

    def _mul_(self,s):
        r'''Compose two similarities.'''
        C = self.__class__
        return C(self._parent,
            s._s+self._s,
            s._t+self._t)

    def __invert__(self):
        r'''Invert a similarity.'''
        C = self.__class__
        return C(self._parent,
            -self._s,
            -self._t)

    def _div_(self,s):
        return self._mul_(s.__invert__())

    def __hash__(self):
        return 13*hash(self._s)+53*hash(self._t)

    def __call__(self,w):
        r'''Return the image of a vector w under the translation.'''
        return vector([ self._s+w[0], self._t+w[1] ])

    def _repr_(self):
        return "Translation by ("+str(self._s)+", "+str(self._t)+")"

    def _cmp_(self, other):
        x=cmp(self._s,other._s)
        if x!=0:
            return x
        return cmp(self._t,other._t)

    __cmp__=_cmp_

    def matrix(self):
        return matrix(self._parent._f,[
            [self._parent._f.one(), self._parent._f.zero(), self._s],
            [self._parent._f.zero(),  self._parent._f.one(), self._t],
            [self._parent._f.zero(), self._parent._f.zero(), self._parent._f.one()]])
            
    def s(self):
        return self._s

    def t(self):
        return self._t

    def sign(self):
        r"""
        Records that this transformation is orientation preserving.
        """
        return 1

class TranslationGroup(UniqueRepresentation,Group):
    r'''Group representing translations in the plane with a multiplicative group operation.
    '''

    Element = Translation

    def _element_constructor_(self, *args, **kwds):
        if len(args)!=1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        return self.element_class(self, x, **kwds)

    def __init__(self, base_field):
        self._f=base_field
        # The vector space of vectors 
        self._vs = VectorSpace(self._f,2)
        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        return "TranslationGroup over field "+str(self._f)

    def one(self):
        return self.element_class(self,self._f.zero(),self._f.zero())

    def an_element(self):
        return self.element_class(self,self._f(ZZ_3),self._f(ZZ_4))

    def is_abelian(self):
        return True

    def gens(self):
        return [
            self.element_class(self._f.one(),self._f.zero()),
            self.element_class(self._f.zero(),self._f.one()) ]
    
    def base_field(self):
        return self._f

