class SurfacePositioner:
    def __init__(self,field):
        self._positions={}
        self._field=field
        self._zero=self.vector_space().zero()

    def get_position(self,index):
        try:
            pos=self._positions[index]
        except KeyError:
            pos=self._zero
        return pos

    def translate_position(self, index, translation):
        curpos=self.getPosition(index)
        try:
            newpos=curpos+translation
        except StandardError:
            from sage.modules.free_module_element import vector
            try:
                newpos=curpos + vector(translation)
            except StandardError:
                try:
                    newpos=curpos+vector(self._field(translation[0]),self._field(translation[1]))
                except:
                    raise ValueError("can not convert translation to a vector in R^2")
        if newpos==self._zero:
            try:
                del self._positions[index]
            except KeyError:
                pass
        else:
            self._positions[index]=newpos

    def set_position(self, index, position):
        try:
            position=self._zero+position
        except StandardError:
            from sage.modules.free_module_element import vector
            try:
                position=self._zero+vector(position)
            except StandardError:
                try:
                    position=vector(self._field(translation[0]),self._field(translation[1]))
                except:
                    raise ValueError("can not convert translation to a vector in R^2")
        if position==self._zero:
            try:
                del self._positions[index]
            except KeyError:
                pass
        else:
            self._positions[index]=position

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self._field, 2)

