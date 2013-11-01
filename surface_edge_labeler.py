class SurfaceEdgeLabeler:
    def __init__(self):
        self._labels={}
        self._next='A'

    def getLabel(self,index,n):
        try:
            l=self._labels[(index,n)]
        except KeyError:
            l=self._next
            self._labels[(index,n)]=l
            """ 
            The iteration below should be improved so as not to loop!
            I'd like it to switch to lower case, and then use strings of length 2...
            """
            if self._next[-1]=='Z':
                self._next='A'
            else:
                self._next=chr(ord(self._next) + 1)
        return l
