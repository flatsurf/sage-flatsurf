r"""
Class to store and manipulate edge identifications between polygons.
"""

class EditableEdgeGluing:
    def __init__(self):
        self._gluings={}
        self._labels={}
        self._label_to_edge={}

    def get_edge(self,label):
        return self._label_to_edge.get(label)

    def get_opposite(self, p1, e1):
        return self._gluings.get( (p1,e1) )

    def remove_label(self, label):
        if self._label_to_edge.has_key(label):
            pair1=self.get_edge(label)
            pair2=self.get_opposite(pair1[0], pair1[1])
            del self._labels[pair1]
            del self._labels[pair2]
            del self._label_to_edge[label]

    def remove_glue(self, p1, e1):
        pair2=self.get_opposite(p1,e1)
        if pair2 is None:
            return
        label=self.get_label(p1,e1)
        if label is not None:
            self.remove_label(label)
        del self._gluings[(p1,e1)]
        del self._gluings[pair2]

    def glue(self, p1, e1, p2, e2):
        if (p1 != p2) or (e1 != e2):
            self.remove_glue(p1,e1)
            self.remove_glue(p2,e2)
            self._gluings[(p1,e1)]=(p2,e2)
            self._gluings[(p2,e2)]=(p1,e1)
        else:
            raise ValueError("The edges being identified must be distinct!")

    def get_label(self,p1,e1):
        return self._labels.get( (p1,e1) )

    def get_unused_label(self):
        start=ord('a')
        for i in range(26):
            label=chr(start + i)
            if not self._label_to_edge.has_key(label):
                return label
        start=ord('A')
        for i in range(26):
            label=chr(start + i)
            if not self._label_to_edge.has_key(label):
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

    def glue_and_label(self, p1, e1, p2, e2, label=None):
        self.glue(p1, e1, p2, e2)
        self.add_label(p1, e1, label=label)


