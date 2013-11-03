class AffineSurfacePlot:
    r"""
    Class for various plotting of translation surfaces

    TODO:

    - different output: standard sage plot, tikz, Tk editor
    - class for plotting saddle connections, cylinders, circles.
    """
    def __init__(self, surface, pos=None):
        self._surface = surface
        self._pos = pos

    def _set_positions(self):
        from sage.misc.cachefunc import cached_function

        surface = self._surface
        V = surface.vector_space()
        def rand_pos(x):
            return 12 * V.random_element()
        self._pos = Family(indices=surface.polygons().keys(),
                           cached_function(rand_pos))

    def plot(self, edge_labels=False, polygon_labels=False):
        from sage.plot.graphics import Graphics
        surface = self._surface
        polygons = surface.polygons()

        if not polygons.is_finite():
            raise NotImplementedError("no implementation for surface built from infinitely many polygons")

        if self._pos is None:
            self._set_positions()
        if self._edge_labels is None:
            self._set_edge_labels()

        G = Graphics()
        for a in polygons.keys():
            p = polygons[a]
            v = p.vertices(translation=self._pos[a])
            G += p.plot(translation=self._pos[a])

            if edge_labels:
                for i in xrange(p.num_edges()):
                    G += text(str(self._edge_labels[a,i]),
                        (v[i]+v[(i+1)%p.num_edges()])/2,
                        color='black')

        # then we possibly plot the lable inside each polygon
        if len(self._index_set) != 1:
            for a in polygons.keys():
                p = polygons[a]
                if polygon_labels:
                    m = sum(p.vertices(translation=self._pos[a])) / p.num_edges()
                    G += text(str(a), m, color='black')
        return G



