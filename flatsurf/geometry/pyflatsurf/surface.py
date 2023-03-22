from flatsurf.geometry.surface import Surface


class Surface_pyflatsurf(Surface):
    r"""
    EXAMPLES::

        sage: from flatsurf import polygons
        sage: from flatsurf.geometry.surface import Surface_dict

        sage: S = Surface_dict(QQ)
        sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 0), (1, 1)]), label=0)
        0
        sage: S.add_polygon(polygons(vertices=[(0, 0), (1, 1), (0, 1)]), label=1)
        1

        sage: S.set_edge_pairing(0, 0, 1, 1)
        sage: S.set_edge_pairing(0, 1, 1, 2)
        sage: S.set_edge_pairing(0, 2, 1, 0)

        sage: T, _ = S._pyflatsurf()

    TESTS::

        sage: from flatsurf.geometry.pyflatsurf.surface import Surface_pyflatsurf

        sage: isinstance(T, Surface_pyflatsurf)
        True
        sage: TestSuite(T).run()

    """

    def __init__(self, flat_triangulation):
        self._flat_triangulation = flat_triangulation

        # TODO: We have to be smarter about the ring bridge here.
        from flatsurf.geometry.pyflatsurf_conversion import sage_ring

        base_ring = sage_ring(flat_triangulation)

        base_label = map(int, next(iter(flat_triangulation.faces())))

        super().__init__(
            base_ring=base_ring, base_label=base_label, finite=True, mutable=False
        )

    def pyflatsurf(self):
        return self, Deformation_pyflatsurf(self, self)

    @classmethod
    def _from_flatsurf(cls, surface):
        if not surface.is_finite():
            raise ValueError("surface must be finite to convert to pyflatsurf")

        if not surface.is_triangulated():
            raise ValueError("surface must be triangulated to convert to pyflatsurf")

        # populate half edges and vectors
        n = sum(surface.polygon(lab).num_edges() for lab in surface.label_iterator())
        half_edge_labels = {}  # map: (face lab, edge num) in flatsurf -> integer
        vec = []  # vectors
        k = 1  # half edge label in {1, ..., n}
        for t0, t1 in surface.edge_gluing_iterator():
            if t0 in half_edge_labels:
                continue

            half_edge_labels[t0] = k
            half_edge_labels[t1] = -k

            f0, e0 = t0
            p = surface.polygon(f0)
            vec.append(p.edge(e0))

            k += 1

        # compute vertex and face permutations
        vp = [None] * (n + 1)  # vertex permutation
        fp = [None] * (n + 1)  # face permutation
        for t in surface.edge_iterator():
            e = half_edge_labels[t]
            j = (t[1] + 1) % surface.polygon(t[0]).num_edges()
            fp[e] = half_edge_labels[(t[0], j)]
            vp[fp[e]] = -e

        def _cycle_decomposition(p):
            n = len(p) - 1
            assert n % 2 == 0
            cycles = []
            unseen = [True] * (n + 1)
            for i in list(range(-n // 2 + 1, 0)) + list(range(1, n // 2)):
                if unseen[i]:
                    j = i
                    cycle = []
                    while unseen[j]:
                        unseen[j] = False
                        cycle.append(j)
                        j = p[j]
                    cycles.append(cycle)
            return cycles

        # convert the vp permutation into cycles
        verts = _cycle_decomposition(vp)

        # find a finite SageMath base ring that contains all the coordinates
        base_ring = surface.base_ring()

        from sage.all import AA

        if base_ring is AA:
            from sage.rings.qqbar import number_field_elements_from_algebraics
            from itertools import chain

            base_ring = number_field_elements_from_algebraics(
                list(chain(*[list(v) for v in vec])), embedded=True
            )[0]

        from flatsurf.features import pyflatsurf_feature

        pyflatsurf_feature.require()
        from pyflatsurf.vector import Vectors

        V = Vectors(base_ring)
        vec = [V(v).vector for v in vec]

        def _check_data(vp, fp, vec):
            r"""
            Check consistency of data

            vp - vector permutation
            fp - face permutation
            vec - vectors of the flat structure
            """
            assert isinstance(vp, list)
            assert isinstance(fp, list)
            assert isinstance(vec, list)

            n = len(vp) - 1

            assert n % 2 == 0, n
            assert len(fp) == n + 1
            assert len(vec) == n // 2

            assert vp[0] is None
            assert fp[0] is None

            for i in range(1, n // 2 + 1):
                # check fp/vp consistency
                assert fp[-vp[i]] == i, i

                # check that each face is a triangle and that vectors sum up to zero
                j = fp[i]
                k = fp[j]
                assert i != j and i != k and fp[k] == i, (i, j, k)
                vi = vec[i - 1] if i >= 1 else -vec[-i - 1]
                vj = vec[j - 1] if j >= 1 else -vec[-j - 1]
                vk = vec[k - 1] if k >= 1 else -vec[-k - 1]
                v = vi + vj + vk
                assert v.x() == 0, v.x()
                assert v.y() == 0, v.y()

        _check_data(vp, fp, vec)

        from pyflatsurf.factory import make_surface

        return Surface_pyflatsurf(make_surface(verts, vec)), None
