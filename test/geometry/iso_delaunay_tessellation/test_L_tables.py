import pytest

from sage.all import ZZ, QQ
import flatsurf

# See Appendix B of Mukamel ["Orbifold points on Teichmüller curves and Jacobians with complex multiplication", 2013]
@pytest.mark.parametrize("discriminant,chi,genus,ncusps,nu2", [
    (5, QQ(-3/10), 0, 1, 1), #chi_top = (2 - 2 * 0) - 1 = 1?
    (44, QQ(-21/2), 1, 9, 3),
    (53, QQ(-21/2), 2, 7, 3)
])

def test_L_tables(discriminant, chi, genus, ncusps, nu2):
    e = 0 if discriminant % 2 == 0 else -1
    w = (discriminant - e**2)/4
    L = flatsurf.translation_surfaces.mcmullen_genus2_prototype(w, 1, 0, e)
    idt = flatsurf.IsoDelaunayTessellation(L)
    idt.explore()
    idt.plot().show()
    assert idt.topological_euler_characteristic() == (2 - 2 * genus) - ncusps
    assert idt.orbifold_euler_characteristic() == chi
    assert idt.genus() == genus
    assert len(idt.cusps()) == ncusps
    assert len(idt.orbifold_points(2)) == nu2
    assert len(idt.orbifold_points(3)) == 0

    # idt.veech_group() TODO: implement?
    # matching elements + orbifold symmetries = gen. set for Veech group
    # TODO: How do we iterate over IDRs?
    # (R, tau) is an IDR w/ R a hyperbolic polygon and tau a triangulation of M_0
