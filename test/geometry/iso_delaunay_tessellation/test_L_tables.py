import pytest

from sage.all import ZZ, QQ
import flatsurf

# See Appendix C of McMullen ["Billiard and Teichmüller Curves", 2022]
@pytest.mark.parametrize(
    "discriminant,chi,genus,ncusps,nu2",
    [
        (5, QQ(-3 / 10), 0, 1, 1),  # chi_top = (2 - 2 * 0) - 1 = 1?
        (44, QQ(-21 / 2), 1, 9, 3),
        (53, QQ(-21 / 2), 2, 7, 3),
    ],
)
def test_L_tables(discriminant, chi, genus, ncusps, nu2):
    e = 0 if discriminant % 2 == 0 else -1
    w = (discriminant - e**2) / 4
    L = flatsurf.translation_surfaces.mcmullen_genus2_prototype(w, 1, 0, e)
    idt = flatsurf.IsoDelaunayTessellation(L)
    idt.explore()
    idt.plot().show()
    assert idt.genus() == genus
    assert len(idt.cusps()) == ncusps
    assert idt.orbifold_euler_characteristic() == chi
    assert len(idt.orbifold_points(2)) == nu2
    assert len(idt.orbifold_points(3)) == 0
