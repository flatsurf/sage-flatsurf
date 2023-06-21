import pytest

from sage.all import ZZ, QQ
import flatsurf

# See Appendix C of McMullen ["Billiard and Teichmüller Curves", 2022]
@pytest.mark.parametrize(
    "discriminant,genus,nu2,nu3,ncusps,chi",
    [
        (8, 0, 0, 1, 1, QQ(-5,12)),
        (12, 0, 0, 0, 2, QQ(-5,6)),
        (17, 0, 0, 1, 3, QQ(-5,3)),
        (20, 0, 1, 0, 4, QQ(-5,2)),
        (24, 0, 1, 0, 4, QQ(-5,2)),
        (28, 0, 0, 2, 4, QQ(-10,3)),
        (32, 0, 0, 0, 7, QQ(-5)),
        (33, 0, 0, 0, 7, QQ(-5)),
        (40, 0, 1, 2, 6, QQ(-35,6)),
        (41, 0, 0, 1, 8, QQ(-20,3)),
        (44, 0, 1, 2, 6, QQ(-35,6)),
        (48, 1, 0, 0, 10, QQ(-10)),
        (52, 1, 1, 0, 12, QQ(-25,2)),
        (56, 1, 2, 2, 6, QQ(-25,3)),
        (57, 1, 0, 1, 11, QQ(-35,3)),
        (60, 2, 0, 0, 8, QQ(-10))
    ],
)
def test_S_tables(discriminant, chi, genus, ncusps, nu2, nu3):
    e = 0 if discriminant % 2 == 0 else -1
    w = (discriminant - e**2) / 4
    L = flatsurf.translation_surfaces.lanneau_nguyen_genus3_prototype(w, 1, 0, e)
    idt = flatsurf.IsoDelaunayTessellation(L)
    idt.explore()
    idt.plot().show()
    assert idt.genus() == genus
    assert len(idt.cusps()) == ncusps
    assert idt.orbifold_euler_characteristic() == chi
    assert len(idt.orbifold_points(2)) == nu2
    assert len(idt.orbifold_points(3)) == 0
