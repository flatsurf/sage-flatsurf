import pytest

from sage.all import ZZ, QQ
import flatsurf

# See Appendix C of McMullen ["Billiard and Teichmüller Curves", 2022]
@pytest.mark.parametrize(
    "discriminant,genus,nu2,nu3,ncusps,chi",
    [
        (5, 0, 0, 1, 1, QQ(-7,15)),
        (8, 0, 1, 1, 2, QQ(-7,6)),
        (12, 0, 1, 0, 3, QQ(-7,3)),
        (13, 0, 0, 2, 3, QQ(-7,3)),
        (17, 0, 0, 1, 6, QQ(-14,3)),
        (20, 0, 2, 1, 5, QQ(-14,3)),
        (21, 1, 0, 1, 4, QQ(-14,3)),
        (24, 1, 2, 0, 6, QQ(-7)),
        (28, 1, 2, 2, 7, QQ(-28,3)),
        (29, 1, 0, 3, 5, QQ(-7)),
        (32, 1, 2, 2, 7, QQ(-28,3)),
        (33, 2, 0, 0, 12, QQ(-14)),
        (37, 1, 0, 4, 9, QQ(-35,3)),
        (40, 2, 2, 2, 12, QQ(-49,3)),
        (41, 3, 0, 1, 14, QQ(-56,3)),
        (44, 3, 4, 2, 9, QQ(-49,3)),
        (45, 4, 0, 0, 8, QQ(-14)),
        (48, 4, 2, 1, 11, QQ(-56,3)),
        (52, 4, 2, 2, 15, QQ(-70,3)),
        (53, 4, 0, 5, 7, QQ(-49,3)),
        (56, 6, 4, 2, 10, QQ(-70,3)),
        (57, 7, 0, 1, 20, QQ(-98,3)),
        (60, 8, 4, 0, 12, QQ(-28))
    ],
)
def test_X_tables(discriminant, chi, genus, ncusps, nu2, nu3):
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
