import pytest

from sage.all import ZZ, QQ
import flatsurf

# See Appendix C of McMullen ["Billiard and Teichmüller Curves", 2022]
@pytest.mark.parametrize(
    "discriminant,genus,nu2,ncusps,chi",
    [
        (5, 0, 1, 1, QQ(-3, 10)),
        (8, 0, 0, 2, QQ(-3, 4)),
        (12, 0, 1, 3, QQ(-3, 2)),
        (13, 0, 1, 3, QQ(-3, 2)),
        (17, 0, 1, 3, QQ(-3, 2)),
        (20, 0, 0, 5, -3),
        (21, 0, 2, 4, -3),
        (24, 0, 1, 6, QQ(-9, 2)),
        (28, 0, 2, 7, -6),
        (29, 0, 3, 5, QQ(-9, 2)),
        (32, 0, 2, 7, -6),
        (33, 0, 1, 6, QQ(-9, 2)),
        (37, 0, 1, 9, QQ(-15, 2)),
        (40, 0, 1, 12, QQ(-21, 2)),
        (41, 0, 2, 7, -6),
        (44, 1, 3, 9, QQ(-21, 2)),
        (45, 1, 2, 8, -9),
        (48, 1, 2, 11, -12),
        (52, 1, 0, 15, -15),
        (53, 2, 3, 7, QQ(-21, 2)),
        (56, 3, 2, 10, -15),
        (57, 1, 1, 10, QQ(-21, 2)),
        (60, 3, 4, 12, -18),
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
