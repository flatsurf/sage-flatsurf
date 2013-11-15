def test1():
    from similarity_surface_generators import SimilaritySurfaceGenerators
    ss=SimilaritySurfaceGenerators.example()
    ss.edit()
    sm,sb = ss.get_bundle()
    sb.set_transform(100,-100,300,300)
    sb.redraw_all()
    hol=(10000000,1000000)
    sb.draw_flow(hol)


def test2():
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    S = SymmetricGroup(3)
    r = S('(1,2)')
    u = S('(1,3)')
    o = translation_surfaces.origami(r,u)
    o.edit()
    se,sb=o.get_bundle()
    sb.set_transform(200,-200,300,300)
    sb.redraw_all()



def test3():
    from similarity_surface_generators import TranslationSurfaceGenerators
    ss=TranslationSurfaceGenerators.regular_octagon()
    ss.edit()
    sm,sb = ss.get_bundle()
    sb.set_transform(100,-100,500,400)
    sb.redraw_all()
    pt=sb.pick_point()
    print str(pt)
    V=ss.vector_space()
    hol=V((100,110))
    segments=pt.flow_segments(hol)
    sb.render_segments(segments)

def test4():
    from similarity_surface_generators import TranslationSurfaceGenerators
    ss=TranslationSurfaceGenerators.infinite_origami_example()
    ss.edit()
    se,sb=ss.get_bundle()
    sb.set_transform(100,-100,300,300)
    sb.redraw_all()

