from similarity_surface import *
ss=SimilaritySurfaceGenerators.example()
ss.edit()
sm,sb = ss.get_bundle()
sb.set_transform(100,-100,300,300)
sb.redraw_all()

from sage.groups.perm_gps.permgroup_named import SymmetricGroup
S = SymmetricGroup(3)
r = S('(1,2)')
u = S('(1,3)')
o = translation_surfaces.origami(r,u)
o.edit()
se,sb=o.get_bundle()
sb.set_transform(200,-200,300,300)
sb.redraw_all()



