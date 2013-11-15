from surface_manipulator import *
from similarity_surface_generators import *

sm=SurfaceManipulator.launch()

def demo1():
    s=TranslationSurfaceGenerators.octagon_and_squares()
    sm,sb=s.get_bundle()
    sm.set_surface(sb)
    globals()['s'] = s
    globals()['sb'] = sb
    globals()['sm'] = sm


def demo2():
    s=InfiniteStaircase()
    sm,sb=s.get_bundle()
    sm.set_surface(sb)
    globals()['s'] = s
    globals()['sb'] = sb
    globals()['sm'] = sm

def demo3():
    s=SimilaritySurfaceGenerators.right_angle_triangle(3,4)
    sm,sb=s.get_bundle()
    sm.set_surface(sb)
    globals()['s'] = s
    globals()['sb'] = sb
    globals()['sm'] = sm

