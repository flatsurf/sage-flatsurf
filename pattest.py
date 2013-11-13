from surface_manipulator import *
from similarity_surface import *

def run():
    root = Tk()
    root.geometry("800x700+10+10")
    app = SurfaceManipulator(root)
    pc=PolygonCreator()
    pc.add_vertex((0,0))
    pc.add_vertex((2,-2))
    pc.add_vertex((2,0))
    p0=pc.get_polygon()
    pc=PolygonCreator()
    pc.add_vertex((0,0))
    pc.add_vertex((2,0))
    pc.add_vertex((1,3))
    p1=pc.get_polygon()
    ps=(p0,p1)
    glue={ (0,2):(1,0), (0,0):(1,1), (0,1):(1,2), (1,0):(0,2), (1,1):(0,0), (1,2):(0,1) }
    ss=SimilaritySurface(ps,glue)
    sb=SimilaritySurfaceBundle("Test Similarity Surface", app,ss)
    app.set_surface(sb)
    sb.zoom(100,-1,-4)
    root.mainloop()

run()
