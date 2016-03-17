def mega_plot(K, x, y, dx, dy, directory="/tmp/", pause_time = 100):
    r"""
    Plots a trajectory on the MegaWollmilchsau which starts at (x,y) in the square and
    moves in direction of the vector (dx,dy). It only plots within three squares
    surrounding the base square. 
    
    The directory given in the parameter will be filled with some data pertaining to the trajectory.
    A text file named orbit_setup.txt will be produced to describe the initial conditions.
    Periodically (after drawing pause_time segments), the program will write to a text file
    orbitN.txt about the segments and produce orbitN.svg depicting the trajectory.
    
    The program will run until you press Control-C (or until the power is cut to your processor :) ).

    EXAMPLE::
    
        sage: from demo.mega_wollmilchsau_plotter import mega_plot
        sage: from sage.calculus.predefined import x
        sage: from sage.rings.integer_ring import ZZ
        sage: from sage.rings.number_field.number_field import NumberField
        sage: K=NumberField(x**2-2,'s',embedding=ZZ(1))
        sage: sqrt2=K.gens()[0]
        sage: x=K.zero()
        sage: y=K.zero()
        sage: dx=K.one()
        sage: dy=sqrt2
        sage: mega_plot(K,x,y,dx,dy, directory="/tmp/", pause_time = 10)
        Segment #0 (plotted 0)
        Segment #1 (plotted 1)...
    """
    # Check the directory:
    from os import path
    if not path.isdir(directory):
        ValueError("Provided directory must be a directory")

    # Construct the surface:
    from geometry.mega_wollmilchsau import MegaWollmilchsau
    s=MegaWollmilchsau()

    # Get the labels for 3 polygons.
    l0=s.base_label()
    l1=s.opposite_edge(l0,1)[0]
    l2=s.opposite_edge(l0,2)[0]
    labels=[l0,l1,l2]

    # Setup the portion of the surface we want to draw
    from graphical.surface import GraphicalSurface
    gs=GraphicalSurface(s)
    gs.make_adjacent_and_visible(l0,1)
    gs.make_adjacent_and_visible(l0,2)

    # Construct a tangent vector
    from geometry.tangent_bundle import SimilaritySurfaceTangentBundle, SimilaritySurfaceTangentVector
    tb = SimilaritySurfaceTangentBundle(s)
    from sage.modules.free_module import VectorSpace
    V=VectorSpace(K,2)
    v=SimilaritySurfaceTangentVector(tb, l0, V((x,y)), V((dx,dy)))

    # Construct a segment using the tangent vector
    from geometry.straight_line_trajectory import SegmentInPolygon
    seg = SegmentInPolygon(v)

    # Plot the surface
    p = gs.plot()

    # Used for drawing segments:
    from graphical.straight_line_trajectory import GraphicalSegmentInPolygon

    # Prepare to loop:
    all=0
    count=0

    # For providing human readable points
    from sage.rings.real_mpfr import RR

    # Record the initial conditions to a file:
    f = open(path.join(directory,"orbit_setup.txt"), 'w')
    f.write("Field: "+str(K)+"\n")
    f.write("x: "+str(x)+"\n")
    f.write("y: "+str(y)+"\n")
    f.write("dx: "+str(dx)+"\n")
    f.write("dy: "+str(dy)+"\n")
    f.close()


    # Move forward through segments. Pay attention when the segment lies within our list of 
    # polygons we care about. Plot them and print some data about them to a file.
    # Every 100 times we do this, we move to a new file.

    while True:
        f = open(path.join(directory, "orbit"+str(count/pause_time)+".txt"), 'w')
        count_start=count
        while (count==count_start) or (count % pause_time != 0):
            if seg.polygon_label() in labels:
                print("Segment #"+str(all)+" (plotted "+str(count)+")")
                f.write("Segment #"+str(all)+" (plotted "+str(count)+")\n")
                f.write("Label: "+str(seg.polygon_label())+"\n")
                f.write("point: "+str(seg.start_point())+"\n")
                f.write("RR point: "+str(RR(seg.start_point()[0]))+", "+\
                    str(RR(seg.start_point()[1]))+"\n")
                gseg = GraphicalSegmentInPolygon(gs, seg)
                p += gseg.plot()
                count += 1
            # Move to next segment under flow:
            seg = seg.next()
            all += 1
        # Save the plot to a file.
        f.close()
        p.save(path.join(directory,"orbit"+str(count/pause_time-1)+".svg"))
        print("Wrote to files named 'orbit"+str(count/pause_time-1)+"' with extension .txt and .svg")



