def straight_line_plot(gs, v, directory="/tmp/", pause_time = 100, width_in_inches=10, file_offset=0):
    r"""
    Produces successively longer plots of a trajectory on the graphical surface gs which starts with the provided
    SimilaritySurfaceTangentVector, v. 
    
    The graphical surface gs is used to decide which segments to draw. Only segments drawn in visible
    polygons will be recorded.
    
    The directory given in the parameter will be filled with some data pertaining to the trajectory.
    A text file named orbit_setup.txt will be produced to record v. Periodically (after drawing pause_time 
    segments), the program will write to a text file, orbitN.txt about the segments and produce orbitN.svg 
    depicting the trajectory.
    
    The program will run until you press Control-C (or until your computer dies from an overheated processor :) ).

    EXAMPLE::
        
        from geometry.chamanara import GraphicalChamanaraSurface
        s=GraphicalChamanaraSurface(QQ(1)/2,8)
        
        from geometry.tangent_bundle import SimilaritySurfaceTangentBundle, SimilaritySurfaceTangentVector
        K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        tb = SimilaritySurfaceTangentBundle(s.get_surface())
        from sage.modules.free_module_element import vector
        v=SimilaritySurfaceTangentVector(tb, 0, vector((0,0)), vector((1,sqrt2)))
        
        from graphical.straight_line_plotter import straight_line_plot
        straight_line_plot(s,v, directory=".", pause_time = 100)
    """
    # Check the directory:
    from os import path
    if not path.isdir(directory):
        ValueError("Provided directory must be a directory")

    # Get a bounding box for the figure:
    xmin,ymin,xmax,ymax=gs.bounding_box()
    width=xmax-xmin
    height=ymax-ymin
    height_in_inches=height*width_in_inches/width
    figsize=[width_in_inches,height_in_inches]
    # Construct a segment using the tangent vector
    from geometry.straight_line_trajectory import SegmentInPolygon
    seg = SegmentInPolygon(v)

    # Used for drawing segments:
    from graphical.straight_line_trajectory import GraphicalSegmentInPolygon

    # Prepare to loop:
    all=0
    count=0

    # For providing human readable points
    from sage.rings.real_mpfr import RR

   # Record the initial conditions to a file:
    f = open(path.join(directory,"orbit_setup.txt"), 'w')
    f.write("v: "+str(v)+"\n")
    f.close()
    
    # Plot the surface
    from sage.plot.graphics import Graphics
    p = gs.plot()
    #p.save(path.join(directory,"orbit0.svg"),xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, \
    #    figsize=figsize,axes=False,aspect_ratio=1,dpi=72)
    plot = p.matplotlib(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, axes=False,aspect_ratio=1)
    plot.set_size_inches(figsize)
    plot.subplots_adjust(left=0,right=1,top=1,bottom=0)
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    plot.set_canvas(FigureCanvasAgg(plot))
    plot.savefig(path.join(directory,"orbit0.svg"))
    
    #plot.savefig(path.join(directory,"orbit0.svg"),bbox_inches=0)
    p=Graphics()
    
    # Move forward through segments. Pay attention when the segment lies within our list of 
    # polygons we care about. Plot them and print some data about them to a file.
    # Every 100 times we do this, we move to a new file.

    while True:
        filename_f = "orbit"+str(count/pause_time+file_offset)+".txt"
        f = open(path.join(directory, filename_f), 'w')
        count_start=count
        while (count==count_start) or (count % pause_time != 0):
            if gs.is_visible(seg.polygon_label()):
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
        #p.save(path.join(directory,"orbit"+str(count/pause_time)+".svg"),xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fig_tight=True, \
        #    figsize=figsize,axes=False,aspect_ratio=1,transparent=True)

        plot = p.matplotlib(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, axes=False,aspect_ratio=1)
        plot.set_size_inches(figsize)
        plot.subplots_adjust(left=0,right=1,top=1,bottom=0)
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        plot.set_canvas(FigureCanvasAgg(plot))
        filename_p = "orbit"+str(count/pause_time+file_offset)+".svg"
        plot.savefig(path.join(directory,"orbit"+str(count/pause_time+file_offset)+".svg"), transparent=True)

        # Write the pile file:
        filename_a = "orbit_pile"+str(count/pause_time+file_offset)+".svg"
        a = open(path.join(directory,filename_a), 'w')
        a.write('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n');
        a.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
        a.write('<svg height="'+str(height_in_inches*72)+'pt" version="1.1" viewBox="0 0 '+ \
            str(width_in_inches*72)+' '+str(height_in_inches*72)+'" width="'+ \
            str(width_in_inches*72)+'pt" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
        for i in range(count/pause_time+1+file_offset):
            a.write('<image x="0" y="0" width="'+str(width_in_inches*72)+'" height="'+str(height_in_inches*72)+\
                '" xlink:href="orbit'+str(i)+'.svg" />\n')
        a.write('</svg>')
        a.close()
        
        p=Graphics()
        print("Wrote to files '"+filename_f+"', '"+filename_p+"' and '"+filename_a+"'.")
