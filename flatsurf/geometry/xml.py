r"""
This module contains functions to output surfaces to a string or a file in a human and computer 
readable format using XML. 

We also have a function to recreate the surface from the string/file.

EXAMPLES::

    sage: from flatsurf import *
    sage: from flatsurf.geometry.xml import *
    sage: s=translation_surfaces.square_torus().underlying_surface()
    sage: ss=surface_from_xml_string(surface_to_xml_string(s))
    sage: print(ss==s)
    True
"""

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

def surface_to_xml_string(s, complain=True):
    r""" 
    Convert the surface s into a XML string.
    
    Note that this only encodes the Surface part of the object and not the SimilaritySurface 
    wrapper.
    
    By default complain=True, and it will print warning messages and raise errors if it doesn't 
    think you will be able to reconstruct the surface using surface_from_xml_string.

    Currently s should be Surface_list (or a wrapped Surface_list). If not and complain=True,
    an error message will be printed and the surface will be converted to Surface_list. If 
    complain=False, the surface will be encoded but you may not be able to recover it.
    
    Also the surface should be defined over the rationals. Otherwise a ValueError is raised, which
    can be disabled by setting complain=False.
    """
    if not s.is_finite():
        raise ValueError("Can only xml encode a finite surface.")
    from flatsurf.geometry.similarity_surface import SimilaritySurface
    if isinstance(s,SimilaritySurface):
        s=s.underlying_surface()
    from flatsurf.geometry.surface import Surface_list
    if complain:
        if not isinstance(s,Surface_list):
            # Convert to surface list
            print("Warning:surface_to_xml_string is converting to Surface_list before encoding, "+\
                  " labels will likely be changed. "+\
                  "If this is not desired, "+
                  "call surface_to_xml_string with complain=False.")
            s=Surface_list(surface=s)
    
    from xml.sax.saxutils import escape
    output=[]
    output.append('<surface>')
    
    # Include the type of surface (in case we care about it in the future)
    output.append('<type>')
    output.append(escape(repr(type(s))))
    output.append('</type>')

    from sage.rings.rational_field import QQ
    if complain and s.base_ring()!=QQ:
        raise ValueError("Refusing to encode surface with base_ring!=QQ, "+\
                         "because we can not guarantee we bring the surface back. "+
                         "If you would like to encode it anyway, "+
                         "call surface_to_xml_string with complain=False.")

    output.append('<base_ring>')
    output.append(escape(repr(s.base_ring())))
    output.append('</base_ring>')
    
    output.append('<base_label>')
    output.append(escape(repr(s.base_label())))
    output.append('</base_label>')
    
    output.append('<polygons>')
    for label in s.label_iterator():
        output.append('<polygon>')
        output.append('<label>')
        output.append(escape(repr(label)))
        output.append('</label>')
        p=s.polygon(label)
        for e in range(p.num_edges()):
            output.append('<edge>')
            v=p.edge(e)
            output.append('<x>')
            output.append(escape(repr(v[0])))
            output.append('</x>')
            output.append('<y>')
            output.append(escape(repr(v[1])))
            output.append('</y>')
            output.append('</edge>')
        output.append('</polygon>')
    output.append('</polygons>')

    output.append('<gluings>')
    for l,e in s.edge_iterator():
        ll,ee=s.opposite_edge(l,e)
        if ll<l or (ll==l and ee<e):
            continue
        output.append('<glue>')
        output.append('<l>')
        output.append(escape(repr(l)))
        output.append('</l>')
        output.append('<e>')
        output.append(escape(repr(e)))
        output.append('</e>')
        output.append('<ll>')
        output.append(escape(repr(ll)))
        output.append('</ll>')
        output.append('<ee>')
        output.append(escape(repr(ee)))
        output.append('</ee>')
        output.append('</glue>')
    output.append('</gluings>')
    
    output.append('</surface>')
    return ''.join(output)

def surface_to_xml_file(s,filename,complain=True):
    r"""
    Convert the surface s to a string using surface_to_xml_string,
    then write an XML header and this string to the file with filename provided.
    
    The complain flag is passed to surface_to_xml_string.
    """
    string = surface_to_xml_string(s,complain=complain)
    f = open(filename, 'w')
    f.write('<?xml version="1.0"?>\n')
    f.write(string)
    f.close()

def surface_from_xml_string(string):
    r"""
    Attempts to reconstruct a Surface from a string storing an XML representation of the surface.
    
    Currently, this works only if the surface stored was a Surface_list and the surface was defined
    over QQ.
    """
    import xml.etree.ElementTree as ET
    tree = ET.fromstring(string)
    return _surface_from_ElementTree(tree)

def surface_from_xml_file(filename):
    r"""
    Attempts to reconstruct a Surface from a string storing an XML representation of the surface.
    
    Currently, this works only if the surface stored was a Surface_list and the surface was defined
    over QQ.
    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(filename)
    return _surface_from_ElementTree(tree)

def _surface_from_ElementTree(tree):
    from flatsurf.geometry.surface import Surface_list
    node=tree.find("type")
    if node is None:
        raise ValueError('Failed to find tag named "node"')
    if node.text != repr(Surface_list):
        raise NotImplementedError("Currently can only reconstruct from Surface_list. "+\
                        "Found type "+node.text)

    node=tree.find("base_ring")
    if node is None:
        raise ValueError('Failed to find tag named "base_ring"')
    from sage.rings.rational_field import QQ
    if node.text == repr(QQ):
        base_ring=QQ
    else:
        raise NotImplementedError("Can only reconstruct when base_ring=QQ. "+\
                                  "Found "+tree.find("base_ring").text)

    s=Surface_list(QQ)
    
    # Import the polygons
    polygons=tree.find("polygons")
    if polygons is None:
        raise ValueError('Failed to find tag named "polygons"')
    from flatsurf.geometry.polygon import ConvexPolygons
    P = ConvexPolygons(base_ring)
    for polygon in polygons.findall("polygon"):
        node=polygon.find("label")
        if node is None:
            raise ValueError('Failed to find tag <label> in <polygon>')
        label=int(node.text)
        edges=[]
        for edge in polygon.findall("edge"):
            node=edge.find("x")
            if node is None:
                raise ValueError('Failed to find tag <x> in <edge> in <polygon> with label='+str(label))
            x=base_ring(node.text)
            node=edge.find("y")
            if node is None:
                raise ValueError('Failed to find tag <y> in <edge> in <polygon> with label='+str(label))
            y=base_ring(node.text)
            edges.append((x,y))
        p=P(edges=edges)
        s.add_polygon(p,label=label)
    if s.num_polygons()==0:
        raise ValueError("Failed to add any polygons.")
    
    # Glue the edges
    gluings=tree.find("gluings")
    if gluings is None:
        raise ValueError('Failed to find tag named "gluings"')
    for glue in gluings.findall("glue"):
        node=glue.find("l")
        if node is None:
            raise ValueError('Failed to find tag <l> in <glue>.')
        l=int(node.text)
        node=glue.find("e")
        if node is None:
            raise ValueError('Failed to find tag <e> in <glue>.')
        e=int(node.text)
        node=glue.find("ll")
        if node is None:
            raise ValueError('Failed to find tag <ll> in <glue>.')
        ll=int(node.text)
        node=glue.find("ee")
        if node is None:
            raise ValueError('Failed to find tag <ee> in <glue>.')
        ee=int(node.text)
        s.change_edge_gluing(l,e,ll,ee)

    node=tree.find("base_label")
    if node is None:
        raise ValueError('Failed to find tag named "base_label"')
    base_label = int(node.text)
    s.change_base_label(base_label)
        
    # Return the surface:
    return s

