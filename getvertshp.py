import fiona
import numpy as np
from shapely.geometry import mapping, shape
from fiona import collection

def getvert(shp, poly, attr='bacia', buffer=False):

    if buffer:
        with collection(shp+'.shp', "r") as input:
            schema = input.schema.copy()
            with collection(
                    "with-shapely.shp", "w", "ESRI Shapefile", schema
                    ) as output:
                for f in input:
        
                    try:
                        # Make a shapely object from the dict.
                        geom = shape(f['geometry'])
                        geom = geom.buffer(buffer)
                        # Make a dict from the shapely object.
                        f['geometry'] = mapping(geom)
                        output.write(f)
        
                    except Exception as e:
                        # Writing uncleanable features to a different shapefile
                        # is another option.
                        print("Error cleaning feature %s:", f['id'])
    
    if buffer:
        shpe = "with-shapely.shp"
    else:
        shpe = shp+'.shp'
    
    with fiona.open(shpe) as lines:
        print(shpe)
        
        if len(lines) == 1:
            vertices = []
            for line in lines:
                for vert in line['geometry']['coordinates'][0]:
                    vertices.append(vert)
        
            vertices = np.array(vertices)
        
        else:
            for x, poligon in enumerate(lines):
                
                if poligon['properties'][attr] == poly:
                    vert = poligon['geometry']['coordinates']
                    vertices = np.asarray(vert)
                    if len(vertices.shape) < 3:
                        verts2 = []
                        for m in vertices:
                            for x in m:
                                for k in x:
                                    verts2.append(k)
                        vertices = np.asarray(verts2)
            vertices = np.squeeze(vertices)
    
    return vertices