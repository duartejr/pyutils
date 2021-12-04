#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:35:52 2018

@author: duarte
"""

import ogr
import os

def buffer(infile, outfile, buffdist):

    #try:
    ds_in=ogr.Open( infile )
    lyr_in=ds_in.GetLayer( 0 )
    drv=ds_in.GetDriver()

    if os.path.exists( outfile ):
        drv.DeleteDataSource(outfile)

    ds_out = drv.CreateDataSource( outfile )

    layer = ds_out.CreateLayer( lyr_in.GetLayerDefn().GetName(), \
                                lyr_in.GetSpatialRef(), ogr.wkbPolygon)
    
    for i in range ( lyr_in.GetLayerDefn().GetFieldCount() ):
        field_in = lyr_in.GetLayerDefn().GetFieldDefn( i )
        fielddef = ogr.FieldDefn( field_in.GetName(), field_in.GetType() )   
        layer.CreateField ( fielddef )
            
    for feat in lyr_in:
        geom = feat.GetGeometryRef()
        feature = feat.Clone()
        feature.SetGeometry(geom.Buffer(float(buffdist))) 
        layer.CreateFeature(feature)
        del geom
    ds_out.Destroy()
    
    print(ds_out)