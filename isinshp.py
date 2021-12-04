#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:27:09 2021

@author: duarte
"""
import geopandas as gpd
import numpy as np

def isinshp(shape, lat, lon, buffer=0.01):
    gpd_shp = gpd.read_file(shape)
    gpd_geom = gpd_shp['geometry'].buffer(buffer).to_crs('epsg:4674')
    
    xi, yi = np.meshgrid(lon, lat)
    coords = xi.flatten(), yi.flatten()
    coords = np.array(coords).T
    pr_df = gpd.GeoDataFrame(coords,
                             geometry=gpd.points_from_xy(coords[:,0],
                                                         coords[:,1]),
                             crs='epsg:4674')
    
    coords_in = []
    for i in range(len(coords)):
        if gpd_geom.contains(pr_df.geometry[i])[0]:
            coords_in.append(coords[i])
    coords_in = np.array(coords_in)
    
    return coords_in
