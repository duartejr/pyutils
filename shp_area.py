#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 21:32:57 2018

@author: duarte
"""
import pyproj
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial


def get_area(geom):
    geom = Polygon(geom.points)

    geom_area = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4674'),
            pyproj.Proj(
                proj='aea',
                lat1=geom.bounds[1],
                lat2=geom.bounds[3])),
        geom)
    # Print the area in m^2
    return geom_area.area
